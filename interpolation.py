import numpy as np
import scipy
import scipy.interpolate

from . import useful
import importlib
importlib.reload(useful)
from .useful import *

# Strategy for interpolation with optimal noise reduction when output has an equally spaced domain:
# 
# <ol>
#     <li>Step along each domain point.
#     <li>Count thern number of points on the left and right sides of the window.
#     <li>If left >= 1 and right >= 1 and left+right >= 3, continue.
#     <li>Otherwise, expand window until conditions are met. Beware of window functions that go to zero at the edges.
#     <li>Take the fit polynomial about that point.
# </ol>

# In[42]:
def auto_decimate(x, dx):
    tf = np.ones(len(x), dtype=bool)
    xlow = x[0]
    for i in range(1,len(x)-1):
        # if our absence would not cause the resulting interval to fall below dx,
        # mark for deletion
        if abs(x[i+1]-xlow) < dx:
            tf[i] = False
        else:
            xlow = x[i]
    return tf


# Bracket the indices of an array surrounding value limits.
class BracketLimitsError(Exception):
    """Basic exception when desired value cannot be bracketed by upper/lower indices"""
    def __init__(self, description=None):
        self.description = description
    def __str__(self):
        return self.description

LowerLimit = BracketLimitsError('LowerLimit')
BelowLowerLimit = BracketLimitsError('BelowLowerLimit')
UpperLimit = BracketLimitsError('UpperLimit')
AboveUpperLimit = BracketLimitsError('AboveUpperLimit')


def search_back_strict(A, y, i):
    while A[i-1] >= y:
        i -= 1
    return i
def search_back_inclusive(A, y, i):
    while A[i] >= y:
        i -= 1
    return i
def search_forward_strict(A, y, i):
    while A[i+1] <= y:
        i += 1
    return i
def search_forward_inclusive(A, y, i):
    while A[i] <= y:
        i += 1
    return i

# a = np.asarray([1,2,3,4,5])
# b = np.asarray([2,3,4,5,6])
# returns limits corresponding to the subarrays [2,3,4,5]
def domain_intersection_limits(a,b):
    a_ll = search_forward_strict(a,b[0],0)
    b_ll = search_forward_strict(b,a[0],0)
    a_ul = search_back_strict(a,b[len(b)-1],len(a)-1)
    b_ul = search_back_strict(b,a[len(a)-1],len(b)-1)
    return a_ll,a_ul+1,b_ll,b_ul+1

# inclusive: get indeces completely surrounding the desired interval (returns one index past)
# allow_endpoints: in inclusive mode, allow an endpoint to be returned as a bracketing value, even if it is on the interval boundary and not beyond it.
# allow_all: do not throw any boundary exceptions. just return boundary value.
# NOTE for numpy users: to select a subarray from the return value, use (lower, upper+1)
def bracket_interval(A, lower, upper, start=None, inclusive=None, allow_endpoints=None, allow_all=None):
    # coerce starting value
    if start is None or start < 0:
        start = 0
    if start >= len(A):
        start = len(A)-1
    
    # defaults
    if inclusive is None:
        inclusive = False
    
    # defaults
    if allow_endpoints is None:
        allow_endpoints = False
        
    # defaults
    if allow_all is None:
        allow_all = False


    # to get the first of a repeated element, search strict
    # to go one past all repeated elements, search inclusive
    #     if you get an index error, consider the last two elements exactly
    # to get the last of a repeated element, search inclusive and subtract one
    #     if you get an index error, consider the last two elements exactly

    if inclusive: # inclusive
        # search forward, backward, then forward.
        try:
            a = search_forward_inclusive(A,lower,start)
        except IndexError:
            # we overran, but last value might be in requested interval
            a = len(A)-1
        try:
            b = search_back_inclusive(A,lower,a)
        except IndexError:
            # overran start of array. ok if allow endpoints.
            if allow_endpoints:
                b = 0
            elif A[0] == lower:
                if allow_all:
                    b = 0
                else:
                    raise LowerLimit
            else:
                if allow_all:
                    b = 0
                else:
                    raise BelowLowerLimit
        try:
            c = search_forward_inclusive(A,upper,b)
        except IndexError:
            if allow_endpoints:
                c = len(A)-1
            elif A[len(A)-1] == upper:
                if allow_all:
                    c = len(A)-1
                else:
                    raise UpperLimit
            else:
                if allow_all:
                    c = len(A)-1
                else:
                    raise AboveUpperLimit
                
        # we have a bracketing strict interval, but in the case of a repeated value, our interval
        # may only include the first or last repeated entry. For fitting, these identical points 
        # provide additional data that should be considered. Endpoints are always allowed in this case.
        try:
            b = search_back_strict(A,A[b],b)
        except IndexError:
            b = 0
        try:
            c = search_forward_strict(A,A[c],c)
        except IndexError:
            c = len(A)-1

        return (b,c)
    else:
        # search forward (inclusive), backward (strict), then forward (inclusive), then backward (strict)
        # to bracket entire set of potentially identical values
        try:
            a = search_forward_inclusive(A,lower,start)
        except IndexError:
            # we overran, but last value might be in requested interval
            a = len(A)-1
        try:
            b = search_back_strict(A,lower,a)
        except IndexError:
            # overran start of array. ok if it matches value.
            if A[0] == lower:
                b = 0
            else:
                if allow_all:
                    b = 0
                else:
                    raise BelowLowerLimit
        try:
            c = search_forward_strict(A,upper,b)
        except IndexError:
            # we overran, but last value might be in requested interval
            if A[len(A)-1] == upper:
                c = len(A)-1
            else:
                if allow_all:
                    c = len(A)-1
                else:
                    raise AboveUpperLimit

        return (b,c)

        

#A = [1,1,1,2,5,7,9,11,12,13,13]
#try:
#    print bracket_interval(A,1.5, 12.5 ,start=-1, inclusive=False, allow_endpoints=False)
#except BracketLimitsError as e:
#    print e


# <h1>Interpolant point picking strategy</h1>
# <ul>
#     <li>Bracket control points at desired sizes and intervals (<tt>dx</tt> and <tt>fixed_x=0.0</tt>)
#     <li>Do we expand interval if we don't have enough data or just duplicate existing control points? (<tt>expand_over_dead_zone=False</tt>)
#     <li>Do we remove control points that are too close together? (<tt>min_control_point_density=0.0</tt>)
#     <li>Fitting window on interval [-0.5:0.5] = <code>fit_window = lambda x: 0.54 + 0.46 \* np.cos(2 \* np.pi \* x)</code>
#     <li>Fit each interval.
# </ul>

# In[16]:


#X=np.asarray([ 1.75 , 3. ,   4.25 , 6.75 , 8.  ])
#print X[-1]
#dx= 1.25
#fixed_x= 0.5
#x_min = int((X[0]-fixed_x)/dx)-1
#x_max = int((X[-1]-fixed_x)/dx)+1
#print range(x_max-x_min+1)
#
#values_at_interval([0,1,2.2],0.25,0.25)
#
#values_at_interval(np.asarray([ 1.75 , 3. ,   4.25 , 6.75 , 8.  ]) , dx= 1.25 , fixed_x= 0.5 )
        
        


# In[17]:


# global default; may be changed by user for preview-style calculations
fast_mode = False

# dx and fixed_x define the control points. do they also define our desired output data?
def interpolate_smoothly(X, Y, dx=None, fixed_x=None, window_size=None, fit_window=None, poly_order=None, expand_over_gap=None, min_control_point_spacing=None, return_cubic_spline=None, fast=None):
    if fast is None:
      fast = fast_mode

    XI = np.argsort(X)
    #print "sorting table:",XI
    #print "len = ",len(XI)
    X = shuffle_array(np.asarray(X), XI)
    Y = shuffle_array(np.asarray(Y), XI)
    
    minX = X[0] #min(X)
    maxX = X[-1] #max(X)
    
    if fixed_x is None:
        fixed_x = 0
    if dx is None:
        # Default point density is equal to original point density 
        dx = (maxX-minX)/(len(X)-1)
    if window_size is None:
        # Default window size is twice the point density
        window_size = dx*2
    if fit_window is None:
        # Hamming window weights
        fit_window = lambda x: 0.54 + 0.46*np.cos(2*np.pi*x)
    if poly_order is None:
        poly_order = 2 # linear
    if expand_over_gap is None:
        expand_over_gap = True
        
    # default minimum spacing is min(dx, min(X spacing))/2
    if min_control_point_spacing is None:
        xprev=X[0]
        min_control_point_spacing = dx
        for x in X[1:]:
            min_control_point_spacing = min(abs(x-xprev),min_control_point_spacing)
        min_control_point_spacing /= 2
    
    if return_cubic_spline is None:
        return_cubic_spline = False
    
    # fast mode for preview
    if fast:
      # first remove duplicate or closely spaced entries
      tf = auto_decimate(X,dx/4)
      csx = list(X[tf])
      csy = list(Y[tf])

    # main routine for real data
    else:
      # Pick a set of control points
      control_points = values_at_interval(X, dx, fixed_x)
      fit_domains = [(x-window_size/2.0,x+window_size/2.0) for x in control_points]
      lost_control_points = [] # starts empty. may end up being unused
      fit_sets_X = []
      fit_sets_Y = []
      fit_sets_W = []
      exact_cs_points = set() # put points we wish to include in cubic spline explicitly
      xy_points_fitted = set()
      
      #print "Chose control points and fit domains:",zip(control_points,fit_domains)
      
      # speed up bracket_interval by keeping track of approximate operating position within data array
      rough_data_index = 0
      
      # do it in reversed order so we can remove points from control_points if we can't use them.
      for i in reversed(list(range(len(control_points)))):
          #print ""
          #print "control_points[",i,"] = ",control_points[i]
                        
          # need a try block here if control point could possibly be outside the range of X
          # count data points covered by window
          #print "Fit domain is",(fit_domains[i][0], fit_domains[i][1])
          (a,b) = bracket_interval(X, fit_domains[i][0], fit_domains[i][1], start=rough_data_index, inclusive=False, allow_all=True)
          points = b-a+1
          #print "Bracket interval is",(a,b),"; points=",points
                        
          # if not enough:
          if poly_order >= points:
              # expand_over_gap is now misnamed. this might become default behavior with no need for
              # a variable name
              if expand_over_gap:
                  
                  # expand window to get enough control points
                  # We need a sane way to symmetrically (except near endpoints) expand the window until we have
                  # enough data to perform the fit.
                  #
                  # Here's a possible strategy:
                  # Call bracket_interval() again, but with inclusive=True. That will get us at least one more index,
                  # but we can allow_all so nothing overflows an endpoint. Since we want to expand symmetrically,
                  # take the largest distance from control_points[i] as half_interval, and call bracket_interval()
                  # yet again with inclusive=False and a=control_points[i]-half_interval,
                  # b=control_points[i]+half_interval and allow_all=True. Repeat until poly_order +1 > points.
                  #
                  #
                  # Include data points explicitly. Duplicates will be removed later.
                  #print "not enough points to fit"
                  for x,y in zip(X[a:b+1],Y[a:b+1]):
                      exact_cs_points.add( (x,y) )
                      #print "Adding to exact_cs_points: ", (x,y)
                  #print "Deleting control point",control_points[i]
                  lost_control_points += [control_points[i]]
                  del(control_points[i])
                  
                  # It would be nice to mark points we have already used for fitting and not use them in 
                  # exact_cs_points
                  continue
              else:
                  print("Not enough points to fit. Don't use expand_over_gap=False")
                  # add to lost_control_points and skip the rest of the loop. Can we actually do anything
                  # with the lost control points? There are gaps in the data here. Might have to discard.
                  lost_control_points += [control_points[i]]
                  del(control_points[i])
                  continue
          # Add subarray on interval (a,b) (inclusive) to lists of x and y arrays for fitting.
          #print "adding to fit set:",(X[a:b+1]-control_points[i],Y[a:b+1],fit_window((X[a:b+1]-control_points[i])/(fit_domains[i][1]- fit_domains[i][0])))
          fit_sets_X.append(X[a:b+1]-control_points[i])
          fit_sets_Y.append(Y[a:b+1])
          fit_sets_W.append(fit_window((X[a:b+1]-control_points[i])/(fit_domains[i][1]- fit_domains[i][0])))
          #print "Computing fit window:"
          #print "X[a:b+1]-control_points[i]=",X[a:b+1]-control_points[i]
          #print "(fit_domains[i][1]- fit_domains[i][0])=",(fit_domains[i][1]- fit_domains[i][0])
          
          # now add the data points we include in the fit into the set xy_points_fitted:
          for xy in zip(X[a:b+1],Y[a:b+1]):
              xy_points_fitted.add(xy)
          #print "xy_points_fitted now contains:",xy_points_fitted
              
      #print ""
      #print "Do the fitting."
      
      #print "control_points=",control_points
      
      csx = reversed(control_points)
      csy = []
      for x,fitX,fitY,fitW in zip(reversed(control_points),fit_sets_X,fit_sets_Y,fit_sets_W):
          #print "Fit point:",x
          #print "X:",fitX,"Y:",fitY,"W:",fitW
          
          # currently not using this bit of kludge code
          if False:
              # if we need an exact solution, force the fitter to find it by adding an irrelevant point to the fit.
              if poly_order == len(fitX):
                  fitX = np.append(fitX, 0.0)
                  fitY = np.append(fitY, sum(fitY)/float(len(fitY)) )
                  fitW = np.append(fitW, min(fitW)*1e-6)
                  #print "Adjusted fit vectors for an exact solution: fitX=",fitX,"fitY=",fitY,"fitW=",fitW
              
          # do the weighted poly fit on its domain
          p=np.polyfit(fitX, fitY, 1, w=fitW)
          csy.append(np.polyval(p,0.0))
          # add to list of cubic splines
          #print "  result:",np.polyval(p,0.0)

      # build the control point arrays
      csx = control_points
      csy = list(reversed(csy))
      
      # zip lists so we can check for duplicates
      cs = list(zip(csx,csy))
          
      #print "exact_cs_points=",exact_cs_points
      #print "xy_points_fitted=",xy_points_fitted
      #print "set(cs)=",set(cs)
      
      # For the purposes of excluding fitted points from being used as exact_cs_points, explicit endpoints
      # are always allowed if the fit failed there.
      xy_points_fitted.discard( (X[0],Y[0]) )
      xy_points_fitted.discard( (X[-1],Y[-1]) )

          
      # Include exact spline handles (except for those that have been used as fit points). Remove 
      # duplicate entries by converting to a set, then unzipping. Note that zip(*x) is inverse of zip.
      points =  list(zip(*list(set(cs).union(exact_cs_points-xy_points_fitted))))
      #print "Duplicates removed. zip(exact_cs_points)=",points
      
      csx = list(points[0])
      csy = list(points[1])

      # we now have points out of order, so sort
      csxi = np.argsort(csx)
      csx = shuffle_list(csx,csxi)
      csy = shuffle_list(csy,csxi)

    # fast mode and high quality mode rejoin here
    if min_control_point_spacing > 0.0:
        # Ideal method:
        # find spacings between each pair of control points
        # throw out the point spaced closest to its peers
        # adjust lists and repeat process until no pair is closer than min_control_point_spacing
        #
        # Implemented method: delete from end to beginning if spacings are too close. average in a 
        # way that becomes crude if we end up deleting more points.
        for i in reversed(list(range(0,len(csx)-1))):
            if abs(csx[i]-csx[i+1]) < min_control_point_spacing:
                csy[i+1] = (csy[i+1]+csy[i])/2
                del(csx[i])
                del(csy[i])
    
    csx = np.asarray(csx)
    csy = np.asarray(csy)
    
    #print ""
    #print "doing spline on csx and csy:"
    #print "csx=",csx
    #print "csy=",csy
    #print "dtype of csy=",csy.dtype
    
    
    if np.any(np.iscomplex(csy)):
      # do cubic spline; default smoothing s is not zero, so must pass explicitly
      iscomplex = True
      cs_real = scipy.interpolate.splrep(csx, np.real(csy), s=0)
      cs_imag = scipy.interpolate.splrep(csx, np.imag(csy), s=0)
    else:
      iscomplex = False
      try:
        cs = scipy.interpolate.splrep(csx, csy, s=0)
      except TypeError:
        f_interp = scipy.interpolate.interp1d(csx, csy, kind='linear')
        print 'Falling back to linear interpolation.'


    if return_cubic_spline:
      if iscomplex:
        return cs_real,cs_imag
      else:
        return cs

    else:
        # figure out the domain. start at fixed_x and proceed backward or forward by dx steps until
        # we find the minimum and maximum possible values.
        
        # find largest domain of fixed_x + i*dx that is supported by data
        xnew = np.asarray( values_at_interval(csx, dx, fixed_x) )
        #print "interpolate_smoothly() csx: {} {}".format(min(csx),max(csx))
        #print "interpolate_smoothly() from values_at_interval: {} {}".format(min(xnew),max(xnew))
        
        # if cs doesn't exist, try the complex version
        if iscomplex:
          ynew_real = scipy.interpolate.splev(xnew, cs_real, der=0)
          ynew_imag = scipy.interpolate.splev(xnew, cs_imag, der=0)
          ynew = ynew_real + 1.j*ynew_imag
        else:
          try:
            ynew = scipy.interpolate.splev(xnew, cs, der=0)
          except:
            ynew = f_interp(xnew)

        # return result
        return xnew,ynew

if __name__ == '__main__':
  import matplotlib.pyplot as plt

  x = [0,1,2,3,4,5,6,7,8,9]
  y = [0,1,4,9,16,25,36,49,64,81]

  x2,y2 = interpolate_smoothly(x,y,dx=0.25)
  plt.scatter(x,y)
  plt.plot(x2,y2)
  plt.show()

# splines X-Y data along X, keeps domain (X-values) that are in-common
# performs operation on Y-arrays
#
# if XX is a list of x-arrays
#  and YY is a list of y-arrays:
# zipXXYY = zip(XX, YY)
#
# returns the common domain X and the reduced range data Y
def common_domain_reduce(zipXXYY, operation, dx=None):
  x_data, y_data = zip(*zipXXYY)
  if dx is None:
      # get dx from first dataset
      dx = (x_data[0][-1] - x_data[0][0])/(len(x_data[0])-1)

  # the notation "sum" is used here, but the operation is performed generally
  x_sum,y_sum = interpolate_smoothly(x_data[0],y_data[0],dx)
  y_sum[:] = 0 

  for x,y in zip(x_data,y_data): 
      x,y = interpolate_smoothly(x,y,dx)
      
      # find out the respective limits of the x 
      a_ll,a_ul,b_ll,b_ul = domain_intersection_limits(x,x_sum) 
      # subset x 
      x_sum = x_sum[b_ll:b_ul] 
      # subset and sum y_data 
      y_sum = operation(y_sum[b_ll:b_ul], y[a_ll:a_ul])
  return x_sum,y_sum

# if XX is a list of x-arrays
#  and YY is a list of y-arrays:
# zipXXYY = zip(XX, YY)
def average_xy_datasets(zipXXYY, dx=None):
  x_sum, y_sum = common_domain_reduce(zipXXYY, np.add, dx=dx)
  return x_sum, y_sum/len(zipXXYY)
