import numpy as np

# shuffle array a using the indices given in indices
def shuffle_array(a, indices):
    if len(a) != len(indices):
        print 'Error: array length mismatch: a='+str(len(a))+', indices='+str(len(indices))
    #print "shuffle array. a =",a
    c = np.zeros(len(indices),dtype=a.dtype)
    for i in range(len(indices)):
        #print "shuffle array. i =",i
        #print "  a[indices[i]]=",a[indices[i]]
        c[i]= a[indices[i]]  
    return c
def shuffle_list(a, indices):
    if len(a) != len(indices):
        print 'Error: array length mismatch: a='+str(len(a))+', indices='+str(len(indices))
    #print "shuffle array. a =",a
    c = [0.0]*(len(indices))
    for i in range(len(indices)):
        c[i]= a[indices[i]]  
    return c
# how to use:
#XI = np.argsort(X)
#X = shuffle_array(X, XI)
#Y = shuffle_array(Y, XI)

def analytic_signal(x):
    from scipy.fftpack import fft,ifft
    N = len(x)
    X = fft(x,N)
    h = np.zeros(N)
    h[0] = 1
    h[1:N//2] = 2*np.ones(N//2-1)
    h[N//2] = 1
    Z = X*h
    z = ifft(Z,N)
    return z

# values_at_interval works a bit like np.linspace, but instead of using
# and explicit start point, we calculate the start point from available
# data. Furthermore, that start point must be a perfect multiple of 
# dx that is offset by fixed_x. So several datasets on different domains
# can be lined up, and all the points will match.
#
# X must be a pre-sorted iterable representing the full domain of the data
# dx is the resulting spacing of returned points with one point fixed at fixed_x
def values_at_interval(X, dx, fixed_x):
    #print "values_at_interval"
    #print "values_at_interval(X=from {} to {}".format(min(X),max(X)), "dx=",dx,", fixed_x=",fixed_x,")"
    # find largest domain of fixed_x + i*dx that is supported by data
    x_min = int((X[0]-fixed_x)/dx)-1
    x_max = int((X[-1]-fixed_x)/dx)+1
    #print "x_min=",x_min
    #print "x_max=",x_max
    x_list = [(x+x_min)*dx+fixed_x for x in range(x_max-x_min+1)]
    #print "min/max (x_list)= {} {}".format(min(x_list),max(x_list))
    while x_list[0] < X[0]:
        #print "deleted point from x_list[0]"
        del x_list[0]
    while x_list[-1] > X[-1]:
        #print "deleted point from x_list[-1]"
        del x_list[-1]
    #print "before returning from values_at_interval, min/max (x_list)= {} {}".format(min(x_list),max(x_list))
    return x_list

def make_array_if_not_already(a,length):
    try:
        if len(a) != length:
            a = np.ones(length)*a[0]
    except TypeError:
        a = np.ones(length)*a
    return a


