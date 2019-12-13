import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

# Rational function aproximation
# Nat Fortune and Scott Hannahs call this Pade approximation, but that's a specific
# procedure to transform a Taylor Series into a rational fraction
import scipy.optimize

# originally from https://stackoverflow.com/questions/29815094/rational-function-curve-fitting-in-python
# but had to patch bug
def rational(x, p, q):
    """
    The general rational function description.
    p is a list with the polynomial coefficients in the numerator
    q is a list with the polynomial coefficients (except the first one)
    in the denominator
    The first coefficient of the denominator polynomial is fixed at 1.
    """
    #r= np.polyval(p, x) / np.polyval(list(q) + [1.0],x)
    #print "rational() with x=",x,"p=",p,"q=",q," = ","n=",n,", d=",d
    #return r
    return np.polyval(p, x) / np.polyval(list(q) + [1.0],x)

# Here's some Python trickery. Original idea from 
# https://stackoverflow.com/questions/10250461/passing-additional-arguments-using-scipy-optimize-curve-fit
#
# We want to pass some fixed parameters into our rational function x-evaluator. 
# However, unlike scipy.optimize.minimize, there is no args pass-through variable.
# So we have to create one by defining a function that returns the objective
# function with the constants baked in.
#
# Here's how we define exactly what rational function to use:
# args:
#    'len_p' : len(p), where len(coeffs) - len_p = len(q)
#    'p0' : value to fix constant for p0. If not present, fit p0.
#
#    Note that len_p must account for any fixed p0. That is, in order to fit
#    (p2*x**2 + p1*x + p0)/(q1*x + q0) with p0=0.0 fixed, {'len_p':3,'p0':0.0}
def rational_objective_func_builder(args):
    len_p = args['len_p']
    if 'p0' in args:
        p0 = args['p0']
        def rational_objective_func(x,*coeffs):
            coeffs = list(coeffs)
            p = coeffs[:len_p]+[p0]
            q = coeffs[len_p:]
            return rational(x,p,q)
    else:
        def rational_objective_func(x,*coeffs):
            coeffs = list(coeffs)
            p = coeffs[:len_p]
            q = coeffs[len_p:]
            return rational(x,p,q)
    
    return rational_objective_func

# config is a dictionary that looks like: {'len_p':3,'len_q':1,'p0':1.0}
# p0 is an optional parameter fixing the x=0 point if desired
#
# p_guess is an optional starting guess in the form (p,q) where p and q are the
# coefficient arrays of the numerator and denominator, respectively.
def fit_rational(x,y,config,p_guess=None):
    if p_guess is None:
        p_guess=[0.0]*(config['len_p']+config['len_q'])
    else:
        p_guess = list(p_guess[0]) + list(p_guess[1])
    popt, pcov = scipy.optimize.curve_fit(rational_objective_func_builder(config), x, y, p0=p_guess, method='dogbox')
    print(popt)
    p,q = popt[:config['len_p']],popt[config['len_p']:]
    if 'p0' in config:
        p = list(p) + [config['p0']]
    return list(p),list(q)

# demo
if False:
    x = np.arange(0,1.5,0.003)
    #y = (3*x**5-x**2-x+3.0)/(1+3.0*x)
    y = np.sqrt(8*x)+1

    rational_config = {'len_p':3,'len_q':1,'p0':1.0}
    p,q = fit_rational(x,y,rational_config)
    plt.plot(x, y, label='original')
    plt.plot(x, rational(x,p,q), label='rational fit')

    d = (min(x),max(x))
    x2 = ((x-d[0])/(d[1]-d[0])-0.5)*2

    cheb = np.polynomial.chebyshev.chebfit(x2,y,4)
    plt.plot(x,np.polynomial.chebyshev.chebval(x2, cheb), label='cheb fit' )

    #plt.xlim((-0.01,0.1))
    plt.legend()
    plt.show()


# list of float. necessary for serializing data that might be a python list or a numpy array
def lof(a):
  return [float(x) for x in a]


# Chebyshev approximation tools for use with thermometry calibration

def to_unity_range(x,d):
    return ((x-d[0])/(d[1]-d[0])-0.5)*2
def from_unity_range(x,d):
    return ((x+1)/2.0*(d[1]-d[0])) + d[0]

class ScaledChebyshev:
    def __init__(self, cheb_order=0):
        self.cheb_order = int(cheb_order) # order of the Chebyshev polynomial mapping a and b
        self.a_limits = None
        self.b_limits = None
        self.a_log = False
        self.b_log = False
        

    def read(self, c):
        self.set_a_limits(c['a_limits'])
        self.set_b_limits(c['b_limits'])
        self.ab_cheb = c['ab_cheb']
        self.ba_cheb = c['ba_cheb']
        if 'a_log' in c:
          self.a_log = c['a_log']
        if 'b_log' in c:
          self.b_log = c['b_log']
    def write(self):
        c = {}
        c['a_limits'] = lof(self.get_a_limits())
        c['b_limits'] = lof(self.get_b_limits())
        c['ab_cheb'] = lof(self.ab_cheb)
        c['ba_cheb'] = lof(self.ba_cheb)
        if self.a_log == True:
          c['a_log'] = True
        if self.b_log == True:
          c['b_log'] = True
        self.cheb_order = len(self.ab_cheb)
        return c
    
    def set_a_limits(self, d):
        if d is None:
            self.a_limits = None
        else:
            self.a_limits = np.log(d) if self.a_log else d
    def set_b_limits(self, d):
        if d is None:
            self.b_limits = None
        else:
            self.b_limits = np.log(d) if self.b_log else d
    def get_a_limits(self):
        if self.a_limits is None:
            return None
        return np.exp(self.a_limits) if self.a_log else self.a_limits
    def get_b_limits(self):
        if self.b_limits is None:
            return None
        return np.exp(self.b_limits) if self.b_log else self.b_limits

            
    def fit(self, a, b, cheb_order=None):
        if cheb_order is None:
          cheb_order = self.cheb_order
        else:
          self.cheb_order = cheb_order
        a = np.log(a) if self.a_log else a
        b = np.log(b) if self.b_log else b

        if self.a_limits is None:
            self.a_limits = (min(a),max(a))
        if self.b_limits is None:
            self.b_limits = (min(b),max(b))

        a = to_unity_range(a,self.a_limits)
        b = to_unity_range(b,self.b_limits)

        self.ab_cheb = np.polynomial.chebyshev.chebfit(a,b,self.cheb_order-1)
        self.ba_cheb = np.polynomial.chebyshev.chebfit(b,a,self.cheb_order-1)

    
    def a_to_b(self, a):
        if self.cheb_order == 0:
            return a
        a = np.log(a) if self.a_log else a
        a = to_unity_range(a,self.a_limits)
        b = np.polynomial.chebyshev.chebval(a, self.ab_cheb, tensor=False)
        b = from_unity_range(b,self.b_limits)
        return np.exp(b) if self.b_log else b

    def b_to_a(self, b):
        if self.cheb_order == 0:
            return b
        b = np.log(b) if self.b_log else b
        b = to_unity_range(b,self.b_limits)
        a = np.polynomial.chebyshev.chebval(b, self.ba_cheb, tensor=False)
        a = from_unity_range(a,self.a_limits)
        return np.exp(a) if self.a_log else a


class ThermometryChebyshev(ScaledChebyshev):
    def __init__(self,cheb_order=0, name=''):
        ScaledChebyshev.__init__(self,cheb_order)
        self.name = name
        self.a_log = self.b_log = True

    def r_to_t(self, R):
        return self.a_to_b(R)
    def t_to_r(self, T):
        return self.b_to_a(T)
    def set_r_limits(self, d):
        self.set_a_limits(d)
    def set_t_limits(self, d):
        self.set_b_limits(d)
    def get_r_limits(self):
        return self.get_a_limits()
    def get_t_limits(self):
        return self.get_b_limits()

    def load_tr_file(self, filename, columns=['T','R'], header=None):
        self.df = pd.read_csv(filename, comment='#', header=header, names=columns, index_col=False, dtype=np.float32, delim_whitespace=True)
        self.fit(self.df['R'].values,self.df['T'].values)

    def read(self, c):
        if 'name' in c:
          self.name = c['name']
        self.set_r_limits(c['resistance_range'])
        self.set_t_limits(c['temperature_range'])
        self.ab_cheb = c['rt_cheb']
        self.ba_cheb = c['tr_cheb']
        self.cheb_order = len(self.ab_cheb)
    def write(self):
        c = {}
        if self.name != '':
          c['name'] = self.name
        c['resistance_range'] = lof(self.get_r_limits())
        c['temperature_range'] = lof(self.get_t_limits())
        c['rt_cheb'] = lof(self.ab_cheb)
        c['tr_cheb'] = lof(self.ba_cheb)
        return c


#base_folder = '../Thermometer Calibrations/NML SystemE/'
#cx_fn = base_folder + 'NML Cernox X63617 SysE Pb1.cal'
#
#th = ThermometryChebyshev(13)
#th.load_thermometer_file(cx_fn)
#print th.t_to_r(4.2)



class ThermometryInFieldCurve:
    def __init__(self, zero_field_calibration=ThermometryChebyshev(), name=''):
        self.name = name
        self.th = zero_field_calibration # default value performs no transformations
        
    def read(self,c):
        # Data structure: (Actual coefficeints taken from Jeonghoon's work. Values will 
        #                  not be compatible with this implementation)
        # c = {'name' : 'X89065',
        #     'temperature_limits' : [0.3, 4.2],
        #     'field_limits' : [0.0, 18.0],
        #     'coeffs' : [
        #         {'p' : [2.831225792834143,1.734827296534389,-0.0007534116260401212], 'q' : [0.6507947410714512]},
        #         {'p' : [-0.8070888352090849,-0.5868359288037652,0.003293056015132634], 'q' : [1.005139513043437]},
        #         {'p' : [0.162849448895288,0.05254755389668063,0.0002296367139975507], 'q' : [1.019006126715248]},
        #         {'p' : [-0.03026666605430712,-14640880872.92773,-171445522.4987663], 'q' : [1269109740654.296]} ]
        #     }
        if 'name' in c:
          self.name = c['name']
        self.read_limits = c['read_limits']  # (-1,1)
        self.actual_limits = c['actual_limits']  # (-1,1)
        if 'field_limits' in c:
          self.field_limits = c['field_limits'] # not enforced, but important to document
        self.th = ThermometryChebyshev()
        self.th.read(c['thermometer'])
    
        # we wish to make the number of coefficients be self-describing, because 
        # we may need an extra fit component on some coefficients but not others.
        # if we keep p and q separate, we can have a variable number of coefficients
        # and keep the total number of computations and amount of redundancy lower.
        self.fwd_p_coeffs = [np.array(list(reversed(a['p']))) for a in c['fwd_coeffs']]
        self.fwd_q_coeffs = [np.array(list(reversed(a['q']))) for a in c['fwd_coeffs']]
        self.rev_p_coeffs = [np.array(list(reversed(a['p']))) for a in c['rev_coeffs']]
        self.rev_q_coeffs = [np.array(list(reversed(a['q']))) for a in c['rev_coeffs']]

    def write(self):
        c = {}
        if self.name != '':
          c['name'] = self.name
        c['read_limits'] = self.read_limits
        c['actual_limits'] = self.actual_limits
        try:
          c['field_limits'] = self.field_limits
        except AttributeError:
          pass
        c['thermometer'] = self.th.write()
        c['fwd_coeffs'] = [{'p' : lof(reversed(p)), 'q' : lof(reversed(q))} for p,q in zip(self.fwd_p_coeffs,self.fwd_q_coeffs)]
        c['rev_coeffs'] = [{'p' : lof(reversed(p)), 'q' : lof(reversed(q))} for p,q in zip(self.rev_p_coeffs,self.rev_q_coeffs)]
        return c

        
        
    # Procedure to convert (B-field, Resistance or Temperature reading) pair to Actual Temperature or equivalent zero-field resistance
    #  For each chebyshev component:
    #    Use rational fraction coefficients to compute Chebyshev coefficient at field B
    #  Do Chebyshev evaluation at temperature T
    #
    # Pass magnetic field in B and EITHER R for measured resistance or T for temperature reading.
    # output is either 'T' to return actual temperature or 'R' to return equivalent zero-field resistance
    def measured_to_actual(self, B, R=None,T=None, output='T'):
        B = np.abs(B)
        if R is None:
          R = self.th.t_to_r(T)
        
        cheb_order = len(self.fwd_p_coeffs)
        field_ch_calculated = ScaledChebyshev(cheb_order=cheb_order)
        field_ch_calculated.log_a = field_ch_calculated.log_b = True

        field_ch_calculated.set_a_limits(self.read_limits)
        field_ch_calculated.set_b_limits(self.actual_limits)

        field_ch_calculated.ab_cheb = [rational(B,p,q) for p,q in zip(self.fwd_p_coeffs,self.fwd_q_coeffs)]

        b = field_ch_calculated.a_to_b(R)
        return b if output == 'R' else self.th.r_to_t(b)

    # Convert actual temperature or equivalent resistance to values that we expect to read given a certain B-field
    def actual_to_measured(self, B, R=None,T=None, output='T'):
        B = np.abs(B)
        if R is None:
          R = self.th.t_to_r(T)
        
        cheb_order = len(self.rev_p_coeffs)
        field_ch_calculated = ScaledChebyshev(cheb_order=cheb_order)
        field_ch_calculated.log_a = field_ch_calculated.log_b = True

        field_ch_calculated.set_a_limits(self.read_limits)
        field_ch_calculated.set_b_limits(self.actual_limits)

        field_ch_calculated.ba_cheb = [rational(B,p,q) for p,q in zip(self.rev_p_coeffs,self.rev_q_coeffs)]

        a = field_ch_calculated.b_to_a(R)
        return a if output == 'R' else self.th.r_to_t(a)

    # Create a resistance-field calibration surface based on:
    #   R_read     the resistance of the sensor
    #   B          magnetic field
    #   R_actual   resistance this sensor would show at zero magnetic field and the same temperature
    #
    #   p_order, q_order    numerator and denominator order for fitting in field
    #   cheb_order          order of the temperature-correcting chebysheb polynomial (does not need to be same
    #                       order as the zero-field calibration)
    def fit_field_ranges_rr(self, B, R_read, R_actual, p_order, q_order, cheb_order, field_partition_points=None):
        B = np.abs(B)
        min_field = float(min(B))
        max_field = float(max(B))
        if field_partition_points is None:
            field_partition_points = np.linspace(min_field,max_field,(p_order + q_order)*3)
        field_lower_points = np.asarray(field_partition_points[:-1])
        field_upper_points = np.asarray(field_partition_points[1:])
        field_ranges = [(l,u) for l,u in zip(field_lower_points,field_upper_points)]
        self.field_limits = (min(field_lower_points),max(field_upper_points))
                                
        self.read_limits = (min(R_read),max(R_read))
        self.actual_limits = (min(R_actual),max(R_actual))

        geomspace_r = np.geomspace(self.read_limits[0],self.read_limits[1])

        fwd_cheb_coeffs = list()
        rev_cheb_coeffs = list()
        print("field_ranges=",field_ranges)
        field_centers = [sum(x)/2 for x in field_ranges]
        for low,high in field_ranges:
            # filter only the points within our field range
            # rba stands for read-bfield,actual
            rba = [rba for rba in zip(R_read,B,R_actual) if rba[1]>=low and rba[1]<high]
            
            # here, divide by the zero-field calibration curve so we only compute
            # a difference curve
            
            field_ch = ScaledChebyshev(cheb_order=cheb_order)
            field_ch.log_a = field_ch.log_b = True
            field_ch.set_a_limits(self.read_limits)
            field_ch.set_b_limits(self.actual_limits)

            r,b,a = list(zip(*rba))

            field_ch.fit(r, a)
            print("field=",(low+high)/2,"ab_cheb: ",field_ch.ab_cheb)
            print("field=",(low+high)/2,"ba_cheb: ",field_ch.ba_cheb)
            
            #x-axis is log(resistance read), y-axis is log(resistance corresponding to actual)
            plt.scatter(np.log(r),np.log(a),label=str((low+high)/2))
            plt.plot(np.log(geomspace_r), np.log(geomspace_r))
            plt.plot(np.log(geomspace_r),np.log([field_ch.a_to_b(r) for r in geomspace_r]))
            plt.show()
            
            fwd_cheb_coeffs.append(field_ch.ab_cheb)
            rev_cheb_coeffs.append(field_ch.ba_cheb)
        # transpose list
        fwd_cheb_coeffs = list(map(list, list(zip(*fwd_cheb_coeffs))))
        rev_cheb_coeffs = list(map(list, list(zip(*rev_cheb_coeffs))))
        
        # find the cheb coeffs for a 1:1 line
        ch_line = ScaledChebyshev(cheb_order=cheb_order)
        ch_line.log_a = ch_line.log_b = True
        ch_line.set_a_limits(self.read_limits)
        ch_line.set_b_limits(self.actual_limits)
        ch_line.fit(geomspace_r, geomspace_r)
        print("straight line cheb:",ch_line.ab_cheb, "and in reverse: ",ch_line.ba_cheb)
            
        self.fwd_p_coeffs = list()
        self.fwd_q_coeffs = list()
        for c_vs_field,zero_point in zip(fwd_cheb_coeffs,ch_line.ab_cheb):
            p,q = fit_rational(field_centers,c_vs_field,{'len_p':p_order,'len_q':q_order,'p0':zero_point}) #
            self.fwd_p_coeffs.append(p)
            self.fwd_q_coeffs.append(q)
            print("PQ(0) must =",zero_point)
            plt.scatter(field_centers,c_vs_field)
            plt.plot(np.linspace(0,self.field_limits[1]),[rational(x,p,q) for x in np.linspace(0,self.field_limits[1])])
            plt.show()
            print(p,q)
        print("fwd p coeffs=",self.fwd_p_coeffs)
        print("fwd q coeffs=",self.fwd_q_coeffs)

        self.rev_p_coeffs = list()
        self.rev_q_coeffs = list()
        for c_vs_field,zero_point in zip(rev_cheb_coeffs,ch_line.ba_cheb):
            p,q = fit_rational(field_centers,c_vs_field,{'len_p':p_order,'len_q':q_order,'p0':zero_point}) #
            self.rev_p_coeffs.append(p)
            self.rev_q_coeffs.append(q)
            print("PQ(0) must =",zero_point)
            plt.scatter(field_centers,c_vs_field)
            plt.plot(np.linspace(0,self.field_limits[1]),[rational(x,p,q) for x in np.linspace(0,self.field_limits[1])])
            plt.show()
            print(p,q)
        print("rev p coeffs=",self.rev_p_coeffs)
        print("rev q coeffs=",self.rev_q_coeffs)


    # convenience functions that let you call the calibration with various resistance or temperature datasets
    def fit_field_ranges_rt(self, B, R_read, T_actual, p_order, q_order, cheb_order, field_partition_points=None):
        self.fit_field_ranges_rr(B, R_read, self.th.t_to_r(T_actual), p_order, q_order, cheb_order, field_partition_points)
    def fit_field_ranges_tr(self, B, T_read, R_actual, p_order, q_order, cheb_order, field_partition_points=None):
        self.fit_field_ranges_rr(B, self.th.t_to_r(T_read), R_actual, p_order, q_order, cheb_order, field_partition_points)
    def fit_field_ranges_tt(self, B, T_read, T_actual, p_order, q_order, cheb_order, field_partition_points=None):
        self.fit_field_ranges_rr(B, self.th.t_to_r(T_read), self.th.t_to_r(T_actual), p_order, q_order, cheb_order, field_partition_points)

            
    # outdated, but let's not lose the algorithm yet in case we can update it to refine our high dimensional fits.
    if False:
      # Too many parameters to fit. Perhaps it could be used to refine parameters. Otherwise, avoid.
      # Given: equal length arrays of (R,B,T) points. Find coefficients.
      def fit_arbitrary_points(self, R, B, T, p_order, q_order, cheb_order):
          rough_curve = ThermometryChebyshev(cheb_order)
          self.cheb_order = cheb_order
          
          # Get approximate chebyshev fit over all field values. This also sets appropriate
          # limits so the fits stay close to the [-1:1] domain
          rough_curve.cheb_thermometer_calibration(R, T)
          print("rough calibration for all B:",rough_curve.fwd_cheb)
          
          self.p_coeffs = np.zeros((self.cheb_order,p_order))
          #print "setting up p array:",self.p_coeffs
          self.q_coeffs = np.zeros((self.cheb_order,q_order))
          #print "setting up q array:",self.q_coeffs
          self.p_coeffs[:,0] = rough_curve.fwd_cheb # put zeroth order guess in coefficients
          #if q_order>0:
          #    self.q_coeffs = np.ones(len(self.q_coeffs))*3.0 
          #else:
          #    self.q_coeffs = np.asarray([])
          
          R = np.log(R)
          T = np.log(T)
          
          # flattened array of variables for curve_fit
          v1 = self.p_coeffs.reshape(self.cheb_order*p_order)
          if q_order>0:
              v2 = self.q_coeffs.reshape(self.cheb_order*q_order)
              v = np.concatenate( (v1, v2) )
          else:
              v = v1
          s = self
                  
          # define callback function with some constants baked in
          def rational_vs_cheb_eval(RB, *v):
              #print "v=",v
              R,B = RB
              P = np.asarray(v[:self.cheb_order*p_order]).reshape((self.cheb_order,p_order))
              Q = np.asarray(v[self.cheb_order*p_order:]).reshape((self.cheb_order,q_order))
              #print "rational_vs_cheb_eval(): (P,Q) =",(P,Q)
              cheb = np.zeros(self.cheb_order)
              return [s.actual_temperature(r,b,P,Q) for r,b in zip(R,B)]

          print("initial guess v=",v)
          popt, pcov = scipy.optimize.curve_fit(rational_vs_cheb_eval,(R,B),T,v,epsfcn=1e-9)
          self.p_coeffs = popt[:self.cheb_order*p_order].reshape((self.cheb_order,p_order))
          self.q_coeffs = popt[self.cheb_order*p_order:].reshape((self.cheb_order,q_order))


# little snippet of code to express 2d chebyshev polynomials using the recursion relation
# commented so as to not make the library require sympy
if False:
    import sympy
    x,y = sympy.symbols('x y')

    # given a list of n chebychev polynomials, generate the n+1 entry
    def cheb_n_plus_one(x, cheb_list):
        n = len(cheb_list)-1
        cheb_list.append(sympy.expand(2*x*cheb_list[n]-cheb_list[n-1]) )
        return cheb_list

    chebx = [1.0+0*x, x]
    for i in range(2):
        chebx = cheb_n_plus_one(x,chebx)
    cheby = [c.subs(x,y) for c in chebx]
    chebxy = {}
    for i,cx in enumerate(chebx):
        for j,cy in enumerate(cheby):
            chebxy[(i,j)] = sympy.simplify(cx*cy)
    display(chebxy)
