import re
import dateutil
import datetime
import pandas as pd
import numpy as np
import scipy
import scipy.interpolate
import string

from .useful import *
from . import interpolation

debug = False



# respline data where duplicate values indicate lack of update, not constant response
def un_stairstep_data(X, Y, tolerance=None):
    if tolerance is None:
      tolerance = 0.0

    # first, ensure sorting in X
    XI = np.argsort(X)
    X = shuffle_list(X, XI)
    Y = shuffle_list(Y, XI)

    # make a list if X is an array
    X = list(X)
    
    current_value = Y[0]
    delete_list = []
    for i in range(1,len(X)-1): # preserve endpoints
        if abs(current_value-Y[i]) < tolerance:
            delete_list += [i]
        else:
            current_value = Y[i]
    for i in reversed(delete_list):
        del(X[i])
        del(Y[i])
    
    return list(X),list(Y)

# delete points that would result in a non-monotonic function
def make_monotonic(X,Y):
    # first, ensure sorting in X
    if Y[0] < Y[-1]:
        reverse_direction = False # Y rises with rising X
        XI = np.argsort(X)
        #print 'forward direction'
    else:
        # Y falls with rising X
        reverse_direction = True
        XI = np.argsort(X[::-1])
        #print 'reverse direction'
    
    X = shuffle_list(X, XI)
    Y = shuffle_list(Y, XI)  

    current_value = Y[0]
    delete_list = []
    for i in range(1,len(X)-1): # preserve endpoints
        if current_value > Y[i]:
            delete_list += [i]
            #print "scheduling ",i,"for deletion"
        else:
            current_value = Y[i]
    for i in reversed(delete_list):
        del(X[i])
        del(Y[i])
    
    if reverse_direction:
        return list(reversed(X)),list(reversed(Y))
    else:
        return X,Y


exp_notation_regex = r'([+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:E[+\-]?\d+)?|NAN)' # find regex

# f is a file handle, not a filename
# require_line_count is how many lines to read into the data to confirm the numerical structure of the file
# returns header list of lines and the pandas dataframe df for all the data in the file
def read_nml_header_and_data(f,require_line_count=20, dtype=np.float32):
    # Regex strategy: Search for columns of float/scientific numbers as a whole-line match.
    # If found, figure out how many columns we have. Insist that number be constant.
    r1 = r'^\s*(?:([+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:E[+\-]?\d+)?|NAN)\s*)+$' # match regex
    count = 0
    columns = 0
    x0 = x1 = f.tell() # Mark our position in the file for seekback

    i = 0 # Line number
    y = 0 # Mark our position in the header buffer list
    header = []

    for l in iter(f.readline, ''):
        header.append(l)
        if count == 0:
            if not re.match(r1,l,flags=re.I) is None:
                newcolumns = len(re.findall(exp_notation_regex,l,flags=re.I))
            else:
                newcolumns = 0
                x0 = x1
                x1 = f.tell()
                y = i
        if newcolumns == columns and newcolumns > 0:
            count += 1
        else:
            columns = newcolumns
            count = 0
        if count == require_line_count:
            break
        i += 1
    f.seek(x0) # Seek back to the start of the column names
    header = header[:y]
    if debug:
      print("Header: ",header)
    if len(header) == 0:
      raise ValueError('No data in file')
    df = pd.read_csv(f, dtype=dtype, delim_whitespace=True)
    return header,df

class NHMFL_DC_Data:
    # instance variables:
    #   all variables are written when the file is written to disk at the end, so everything reflects
    #   the state of the data at that time, including timestamp, trigger information, and so on.
    #header = []
    #filename # read from header, not from NHMFL_DC_Data(fn)
    #end_time # datetime object when the file was written (usually close to the last point)
    #magnet   # name of magnet if that was set
    #trigger_spacing # when file was written
    #trigger_column  # when file was written
    #scale = []  # list of column scalers
    #df   # pandas DataFrame with the actual data in it. Access using df['Field'].values
          # df['Datetime'] column has been added so data may be examined vs time of day
    
    def __init__(self,fn):
        # 'rU' mode reads "universal" line endings
        with open(fn,'rU') as f:
        #with open(fn,'r') as f:
            header,self.df = read_nml_header_and_data(f)

            # Some old files have no txt extension
            #r3 = r'^(.*\.(?:txt|text)) ((?:Mon|Tue|Wed|Thu|Fri|Sat|Sun).*(?:am|pm)) (.*)'
            r3 = r'^(.*) ((?:Mon|Tue|Wed|Thu|Fri|Sat|Sun).*(?:am|pm)) (.*)'
            self.filename,datestring,self.user = re.match(r3,header[0],flags=re.I).groups()       
            self.scale = [float(x) for x in re.findall(exp_notation_regex,header[3],flags=re.I)]
            self.magnet = re.findall(r'^Magnet.*: (.*)$',header[1])[0]
            trigger_spacing,trigger_column = re.match(r'Trigger spacing of '+exp_notation_regex+' on column #([0-9]+)',header[2]).groups()
            self.trigger_spacing = float(trigger_spacing)
            self.trigger_column = int(trigger_column)
            self.end_time = dateutil.parser.parse(datestring)

            # drop the file number from column names
            self.df = self.df.rename(columns=lambda x: re.match('^(.*)_\d{3}',x).groups()[0])
            self.df['Datetime'] = pd.Series(pd.to_datetime(self.end_time) + pd.to_timedelta(self.df['Timestamp'].values - self.df['Timestamp'].values[-1], unit='s'))

    # shortcut to the numpy array of a data field
    def __getitem__(self, x):
      return self.df[x].values
    
    def clean_field(self,smoothing_time=5.0,remove_stairsteps=True,ensure_monotonic=False):
      time = self.df['Timestamp'].values
      field = self.df['Field'].values

      # remove NaN
      time = time[~np.isnan(field)]
      field = field[~np.isnan(field)]

      if remove_stairsteps:
        time,field = un_stairstep_data(time,field)
      if ensure_monotonic:
        time,field = make_monotonic(time,field)

      # remove-stairstep and ensure-monotonic preserve endpoints, which might cause artifacts.
      # let's just remove the endpoints ourselves now to avoid that problem.
      cs = interpolation.interpolate_smoothly(time,field,dx=smoothing_time,return_cubic_spline=True)

      # Cleaning can result in loss of points at the start or end of dataset. Truncate
      # the dataframe to remove these points.
      tmin = min(time[1:-2])
      tmax = max(time[1:-2])
      imin,imax = interpolation.bracket_interval(self.df['Timestamp'].values,tmin,tmax)
      self.df = self.df.truncate(before=imin,after=imax)

      self.df['Field'] = scipy.interpolate.splev(self.df['Timestamp'].values,cs)




# current and temperature are arrays.
# r_coeff is a list of polynomial coefficients in order 0, P(1), P(2) in mOhms
# t_coeff is the temperature coefficent of resistance in Celsius
# inductance is in mH
# dt is time between points
def rmps_compute_expected_voltage(current, temperature, r_coeff, t_coeff, inductance, dt):
    dIdt = np.zeros(len(current))
    dIdt[1:] = np.diff(current)/dt

    r_coeff = np.asarray([0.0, r_coeff[0], 0.0, r_coeff[1], 0.0, r_coeff[2]])[::-1]
    return np.polyval(r_coeff,current) * (1.0+t_coeff*temperature) + inductance*dIdt


class NHMFL_RMPS_Data:
    # instance variables:
    #   all variables are written when the file is written to disk at the end, so everything reflects
    #   the state of the data at that time, including timestamp, trigger information, and so on.
    #header = []
    #filename # read from NHMFL_RMPS_Data(fn)
    #hostname # RMPSxMac
    #cell # cell number
    #dt # delta time per point
    #end_time # datetime object when the file was written (usually close to the last point)
    #coil_inductances[M] # list of inductances in mH
    #N_coils # number of coils in this magnet
    #max_coil_currents[M] # list of max currents in kA
    #max_voltage_deviations[M] # list of trip levels for each coil
    #max_coil_temperatures[M] # list of maximum allowed coil temperatures
    #pressure_drop_range # (min,max) pair of expected water pressure drop across magnet
    #coil_coeffs[M][N] # where M is the number of coils and N is a polynomial in mOhms (R0 = 0)
    #coil_temp_coeffs[M] # temperature coeffients for each coil in change per degree C
    #cell_calibration_scales[P] # values used to calibrate cell to RMPS
    #cell_calibration_offsets[P] # values used to calibrate cell to RMPS

    #df   # pandas DataFrame with the actual data in it. Access using df['Column'].values
          # df['Datetime'] column has been added so data may be examined vs time of day
    
    def __init__(self,fn):
        # 'rU' mode reads "universal" line endings
        with open(fn,'rU') as f:
            header,self.df = read_nml_header_and_data(f)
            self.filename = fn

            # if there is no hostname, then we should accept anything before the comma as a hostname
            #r3 = r'^Normal Data Dump,?\s+([A-Za-z0-9\-]*[A-Za-z0-9]),?\s+for\s+Cell\s+([0-9]+)\s+((?:Mon|Tue|Wed|Thu|Fri|Sat|Sun).*(?:am|pm))\s*'
            r3 = r'^Normal Data Dump,?\s+([^,]*),?\s+for\s+Cell\s+([0-9]+)\s+((?:Mon|Tue|Wed|Thu|Fri|Sat|Sun).*(?:am|pm))\s*'
            self.hostname,self.cell,datestring = re.match(r3,header[0],flags=re.I).groups()       
            self.dt = float(re.findall(exp_notation_regex,header[1],flags=re.I)[0])
            self.end_time = dateutil.parser.parse(datestring)
            self.coil_inductances = [float(x) for x in re.findall(exp_notation_regex,header[3],flags=re.I)]
            self.N_coils = len(self.coil_inductances)
            self.coil_names = list(string.ascii_uppercase)[:self.N_coils]
            self.max_coil_currents = [float(x) for x in re.findall(exp_notation_regex,header[4],flags=re.I)]
            self.max_coil_temperatures = [float(x) for x in re.findall(exp_notation_regex,header[5],flags=re.I)]
            self.pressure_drop_range = [float(x) for x in re.findall(exp_notation_regex,header[6],flags=re.I)][0:2]
            self.coil_coeffs = [False]*self.N_coils # initalize list so we may access it by indexing
            for i,c in enumerate(self.coil_coeffs):
                self.coil_coeffs[i] = [float(x) for x in re.findall(exp_notation_regex,header[9+i],flags=re.I)]
            self.cell_temp_coeffs = [float(x) for x in re.findall(exp_notation_regex,header[9+self.N_coils],flags=re.I)]
            self.cell_calibration_scales = [float(x) for x in re.findall(exp_notation_regex,header[10+self.N_coils],flags=re.I)]
            self.cell_calibration_offsets = [float(x) for x in re.findall(exp_notation_regex,header[11+self.N_coils],flags=re.I)]

            # current comes without sign. use the sign of the A-coil voltage to determine sign of current
            self.df['Current'] = self.df['Current'].values * np.sign(self.df['VA'].values)

            self.df['Datetime'] = pd.Series(pd.to_datetime(self.end_time) + pd.to_timedelta(self.dt*(self.df.index.values - self.df.index.values[-1]), unit='s'))
            self.df.set_index('Datetime', inplace=True)

    def deviations(self):
      for i,name in enumerate(self.coil_names):
          self.df['V'+name+'_expect'] = rmps_compute_expected_voltage(self.df['Current'].values,self.df['Tin'].values,self.coil_coeffs[i],self.cell_temp_coeffs[i],self.coil_inductances[i], self.dt)
          self.df['V'+name+'_deviation'] = self.df['V'+name].values - self.df['V'+name+'_expect'].values


if __name__ == '__main__':
  import matplotlib.pyplot as plt

  dc = NHMFL_DC_Data('bets2018.158.txt')
  #plt.plot(dc.df['Field'],dc.df['PCC'])
  #plt.show()

  print("dc field:",dc.df['Field'].values[-7:])

  t,b = make_monotonic(dc.df['Timestamp'].values,dc.df['Field'].values)
  t,b = un_stairstep_data(t,b)
  
  cs = interpolation.interpolate_smoothly(t,b,dx=1.5,return_cubic_spline=True)
  b2 = scipy.interpolate.splev(dc.df['Timestamp'].values,cs,der=0)

  plt.plot(dc.df['Timestamp'].values,dc.df['Field'].values,label='original')
  plt.plot(dc.df['Timestamp'].values,b2,label='smoothed')
  plt.legend()
  plt.show()

  orig_cs = scipy.interpolate.splrep(dc.df['Timestamp'].values,dc.df['Field'].values)

  b2diff = scipy.interpolate.splev(dc.df['Timestamp'].values,cs,der=1)
  plt.plot(dc.df['Timestamp'].values,scipy.interpolate.splev(dc.df['Timestamp'].values,orig_cs,der=1),label='original')
  plt.plot(dc.df['Timestamp'].values,b2diff,label='smooth')
  plt.legend()
  plt.show()
