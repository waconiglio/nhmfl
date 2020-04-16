import dateutil.parser as dt
import copy
import math
from enum import Enum
import collections


# perhaps this would be better described as the sample's current state
class Sample:
    # name is the name of the sample
    def __init__(self, name='', **kwargs):
        self.name = name
        self.variables = kwargs
        try:
          del(self.variables['name'])
        except KeyError:
          pass
    def __getitem__(self,v):
        if v == 'name':
            return self.name
        return self.variables[v]

    def __setitem__(self,v,value):
        if v == 'name':
            raise KeyError('\'name\' is reserved in class Sample')
        self.variables[v] = value

    def __delitem__(self,v):
        if v == 'name':
            raise KeyError('\'name\' is reserved in class Sample')
        del(self.variables[v])

    def __contains__(self, v):
        return v in self.variables

    def __eq__(self, other):
        try:
            return (self.name == other.name and self.variables == other.variables)
        except:
            return hash(self) == hash(other)

    def __str__(self):
        v = self.variables.copy()
        v['name'] = self.name
        return str(v)
    
    def __repr__(self):
        return str(self)

    # if we have a mapper function, return it
    def __getattr__(self, name):
        # avoid recursion by checking for this one explicitly
        if name == 'variables':
            return []

        if name not in self.variables:
            raise AttributeError('{:} is not a method of Sample, nor is it in the variable list.'.format(name))
        if not callable(self.variables[name]):
            raise AttributeError('{:} is not callable, as required to be returned from Sample.'.format(name))
        def run_map():
            self.merge(self.variables[name](self))
            return self
        return run_map

    def __hash__(self):
        return hash(self.name + repr(self.variables))

    def __len__(self):
        return len(self.variables)

    # just like dict.update(), Sample.update can take a dictionary-like via kargs or a list of keys
    # and values via kwargs
    def update(self, *kargs, **kwargs):
        for k in kargs:
            self.variables.update(k)
        self.variables.update(kwargs)
    def up(self, **kwargs):
        return self.update(**kwargs)
    def remove(self, *var_list):
        for v in var_list:
            try:
                del(self.variables[v])
            except KeyError:
                pass

    
    def load(self,d):
        if type(d) == type(self):   # load from another Sample
            self.name = d['name']
            self.variables = {}
            self.merge(d)
        else:   # load from a dictionary
            self.name = d['name']
            self.variables = {}
            self.variables.update(d)
            del(self.variables['name'])
        return self
    def save(self):
        d = self.variables.copy()
        d.update({'name':self.name})
        return d
    def merge(self, other): # merge variables (except name) from another sample into ours. other overwrites self.
        for v in other.variables:
          self.variables[v] = other.variables[v]
        return self
    def copy(self):
        c = Sample(self.name)
        c.merge(self)
        return c
        
class ExperimentIterator(collections.Iterator):
    def __init__(self, e, rev=False):
        if rev:
            self.it = reversed(e.snapshots)
        else:
            self.it = iter(e.snapshots)

    def __next__(self):
        return Sample().load(next(self.it))

class Experiment:
    class Mode(Enum):
        snapping = 1
        analyzing = 2


    def __init__(self, snapshots=[]):
        self.default_exp = Sample()
        self.load(snapshots)
        self.primary_key('{:d}', ['sequence'])
        if len(snapshots) == 0:
            self.mode = Experiment.Mode.snapping
        else:
            self.mode = Experiment.Mode.analyzing

    def primary_key(self, prikey_format='[{:}_{:03d}]', variables=['name','n'], functions=None):
        self.prikey_format = prikey_format
        if functions is not None:
            self.prikey_kargs = functions
        else:
            self.prikey_kargs = [lambda x: self.snapshots[x] for x in variables]


    # When documenting the experiment, we often have need to index based on the sample name 
    # to set unique variables or such. When analyzing data, we need access to the snapshot rows.
    # 
    # Overload __getitem__ to work in either mode, depending on what we are doing at the time.
    # Set the mode enum, or let it get set to Mode.snapping automatically at creation, then change
    # to Mode.analyzing when we call end()
    def __getitem__(self,key):
        if self.mode == Experiment.Mode.snapping:
            return [s for s in self.samples if s['name'] == key][0]
        else:
            return [s for s in self.snapshots if s['primary_key'] == key][0]

    def __len__(self):
        return len(self.snapshots)

    def __iter__(self):
        return ExperimentIterator(self)

    def __reversed__(self):
        return ExperimentIterator(self, rev=True)


    # operates on primary key
    def __contains__(self, a):
        return any([x['primary_key'] == a for x in self.snapshots])

    # i can't think of a time to use this, but here it is. maybe if we import some lists that
    # are missing private keys, we need to run it to get indexing capability.
    def gen_prikey(self, overwrite=False):
        for s in self.snapshots:
            if 'primary_key' not in s or overwrite:
                self.snapshot_gen_prikey(s)

    # generate private key for snapshot s
    # format follows string.format().
    # further arguments refer to {:} portions of the format specifier, where the functions take the 
    #   sample data as input.
    def snapshot_gen_prikey(self, s):
        output_key = 'primary_key'
        format_args = [k(s) for k in self.prikey_kargs]
        prikey = self.prikey_format.format(*format_args)
        extension = ord('a')

        # at some point, we should lift the 26-limit here and create aa, ab, ac... as needed
        try:
            while any(x[output_key] == prikey if output_key in x else False for x in self.snapshots):
                prikey = self.prikey_format.format(*format_args) + chr(extension)
                extension += 1
                if extension > ord('z'):
                    raise KeyError
        except StopIteration:
            pass
        s[output_key] = prikey

    
    # perhaps these should be named include and exclude, so as to differentiate between add and remove for variables.
    def include(self,*samples):
        # Add to sample list. Samples are indexed by name.
        for e in samples:
            e = copy.copy(e) # shallow copy, so that the original may be reused later
            self.samples[e.name] = e.merge(self.default_exp)
            
        
    def inc(self,*samples):
        return self.include(*samples)

    def exclude(self,*sample_names):
        for e in sample_names:
            del(self.samples[e])
    def exc(self,*samples):
        return self.exclude(*samples)

    # clear all samples and the internal state off the experiment and start fresh
    # keep all the existing snapshots
    def clear(self):
        self.samples = {}
        self.default_exp = Sample()
        
    # update child samples with new or changed variables
    def update(self, **kwargs):
        # update child samples
        for e in self.samples:
            self.samples[e].update(**kwargs)
        self.default_exp.update(**kwargs)
    def up(self, **kwargs):
        return self.update(**kwargs)

    # we update the date all the time, so let's simplify that line
    # call as p.date('29 april 2019')
    def date(self, d):
        self.update(date=dt.parse(d))
            
    # update and snapshot in one step
    def snap(self, **kwargs):
        self.update(**kwargs)
        for e in self.samples:
            exp_copy = self.samples[e].copy()
            exp_copy.update(sequence=self.sequence)
            self.snapshot_gen_prikey(exp_copy)
            self.snapshots.append(exp_copy)
        self.sequence += 1

    # adjust something on the most recent snapshot temporarily withou updating the ongoing state of the samples
    def amend(self, **kwargs):
        for sample,snapshot in zip(reversed(self.samples),reversed(self.snapshots)):
            snapshot.update(**kwargs)

    # call at end of experiment logging when you are ready to analyze data
    def end(self):
        self.mode = Experiment.Mode.analyzing

        
    # remove a variable from child samples
    def remove(self, *var_list):
        for e in self.samples:
            self.samples[e].remove(*var_list)
        self.default_exp.remove(*var_list)
    # alias for remove
    def rm(self, *var_list):
        return self.remove(*var_list)
    
    # save the entire history
    def save(self):
        return self.snapshots
    
    # load the entire history; if blank, initialize the object
    def load(self, snapshots):
        self.samples = {}
        self.snapshots = list(snapshots)
        self.sequence = 0 # increments with each snapshot

        if len(self.snapshots) == 0:
            return

        # find all unique samples as of final snapshot
        self.sequence = self.snapshots[-1]['sequence'] + 1
        
        # create samples as of the final sequence snapshot
        samples = [s for s in self.snapshots if s['sequence'] == self.sequence - 1]
        for e in samples:
            e = e.copy()
            try:
                del(e['sequence'])
            except:
                pass
            ee = Sample()
            ee.load(e)
            self.include(ee)

        return self
    
    # filter data with a function (perhaps a lambda), but protect against KeyError exceptions
    # Returns a generator object. Make a list with list() or iterate with a comprehension.
    # example:
    # p.filter(lambda e : e['Temperature'] < 1.5)
    #  output: a list we can load with Experiment.load() that contains all the datasets < 1.5 kelvin
    def filter(self, condition):
        l = []
        for e in self.save():
            try:
                if condition(e):
                    l.append(e)
            except:
                pass
        p = copy.copy(self)
        return p.load(l)

    # simple version of filter that matches an exact set of values
    def match(self, **kwargs):
        l = []
        for e in self.save():
          match = True
          try:
            for k in kwargs:
              # check for a date string and convert it if necessary
              if k == 'date' and isinstance(kwargs[k], type('')):
                if e[k] != dt.parse(kwargs[k]):
                  match = False
                  break
              # not a date in string format; carry on.
              elif e[k] != kwargs[k]:
                match = False
                break
          except:
            match = False

          if match:
            l.append(e)

        p = copy.copy(self)
        return p.load(l)

    def map(self, function):
        e = Experiment().load( map(function,self) )
        return e

# class Angle transforms angle data (and makes sure you don't do it twice), generally for use on a Experiment
# object that keeps track of a bunch of Sample objects. Assume default key names to keep the calling
# code short. If you need different ones, see AngleTransformer below.
#
# p = Experiment
#  ... populate experiment ...
# p.update(angle_mapper=Angle.shifted(-10))
# p.snap(angle=45)
# p = Experiment.map(Angle.transform)
class Angle(float):
    def __new__(cls, val, mapper=lambda x : x):
        if val is None:             # converts Angle(None) to legal float representation NaN
            return super(Angle, cls).__new__(cls, 'NaN')
        if isinstance(val,Angle):   # covers case where mapper operation has already been applied
            return val
        return super(Angle, cls).__new__(cls, mapper(val))
    def __str__(self):
        if math.isnan(self):
            return '(no angle)'
        return '{:f} deg'.format(float(self))
    def __repr__(self):
        return self.__str__()
    
    identity = staticmethod(lambda x : x)
    none = staticmethod(lambda x : Angle('NaN'))
    @staticmethod
    def fixed(x):
        return lambda y : x
    @staticmethod
    def shifted(d):
        return lambda x : x + d

    # Apply angle transformations to an sample dictionary. Use with Experiment.map()
    @staticmethod
    def transform(dict_like, angle_variable='angle', angle_mapper='angle_mapper', mask_key_error=True):
      try:
        if angle_variable in dict_like and angle_mapper in dict_like:
            dict_like[angle_variable] = Angle(dict_like[angle_variable], dict_like[angle_mapper])
        if dict_like[angle_variable] is None:
            del dict_like[angle_variable]
        return dict_like
      except KeyError:
        if mask_key_error:
          return dict_like
        else:
          raise

# if you need to prepare a function for a call to map, but you are using different names for your keys,
# p = Experiment
#  ... populate experiment ...
# t = AngleTransformer(angle_variable='angle_1')
# p = Experiment.map(t)
class AngleTransformer:
    def __init__(self,angle_variable='angle', angle_mapper='angle_mapper', mask_key_error=True):
        self.angle_variable = angle_variable
        self.angle_mapper = angle_mapper
        self.mask_key_error = mask_key_error
    def __call__(self,dict_like):
        return Angle.transform(dict_like, self.angle_variable, self.angle_mapper, self.mask_key_error)

# call a function saved in a dictionary with the contents of the entire dictionary as a parameter. save the result
# back to that dictionary.
#
# CODE:
# A = [ {'a' : 1, 'b' : 2, 'op' : lambda D : D['a']+D['b'] } ,
#       {'a' : 5, 'b' : 1, 'op' : lambda D : D['a']-D['b'] } ]
# map(ExperimentVariableMapper('op'), A)
#
# OUTPUT:
# [3, 4]
class ExperimentVariableMapper:
  def __init__(self, key, mask_key_error=True):
    self.key = key
    self.mask_key_error = mask_key_error
  def __call__(self, dict_like):
    try:
      return dict_like[self.key](dict_like)
    except KeyError:
      if self.mask_key_error:
        return dict_like
      else:
        raise

#class DatasetCacheFile:
    # opens and closes hdf5 file as needed
    # creates groups necessary to return the requested node
        
#class CachedDatasetGroup:
    # caches a bunch of related datasets in an hdf5 group with a hash
    #
    # methods: get() tries to load the data from the cache
    #          set() stores the hash and data, replacing everything in the group
