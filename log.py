import dateutil.parser as dt
import copy
import math


# perhaps this would be better described as the sample's current state
class Sample:
    # name is the name of the sample
    def __init__(self, name='', **kwargs):
        self.name = name
        self.variables = kwargs
    def __getitem__(self,v):
        return self.variables[v]
    def update(self, **kwargs):
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
        self.name = d['name']
        self.variables = {}
        self.variables.update(d)
        del(self.variables['name'])
    def save(self):
        d = self.variables.copy()
        d.update({'name':self.name})
        return d
    def merge(self, other): # merge variables (except name) from another sample into ours. other overwrites self.
        for v in other.variables:
          self.variables[v] = other.variables[v]
        return self
            
        
class Experiment:
    def __init__(self, snapshots=[]):
        self.default_exp = Sample()
        self.load(snapshots)
        
    # p['var_name']
    def __getitem__(self,sample_name):
        return self.samples[sample_name]
    
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
            exp_copy = self.samples[e].save()
            exp_copy.update({'sequence':self.sequence})
            self.snapshots.append(exp_copy)
        self.sequence += 1
        
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
        samples = filter(lambda s: s['sequence'] == self.sequence - 1, self.snapshots)
        for e in samples:
            e = e.copy()
            del(e['sequence'])
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
              if k == 'date' and type(kwargs[k]) == type(''):
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
        p = copy.copy(self)
        return p.load(map(function,self.save()))


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
    def __new__(cls, val, mapper=Angle.identity):
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
    def transform(dict_like, angle_variable='angle', angle_mapper='angle_mapper'):
        if angle_variable in dict_like and angle_mapper in dict_like:
            dict_like[angle_variable] = Angle(dict_like[angle_variable], dict_like[angle_mapper])
        if dict_like[angle_variable] is None:
            del dict_like[angle_variable]
        return dict_like

# if you need to prepare a function for a call to map, but you are using different names for your keys,
# p = Experiment
#  ... populate experiment ...
# t = AngleTransformer(angle_variable='angle_1')
# p = Experiment.map(t)
class AngleTransformer:
    def __init__(self,angle_variable='angle', angle_mapper='angle_mapper'):
        self.angle_variable = angle_variable
        self.angle_mapper = angle_mapper
    def __call__(self,dict_like):
        return Angle.transform(dict_like, self.angle_variable, self.angle_mapper)

