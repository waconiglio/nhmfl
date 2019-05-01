import dateutil.parser as dt
import copy


# perhaps this would be better described as the experiment's current state
class Experiment:
    # name is the name of the experiment
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
    def merge(self, other): # merge variables (except name) from another experiment into ours. other overwrites self.
        for v in other.variables:
          self.variables[v] = other.variables[v]
        return self
            
        
class Probe:
    def __init__(self, snapshots=[]):
        self.default_exp = Experiment()
        self.load(snapshots)
        
    # p['var_name']
    def __getitem__(self,experiment_name):
        return self.experiments[experiment_name]
    
    # perhaps these should be named include and exclude, so as to differentiate between add and remove for variables.
    def include(self,*experiments):
        # Add to experiment list. Experiments are indexed by name.
        for e in experiments:
            e = copy.copy(e) # shallow copy, so that the original may be reused later
            self.experiments[e.name] = e.merge(self.default_exp)
            
        
    def inc(self,*experiments):
        return self.include(*experiments)

    def exclude(self,*experiment_names):
        for e in experiment_names:
            del(self.experiments[e])
    def exc(self,*experiments):
        return self.exclude(*experiments)

    # clear all experiments and the internal state off the probe and start fresh
    # keep all the existing snapshots
    def clear(self):
        self.experiments = {}
        self.default_exp = Experiment()
        
    # update child experiments with new or changed variables
    def update(self, **kwargs):
        # update child experiments
        for e in self.experiments:
            self.experiments[e].update(**kwargs)
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
        for e in self.experiments:
            exp_copy = self.experiments[e].save()
            exp_copy.update({'sequence':self.sequence})
            self.snapshots.append(exp_copy)
        self.sequence += 1
        
    # remove a variable from child experiments
    def remove(self, *var_list):
        for e in self.experiments:
            self.experiments[e].remove(*var_list)
        self.default_exp.remove(*var_list)
    # alias for remove
    def rm(self, *var_list):
        return self.remove(*var_list)
    
    # save the entire history
    def save(self):
        return self.snapshots
    
    # load the entire history; if blank, initialize the object
    def load(self, snapshots):
        self.experiments = {}
        self.snapshots = list(snapshots)
        self.sequence = 0 # increments with each snapshot

        if len(self.snapshots) == 0:
            return

        # find all unique experiments as of final snapshot
        self.sequence = self.snapshots[-1]['sequence'] + 1
        
        # create experiments as of the final sequence snapshot
        experiments = filter(lambda s: s['sequence'] == self.sequence - 1, self.snapshots)
        for e in experiments:
            e = e.copy()
            del(e['sequence'])
            ee = Experiment()
            ee.load(e)
            self.include(ee)

        return self
    
    # filter data with a function (perhaps a lambda), but protect against KeyError exceptions
    # Returns a generator object. Make a list with list() or iterate with a comprehension.
    # example:
    # p.filter(lambda e : e['Temperature'] < 1.5)
    #  output: a list we can load with Probe.load() that contains all the datasets < 1.5 kelvin
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
