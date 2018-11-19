class Experiment:
    # name is the name of the experiment
    def __init__(self, name='', **kwargs):
        self.name = name
        self.vars = kwargs
    def __getitem__(self,v):
        return self.vars[v]
    def update(self, **kwargs):
        self.vars.update(kwargs)
    def remove(self, *var_list):
        for v in var_list:
            try:
                del(self.vars[v])
            except KeyError:
                pass

    
    def load(self,d):
        self.name = d['name']
        self.vars = {}
        self.vars.update(d)
        del(self.vars['name'])
    def save(self):
        d = self.vars.copy()
        d.update({'name':self.name})
        return d
            
        
class Probe:
    def __init__(self, snapshots=[]):
        self.load(snapshots)
        
    def add(self,*experiments):
        # Add to experiment list. Experiments are indexed by name.
        for e in experiments:
            self.experiments[e.name] = e
    def remove(self,*experiment_names):
        for e in experiment_names:
            del(self.experiments[e])
        
    def __getitem__(self,experiment_name):
        return self.experiments[experiment_name]
    
    # update child experiments with new or changed variables
    def update(self, **kwargs):
        # update child experiments
        for e in self.experiments:
            self.experiments[e].update(**kwargs)
            
    # update and snapshot in one step
    def snap(self, **kwargs):
        self.update(**kwargs)
        for e in self.experiments:
            exp_copy = self.experiments[e].save()
            exp_copy.update({'sequence':self.sequence})
            self.snapshots.append(exp_copy)
        self.sequence += 1
        
    # remove a variable from child experiments
    def rm_var(self, *var_list):
        for e in self.experiments:
            self.experiments[e].remove(*var_list)
    
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
            self.add(ee)
    
    # filter data with a function (perhaps a lambda), but protect against KeyError exceptions
    # Returns a generator object. Make a list with list() or iterate with a comprehension.
    # example:
    # p.filter(lambda e : e['Temperature'] < 1.5)
    #  output: a list we can load with Probe.load() that contains all the datasets < 1.5 kelvin
    def filter(self, condition):
        snaps = []
        for e in self.save():
            try:
                if condition(e):
                    yield e
            except:
                pass

