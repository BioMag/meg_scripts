# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 09:34:26 2016

Get averaging info (events and categories) defined in the acquisition parameters for 
Elekta TRIUX/Vectorview systems.

@author: jussi
"""


class Elekta_event(object):
    """ Represents trigger event in Elekta system.  """
    
    vars = ['Name', 'Channel', 'NewBits', 'OldBits', 'NewMask', 'OldMask', 'Delay', 'Comment']
    
    def __init__(self):
        pass

    def __repr__(self):
        return '<Elekta_event | Name: {} Comment: {} NewBits: {} OldBits: {} NewMask: {} OldMask: {} Delay: {}>'.format(
        self.Name, self.Comment, self.NewBits, self.OldBits, self.NewMask, self.OldMask, self.Delay)


class Elekta_category(object):
    """ Represents averaging category in Elekta system. """
    
    vars =  ['Comment', 'Display', 'Start', 'State', 'End', 'Event', 'Nave', 'ReqEvent', 
    'ReqWhen', 'ReqWithin',  'SubAve']

    def __init__(self):
        pass
    
    def __repr__(self):
        return '<Elekta_category | Comment: {} Event: {} ReqEvent: {} Start: {} End: {}>'.format(
        self.Comment, self.Event, self.ReqEvent, self.Start, self.End)
    

class Elekta_averager(object):
    """ Represents averager settings of Elekta TRIUX/Vectorview systems, including all defined
    trigger events and categories. 
    """

    # averager related variable names in dacq parameters
    vars = ['magMax', 'magMin', 'magNoise', 'magSlope', 'magSpike', 'megMax', 'megMin', 'megNoise',
            'megSlope', 'megSpike', 'ncateg', 'nevent', 'stimSource', 'triggerMap', 'update', 'version',
            'artefIgnore', 'averUpdate']    

    def __init__(self, acq_pars):
        """ acq_pars usually is obtained as Data.info['acq_pars'], where Data can be 
        instance of Raw, Epochs or Evoked. """
        self.acq_dict = Elekta_averager._acqpars_dict(acq_pars)
        for var in Elekta_averager.vars:
            setattr(self, var, self.acq_dict['ERF'+var])
        self.ncateg = int(self.ncateg)
        self.nevent = int(self.nevent)
        self.events = self._events_from_acq_pars()
        self.categories = self._categories_from_acq_pars()

    def __repr__(self):
        return '<Elekta_averager | Version: {} Categories: {} Events: {} StimSource: {}>'.format(
        self.version, self.ncateg, self.nevent, self.stimSource)

    @staticmethod
    def _acqpars_gen(acq_pars):
        """ Yields key/value pairs from a string of acquisition parameters. """
        acq_keys = ['ERF', 'DEF', 'ACQ', 'TCP']  # keys start with one of these
        for line in acq_pars.split():
            if any([line.startswith(x) for x in acq_keys]):
                key = line
                val = ''
            else:
                val += ' ' + line if val else line
            yield key, val

    @staticmethod
    def _acqpars_dict(acq_pars):
        """ Makes a dict from a string of acquisition parameters, which is info['acq_pars'] 
        for Elekta Vectorview/TRIUX systems. """
        return dict(Elekta_averager._acqpars_gen(acq_pars))

    def _events_from_acq_pars(self):
        events = {}
        for evnum in [str(x).zfill(2) for x in range(1,self.ncateg+1)]:  # '01', '02', etc.
            event = Elekta_event()
            for var in Elekta_event.vars:
                key = 'ERFevent'+var+evnum
                if key in self.acq_dict:
                    setattr(event, var, self.acq_dict[key])
                else:
                    raise Exception('Required key not in data acquisition parameters: '+key)
            events[int(evnum)] = event  # key dict by event number, that's how cats reference them
        return events

    def _categories_from_acq_pars(self, all_categories=False):
        cats = {}
        for catnum in [str(x).zfill(2) for x in range(1,self.nevent+1)]:  # '01', '02', etc.
            cat = Elekta_category()
            for var in Elekta_category.vars:
                key = 'ERFcat'+var+catnum
                if key in self.acq_dict:
                    setattr(cat, var, self.acq_dict[key])
                else:
                    raise Exception('Required key not in data acquisition parameters: '+key)
            if int(cat.State) == 1 or all_categories:  # category enabled
                cats[cat.Comment] = cat
        return cats
        
    def _events_in_use(self):
        eset = set()
        for cat in self.categories.values():
            eset.add(cat.Event)
            eset.add(cat.ReqEvent)
        eset.discard(u'0')  # indicates that no event for referred to
        return eset
            
    def event_in_use(self, event):
        return event.Name in self._events_in_use()
        
        




