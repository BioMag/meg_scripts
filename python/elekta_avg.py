# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 09:34:26 2016

Get averaging info (events and categories) defined in the acquisition parameters for 
Elekta TRIUX/Vectorview systems.

TODO: adapt event & cat creation code for new objects

@author: jussi
"""


class Elekta_event(object):
    """ Represents trigger event in Elekta system.  """
    
    # original dacq variable names
    vars = ['Name', 'Channel', 'NewBits', 'OldBits', 'NewMask', 'OldMask', 'Delay', 'Comment']
    
    def __init__(self, name, channel, newbits, oldbits, newmask, oldmask, delay, comment):
        self.name = name
        self.number = int(name)     # for convenience; dacq variable 'Name' actually represents event number 1..32
        self.channel = channel
        self.newbits = int(newbits)
        self.oldbits = int(oldbits)
        self.newmask = int(newmask)
        self.oldmask = int(oldmask)
        self.delay = int(delay)
        self.comment = comment

    def __repr__(self):
        return '<Elekta_event | Name: {} Comment: {} NewBits: {} OldBits: {} NewMask: {} OldMask: {} Delay: {}>'.format(
        self.name, self.comment, self.newbits, self.oldbits, self.newmask, self.oldmask, self.delay)


class Elekta_category(object):
    """ Represents averaging category in Elekta system. """
    
    vars =  ['Comment', 'Display', 'Start', 'State', 'End', 'Event', 'Nave', 'ReqEvent', 
    'ReqWhen', 'ReqWithin',  'SubAve']

    def __init__(self, comment, display, start, state, end, event, nave, reqevent, reqwhen, reqwithin, subave):
        # TODO: some typecasts
        self.comment = comment
        self.display = display
        self.start = start
        self.state = state
        self.end = end
        self.event = event
        self.nave = nave
        self.reqevent = reqevent
        self.reqwhen = reqwhen
        self.reqwithin = reqwithin
        self.subave = subave
    
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
            for var in Elekta_event.vars:
                tdi = {}
                acq_key = 'ERFevent'+var+evnum
                if key in self.acq_dict:
                    clsvar = var.lower()
                    
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
        
        




