# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 09:34:26 2016

Get averaging info (events and categories) defined in the acquisition parameters for 
Elekta TRIUX/Vectorview systems.

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
        return '<Elekta_event | name: {} comment: {} new bits: {} old bits: {} new mask: {} old mask: {} delay: {}>'.format(
        self.name, self.comment, self.newbits, self.oldbits, self.newmask, self.oldmask, self.delay)


class Elekta_category(object):
    """ Represents averaging category in Elekta system. """

    # original dacq variable names    
    vars =  ['Comment', 'Display', 'Start', 'State', 'End', 'Event', 'Nave', 'ReqEvent', 
    'ReqWhen', 'ReqWithin',  'SubAve']

    def __init__(self, comment, display, start, state, end, event, nave, reqevent, reqwhen, reqwithin, subave):
        self.comment = comment
        self.display = True if display==u'1' else False
        self.state = True if state==u'1' else False
        self.start = float(start)
        self.end = float(end)
        self.nave = int(nave)
        self.event = int(event)  # categories are referred to by numbers
        self.reqevent = int(reqevent)
        self.reqwhen = float(reqwhen)
        self.reqwithin = float(reqwithin)
        self.subave = True if subave==u'1' else False
    
    def __repr__(self):
        return '<Elekta_category | comment: {} event: {} reqevent: {} reqwhen: {} reqwithin: {} start: {} end: {}>'.format(
        self.comment, self.event, self.reqevent, self.reqwhen, self.reqwithin, self.start, self.end)
    

class Elekta_averager(object):
    """ Represents averager settings of Elekta TRIUX/Vectorview systems, including all defined
    trigger events and categories. 
    """

    # averager related variable names in dacq parameters
    vars = ['magMax', 'magMin', 'magNoise', 'magSlope', 'magSpike', 'megMax', 'megMin', 'megNoise',
            'megSlope', 'megSpike', 'ncateg', 'nevent', 'stimSource', 'triggerMap', 'update', 'version',
            'artefIgnore', 'averUpdate']    

    acq_var_magic = ['ERF', 'DEF', 'ACQ', 'TCP']  # dacq variable names start with one of these

    def __init__(self, acq_pars):
        """ acq_pars usually is obtained as Data.info['acq_pars'], where Data can be 
        instance of Raw, Epochs or Evoked. """
        self.acq_dict = Elekta_averager._acqpars_dict(acq_pars)
        # sets instance variables (lowercase versions of dacq variable names), to avoid a lot of boilerplate code
        for var in Elekta_averager.vars:
            val = self.acq_dict['ERF'+var]
            if var[:3] in ['mag','meg']:
                val = float(val)
            elif var in ['ncateg','nevent']:
                val = int(val)
            setattr(self, var.lower(), val)
        self.events = self._events_from_acq_pars()
        self.categories = self._categories_from_acq_pars()

    def __repr__(self):
        return '<Elekta_averager | Version: {} Categories: {} Events: {} StimSource: {}>'.format(
        self.version, self.ncateg, self.nevent, self.stimSource)

    @staticmethod
    def _acqpars_gen(acq_pars):
        """ Yields key/value pairs from a string of acquisition parameters. """
        for line in acq_pars.split():
            if any([line.startswith(x) for x in Elekta_averager.acq_var_magic]):
                key = line
                val = ''
            else:
                val += ' ' + line if val else line  # dacq splits items with spaces into multiple lines
            yield key, val

    @staticmethod
    def _acqpars_dict(acq_pars):
        """ Makes a dict from a string of acquisition parameters, which is info['acq_pars'] 
        for Elekta Vectorview/TRIUX systems. """
        return dict(Elekta_averager._acqpars_gen(acq_pars))

    def _events_from_acq_pars(self):
        events = {}
        for evnum in [str(x).zfill(2) for x in range(1,self.ncateg+1)]:  # '01', '02', etc.
            evdi = {}
            for var in Elekta_event.vars:
                acq_key = 'ERFevent'+var+evnum  # name of dacq variable, e.g. 'ERFeventNewBits01'
                class_key = var.lower()         # corresponding instance variable, e.g. 'newbits'
                evdi[class_key] = self.acq_dict[acq_key]
            events[int(evdi['name'])] = Elekta_event(**evdi)  # events are keyed by number starting from 1
        return events

    def _categories_from_acq_pars(self, all_categories=False):
        cats = {}
        for catnum in [str(x).zfill(2) for x in range(1,self.nevent+1)]:  # '01', '02', etc.
            catdi = {}
            for var in Elekta_category.vars:
                acq_key = 'ERFcat'+var+catnum
                class_key = var.lower()
                catdi[class_key] = self.acq_dict[acq_key]
            if int(catdi['state']) == 1 or all_categories:  # category enabled
                cats[catdi['comment']] = Elekta_category(**catdi)
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
        
        




