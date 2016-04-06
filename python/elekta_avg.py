# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 09:34:26 2016

Defines functions for getting averaging info (events and categories)
for Elekta TRIUX/Vectorview(?) systems.

@author: jussi
"""


class Elekta_event(object):
    """ Represents trigger event in Elekta system, as defined in acquisition settings.  """
    
    vars = ['Name', 'Channel', 'NewBits', 'OldBits', 'NewMask', 'OldMask', 'Delay', 'Comment']
    
    def __init__(self):
        pass

class Elekta_category(object):
    """ Represents averaging category in Elekta system, as defined in acquisition settings. """
    
    vars =  ['Comment', 'Display', 'Start', 'State', 'End', 'Event', 'Nave', 'ReqEvent', 
    'ReqWhen', 'ReqWithin',  'SubAve']

    def __init__(self):
        pass
    
class Elekta_averaging_params(object):
    """ Averaging parameters for Elekta system, as defined in acquisition settings. """

    vars = ['magMax', 'magMin', 'magNoise', 'magSlope', 'magSpike', 'megMax', 'megMin', 'megNoise',
            'megSlope', 'megSpike', 'ncateg', 'nevent', 'stimSource', 'triggerMap', 'update', 'version',
            'artefIgnore', 'averUpdate']    
    
    def __init__(self, acq_pars):
        acqdi = _acqpars_dict(acq_pars)
        for var in Elekta_averaging_params.vars:
            setattr(self, var, acqdi['ERF'+var])

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

def _acqpars_dict(acq_pars):
    """ Makes a dict from a string of acquisition parameters, which is info['acq_pars'] 
    for Elekta Vectorview/TRIUX systems. """
    return dict(_acqpars_gen(acq_pars))

def _events_from_acq_pars(acq_pars):
    acqdi = _acqpars_dict(acq_pars)
    ncateg = int(acqdi['ERFncateg'])
    events = {}
    for evnum in [str(x).zfill(2) for x in range(1,ncateg)]:  # '01', '02', etc.
        event = Elekta_event()
        for var in Elekta_event.vars:
            eventvar = lambda var, evnum=evnum: 'ERFevent'+var+evnum  # e.g. 'ERFEventName01'
            setattr(event, var, acqdi[eventvar(var)])
        events[evnum] = event
    return events

def _categories_from_acq_pars(acq_pars):
    acqdi = _acqpars_dict(acq_pars)
    nevent = int(acqdi['ERFnevent'])
    cats = {}
    for catnum in [str(x).zfill(2) for x in range(1,nevent)]:  # '01', '02', etc.
        cat = Elekta_category()
        for var in Elekta_category.vars:
            catvar = lambda var, catnum=catnum: 'ERFcat'+var+catnum  # e.g. 'ERFcatName01'
            setattr(cat, var, acqdi[catvar(var)])
        cats[catnum] = cat
    return cats

def elekta_averaging_info(acq_pars):
    """ Returns event and category dicts, given info.['acq_pars'] """
    return _events_from_acq_pars(acq_pars), _categories_from_acq_pars(acq_pars), Elekta_averaging_params(acq_pars)
    



