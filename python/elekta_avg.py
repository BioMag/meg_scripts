# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 09:34:26 2016

Get averaging info (events and categories) defined in the acquisition parameters for 
Elekta TRIUX/Vectorview systems.

@author: jussi
"""

import numpy as np
import mne


class Elekta_event(object):
    """ Represents trigger event in Elekta/Neuromag system.  """
    
    # original dacq variable names
    vars = ['Name', 'Channel', 'NewBits', 'OldBits', 'NewMask', 'OldMask', 'Delay', 'Comment']
    
    def __init__(self, name, channel, newbits, oldbits, newmask, oldmask, delay, comment):
        self.name = name  # short name for event
        self.channel = channel  # stimulus channel
        self.newbits = int(newbits)  # state after trigger transition
        self.oldbits = int(oldbits)  # state before trigger transition
        self.newmask = int(newmask)  # bitmask for post-transition state. 0-bits = ignore
        self.oldmask = int(oldmask)  # bitmask for pre-transition state
        self.delay = float(delay)  # delay to stimulus (s)
        self.comment = comment  # verbose comment for the event
        # non-dacq vars
        self.in_use = False  # whether event is in use, i.e. referred to by a category

    def __repr__(self):
        return '<Elekta_event | name: {} comment: {} new bits: {} old bits: {} new mask: {} old mask: {} delay: {}>'.format(
        self.name, self.comment, self.newbits, self.oldbits, self.newmask, self.oldmask, self.delay)


class Elekta_category(object):
    """ Represents averaging category in Elekta/Neuromag system. """

    # original dacq variable names    
    vars =  ['Comment', 'Display', 'Start', 'State', 'End', 'Event', 'Nave', 'ReqEvent', 
    'ReqWhen', 'ReqWithin',  'SubAve']

    def __init__(self, comment, display, start, state, end, event, nave, reqevent, reqwhen, reqwithin, subave):
        self.comment = comment
        self.display = True if display==u'1' else False  # whether do display online in dacq
        self.state = True if state==u'1' else False  # enabled or not
        self.start = float(start)  # epoch start
        self.end = float(end)  # epoch end
        self.nave = int(nave)  # desired n of averages
        self.event = int(event)  # the reference event (index to event list)
        self.reqevent = int(reqevent)  # additional required event (index to event list)
        self.reqwhen = int(reqwhen)  # 1=before, 0=after ref. event
        self.reqwithin = float(reqwithin)  # time window before or after ref. event (s)
        self.subave = int(subave)  # subaverage size
        # non-dacq vars
        self.times = None  # list of t=0 times (in samples) for corresponding epochs
    
    def __repr__(self):
        return '<Elekta_category | comment: \'{}\' event: {} reqevent: {} reqwhen: {} reqwithin: {} start: {} end: {}>'.format(
        self.comment, self.event, self.reqevent, self.reqwhen, self.reqwithin, self.start, self.end)
    

class Elekta_averager(object):
    """ Represents averager settings of Elekta TRIUX/Vectorview systems, including all defined
    trigger events and categories. 
    """

    # these are DACQ averager related variable names without preceding 'ERF'
    vars = ['magMax', 'magMin', 'magNoise', 'magSlope', 'magSpike', 'megMax', 'megMin', 'megNoise',
            'megSlope', 'megSpike', 'eegMax', 'eegMin', 'eegNoise', 'eegSlope', 'eegSpike', 'eogMax',
            'ecgMax', 'ncateg', 'nevent', 'stimSource', 'triggerMap', 'update', 'version',
            'artefIgnore', 'averUpdate']    

    acq_var_magic = ['ERF', 'DEF', 'ACQ', 'TCP']  # dacq variable names start with one of these

    def __init__(self, acq_pars):
        """ acq_pars usually is obtained as Data.info['acq_pars'], where Data can be 
        instance of Raw, Epochs or Evoked. """
        self.acq_dict = Elekta_averager._acqpars_dict(acq_pars)
        # sets instance variables (lowercase versions of dacq variable names), to avoid a lot of boilerplate code
        for var in Elekta_averager.vars:
            val = self.acq_dict['ERF'+var]
            if var[:3] in ['mag','meg','eeg','eog','ecg']:  # rejection criteria are floats
                val = float(val)
            elif var in ['ncateg','nevent']:
                val = int(val)
            setattr(self, var.lower(), val)
        self.stimsource = u'Internal' if self.stimsource == u'1' else u'External'
        # collect events and categories as instance dicts
        self.events = self._events_from_acq_pars()
        self.categories = self._categories_from_acq_pars()
        # tag the events that are actually in use
        for cat in self.categories.values():
            if cat.event:
                self.events[cat.event].in_use = True
            if cat.reqevent:
                self.events[cat.reqevent].in_use = True

    def __repr__(self):
        return '<Elekta_averager | Version: {} Categories: {} Events: {} Stim source: {}>'.format(
        self.version, self.ncateg, self.nevent, self.stimsource)

    @staticmethod
    def _acqpars_gen(acq_pars):
        """ Yields key/value pairs from a string of acquisition parameters. """
        for line in acq_pars.split():
            if any([line.startswith(x) for x in Elekta_averager.acq_var_magic]):
                key = line
                val = ''
            else:
                val += ' ' + line if val else line  # DACQ splits items with spaces into multiple lines
            yield key, val

    @staticmethod
    def _acqpars_dict(acq_pars):
        """ Makes a dict from a string of acquisition parameters, which is info['acq_pars'] 
        for Elekta Vectorview/TRIUX systems. """
        return dict(Elekta_averager._acqpars_gen(acq_pars))

    def _events_from_acq_pars(self):
        """ Collects DACQ defined events into a dict. """
        events = {}
        for evnum in [str(x).zfill(2) for x in range(1,self.ncateg+1)]:  # '01', '02', etc.
            evdi = {}
            for var in Elekta_event.vars:
                acq_key = 'ERFevent'+var+evnum  # name of dacq variable, e.g. 'ERFeventNewBits01'
                class_key = var.lower()         # corresponding instance variable, e.g. 'newbits'
                evdi[class_key] = self.acq_dict[acq_key]
            events[int(evnum)] = Elekta_event(**evdi)  # events are keyed by number starting from 1
        return events

    def _categories_from_acq_pars(self, all_categories=False):
        """ Collects DACQ averaging categories into a dict. """
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
        
    def _mne_events_to_dacq(self, mne_events):
        """ Creates list of dacq events based on mne trigger transitions list.
        mne_events is typically given by mne.find_events (use consecutive=True
        to get all transitions). Output consists of rows in the form [t, 0, event_code]
        where t is time in samples and event_code is all events compatible with the transition
        bitwise anded together. """
        events_ = mne_events.copy()
        events_[:,1] = 0
        events_[:,2] = 0
        for n,ev in self.events.iteritems():
            if ev.in_use:
                pre_ok = np.bitwise_and(ev.oldmask, mne_events[:,1]) == ev.oldbits
                post_ok = np.bitwise_and(ev.newbits, mne_events[:,2]) == ev.newbits
                ok_ind = np.where(pre_ok & post_ok)
                if np.all(events_[ok_ind,2] == 0):
                    events_[ok_ind,2] |= 1 << (ev.number - 1)  # switch on the bit corresponding to event number
                else:
                    raise Exception('Multiple dacq events match trigger transition')
        return events_

    def _mne_events_to_category_times(self, cat, mne_events, sfreq):
        """ Translate mne events to times when epochs for a given category should
        be collected. """
        events = self._mne_events_to_dacq(mne_events)
        times = events[:,0]
        refEvents_inds = np.where(events[:,2] & (1 << cat.event-1))[0]  # indices of times where ref. event occurs
        refEvents_t = times[refEvents_inds]
        catEvents_inds = refEvents_inds
        catEvents_t = times[catEvents_inds]
        if cat.reqevent:
            reqEvents_inds = np.where(events[:,2] & (1 << cat.reqevent-1))[0]  # indices of times where req. event occurs
            reqEvents_t = times[reqEvents_inds]
            # relative (to refevent) time window (in samples) where req. event must occur (e.g. [0 200])
            win = np.round(np.array(sorted([0, (-1)**(cat.reqwhen)*cat.reqwithin]))*sfreq)
            refEvents_wins = refEvents_t[:,None] + win
            req_acc = np.zeros(refEvents_inds.shape, dtype=bool)
            for t in reqEvents_t:
                # true for windows which have the given reqEvent in them
                reqEvent_in_win = np.logical_and(t >= refEvents_wins[:,0], t <= refEvents_wins[:,1])
                req_acc |= reqEvent_in_win
            # leave only ref. events where req. event condition is satisfied
            catEvents_inds = catEvents_inds[np.where(req_acc)]
            catEvents_t = times[catEvents_inds]            
        # adjust for trigger-stimulus delay by delaying the ref. event
        catEvents_t += np.round(self.events[cat.event].delay * sfreq)
        return catEvents_t

    def get_mne_rejection_dict(self):
        """ Makes a mne rejection dict based on the averager parameters. Result
        can be used e.g. with mne.Epochs """
        return {'grad':self.megmax, 'mag':self.magmax, 'eeg':self.eegmax,
                'eog':self.eogmax, 'ecg':self.ecgmax}

    def get_epochs(self, raw, catname, picks=None, reject=None, stim_channel=None, mask=0):
        """ Get mne.Epochs instance corresponding to the given category. """
        cat = self.categories[catname]
        mne_events = mne.find_events(raw, stim_channel=stim_channel, mask=mask, consecutive=True)
        sfreq = raw.info['sfreq']
        cat_t = self._mne_events_to_category_times(cat, mne_events, sfreq)
        catev = np.c_[cat_t, np.zeros(cat_t.shape), np. ones(cat_t.shape)].astype(np.uint32)
        id = {cat.comment: 1}
        return mne.Epochs(raw, catev, event_id=id, reject=reject, tmin=cat.start, tmax=cat.end, baseline=None,
                          picks=picks, preload=True)

    def _event_id_dict(self):
        """ Returns a simple event id dict that can be used with mne.find_events.
        Keys are category comments and values are the corresponding ref. event numbers. """
        catdi = {}
        for cat in self.categories.values():
            if cat.state == 1:
                catdi[cat.comment] = cat.event
            if cat.reqevent:
                raise ValueError('Req. events not supported yet')
        return catdi
        
        




