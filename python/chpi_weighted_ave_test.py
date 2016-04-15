# -*- coding: utf-8 -*-

"""
Testbench for chpi_weighted_average functions

@author: jussi
"""


import numpy as np
import mne
from mne.datasets import sample
import sys
import matplotlib.pyplot as plt
from elekta_avg import Elekta_averager
#from chpi_weighted_average import chpi_snr_epochs, weighted_average_epochs



if sys.platform == 'win32':
    DATA_ROOT = 'C:/Users/HUS20664877/Dropbox/'
else:  # Linux
    DATA_ROOT = '/home/jussi/Dropbox/'


chpi_raw_fname = DATA_ROOT+'bad_203_am_raw.fif'
raw = mne.io.Raw(chpi_raw_fname, allow_maxshield=True)
#raw = mne.io.Raw('avg_setup_test.fif')

picks = mne.pick_types(raw.info, meg=True)
event_id=[1,4,8]
tmin, tmax = -0.2, 0.8
tmin, tmax = 0, 1
sfreq = raw.info['sfreq']

events = mne.find_events(raw, stim_channel='STI101', consecutive=True)
eav = Elekta_averager(raw.info['acq_pars'])

evs, ev_id = eav.events_from_mne_triggers(events)


""" Replace each trigger transition by the corresponding Elekta event(s) (if any). """
eav = Elekta_averager(raw.info['acq_pars'])
events = mne.find_events(raw, stim_channel='STI101', consecutive=True)
events_ = events.copy()
events_[:,2] = 0
for n,ev in eav.events.iteritems():
    pre_ok = np.bitwise_and(ev.oldmask, events[:,1]) == ev.oldbits
    post_ok = np.bitwise_and(ev.newbits, events[:,2]) == ev.newbits
    ok_ind = np.where(pre_ok & post_ok)
    if np.all(events_[ok_ind,2] == 0):
        events_[ok_ind,2] |= 1 << (n - 1)  # switch on the bit corresponding to event number
    """ Adjust for trigger-stimulus delay by delaying the event. """
    events[ok_ind,0] += ev.delay
        

""" Replace Elekta events by 'category triggers', i.e. times where conditions for averaging a given
category are fulfilled. This requires consideration of both ref. and req. events.  """
cat_triggers = events.copy()
times = events[:,0]
for n,cat in eav.categories.iteritems():
    refEvents_inds = np.where(events_[:,2] & (1 << cat.event-1))[0]  # indices of times where ref. event occurs
    refEvents_t = times[refEvents_inds]
    if cat.reqevent:
        reqEvents_inds = np.where(events_[:,2] & (1 << cat.reqevent-1))[0]  # indices of times where req. event occurs
        reqEvents_t = times[reqEvents_inds]
        # relative (to refevent) time window (in samples) where req. event must occur (e.g. [0 200])
        win = np.round(np.array(sorted([0, (-1)**(0-cat.reqwhen)*cat.reqwithin]))*sfreq)
        refEvents_wins = refEvents_t[:,None] + win
        req_acc = np.zeros(refEvents_inds.shape, dtype=bool)
        for t in reqEvents_t:
            # true for windows which have the given reqEvent in them
            reqEvent_in_win = np.logical_and(t >= refEvents_wins[:,0], t <= refEvents_wins[:,1])
            req_acc |= reqEvent_in_win
        # leave only ref. events where req. event condition is satisfied
        refEvents_inds = refEvents_inds[np.where(req_acc)]
        refEvents_t = times[refEvents_inds]            
    # adjust for trigger-stimulus delay by delaying the ref. event
    refEvents_t += eav.events[cat.event].delay
    



    
    

            
            





    if cat.reqevent:
        # indices where the req. events occur        
        reqEvents_inds = np.where(events_[:,2] & (1 << cat.reqevent-1))
        # relative (to refevent) time window (in samples) where req. event must occur (e.g. [0 200])
        win = np.round(np.array(sorted([0, (-1)**(0-cat.reqwhen)*cat.reqwithin]))*sfreq)
        # time window for each ref. event
        refEvents_wins = refEvents_t[:,None] + win[None,t]
        req_in_window = np.zeros(refEvents_t.shape)
        

"""        for t in refEvents_wins:


        reqEvent_ok = np.logical_and(
        
        -lines where single given req. event time x is in the window:
        np.logical_and(x >= refEvents_wins[:,0], x <= refEvents_wins[:,1])

    #for t in refEvents_t:
    #    if t-cat.reqwhen
        
    #cat_triggers[refEvents_t,2] = 1
   """ 
    
    






"""
for each cat:
    find where reqEvents and refEvents occur
    take times = t where refEvents occur and reqEvent occurs in given window
    write t and category code into array
produce event_id with category codes and comments


"""
    
    
            
        
    
    
    





chpi_epochs = mne.Epochs(raw, events, id, tmin=-.2, tmax=.8, baseline=None, picks=picks, preload=True)





""""
events array is of form:

[t0 x_pre x_post
t1 x_pre x_post
...
]

For each Elekta event, find array rows where:

OldMask & x_pre == OldBits
NewMask & x_post == NewBits





"""

def events_from_dacq_pars(data):
    """ Create events and event_id from Elekta data acquisition parameters. """
    
    events = mne.find_events(raw, stim_channel='STI101', consecutive=True)
    for ev in evs:
        premask = int(ev.OldMask)
        postmask = int(ev.NewMask)
        
    
    # -apply bitmasking from Elekta event categories
    # -compare with defined events
    # -produce array with [t 0 Eventcode], where Eventcode is the Elekta event number
    # (check whether this works with mne.Epochs)
    # -write corresponding event_id, taking dict keys from avg categories
    
    
    
    
    
    
    
    
    


#for id in event_id:
#    chpi_epochs = mne.Epochs(raw, events, id, tmin, tmax, baseline=(None, 0), picks=picks, preload=True)
#    snr = chpi_snr_epochs(chpi_epochs)


sys.exit()



""" Average according to per-epoch SNRs """
sss_raw_fname = '/home/jussi/Dropbox/bad_203_am_raw_sss.fif'
raw = mne.io.Raw(sss_raw_fname, allow_maxshield=True)
events = mne.find_events(raw, stim_channel='STI101')
picks = mne.pick_types(raw.info, meg=True)
event_id={'Eka': 1}
tmin, tmax = -0.2, 0.8
sss_epochs = mne.Epochs(raw, events, event_id, tmin, tmax, baseline=(None, 0), picks=picks, preload=True)
weights = snr_grad
ev = weighted_average_epochs(sss_epochs, weights)
ev.savgol_filter(60)
ev.plot()
plt.figure()
ev0 = sss_epochs.average()
ev0.savgol_filter(60)
ev0.plot()







