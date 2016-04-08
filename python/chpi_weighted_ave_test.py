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

picks = mne.pick_types(raw.info, meg=True)
event_id=[1,4,8]
tmin, tmax = -0.2, 0.8
tmin, tmax = 0, 1

eav = Elekta_averager(raw.info['acq_pars'])
events = mne.find_events(raw, stim_channel='STI101', consecutive=True)

""" Transform events into corresponding list of Elekta events. """
events_ = events.copy()
events_[:,2] = 0
for n,ev in eav.events.iteritems():
    pre_ok = np.bitwise_and(int(ev.OldMask), events[:,1]) == int(ev.OldBits)
    post_ok = np.bitwise_and(int(ev.NewMask), events[:,2]) == int(ev.NewBits)
    ok_ind = np.where(pre_ok & post_ok)
    if np.all(events_[ok_ind,2] == 0):
        events_[ok_ind,2] |= 1 << (int(ev.Name) - 1)  # switch on the bit corresponding to event number

""" Transform Elekta events into averaging categories. """
for n,cat in eav.categories.iteritems():
    refEvents_t = np.where(events_[:,2] & (1 << int(cat.Event)-1))
    reqEvents_t = np.where(events_[:,2] & (1 << int(cat.ReqEvent)-1))
    






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







