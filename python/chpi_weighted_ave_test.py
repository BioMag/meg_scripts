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
from elekta_avg import elekta_averaging_info
#from chpi_weighted_average import chpi_snr_epochs, weighted_average_epochs




""" Get the per-epoch SNRs """
chpi_raw_fname = '/home/jussi/Dropbox/bad_203_am_raw.fif'
raw = mne.io.Raw(chpi_raw_fname, allow_maxshield=True)

picks = mne.pick_types(raw.info, meg=True)
event_id=[1,4,8]
tmin, tmax = -0.2, 0.8
tmin, tmax = 0, 1

evs, cats, av_info = elekta_averaging_info(raw)

events = mne.find_events(raw, stim_channel='STI101', consecutive=False)

chpi_epochs = mne.Epochs(raw, events, id, tmin=args.epoch_start, tmax=args.epoch_end, baseline=None, picks=picks, preload=True, verbose=verbose)


""""
[t x_pre x_post]
"""

def events_from_dacq_pars(data):
    """ Create events and event_id from Elekta data acquisition parameters. """
    evs, cats, av_info = elekta_averaging_info(data)
    events = mne.find_events(raw, stim_channel='STI101', consecutive=True)
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







