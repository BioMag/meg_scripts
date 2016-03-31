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

from chpi_weighted_average import chpi_snr_epochs, weighted_average_epochs




""" Get the per-epoch SNRs """
chpi_raw_fname = '/home/jussi/Dropbox/bad_203_am_raw.fif'
raw = mne.io.Raw(chpi_raw_fname, allow_maxshield=True)
events = mne.find_events(raw, stim_channel='STI101', consecutive=False)  # only find 0-> transitions
picks = mne.pick_types(raw.info, meg=True)
event_id=[1,4,8]
tmin, tmax = -0.2, 0.8
tmin, tmax = 0, 1

for id in event_id:
    chpi_epochs = mne.Epochs(raw, events, id, tmin, tmax, baseline=(None, 0), picks=picks, preload=True)
    snr = chpi_snr_epochs(chpi_epochs)


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







