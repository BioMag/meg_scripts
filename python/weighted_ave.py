# -*- coding: utf-8 -*-

"""
Weighted averager of epochs.

@author: jussi
"""


import numpy as np
import mne
from mne import io
from mne.io.constants import FIFF
from mne.datasets import sample
import sys
import matplotlib.pyplot as plt




def chpi_freqs(info):
    """ Get chpi frequencies from info dict (e.g. raw.info) """
    for coil in info['hpi_meas'][0]['hpi_coils']:
        yield coil['coil_freq'][0]


def chpi_snr_epochs(epochs, n_lineharm=2, channels='grad', hpi_coil='median'):
    """ Return estimated continuous HPI SNR for each epoch in epochs
    (mne.Epochs object). SNR estimation is done by fitting a GLM to the data,
    comparing to residual (unexplained) variance and averaging the resulting
    SNR across channels.
    Parameters:
    epochs: mne.epochs instance with a single condition
    n_lineharm: number of line frequency harmonics to use in the GLM
    channels: 'grad', 'mag': which channels to consider in SNR averaging
    hpi_coil: 'best', 'worst', 'median': which hpi coil to track SNR for.
    'best' selects best coil at each epoch. 'median' selects the coil with
    the median SNR value at each epoch. 'worst' selects the coil with lowest
    SNR at each epoch.
    
    TODO:
    handling of bad channels?
    """

    if len(epochs.event_id) > 1:
        raise ValueError('Epochs object should contain one category only')

    alldata = epochs.pick_types(meg=True).get_data()    
    nepochs, nchan, buflen = alldata.shape
    sfreq = epochs.info['sfreq']
    t = np.linspace(0, buflen/sfreq, endpoint=False, num=buflen)
    cfreqs = list(chpi_freqs(epochs.info))
    ncoils = len(cfreqs)
    linefreq = raw.info['line_freq']
    linefreqs = (np.arange(n_lineharm+1)+1) * linefreq
    
    # gradiometer and magmeter indices
    pick_meg = mne.pick_types(epochs.info, meg=True)
    pick_mag = mne.pick_types(epochs.info, meg='mag')
    pick_grad = mne.pick_types(epochs.info, meg='grad')

    # create linear model    
    model = np.c_[t, np.ones(t.shape)]  # model slope and DC
    for f in list(linefreqs)+cfreqs:  # add sine and cosine term for each freq
        model = np.c_[model, np.cos(2*np.pi*f*t), np.sin(2*np.pi*f*t)]
    inv_model = np.linalg.pinv(model)

    # loop through epochs
    snr_avg_grad = np.zeros([ncoils, nepochs])
    snr_avg_mag = np.zeros([ncoils, nepochs])
    resid_vars = np.zeros([nchan, nepochs])

    for epn in range(nepochs):
        epdata = alldata[epn,:,:].transpose()
        coeffs = np.dot(inv_model, epdata)
        print(coeffs)
        coeffs_hpi = coeffs[2 + 2*len(linefreqs):]
        resid_vars[:,epn] = np.var(epdata - np.dot(model, coeffs), 0)
        # get total hpi amplitudes by combining sine and cosine terms
        hpi_amps = np.sqrt(coeffs_hpi[0::2,:]**2 + coeffs_hpi[1::2,:]**2)
        # divide average HPI power by average variance
        snr_avg_grad[:,epn] = np.divide((hpi_amps**2/2)[:,pick_grad].mean(1),resid_vars[pick_grad,epn].mean())
        snr_avg_mag[:,epn] = np.divide((hpi_amps**2/2)[:,pick_mag].mean(1),resid_vars[pick_mag,epn].mean())

    if hpi_coil == 'median':
        fun = np.median
    elif hpi_coil == 'best':
        fun = np.max
    elif hpi_coil == 'worst':
        fun = np.min        
    else:
        raise ValueError('Invalid HPI coil selection parameter')        

    if channels == 'grad':
        snr = fun(snr_avg_grad, axis=0)
    elif channels == 'mag':
        snr = fun(snr_avg_mag, axis=0)
    
    return snr


def weighted_average_epochs(epochs, weights):
    """  Compute weighted average of epochs. epochs is a mne.Epochs object.
    weights is a list or 1-d numpy array with leading dim of n_epochs.
    The weighted average is divided by the sum of weights. """
  
    weights = np.array(weights)  # will accept either list or numpy array
    n_epochs = len(epochs)
    if not len(weights) == n_epochs:
        raise ValueError('Need as many weights as epochs')
    w_ = weights.squeeze()[:,np.newaxis,np.newaxis]  # reshape for broadcasting
    epw = epochs.get_data() * w_ / np.sum(w_) # normalize
    epw_av = np.sum(epw, axis=0)
    return epochs._evoked_from_epoch_data(epw_av, epochs.info, None, n_epochs, FIFF.FIFFV_ASPECT_AVERAGE)
    return mne.EpochsArray


""" Get the per-epoch SNRs """
chpi_raw_fname = '/home/jussi/Dropbox/bad_203_am_raw.fif'
raw = io.Raw(chpi_raw_fname, allow_maxshield=True)
events = mne.find_events(raw, stim_channel='STI101', consecutive=False)  # only find 0-> transitions
picks = mne.pick_types(raw.info, meg=True)
event_id=1
tmin, tmax = -0.2, 0.8
tmin, tmax = 0, 1
chpi_epochs = mne.Epochs(raw, events, event_id, tmin, tmax, baseline=(None, 0), picks=picks, preload=True)

c

snr_grad, snr_mag = chpi_snr_epochs(chpi_epochs)



""" Average according to per-epoch SNRs """
sss_raw_fname = '/home/jussi/Dropbox/bad_203_am_raw_sss.fif'
raw = io.Raw(sss_raw_fname, allow_maxshield=True)
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







