# -*- coding: utf-8 -*-
"""
Weighted averaging of epochs according to SNR of continuous HPI (cHPI) signal.


@author: Jussi Nurminen (jnu@iki.fi)
"""


from __future__ import print_function
import matplotlib.pyplot as plt
import mne
import numpy as np
import argparse
import os
from mne import AcqParserFIF
from mne.chpi import _get_hpi_info


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

    """

    if len(epochs.event_id) > 1:
        raise ValueError('Epochs object should contain one category only')

    epochs.load_data()
    alldata = epochs.pick_types(meg=True).get_data()
    nepochs, nchan, buflen = alldata.shape
    sfreq = epochs.info['sfreq']
    t = np.linspace(0, buflen/sfreq, endpoint=False, num=buflen)
    cfreqs = list(chpi_freqs(epochs.info))
    # this cannot be used yet, see
    # https://github.com/mne-tools/mne-python/pull/3670
    # (cfreqs, _, _, _, _) = _get_hpi_info(epochs.info)
    ncoils = len(cfreqs)
    linefreq = epochs.info['line_freq']
    linefreqs = (np.arange(n_lineharm+1)+1) * linefreq

    # gradiometer and magmeter indices
    pick_meg = mne.pick_types(epochs.info, meg=True)
    pick_mag = mne.pick_types(epochs.info, meg='mag')
    pick_grad = mne.pick_types(epochs.info, meg='grad')

    # create linear model
    model = np.c_[t, np.ones(t.shape)]  # model slope and DC
    # add sine and cosine term for each freq
    for f in list(linefreqs)+list(cfreqs):
        model = np.c_[model, np.cos(2*np.pi*f*t), np.sin(2*np.pi*f*t)]
    inv_model = np.linalg.pinv(model)

    # loop through epochs
    snr_avg_grad = np.zeros([ncoils, nepochs])
    snr_avg_mag = np.zeros([ncoils, nepochs])
    resid_vars = np.zeros([nchan, nepochs])

    for epn in range(nepochs):
        epdata = alldata[epn, :, :].transpose()
        coeffs = np.dot(inv_model, epdata)
        coeffs_hpi = coeffs[2 + 2*len(linefreqs):]
        resid_vars[:, epn] = np.var(epdata - np.dot(model, coeffs), 0)
        # get total hpi amplitudes by combining sine and cosine terms
        hpi_amps = np.sqrt(coeffs_hpi[0::2, :]**2 + coeffs_hpi[1::2, :]**2)
        # divide average HPI power by average variance
        snr_avg_grad[:, epn] = np.divide((hpi_amps**2/2)[:, pick_grad].mean(1),
                                         resid_vars[pick_grad, epn].mean())
        snr_avg_mag[:, epn] = np.divide((hpi_amps**2/2)[:, pick_mag].mean(1),
                                        resid_vars[pick_mag, epn].mean())

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


def weigh_epochs(epochs, weights):
    """ Weigh epochs in mne.Epochs object (elementwise multiply by
    weights vector). """
    epochs.load_data()
    weights = np.array(weights)  # will accept either list or numpy array
    n_epochs = len(epochs)
    if not len(weights) == n_epochs:
        raise ValueError('Need as many weights as epochs')
    # reshape for broadcasting
    w_ = weights.squeeze()[:, np.newaxis, np.newaxis]
    # normalize by n_epochs/sum(weights), so that averaging
    # will yield a correct result
    epochs._data *= n_epochs * w_ / np.sum(w_)


if __name__ == '__main__':

    """ Default parameters """
    default_nharm = 2  # number of line harmonics to include
    default_stim_channel = 'STI101'

    """ Parse command line """
    parser = argparse.ArgumentParser(description='Weighted averaging of '
                                                 'epochs according to cHPI '
                                                 'SNR.')
    help_snr_file = ('Name of raw fiff file. cHPI SNR will be computed from '
                     'this file. Epochs for averaging will be taken from '
                     'the same file, unless --epochs_file is specified.')
    help_epochs_file = ('Raw fiff file to get the epochs from. It must '
                        'have the same epochs/categories as snr_file.')
    help_reject = ('Whether to use the rejection parameters defined in '
                   'the data acquisition when averaging epochs.')
    help_nharm = 'Number of line frequency harmonics to include.'
    help_epoch_start = 'Epoch start relative to trigger (s)'
    help_epoch_end = 'Epoch end relative to trigger (s)'
    help_plot_snr = 'Whether to plot SNR or not'
    help_sti_mask = 'Mask to apply to the stim channel'
    parser.add_argument('snr_file', help=help_snr_file)
    parser.add_argument('--epochs_file', type=str, default=None,
                        help=help_epochs_file)
    parser.add_argument('--reject', help=help_reject, action='store_true')
    parser.add_argument('--nharm', type=int, default=default_nharm,
                        choices=[0, 1, 2, 3, 4], help=help_nharm)
    parser.add_argument('--epoch_start', type=float, default=None,
                        help=help_epoch_start)
    parser.add_argument('--epoch_end', type=float, default=None,
                        help=help_epoch_end)
    parser.add_argument('--plot_snr', help=help_plot_snr, action='store_true')
    parser.add_argument('--sti_mask', type=int, default=0,
                        help=help_sti_mask)

    args = parser.parse_args()

    if args.epochs_file:
        fnbase = os.path.basename(os.path.splitext(args.epochs_file)[0])
    else:
        fnbase = os.path.basename(os.path.splitext(args.snr_file)[0])
    verbose = False

    """ Load cHPI SNR file. It is typically not maxfiltered, so ignore
    MaxShield warnings. """
    raw_chpi = mne.io.Raw(args.snr_file, allow_maxshield='yes',
                          verbose=verbose)
    picks = mne.pick_types(raw_chpi.info, meg=True)

    """ If using a separate file for the actual data epochs, load it. """
    if args.epochs_file:
        raw_epochs = mne.io.Raw(args.epochs_file, allow_maxshield=True,
                                verbose=verbose)

    """ Get averaging parameters. These should be identical for the SNR and
    epochs files. """
    ap = AcqParserFIF(raw_chpi.info)

    """ For each category, compute the SNR weights and the weighted
    average. """
    evokeds = []
    for cat in ap.categories:
        print('\nProcessing category: %s' % cat['comment'])
        print('Loading epochs for cHPI SNR...')
        cond = ap.get_condition(raw_chpi, cat)
        if args.epoch_start:
            epoch_start = args.epoch_start
            epoch_end = args.epoch_end
        else:
            epoch_start = cond['tmin']
            epoch_end = cond['tmax']
        if args.reject:
            print('Warning: epoch rejection not implemented yet')
            # reject = ap.reject if args.reject else None
            # flat = ap.flat if args.reject else None
            reject, flat = None, None
        chpi_epochs = mne.Epochs(raw_chpi, reject=reject, flat=flat, **cond)
        print('Computing SNR...')
        w_snr = chpi_snr_epochs(chpi_epochs, n_lineharm=args.nharm,
                                channels='grad', hpi_coil='median')
        if not args.epochs_file:
            data_epochs = chpi_epochs
        else:
            print('Loading epochs for averaging...')
            data_epochs = mne.Epochs(raw_epochs, reject=reject, flat=flat,
                                     **cond)
        print('Computing weighted average...')
        weigh_epochs(data_epochs, w_snr)
        evoked = data_epochs.average()
        evoked.comment = cat['comment']
        evokeds.append(evoked)

    if args.plot_snr:
        plt.figure()
        plt.plot(20*np.log10(w_snr))
        plt.xlabel('Epoch n')
        plt.ylabel('SNR (dB)')
        plt.show()

    """ Write all resulting evoked objects to a fiff file. """
    fn = fnbase + '_chpi_weighted-ave.fif'
    print('Saving', fn)
    mne.write_evokeds(fn, evokeds)
