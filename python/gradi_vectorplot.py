# -*- coding: utf-8 -*-
"""
gradi_vectorplot.py

Combine gradiometer data from evoked responses into a pairwise vector sum. 
Plot combined responses in topographic plot (XPlotter style).

NOTE: run with ipython --matplotlib -i -- gradi_vectorplot_py [args]
so that single-channel plot (left click on topographic display) works properly.



@author: jussi (jnu@iki.fi)
"""


from __future__ import print_function


import sys
import mne
import argparse
import mne.viz
import numpy as np
from scipy import signal


# params
BORD = 5

# parse command line
parser = argparse.ArgumentParser()
parser.add_argument('files', metavar='evoked_file',
                    help='Name of evoked fiff file', nargs='+')
parser.add_argument('--lowpass', type=float, metavar='f',
                    default=None, help='Lowpass frequency (Hz)')
args = parser.parse_args()


# read evoked data
evokeds = []
for file in args.files:
    evokeds_ = mne.read_evokeds(file)
    evokeds.extend(evokeds_)

# preprocess
for evoked in evokeds:

    # filter
    if args.lowpass:
        data = evoked.data
        sfreq = evoked.info['sfreq']
        lpfreqn = 2 * np.array(args.lowpass) / sfreq
        b, a = signal.butter(BORD, lpfreqn)
        evoked.data = signal.filtfilt(b, a, data)

    # get peak amplitude - TODO: does not work
    pch, plat = evoked.get_peak(ch_type='grad')
    print('%s peak amplitude: channel pair %s, latency %.2f ms' %
          (evoked.comment, pch, plat*1e3))

mne.viz.plot_evoked_topo(evokeds, merge_grads=True)


