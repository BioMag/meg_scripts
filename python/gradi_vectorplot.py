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
import matplotlib.pyplot as plt
import os.path as op
from scipy import signal
from mne.channels.layout import _merge_grad_data, _pair_grad_sensors


# order of lowpass IIR filter
BORD = 5
# colors for evoked sets
colors = [(166, 206, 227), (31, 120, 180), (178, 223, 138), (51, 160, 44),
          (251, 154, 153), (227, 26, 28), (253, 191, 111), (255, 127, 0),
          (202, 178, 214), (106, 61, 154)]
colors = [(x[0]/255., x[1]/255., x[2]/255.) for x in colors]

# parse command line
parser = argparse.ArgumentParser()
parser.add_argument('files', metavar='evoked_file',
                    help='Name of evoked fiff file', nargs='+')
parser.add_argument('--lowpass', type=float, metavar='f',
                    default=None, help='Lowpass frequency (Hz)')
parser.add_argument('--baseline', nargs=2, metavar=('t0', 't1'),
                    help='Start and end of baseline period (s)')
args = parser.parse_args()

# read evoked data
evokeds = []
filenames = []

for file in args.files:
    evokeds_ = mne.read_evokeds(file)
    evokeds.extend(evokeds_)
    filenames.extend([op.split(file)[1]]*len(evokeds_))

nev = len(evokeds)
if nev > len(colors):
    raise ValueError('Too many evoked sets')

# preprocess

if args.baseline:
    bl0 = float(args.baseline[0])
    bl1 = float(args.baseline[1])
    if bl0 >= bl1:
        raise ValueError('Baseline end should be larger than start')


for evoked, fn in zip(evokeds, filenames):
    
    evoked.comment = fn + ': ' + evoked.comment

    # filter
    if args.lowpass:
        data = evoked.data
        sfreq = evoked.info['sfreq']
        lpfreqn = 2 * np.array(args.lowpass) / sfreq
        b, a = signal.butter(BORD, lpfreqn)
        evoked.data = signal.filtfilt(b, a, data)

    # get peak amplitude - TODO: does not combine grads
    pch, plat = evoked.get_peak(ch_type='grad', merge_grads=True)
    print('%s peak amplitude: channel %s, latency %.2f ms' %
          (evoked.comment, pch, plat*1e3))

    # compute std on baseline
    if args.baseline:
        print('Baseline std. dev for each channel pair:')
        i0 = evoked.time_as_index(bl0)[0]
        i1 = evoked.time_as_index(bl1)[0]

        picks = _pair_grad_sensors(evoked.info, topomap_coords=False)
        data = evoked.data
        ch_names = evoked.ch_names

        data = data[picks]
        ch_names = [ch_names[k] for k in picks]

        data = _merge_grad_data(data)
        ch_names = [ch_name[:-1] + 'X' for ch_name in ch_names[::2]]

        bl_stds = data[:, i0:i1].std(axis=1)
        for std, ch in zip(bl_stds, ch_names):
            print('%s: %.4f fT/cm' % (ch, std*1e13))


colors_ = colors[:nev]

# we do our own legend with filenames
mne.viz.plot_evoked_topo(evokeds, color=colors_, merge_grads=True)


plt.show()

