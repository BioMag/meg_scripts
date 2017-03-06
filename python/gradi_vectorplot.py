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


import mne
import argparse
import mne.viz
import numpy as np
from scipy import signal
from mne.channels.layout import _merge_grad_data, _pair_grad_sensors


# parse command line
parser = argparse.ArgumentParser()
parser.add_argument('evoked_file', help='Name of evoked fiff file')
parser.add_argument('category', help='Name of evoked category to plot')
parser.add_argument('--lowpass', type=float, metavar='f', default=None, help='Lowpass frequency (Hz)')
args = parser.parse_args()

# read evoked data
evoked = mne.read_evokeds(args.evoked_file, condition=args.category)

# filter and replace data in-place
BUTTORD = 5  # order of Butterworth IIR lowpass filter
if args.lowpass:
    data = evoked.data
    sfreq = evoked.info['sfreq']
    lpfreqn = 2 * np.array(args.lowpass) / sfreq
    b, a = signal.butter(BUTTORD, lpfreqn)
    evoked.data = signal.filtfilt(b, a, data)

# compute gradiometer RMS data
picks_gradc = _pair_grad_sensors(evoked.info, topomap_coords=False)
gradc_data = _merge_grad_data(evoked.data[picks_gradc])
# create list of combined channel names (e.g. MEG 204x)
ch_names_gradc = list()
for ind in picks_gradc:
    ch = evoked.ch_names[ind]
    if ch[-1] == '3':
        ch_names_gradc.append(ch[:-1]+'x')
    
# create evoked set with magnetometer channels only
evoked_mag = evoked.pick_types(meg='mag')

# change channel names of evoked set to represent gradiometer pairs (e.g. MEG 2041 -> MEG 204x)
for i,nm in enumerate(evoked_mag.ch_names):
    evoked_mag.ch_names[i] = nm[:-1]+'x'
for i,ch in enumerate(evoked_mag.info['chs']):
    evoked_mag.info['chs'][i]['ch_name'] = evoked_mag.info['chs'][i]['ch_name'][:-1]+'x'
    evoked_mag.info['chs'][i]['unit'] = 201  # change unit to T/m
   
# get magnetometer layout for platting, replace magnetometer names as above
laym = mne.channels.read_layout('Vectorview-mag')
for i,nm in enumerate(laym.names):
    laym.names[i] = nm[:-1]+'x'
    
# replace magnetometer data array with combined gradiometer data
evoked_mag.data = np.zeros(evoked_mag.data.shape)
for j,ch in enumerate(ch_names_gradc):
    for i,ch1 in enumerate(evoked_mag.ch_names):
        if ch1 == ch:
            #print(ch1,'matches',ch,'at index',i)
            evoked_mag.data[i,:] = gradc_data[j,:]

# get peak
pch, plat = evoked_mag.get_peak()
print('Peak amplitude: channel pair',pch,'at latency',plat*1e3,'ms')

mne.viz.plot_evoked_topo(evoked_mag, layout=laym, title=args.evoked_file)


