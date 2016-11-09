#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 10:49:06 2016

@author: jussi
"""

from __future__ import print_function

import mne
import matplotlib.pyplot as plt
import numpy as np
import argparse
from mne.transforms import quat_to_rot, rotation_angles, _angle_between_quats
from scipy.signal import medfilt, filtfilt, butter
from scipy.linalg import inv


# lowpass filter
def lowpass(y, sfrate, corner, axis=-1):
    """ lowpass filter given data """
    ord = 5
    if corner is None:
        return y
    cornern = 2 * corner / sfrate
    b, a = butter(ord, cornern)
    yfilt = filtfilt(b, a, y, axis=axis)
    return yfilt

# parse command line
parser = argparse.ArgumentParser(description='Get head movement info from raw '
                                             'fiff file.')
parser.add_argument('fiff_file', help='Name of raw fiff file')

parser.add_argument('--plot', help='Plot HPI data', action='store_true')
parser.add_argument('--lpcorner', type=int, default=10,
                    help='Lowpass corner frequency for filtering position '
                         'and angle data')

args = parser.parse_args()

print(args.plot)

GOF_LIMIT = .98  # goodness of fit limit where chpi fit is still ok
DS_FREQ = 100  # frequency to downsample to 

raw = mne.io.Raw(args.fiff_file)
if 'CHPI001' not in raw.info['ch_names']:
    raise ValueError('Cannot find movement data in file. File needs to be '
                     'processed with MaxFilter first (movement estimation '
                     'or compensation).')
raw.load_data()
raw.pick_types(meg=False, chpi=True)
raw.resample(DS_FREQ)

picks_chpi_r = mne.pick_types(raw.info, meg=False, include=['CHPI001', 'CHPI002',
                                                          'CHPI003'])
picks_chpi_t = mne.pick_types(raw.info, meg=False, include=['CHPI004', 'CHPI005',
                                                          'CHPI006'])
picks_chpi_g = mne.pick_types(raw.info, meg=False, include=['CHPI007'])

# goodness of fit
g, times_all = raw[picks_chpi_g, :]
g = g.squeeze()
ind_ok = np.where(g > GOF_LIMIT)[0]
times_ok = times_all[ind_ok]

"""
quats give device -> head transformation: y = Rx + T 
where y is in head coords and x in device coords.
to get pos of head origin in device coords: y=0 -> Rx = -T -> x = -R.T * T 
(symmetric rot matrix)
result agrees with maxfilter (note that maxfilter plots coords of its sss
expansion origin, not of the head origin!)
rotation matrix R is negative of maxfilter (opposite rotations)
"""

# rotation quaternions
Rq, _ = raw[picks_chpi_r, :]
Rq = Rq[:, ind_ok].T

# convert to rotations around x, y and z axes
R = quat_to_rot(Rq)
nrot = R.shape[0]
A = np.zeros((nrot, 3))
for k in np.arange(nrot):
    A[k, :] = np.array(rotation_angles(R[k, :, :]))/np.pi * 180

# filter quaternion data and get difference angles directly
Rq = lowpass(Rq, raw.info['sfreq'], args.lpcorner, axis=0)
dA = _angle_between_quats(Rq[1:, :], Rq[:-1, :]) / np.pi * 180

# translation
T, _ = raw[picks_chpi_t, :]
T = T[:, ind_ok].T
T = lowpass(T, raw.info['sfreq'], args.lpcorner, axis=0)
dT = np.diff(T, axis=0)
dTv = np.sqrt(np.sum(dT**2, axis=1))  # magnitude of movement at each time point
tlen = times_ok[-1] - times_ok[0]

dTtot = np.sum(dTv)
dAtot = np.sum(dA)
print('\nTotal head translation during recording: %.2f mm, average %.2f mm/s' % 
      (dTtot*1e3, dTtot*1e3 / tlen))
print('Total head rotation during recording: %.2f deg, average %.2f deg/s\n' % 
      (dAtot, dAtot / tlen))

if args.plot:
    plt.figure()
    ax = plt.subplot(4, 1, 1)
    ax.plot(times_ok[1:], 1e3*dTv)
    ax.set_title('Head translations (where GOF > %g)' % GOF_LIMIT)
    ax.set(ylabel='Translation (mm)')

    ax = plt.subplot(4, 1, 2)
    ax.plot(times_ok[1:], dA)
    ax.set_title('Head rotations (where GOF > %g)' % GOF_LIMIT)
    ax.set(ylabel='Rotation (deg)')

    ax = plt.subplot(4, 1, 3)
    ax.plot(times_ok, T*1e3)
    ax.set_title('Position of head coordinate origin (whole file)')
    ax.legend(['x','y','z'])
    ax.set(ylabel='Position (mm)')

    ax = plt.subplot(4, 1, 4)
    ax.plot(times_all, g)
    ax.set_title('Goodness of fit (whole recording)')
    ax.set(ylabel='Translation (mm)', xlabel='Time (s)')
    ax.set_ylim([.8, 1.05])
    
    plt.tight_layout()
   
    plt.show()

