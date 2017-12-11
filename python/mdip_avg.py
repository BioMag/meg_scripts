#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read and analyze Elekta/Neuromag .mdip file (exported source waveform)

mdip lines are as follows (excl. comments):

    ndip nsamp
    t0 t1 dt

    (repeat following 3 lines for each dipole)
    1 (?)
    x y z Qx Qy Qz
    waveform (Q as func of time)

    GOF as func of time

@author: jnu@iki.fi
"""

from __future__ import print_function
import argparse
import numpy as np


def _datalines(fn):
    """Yield non-comment lines from text file fn"""
    with open(fn) as f:
        for line in f:
            if line[0] != '#':
                yield line


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get average source '
                                                 'amplitudes from xfit mdip '
                                                 'file')
    parser.add_argument('fn', help='Name of mdip file')
    parser.add_argument('winl', type=int, help='Window half length (ms)')
    parser.add_argument('--pos', help='Positive peaks only',
                        action='store_true')

    args = parser.parse_args()
    lines = [np.fromstring(l, sep=' ') for l in _datalines(args.fn)]

    print('Analyzing %s with averaging window of %d ms'
          % (args.fn, 2*args.winl+1))

    # read header
    ndip, nsamp = lines[0]
    t0, t1, dt = lines[1]
    t = np.arange(t0, t1, dt)
    print('%d dipole(s), %d samples, time axis %.1f..%.1f ms'
          % (ndip, nsamp, t0, t1))

    # read the waveform for each dipole
    for k in range(int(ndip)):
        dip = lines[3 + k * 3]
        wave = lines[4 + k * 3]
        waveabs = wave if args.pos else abs(wave)
        maxind = np.argmax(waveabs)
        maxt = t[maxind]
        win0, win1 = maxt - args.winl, maxt + args.winl
        ind0, ind1 = np.argmin(abs(t - win0)), np.argmin(abs(t - win1))
        avg = wave[ind0:ind1].mean()
        print('dipole %d: Qavg = %.2f nAm on period %.1f..%.1f ms'
              ' (max %.2f nAm at %.1f ms)'
              % (k+1, avg, win0, win1, wave[maxind], maxt))
