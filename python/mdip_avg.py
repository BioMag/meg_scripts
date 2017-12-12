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


def _nearest(A, x):
    """Find index of value in A nearest to x"""
    return np.argmin(abs(A - x))


if __name__ == '__main__':

    desc = """
    Get peak amplitudes from an .mdip file produced by XFit
    (source modelling). The .mdip contains source waveforms for time varying
    dipoles.
    The amplitudes are determined by averaging around the peak using a given
    half window length.
    """

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('fn', help='Name of mdip file')
    parser.add_argument('winl', type=int, help='Window half length (ms)')
    parser.add_argument('--search', type=float, nargs=2, metavar=('t0', 't1'),
                        help='Restrict peak search to interval t0..t1 (ms)')
    parser.add_argument('--pos', help='Find positive peaks only',
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

    if args.search:
        tstart, tend = args.search
        if tstart < t0 or tend > t1:
            raise ValueError('Time range must be inside dipole time course')
        indstart, indend = _nearest(t, tstart), _nearest(t, tend)
        t = t[indstart:indend+1]

    # read the waveform for each dipole
    for k in range(int(ndip)):
        dip = lines[3 + k * 3]
        wave = lines[4 + k * 3]
        wave = wave[indstart:indend] if args.search else wave
        waveabs = wave if args.pos else abs(wave)
        maxind = np.argmax(waveabs)
        maxt = t[maxind]
        win0, win1 = maxt - args.winl, maxt + args.winl
        ind0, ind1 = _nearest(t, win0), _nearest(t, win1)
        avg = wave[ind0:ind1].mean()
        print('dipole %d: Qavg = %.2f nAm at %.1f..%.1f ms'
              ' (max = %.2f nAm at %.1f ms)'
              % (k+1, avg, t[ind0], t[ind1], wave[maxind], maxt))
