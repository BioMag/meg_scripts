#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read and analyze Elekta/Neuromag .mdip file (exported source waveform)

lines:
    ndip nsamp
    t0 t1 dt

    (following 3 lines for each dipole)
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
        s = f.readline()
        while s:
            if s[0] != '#':
                yield s
            s = f.readline()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get average source '
                                                 'amplitudes from xfit mdip '
                                                 'file')
    parser.add_argument('fn', help='Name of mdip file')
    parser.add_argument('winl', type=int, help='Window half length')
    args = parser.parse_args()
    lines = [np.fromstring(l, sep=' ') for l in _datalines(args.fn)]

    print('Analyzing %s with averaging window of %d samples'
          % (args.fn, args.winl))

    # read header
    ndip = int(lines[0][0])
    nsamp = int(lines[0][1])
    t0 = lines[1][0]
    t1 = lines[1][1]
    dt = lines[1][2]
    t = np.arange(t0, t1, dt)
    print('%d dipole(s), %d samples, time axis %.1f..%.1f ms'
          % (ndip, nsamp, t0, t1))

    # read the waveform for each dipole
    for k in range(ndip):
        dip = lines[3 + k * 3]
        wave = lines[4 + k * 3]
        maxind = np.argmax(wave)
        maxt = t[maxind]
        avg = wave[maxind-args.winl:maxind+args.winl].mean()
        print('dipole %d: Qavg = %.2f nAm around the maximum'
              ' (%.2f nAm at %.2f ms)' % (k+1, avg, wave.max(), maxt))
