#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot waveforms from Elekta/Neuromag .mdip files (exported source waveform)

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
import matplotlib.pyplot as plt


def _datalines(fn):
    """Yield non-comment lines from text file fn"""
    with open(fn) as f:
        for line in f:
            if line[0] != '#':
                yield line

def header(v):
    """Return header info from array v"""
    ndip, nsamp = v[0]
    t0, t1, dt = v[1]
    t = np.arange(t0, t1, dt)
    return ndip, t


def waveforms(v):
    """Return all dipole waveforms from array v"""
    ndip, nsamp = v[0]
    for k in range(int(ndip)):
       yield v[4 + k * 3]
        #dip = lines[3 + k * 3]


if __name__ == '__main__':

    desc = """
    Superpose & average mdip waveforms
    """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('fns', nargs='+', help='Names of mdip files',
                        metavar='filename')
    parser.add_argument('--dipn', type=int, default=1,
                        help='index of dipole (0 for the last one)')
    args = parser.parse_args()
    dipn = args.dipn - 1  # 1-based index, 0 is last (becomes -1)
    waves_all = list()
    
    for fn in args.fns:
        print(fn)
        v = [np.fromstring(l, sep=' ') for l in _datalines(fn)]    
        ndip, t = header(v)
        waves = list(waveforms(v))
        if args.dipn > ndip:
            raise ValueError('Dipole #%d does not exist in %s' % (args.dipn, fn))
        wave = waves[dipn]
        plt.plot(t, wave)
        waves_all.append(wave)

    plt.legend(args.fns)
    plt.title('Superposed waveforms (dipole %d)' % dipn)
    plt.xlabel('Time (ms)')
    plt.ylabel('Dipole amplitude (nAm)')
    wall = np.array(waves_all)
    wall_avg = wall.mean(axis=0)
    plt.figure()
    plt.plot(wall_avg)
    plt.title('Averaged waveform (dipole %d)' % dipn)
    plt.xlabel('Time (ms)')
    plt.ylabel('Dipole amplitude (nAm)')














