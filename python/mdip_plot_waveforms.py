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
    return ndip, nsamp, t0, t1, dt


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
    parser.add_argument('--autoscale', help='autoscale average plot',
                        action='store_true')
    parser.add_argument('--legend', help='create a legend',
                        action='store_true')
    parser.add_argument('--average', help='plot the average of all waveforms',
                        action='store_true')
    args = parser.parse_args()
    dipn = args.dipn - 1  # 1-based index, 0 is last (becomes -1)
    waves_all = list()

    for fn in args.fns:
        print(fn)
        v = [np.fromstring(l, sep=' ') for l in _datalines(fn)]
        ndip, nsamp, t0, t1, dt = header(v)
        t = np.linspace(t0, t1, int(nsamp))
        waves = list(waveforms(v))
        if args.dipn > ndip:
            raise ValueError('Dipole #%d does not exist in %s' %
                             (args.dipn, fn))
        wave = waves[dipn]
        plt.plot(t, wave)
        waves_all.append(wave)

    if args.legend:
        plt.legend(args.fns)

    plt.title('Superposed waveforms (dipole %d)' % args.dipn)
    plt.xlabel('Time (ms)')
    plt.ylabel('Dipole amplitude (nAm)')
    ylim = plt.gca().get_ylim()

    if args.average:
        plt.figure()
        wall = np.array(waves_all)
        wall_avg = wall.mean(axis=0)
        plt.plot(t, wall_avg)
        plt.title('Averaged waveform (dipole %d)' % args.dipn)
        plt.xlabel('Time (ms)')
        plt.ylabel('Dipole amplitude (nAm)')
        if not args.autoscale:
            plt.gca().set_ylim(ylim[0], ylim[1])













