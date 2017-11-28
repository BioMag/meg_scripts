#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read Elekta/Neuromag .bdip file and write the dipole data into Excel file.

@author: jnu@iki.fi
"""

import os.path as op
import argparse
import numpy as np
from openpyxl import Workbook
import struct


def write_workbook(data, filename, title='Data', first_col=1, first_row=1):
    """Write data into .xlsx file (filename). data must be a iterable of
    iterables which represents the rows to write. first_col and
    first_row specify the column and row where to start writing (1-based).
    To write data as columns, pass zip(*data) as arg instead of data.
    """
    wb = Workbook()
    ws1 = wb.active
    ws1.title = title
    for row_ind, row in enumerate(data):
        for col_ind, val in enumerate(row):
            ws1.cell(column=col_ind+first_col+1, row=row_ind+first_row+1,
                     value=val)
    wb.save(filename=filename)


def decode_bdip(fn):
    """Read bdip as 196 byte chunks and decode each chunk. See XFit manual,
    Appendix B for more details. """
    BLOCKSIZE = 196
    fmt = '>i' + 12*'f' + 'i' + 35*'f'  # sequence of IEEE floats and ints
    with open(fn) as f:
        b = f.read(BLOCKSIZE)
        while b:
            yield struct.unpack(fmt, b)
            b = f.read(BLOCKSIZE)


def _bdip_vars(add_nice_units=True):
    """bdip variable names, in correct order"""
    vars = ['dipole', 'begin', 'end', 'r0x', 'r0y', 'r0z', 'x', 'y', 'z',
            'Qx', 'Qy', 'Qz', 'goodness', 'errors_computed', 'noise_level']
    # single parameter confidence limits
    vars += ['single_errors_%d' % k for k in range(5)]
    # error covariance matrix
    vars += ['error_matrix_%d%d' % (i, j) for i in range(5) for j in range(5)]
    # stats
    vars += ['conf_vol', 'khi2', 'prob', 'noise_est']
    # add (non-SI) units, see _nice_units()
    if add_nice_units:
        for ind in (1, 2):
            vars[ind] += ' (ms)'
        for ind in range(3, 9):
            vars[ind] += ' (mm)'
        for ind in range(9, 12):
            vars[ind] += ' (nAm)'
    return vars


def _nice_units(vals):
    """Convert vals returned by _decode_bdip() into nicer (non-SI) units"""
    vals_ = np.array(vals)
    vals_[1:3] *= 1e3  # s -> ms
    vals_[3:9] *= 1e3  # m -> mm
    vals_[9:12] *= 1e9  # Am -> nAm
    return vals_


def _filter_vars(vars):
    """Return indices of desired variables"""
    _exclude = 'r0', 'error_matrix', 'single_errors'
    inds_vars = [(i, var) for i, var in enumerate(vars)
                 if not any([s in var for s in _exclude])]
    return map(list, zip(*inds_vars))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('fname', help='Name of bdip file')
    args = parser.parse_args()
    fname_out = op.splitext(args.fname)[0] + '.xlsx'

    inds_ok, vars = _filter_vars(_bdip_vars())
    data = [[''] + vars]  # add empty cell to align with other rows

    data += [[args.fname] + list(_nice_units(data_dip)[inds_ok])
             for data_dip in decode_bdip(args.fname)]

    write_workbook(data, fname_out)
    print 'saved %s' % fname_out
