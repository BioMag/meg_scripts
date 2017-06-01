#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 15:52:33 2017

@author: jussi
"""

from __future__ import print_function
import mne
import os.path as op
import argparse
import sys

# parse command line
parser = argparse.ArgumentParser()
parser.add_argument('raw_fname', help='Name of raw fiff file')
parser.add_argument('--magrej', type=float, default=None,
                    help='Magnetometer rejection limit (T)')
parser.add_argument('--gradrej', type=float, default=None,
                    help='Gradiometer rejection limit (T/m)')
parser.add_argument('--magflat', type=float, default=None,
                    help='Magnetometer flatness limit (T)')
parser.add_argument('--gradflat', type=float, default=None,
                    help='Gradiometer flatness limit (T/m)')
parser.add_argument('--eogrej', type=float, default=None,
                    help='EOG rejection limit (V)')
args = parser.parse_args()


# raw_fname = u'/net/tera2/data/neuro-data/kekeke_kuu/Jussille_silmänliikkeet/5131/IE_taktilli_OB_Lhand_raw_st16_mc_b.fif'
# raw_fname = u'/net/tera2/data/neuro-data/kekeke_kuu/Jussille_silmänliikkeet/3282/OV_taktiiliOB_RHand_raw_st16_mc_bb.fif'


raw = mne.io.read_raw_fif(args.raw_fname, preload=True)

raw_fname_ = op.split(args.raw_fname)[1]
out_fname = op.splitext(raw_fname_)[0] + '_responses-ave.fif'

# mask high trigger lines to get rid of videoMEG trigger etc.
cnd = raw.acqparser.get_condition(raw, 'Pikkurilli', uint_cast=True,
                                  mask=2**8-1)
# set rejection criteria
reject = {'mag': args.magrej, 'grad': args.gradrej, 'eog': args.eogrej}
flat = {'mag': args.magflat, 'grad': args.gradflat}
# filter None (unset) values away
reject = {k: v for k, v in reject.items() if v is not None}
flat = {k: v for k, v in flat.items() if v is not None}
# gather epochs
eps = mne.Epochs(raw, reject=reject, flat=flat, preload=True, **cnd)
print('\n')

if len(eps) == 0:
    print('No epochs left after rejection!')
    sys.exit()

# pick EMG chs only
emg_picks = mne.pick_types(eps.info, meg=False, emg=True)
# start gui for rejecting epochs
# at this point, bad epochs have been dropped already
eps.plot(picks=emg_picks, block=True, title='Please select bad epochs')
# print stats
user_dropped = [k for k, reason in enumerate(eps.drop_log) if 'USER' in reason]
other_dropped = [k for k, reason in enumerate(eps.drop_log) if reason and
                 'USER'not in reason]

if user_dropped:
    print('%d epochs rejected by user:' % len(user_dropped), *user_dropped)
if other_dropped:
    print('%d epochs rejected by amplitude/flatness:' % len(other_dropped),
          *other_dropped)
print('%d epochs remain' % len(eps))
evs = eps.average()
print('Saving %s' % out_fname)
evs.save(out_fname)






            
            
            
        
        





