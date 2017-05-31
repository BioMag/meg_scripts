#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 15:52:33 2017

@author: jussi
"""

from __future__ import print_function
import mne
import os.path as op


raw_fname = u'/net/tera2/data/neuro-data/kekeke_kuu/Jussille_silmänliikkeet/5131/IE_taktilli_OB_Lhand_raw_st16_mc_b.fif'
raw = mne.io.read_raw_fif(raw_fname, preload=True)

raw_fname_ = op.split(raw_fname)[1]
out_fname = op.splitext(raw_fname_)[0] + '_responses-ave.fif'

# get averaging category
ap = raw.acqparser
cnd = ap.get_condition(raw, 'Pikkurilli', uint_cast=True, mask=2**15,
                       mask_type='not_and')
# set rejection criteria
reject = {'eog': 200e-6, 'mag': 5e-12, 'grad': 5e-10}
flat = {'mag': 1e-15, 'grad': 1e-13}
# gather epochs
eps = mne.Epochs(raw, reject=reject, flat=flat, preload=True, **cnd)
# pick EMG chs only
emg_picks = mne.pick_types(eps.info, meg=False, emg=True)
# start gui for rejecting epochs
eps.plot(picks=emg_picks, block=True)
# print stats
user_dropped = [k for k, reason in enumerate(eps.drop_log) if 'USER' in reason]
print('\nUser rejected %d epochs:' % len(user_dropped), *user_dropped)
evs = eps.average()
evs.save(out_fname)






            
            
            
        
        





