# -*- coding: utf-8 -*-

"""
SHow off elekta_averager

@author: jussi
"""


import numpy as np
import mne
from mne.datasets import sample
import sys
import matplotlib.pyplot as plt
from elekta_avg import Elekta_averager
#from chpi_weighted_average import chpi_snr_epochs, weighted_average_epochs



if sys.platform == 'win32':
    DATA_ROOT = 'C:/Users/HUS20664877/Dropbox/'
else:  # Linux
    DATA_ROOT = '/home/jussi/Dropbox/'

raw = mne.io.Raw('/home/jussi/projects/meg_scripts_git/test_data/aud_hilo_head_still_raw_tsss_mc.fif')
picks = mne.pick_types(raw.info, meg=True)

catname = 'auditory high'
eav = Elekta_averager(raw.info['acq_pars'])
eps = eav.get_epochs(raw, catname, picks,  mask=16384)

eva = eps.average()
eva.comment = catname
eva.save('aud_test_uusi.fif')


