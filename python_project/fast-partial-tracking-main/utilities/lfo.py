"""
Documentation
============

This code provides functions to generate LFO

"""


import numpy as np
import matplotlib.pyplot as plt

def lfo_sin(rate,amp,length,fs):
    t = np.arange(0,length/fs,1/fs)
    wave = amp*np.sin(2*np.pi*rate*t)
    return(wave)