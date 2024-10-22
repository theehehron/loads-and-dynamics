#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 18:15:36 2024

@author: aarongreisen
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, fft
import dynamics.vibration as vibe
# import dynamics.vibration as dy

acceptance_level = np.array([[  20, 0.0053],
                            [ 150, 0.04],
                            [ 600, 0.04],
                            [2000, 0.0036]])

# 1. Generate low passed white noise time history.
fs = 4000 # (Hz)
T = 30 # (sec)
N = int(fs*T)

mean = 0
std_dev = 19
white_noise = np.random.normal(mean, std_dev, N)

lowpass_filt_2000Hz = signal.butter(10, fs/2-1, 'lp', fs=fs, output='sos')
white_noise_filtered = signal.sosfilt(lowpass_filt_2000Hz, white_noise)
time = np.arange(N)/fs


# 2. Compute the fft of the white noise time history.
white_fft = fft.fft(white_noise_filtered)
fft_freqs = fft.fftfreq(N, 1/fs)
PSD_target, PSD_target_freqs = vibe.spectrum_points(acceptance_level, np.abs(fft_freqs))


# 3. scale the white noise signal in the frequency domain to match the target PSD
scale_fft = np.sqrt(PSD_target*fs/(2*np.power(std_dev, 2)))
synth_fft = scale_fft*white_fft


# 4. Compute the inverse fft to get the time history
synth_time_history = np.real(fft.ifft(synth_fft))


# 5. Time history analysis: evaluate the PSD, kurtosis, & histogram
f_synth, psd_synth = signal.welch(synth_time_history, fs)
kurtosis = np.sum(np.power(synth_time_history-np.mean(synth_time_history), 4))/(synth_time_history.shape[0]*np.power(np.std(synth_time_history), 4))

fig0, ax0 = plt.subplots()
ax0.loglog(acceptance_level[:,0], acceptance_level[:, 1])
ax0.loglog(f_synth, psd_synth)
ax0.set_ylim(0.001, 0.1)
ax0.set_xlim(np.min(acceptance_level[:,0]), np.max(acceptance_level[:,0]))
ax0.grid()
ax0.grid(which='minor', linestyle=':', linewidth='0.5', color='gray') 
plt.xlabel("Frequency (Hz)")
plt.ylabel("Acceleration PSD ($g^2$/Hz)")
plt.title('PSD of Synthesized Random Vibration Signal')

fig1, ax1 = plt.subplots()
ax1.hist(synth_time_history, bins=40)

