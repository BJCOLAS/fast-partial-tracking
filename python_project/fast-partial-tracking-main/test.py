# # ## My own STFT : 

# from scipy.io.wavfile import read
# import numpy as np
# from utilities import functions 
# import matplotlib.pyplot as plt
# from scipy.signal import stft
# from scipy.io import loadmat, savemat
# import pandas as pd
# # # Audio file input : 

# # input_filename = "/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/sound/sin8_A220Hz.wav"
# # fs, y = read(input_filename)

# # # Convert to float32 (normalize between -1 and 1)
# # y = y[1,:]
# # y = y.astype('float32') / 32768.0
# # # Normalisation 


# # # Analysis window length : 

# # N = 2 ** (int(np.ceil(np.log2(0.025 * fs))))-1 

# # # Oversampling factor : 

# # OverSample = 2 

# # # Hop size factor (HopSize = N / HopFactor )  ==> Indice de recouvrement 

# # HopFactor = 4

# # # Checking the input values : 

# # N = max(4,round(N))
# # HopFactor = 2 ** (int(np.ceil(np.log2(max(1,HopFactor)))))
# # OverSample = 2 ** (int(np.ceil(np.log2(max(1,OverSample)))))

# # L = len(y)

# # # Window selection 

# # if OverSample > 1 :
# #     win = functions.bh_window(N,1)
# #     HopFactor = max(8,HopFactor)
# # else : 
# #     win = functions.hann_window(N,1)
# #     HopFactor = max(4,HopFactor)

# # # Hopsize calculation 

# # H = round((N-N%2)/HopFactor)

# # # Size of the DFT using the oversampling factor : 

# # Ndft = 2 ** (round(np.log2(OverSample*2**round(np.log2(N)))))

# # # Check if Ndft is greater than N , else apply an OverSampling factor of 2 

# # if Ndft < N :
# #     Ndft = 2 ** round(np.log2(N)+1)

# # # Number of short frames calculation + padding !!! 

# # padding = np.zeros(2)
# # padding[0] = Ndft-H
# # M = int(np.ceil((Ndft + L -2*H -1)/H + 1))
# # padding[1] = ((M-1)*H+Ndft)-(padding[0]+L)

# # # Zero padding of the signal :
# # ypad = np.concatenate((np.zeros(int(padding[0])),y,np.zeros(int(padding[1]))))


# # # Time, frequency and Stft array

# # ndft = np.arange(Ndft)
# # time = np.zeros(M)
# # f = ndft / Ndft
# # S = np.zeros((Ndft,M),dtype = complex)

# # # Center of the analysis window ; 

# # n_center = int((N-N%2)/2 + 1) 

# # for m in range (M):

# #     # Time index at the center of the analysis window : 
# #     time[m] = m*H + n_center

# #     # Short term signal : 
# #     yst = ypad[int(time[m]-n_center): int(time[m]+N-n_center)]
# #     i = np.nonzero(yst) 
# #     y_win = win * yst

# #     S[:,m] = np.fft.fft(y_win,n = Ndft)

# #     # if m ==  202 :
# #     #     data = loadmat('/Users/colas/Documents/Programmation/spectro1.mat')
# #     #     A = data['S_9']  # Extract matrix
# #     #     plt.plot(A)
# #     #     plt.plot(np.fft.fft(y_win,n = Ndft),'--')
# #     #     #plt.plot(yst,'--')
# #     #     plt.show()

# # savemat("spectro.mat",{"S":S})

# # # Compute frequency axis
# # Ndft = S.shape[0]
# # f = fs * np.arange(Ndft) / Ndft
# # f = f[:Ndft // 2]

# # # Compute magnitude spectrogram and normalize
# # Smag = np.abs(S[:Ndft // 2, :])
# # Smag = Smag / Smag.max()

# # # Plot
# # plt.figure(figsize=(8, 5))
# # extent = [(time[0] - padding[0]) / fs, (time[-1] - padding[0]) / fs, f[0], f[-1]]
# # plt.imshow(20 * np.log10(Smag), aspect='auto', cmap='magma', extent=extent, origin='lower')

# # # Set color limits
# # plt.clim([-100, np.inf])

# # # Labels and formatting
# # plt.xlabel("Time (seconds)")
# # plt.ylabel("Frequency (Hz)")
# # plt.title("Spectrogram")
# # plt.colorbar(label="Magnitude (dB)")
# # plt.grid(True)
# # plt.axis([0, L / fs, 0, 5000])  # Adjusting frequency axis

# # plt.show()


# # df_loudness = pd.read_csv('/Users/colas/Documents/Programmation/Python/Snail_estimator/sin1/sine_wave_440_loudness.mat.csv')
# # loudness = df_loudness.to_numpy()
# # #logAmp = loudness*np.log(10)/20 # Conversion of the loudness to log-amplitude
# # df_deltaPhi = pd.read_csv('/Users/colas/Documents/Programmation/Python/Snail_estimator/sin1/sine_wave_440_deltaPhi.mat.csv')
# # deltaPhi = df_deltaPhi.to_numpy()
# # df_phiDem = pd.read_csv('/Users/colas/Documents/Programmation/Python/Snail_estimator/sin1/sine_wave_440_phiDem.mat.csv')
# # phiDem = df_phiDem.to_numpy()
# # Variables of the snail analyseur for sine 1

# # audio_frame = 48000
# # window_size = 50 #ms
# # step_size = 1 #ms
# # fs = 48000 # Hz
# # N = int(fs * window_size / 1000) # window size in samples
# # HopSize = int(fs * step_size / 1000) # step size in samples


# df_loudness = pd.read_csv('/Users/colas/Documents/Programmation/Python/Snail3/Sin8_A220Hz_loudness.mat.csv')
# loudness = df_loudness.to_numpy()    
# # Normalization of the loudness to be between -96 and 0 dB
# df_deltaPhi = pd.read_csv('/Users/colas/Documents/Programmation/Python/Snail3/Sin8_A220Hz_deltaPhi.mat.csv')
# deltaPhi = df_deltaPhi.to_numpy()   
# df_phiDem = pd.read_csv('/Users/colas/Documents/Programmation/Python/Snail3/Sin8_A220Hz_phiDem.mat.csv')
# phiDem = df_phiDem.to_numpy()  

# audio_frame = 144000
# fs = 48000 # Hz
# win_size = 50 #ms
# step_size = 5 #ms
# N = int(fs * win_size / 1000) # window size in samples
# HopSize = int(fs * step_size / 1000) # step size in samples
# HopFactor = int(N / HopSize) # HopFactor

# trame_deltaPhi = deltaPhi[5,:]
# trame_Loudness = loudness[5,:]
# trame_PhiDem = phiDem[5,:]
# # period
# T = 1/48000
# # peak threshold
# G_g = -55 #dB
# # Frequencies computation
# frequencies = functions.compute_frequencies()

# # f_true and f_true_index computation for one trame
# true_f_, true_f_index= functions.find_f_true(trame_deltaPhi,frequencies)
# # Amplitude and derivative of amplitude computation for one trame
# amp_true = functions.interp_amp(trame_Loudness,true_f_,true_f_index,frequencies)

# amp_true, true_f_, true_f_index = functions.amplitude_thresholding(G_g, amp_true, true_f_, true_f_index)
# log_amp_true = 1/2*amp_true*np.log(10)/20
# damp, damp_true = functions.interp_damp(trame_Loudness,T,true_f_,true_f_index,frequencies)
# # Phase computation for one trame
# phase_true = functions.interp_phase(trame_PhiDem,true_f_,true_f_index,frequencies)
# # Number of peaks estimated
# num_peaks = len(true_f_)



import numpy as np
import matplotlib.pyplot as plt

# Paramètres
fs = 48000                   # Fréquence d'échantillonnage
f0 = 480                    # Fréquence du signal
duration = 1.0              # Durée en secondes
dft_size = 4096
win_size = 2047
hop_size = 256

# Signal sinusoïdal
t = np.arange(0, duration, 1/fs)
x = np.sin(2 * np.pi * f0 * t)

# Fenêtre de Hanning
window = np.hanning(win_size)

# Nombre de trames
num_frames = (len(x) - win_size) // hop_size + 1

# STFT manuelle avec np.fft.fft
phases = []
bin_freqs = np.fft.fftfreq(dft_size, d=1/fs)
bin_index = np.argmin(np.abs(bin_freqs - f0))  # index du bin le plus proche de 480 Hz

for i in range(num_frames):
    start = i * hop_size
    frame = x[start:start + win_size] * window
    frame_padded = np.zeros(dft_size)
    frame_padded[:win_size] = frame
    spectrum = np.fft.fft(frame_padded)
    phase = np.angle(spectrum[bin_index])
    phases.append(phase)

# Unwrap pour suivi de phase lisse
#phases = np.unwrap(phases)

# Axe temporel
time_axis = np.arange(num_frames) * hop_size / fs

# Affichage
plt.figure(figsize=(10, 4))
plt.plot(time_axis, phases, label="Phase à 480 Hz")
plt.xlabel("Temps (s)")
plt.ylabel("Phase (rad)")
plt.title("Évolution de la phase à 480 Hz")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


    
