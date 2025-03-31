# ## My own STFT : 

from scipy.io.wavfile import read
import numpy as np
from utilities import functions 
import matplotlib.pyplot as plt
from scipy.signal import stft
from scipy.io import loadmat, savemat

# Audio file input : 

input_filename = "/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/sound/sin8_A220Hz.wav"
fs, y = read(input_filename)

# Convert to float32 (normalize between -1 and 1)
y = y[1,:]
y = y.astype('float32') / 32768.0
# Normalisation 


# Analysis window length : 

N = 2 ** (int(np.ceil(np.log2(0.025 * fs))))-1 

# Oversampling factor : 

OverSample = 2 

# Hop size factor (HopSize = N / HopFactor )  ==> Indice de recouvrement 

HopFactor = 4

# Checking the input values : 

N = max(4,round(N))
HopFactor = 2 ** (int(np.ceil(np.log2(max(1,HopFactor)))))
OverSample = 2 ** (int(np.ceil(np.log2(max(1,OverSample)))))

L = len(y)

# Window selection 

if OverSample > 1 :
    win = functions.bh_window(N,1)
    HopFactor = max(8,HopFactor)
else : 
    win = functions.hann_window(N,1)
    HopFactor = max(4,HopFactor)

# Hopsize calculation 

H = round((N-N%2)/HopFactor)

# Size of the DFT using the oversampling factor : 

Ndft = 2 ** (round(np.log2(OverSample*2**round(np.log2(N)))))

# Check if Ndft is greater than N , else apply an OverSampling factor of 2 

if Ndft < N :
    Ndft = 2 ** round(np.log2(N)+1)

# Number of short frames calculation + padding !!! 

padding = np.zeros(2)
padding[0] = Ndft-H
M = int(np.ceil((Ndft + L -2*H -1)/H + 1))
padding[1] = ((M-1)*H+Ndft)-(padding[0]+L)

# Zero padding of the signal :
ypad = np.concatenate((np.zeros(int(padding[0])),y,np.zeros(int(padding[1]))))


# Time, frequency and Stft array

ndft = np.arange(Ndft)
time = np.zeros(M)
f = ndft / Ndft
S = np.zeros((Ndft,M),dtype = complex)

# Center of the analysis window ; 

n_center = int((N-N%2)/2 + 1) 

for m in range (M):

    # Time index at the center of the analysis window : 
    time[m] = m*H + n_center

    # Short term signal : 
    yst = ypad[int(time[m]-n_center): int(time[m]+N-n_center)]
    i = np.nonzero(yst) 
    y_win = win * yst

    S[:,m] = np.fft.fft(y_win,n = Ndft)

    # if m ==  202 :
    #     data = loadmat('/Users/colas/Documents/Programmation/spectro1.mat')
    #     A = data['S_9']  # Extract matrix
    #     plt.plot(A)
    #     plt.plot(np.fft.fft(y_win,n = Ndft),'--')
    #     #plt.plot(yst,'--')
    #     plt.show()

savemat("spectro.mat",{"S":S})

# Compute frequency axis
Ndft = S.shape[0]
f = fs * np.arange(Ndft) / Ndft
f = f[:Ndft // 2]

# Compute magnitude spectrogram and normalize
Smag = np.abs(S[:Ndft // 2, :])
Smag = Smag / Smag.max()

# Plot
plt.figure(figsize=(8, 5))
extent = [(time[0] - padding[0]) / fs, (time[-1] - padding[0]) / fs, f[0], f[-1]]
plt.imshow(20 * np.log10(Smag), aspect='auto', cmap='magma', extent=extent, origin='lower')

# Set color limits
plt.clim([-100, np.inf])

# Labels and formatting
plt.xlabel("Time (seconds)")
plt.ylabel("Frequency (Hz)")
plt.title("Spectrogram")
plt.colorbar(label="Magnitude (dB)")
plt.grid(True)
plt.axis([0, L / fs, 0, 5000])  # Adjusting frequency axis

plt.show()
