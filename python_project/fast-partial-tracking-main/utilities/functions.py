"""
Documentation
============

This code provides functions usefull for both the state of the art analysis/synthesis method from Depalle and Neri and our own method developped for INTIM CNRS project

"""

import numpy as np
import scipy.io.wavfile as wav
from scipy.signal.windows import hann
import pandas as pd
import os 
import pickle

def time_vector_computation(n,Q) : 
    """
    Compute the time vector 

    Input : n is the discrete time instant
            Q is the polynomial order 

    Output : the time vector

    """
    return n ** np.arange(Q+1)
    

def derivative_time_vector_computation(n,Q):
    """
    Compute the derivative of the time vector 

    Input : n is the discrete time instant
            Q is the polynomial order 

    Output : the derivative time vector

    """
    return n ** np.arange(0,Q) * np.arange(1, Q+1)

def FreqAmpPha(fs, alpha, p, pPrime) :
    """
    Compute the log-amplitude, the frequency and the phase from alpha

    Input : fs is the sampling frequency
            alpha are the paramater estimated by ddm method
            p and pPrime are respectively time index and its derivative


    Output : Frequency, amplitude and phase vector calculated from alpha parameters

    """
    # Eq. (7)

    freq = np.array(np.expand_dims( np.imag(pPrime @ alpha[1:, :]),axis=1))

    # Eq. (6)

    Pha = np.array(np.expand_dims(np.imag(p @ alpha),axis=1))

    # Eq. (5)

    amp  = np.array(np.expand_dims(np.real(p @ alpha),axis=1))

    return np.array([freq, amp, Pha]).T[:]

def FreqAmpPhadAmp(fs, alpha, p, pPrime) :
    """
    Compute the log-amplitude, the frequency and the phase and the derivative of the log-amplitude from alpha

    Input : fs is the sampling frequency
            alpha are the paramater estimated by ddm method
            p and pPrime are respectively time index and its derivative


    Output : Frequency, amplitude, phase vector and derivative of the amplitude calculated from alpha parameters

    """
    # Eq. (7)

    freq = np.array(np.expand_dims( np.imag(pPrime @ alpha[1:, :]),axis=1))

    # Eq. (6)

    Pha = np.array(np.expand_dims(np.imag(p @ alpha),axis=1))

    # Eq. (5)

    amp  = np.array(np.expand_dims(np.real(p @ alpha),axis=1))

    # derivative of the log-amplitude 

    damp = np.array(np.expand_dims( np.real(pPrime @ alpha[1:, :]),axis=1))

    return np.array([freq, amp, Pha, damp]).T[:]

def match_sets(A, B):
        """
        Matches the elements of A and B.
        
        Parameters:
        A (list or array): First set of unique elements.
        B (list or array): Second set of unique elements.
        
        Returns:
        tuple: (match, nomatch_A, nomatch_B)
            match (numpy array): Pairs of matching indices.
            nomatch_A (numpy array): Indices of elements in A that did not match.
            nomatch_B (numpy array): Indices of elements in B that did not match.
        """

        nA = len(A)
        nB = len(B)
        
        # Ensure elements in A and B are unique
        if len(set(A)) != nA:
            print("A must be unique")
            return None, None, None
        elif len(set(B)) != nB:
            print("B must be unique")
            return None, None, None
        
        validA = np.ones(nA, dtype=bool)  # Boolean mask for valid elements in A
        validB = np.ones(nB, dtype=bool)  # Boolean mask for valid elements in B
        
        match = np.zeros((min(nA, nB), 2), dtype=int)  # Array to store matched indices
        counter = 0  # Counter for matched pairs
        
        for a in range(nA):
            if validA[a]:
                # Find indices in B that match A[a] and are still valid
                indices = np.where(B[validB] == A[a])[0]
                
                if indices.size > 0:
                    counter += 1
                    bs = np.where(validB)[0]  # Get valid indices in B
                    b = bs[indices[0]]  # Select the first match
                    
                    match[counter - 1, :] = [a, b]
                    validA[a] = False
                    validB[b] = False
        
        match = match[:counter, :]  # Trim to the actual number of matches
        nomatch_A = np.where(validA)[0]  # Indices of unmatched elements in A
        nomatch_B = np.where(validB)[0]  # Indices of unmatched elements in B
        
        return match, nomatch_A, nomatch_B

def bh_window(N, d=1):

    """
    Blackman-Harris Window
    
    Parameters:
    N : int
        Length of window
    d : int, optional
        Order of derivative (default is 1)
        
    Returns:
    win : ndarray
        The computed window
    n : ndarray
        The time indices
    """
    if N % 2:
        Nwinhf = (N - 1) // 2
        n = np.arange(-Nwinhf, Nwinhf + 1)
    else:
        Nwinhf = N // 2 - 1
        n = np.arange(-Nwinhf - 1, Nwinhf + 1)
    
    a = np.array([0.35875, 0.48829, 0.14128, 0.01168])
    in_vals = np.array([2, 4, 6]) * np.pi / N
    
    if d == 1:
        win = a[0] + a[1] * np.cos(in_vals[0] * n) + a[2] * np.cos(in_vals[1] * n) + a[3] * np.cos(in_vals[2] * n)
    elif d == 2:
        win = -a[1] * in_vals[0] * np.sin(in_vals[0] * n) \
              -a[2] * in_vals[1] * np.sin(in_vals[1] * n) \
              -a[3] * in_vals[2] * np.sin(in_vals[2] * n)
    elif d == 3:
        win = -a[1] * in_vals[0]**2 * np.cos(in_vals[0] * n) \
              -a[2] * in_vals[1]**2 * np.cos(in_vals[1] * n) \
              -a[3] * in_vals[2]**2 * np.cos(in_vals[2] * n)
    else:
        raise ValueError("Derivative order d must be 1, 2, or 3")
    
    return win

def hann_window(N, derivative=1):
    """
    Hann Window
    
    Parameters:
    N : int
        Length of window
    derivative : int, optional
        Order of derivative (default is 1)
        
    Returns:
    win : ndarray
        The computed window
    """
    if N % 2 == 1:
        Nhf = (N - 1) // 2
        n = np.arange(0, N) - Nhf
        in_val = 2 * np.pi / (N - 1)
        
        if derivative == 1:
            win = 0.5 + 0.5 * np.cos(in_val * n)
        elif derivative == 2:
            win = -0.5 * in_val * np.sin(in_val * n)
        elif derivative == 3:
            win = -0.5 * in_val**2 * np.cos(in_val * n)
        else:
            raise ValueError("Derivative order must be 1, 2, or 3")
    else:
        n = np.arange(0, N)
        in_val = 2 * np.pi / (N - 1)
        
        if derivative == 1:
            win = 0.5 - 0.5 * np.cos(in_val * n)
        elif derivative == 2:
            win = 0.5 * in_val * np.sin(in_val * n)
        elif derivative == 3:
            win = 0.5 * in_val**2 * np.cos(in_val * n)
        else:
            raise ValueError("Derivative order must be 1, 2, or 3")
    
    return win


def compute_frequencies_fourier():
    """
    Compute the vector of frequencies of Fourier transform for our specific implementation of Ndft = 4096 and fs = 48000

        
    Returns:
    freq : vector of frequencies 
    """

    Ndft = 4096
    fs = 48000
    frequencies = np.fft.rfftfreq(Ndft, d=1/fs)
    return frequencies


def find_f_true(trame_deltaPhi, frequencies):
    """
    Compute the vector of estimated frequencies using the Snail-Analyzer descriptor and our method designed during the project

    Input : deltaPhi descripor for the current frame, vector of frequencies
        
    Returns:
    true_f : vector of estimated peak frequencies 
    true_f_index : index aroud true_f in the vector for a frae 
    """

    N = len(trame_deltaPhi)
    true_f_ = []
    true_f_index = []
    for k in range(N - 1):
        if trame_deltaPhi[k] >= 0 and trame_deltaPhi[k + 1] < 0:
            true_f_index.append([k, k + 1])
            x = np.array([trame_deltaPhi[k], trame_deltaPhi[k + 1]])
            y = np.array([frequencies[k], frequencies[k + 1]])
            f_true = y[0] + (y[1] - y[0]) * (0 - x[0]) / (x[1] - x[0])
            true_f_.append(f_true)
    return true_f_, np.array(true_f_index)



def interp_amp(trame_amp, true_f_, true_f_index,frequencies):
    """
    Compute the vector of estimated amplitude using the Snail-Analyzer descriptor and our method designed during the project : linear interpolation

    Input : Amplitude descripor for the current frame, peak frequencies true_f, frequencies index around the peak, vector of frequencies
        
    Returns:
    true_amp : vector of estimated peak amplitude
    true_f_index : index aroud true_f in the vector for a frae 
    """

    amp_true_ = []
    for k in range(len(true_f_)):
        x = np.array([trame_amp[true_f_index[k, 0]], trame_amp[true_f_index[k, 1]]])
        y = np.array([frequencies[true_f_index[k, 0]], frequencies[true_f_index[k, 1]]])
        amp_true = np.interp(true_f_[k], y, x)
        amp_true_.append(amp_true)
    return amp_true_


def interp_damp(tram_amp,trame_amp_last_raw,true_f_,frequencies,true_f_index,H):
    """
    Compute the vector of estimated derivative amplitude using the Snail-Analyzer descriptor and our method designed during the project : finite difference, WARNING : the implemented method is not correct 

    Input : Amplitude descripor for the current frame and the previous frame, peak frequencies true_f, frequencies index around the peak, vector of frequencies, step size
        
    Returns:
    damp : vector of estimated peak derivative of amplitude

    """
    N = len(true_f_)
    damp = np.zeros(N)
    tram_diff = np.zeros(N)
    for k in range (N):
        x = np.array([trame_amp_last_raw[true_f_index[k,0]] ,trame_amp_last_raw[true_f_index[k,1]]])
        y = np.array([frequencies[true_f_index[k, 0]], frequencies[true_f_index[k, 1]]])
        tram_diff[k] =  np.interp(true_f_[k], y, x)
        damp[k] = (tram_amp[k] - tram_diff[k]) / H
    return np.array(damp)

def interp_phase(trame_phase,true_f_,true_f_index,frequencies):
    """
    Compute the vector of estimated phase using the Snail-Analyzer descriptor and our method designed during the project : basic linear interpolation without unwrapping 

    Input : Fourier phase descripor for the current frame, peak frequencies true_f, frequencies index around the peak, vector of frequencies, step size
        
    Returns:
    phase_true : vector of estimated peak derivative of amplitude

    """
    N = len(trame_phase)
    phase_true_ = [] 
    for k in range(len(true_f_)):
        x = np.array([trame_phase[true_f_index[k,0]] ,trame_phase[true_f_index[k,1]]])
        y = np.array([frequencies[k], frequencies[k+1]])
        phase_true = np.angle(np.exp(1j*x[0])+ (true_f_[k]-frequencies[true_f_index[k,0]])/(frequencies[true_f_index[k,1]]-frequencies[true_f_index[k,0]])*(np.exp(1j*x[1])-np.exp(1j*x[0])))
        phase_true_.append(phase_true)
    return np.array(phase_true_)

def interp_phase2(trame_phase, true_f_, true_f_index, frequencies):

    """
    Compute the vector of estimated phase using the Snail-Analyzer descriptor and our method designed during the project : linear interpolation with unwrapping 

    Input : Fourier phase descripor for the current frame, peak frequencies true_f, frequencies index around the peak, vector of frequencies
        
    Returns:
    phase_true : vector of estimated peak derivative of amplitude

    """
    phase_true_ = []
    for k in range(len(true_f_)):
        # Indices des fréquences encadrantes
        idx_low, idx_high = true_f_index[k]

        # Phases correspondantes
        phi0 = trame_phase[idx_low]
        phi1 = trame_phase[idx_high]

        # Dépliage de la phase pour éviter les discontinuités 
        phi1_unwrapped = phi0 + np.angle(np.exp(1j * (phi1 - phi0)))

        # Interpolation linéaire dans l’espace angulaire
        f0 = frequencies[idx_low]
        f1 = frequencies[idx_high]
        f_true = true_f_[k]
        w = (f_true - f0) / (f1 - f0) if f1 != f0 else 0.0

        phase_interp = (1 - w) * phi0 + w * phi1_unwrapped
        phase_true_.append(phase_interp)

    return np.array(phase_true_)

def amplitude_thresholding(threshold, amp_true, f_true, true_f_index):
    """
    Compute the vector of estimated ampltiude and frequncies using the Snail-Analyzer descriptor and our method designed during the project : linear interpolation with unwrapping 

    Input : treshold value in dB, peak frequencies true_f, frequencies index around the peak
        
    Returns:
    amp_treshold : vector of estimated amplitude with applied treshold
    f_true_treshold and true_f_index_treshold : associated frequencies and index around it for the threshold values of frequencies

    """
    f_true_threshold = []
    true_f_index_threshold = []
    amp_true_threshold = []
    for k in range(len(amp_true)):
        if amp_true[k] > threshold:
            amp_true_threshold.append(amp_true[k])
            f_true_threshold.append(f_true[k])
            true_f_index_threshold.append(true_f_index[k,:])

    return np.array(amp_true_threshold),np.array(f_true_threshold), np.array(true_f_index_threshold)

def window_amplitude_interpolation(trame_amp, f_true_th, f_true_index_th,frequencies,window_hat,window_freq):
    """
    Compute the vector of estimated amplitude using the Snail-Analyzer descriptor and our method designed during the project : window interpolation with 2 points

    Input : Amplitude descripor for the current frame and the previous frame, peak frequencies true_f, frequencies index around the peak, interpolated value of the fourier transform of the window, values of frequencies of this window
        
    Returns:
    damp : vector of estimated peak derivative of amplitude

    """

    N = len(f_true_th)
    amp_ = np.zeros(N)
    for k in range(N): 

        # point 1 
        A1_dB = trame_amp[f_true_index_th[k,0]]
        A1_lin = 10**(A1_dB/20)
        f1 = frequencies[f_true_index_th[k,0]]
        # point 2
        A2_dB = trame_amp[f_true_index_th[k,1]]
        A2_lin = 10**(A2_dB/20)
        f2 = frequencies[f_true_index_th[k,1]]

        delta_f1 = np.abs(f_true_th[k]-f1)
        delta_f2 = np.abs(f_true_th[k]-f2) 

       # W = np.interp(delta_final,window_freq,window_hat)
        W1 = np.interp(delta_f1,window_freq,window_hat) #cx
        W2 = np.interp(delta_f2,window_freq,window_hat) #cx
        #A_star = A_final/np.abs(W)
        A_star = 1/(np.abs(W1)**2+np.abs(W2)**2) * (np.abs(W1)*A1_lin+np.abs(W2)*A2_lin)

        amp_[k] = 20*np.log10(np.abs(A_star))

    return(amp_)


def compute_residual(original,synthesis):
    """
    Compute residual signal from two signals 

    Input : two signals of same size, warning is printed if not

    Returns : residual signal

    """

    if len(original)!=len(synthesis) : 
        print('Issue : size difference')
    else : 
        residual= original - synthesis
    return(residual)


def load_audio(file_path) : 
    """
    fonction to load audio from a wav file, convert to mono, normalize between [-1,1]

    """
    fs, signal = wav.read(file_path)

    # convert to mono if stereo
    if len(np.shape(signal))>1:
        signal = 1/2 * (signal[:,0]+signal[:,1])

    #Convert to float32 (normalize between -1 and 1)
    if np.max(signal)>1 :
        signal = signal.astype('float32') / 32768.0
    return fs, signal


def write_audio(signal,fs,config_name,method, save_dir) : 
    """
    fonction to write audio to a wav file, automated procedure of naming from config file

    """
    os.makedirs(save_dir,exist_ok=True)
    filename = f"{method}_partials_{config_name}.wav"
    path = os.path.join(save_dir,filename)
    wav.write(path, rate = fs , data = signal)
    print("Signal have been exported as .wav")


def load_snail_estimators(paths_dict):
    """
    function to load snail estimators from a CSV file using the config file

    """
    loudness = pd.read_csv(paths_dict["loudness"], header=None).to_numpy()
    deltaPhi = pd.read_csv(paths_dict["deltaPhi"], header=None).to_numpy()
    phiDem = pd.read_csv(paths_dict["phiDem"], header=None).to_numpy()
    return deltaPhi, loudness, phiDem

def save_partials(Partials,config_name,method, save_dir):
    """
    function to to save tracked partials to a .pkl file using the config file

    """
    os.makedirs(save_dir,exist_ok=True)
    filename = f"{method}_partials_{config_name}.pkl"
    path = os.path.join(save_dir,filename)
    with open(path, "wb") as f:
        pickle.dump(Partials, f)

    print(f"Partiels sauvegardés dans : {path}")

def save_spectro(spectro,config_name,method, save_dir):
    """
    function to to save spectrogramm to a .npy file file using the config file

    """
    os.makedirs(save_dir,exist_ok=True)
    filename = f"{method}_spectro_{config_name}.npy"
    path = os.path.join(save_dir,filename)
    np.save(path,spectro)

    print(f"spectrogramme sauvegardé dans : {path}")



def snr(x, y):
    """
    Compute the signal to noise ratio of two signals

    Parameters:
    x,y : input signals

    Returns:
    
    snr : signal to noise ration in dB
    """

    signal_power = np.sum(x ** 2)
    noise_power = np.sum((y - x) ** 2)
    snr_value = 10 * np.log10(signal_power / noise_power)
    
    return snr_value



def analyze_partials_from_dict(Partials):
    """
    Analyse le dictionnaire Partials pour extraire :
    - le nombre max de partiels actifs sur une frame
    - le nombre moyen de partiels actifs
    - le nombre de partiels actifs par frame

    Parameters:
    Partials : dict
        Dictionnaire indexé par frame, chaque valeur est un tableau (N, 5)

    Returns:
    stats : dict
        Statistiques sur l'activité des partiels
    """
    frames = sorted(Partials.keys())
    active_counts = []

    for m in frames:
        frame_data = Partials[m]
        if frame_data is not None and frame_data.shape[0] > 0:
            active_counts.append(frame_data.shape[0])
        else:
            active_counts.append(0)

    active_counts = np.array(active_counts)
    max_active = int(np.max(active_counts))
    mean_active = float(np.mean(active_counts))

    stats = {
        'max_active_partials': max_active,
        'mean_active_partials': mean_active,
        'active_counts_per_frame': active_counts,
        'frames': frames
    }

    return stats