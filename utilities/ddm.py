"""
Documentation
============

This code provides the ddm method for estimating the parameter of a given short term signal based on DDM method from M. Betser, "Sinusoidal Polynomial Parameter Estimation Using the Distribution Derivative", IEEE Transactions on Signal Processing, vol.
57, no. 12, pp. 4633-4645, Dec. 2009.

Inspired by the matlab code from : J. Neri and P. Depalle. "Fast partial tracking of audio with real-time capability through linear programming." Proceedings of the International Conference on Digital Audio Effects (DAFx). 2018.

"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# S
# ---------------------------------------------------------------------------


def pickPeaks(X, G_g):
    """
    Pick-peaking
    INPUTS:
    X : spectrum
    G_g : Peaks below this magnitude threshold (in dB) are not considered.
    OUTPUTS:
    peaks_inds :
    n_peaks :
    lefts_inds :
    rights_inds :
    """
    N, M = np.shape(X)
    idx_peaks = np.zeros([N, M])
    n_peaks = np.zeros([1, M])
    leftBin = np.zeros([N, M])
    rightBin = np.zeros([N, M])
    X = 20 * np.log10(np.abs(X))
    for m in range(M):  # loop through each frame ( however only one frame each time it is called ... )
        num = 0
        num_temp = 0
        for n in range(2, N-1):  # loop through each frequency bin
            if X[n, m] > G_g:  # check that peak are above the magnitude G_g
                if X[n, m] > X[n-1, m] and X[n, m] > X[n+1, m]:
                    num_temp += 1
                    left = n - 1
                    while (left > 0) and (X[left, m] > X[left - 1, m]):
                        left -= 1
                    right = n + 1
                    while (right < N-1) and (X[right, m] > X[right + 1, m]):
                        right += 1
                    num += 1
                    idx_peaks[num, m] = n
                    leftBin[num, m] = left
                    rightBin[num, m] = right
        n_peaks[m] = num  # number of peaks per one frame
    maxpeaks = int(n_peaks.item())
    peaks_inds = idx_peaks[1:maxpeaks+1, :]
    left_inds = leftBin[1:maxpeaks+1, :]
    right_inds = rightBin[1:maxpeaks+1, :]
    return peaks_inds, n_peaks, left_inds, right_inds


def FreqAmpPha(alpha,p,p_prime) :
    """
    Calculates the frequency, amplitude, and phase from alpha.
    
    Parameters:
    alpha : numpy array
        Input array containing alpha values.
    p : numpy array
        Array representing n^i for i=0:Q.
    p_prime : numpy array
        Derivative of p with respect to n.
    
    Returns:
    numpy array
        A concatenated array of frequency, amplitude, and phase values.
    """
    # Equation (7)
    frequency = np.imag(np.dot(p_prime, alpha[1:,:]))
    # Equation (5)
    amplitude = np.real(np.dot(p, alpha))
    # Equation (6)
    phase = np.imag(np.dot(p, alpha))
    
    return np.column_stack((frequency, amplitude, phase))



def ddm(y, Q, R, G_g, win, winD, p, pD, centering, Ndft, ft_mat, omega):
    """
    Distribution Derivative Method (DDM) for estimating sinusoidal model parameters
    
    Reference:
    M. Betser, "Sinusoidal Polynomial Parameter Estimation Using the Distribution Derivative",
    IEEE Transactions on Signal Processing, vol. 57, no. 12, pp. 4633-4645, Dec. 2009.
    
    Parameters:
    y: input signal
    Q: polynomial order
    R: number of peaks to use for estimation
    G_g: peak amplitude threshold (dB)
    win: window
    winD: derivative of window
    p: time vector n^i, where i = 0:Q
    pD: time derivative of p
    centering: centering vector for frequency domain operations
    Ndft: number of points for DFT
    ft_mat: DFT matrix for alpha0 estimation
    omega: frequency bin indices
    
    Returns:
    Alpha: estimated parameters
    num_peaks: number of detected peaks
    S: FFT of windowed signal
    """


    N = len(y)
    max_y = 2 * np.max(np.abs(y))
    
    # Windowed signal
    yst = win * y

    # DFT (one frame of the signal's STFT)
    S = np.array(np.expand_dims(np.fft.fft(yst, Ndft),axis = 1))
    
    # Pick Peaks (Placeholder for peak-picking function)
    peakBin, num_peaks, LeftBin, RightBin = pickPeaks(S[:Ndft // 2], G_g)

    Alpha = np.zeros((Q+1, int(num_peaks.item()))) if int(num_peaks.item()) > 0 else np.zeros((Q+1, 0))
    
    if num_peaks > 0:
        # Stores the various DFTs involved with the DDM method

        Sp = np.zeros((Ndft, Q), dtype=complex)
        Sp[:, 0] = np.array(S[:,0])
        Sp[:, 0] =  Sp[:, 0] * centering
        for i in range(1, Q):
            Sp[:, i] = np.fft.fft(yst * pD[:, i], Ndft) * centering
        
        # Derivative Windowed Signal 
        yDst = winD * y
        SD = np.fft.fft(yDst, Ndft) * centering + Sp[:, 0] * (-1j * omega)
        
        alpha_hat = np.zeros((Q+1, int(num_peaks.item())), dtype=complex)
        useful = 0
        
        for jj in range(int(num_peaks.item())):

            # Use the highest energy bins around the peak for DDM
            pbl = max(LeftBin[jj], peakBin[jj] - (R - 1))
            pbr = min([RightBin[jj], peakBin[jj] + (R - 1), Ndft // 2])
            pb_sides = np.sort(np.arange(pbl, pbr+1))[::-1]
            pbs = np.concatenate(([peakBin[jj].item()], pb_sides))
            
            # Define matrix A and vector b
            pbs = pbs.astype(int)
            A = Sp[pbs, :Q]
            b = -SD[pbs]
            
            # Solve for alpha using least squares
            alpha_temp = np.zeros(Q+1, dtype=complex)
            alpha_temp[1:] = np.linalg.lstsq(A, b, rcond=None)[0]
            
            # Check edges
            if np.any(np.isnan(alpha_temp)):
                print(f'NAN estimate, skipping peak {jj}')
            elif np.any(np.abs(np.dot(p[[0, N-1], 1:], alpha_temp[1:])) > 1e100):
                print(f'Infinite edges, skipping peak {jj}')
            else:
            
                # Solve for alpha0 (approximate) # see IV. AMPLITUDE AND PHASE ESTIMATION
                gam = win * np.exp(p[:, 1:] @ alpha_temp[1:])

                T_gam = np.dot(ft_mat[:N, int(peakBin[jj].item())].conj().T, gam)               
                alpha_temp[0] = np.log(Sp[int(peakBin[jj].item()), 0]) - np.log(T_gam)
                
                # Final check and save
                if np.all(np.dot(p[[0, N//2, N-1], :], alpha_temp.real) <= max_y):
                    alpha_hat[:, useful] = alpha_temp
                    useful += 1
        
        # Save the estimates
        num_peaks = useful
        Alpha = alpha_hat[:, :num_peaks]
    
    return Alpha, num_peaks, S


# ---------------------------------------------------------------------------
# tests
# --------------------------------------------------------------------------

# from scipy.io.wavfile import read
# import functions

# def test_ddm():
    
#     #fs = 44100  # Sampling rate
#     fs, signal = read("/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/sound/Sin8_A220Hz.wav")
    
#     N = 1024  # Frame size
#     t = np.arange(N) / fs
    
#     signal = signal[0:N,0]

#     win = np.hanning(N)
#     winD = np.gradient(win)
#     Ndft = 2048
#     ndft = np.arange(0,Ndft,1)
#     omega = 2 * np.pi * ndft/Ndft
    
#     # Fixing polynomial matrices to ensure correct shape
#     Q = 2
#     n_center = (N - (N % 2)) // 2
#     n = np.expand_dims(np.arange(N),axis=1) - n_center
#     p = functions.time_vector_computation(n,Q)
#     pD = functions.derivative_time_vector_computation(n,Q)

#     n_center = (N - (N % 2)) // 2
#     centering = np.exp(1j * omega * n_center)
#     ft_mat = np.exp(1j*np.outer(n , omega[:Ndft//2].T.conj()))
    
#     G_g = -60

#     Alpha, num_peaks, S = ddm(signal, Q, 4, G_g, win, winD, p, pD, centering, Ndft, ft_mat, omega)

#     freq = np.arange((Ndft//2))*fs/(Ndft//2)

#     # plt.figure()
#     # plt.plot(freq, np.abs(S[:Ndft//2]))
#     # plt.title("FFT Magnitude Spectrum")
#     # plt.show()
#     print(Alpha.shape)

#     results = np.zeros((len(Alpha),3))
#     amp = np.zeros(len(Alpha))
#     Pha = np.zeros(len(Alpha))

#     for k in range (len(Alpha)): 
#         results = FreqAmpPha(Alpha,p[n_center],pD[n_center])
    

# test_ddm()

