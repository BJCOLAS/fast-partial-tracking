"""
Documentation
============

This class provides the tracking of the partial coming from an audio signal. 

Inspired by the matlab code from : J. Neri and P. Depalle. "Fast partial tracking of audio with real-time capability through linear programming." Proceedings of the International Conference on Digital Audio Effects (DAFx). 2018.

"""



# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------


import numpy as np
from utilities import ddm, functions, munkres 
from pprint import pprint
from scipy.optimize import linear_sum_assignment
from scipy.signal.windows import blackmanharris
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Classes
# ---------------------------------------------------------------------------


class trackingPartials:

    """
    Implement a tracking method of partials : 

        Involves the following two main procedures for each short-term frame
    
     1. Short-term sinusoidal model parameters are estimated using the 
        Distribution Derivative Method (DDM).
     2. Peaks are connected over consecutive analysis frames by solving an
        assignment problem, using the Hungarian algorithm (munkres.m).
    
    
     INPUTS
     y: input audio signal
     fs: sampling rate (Hz)
     Variable inputs:
     N: length of short-term analysis window
     HopFactor: Hop Size = N/HopFactor
     OverSample: Oversampling factor
     G_g: Peaks below this magnitude threshold (in dB) are not considered.
     Q: short-term model polynomial order (DDM)
     delta: controls the relative preference towards useful assignments.
            larger values of delta promote useful assignments (0<delta<1)
     zetaF, zetaA: the approximate cross-over values from a useful
            to spurious assignment. 
    
     OUTPUTS
     Partials: a cell array, one cell for each frame containing the
     instantaneous estimates of frequency, amplitude, and phase, and index of
     the partials active in that frame.
     Time: the locations in time of the short-term frames
     Padding: the amount of zeros padding the signal
     L: the signal length
     S: the STFT (for plotting later on).
    """

    def __init__(self, y , fs , N = 2^11-1 , hopFactor=4 , overSample=1 , G_g=-40 , Q=2 , delta=0.2 , zetaF = 50 , zetaA = 15):
        """Create a new instance"""

        self.y = y 
        self.fs = fs
        self.N = max(4,round(N))
        self.hopFactor = max(2, int(2**np.ceil(np.log2(hopFactor))))
        self.overSample = max(2, int(2**np.ceil(np.log2(overSample))))
        self.G_g = np.abs(G_g)
        self.Q = max(1,round(Q))
        self.zetaF = max(1e-1, zetaF)
        self.zetaA = max(1e-1, zetaA)/20 * np.log(10)
        self.delta = delta 

        # Initialization 
        self.partials = []
        self.time = []
        self.padding = None
        self.M = 0
        self.omega = 0
        self.ndft = []
        self.Alpha = []

    def setup(self) : 
        """
        Setup method of the main variables use of the parameter estimation
        """
        # converting to mono
        self.y = self.y[:,1] if np.ndim(self.y)>1 else self.y

        # short term signal length
        self.L = len(self.y)

        # Compute variance for tracking
        self.varF = -self.zetaF**2 * np.log((self.delta - 1) / (self.delta - 2))
        self.varA = -self.zetaA**2 * np.log((self.delta - 1) / (self.delta - 2))

        # Window selection

        if self.overSample > 1 :
            self.win = functions.bh_window(self.N)
            self.winD = functions.bh_window(self.N,d=2)
            self.hopFactor = max(8,self.hopFactor)
        else :
            self.win = functions.hann_window(self.N)
            self.winD = functions.hann_window(self.N,d=2)
            self.hopFactor = max(4,self.hopFactor)
        
        # Hop size & DFT size
        self.H = round((self.N - self.N % 2) / self.hopFactor)
        self.Ndft = 2 ** round(np.log2(self.overSample * 2**round(np.log2(self.N))))
        if self.Ndft < self.N:
            self.Ndft = 2**round(np.log2(self.N) + 1)
        
        self.ndft = np.arange(0,self.Ndft,step=1).T
        
        # Compute zero-padding
        self.padding = [self.Ndft - self.H, 0]
        self.M = int(np.ceil((self.Ndft + self.L - 2 * self.H - 1) / self.H) + 1)
        self.padding[1] = ((self.M - 1) * self.H + self.Ndft) - (self.padding[0] + self.L)
        self.ypad = np.concatenate((np.zeros(self.padding[0]), self.y, np.zeros(self.padding[1])))

        # Time vector for each frame
        self.time = np.zeros(self.M)
        self.omega = 2*np.pi*self.ndft*1/self.Ndft
        self.S = np.zeros((self.Ndft, self.M), dtype=complex)

        # Threshold : 

        max_sample = np.argmax(np.abs(self.ypad))  # Find index of max absolute value in ypad

        n_center = (self.N - (self.N % 2)) // 2

        # Extract the segment of ypad centered at max_sample
        segment = self.ypad[max_sample + np.arange(self.N) - n_center]

        # Compute FFT on the windowed segment
        fft_result = np.fft.fft( np.hanning(self.N)* segment, self.Ndft)
        
       
        # Compute maximum spectral magnitude in dB
        maxS = 20 * np.log10(np.max(np.abs(fft_result)))
        # Adjust relative amplitude threshold
        self.G_g = maxS - self.G_g



    def _costMatrixComputation(self,Alpha_last,Alpha):

        """
        Compute the cost matrix used in the assignement problem 

        Input : estimated complex paramater of frame k and frame k-1 

        Output : Cost matrix of size : ()
        """

        # Amplitudes/Frequencies at midpoint of analysis frames 

        n_center = (self.N-self.N%2) /2 +1
        n_midpoint = (np.array([-np.floor(self.H/2),np.ceil(self.H/2)]) + n_center).astype(np.int32)
  
        time1 = functions.time_vector_computation(n_midpoint[1],self.Q)   
        time2 = functions.time_vector_computation(n_midpoint[0],self.Q) 

        dTime1 = functions.derivative_time_vector_computation(n_midpoint[1],self.Q)
        dTime2 = functions.derivative_time_vector_computation(n_midpoint[0],self.Q)

        # a(+H/2) [k-1]
        mA1 = np.array(np.expand_dims(np.real(time1 @ Alpha_last),axis=1))
        # a(-H/2) [k]
        mA2 = np.array(np.expand_dims(np.real(time2 @ Alpha),axis=1))

        # f(+H/2)
        mF1 = np.array(np.expand_dims(self.fs/(2*np.pi) * np.imag(dTime1 @ Alpha_last[1:, :]),axis=1))
        # f(-H/2)
        mF2 = np.array(np.expand_dims(self.fs/(2*np.pi) * np.imag(dTime2 @ Alpha[1:, :]),axis=1))


        #eq. (10)
        deltaA = mA1.T - mA2
        #eq. (11)
        deltaF = mF1.T - mF2 

        #eq. (8)
        Type = "A"
        Type = "B"
        A_useful = 1 - np.exp(-deltaF**2/(2*self.varF**2)-deltaA**2/(2*self.varA**2))
        B_spurious = 1 - (1-self.delta)*A_useful 
        
        #eq. (14)
        # Compute the element-wise minimum values
        min_values = np.minimum(A_useful, B_spurious)

        # Determine which matrix contributed the minimum value
        labels = np.where(A_useful <= B_spurious, "A", "B")

        # Combine min values and labels into final structured format
        costMatrix = np.stack((min_values, labels), axis=-1)

        return costMatrix
    
    def _parameterEstimation(self,signal):

        """
        Compute the estimated complex parameters alphas for each frame

        Input : 

        Output : vector Alpha of estimated complex parameter
        """

        #win = np.hanning(self.N)
        #winD = np.gradient(win)
        # win =np.array(blackmanharris(self.N,sym = False))
        # winD = np.gradient(win)

        Ndft = self.Ndft
        ndft = np.arange(0,Ndft,1)
        omega = 2 * np.pi * ndft/Ndft
        t = np.arange(self.N) / self.fs

        n_center = (self.N - (self.N % 2)) // 2
        centering = np.exp(1j * omega * n_center)
        

        n = np.expand_dims(np.arange(self.N),axis=1) - n_center
        ft_mat = np.exp(1j*np.outer(n , omega[:Ndft//2].T.conj()))

        p = functions.time_vector_computation(n,self.Q)
        pD = functions.derivative_time_vector_computation(n,self.Q)
        
        Rddm = self.Q*self.overSample
        Alpha, num_peaks, S = ddm.ddm(signal, self.Q, Rddm, self.G_g, self.win, self.winD, p, pD, centering, self.Ndft, ft_mat, omega)
        
        return Alpha, num_peaks, S[:,0]

    

    def partialTracking(self) :
    
        # define center of the analysis window 
        n_center = (self.N - self.N%2)/2 + 1 

        # Intilialize Alpha 
        Alpha = [] 
        num_peaks = []

        # Initialization for trajectory labelling : 

        Partials = {}
        pairs_ = []
        tracks_ = []
        num_tracks = 0
        num_active = 0
        Spectro = np.zeros((self.M,self.Ndft),dtype=complex)
        for m in range(self.M) :

            # SHORT-TERM PARAMETER ESTIMATION
            
            # time instant for frame m 
            self.time[m] = m*self.H + n_center

            # Short-term signal to estimate parameters 
            yst = self.ypad[int(self.time[m]-n_center):int(self.time[m]+self.N-n_center)]
            
            # Saving the estimate of the previous frame : 
            Alpha_last = Alpha  
            num_peaks_last = num_peaks
            
            # Estimate short-term model parameters of each peak using DDM
            
            Alpha, num_peaks, Spectro[m,:] = self._parameterEstimation(yst)
            
            # plt.close()
            # plt.figure()
            # plt.plot(abs(S[:self.Ndft//2]))
            # plt.show()


            # PARTIAL TRACKING USING HUNGARIAN ALGORITHM

            if ( m > 0 ) : # we cannot track the first frame

                ## TODO : Find the assignement problem index solution
                # Identify in the solution the non zeros (not necessary with the scipy solver)

                num_assignement = 0

                if num_peaks > 0 and num_peaks_last > 0 :

                    # Cost matrix computation for the frame m
                    cost_matrix = self._costMatrixComputation(Alpha_last,Alpha)

                    if m == 40 : 
                        np.save("cost_matrix_sin8",cost_matrix)

                    np.save("cost",cost_matrix)
                    # Scipy linear_solve to solve the assignement problem
                
                    column, row = linear_sum_assignment(cost_matrix[:,:,0].astype(float))
                    
                    final_indx = []

                    for k in range(len(column)) : 
                        if (cost_matrix[column[k],row[k],1]=="A"):
                            final_indx.append([row[k],column[k]])

                    Assignments = np.array(final_indx)
                    num_assignement = len(final_indx)


                    # Trajectory labelling

                    num_active_last = num_active
                    pairs_last = pairs_
                    tracks_last = tracks_ 

                    num_active = 0
                    pairs_ = np.zeros((num_assignement,1))  # paires => peak assignment label
                    tracks_ = np.zeros((num_assignement,1)) # tracks => trajectory label
                    

                    if num_assignement > 0 : 

                        # An assignement either continues an existing partials or create a new one (birth)

                        toContinue, _, toBirth = functions.match_sets(pairs_last, Assignments[:, 0])

                        num_toContinue = len(toContinue)
                        num_toBirth = len(toBirth)
                        num_active = num_toContinue + num_toBirth

                        pairs_ = Assignments[np.concatenate((toContinue[:, 1], toBirth)), 1]
                        #tracks_ = np.concatenate((tracks_last[toContinue[:, 0]], num_tracks + np.arange(1, num_toBirth + 1)))
                        
                        # Handle empty toContinue separately
                        if toContinue.size > 0:
                            toContinue_first_col = toContinue[:, 0]  # Extract indices safely
                            continued_tracks = tracks_last[toContinue_first_col]
                        else:
                            continued_tracks = np.array([], dtype=int)  # Ensure an empty integer array

                        if num_toBirth > 0:
                            new_tracks = num_tracks + np.arange(num_toBirth)
                        else:
                            new_tracks = np.array([], dtype=int)  # Ensure an empty integer array

                        # Concatenate only when we have valid values

                        tracks_ = np.concatenate((continued_tracks, new_tracks))
                        ii_active_last = num_active_last + np.arange(num_toBirth)
                        
                        # Save partials for frame m-1 and m
                        
                        # polynomiam time vectors
                        p = functions.time_vector_computation(n_center,self.Q)
                        pD = functions.derivative_time_vector_computation(n_center,self.Q)

                        # Initialization of the entry of the dictionnary

                        if (m-1) not in Partials:
                            size_m_1 = num_active_last + num_toBirth  
                            Partials[m-1] = np.zeros((size_m_1, 4))  

                        # add zeros at the end of partials at frame k-1 to get the new birth of partials ...
                        Partials[m-1] = np.vstack((Partials[m-1],np.zeros((int(num_toBirth),4))))

                        if m not in Partials:
                            size_m = num_active  
                            Partials[m] = np.zeros((size_m, 4))  
                        
                        Partials[m-1][ii_active_last, 0:3] = functions.FreqAmpPha(self.fs,Alpha_last[:, Assignments[toBirth, 0]], p, pD)
                        Partials[m-1][ii_active_last, 3] = tracks_[num_toContinue:num_active]

                        Partials[m][0:num_active, 0:3] = functions.FreqAmpPha(self.fs, Alpha[:, pairs_], p, pD)
                        Partials[m][0:num_active, 3] = tracks_
                    
                        num_active_last += num_toBirth
                        num_tracks += num_toBirth

            if m not in Partials :
                Partials[m] = np.zeros((1,4))           
        return Spectro, Partials
