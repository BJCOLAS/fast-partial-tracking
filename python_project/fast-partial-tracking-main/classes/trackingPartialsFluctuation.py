"""
Documentation
============

This class provides the tracking of the partial using the descriptor from the snail analyseur

Inspired by the matlab code from : J. Neri and P. Depalle. "Fast partial tracking of audio with real-time capability through linear programming." Proceedings of the International Conference on Digital Audio Effects (DAFx). 2018.

"""


# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------


import numpy as np
from utilities import ddm, functions, munkres 
from pprint import pprint
from scipy.optimize import linear_sum_assignment
from scipy.signal.windows import blackmanharris, hann
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d 


# ---------------------------------------------------------------------------
# Classes
# ---------------------------------------------------------------------------


class trackingPartialsSnailFluctuation:

    """
    Implement a tracking method of partials : 

        Involves the following two main procedures for each short-term frame
    
     1. Read the values from the snail parameters estimations and compute the cost matrix / cost function
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

    def __init__(self, deltaPhi, Loudness, PhiDem, fs , N , hopSize, G_g=-40 , Q=1 , delta=0.2 , zetaF = 50 , zetaA = 15, WINDOW_INTERP=False):
        """Create a new instance"""

        # Sampling rate
        self.fs = fs 
        # Window size 
        self.N = N
        # Hopsize 
        self.Hopsize = hopSize

        self.deltaPhi = deltaPhi
        self.Loudness = Loudness
        self.PhiDem = 2*np.pi*PhiDem   
        #self.PhiDem = PhiDem

        self.G_g = np.abs(G_g)
        self.Q = max(1,round(Q))
        self.zetaF = max(1e-1, zetaF)
        self.zetaA = max(1e-1, zetaA)/20 * np.log(10)
        self.delta = delta 

        # Raffinement 

        self.window_interp = WINDOW_INTERP
        #self.window_interp = False


    def setup(self) : 
        """
        Setup method of the main variables use of the parameter estimation
        """
        # Compute variance for tracking
        self.varF = -self.zetaF**2 * np.log((self.delta - 1) / (self.delta - 2))
        self.varA = -self.zetaA**2 * np.log((self.delta - 1) / (self.delta - 2))

        # Compute the frequency vector 
        #self.frequencies = functions.compute_frequencies()
        self.frequencies = functions.compute_frequencies_fourier()

        # Number of frame 

        self.M = self.Loudness.shape[0]

        # Sampling Period 
        self.T = 1/self.fs

        # time vector 
        self.time = np.zeros(self.M)

        # store last loudness trame for dervivative computation
        self.trame_log_amp_last = np.zeros(len(self.frequencies))

        # center of a frame 

        #self.n_center = int((self.N-self.N%2)/2 + 1) 
        #self.n_center = self.N//2
        self.n_center = 0

        # max value for thresholding 

        self.G_g =  np.max(self.Loudness) - self.G_g
        #self.G_g = - self.G_g

        ## window for amplitude interpolation ###
        window = hann(self.N)
        Ndft_interp = 4096 * 20 
        nf = 1/np.sum(window)
        self.window_hat = nf*np.fft.fftshift(np.fft.fft(window,Ndft_interp))
        self.window_freq = np.fft.fftshift(np.fft.fftfreq(Ndft_interp,1/self.fs))

        ## jitter
        # Longueur totale en samples
        total_length = self.M  # ou self.M * self.Hopsize si tu préfères

        # Jitter global en omega normalisé (à adapter selon ton domaine)
        sigma = 30 # très petit pour éviter le chaos
        t_ms = 0.15
        g_poisson = 1
        fs = 48000
        lambda_poisson = 1 / (t_ms * fs / 1000) * g_poisson

        # Instants du processus de Poisson
        inter_arrivals = np.random.exponential(scale=1 / lambda_poisson, size=total_length)
        instants = np.cumsum(inter_arrivals)
        instants = instants[instants < total_length]
        instants = np.insert(instants, 0, 0)
        instants = np.append(instants, total_length - 1)

        # Valeurs de jitter gaussien
        jitter_values = np.random.normal(loc=0, scale=sigma, size=len(instants))

        # Interpolation linéaire
        interp_fn = interp1d(instants, jitter_values, kind='linear')
        tvec = np.arange(total_length)
        self.global_jitter = interp_fn(tvec)


        ## shimmer
        # Longueur totale 
        total_length = self.M 

        # Parametres du shimmer
        sigma = 0.2 # très petit pour éviter le chaos
        t_ms = 0.1
        g_poisson = 1
        fs = 48000
        lambda_poisson = 1 / (t_ms * fs / 1000) * g_poisson

        # Processus de poisson
        inter_arrivals = np.random.exponential(scale=1 / lambda_poisson, size=total_length)
        instants = np.cumsum(inter_arrivals)
        instants = instants[instants < total_length]
        instants = np.insert(instants, 0, 0)
        instants = np.append(instants, total_length - 1)

        # Valeurs jitter loi Gausienne
        shimmer_values = np.random.normal(loc=0, scale=sigma, size=len(instants))

        # Interpolation linéaire
        interp_fn = interp1d(instants, shimmer_values, kind='linear')
        tvec = np.arange(total_length)
        self.global_shimmer = interp_fn(tvec)




    def _parameterEstimation(self,trame_deltaPhi, trame_Loudness, trame_PhiDem,trame_log_amp_last):

        """
        Compute the estimated parameters for each frame using the Snail Analyser estimator

        Input : 

        Output : vector Alpha of estimated complex parameter
        """
        
        # f_true and f_true_index computation for one trame
        true_f_, true_f_index= functions.find_f_true(trame_deltaPhi,self.frequencies)

        # Amplitude computation for one trame
        amp_true = functions.interp_amp(trame_Loudness,true_f_,true_f_index,self.frequencies)

        # Amplitude thresholding
        amp_true, true_f_, true_f_index = functions.amplitude_thresholding(self.G_g, amp_true, true_f_, true_f_index)

        # amplitude interpolation with fourier window

        if self.window_interp == True : 
            amp_true = functions.window_amplitude_interpolation(trame_Loudness, true_f_, true_f_index,self.frequencies,self.window_hat,self.window_freq)

        # Amplitude conversion dB to logamp
        
        log_amp_true = amp_true * np.log(10) / 20

        # Derivative of amplitude computation for one frame
        damp_true = functions.interp_damp(log_amp_true,trame_log_amp_last,true_f_,self.frequencies,true_f_index,self.Hopsize*self.T)
        

        # Phase computation for one frame
        phase_true = functions.interp_phase2(trame_PhiDem,true_f_,true_f_index,self.frequencies)

        # Number of peaks estimated for one frame
        num_peaks = len(true_f_)

        return true_f_,log_amp_true,damp_true,phase_true, num_peaks
    

    def _costMatrixComputation(self,parameter_last,parameter):

        """
        Compute the cost matrix used in the assignement problem 

        Input : estimated complex paramater of frame k and frame k-1 

        Output : Cost matrix of size : ()
        """

        true_f_1,amp_true1,damp_true1,_,_= parameter_last
        true_f_2,amp_true2,damp_true2,_,_ = parameter

        # midpoint amplitude computation

        # mA1 = np.array(np.expand_dims(amp_true1 + damp_true1*(self.Hopsize/(2*self.fs)),axis=1))
        # mA2 = np.array(np.expand_dims(amp_true2 + damp_true2*(-self.Hopsize/(2*self.fs)),axis=1))

        mA1 = np.array(np.expand_dims(amp_true1,axis=1))
        mA2 = np.array(np.expand_dims(amp_true2,axis=1))

        # mid frame frequency computation (we don't have access to the derivative of the frequency)

        mF1 = np.array(np.expand_dims(true_f_1,axis=1))
        mF2 = np.array(np.expand_dims(true_f_2,axis=1))

        # mid frame dA computation 

        deltaA = (mA1.T - mA2).T
        
        # mid frame dF computation

        deltaF = (mF1.T - mF2).T
        
        #eq. (8)
        Type = "A"
        Type = "B"
        A_useful = 1 - np.exp(-(deltaF**2)/(self.varF)-(deltaA**2)/(self.varA))
        B_spurious = 1 - (1-self.delta)*A_useful 
        
        #eq. (14)
        # Compute the element-wise minimum values
        min_values = np.minimum(A_useful, B_spurious)

        # Determine which matrix contributed the minimum value
        labels = np.where(A_useful <= B_spurious, "A", "B")

        # Combine min values and labels into final structured format
        costMatrix = np.stack((min_values, labels), axis=-1)

        return costMatrix
    

    def partialTracking(self) : 
        
        # initialisation 
        parameter = [] 
        num_peaks = []
        Partials = {}
        pairs_ = []
        tracks_ = []
        num_tracks = 0
        num_active = 0

        for m in range(self.M) :

            # SHORT-TERM PARAMETER ESTIMATION

            # Saving the estimate of the previous frame : 
            parameter_last = parameter  
            num_peaks_last = num_peaks

            # Time instant : 
            #self.time[m] = m * self.Hopsize / self.fs + self.Hopsize / (2 * self.fs) # middle of the frame

            self.time[m] = self.n_center + (m)* self.Hopsize 

            # Compute the estimated parameters for each frame using the Snail Analyser estimator
            
            parameter = self._parameterEstimation(self.deltaPhi[m,:],self.Loudness[m,:],self.PhiDem[m,:],self.trame_log_amp_last)

            ## Ajout du jitter ##


            freqs = parameter[0]
            freqs_jittered = freqs + self.global_jitter[m]
            #freqs_jittered = freqs 

            amps = parameter[1]
            #amps_shimmered = amps + self.global_shimmer[m]
            amps_shimmered = amps
            # Mettre à jour
            parameter = (freqs_jittered, amps_shimmered, parameter[2], parameter[3], parameter[4])



                    #####################


            # Store the last loudness converted as log_amp trame for derivative computation
            self.trame_log_amp_last = self.Loudness[m,:]*np.log(10)/20

    
            num_peaks = parameter[4]

            # PARTIAL TRACKING USING HUNGARIAN ALGORITHM

            if ( m > 0 ) : # we cannot track the first frame

                num_assignement = 0

                if num_peaks > 0 and num_peaks_last > 0 :

                    # Cost matrix computation for the frame m
                    cost_matrix = self._costMatrixComputation(parameter,parameter_last)
                    cm = cost_matrix[:,:,0].astype(np.float64)  # ou simplement float
                    #np.savetxt("cost_matrix.csv", cm, delimiter=",", fmt="%.18e")
                    #np.savetxt("cost_matrix.csv", cost_matrix[:,:,0], delimiter=",")

                    # Scipy linear_solve to solve the assignement problem
                
                    column, row = linear_sum_assignment(cost_matrix[:,:,0].astype(float))

                    final_indx = []

                    for k in range(len(column)) : 
                        if (cost_matrix[column[k],row[k],1]=="A"):
                            final_indx.append([row[k],column[k]])

                    Assignments = np.array(final_indx)
                    num_assignement = len(final_indx)

                    num_active_last = num_active
                    pairs_last = pairs_
                    tracks_last = tracks_ 

                    num_active = 0
                    pairs_ = np.zeros((num_assignement,1))  # paires => peak assignment label
                    tracks_ = np.zeros((num_assignement,1)) # tracks => trajectory label
                    
                    # trajectory labelling : 

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
                        
                        # Initialization of the entry of the dictionnary

                        # if (m-1) not in Partials:
                        #     size_m_1 = num_active_last + num_toBirth  
                        #     Partials[m-1] = np.zeros((size_m_1, 5))  

                        # # add zeros at the end of partials at frame k-1 to get the new birth of partials ...
                        # Partials[m-1] = np.vstack((Partials[m-1],np.zeros((int(num_toBirth),5))))
                        size_needed = num_active_last + num_toBirth

                        if (m-1) not in Partials:
                            Partials[m-1] = np.zeros((size_needed, 5))
                        else:
                            # Si le buffer existe mais est trop petit, on l’agrandit
                            current_size = Partials[m-1].shape[0]
                            if current_size < size_needed:
                                missing = size_needed - current_size
                                Partials[m-1] = np.vstack((Partials[m-1], np.zeros((missing, 5))))

                        if m not in Partials:
                            size_m = num_active  
                            Partials[m] = np.zeros((size_m, 5))  
                        
                        # store the values of the birth of partials in the dictionnary of the previous frame
                        true_f_1,amp_true1,damp1,pha_true1,_= parameter_last
                        omega_1 = 2*np.pi*true_f_1/self.fs

                        true_f_2,amp_true2,damp2,pha_true2,_ = parameter
                        omega_2 = 2*np.pi*true_f_2/self.fs
         
                        Partials[m-1][ii_active_last, 0:4] = np.array([omega_1[Assignments[toBirth, 0]], amp_true1[Assignments[toBirth, 0]], pha_true1[Assignments[toBirth, 0]],damp1[Assignments[toBirth, 0]]]).T
                        Partials[m-1][ii_active_last, 4] = tracks_[num_toContinue:num_active]

                        Partials[m][0:num_active, 0:4] = np.array([omega_2[pairs_], amp_true2[pairs_], pha_true2[pairs_],damp2[pairs_]]).T
                        Partials[m][0:num_active, 4] = tracks_
                    
                        num_active_last += num_toBirth
                        num_tracks += num_toBirth

            if m not in Partials :
                Partials[m] = np.zeros((1,5))
        return Partials
