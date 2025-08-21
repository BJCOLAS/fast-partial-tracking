"""
Documentation
============

"""


# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------


import numpy as np
from utilities.functions import match_sets
from classes.trackingPartialsSnail import trackingPartialsSnail
from scipy.interpolate import interp1d

# ---------------------------------------------------------------------------
# Classe
# ---------------------------------------------------------------------------


class synthesisDecoder:

    def __init__(self,parameters,tracker: trackingPartialsSnail):
        self.parameters = parameters
        self.tracker = tracker
    
    def _compute_amp_per_id(self):
        amp = {}
        time = self.tracker.time
        
        for frame in self.parameters : 
            if frame not in amp:
                amp[frame] = []
            if frame == 0 : 
                continue 
            for track in self.parameters[frame] :

                time_m = [time[frame-1], time[frame]]
                H = abs(np.diff(time_m))[0]
                t_span = [0, H]
                t = np.arange(t_span[0], t_span[1])
                Amplitude = track["Amp_m"]
                Amp_ = np.log(interp1d(t_span, np.exp(Amplitude), kind='linear')(t)+1e-10)

                amp[frame].append({"track_id":track["track_id"],
                                   "amplitudes":Amp_})
                
        return amp
        
    def _compute_phase_per_id(self):

        time = self.tracker.time
        frame_phase = {} 
        for frame in self.parameters :

            if frame not in frame_phase:
                frame_phase[frame] = []
            if frame == 0 : 
                continue 
            for track in self.parameters[frame] :

                time_m = [time[frame-1], time[frame]]
                H = abs(np.diff(time_m))[0]
                t_span = [0, H]
                t = np.arange(t_span[0], t_span[1])

                Omega = track["Omega_m"] 
  
                Phase = track["Phase_m"]

                M_star = (H*Omega[0]+H*Omega[1] + 2 *(Phase[0]-Phase[1]))/(4*np.pi)
                M_star = np.floor(M_star+0.5)

                gamma_0  = 1j*Phase[0]
                gamma_prime_0 = 1j*Omega[0]
                gamma_1 = 1j*Phase[1]
                gamma_prime_1 = 1j*Omega[1]

                P_bar = np.array([[gamma_0],
                                  [gamma_prime_0],
                                  [6*M_star*1j*np.pi/H**2 - 2*gamma_prime_0/H - 3*gamma_0/H**2- gamma_prime_1/H + 3*gamma_1/H**2],
                                  [-4*M_star*1j*np.pi/H**3 + gamma_prime_0/H**2 + 2*gamma_0/H**3 + gamma_prime_1/H**2 - 2*gamma_1/H**3]])

                Phase_ = np.imag(P_bar[0,0]) + np.imag(P_bar[1,0]) * t + np.imag(P_bar[2,0]) * t**2 + np.imag(P_bar[3,0]) * t**3

                frame_phase[frame].append({"track_id" : track["track_id"],
                                           "phase" : Phase_})
        return frame_phase
    
    def _compute_freq_per_id(self):

        time = self.tracker.time
        frame_freq = {} 
        for frame in self.parameters :

            if frame not in frame_freq:
                frame_freq[frame] = []
            if frame == 0 : 
                continue 
            for track in self.parameters[frame] :

                time_m = [time[frame-1], time[frame]]
                H = abs(np.diff(time_m))[0]
                t_span = [0, H]
                t = np.arange(t_span[0], t_span[1])

                Omega = track["Omega_m"] 

                Phase = track["Phase_m"]

                M_star = (H*Omega[0]+H*Omega[1] + 2 *(Phase[0]-Phase[1]))/(4*np.pi)
                M_star = np.floor(M_star+0.5)

                gamma_0  = 1j*Phase[0]
                gamma_prime_0 = 1j*Omega[0]
                gamma_1 = 1j*Phase[1]
                gamma_prime_1 = 1j*Omega[1]

                P_bar = np.array([[gamma_0],
                                    [gamma_prime_0],
                                    [6*M_star*1j*np.pi/H**2 - 2*gamma_prime_0/H - 3*gamma_0/H**2- gamma_prime_1/H + 3*gamma_1/H**2],
                                    [-4*M_star*1j*np.pi/H**3 + gamma_prime_0/H**2 + 2*gamma_0/H**3 + gamma_prime_1/H**2 - 2*gamma_1/H**3]])

                Freq_ = self.tracker.fs/(2*np.pi) * (np.imag(P_bar[1,0]) + 2 *np.imag(P_bar[2,0]) * t + 3 * np.imag(P_bar[3,0]) * t**2)

                frame_freq[frame].append({"track_id" : track["track_id"],
                                          "frequency" : Freq_})
        return frame_freq
    
    
    def compute_output(self):
        phases = self._compute_phase_per_id()
        amplitudes = self._compute_amp_per_id()
        frequencies = self._compute_freq_per_id()
        return phases, amplitudes, frequencies

