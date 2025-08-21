"""
Documentation
============

Generate the usefull parameter per track id for the interpolation with the management of birth, continuity and death of partials
"""


# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------


import numpy as np
from utilities.functions import match_sets

# ---------------------------------------------------------------------------
# Classes
# ---------------------------------------------------------------------------


class trackedParameterSynthesis:

    """
    Implement the generation of the right parameter struct

    Input :  
    Analysis Partials Dictionnary
    Time
    Global parameters 

    Output Synthesis Parameter : 
    - n° de trame de synthèse
    - liste ID track/partiels
    - Pour chaque ID : Param_phase, Param_A 

    Param_phase : Phi(k-1),Phi(k),Freq(k-1),Freq(k)
    Param_A : log_amp(k-1),log_amp(k)
    """

    def __init__(self, tracker, partials):
        """Create a new instance"""

        self.tracker = tracker
        self.partials = partials 
        self.time = tracker.time 
        self.fs = tracker.fs
        self.M = tracker.M
        self.Nw = tracker.Hopsize
    

    def transform(self):
        output = {}
        partials = self.partials
        time = self.time
        M = len(partials)
        num_active = 0 
        Nw = 2400 # window size
         # Centroïde temporel
        Tc = Nw / 2
        time = time + int(Tc)

        for m in range(1,M) : 
            if m >= 1360 :
            #test break
                a = time[m] 
            if time[m]< (self.tracker.M*self.tracker.Hopsize) :
                output[m] = []

                num_active_last = num_active
                num_active = np.count_nonzero(partials[m][:, 1])

                # naissance : 
                if num_active > 0 and num_active_last == 0:
                    time_m = [time[m-1], time[m]]
                    H = abs(np.diff(time_m))[0]
                    n_unmatched = num_active

                    for idx in range(n_unmatched):
                        
                        trackId = partials[m][idx, 4]

                        Omega_m = [partials[m][idx, 0], partials[m][idx, 0]]
                        Amp_m = [-1e10, partials[m][idx, 1]]

                        Phase_1_bary = partials[m][idx, 2] +  partials[m][idx, 0] * Tc
                        Phase_0_bary = Phase_1_bary - partials[m][idx, 0] * H
                        Phase_m = [Phase_0_bary , Phase_1_bary]

                        # stockage
                        output[m].append({
                            "track_id": trackId,
                            "Omega_m": Omega_m,
                            "Amp_m": Amp_m,
                            "Phase_m": Phase_m
                        })
                # mort : 


                # death of all partials
                elif num_active == 0 and num_active_last > 0:
                    time_m = [time[m-1], time[m]]
                    H = abs(np.diff(time_m))[0]
                    n_unmatched_last = num_active_last

                    for idx in range(n_unmatched_last):
                    # probleme : extinction 1 partiel => n_unmatched_last = 2 donc idx =1 pose probleme
                        trackId = partials[m-1][idx, 4]

                        Omega_m = [partials[m-1][idx, 0], partials[m-1][idx, 0]]
                        Amp_m = [partials[m-1][idx, 1], -1e10]
                        
                        Phase_0_bary = partials[m-1][idx, 2] +  partials[m-1][idx, 0] * Tc
                        Phase_1_bary = Phase_0_bary + partials[m-1][idx, 0]*H
                        Phase_m = Phase_m = [Phase_0_bary , Phase_1_bary]
                        
                        output[m].append({
                            "track_id": trackId,
                            "Omega_m": Omega_m,
                            "Amp_m": Amp_m,
                            "Phase_m": Phase_m
                        })
                
                # general
                elif num_active > 0 and num_active_last > 0:
                    matches, unmatched_last, unmatched = match_sets(partials[m-1][:num_active_last, 4], partials[m][:num_active, 4])
                    n_matches = len(matches)
                    n_unmatched_last = len(unmatched_last)
                    n_unmatched = len(unmatched)

                    time_m = [time[m-1], time[m]]
                    H = abs(np.diff(time_m))[0]

                    for idx in range(n_unmatched):

                        trackId = int(partials[m][unmatched[idx], 4])
                        Omega_m = [partials[m][unmatched[idx], 0], partials[m][unmatched[idx], 0]]
                        Amp_m = [-1e10, partials[m][unmatched[idx], 1]]
                        
                        Phase_1_bary = partials[m][unmatched[idx], 2] +  partials[m][unmatched[idx], 0] * Tc
                        Phase_0_bary = partials[m][unmatched[idx], 2] - partials[m][unmatched[idx], 0] * H
                        Phase_m = [Phase_0_bary , Phase_1_bary]
                        output[m].append({
                            "track_id": trackId,
                            "Omega_m": Omega_m,
                            "Amp_m": Amp_m,
                            "Phase_m": Phase_m
                        })

                    for idx in range(n_unmatched_last):

                        trackId = partials[m-1][unmatched_last[idx], 4]
                        Omega_m = [partials[m-1][unmatched_last[idx], 0], partials[m-1][unmatched_last[idx], 0]]
                        Amp_m = [partials[m-1][unmatched_last[idx], 1], -1e10]
                        
                        Phase_0_bary = partials[m-1][unmatched_last[idx], 2] +  partials[m-1][unmatched_last[idx], 0] * Tc
                        Phase_1_bary = Phase_0_bary + partials[m-1][unmatched_last[idx], 0]*H
                        Phase_m = [Phase_0_bary , Phase_1_bary]
                        output[m].append({
                            "track_id": trackId,
                            "Omega_m": Omega_m,
                            "Amp_m": Amp_m,
                            "Phase_m": Phase_m
                        })

                    for idx in range(n_matches):
                        trackId = partials[m][matches[idx][1], 4]
                        Omega_m = [partials[m-1][matches[idx][0], 0], partials[m][matches[idx][1], 0]]
                        Amp_m = [partials[m-1][matches[idx][0], 1], partials[m][matches[idx][1], 1]]

                        Phase_0_bary = partials[m-1][matches[idx][0], 2] +  partials[m-1][matches[idx][0], 0] * Tc
                        Phase_1_bary = partials[m][matches[idx][1], 2] + partials[m][matches[idx][1], 0] * Tc
                        Phase_m = [Phase_0_bary , Phase_1_bary]
                        output[m].append({
                            "track_id": trackId,
                            "Omega_m": Omega_m,
                            "Amp_m": Amp_m,
                            "Phase_m": Phase_m
                        })
        return output 
