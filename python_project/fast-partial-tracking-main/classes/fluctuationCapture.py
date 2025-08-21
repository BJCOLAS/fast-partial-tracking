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
import matplotlib.pyplot as plt
from scipy.signal import firwin, lfilter

# ---------------------------------------------------------------------------
# Classe
# ---------------------------------------------------------------------------


class fluctuationCapture:

    def __init__(self,frequencies):
        self.frequencies = frequencies


    def applyFiltering(self,H,L):

    # 1. Récupérer tous les track_id uniques
        track_ids = sorted({
            entry["track_id"]
            for frame in list(self.frequencies.values())
            for entry in frame
        })

        # 2. Préparer les buffers (phase, amplitude, signal)
        frequency_buffer  = {tid: np.full(L, np.nan) for tid in track_ids}


        # 3. Remplir les buffers frame par frame
        for frame, p_tracks in self.frequencies.items():
            start = frame * H

            if start < L :
                for p_track in p_tracks:
                    tid = p_track["track_id"]
                    fh  = p_track["frequency"]


                    # Vérifications
                    if fh is None or fh.size == 0:
                        continue

                    # Tronquer à la taille minimale
                    n   = fh.size
                    end = min(start + n, L)
                    f_seg = fh[:n]

                    # Remplissage des buffers
                    frequency_buffer[tid][start:end]  = f_seg

        signal = frequency_buffer[0]
        signal = np.nan_to_num(signal)

        plt.plot(signal)
        plt.show()
        # Filtres 

        fs = 48000 # Hz
        fc_hz = 0.5  # Hz (fréquence de coupure réelle)
        fc = fc_hz / (fs / 2)  # Normalisation pour firwin
        N = 601  # Filtre plus large pour meilleure séparation

        h_lp = firwin(N, fc, window='hamming')  # Passe-bas
        delta = np.zeros(N)
        delta[N // 2] = 1
        h_hp = delta - h_lp  # Passe-haut apparié

        slow_fluct = lfilter(h_lp, 1.0, signal)
        fast_fluct = lfilter(h_hp, 1.0, signal)
        
        time_vec = np.arange(len(slow_fluct))/fs
        # 6. Affichage
        plt.figure(figsize=(10, 5))
        plt.plot(time_vec,signal, label='Fréquences du premier partiel', alpha=0.5)
        plt.plot(time_vec,slow_fluct, label='Fluctuations lentes (passe-bas)', linewidth=2)
        plt.plot(time_vec,fast_fluct, label='Fluctuations rapides (passe-haut)', linewidth=2)
        plt.legend()
        plt.title("Séparation des fluctuations")
        plt.xlabel("temps(s)")
        plt.xlim([0,0.1])
        plt.ylabel("fréquences(Hz)")
        plt.grid()
        plt.show()


        #plt.figure()
        #plt.plot(signal)
        #plt.show()

