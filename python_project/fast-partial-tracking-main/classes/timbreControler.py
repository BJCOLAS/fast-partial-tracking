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
from utilities.lfo import lfo_sin
import matplotlib.pyplot as plt
# ---------------------------------------------------------------------------
# Classe
# ---------------------------------------------------------------------------

class timbreController:
    def __init__(self, amplitudes: dict, phases: dict, freqs: dict):
        self.amplitudes = amplitudes
        self.phases     = phases
        self.freqs      = freqs

    def apply_LFO_global(self,
                         lfo_type: str,
                         rate: float,
                         sigma_hz: float,
                         fs: float) -> dict:

        if lfo_type != "sin":
            raise ValueError(f"LFO '{lfo_type}' non implémenté")


        frames       = sorted(self.phases.keys())
        H            = len(self.phases[frames[0]][0]["phase"])
        total_length = frames[-1] * H


        lfo = lfo_sin(rate, sigma_hz, total_length, fs)


        track_ids = {
            entry["track_id"]
            for frame in frames
            for entry in self.phases[frame]
        }


        phases_out = { frame: [] for frame in frames }

        for tid in sorted(track_ids):

            f_full   = np.zeros(total_length)
            phi_orig = np.zeros(total_length)
            for frame in frames:
                start = frame * H
                end   = start + H

                for e in self.freqs.get(frame, []):
                    if e["track_id"] == tid:
                        if end <= total_length:
                            f_full[start:end] = e["frequency"]
                        break

                for e in self.phases[frame]:
                    if e["track_id"] == tid:
                        if end <= total_length:
                            phi_orig[start:end] = e["phase"]
                        break

            f_mod = f_full + lfo

            phi_full = np.zeros(total_length)
            acc       = 0.0
            prev_zero = True

            for i, f in enumerate(f_mod):
                if f_full[i] == 0:
                    acc = 0.0
                    prev_zero = True
                else:
                    if prev_zero:

                        acc = phi_orig[i]
                        prev_zero = False

                    acc += 2 * np.pi * f / fs
                phi_full[i] = acc


            for frame in frames:
                start = frame * H
                end   = start + H
                if any(e["track_id"] == tid for e in self.freqs.get(frame, [])):
                    phases_out[frame].append({
                        "track_id": tid,
                        "phase": phi_full[start:end]
                    })

        return phases_out

    def apply_LFO_test(self,fs):

        phases_out = []

        frames = sorted(self.phases.keys())
        H = len(self.phases[frames[1]][0]["phase"])  
        total_length = frames[-1]* H

        track_ids = {entry["track_id"] for frame in frames for entry in self.phases[frame]}

        for tid in sorted(track_ids):

            f_full     = np.zeros(total_length)
            phi_orig   = np.zeros(total_length)

            for frame in frames:

                for entry in self.freqs[frame]:
                    if entry["track_id"] == tid:
                        start, end = frame * H, frame * H + H
                        if end <= total_length:
                            f_full[start:end] = entry["frequency"]
                        break

                for entry in self.phases[frame]:
                    if entry["track_id"] == tid:
                        start, end = frame * H, frame * H + H
                        if end <= total_length:
                            phi_orig[start:end] = entry["phase"]
                        break


            phi_integrated = np.zeros_like(f_full)
            acc = 0.0

            for i, f in enumerate(f_full):
                    if f == 0:
                        acc = 0.0
                        prev_zero = True
                    else:
                        if prev_zero:
                            phi0 = phi_orig[i]
                            acc  = phi0
                            prev_zero = False

                        acc += 2 * np.pi * f / fs

                    phi_integrated[i] = acc


            # Figure : 
            plt.figure(figsize=(10, 3))
            plt.plot(phi_orig,     '--', label=f"Phase originale (tid={tid})")
            plt.plot(phi_integrated, '-', label=f"Phase intégrée (tid={tid})")
            plt.title("Phase originale vs phase intégrée")
            plt.xlabel("Échantillons")
            plt.ylabel("Phase (radians)")
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.show()

            signal_orig = np.exp(1j * phi_orig)
            signal_int  = np.exp(1j * phi_integrated)

            plt.figure(figsize=(10, 3))
            plt.plot(np.real(signal_orig), '--', label=f"Re(e^{'{j·φ_orig}'}) tid={tid}")
            plt.plot(np.real(signal_int),  '-', label=f"Re(e^{'{j·φ_int}'}) tid={tid}")
            plt.title("Oscillations — Partie réelle des exponentielles")
            plt.xlabel("Échantillons")
            plt.ylabel("Amplitude")
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.show()

        return phases_out


    def apply_LFO_jitter_per_track2(self,
                                    sigmacent: float,
                                    t_ms: float,
                                    poisson_gain: float,
                                    fs: float,
                                    total_length: int) -> dict:
        from scipy.interpolate import interp1d

        ln2_over_1200 = np.log(2) / 1200
        lambda_poisson = (1 / (t_ms * fs / 1000)) * poisson_gain

        frames = sorted(f for f in self.freqs if f != 0)
        H = len(self.freqs[frames[0]][0]["frequency"])

        track_ids = {entry["track_id"] for f in frames for entry in self.freqs[f]}
        phases_out = {frame: [] for frame in frames}

        for tid in track_ids:
            # -- reconstruct freq signal
            f_full = np.zeros(total_length)
            phi0 = None
            for frame in frames:
                for entry in self.freqs[frame]:
                    if entry["track_id"] == tid:
                        start = frame * H
                        end = start + H
                        f_full[start:end] = entry["frequency"]
                        if phi0 is None:
                            pe = next(e for e in self.phases[frame] if e["track_id"] == tid)
                            phi0 = pe["phase"][0]
                        break

            sigma_local = sigmacent * ln2_over_1200 * f_full

            inter_arrivals = np.random.exponential(scale=1 / lambda_poisson, size=total_length)
            instants = np.cumsum(inter_arrivals)
            instants = instants[instants < total_length]
            instants = np.insert(instants, 0, 0)
            instants = np.append(instants, total_length - 1)

            jitter_vals = np.random.normal(loc=0, scale=sigma_local[instants.astype(int)], size=len(instants))
            jitter_fn = interp1d(instants, jitter_vals, kind='linear')
            jitter = jitter_fn(np.arange(total_length))

            f_mod = f_full + jitter
            phi_full = phi0 + 2 * np.pi * np.cumsum(f_mod) / fs

            for frame in frames:
                start = frame * H
                end = start + H
                phi_seg = phi_full[start:end]
                phases_out[frame].append({
                    "track_id": tid,
                    "phase": phi_seg
                })

        return phases_out
    

    def apply_LFO_jitter_per_track(self, sigmacent: float, t_ms: float, poisson_gain: float, fs: float, total_length: int) -> dict: # on informe les types pour faciliter l'intégraton et la complétion plus tard
        """
        Applique un jitter de fréquence à chaque partiel selon un
        processus Poisson/exponentiel, avec intégration de phase

        """

        ln2_over_1200 = np.log(2) / 1200.0
 
        frames = sorted(f for f in self.freqs if f != 0)

        valid_frames = [f for f in frames if len(self.freqs[f]) > 0]
        if not valid_frames:
            raise ValueError("Aucune frame non-vide dans self.freqs")


        H = len(self.freqs[valid_frames[0]][0]["frequency"])


        lambda_poisson = (1.0 / (t_ms * fs / 1000.0)) * poisson_gain

        track_ids  = {
            e["track_id"]
            for f in frames
            for e in self.freqs[f]
        }
        phases_out = {frame: [] for frame in frames}

        for tid in sorted(track_ids):

            f_full   = np.zeros(total_length)
            phi_orig = np.zeros(total_length)
            phi0     = None

            for frame in frames:
                start = frame * H
                end   = start + H

                for e in self.freqs.get(frame, []):
                    if e["track_id"] == tid:
                        f_full[start:end] = e["frequency"]
                        break

                if phi0 is None:
                    pe = next((e for e in self.phases[frame]
                               if e["track_id"] == tid), None)
                    if pe is not None:
                        phi0 = pe["phase"][0]
                # on stocke la phase d'origine 
                for e in self.phases.get(frame, []):
                    if e["track_id"] == tid:
                        phi_orig[start:end] = e["phase"]
                        break

            if phi0 is None:
                continue

            sigma_local = sigmacent * ln2_over_1200 * f_full

            # instants de déviation Poisson/expo
            inter_arrivals = np.random.exponential(
                scale=1.0 / lambda_poisson, size=total_length)
            instants = np.cumsum(inter_arrivals)
            instants = instants[instants < total_length].astype(int)
            instants = np.insert(instants, 0, 0)
            instants = np.append(instants, total_length - 1)

            # valeurs aléatoires de jitter aux instants 
            jitter_vals = np.random.normal(
                loc=0.0,
                scale=sigma_local[instants],
                size=instants.shape[0]
            )
            jitter_fn = interp1d(instants, jitter_vals,
                                 kind="linear",
                                 bounds_error=False,
                                 fill_value=0.0)
            jitter = jitter_fn(np.arange(total_length))

            # Intégration
            phi_full = np.zeros(total_length)
            acc       = 0.0
            prev_zero = True

            # f_mod = f_full + jitter 
            f_mod = f_full + jitter

            for i in range(total_length):
                if f_full[i] == 0:
                    acc = 0.0
                    prev_zero = True
                else:
                    if prev_zero:
                        acc = phi_orig[i] if phi_orig[i] != 0 else phi0
                        prev_zero = False
                    acc += 2.0 * np.pi * f_mod[i] / fs
                phi_full[i] = acc

            # découpage en frames pour la sortie
            for frame in frames:
                start = frame * H
                end   = start + H

                if any(e["track_id"] == tid for e in self.freqs.get(frame, [])):
                    phases_out[frame].append({
                        "track_id": tid,
                        "phase": phi_full[start:end]
                    })

        return phases_out
    