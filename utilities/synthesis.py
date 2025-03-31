import numpy as np
from utilities.functions import match_sets
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def jun_synthesize_partials(Partials, time, L):
    out = np.zeros(L)
    M = len(Partials)
    num_active = 0

    for m in range(M):
        num_active_last = num_active
        num_active = np.count_nonzero(Partials[m][:, 1])

        if num_active > 0 and num_active_last == 0:
            # Fade in all
            if m == 0:
                time_m = [1, time[m]]
            else:
                time_m = [time[m-1], time[m]]
            H = abs(np.diff(time_m))[0]
            n_unmatched = num_active
            for idx in range(n_unmatched):
                Omega_m = [Partials[m][idx, 0], Partials[m][idx, 0]]
                Amp_m = [-1e10, Partials[m][idx, 1]]
                Phase_m = [Partials[m][idx, 2] - Omega_m[1] * H, Partials[m][idx, 2]]

                ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                out[time_m[0]:time_m[1]] += ytemp[:-1]

        elif num_active == 0 and num_active_last > 0:
            time_m = [time[m-1], time[m]]
            H = abs(np.diff(time_m))[0]
            n_unmatched_last = num_active_last
            for idx in range(n_unmatched_last):
                Omega_m = [Partials[m-1][idx, 0], Partials[m-1][idx, 0]]
                Amp_m = [Partials[m-1][idx, 1], -np.inf]
                Phase_m = [Partials[m-1][idx, 2], Partials[m-1][idx, 2] + Omega_m[0] * H]

                ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                out[time_m[0]:time_m[1]+1] += ytemp

        elif num_active > 0 and num_active_last > 0:
            matches, unmatched_last, unmatched = match_sets(Partials[m-1][:num_active_last, 3], Partials[m][:num_active, 3])
            n_matches = len(matches)
            n_unmatched_last = len(unmatched_last)
            n_unmatched = len(unmatched)

            time_m = [time[m-1], time[m]]
            H = abs(np.diff(time_m))[0]

            for idx in range(n_unmatched):
                Omega_m = [Partials[m][unmatched[idx], 0], Partials[m][unmatched[idx], 0]]
                Amp_m = [-np.inf, Partials[m][unmatched[idx], 1]]
                Phase_m = [Partials[m][unmatched[idx], 2] - Omega_m[1] * H, Partials[m][unmatched[idx], 2]]

                ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                out[time_m[0]:time_m[1]] += ytemp[:-1]

            for idx in range(n_unmatched_last):
                Omega_m = [Partials[m-1][unmatched_last[idx], 0], Partials[m-1][unmatched_last[idx], 0]]
                Amp_m = [Partials[m-1][unmatched_last[idx], 1], -1e10]
                Phase_m = [Partials[m-1][unmatched_last[idx], 2], Partials[m-1][unmatched_last[idx], 2] + Omega_m[0] * H]

                ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                out[time_m[0]:time_m[1]+1] += ytemp

            for idx in range(n_matches):
                Omega_m = [Partials[m-1][matches[idx][0], 0], Partials[m][matches[idx][1], 0]]
                Amp_m = [Partials[m-1][matches[idx][0], 1], Partials[m][matches[idx][1], 1]]
                Phase_m = [Partials[m-1][matches[idx][0], 2], Partials[m][matches[idx][1], 2]]

                ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                out[time_m[0]:time_m[1]] += ytemp[:-1]
    return out

def jun_synthesize_one_partial(Omega, Amplitude, Phase, H):
    t_span = [0, H]
    t = np.arange(t_span[0], t_span[1] + 1)

    # Linear Amplitude Interpolation
    Amp_ = np.exp(interp1d(t_span, Amplitude, kind='linear', fill_value="extrapolate")(t))

    # Phase Interpolation
    T_mx = np.array([[3/H**2, -1/H], [-2/H**3, 1/H**2]])
    M_star = (Phase[0] + Omega[0] * H - Phase[1]) + (Omega[1] - Omega[0]) * H / 2
    M_star = round(1 / (2 * np.pi) * M_star)
    rowv = np.array([Phase[1] - Phase[0] - Omega[0] * H + 2 * np.pi * M_star, Omega[1] - Omega[0]])
    ab = T_mx @ rowv
    alpha, beta = ab

    Phase_ = Phase[0] + Omega[0] * t + alpha * t**2 + beta * t**3

    # Synthesize the Partial
    out = 2 * np.real(np.exp(np.log(Amp_ + 1e-10) + 1j * Phase_))

    return out

