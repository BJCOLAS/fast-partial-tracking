import numpy as np
from utilities.functions import match_sets
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# def jun_synthesize_partials(Partials, time, L,barycentric,Nw):
#     out = np.zeros(L)
#     M = len(Partials)
#     num_active = 0
#     T_bary = Nw / 2
#     for m in range(M): 
#         if time[m]<L :     
#             #plt.plot(out)
#             num_active_last = num_active
#             num_active = np.count_nonzero(Partials[m][:, 1])

#             #birth of partials
#             if num_active > 0 and num_active_last == 0:
#                 # Fade in all
#                 if m == 0:
#                     time_m = [0, time[m]]
#                     #time_m = [0, N]
#                 else:
#                     time_m = [time[m-1], time[m]]
#                 H = abs(np.diff(time_m))[0]
#                 n_unmatched = num_active
#                 for idx in range(n_unmatched):
#                     Omega_m = [Partials[m][idx, 0], Partials[m][idx, 0]]
#                     Amp_m = [-1e10, Partials[m][idx, 1]]
#                     #Phase_m = [Partials[m][idx, 2] + Omega_m[0]*time_m[0]- Omega_m[1] * H ,Partials[m][idx, 2]+Omega_m[1]*time_m[1]]
#                     Phase_m = [Partials[m][idx, 2] - Omega_m[1] * H ,Partials[m][idx, 2]]
#                     ytemp,Freq_, M_star_pi, M_star_int, diff= jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H, barycentric,Nw,time_m)
#                     out[time_m[0]:time_m[1]] += ytemp

#             elif num_active == 0 and num_active_last > 0:
#                 time_m = [time[m-1], time[m]]
#                 H = abs(np.diff(time_m))[0]
#                 n_unmatched_last = num_active_last
#                 for idx in range(n_unmatched_last):
#                     Omega_m = [Partials[m-1][idx, 0], Partials[m-1][idx, 0]]
#                     Amp_m = [Partials[m-1][idx, 1], -1e10]
#                     #Phase_m = [Partials[m-1][idx, 2]+Omega_m[0]*time_m[0], Partials[m-1][idx, 2] + Omega_m[0] * H + Omega_m[1]*time_m[1]]
#                     Phase_m = [Partials[m-1][idx, 2], Partials[m-1][idx, 2] + Omega_m[0] * H]
#                     ytemp,Freq_, M_star_pi, M_star_int, diff = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H, barycentric,Nw,time_m)
#                     out[time_m[0]:time_m[1]] += ytemp

#             elif num_active > 0 and num_active_last > 0:
#                 matches, unmatched_last, unmatched = match_sets(Partials[m-1][:num_active_last, 3], Partials[m][:num_active, 3])
#                 n_matches = len(matches)
#                 n_unmatched_last = len(unmatched_last)
#                 n_unmatched = len(unmatched)

#                 time_m = [time[m-1], time[m]]
#                 H = abs(np.diff(time_m))[0]

#                 for idx in range(n_unmatched):
#                     Omega_m = [Partials[m][unmatched[idx], 0], Partials[m][unmatched[idx], 0]]
#                     Amp_m = [-1e10, Partials[m][unmatched[idx], 1]]
#                     #Phase_m = [Partials[m][unmatched[idx], 2] - Omega_m[1] * H +Omega_m[0]*time_m[0], Partials[m][unmatched[idx], 2] +Omega_m[1]*time_m[1]]
#                     Phase_m = [Partials[m][unmatched[idx], 2] - Omega_m[1] * H, Partials[m][unmatched[idx], 2]]
#                     ytemp,Freq_, M_star_pi, M_star_int, diff = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H, barycentric,Nw,time_m)
#                     out[time_m[0]:time_m[1]] += ytemp

#                 for idx in range(n_unmatched_last):
#                     Omega_m = [Partials[m-1][unmatched_last[idx], 0], Partials[m-1][unmatched_last[idx], 0]]
#                     Amp_m = [Partials[m-1][unmatched_last[idx], 1], -1e10]
#                     #Phase_m = [Partials[m-1][unmatched_last[idx], 2]+Omega_m[0]*time_m[0], Partials[m-1][unmatched_last[idx], 2] + Omega_m[0] * H +Omega_m[1]*time_m[1]]
#                     Phase_m = [Partials[m-1][unmatched_last[idx], 2], Partials[m-1][unmatched_last[idx], 2] + Omega_m[0] * H]
#                     ytemp,Freq_, M_star_pi, M_star_int, diff = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H, barycentric,Nw,time_m)
#                     out[time_m[0]:time_m[1]] += ytemp

#                 for idx in range(n_matches):
#                     Omega_m = [Partials[m-1][matches[idx][0], 0], Partials[m][matches[idx][1], 0]]
#                     Amp_m = [Partials[m-1][matches[idx][0], 1], Partials[m][matches[idx][1], 1]]
#                     #Phase_m = [Partials[m-1][matches[idx][0], 2]+Omega_m[0]*time_m[0], Partials[m][matches[idx][1], 2]+Omega_m[1]*time_m[1]]
#                     Phase_m = [Partials[m-1][matches[idx][0], 2], Partials[m][matches[idx][1], 2]]
#                     ytemp,Freq_, M_star_pi, M_star_int, diff = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H, barycentric,Nw,time_m)
#                     out[time_m[0]:time_m[1]] += ytemp
        

#                 #plt.plot(np.arange(time_m[0],time_m[1])/48000*1000,out[time_m[0]:time_m[1]])
#                 #plt.ylim(439,441)

#                 #plt.scatter(time_m[1]/48000*1000,diff,marker = 'o',color = 'red',label='Phase Difference')
#                 #plt.scatter(time_m[1]/48000*1000,M_star_pi,marker = 'o',color = 'blue',label='M_star_pi')
#                 #plt.scatter(time_m[1]/48000*1000,M_star_int,marker = '+', color = 'blue',label='M_star_int')
#                 #plt.legend()

#     return out

def synthesize_partials_bary(Partials, time, L,Nw):
    
    out = np.zeros(L)
    M = len(Partials)
    num_active = 0
    T_bary = Nw / 2
    time = time + int(T_bary)
    for m in range(1,M): 
        if m == 1360 :
            a = time[m]
        if time[m]<L :     
            num_active_last = num_active
            num_active = np.count_nonzero(Partials[m][:, 1])

            #birth of all partials

            if num_active > 0 and num_active_last == 0:
                
                time_m = [time[m-1], time[m]]
                H = abs(np.diff(time_m))[0]
                n_unmatched = num_active

                for idx in range(n_unmatched):

                    Omega_m = [Partials[m][idx, 0], Partials[m][idx, 0]]
                    Amp_m = [-1e10, Partials[m][idx, 1]]

                    Phase_1_bary = Partials[m][idx, 2] +  Partials[m][idx, 0] * T_bary
                    Phase_0_bary = Phase_1_bary - Partials[m][idx, 0] * H
                    Phase_m = [Phase_0_bary , Phase_1_bary]

                    ytemp= jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                    out[time_m[0]:time_m[1]] += ytemp

            # death of all partials
            elif num_active == 0 and num_active_last > 0:
                time_m = [time[m-1], time[m]]
                H = abs(np.diff(time_m))[0]
                n_unmatched_last = num_active_last
                for idx in range(n_unmatched_last):
                    Omega_m = [Partials[m-1][idx, 0], Partials[m-1][idx, 0]]
                    Amp_m = [Partials[m-1][idx, 1], -1e10]
                    
                    Phase_0_bary = Partials[m-1][idx, 2] +  Partials[m-1][idx, 0] * T_bary
                    Phase_1_bary = Phase_0_bary + Partials[m-1][idx, 0]*H
                    Phase_m = Phase_m = [Phase_0_bary , Phase_1_bary]

                    ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                    out[time_m[0]:time_m[1]] += ytemp

            # general
            elif num_active > 0 and num_active_last > 0:
                matches, unmatched_last, unmatched = match_sets(Partials[m-1][:num_active_last, 4], Partials[m][:num_active, 4])
                n_matches = len(matches)
                n_unmatched_last = len(unmatched_last)
                n_unmatched = len(unmatched)

                time_m = [time[m-1], time[m]]
                H = abs(np.diff(time_m))[0]

                for idx in range(n_unmatched):
                    Omega_m = [Partials[m][unmatched[idx], 0], Partials[m][unmatched[idx], 0]]
                    Amp_m = [-1e10, Partials[m][unmatched[idx], 1]]
                    
                    Phase_1_bary = Partials[m][unmatched[idx], 2] +  Partials[m][unmatched[idx], 0] * T_bary
                    Phase_0_bary = Partials[m][unmatched[idx], 2] - Partials[m][unmatched[idx], 0] * H
                    Phase_m = [Phase_0_bary , Phase_1_bary]
                    ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                    out[time_m[0]:time_m[1]] += ytemp

                for idx in range(n_unmatched_last):
                    Omega_m = [Partials[m-1][unmatched_last[idx], 0], Partials[m-1][unmatched_last[idx], 0]]
                    Amp_m = [Partials[m-1][unmatched_last[idx], 1], -1e10]
                    
                    Phase_0_bary = Partials[m-1][unmatched_last[idx], 2] +  Partials[m-1][unmatched_last[idx], 0] * T_bary
                    Phase_1_bary = Phase_0_bary + Partials[m-1][unmatched_last[idx], 0]*H
                    Phase_m = Phase_m = [Phase_0_bary , Phase_1_bary]
                    ytemp= jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                    out[time_m[0]:time_m[1]] += ytemp

                for idx in range(n_matches):
                    Omega_m = [Partials[m-1][matches[idx][0], 0], Partials[m][matches[idx][1], 0]]
                    Amp_m = [Partials[m-1][matches[idx][0], 1], Partials[m][matches[idx][1], 1]]

                    Phase_0_bary = Partials[m-1][matches[idx][0], 2] +  Partials[m-1][matches[idx][0], 0] * T_bary
                    Phase_1_bary = Partials[m][matches[idx][1], 2] + Partials[m][matches[idx][1], 0] * T_bary
                    Phase_m = [Phase_0_bary , Phase_1_bary]
                    ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                    out[time_m[0]:time_m[1]] += ytemp

    return out

def jun_synthesize_partials(Partials, time, L):
    
    out = np.zeros(L)
    M = len(Partials)
    num_active = 0

    for m in range(1,M):    
        num_active_last = num_active
        num_active = np.count_nonzero(Partials[m][:, 1])

        #birth of all partials

        if num_active > 0 and num_active_last == 0:
            
            if m == 0:
                #time_m = [1, time[m]]
                time_m = [0, time[m]]
            else:
                time_m = [time[m-1], time[m]]

            H = abs(np.diff(time_m))[0]
            n_unmatched = num_active

            for idx in range(n_unmatched):

                Omega_m = [Partials[m][idx, 0], Partials[m][idx, 0]]
                Amp_m = [-1e10, Partials[m][idx, 1]]

                Phase_m = [Partials[m][idx, 2] - Omega_m[1] * H, Partials[m][idx, 2]]  

                ytemp= jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                out[time_m[0]:time_m[1]] += 2*ytemp

        # death of all partials
        elif num_active == 0 and num_active_last > 0:
            time_m = [time[m-1], time[m]]
            H = abs(np.diff(time_m))[0]
            n_unmatched_last = num_active_last
            for idx in range(n_unmatched_last):
                Omega_m = [Partials[m-1][idx, 0], Partials[m-1][idx, 0]]
                Amp_m = [Partials[m-1][idx, 1], -1e10]
                

                Phase_m = [Partials[m-1][idx, 2], Partials[m-1][idx, 2] + Omega_m[0] * H]

                ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                out[time_m[0]:time_m[1]] += 2*ytemp

        # general
        elif num_active > 0 and num_active_last > 0:
            matches, unmatched_last, unmatched = match_sets(Partials[m-1][:num_active_last, 4], Partials[m][:num_active, 4])
            n_matches = len(matches)
            n_unmatched_last = len(unmatched_last)
            n_unmatched = len(unmatched)

            time_m = [time[m-1], time[m]]
            H = abs(np.diff(time_m))[0]

            for idx in range(n_unmatched):
                Omega_m = [Partials[m][unmatched[idx], 0], Partials[m][unmatched[idx], 0]]
                Amp_m = [-1e10, Partials[m][unmatched[idx], 1]]

                Phase_m = [Partials[m][unmatched[idx], 2] - Omega_m[1] * H, Partials[m][unmatched[idx], 2]]
                ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                out[time_m[0]:time_m[1]] += 2*ytemp

            for idx in range(n_unmatched_last):
                Omega_m = [Partials[m-1][unmatched_last[idx], 0], Partials[m-1][unmatched_last[idx], 0]]
                Amp_m = [Partials[m-1][unmatched_last[idx], 1], -1e10]
                
                Phase_m = [Partials[m-1][unmatched_last[idx], 2], Partials[m-1][unmatched_last[idx], 2] + Omega_m[0] * H]
                ytemp= jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                out[time_m[0]:time_m[1]] += 2*ytemp

            for idx in range(n_matches):
                Omega_m = [Partials[m-1][matches[idx][0], 0], Partials[m][matches[idx][1], 0]]
                Amp_m = [Partials[m-1][matches[idx][0], 1], Partials[m][matches[idx][1], 1]]


                Phase_m = [Partials[m-1][matches[idx][0], 2], Partials[m][matches[idx][1], 2]]

                ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H)
                out[time_m[0]:time_m[1]] += 2*ytemp

    return out

def jun_synthesize_one_partial(Omega, Amplitude, Phase, H):

    t_span = [0, H]
    t = np.arange(t_span[0], t_span[1])

    # Linear Amplitude Interpolation
    Amp_ = interp1d(t_span, np.exp(Amplitude), kind='linear')(t)

    dif = Phase[1] - ( Phase[0] + Omega[0]*H)

    # Phase Interpolation
    T_mx = np.array([[3/H**2, -1/H], [-2/H**3, 1/H**2]])
    M_star = (Phase[0] + Omega[0] * H - Phase[1]) + (Omega[1] - Omega[0]) * H / 2

    M_star_pi = M_star / (2 * np.pi)
    M_star_int = np.floor(M_star_pi+0.5)

    #M_star = round(1 / (2 * np.pi) * M_star)
    M_star = M_star_int
    rowv = np.array([Phase[1] - Phase[0] - Omega[0] * H + 2 * np.pi * M_star, Omega[1] - Omega[0]])
    ab = T_mx @ rowv
    alpha, beta = ab

    Phase_ = Phase[0] + Omega[0] * t + alpha * t**2 + beta * t**3

    Freq_ = 48000/(2*np.pi) * (Omega[0] + 2 * alpha * t + 3 * beta * t**2)

    # Synthesize the Partial
    out = np.real(np.exp(np.log(Amp_ + 1e-10) + 1j * Phase_))

    return out

def synthesize_partials_bary_fluctuation(Partials, time, L,Nw,sigma_cent):
    
    out = np.zeros(L)
    M = len(Partials)
    num_active = 0
    T_bary = Nw / 2
    time = time + int(T_bary)

    ## fluctutation ##

    # Longueur totale en samples

    total_length = L
    # Jitter global en omega normalisé (à adapter selon ton domaine)
    sigma = 1.5  # très petit pour éviter le chaos
    t_ms = 100
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
    tvec_ech = np.arange(total_length)
    global_jitter = interp_fn(tvec_ech)
    
    for m in range(1,M): 
        if time[m]<L :     
            num_active_last = num_active
            num_active = np.count_nonzero(Partials[m][:, 1])

            #birth of all partials

            if num_active > 0 and num_active_last == 0:
                
                time_m = [time[m-1], time[m]]
                H = abs(np.diff(time_m))[0]
                n_unmatched = num_active

                for idx in range(n_unmatched):

                    Omega_m = [Partials[m][idx, 0], Partials[m][idx, 0]]
                    Amp_m = [-1e10, Partials[m][idx, 1]]

                    Phase_1_bary = Partials[m][idx, 2] +  Partials[m][idx, 0] * T_bary
                    Phase_0_bary = Phase_1_bary - Partials[m][idx, 0] * H
                    Phase_m = [Phase_0_bary , Phase_1_bary]

                    ytemp= jun_synthesize_one_partial_fluctuation(Omega_m, Amp_m, Phase_m, H, sigma_cent, fs=48000, mode='sample')
                    out[time_m[0]:time_m[1]] += ytemp

            # death of all partials
            elif num_active == 0 and num_active_last > 0:
                time_m = [time[m-1], time[m]]
                H = abs(np.diff(time_m))[0]
                n_unmatched_last = num_active_last
                for idx in range(n_unmatched_last):
                    Omega_m = [Partials[m-1][idx, 0], Partials[m-1][idx, 0]]
                    Amp_m = [Partials[m-1][idx, 1], -1e10]
                    
                    Phase_0_bary = Partials[m-1][idx, 2] +  Partials[m-1][idx, 0] * T_bary
                    Phase_1_bary = Phase_0_bary + Partials[m-1][idx, 0]*H
                    Phase_m = Phase_m = [Phase_0_bary , Phase_1_bary]

                    ytemp= jun_synthesize_one_partial_fluctuation(Omega_m, Amp_m, Phase_m, H, sigma_cent, fs=48000, mode='sample')
                    out[time_m[0]:time_m[1]] += ytemp

            # general
            elif num_active > 0 and num_active_last > 0:
                matches, unmatched_last, unmatched = match_sets(Partials[m-1][:num_active_last, 4], Partials[m][:num_active, 4])
                n_matches = len(matches)
                n_unmatched_last = len(unmatched_last)
                n_unmatched = len(unmatched)

                time_m = [time[m-1], time[m]]
                H = abs(np.diff(time_m))[0]

                for idx in range(n_unmatched):
                    Omega_m = [Partials[m][unmatched[idx], 0], Partials[m][unmatched[idx], 0]]
                    Amp_m = [-1e10, Partials[m][unmatched[idx], 1]]
                    
                    Phase_1_bary = Partials[m][unmatched[idx], 2] +  Partials[m][unmatched[idx], 0] * T_bary
                    Phase_0_bary = Partials[m][unmatched[idx], 2] - Partials[m][unmatched[idx], 0] * H
                    Phase_m = [Phase_0_bary , Phase_1_bary]
                    ytemp= jun_synthesize_one_partial_fluctuation(Omega_m, Amp_m, Phase_m, H, sigma_cent, fs=48000, mode='sample')
                    out[time_m[0]:time_m[1]] += ytemp

                for idx in range(n_unmatched_last):
                    Omega_m = [Partials[m-1][unmatched_last[idx], 0], Partials[m-1][unmatched_last[idx], 0]]
                    Amp_m = [Partials[m-1][unmatched_last[idx], 1], -1e10]
                    
                    Phase_0_bary = Partials[m-1][unmatched_last[idx], 2] +  Partials[m-1][unmatched_last[idx], 0] * T_bary
                    Phase_1_bary = Phase_0_bary + Partials[m-1][unmatched_last[idx], 0]*H
                    Phase_m = Phase_m = [Phase_0_bary , Phase_1_bary]
                    ytemp= jun_synthesize_one_partial_fluctuation(Omega_m, Amp_m, Phase_m, H, sigma_cent, fs=48000, mode='sample')
                    out[time_m[0]:time_m[1]] += ytemp

                for idx in range(n_matches):
                    Omega_m = [Partials[m-1][matches[idx][0], 0], Partials[m][matches[idx][1], 0]]
                    Amp_m = [Partials[m-1][matches[idx][0], 1], Partials[m][matches[idx][1], 1]]

                    Phase_0_bary = Partials[m-1][matches[idx][0], 2] +  Partials[m-1][matches[idx][0], 0] * T_bary
                    Phase_1_bary = Partials[m][matches[idx][1], 2] + Partials[m][matches[idx][1], 0] * T_bary
                    Phase_m = [Phase_0_bary , Phase_1_bary]
                    ytemp= jun_synthesize_one_partial_fluctuation(Omega_m, Amp_m, Phase_m, H, sigma_cent, fs=48000, mode='sample')
                    out[time_m[0]:time_m[1]] += ytemp
                
                # plt.plot(out[0:H*m])
                # plt.show()

    return out


def jun_synthesize_one_partial_fluctuation(Omega, Amplitude, Phase, H, sigma_cent, fs=48000, mode='sample'):

    t_span = [0, H]
    t = np.arange(t_span[0], t_span[1])

    # Linear Amplitude Interpolation
    Amp_ = interp1d(t_span, np.exp(Amplitude), kind='linear')(t)

    dif = Phase[1] - ( Phase[0] + Omega[0]*H)

    # Phase Interpolation
    T_mx = np.array([[3/H**2, -1/H], [-2/H**3, 1/H**2]])
    M_star = (Phase[0] + Omega[0] * H - Phase[1]) + (Omega[1] - Omega[0]) * H / 2

    M_star_pi = M_star / (2 * np.pi)
    M_star_int = np.floor(M_star_pi+0.5)

    #M_star = round(1 / (2 * np.pi) * M_star)
    M_star = M_star_int
    rowv = np.array([Phase[1] - Phase[0] - Omega[0] * H + 2 * np.pi * M_star, Omega[1] - Omega[0]])
    ab = T_mx @ rowv
    alpha, beta = ab

    Phase_ = (Phase[0] + Omega[0] * t + alpha * t**2 + beta * t**3 )

    out_test = np.real(np.exp(np.log(Amp_ + 1e-10) + 1j * Phase_))

    
    Freq_ = 48000/(2*np.pi) * (Omega[0] + 2 * alpha * t + 3 * beta * t**2)

    length = 256
    fs = 48000
    sigma = 10.0  # Écart-type de la loi normale (tu peux ajuster)
    t_ms = 100 #ms
    g_poisson = 1
    lambda_poisson = 1/(t*fs/1000)*g_poisson   # Taux de Poisson (plus bas = moins d'événements)

    # Étape 1: Générer les instants du processus de Poisson
    inter_arrivals = np.random.exponential(scale=1/lambda_poisson, size=length)
    instants = np.cumsum(inter_arrivals)
    instants = instants[instants < length]  # On garde que dans la fenêtre de 0 à 255
    instants = np.insert(instants, 0, 0)  # Ajoute 0 au début si besoin
    instants = np.append(instants, length - 1)  # Fin aussi

    # Étape 2: Générer les valeurs gaussiennes (centrées sur une fréquence cible, disons 440 Hz)
    base_freq = 0
    freq_values = np.random.normal(loc=base_freq, scale=sigma, size=len(instants))

    

    Freq_fluct = Freq_ + freq_values
    Omega_fluct = 2 * np.pi * Freq_fluct /fs
    Phase_fluct = np.cumsum(Omega_fluct)

    # Phase_ = Phase[0] + Phase_fluct

    out = np.real(np.exp(np.log(Amp_ + 1e-10) + 1j * Phase_))

    # plt.plot(Freq_)
    # plt.plot(Freq_fluct)
    # #plt.plot(out_test,label='original')
    # #plt.plot(out,label='fluct')
    # #plt.plot(lfo)
    # plt.legend()
    # plt.show()

    return out

def additive_synthesis(phases, amplitudes, L,H):
    """
    Calcule le signal synthétisé en sommant exp(amp + j * phase) pour chaque track, par frame.
    Retourne la partie réelle du signal total.
    
    Args:
        phases (dict): dictionnaire {frame: [ {"phase": np.array, "track_id": id }, ... ]}
        amplitudes (dict): même structure que phases mais avec "amp" à la place de "phase"
        L (int): longueur totale du signal à générer

    Returns:
        np.ndarray: signal synthétisé (partie réelle uniquement)
    """
    output = np.zeros(L)
    phase_ = np.zeros(L)
    for frame in phases:
        # Vérifier la cohérence avec les amplitudes
        if frame not in amplitudes:
            continue
        
        tracks_phase = phases[frame]
        tracks_amp = amplitudes[frame]

        for p_track, a_track in zip(tracks_phase, tracks_amp):
            phase = p_track["phase"]
            amp = a_track["amplitudes"]

            
            # Condition: on ne calcule le signal que si phase et amp sont non vides
            if phase is None or amp is None:
                continue
            if phase.size == 0 or amp.size == 0:
                continue

            # On assure également que le signal et la phase sont de même longueur
            n = min(phase.size, amp.size)
            phase = phase[:n]
            amp   = amp[:n]

            # Calcul du signal complexe et ajout de sa partie réelle
            signal = np.exp(amp + 1j * phase)

            start = frame * H
            end   = min(start + n, L)
            output[start:end] += np.real(signal[: end - start])

            phase_[start:end] = phase



    return output

def additive_synthesis_2(phases,amplitudes,L,H):

    output = np.zeros(L)

    # 1. Récupérer tous les track_id uniques
    track_ids = sorted({
        entry["track_id"]
        for frame in list(phases.values()) + list(amplitudes.values())
        for entry in frame
    })

    # 2. Préparer les buffers (phase, amplitude, signal)
    phase_buffer  = {tid: np.full(L, np.nan) for tid in track_ids}
    amp_buffer    = {tid: np.zeros(L)        for tid in track_ids}
    signal_buffer = {tid: np.zeros(L)        for tid in track_ids}

    # 3. Remplir les buffers frame par frame
    for frame, p_tracks in phases.items():
        start = frame * H

        if start<L:
            # Créer un mapping track_id → tableau d'amplitude
            amp_map = {a_track["track_id"]: a_track["amplitudes"]
                    for a_track in amplitudes.get(frame, [])}

            for p_track in p_tracks:
                tid = p_track["track_id"]
                ph  = p_track["phase"]
                am  = amp_map.get(tid)

                # Vérifications
                if ph is None or am is None or ph.size == 0 or am.size == 0:
                    continue

                # Tronquer à la taille minimale
                n   = min(ph.size, am.size)
                end = min(start + n, L)
                ph_seg = ph[:n]
                am_seg = am[:n]

                # Remplissage des buffers
                phase_buffer[tid][start:end]  = ph_seg
                amp_buffer[tid][start:end]    = am_seg
                signal_buffer[tid][start:end] = np.real(np.exp(am_seg + 1j * ph_seg))


    for tid in track_ids :
        output += signal_buffer[tid]
    return output 

