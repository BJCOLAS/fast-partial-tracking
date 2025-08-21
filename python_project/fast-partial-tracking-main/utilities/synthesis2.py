import numpy as np
from utilities.functions import match_sets
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def time_derivative(time):
    N = len(time)
    for k in range(1,N):
        tm_d = (time[k]-time[k-1])/2
    return(tm_d)

def synthesize_partials_bary(Partials, time, L, Nw):
    out = np.zeros(L)
    M = len(Partials)
    num_active = 0

    T_bary = Nw / 2
    time = time + int(T_bary)

    for m in range(1,M):
        if time[m]<L :
        #plt.plot(out)
            num_active_last = num_active
            num_active = np.count_nonzero(Partials[m][:, 1])


            ## Naissance  de trajectoire
            if num_active > 0 and num_active_last == 0:
                # Fade in all
                # if m == 0:
                #     #time_m = [1, time[m]]
                #     time_m = [0, time[m]+(time[m+1]-time[m])]
                # else:
                time_m = [time[m-1], time[m]]
                H = abs(np.diff(time_m))[0]
                n_unmatched = num_active
                for idx in range(n_unmatched):
                    Omega_m = [Partials[m][idx, 0], Partials[m][idx, 0]]
                    Amp_m = [-1e10, Partials[m][idx, 1]]

                    # Phase Barycentric
                    Phase_1_bary = Partials[m][idx, 2] +  Partials[m][idx, 0] * T_bary
                    Phase_0_bary = Phase_1_bary - Partials[m][idx, 0] * H
                    Phase_m = [Phase_0_bary , Phase_1_bary]

                    # New
                    dAmp_m = [0, Partials[m][idx, 3]] # 0 useless 
                    ytemp, Freq_ = synthesize_one_partial_birth(Omega_m, Amp_m, Phase_m,dAmp_m, H)
                    out[time_m[0]:time_m[1]] += ytemp

            ## Mort de trajectoire
            elif num_active == 0 and num_active_last > 0:
                time_m = [time[m-1], time[m]]
                H = abs(np.diff(time_m))[0]
                n_unmatched_last = num_active_last
                for idx in range(n_unmatched_last):
                    Omega_m = [Partials[m-1][idx, 0], Partials[m-1][idx, 0]]
                    Amp_m = [Partials[m-1][idx, 1], -1e10]

                    # Phase Barycentric
                    Phase_0_bary = Partials[m-1][idx, 2] +  Partials[m-1][idx, 0] * T_bary
                    Phase_1_bary = Phase_0_bary + Partials[m-1][idx, 0]*H
                    Phase_m = Phase_m = [Phase_0_bary , Phase_1_bary]

                    # New 
                    dAmp_m = [Partials[m-1][idx, 3],0]

                    ytemp, Freq_ = synthesize_one_partial_death(Omega_m, Amp_m, Phase_m,dAmp_m, H)
                    out[time_m[0]:time_m[1]] += ytemp

            ## Trajectoire active 
            elif num_active > 0 and num_active_last > 0:
                matches, unmatched_last, unmatched = match_sets(Partials[m-1][:num_active_last, 4], Partials[m][:num_active, 4])
                n_matches = len(matches)
                n_unmatched_last = len(unmatched_last)
                n_unmatched = len(unmatched)

                time_m = [time[m-1], time[m]]
                H = abs(np.diff(time_m))[0]
                # Naissance dans les actives
                for idx in range(n_unmatched):
                    Omega_m = [Partials[m][unmatched[idx], 0], Partials[m][unmatched[idx], 0]]
                    Amp_m = [-1e10, Partials[m][unmatched[idx], 1]]

                    # Phase barycentric
                    Phase_1_bary = Partials[m][unmatched[idx], 2] +  Partials[m][unmatched[idx], 0] * T_bary
                    Phase_0_bary = Partials[m][unmatched[idx], 2] - Partials[m][unmatched[idx], 0] * H
                    Phase_m = [Phase_0_bary , Phase_1_bary]
                    
                    # New
                    dAmp_m = [0, Partials[m][unmatched[idx], 3]]
                    ytemp, Freq_ = synthesize_one_partial_birth(Omega_m, Amp_m, Phase_m,dAmp_m, H)
                    out[time_m[0]:time_m[1]] += ytemp

                # Mort dans les actives
                for idx in range(n_unmatched_last):
                    Omega_m = [Partials[m-1][unmatched_last[idx], 0], Partials[m-1][unmatched_last[idx], 0]]
                    Amp_m = [Partials[m-1][unmatched_last[idx], 1], -1e10]

                    # Phase barycentric
                    Phase_0_bary = Partials[m-1][unmatched_last[idx], 2] +  Partials[m-1][unmatched_last[idx], 0] * T_bary
                    Phase_1_bary = Phase_0_bary + Partials[m-1][unmatched_last[idx], 0]*H
                    Phase_m = Phase_m = [Phase_0_bary , Phase_1_bary]

                    # New
                    dAmp_m = [Partials[m-1][unmatched_last[idx], 3],0]
                    ytemp, Freq_ = synthesize_one_partial_death(Omega_m, Amp_m, Phase_m,dAmp_m, H)
                    out[time_m[0]:time_m[1]] += ytemp

                # Actives
                for idx in range(n_matches):
                    Omega_m = [Partials[m-1][matches[idx][0], 0], Partials[m][matches[idx][1], 0]]
                    Amp_m = [Partials[m-1][matches[idx][0], 1], Partials[m][matches[idx][1], 1]]

                    # Phase barycentric
                    Phase_0_bary = Partials[m-1][matches[idx][0], 2] +  Partials[m-1][matches[idx][0], 0] * T_bary
                    Phase_1_bary = Partials[m][matches[idx][1], 2] + Partials[m][matches[idx][1], 0] * T_bary
                    Phase_m = [Phase_0_bary , Phase_1_bary]


                    dAmp_m = [Partials[m-1][matches[idx][0], 3], Partials[m][matches[idx][1], 3]]
                    ytemp, Freq_ = synthesize_one_partial_active(Omega_m, Amp_m, Phase_m,dAmp_m, H)
                    out[time_m[0]:time_m[1]] += ytemp
                    
    return out

def jun_synthesize_partials(Partials, time, L):
    out = np.zeros(L)
    M = len(Partials)
    num_active = 0

    for m in range(M):
        #plt.plot(out)
        num_active_last = num_active
        num_active = np.count_nonzero(Partials[m][:, 1])


        ## Naissance  de trajectoire
        if num_active > 0 and num_active_last == 0:
            # Fade in all
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

                # New
                dAmp_m = [0, Partials[m][idx, 3]] # 0 useless 
                ytemp, Freq_ = synthesize_one_partial_birth(Omega_m, Amp_m, Phase_m,dAmp_m, H)
                #out[time_m[0]:time_m[1]] += ytemp[:-1]
                out[time_m[0]:time_m[1]] += 2*ytemp

        ## Mort de trajectoire
        elif num_active == 0 and num_active_last > 0:
            time_m = [time[m-1], time[m]]
            H = abs(np.diff(time_m))[0]
            n_unmatched_last = num_active_last
            for idx in range(n_unmatched_last):
                Omega_m = [Partials[m-1][idx, 0], Partials[m-1][idx, 0]]
                Amp_m = [Partials[m-1][idx, 1], -1e10]
                Phase_m = [Partials[m-1][idx, 2], Partials[m-1][idx, 2] + Omega_m[0] * H]

                # New 
                dAmp_m = [Partials[m-1][idx, 3],0]

                ytemp, Freq_ = synthesize_one_partial_death(Omega_m, Amp_m, Phase_m,dAmp_m, H)
                #out[time_m[0]:time_m[1]+1] += ytemp
                out[time_m[0]:time_m[1]] += 2*ytemp

        ## Trajectoire active 
        elif num_active > 0 and num_active_last > 0:
            matches, unmatched_last, unmatched = match_sets(Partials[m-1][:num_active_last, 4], Partials[m][:num_active, 4])
            n_matches = len(matches)
            n_unmatched_last = len(unmatched_last)
            n_unmatched = len(unmatched)

            time_m = [time[m-1], time[m]]
            H = abs(np.diff(time_m))[0]
            # Naissance dans les actives
            for idx in range(n_unmatched):
                Omega_m = [Partials[m][unmatched[idx], 0], Partials[m][unmatched[idx], 0]]
                Amp_m = [-1e10, Partials[m][unmatched[idx], 1]]
                Phase_m = [Partials[m][unmatched[idx], 2] - Omega_m[1] * H, Partials[m][unmatched[idx], 2]]
                # New
                dAmp_m = [0, Partials[m][unmatched[idx], 3]]
                ytemp, Freq_ = synthesize_one_partial_birth(Omega_m, Amp_m, Phase_m,dAmp_m, H)
                #out[time_m[0]:time_m[1]] += ytemp[:-1]
                out[time_m[0]:time_m[1]] += 2*ytemp
            # Mort dans les actives
            for idx in range(n_unmatched_last):
                Omega_m = [Partials[m-1][unmatched_last[idx], 0], Partials[m-1][unmatched_last[idx], 0]]
                Amp_m = [Partials[m-1][unmatched_last[idx], 1], -1e10]
                Phase_m = [Partials[m-1][unmatched_last[idx], 2], Partials[m-1][unmatched_last[idx], 2] + Omega_m[0] * H]
                # New
                dAmp_m = [Partials[m-1][unmatched_last[idx], 3],0]
                ytemp, Freq_ = synthesize_one_partial_death(Omega_m, Amp_m, Phase_m,dAmp_m, H)
                #out[time_m[0]:time_m[1]+1] += ytemp
                out[time_m[0]:time_m[1]] += 2*ytemp
            # Actives
            for idx in range(n_matches):
                Omega_m = [Partials[m-1][matches[idx][0], 0], Partials[m][matches[idx][1], 0]]
                Amp_m = [Partials[m-1][matches[idx][0], 1], Partials[m][matches[idx][1], 1]]
                Phase_m = [Partials[m-1][matches[idx][0], 2], Partials[m][matches[idx][1], 2]]
                dAmp_m = [Partials[m-1][matches[idx][0], 3], Partials[m][matches[idx][1], 3]]
                ytemp, Freq_ = synthesize_one_partial_active(Omega_m, Amp_m, Phase_m,dAmp_m, H)
                #out[time_m[0]:time_m[1]] += ytemp[:-1]
                out[time_m[0]:time_m[1]] += 2*ytemp
                
            #plt.plot(np.arange(time_m[0],time_m[1])/48000*1000,out[time_m[0]:time_m[1]])
    return out


def synthesize_one_partial_birth(Omega, Amplitude, Phase,dAmp, H):
    t_span = [0, H]
    t = np.arange(t_span[0], t_span[1])

    gamma_0 = Amplitude[0]+1j*Phase[0]
    gamma_1 = Amplitude[1]+1j*Phase[1]
    gamma_prime_0 = dAmp[0]+1j*Omega[0]
    gamma_prime_1 = dAmp[1]+1j*Omega[1] 


    #M_star = (H*gamma_prime_0+H*gamma_prime_1+2*gamma_0-2*gamma_1)/(4*np.pi*1j)
    M_star = (H*Omega[0]+H*Omega[1] + 2 *(Phase[0]-Phase[1]))/(4*np.pi)
    M_star = round(M_star)


    P_bar = np.array([[gamma_0],
                     [gamma_prime_0],
                     [6*M_star*1j*np.pi/H**2 - 2*gamma_prime_0/H - 3*gamma_0/H**2- gamma_prime_1/H + 3*gamma_1/H**2],
                     [-4*M_star*1j*np.pi/H**3 + gamma_prime_0/H**2 + 2*gamma_0/H**3 + gamma_prime_1/H**2 - 2*gamma_1/H**3]])

    Gamma_ = P_bar[0,0] + P_bar[1,0] * t + P_bar[2,0] * t**2 + P_bar[3,0] * t**3

    Phase_ = np.imag(P_bar[0,0]) + np.imag(P_bar[1,0]) * t + np.imag(P_bar[2,0]) * t**2 + np.imag(P_bar[3,0]) * t**3
    Freq_ = 48000/(2*np.pi) * (np.imag(P_bar[1,0]) + 2 *np.imag(P_bar[2,0]) * t + 3 * np.imag(P_bar[3,0]) * t**2)

    Amp_ = interp1d(t_span, np.exp(Amplitude), kind='linear')(t)
    

    out = np.real(np.exp(np.log(Amp_+1e-10) +1j*Phase_)) 

    return out, Freq_


def synthesize_one_partial_death(Omega, Amplitude, Phase,dAmp, H):
    t_span = [0, H]
    t = np.arange(t_span[0], t_span[1])

    gamma_0 = Amplitude[0]+1j*Phase[0]
    gamma_1 = Amplitude[1]+1j*Phase[1]
    gamma_prime_0 = dAmp[0]+1j*Omega[0]
    gamma_prime_1 = dAmp[1]+1j*Omega[1] 

    M_star = (H*Omega[0]+H*Omega[1] + 2 *(Phase[0]-Phase[1]))/(4*np.pi)
    M_star = round(M_star)

    P_bar = np.array([[gamma_0],
                     [gamma_prime_0],
                     [6*M_star*1j*np.pi/H**2 - 2*gamma_prime_0/H - 3*gamma_0/H**2- gamma_prime_1/H + 3*gamma_1/H**2],
                     [-4*M_star*1j*np.pi/H**3 + gamma_prime_0/H**2 + 2*gamma_0/H**3 + gamma_prime_1/H**2 - 2*gamma_1/H**3]])

    Gamma_ = P_bar[0,0] + P_bar[1,0] * t + P_bar[2,0] * t**2 + P_bar[3,0] * t**3
    Phase_ = np.imag(P_bar[0,0]) + np.imag(P_bar[1,0]) * t + np.imag(P_bar[2,0]) * t**2 + np.imag(P_bar[3,0]) * t**3
    Freq_ = 48000/(2*np.pi) * (np.imag(P_bar[1,0]) + 2 *np.imag(P_bar[2,0]) * t + 3 * np.imag(P_bar[3,0]) * t**2)

    Amp_ = interp1d(t_span, np.exp(Amplitude), kind='linear')(t)
    out = np.real(np.exp(np.log(Amp_+1e-10) +1j*Phase_))

    return out, Freq_

def synthesize_one_partial_active(Omega, Amplitude, Phase,dAmp, H):

    t_span = [0, H]
    t = np.arange(t_span[0], t_span[1])

    gamma_0 = Amplitude[0]+1j*Phase[0]
    gamma_1 = Amplitude[1]+1j*Phase[1]
    gamma_prime_0 = dAmp[0]+1j*Omega[0]
    gamma_prime_1 = dAmp[1]+1j*Omega[1] 

    #M_star = (H*gamma_prime_0+H*gamma_prime_1+2*gamma_0-2*gamma_1)/(4*np.pi*1j)
    M_star = (H*Omega[0]+H*Omega[1] + 2 *(Phase[0]-Phase[1]))/(4*np.pi)
    M_star = round(M_star)

    P_bar = np.array([[gamma_0],
                     [gamma_prime_0],
                     [6*M_star*1j*np.pi/H**2 - 2*gamma_prime_0/H - 3*gamma_0/H**2- gamma_prime_1/H + 3*gamma_1/H**2],
                     [-4*M_star*1j*np.pi/H**3 + gamma_prime_0/H**2 + 2*gamma_0/H**3 + gamma_prime_1/H**2 - 2*gamma_1/H**3]])

    Gamma_ = P_bar[0,0] + P_bar[1,0] * t + P_bar[2,0] * t**2 + P_bar[3,0] * t**3
    Phase_ = np.imag(P_bar[0,0]) + np.imag(P_bar[1,0]) * t + np.imag(P_bar[2,0]) * t**2 + np.imag(P_bar[3,0]) * t**3
    Amp_ = np.real(P_bar[0,0]) + np.real(P_bar[1,0]) * t + np.real(P_bar[2,0]) * t**2 + np.real(P_bar[3,0]) * t**3
    Freq_ = 48000/(2*np.pi) * (np.imag(P_bar[1,0]) + 2 *np.imag(P_bar[2,0]) * t + 3 * np.imag(P_bar[3,0]) * t**2)


    out = np.real(np.exp(Gamma_+1e-10))
    return out, Freq_

