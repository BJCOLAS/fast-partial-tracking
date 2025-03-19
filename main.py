import numpy as np
from classes.trackingPartials import trackingPartials  # Importing your class
from utilities import ddm, functions, munkres
from pprint import pprint
from scipy.optimize import linear_sum_assignment
from scipy.io.wavfile import read, write

def main():

    def F_nextpow2(i):
        N = 1
        while N < i:
            N *= 2
        return N
    

    # Parameters 
    fs = 44100 
    # Analysis window length (5ms in the snail analyser)
    #N = 2^ F_nextpow2(.025*fs)-1; 
    # Oversampling facotr  
    OverSample = 1
    #  HopSize Factor (HopSize = N / HopFactor)
    HopFactor = 4 
    # Magnitude Thresold for peak picking 
    Peak_dB  = -50
    # Polynomial order Q of the short term model : Q = 1  Frequency, Damping
    # Q = 2 Frequency Derivative, Damping Derivative
    # Q = 3 Frequency 2nd Derivative, Damping 2nd Derivative, and so on.
    Q = 2
    # Costs Function parameters 
    delta = 0.2
    zetaF = 50 # Hz
    zetaA = 15 # dB

    #fs = 44100  # Sampling rate
    N = 2^F_nextpow2(.025*fs)-1;   # Frame size
    
    # Test signal 
    # t = np.arange(P) / fs
    # f0 = 500  # Base frequency
    # f1 = 700  # Base frequency

    # InputSignalFrameLast = np.sin(2 * np.pi * f0 * t) + 0.5 * np.sin(2 * np.pi * 2*f0 * t)
    # InputSignalFrame = np.sin(2 * np.pi * f1 * t) + 0.5 * np.sin(2 * np.pi * 2*f1 * t)

    file_descriptor = "/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/sound/Sin8_A220Hz.wav"
    fs, signal = read(file_descriptor)

    # Create an instance of TrackPartials
    tracker = trackingPartials(signal, fs , N, HopFactor, OverSample, Peak_dB, Q, delta, zetaF, zetaA)

    # Setup method to initialize the variables used in the tracking
    trackingPartials.setup(tracker)

    # Tracking method 
    trackingPartials.partialTracking(tracker)


    # win = np.hanning(N)
    # winD = np.gradient(win)
    # Ndft = 2048
    # ndft = np.arange(0,Ndft,1)
    # omega = 2 * np.pi * ndft/Ndft
    
    # n_center = (N - (N % 2)) // 2
    # centering = np.exp(1j * omega * n_center)
    # ft_mat = np.exp(1j * np.outer(t - np.mean(t), omega[:Ndft//2]))

    # n = np.expand_dims(np.arange(N),axis=1) - n_center
    
    # p = functions.time_vector_computation(n,Q)
    # pD = functions.derivative_time_vector_computation(n,Q)

    # Alpha_last, num_peaks, S_last = ddm.ddm(InputSignalFrameLast, Q, 10, Peak_dB, win, winD, p, pD, centering, Ndft, ft_mat, omega)
    # Alpha, num_peaks, S = ddm.ddm(InputSignalFrame, Q, 10, Peak_dB, win, winD, p, pD, centering, Ndft, ft_mat, omega)
    
    # matrix_cost = tracker._costMatrixComputation(Alpha_last,Alpha)

    # np.save("matrix_cost",matrix_cost)

    
    # #Solution = linear_sum_assignment(matrix_cost)

    # #print(Solution)



if __name__ == "__main__":
    main()
