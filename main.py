import numpy as np
import matplotlib.pyplot as plt
from classes.trackingPartials import trackingPartials  # Importing your class
from utilities import ddm, functions, munkres, plotting
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
    #fs = 44100 
    # Analysis window length (5ms in the snail analyser)
    #N = 2^ F_nextpow2(.025*fs)-1; 

    #file_descriptor = "/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/sound/Sin8_A220Hz.wav"
    file_descriptor = "/Users/colas/Documents/Programmation/Matlab_Projects/Fast-Partial-Tracking/demo_sound.wav"
    #file_descriptor = "/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/sound/sine_wave_variating.wav"
    #file_descriptor = "/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/sound/sine_wave_1_440.wav"


    fs, signal = read(file_descriptor)

    # Oversampling facotr  
    OverSample = 2
    #  HopSize Factor (HopSize = N / HopFactor)
    HopFactor = 4 
    # Magnitude Thresold for peak picking 
    Peak_dB  = -40
    # Polynomial order Q of the short term model : Q = 1  Frequency, Damping
    # Q = 2 Frequency Derivative, Damping Derivative
    # Q = 3 Frequency 2nd Derivative, Damping 2nd Derivative, and so on.
    Q = 2
    # Costs Function parameters 
    delta = 0.2
    zetaF = 50 # Hz
    zetaA = 15 # dB

    #fs = 44100  # Sampling rate
    N = 2^F_nextpow2(.025*fs) - 1;   # Frame size



    # Create an instance of TrackPartials
    tracker = trackingPartials(signal, fs , N, HopFactor, OverSample, Peak_dB, Q, delta, zetaF, zetaA)

    # # Setup method to initialize the variables used in the tracking
    trackingPartials.setup(tracker)
    # Tracking method 
    Spectro, Partials = trackingPartials.partialTracking(tracker)

    # np.save("Spectro_Sin8",Spectro)
    # np.save("Tracker_time",tracker.time)

    plotting.jun_plot_partials(Partials,tracker.time,fs,1000,Spectro,tracker.Ndft)


if __name__ == "__main__":
    main()
