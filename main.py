import numpy as np
import matplotlib.pyplot as plt
from classes.trackingPartials import trackingPartials  # Importing your class
from utilities import ddm, functions, munkres, plotting, synthesis
from pprint import pprint
from scipy.optimize import linear_sum_assignment
from scipy.io.wavfile import read, write
import pandas as pd
import time
import pickle 

def main():

    #file_descriptor = "/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/sound/Sin8_A220Hz.wav"
    file_descriptor = "/Users/colas/Documents/Programmation/Matlab_Projects/Fast-Partial-Tracking/demo_sound.wav"
    #file_descriptor = "/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/sound/sine_wave_variating.wav"
    #file_descriptor = "/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/sound/sine_wave_1_440.wav"
    #file_descriptor = "/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/sound/jun_piazolla_original.wav"

    # open signal
    fs, signal = read(file_descriptor)

    # Convert to float32 (normalize between -1 and 1)
    signal = signal.astype('float32') / 32768.0

    # Oversampling factor  
    OverSample = 2
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
    # Analysis window length : 

    N = 2 ** (int(np.ceil(np.log2(0.025 * fs))))-1 



    ############################# Tracking method ##########################
    
    # Create an instance of TrackPartials
    tracker = trackingPartials(signal, fs , N, HopFactor, OverSample, Peak_dB, Q, delta, zetaF, zetaA)

    # # Setup method to initialize the variables used in the tracking
    trackingPartials.setup(tracker)

    # # t_start = time.time()

    Spectro, Partials = trackingPartials.partialTracking(tracker)

    # # t_stop = time.time()-t_start
    # # print("time enabled to perform the tracking :", t_stop)

    # with open("partials.pkl", "wb") as f:
    #     pickle.dump(Partials, f)

    # np.save("spectro.npy",Spectro)

    ############################# Plotting the results ##########################
    
    # Spectro_loaded = np.load("spectro.npy")

    # print(Spectro_loaded)
    # with open("partials.pkl", "rb") as f:
    #     Partials_loaded = pickle.load(f)

    num_tracks = 1000
    plotting.jun_plot_partials(Partials,tracker.time,fs,num_tracks,Spectro,tracker.Ndft)

    ############################# Synthesis method ##########################

    output = synthesis.jun_synthesize_partials(Partials,tracker.time.astype(int),tracker.L+np.sum(tracker.padding))
    
    fig, axs = plt.subplots(2)
    fig.suptitle('Original signal & synthesized signal')
    axs[0].plot(signal)
    axs[1].plot(output/max(output))
    plt.show()
    
    write("demo.wav",rate = fs , data = output/max(output))
    
if __name__ == "__main__":
    main()
