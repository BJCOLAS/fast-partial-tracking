import numpy as np
import matplotlib.pyplot as plt
from classes.trackingPartials import trackingPartials  # Importing your class
from classes.trackingPartialsSnail import trackingPartialsSnail
from utilities import ddm, functions, munkres, plotting, synthesis, synthesis2
from pprint import pprint
from scipy.optimize import linear_sum_assignment
from scipy.io.wavfile import read, write
from scipy.signal import stft
import pandas as pd
import time as tm
import pickle 
from scipy.signal import correlate 
from config import AUDIO_FILES, SNAIL_ESTIMATORS
from utilities.tracking_runner import run_tracking_snail, run_tracking_depalle_neri, run_tracking_snail_fluctuation,run_decoder_snail, apply_Control, capture_Fluct
from scipy.io import loadmat



def main(): 

    signal_name = "partials_voice"

    fs,signal = functions.load_audio(AUDIO_FILES[signal_name])

    RUN_SNAIL = False
    RUN_DEPALLE_NERI = True
    RUN_FLUCT = False
    WINDOW_INTERP = True
    AMP_DERIVATIVE = True 
    RUN_COMPARAISON = False

    ##################### Snail Tracking #######################

    if RUN_SNAIL == True :

        partials, tracker, audio_frame = run_tracking_snail(signal_name,fs,WINDOW_INTERP)

        time_vec = tracker.time.astype(int)

        # # 1) Charge le CSV produit par le C++
        # df = pd.read_csv("/Users/colas/Documents/Programmation/Cpp/Real-Time-Fast-Partial-Tracking/data/Partials/partials_test.csv")
        # # on s’assure des bons noms
        # df.columns = ["frame","partial","omega","amp","phase","trackID"]

        # # 2) Construit un dict { frame → ndarray (N×5) }
        # PartialsCPP = {}
        # for frame, sub in df.groupby("frame"):
        #     # Optionnel : réordonne par "partial" (l’indice local) ou laisse l’ordre CSV
        #     sub = sub.sort_values("partial")

        #     omega   = sub["omega"].to_numpy()[:,None]
        #     amp     = sub["amp"  ].to_numpy()[:,None]
        #     phase   = sub["phase"].to_numpy()[:,None]
        #     damp    = np.zeros_like(omega)          # on met des zéros pour le damp
        #     trackID = sub["trackID"].to_numpy()[:,None]

        #     PartialsCPP[int(frame)] = np.hstack([omega, amp, phase, damp, trackID])

        # plotting.snail_plot_partials(PartialsCPP, time_vec[:1397], fs, tracker.Loudness[:1397,:], tracker.frequencies, num_tracks= 10000)

        #plotting.snail_plot_partials(partials, time_vec, fs, tracker.Loudness, tracker.frequencies, num_tracks= 10000)

        output_snail = synthesis.synthesize_partials_bary(partials, time_vec, audio_frame, tracker.N)
        #output_snail_CPP = synthesis.synthesize_partials_bary(PartialsCPP, time_vec, audio_frame, tracker.N)
        #output_snail_fluctuation = synthesis.synthesize_partials_bary_fluctuation(partials, time_vec, audio_frame, tracker.N,sigma_cent=20)

        residual = functions.compute_residual(signal,output_snail)
        #residual_CPP = functions.compute_residual(signal,output_snail_CPP)
        #plotting.plot_outputs(signal,output_snail,residual,fs)

        

        functions.write_audio(output_snail,fs,config_name=signal_name,method="Snail_output_50dB", save_dir = (SNAIL_ESTIMATORS[signal_name])["save_dir"])
        functions.write_audio(residual,fs,config_name=signal_name,method="Snail_residual_70dB", save_dir = (SNAIL_ESTIMATORS[signal_name])["save_dir"])


    ################### Soa Depalle and Neri tracking ################

    if RUN_DEPALLE_NERI == True :

        partials, spectro, tracker = run_tracking_depalle_neri(signal_name)

        ## With amp derivative on synthesis
        #output_neri_damp = synthesis2.jun_synthesize_partials(partials,tracker.time.astype(int),tracker.L+np.sum(tracker.padding))[:tracker.L]
        #residual_damp = functions.compute_residual(signal,output_neri_damp)
        


        # # Charger ce qui vient de MATLAB
        # mat = loadmat('/Users/colas/Documents/Programmation/matlab_export_tracking.mat', squeeze_me=True, struct_as_record=False)
        # partials_matlab = mat['Partials']  # à adapter si structure
        # t_matlab = mat['t_sec'] * fs  # retransformer en échantillons si nécessaire
        # fs_matlab = float(mat['fs'])

        # # Python : tu as déjà partials, tracker.time, fs, spectro
        # plotting.plot_partials_side_by_side(
        #     partials_matlab=partials_matlab,
        #     time_matlab=t_matlab,
        #     fs_matlab=fs_matlab,
        #     partials_python=partials,
        #     time_python=tracker.time,
        #     fs_python=fs,
        #     spectro=spectro,
        #     ndft=tracker.Ndft,
        #     num_tracks=1000
        # )

        plotting.jun_plot_spectro(partials, tracker.time,fs,num_tracks=1000,Spectro = spectro, Ndft = tracker.Ndft)

        #plotting.jun_plot_partials(partials, tracker.time,fs,num_tracks=10000,Spectro = spectro, Ndft = tracker.Ndft)

        #plotting.plot_outputs(signal,output_neri_damp,residual_damp,fs)

        #functions.write_audio(output_neri_damp,fs,config_name=signal_name,method="Depalle_Neri_output", save_dir = (SNAIL_ESTIMATORS[signal_name])["save_dir"])
        #functions.write_audio(residual_damp,fs,config_name=signal_name,method="Depalle_Neri_residual", save_dir = (SNAIL_ESTIMATORS[signal_name])["save_dir"])
    
        ## Without amp derivative on synthesis

        #output_neri_amp = synthesis.jun_synthesize_partials(partials,tracker.time.astype(int),tracker.L+np.sum(tracker.padding))[:tracker.L]
        #residual_amp = functions.compute_residual(signal,output_neri_amp)

        #plotting.plot_outputs(signal,output_neri_amp,residual_amp,fs)

        #functions.write_audio(output_neri_amp,fs,config_name=signal_name,method="Depalle_Neri_output_amp", save_dir = (SNAIL_ESTIMATORS[signal_name])["save_dir"])
        #functions.write_audio(residual_amp,fs,config_name=signal_name,method="Depalle_Neri_residual_amp", save_dir = (SNAIL_ESTIMATORS[signal_name])["save_dir"])


     ################### Fluctuation tracking ################
    if RUN_FLUCT == True :
        

        ## Decoder without modification
        
        #phases, amplitudes, audio_frame,freq, decoder = run_decoder_snail(signal_name,fs,WINDOW_INTERP)
        

        # Decoder + Timbre control 

       # phases_mod , amplitudes_mod,audio_frame,freq, decoder = apply_Control(signal_name,fs,WINDOW_INTERP)
        
        slow_fluct, fast_fluct = capture_Fluct(signal_name,fs,WINDOW_INTERP)

        # Synthesis

        #output_fluct = synthesis.additive_synthesis_2(phases_mod,amplitudes_mod,audio_frame,decoder.tracker.Hopsize)

        #output_fluct = synthesis.additive_synthesis_2(phases,amplitudes,audio_frame,decoder.tracker.Hopsize)
        
        #partials, tracker, audio_frame = run_tracking_snail_fluctuation(signal_name,fs,WINDOW_INTERP)

        #time_vec = decoder.tracker.time.astype(int)


        #plotting.snail_plot_partials(partials, time_vec, fs, tracker.Loudness, tracker.frequencies, num_tracks= 10000)

        #output_snail = synthesis.synthesize_partials_bary_fluctuation(partials, time_vec, audio_frame, tracker.N,0.01)

        #residual = functions.compute_residual(signal,output_snail)


        #functions.write_audio(output_fluct,fs,config_name=signal_name,method="snail_output_fluct_new_struct_test", save_dir = (SNAIL_ESTIMATORS[signal_name])["save_dir"])
        #functions.write_audio(residual,fs,config_name=signal_name,method="Snail_residual_fluct", save_dir = (SNAIL_ESTIMATORS[signal_name])["save_dir"])

    if RUN_COMPARAISON == True :
        

        partials_snail, tracker, audio_frame = run_tracking_snail(signal_name,fs,WINDOW_INTERP)

        time_vec = tracker.time.astype(int)

        partials, spectro, tracker = run_tracking_depalle_neri(signal_name)

        output_neri_amp = synthesis.jun_synthesize_partials(partials,tracker.time.astype(int),tracker.L+np.sum(tracker.padding))[:tracker.L]
        residual_amp = functions.compute_residual(signal,output_neri_amp)


        output_snail = synthesis.synthesize_partials_bary(partials_snail, time_vec, audio_frame, tracker.N)


        residual = functions.compute_residual(signal,output_snail)

        snr_depalle = functions.snr(output_neri_amp, residual_amp)

        snr_snail = functions.snr(output_snail,residual)


        plotting.plot_outputs(signal,output_snail,residual,fs)


        plotting.plot_outputs(signal,output_neri_amp,residual_amp,fs)

        stats_snail = functions.analyze_partials_from_dict(partials_snail)
        stats_neri = functions.analyze_partials_from_dict(partials)

        print(stats_snail)
        print(stats_neri)

if __name__ == "__main__":
    main()
  