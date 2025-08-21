"""
Documentation
============

"""

import matplotlib.pyplot as plt
import numpy as np


def run_tracking_snail(config_id,fs,WINDOW_INTERP):
    from config import SNAIL_ESTIMATORS, GENERAL_PARAMS
    from classes.trackingPartialsSnail import trackingPartialsSnail
    from utilities.functions import load_snail_estimators, save_partials, save_spectro
    import os, pickle

    snail_cfg = SNAIL_ESTIMATORS[config_id]
    save_dir = snail_cfg["save_dir"]
    # === Try to load saved partials ===
    part_file = os.path.join(save_dir, f"snail_partials_{config_id}.pkl")

    deltaPhi, loudness, phiDem = load_snail_estimators(snail_cfg)

    N = int(fs * GENERAL_PARAMS["window_size_ms"] / 1000)
    HopSize = int(fs * GENERAL_PARAMS["step_size_ms"] / 1000)

    tracker = trackingPartialsSnail(
        deltaPhi, loudness, phiDem, fs, N, HopSize,
        GENERAL_PARAMS["Peak_dB"], GENERAL_PARAMS["Q"],
        GENERAL_PARAMS["delta"], GENERAL_PARAMS["zetaF"], GENERAL_PARAMS["zetaA"],WINDOW_INTERP
    )

    trackingPartialsSnail.setup(tracker)

    import time
    t0 = time.time()
    partials = trackingPartialsSnail.partialTracking(tracker)  

    # sauvegarde
    save_partials(partials, config_id, "snail", snail_cfg["save_dir"])
    save_spectro(tracker.Loudness,config_id,"snail",snail_cfg["save_dir"])

    print("Snail tracking time:", time.time() - t0)

    return partials, tracker, snail_cfg["audio_frame"]

def run_tracking_depalle_neri(config_id):
    from classes.trackingPartials import trackingPartials
    from config import AUDIO_FILES,SNAIL_ESTIMATORS, GENERAL_PARAMS
    from utilities.functions import load_audio
    from utilities.functions import load_snail_estimators, save_partials

    snail_cfg = SNAIL_ESTIMATORS[config_id] # to get the save directory
    fs,signal = load_audio(AUDIO_FILES[config_id])

    N = int(fs * GENERAL_PARAMS["window_size_ms"] / 1000)
    HopSize = int(fs * GENERAL_PARAMS["step_size_ms"] / 1000)
    HopFactor = int(N/HopSize) 

    N = 2047 
    #Hopsize = 256
    HopFactor = 4

    tracker = trackingPartials(signal, fs , N, HopFactor, GENERAL_PARAMS["oversample"], 
                               GENERAL_PARAMS["Peak_dB"], GENERAL_PARAMS["Q"], GENERAL_PARAMS["delta"], 
                               GENERAL_PARAMS["zetaF"], GENERAL_PARAMS["zetaA"])


    trackingPartials.setup(tracker)

    import time
    t0 = time.time()
    Spectro, Partials = trackingPartials.partialTracking(tracker)
    print("Neri tracking time:", time.time() - t0)

    # sauvegarde
    save_partials(Partials, config_id, "depalle_neri", snail_cfg["save_dir"])

    return Partials, Spectro, tracker


def run_tracking_snail_fluctuation(config_id,fs,WINDOW_INTERP):
    from config import SNAIL_ESTIMATORS, GENERAL_PARAMS
    from classes.trackedParameterSynthesis import  trackedParameterSynthesis
    from classes.trackingPartialsSnail import trackingPartialsSnail
    from utilities.functions import load_snail_estimators, save_partials, save_spectro
    import os, pickle

    snail_cfg = SNAIL_ESTIMATORS[config_id]
    save_dir = snail_cfg["save_dir"]
    # === Try to load saved partials ===
    part_file = os.path.join(save_dir, f"snail_partials_{config_id}.pkl")

    deltaPhi, loudness, phiDem = load_snail_estimators(snail_cfg)

    N = int(fs * GENERAL_PARAMS["window_size_ms"] / 1000)
    HopSize = int(fs * GENERAL_PARAMS["step_size_ms"] / 1000)

    tracker = trackingPartialsSnail(
        deltaPhi, loudness, phiDem, fs, N, HopSize,
        GENERAL_PARAMS["Peak_dB"], GENERAL_PARAMS["Q"],
        GENERAL_PARAMS["delta"], GENERAL_PARAMS["zetaF"], GENERAL_PARAMS["zetaA"],WINDOW_INTERP
    )

    trackingPartialsSnail.setup(tracker)

    import time
    t0 = time.time()
    partials = trackingPartialsSnail.partialTracking(tracker)  
    
    # transformation

    parameterSynthesis = trackedParameterSynthesis(tracker,partials)
    print("Snail tracking time:", time.time() - t0)
    parameters =  parameterSynthesis.transform()

    return parameters, tracker, snail_cfg["audio_frame"]


def run_decoder_snail(config_id,fs,WINDOW_INTERP):
    from classes.synthesisDecoder import synthesisDecoder

    parameters , tracker , audio_frame = run_tracking_snail_fluctuation(config_id,fs,WINDOW_INTERP)
    ## Transformation de la structure de donnée
    decoder = synthesisDecoder(parameters,tracker)
    phases,amplitudes, out_complementary_param = decoder.compute_output()

    return phases , amplitudes,audio_frame,out_complementary_param, decoder


def apply_Control(config_id,fs,WINDOW_INTERP):
    from classes.timbreControler import timbreController

    phases, amplitudes, audio_frame,freq, decoder = run_decoder_snail(config_id,fs,WINDOW_INTERP)

    lfo_control = timbreController(amplitudes,phases,freq)

    #phases = lfo_control.apply_LFO_test(fs)

    #phases = lfo_control.apply_LFO_global(lfo_type="sin",rate = 0.8 , sigma_hz=30,fs=fs)


    lfo_params = {
    1: (0.3,  5.0),   # track 1 : LFO 0.3 Hz, ±5 Hz
    2: (1.2, 10.0)    # track 2 : LFO 1.2 Hz, ±10 Hz
        }

    #phases = lfo_control.apply_LFO_2(lfo_type="sin",rate = 0.5 , sigma_hz=10,total_length=audio_frame,fs=fs)

    ## mass effect
    sigmacent = 35
    t_ms = 100
    poisson_gain = 1
    total_length = audio_frame 

    phases =  lfo_control.apply_LFO_jitter_per_track(
                                sigmacent,
                                t_ms,
                                poisson_gain,
                                fs,
                                total_length)
    

    #sigmapercent = 0.05
    #t_ms = 80
    
    #amplitudes = lfo_control.apply_amplitude_jitter_per_track(sigmapercent,t_ms,poisson_gain,fs,total_length)

    return phases , amplitudes,audio_frame,freq, decoder


def capture_Fluct(config_id,fs,WINDOW_INTERP):
    from classes.fluctuationCapture import fluctuationCapture

    phases, amplitudes, audio_frame,freq, decoder = run_decoder_snail(config_id,fs,WINDOW_INTERP)


    capture = fluctuationCapture(freq)
    capture.applyFiltering(H=decoder.tracker.Hopsize,L=audio_frame)
    