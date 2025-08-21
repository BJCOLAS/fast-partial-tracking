from pathlib import Path 


AUDIO_FILES = {
    "oboe": Path("data/sound/oboe_impro.wav"),
    "sin_480": Path("data/sound/sine_wave_480.wav"),
    "sin_440": Path("data/sound/sine_wave_440.wav"),
    "buchla_fm": Path("data/sound/fm_buchla.wav"),
    "voice" : Path("data/sound/modulated_voice.wav"),
    "piano" : Path("data/sound/piano.wav"), 
    "vocals" : Path("/Users/colas/Documents/Programmation/Python/Snail_estimator/out2/out2.wav"), 
    "oboe_long" : Path("data/sound/oboe_long.wav"),
    "2sine" : Path("data/sound/sine_wave_480_930.wav"),
    "3sine" : Path("data/sound/sine3.wav"),
    "1sine" : Path("data/sound/sine1.wav"),
    "female_vocal":Path("data/sound/female_vocal.wav"),
    "demo_sound":Path("data/sound/demo_sound.wav"),
    "strum":Path("data/sound/strum.wav"),
    "partials_oboe":Path("data/tracking/oboe/Snail_residual_70dB_partials_oboe.wav"),
    "partials_voice":Path("ddata/tracking/3sine/snail_output_fluct_new_struct_test_partials_3sine.wav"),

}


SNAIL_ESTIMATORS = {
    "oboe": 
    {
        "loudness": "/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe/oboe_impro_loudness.mat.csv",
        "deltaPhi": "/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe/oboe_impro_deltaPhi.mat.csv",
        "phiDem": "/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe/oboe_impro_phiDem.mat.csv",
        "audio_frame": 337244,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/oboe"
    },
    "sin_480":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/sin480/sine_wave_480_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/sin480/sine_wave_480_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/sin480/sine_wave_480_phiDem.mat.csv",
        "audio_frame": 48000,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/sin_480"
    },
    "demo_sound":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/demo_sound/demo_sound_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/demo_sound/demo_sound_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/demo_sound/demo_sound_phiDem.mat.csv",
        "audio_frame":23479,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/demo_sound"
    },
    "sin_440":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/sin1/sine_wave_440_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/sin1/sine_wave_440_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/sin1/sine_wave_440_phiDem.mat.csv",
        "audio_frame":48000,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/sin_440"
    },
    "buchla_fm":
    {    
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/fm_buchla/fm_buchla_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/fm_buchla/fm_buchla_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/fm_buchla/fm_buchla_phiDem.mat.csv",
        "audio_frame":1009843,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/buchla_fm"
    },
    "voice":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/voice/modulated_voice_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/voice/modulated_voice_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/voice/modulated_voice_phiDem.mat.csv",
        "audio_frame":564480,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/voice"
    },
    "piano":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/piano/piano_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/piano/piano_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/piano/piano_phiDem.mat.csv",
        "audio_frame":176120,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/piano"

    },
    "vocals":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/out2/out2_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/out2/out2_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/out2/out2_phiDem.mat.csv",
        "audio_frame":441001,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/vocals"
    },
    "oboe_long":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe_long/oboe_long_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe_long/oboe_long_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe_long/oboe_long_phiDem.mat.csv",
        "audio_frame":1440000,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/oboe_long"
    },
    "2sine":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/2sine/sine_wave_480_930_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/2sine/sine_wave_480_930_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/2sine/sine_wave_480_930_phiDem.mat.csv",
        "audio_frame":144000,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/2sine"
    },
    "3sine":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/3sine/sine3_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/3sine/sine3_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/3sine/sine3_phiDem.mat.csv",
        "audio_frame":144000,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/3sine"
    },
        "1sine":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/1sine/sine1_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/1sine/sine1_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/1sine/sine1_phiDem.mat.csv",
        "audio_frame":144000,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/1sine"
    },
    "female_vocal":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/female_vocal/female_vocal_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/female_vocal/female_vocal_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/female_vocal/female_vocal_phiDem.mat.csv",
        "audio_frame":447233,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/female_vocal"
    },
    "strum":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/strum/strum_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/strum/strum_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/strum/strum_phiDem.mat.csv",
        "audio_frame":48248,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/strum"
    },
        "partials_oboe":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe_long/oboe_long_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe_long/oboe_long_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe_long/oboe_long_phiDem.mat.csv",
        "audio_frame":337244,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/oboe_long"
    },
        "partials_voice":
    {
        "loudness":"/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe_long/oboe_long_loudness.mat.csv",
        "deltaPhi":"/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe_long/oboe_long_deltaPhi.mat.csv",
        "phiDem":"/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe_long/oboe_long_phiDem.mat.csv",
        "audio_frame":447233,
        "save_dir":"/Users/colas/Documents/Programmation/Python/Fast_Partial_Tracking_Python/data/tracking/oboe_long"
    },
}

GENERAL_PARAMS = {
    "window_size_ms": 50,
    #"window_size_ms":42.645833333333334,
    "step_size_ms": 5,
    #"step_size_ms":5.333333333333333,
    "Peak_dB": -10,
    #     # Polynomial order Q of the short term model : Q = 1  Frequency, Damping
    #     # Q = 2 Frequency Derivative, Damping Derivative
    #     # Q = 3 Frequency 2nd Derivative, Damping 2nd Derivative, and so on.
    "Q": 1,
    "delta": 0.2,
    "zetaF": 50,
    "zetaA": 15,
    "oversample": 2
}

