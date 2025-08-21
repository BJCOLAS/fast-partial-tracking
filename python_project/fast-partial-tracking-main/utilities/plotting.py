"""
Documentation
============

This code provides functions to plot partials, signals, spectrogram ... 
Some test plot are provided also but not useful for the global algorithm 

"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

def jun_plot_partials(Partials, time, fs, num_tracks, Spectro=None, Ndft=None):
    """
    Plots the partials along with a spectrogram.

    Parameters:
    Partials (list of numpy arrays): List of partials for each time frame.
    time (numpy array): Time vector.
    fs (float): Sampling frequency.
    num_tracks (int): Number of tracks to plot.
    Spectro (numpy array, optional): Spectrogram of shape (M, Ndft).
    Ndft (int, optional): FFT size, required for frequency scaling of Spectro.
    """
    M = len(time)

    track_F = -np.inf * np.ones((M, num_tracks))
    track_A = -np.inf * np.ones((M, num_tracks))

    for m in range(M):
        active = Partials[m][:, 4].astype(int)
        track_F[m, active] = fs/(2*np.pi)*Partials[m][:, 0]
        track_A[m, active] = Partials[m][:, 1]


    plt.figure(figsize=(10, 6))

    # Plot partials
    for i in range(num_tracks):
        plt.plot(time / fs, track_F[:, i], linewidth=1.5)

    # Plot spectrogram
    
    if Spectro is not None and Ndft is not None:
        # Spectrogram Normalization : 
        Spectro_norm = Spectro / np.max(Spectro)
        # Convert to dB scale (log magnitude) : 
        freqs = np.fft.rfftfreq(Ndft, d=1/fs)  # Compute correct frequency bins
        extent = [time[0] / fs, time[-1] / fs, freqs[0], freqs[-1]]  # Update extent
        plt.imshow(20 * np.log10(np.abs(np.transpose(Spectro_norm[:,:Ndft//2]))+1e-10), aspect='auto', extent = extent , origin='lower', cmap='grey')


    plt.axis('tight')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.title('Partials Plot with Spectrogram')
    plt.colorbar(label="Power (dB)")
    plt.clim([-80, 0])
    plt.ylim((0, 6000))
    plt.show()

def jun_plot_spectro(Partials, time, fs, num_tracks, Spectro=None, Ndft=None):
    """
    Plots the spectrogram only.

    Parameters:
    Partials (list of numpy arrays): List of partials for each time frame.
    time (numpy array): Time vector.
    fs (float): Sampling frequency.
    num_tracks (int): Number of tracks to plot.
    Spectro (numpy array, optional): Spectrogram of shape (M, Ndft).
    Ndft (int, optional): FFT size, required for frequency scaling of Spectro.
    """
    M = len(time)

    # track_F = -np.inf * np.ones((M, num_tracks))
    # track_A = -np.inf * np.ones((M, num_tracks))

    # for m in range(M):
    #     active = Partials[m][:, 4].astype(int)
    #     track_F[m, active] = fs/(2*np.pi)*Partials[m][:, 0]
    #     track_A[m, active] = Partials[m][:, 1]


    # plt.figure(figsize=(10, 6))

    # # Plot partials
    # for i in range(num_tracks):
    #     plt.plot(time / fs, track_F[:, i], linewidth=1.5)

    # Plot spectrogram
    
    if Spectro is not None and Ndft is not None:
        # Spectrogram Normalization : 
        Spectro_norm = Spectro / np.max(Spectro)
        # Convert to dB scale (log magnitude) : 
        freqs = np.fft.rfftfreq(Ndft, d=1/fs)  # Compute correct frequency bins
        extent = [time[0] / fs, time[-1] / fs, freqs[0], freqs[-1]]  # Update extent
        plt.imshow(20 * np.log10(np.abs(np.transpose(Spectro_norm[:,:Ndft//2]))+1e-10), aspect='auto', extent = extent , origin='lower', cmap='magma')


    plt.axis('tight')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.title('Spectrogram')
    plt.colorbar(label="Power (dB)")
    plt.clim([-80, 0])
    plt.ylim((0, 8000))
    plt.show()


def snail_plot_partials(Partials, time, fs, loudness, frequencies, num_tracks):
    """
    Plots the partials along with a spectrogram.

    Parameters:
    Partials (list of numpy arrays): List of partials for each time frame.
    time (numpy array): Time vector.
    fs (float): Sampling frequency.
    num_tracks (int): Number of tracks to plot.
    Spectro (numpy array, optional): Spectrogram of shape (M, Ndft).
    Ndft (int, optional): FFT size, required for frequency scaling of Spectro.
    """
    M = len(time)

    track_F = -np.inf * np.ones((M, num_tracks))
    track_A = -np.inf * np.ones((M, num_tracks))

    for m in range(M):
        active = Partials[m][:, 4].astype(int)
        track_F[m, active] = 1/(2*np.pi)*fs*Partials[m][:, 0]
        track_A[m, active] = Partials[m][:, 1]


    plt.figure(figsize=(10, 6))

    # Plot partials
    for i in range(num_tracks):
        plt.plot(time/fs, track_F[:, i], linewidth=1.5)
    

    # Plot spectrogram
    
    # Meshgrid-style plotting

    T, F = np.meshgrid(time/fs, frequencies)

    plt.pcolormesh(T, F, loudness.T, shading='auto', cmap='gray')



    plt.axis('tight')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.title('Partials Plot with Spectrogram')
    plt.colorbar(label="Power (dB)")
    plt.clim([-100, 0])
    plt.ylim((0, 6000))
    plt.show()
    

def snail_plot_partials_phase(Partials, time, fs, num_tracks, phiDem=None, frequencies=None):
    """
    Plots the phase partials along with a spectrogram.

    Parameters:
    Partials (list of numpy arrays): List of partials for each time frame.
    time (numpy array): Time vector.
    fs (float): Sampling frequency.
    num_tracks (int): Number of tracks to plot.
    Spectro (numpy array, optional): Spectrogram of shape (M, Ndft).
    Ndft (int, optional): FFT size, required for frequency scaling of Spectro.
    """
    M = len(time)

    track_F = -np.inf * np.ones((M, num_tracks))
    track_Phi = -np.inf * np.ones((M, num_tracks))

    for m in range(M):
        active = Partials[m][:, 3].astype(int)
        track_F[m, active] = 1/(2*np.pi)*fs*Partials[m][:, 0]
        track_Phi[m, active] = Partials[m][:, 2]


    plt.figure(figsize=(10, 6))

    # Plot partials
    for i in range(num_tracks):
        plt.plot(time/fs, track_Phi[:, i], linewidth=1.5)
    

    plt.axis('tight')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.title('Partials Plot with Spectrogram')
    plt.ylim((0, 6000))
    plt.show()

def plot_dAmp(Partials, time, fs, num_tracks) : 
    M = len(time)

    track_F = -np.inf * np.ones((M, num_tracks))
    track_dAmp = -np.inf * np.ones((M, num_tracks))

    for m in range(M):
        active = Partials[m][:, 4].astype(int)
        track_F[m, active] = 1/(2*np.pi)*fs*Partials[m][:, 0]
        track_dAmp[m, active] = Partials[m][:, 3]


    plt.figure(figsize=(10, 6))

    # Plot partials
    for i in range(num_tracks):
        plt.plot(time/fs, track_dAmp[:, i], linewidth=1.5)


    plt.axis('tight')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.title('Partials Plot with Spectrogram')
    #plt.colorbar(label="Power (dB)")
    #plt.clim([-100, 0])
    #plt.ylim((0, 5000))
    plt.show()


def plot_outputs(signal,output,residual,fs):

    t_vec = np.linspace(0,len(signal)/fs,len(signal)) 
    fig, axs = plt.subplots(2)

    fig.suptitle('Comparaison orignal et synthèse & résiduels')
    axs[0].plot(t_vec,signal)
    axs[0].plot(t_vec,signal,label='singal orginal')
    axs[0].plot(t_vec,output,'--',label='signal synthétisé')
    axs[0].legend()
    #axs[0].set_xlabel('time')
    axs[0].set_ylabel('amplitude')

    axs[1].plot(t_vec,residual,label='résiduels')
    axs[1].legend()
    axs[1].set_xlabel('temps(s)')
    axs[1].set_ylabel('amplitude')

    plt.show()


def plot_partials_comparison(Partials_depalle,Partials_snail, time_dePalle,time_snail, fs, num_tracks, phiDem=None, frequencies=None):

    M = len(time_dePalle)

    track_F = -np.inf * np.ones((M, num_tracks))
    track_Phi = -np.inf * np.ones((M, num_tracks))

    for m in range(M):
        active = Partials_depalle[m][:, 3].astype(int)
        track_F[m, active] = 1/(2*np.pi)*fs*Partials_depalle[m][:, 0]
        track_Phi[m, active] = Partials_depalle[m][:, 2]


    plt.figure(figsize=(10, 6))

    # Plot partials
    for i in range(num_tracks):
        plt.plot(time_dePalle/fs, track_Phi[:, i], linewidth=1.5,label='neri')


    M = len(time_snail)

    track_F = -np.inf * np.ones((M, num_tracks))
    track_Phi = -np.inf * np.ones((M, num_tracks))

    for m in range(M):
        active = Partials_snail[m][:, 3].astype(int)
        track_F[m, active] = 1/(2*np.pi)*fs*Partials_snail[m][:, 0]
        track_Phi[m, active] = Partials_snail[m][:, 2]


    # Plot partials
    for i in range(num_tracks):
        plt.plot(time_snail/fs, track_Phi[:, i],'--', c ='red',linewidth=1.5,label='snail')
    
    # Paramètres
    fs = 48000                   # Fréquence d'échantillonnage
    f0 = 585.9375                    # Fréquence du signal
    
    duration = 1.0              # Durée en secondes
    dft_size = 4096
    win_size = 2400
    hop_size = 240

    # Signal sinusoïdal
    t = np.arange(0, duration, 1/fs)
    x = np.sin(2 * np.pi * f0 * t)

    # Fenêtre de Hanning
    window = sig.windows.hann(win_size,False)

    # Nombre de trames
    num_frames = (len(x) - win_size) // hop_size + 1

    # STFT manuelle avec np.fft.fft
    phases = []
    bin_freqs = np.fft.fftfreq(dft_size, d=1/fs)
    bin_index = 50  # index du bin le plus proche de 480 Hz

    for i in range(num_frames):
        start = i * hop_size
        frame = x[start:start + win_size] * window
        frame_padded = np.zeros(dft_size)
        frame_padded[:win_size] = frame
        spectrum = np.fft.fft(frame_padded)
        phase = np.angle(spectrum[bin_index])
        phases.append(phase)

    # Unwrap pour suivi de phase lisse
    #phases = np.unwrap(phases)

    # Axe temporel
    time_axis = np.arange(num_frames) * hop_size / fs + win_size//2/fs
    plt.plot(time_axis, phases, c= 'green', label="Phase à 480 Hz")
    # Plot spectrogram
    
    # Meshgrid-style plotting

    #T, F = np.meshgrid(time, frequencies)

    #plt.pcolormesh(T, F, phiDem.T, shading='auto', cmap='gray')

    plt.axis('tight')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.title('Partials Plot with Spectrogram')
    #plt.legend()
    #plt.colorbar(label="Power (dB)")
    #plt.clim([-100, 0])
    #plt.ylim((0, 5000))
    plt.show()


def plot_phase_and_signal_by_track(phases,phases_mod, L, H):
    """
    Pour chaque track_id présent dans 'phases', on génère :
      - un vecteur de longueur L contenant la phase (nan si pas de données)
      - un signal synthétisé uniquement à partir de la phase
    puis on trace deux sous-graphes : phase & signal.
    """

    # pour la phase normale
    # 1. Récupérer tous les track_id uniques
    track_ids = sorted({
        p_track["track_id"]
        for tracks in phases.values()
        for p_track in tracks
    })

    # 2. Préparer les buffers par piste
    phase_buffer  = {tid: np.full(L, np.nan) for tid in track_ids}
    signal_buffer = {tid: np.zeros(L)        for tid in track_ids}

    # 3. Remplir les buffers frame par frame
    for frame, tracks in phases.items():
        start = frame * H
        for p_track in tracks:
            tid   = p_track["track_id"]
            ph    = p_track["phase"]
            n     = ph.size
            end   = min(start + n, L)

            # Stocker la phase
            phase_buffer[tid][start:end] = ph[:end - start]

            # Générer le signal à partir de la phase seule
            sig = np.real(np.exp(1j * ph[:end - start]))
            signal_buffer[tid][start:end] = sig
    
    # pour la phase mod
        # 1. Récupérer tous les track_id uniques
    track_ids = sorted({
        p_track["track_id"]
        for tracks in phases_mod.values()
        for p_track in tracks
    })

    # 2. Préparer les buffers par piste
    phase_buffer_mod = {tid: np.full(L, np.nan) for tid in track_ids}
    signal_buffer_mod = {tid: np.zeros(L)        for tid in track_ids}

    # 3. Remplir les buffers frame par frame
    for frame, tracks in phases_mod.items():
        start = frame * H
        for p_track in tracks:
            tid   = p_track["track_id"]
            ph    = p_track["phase"]
            n     = ph.size
            end   = min(start + n, L)

            # Stocker la phase
            phase_buffer_mod[tid][start:end] = ph[:end - start]

            # Générer le signal à partir de la phase seule
            sig = np.real(np.exp(1j * ph[:end - start]))
            signal_buffer_mod[tid][start:end] = sig

    # 4. Tracer pour chaque piste
    for tid in track_ids:
        fig, (ax1, ax2) = plt.subplots(
            2, 1, sharex=True, figsize=(10, 5),
            gridspec_kw={"height_ratios": [1, 1]}
        )

        x = np.arange(L)

        ax1.plot(x, phase_buffer[tid], color="C0")
        ax1.plot(x,phase_buffer_mod[tid], '--')
        ax1.set_ylabel("Phase (rad)")
        ax1.set_title(f"Track {tid} – Phase")

        ax2.plot(x, signal_buffer[tid], color="C1")
        ax2.plot(x,signal_buffer_mod[tid],'--')
        ax2.set_ylabel("Amplitude")
        ax2.set_xlabel("Échantillon")
        ax2.set_title(f"Track {tid} – Signal synthétisé (phase seule)")

        plt.tight_layout()
        plt.show()


def plot_phase_amp_and_signal_by_track(phases, amplitudes, phases_mod,amplitudes_mod, L, H):
    """
    Pour chaque track_id, trace :
      - la phase en fonction des échantillons
      - l'amplitude en fonction des échantillons
      - le signal synthétisé : Re{exp(amp + j·phase)}
    """

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

        # 1. Récupérer tous les track_id uniques
    track_ids = sorted({
        entry["track_id"]
        for frame in list(phases.values()) + list(amplitudes.values())
        for entry in frame
    })

    # 2. Préparer les buffers (phase, amplitude, signal)
    phase_buffer_mod  = {tid: np.full(L, np.nan) for tid in track_ids}
    amp_buffer_mod    = {tid: np.zeros(L)        for tid in track_ids}
    signal_buffer_mod = {tid: np.zeros(L)        for tid in track_ids}

    # 3. Remplir les buffers frame par frame
    for frame, p_tracks in phases_mod.items():
        start = frame * H
        # Créer un mapping track_id → tableau d'amplitude
        amp_map = {a_track["track_id"]: a_track["amplitudes"]
                   for a_track in amplitudes_mod.get(frame, [])}

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
            phase_buffer_mod[tid][start:end]  = ph_seg
            amp_buffer_mod[tid][start:end]    = am_seg
            signal_buffer_mod[tid][start:end] = np.real(np.exp(am_seg + 1j * ph_seg))

    # 4. Tracer pour chaque piste
    for tid in track_ids:
        x = np.arange(L)

        fig, (ax1, ax2, ax3) = plt.subplots(
            3, 1, sharex=True, figsize=(10, 7),
            gridspec_kw={"height_ratios": [1, 1, 1]}
        )

        ax1.plot(x, phase_buffer[tid], color="C0")
        ax1.plot(x, phase_buffer_mod[tid], '--')
        ax1.set_ylabel("Phase (rad)")
        ax1.set_title(f"Track {tid} – Phase")

        ax2.plot(x, amp_buffer[tid], color="C1")
        ax2.plot(x, amp_buffer_mod[tid], '--')
        ax2.set_ylabel("Amplitude")
        ax2.set_title(f"Track {tid} – Amplitude")

        ax3.plot(x, signal_buffer[tid], color="C2")
        ax3.plot(x, signal_buffer_mod[tid], '--')
        ax3.set_ylabel("Signal")
        ax3.set_xlabel("Échantillon")
        ax3.set_title(f"Track {tid} – Signal synthétisé")

        plt.tight_layout()
        plt.show()
    
    output = np.zeros(L)

    for tid in track_ids :
        output += signal_buffer_mod[tid]
    
    return output


def plot_partials_comparison(
    partials_matlab, time_matlab, fs_matlab,
    partials_python, time_python, fs_python,
    spectro=None, ndft=None, num_tracks=100,
    freq_max=6000, clim_db=(-80, 0)
):
    """
    Compare et superpose les partiels MATLAB et Python sur un même spectrogramme.

    Parameters:
    ----------
    partials_matlab : list of ndarray
        Liste de frames (M_matlab,) contenant les paramètres des partiels MATLAB.
    time_matlab : ndarray
        Vecteur temps MATLAB (en échantillons).
    fs_matlab : float
        Fréquence d'échantillonnage MATLAB.
    partials_python : list of ndarray
        Liste de frames (M_python,) contenant les paramètres des partiels Python.
    time_python : ndarray
        Vecteur temps Python (en échantillons).
    fs_python : float
        Fréquence d'échantillonnage Python.
    spectro : ndarray, optional
        Spectrogramme (T_frames x Ndft_bins) ou (Ndft_bins x T_frames).
    ndft : int, optional
        Taille FFT, nécessaire pour échelle en fréquence.
    num_tracks : int
        Nombre max de trajectoires à afficher.
    freq_max : float
        Fréquence max affichée (Hz).
    clim_db : tuple
        Limites en dB pour l'échelle du spectrogramme.
    """

    # ---- Helper pour transformer les partiels en matrices F/A ----
    def partials_to_tracks_matlab(partials, time_vec, fs, num_tracks):
        M = len(time_vec)
        track_F = -np.inf * np.ones((M, num_tracks))  # -inf = pas de partiel

        for m in range(M):
            frame = np.array(partials[m])

            # Cas frame vide → on saute, mais on garde l'index m
            if frame.size == 0:
                continue

            # Cas frame 1D → forcer en 2D
            if frame.ndim == 1:
                frame = frame[np.newaxis, :]

            active = frame[:, 3].astype(int)
            track_F[m, active] = fs / (2 * np.pi) * frame[:, 0]

        return track_F
    
    def partials_to_tracks_python(partials, time_vec, fs, num_tracks):
        M = len(time_vec)
        track_F = -np.inf * np.ones((M, num_tracks))  # -inf = pas de partiel

        for m in range(M):
            frame = np.array(partials[m])

            # Cas frame vide → on saute, mais on garde l'index m
            if frame.size == 0:
                continue

            # Cas frame 1D → forcer en 2D
            if frame.ndim == 1:
                frame = frame[np.newaxis, :]

            active = frame[:, 4].astype(int)
            track_F[m, active] = fs / (2 * np.pi) * frame[:, 0]

        return track_F

    # Convertir les deux ensembles de partiels en matrices fréquence
    track_F_matlab = partials_to_tracks_matlab(partials_matlab, time_matlab, fs_matlab, num_tracks)
    track_F_python = partials_to_tracks_python(partials_python, time_python, fs_python, num_tracks)

    plt.figure(figsize=(10, 6))

    # ---- Spectrogramme en fond ----
    if spectro is not None and ndft is not None:
        # Normalisation
        spectro_norm = spectro / np.max(np.abs(spectro))
        freqs = np.fft.rfftfreq(ndft, d=1/fs_python)  # On prend fs_python pour l'affichage
        extent = [time_python[0] / fs_python, time_python[-1] / fs_python, freqs[0], freqs[-1]]
        plt.imshow(
            20 * np.log10(np.abs(spectro_norm[:, :ndft//2].T) + 1e-10),
            aspect='auto', extent=extent, origin='lower', cmap='gray'
        )
        plt.clim(clim_db)

    # ---- Overlay Python (plein) ----
    for i in range(num_tracks):
        plt.plot(time_python / fs_python, track_F_python[:, i], linestyle='-', linewidth=1.0, color='blue', alpha=0.8)

        # ---- Overlay MATLAB (pointillés) ----
    for i in range(num_tracks):
        plt.plot(time_matlab / fs_matlab, track_F_matlab[:, i], linestyle='--', linewidth=1.0, color='red', alpha=0.8)

    plt.xlabel("Time (s)")
    plt.ylabel("Frequency (Hz)")
    plt.ylim(0, freq_max)
    plt.title("Comparaison des partiels MATLAB (-- rouge) et Python (— cyan)")
    plt.colorbar(label="Power (dB)")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

def plot_partials_side_by_side(
    partials_matlab, time_matlab, fs_matlab,
    partials_python, time_python, fs_python,
    spectro=None, ndft=None, num_tracks=100,
    freq_max=4000, clim_db=(-80, 0)
):
    """
    Affiche côte à côte les partiels MATLAB et Python sur deux sous-figures.
    """

    # ---- Helper pour transformer les partiels en matrices F/A ----
    def partials_to_tracks(partials, time_vec, fs, num_tracks, active_col):
        M = len(time_vec)
        track_F = -np.inf * np.ones((M, num_tracks))  # -inf = pas de partiel

        for m in range(M):
            frame = np.array(partials[m])
            if frame.size == 0:
                continue
            if frame.ndim == 1:
                frame = frame[np.newaxis, :]
            active = frame[:, active_col].astype(int)
            track_F[m, active] = fs / (2 * np.pi) * frame[:, 0]
        return track_F

    # Conversion
    track_F_matlab = partials_to_tracks(partials_matlab, time_matlab, fs_matlab, num_tracks, active_col=3)
    track_F_python = partials_to_tracks(partials_python, time_python, fs_python, num_tracks, active_col=4)

    # ---- Création des figures côte à côte ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    # Spectrogramme (si fourni)
    def plot_spectro(ax, time_vec, fs, track_F, title,color):
        if spectro is not None and ndft is not None:
            spectro_norm = spectro / np.max(np.abs(spectro))
            freqs = np.fft.rfftfreq(ndft, d=1/fs)
            extent = [time_vec[0] / fs, time_vec[-1] / fs, freqs[0], freqs[-1]]
            img = ax.imshow(
                20 * np.log10(np.abs(spectro_norm[:, :ndft//2].T) + 1e-10),
                aspect='auto', extent=extent, origin='lower', cmap='gray'
            )
            img.set_clim(clim_db)  # correction ici
        for i in range(num_tracks):
            ax.plot(time_vec / fs, track_F[:, i], linestyle='-', linewidth=1.0, color=color, alpha=0.8)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Frequency (Hz)")
        ax.set_ylim(0, freq_max)
        ax.set_title(title)
        ax.grid(True, linestyle='--', alpha=0.5)

    # MATLAB à gauche
    plot_spectro(axes[0], time_matlab, fs_matlab, track_F_matlab, "Partiels MATLAB",color='blue')

    # Python à droite
    plot_spectro(axes[1], time_python, fs_python, track_F_python, "Partiels Python",color='red')

    # Barre de couleur commune si spectrogramme présent
    # if spectro is not None and ndft is not None:
    #     fig.colorbar(img, ax=axes, label="Power (dB)")

    plt.tight_layout()
    plt.show()