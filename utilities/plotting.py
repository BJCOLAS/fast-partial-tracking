import numpy as np
import matplotlib.pyplot as plt

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
        active = Partials[m][:, 3].astype(int)
        track_F[m, active] = Partials[m][:, 0]
        track_A[m, active] = Partials[m][:, 1]


    plt.figure(figsize=(10, 6))

    # Plot partials
    for i in range(num_tracks):
        plt.plot(time / fs, track_F[:, i], linewidth=1.5)

    # Plot spectrogram
    
    if Spectro is not None and Ndft is not None:
        freqs = np.fft.rfftfreq(Ndft, d=1/fs)  # Compute correct frequency bins
        extent = [time[0] / fs, time[-1] / fs, freqs[0], freqs[-1]]  # Update extent
        plt.imshow(20 * np.log10(np.abs(np.transpose(Spectro[:,:Ndft//2])+1e-10)), aspect='auto', extent = extent , origin='lower', cmap='grey')
    

    plt.axis('tight')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.title('Partials Plot with Spectrogram')
    plt.colorbar(label="Power (dB)")
    plt.ylim((0, 5000))
    plt.show()

