import numpy as np
import matplotlib.pyplot as plt

def jun_plot_partials(Partials, time, fs, num_tracks):
    """
    Plots the partials.

    Parameters:
    Partials (list of numpy arrays): List of partials for each time frame.
    time (numpy array): Time vector.
    fs (float): Sampling frequency.
    num_tracks (int): Number of tracks to plot.
    """
    M = len(time)

    track_F = -np.inf * np.ones((M, num_tracks))
    track_A = -np.inf * np.ones((M, num_tracks))

    for m in range(M):
        active = Partials[m][:, 3].astype(int)
        track_F[m, active] = Partials[m][:, 0]
        track_A[m, active] = Partials[m][:, 1]

    plt.figure()
    for i in range(num_tracks):
        plt.plot(time / fs, fs / (2 * np.pi) * track_F[:, i], linewidth=1.5)

    plt.axis('tight')
    plt.grid(True)
    plt.box(True)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.title('Partials Plot')
    plt.show()

# Example usage:
# Partials = [np.array([[f1, a1, phi1, active1], [f2, a2, phi2, active2], ...]), ...]
# time = np.array([t1, t2, ...])
# fs = 44100  # Sampling frequency
# num_tracks = 5  # Number of tracks
# jun_plot_partials(Partials, time, fs, num_tracks)
