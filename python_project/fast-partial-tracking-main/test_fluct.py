import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import matplotlib.cm as cm
import pickle
import matplotlib

# 🔧 Paramètres
fs = 48000
hopsize = 0.005
persistence_frames = 40
max_ids = 256
interval_ms = 20  # intervalle du timer

colors = cm.get_cmap('inferno', max_ids)


# 📥 Chargement spectro & partiels

spec = np.load('data/tracking/oboe/snail_spectro_oboe.npy').T
with open('data/tracking/oboe/snail_partials_oboe.pkl', 'rb') as f:
    data = pickle.load(f)

# spec = np.load('data/tracking/voice/snail_spectro_voice.npy').T
# with open('data/tracking/voice/snail_partials_voice.pkl', 'rb') as f:
#     data = pickle.load(f)

frame_keys = sorted(data.keys())
history = []

# 🕒 Prépa axes spectro
num_freq_bins, num_frames = spec.shape
time = np.arange(num_frames) * hopsize
frequencies = np.linspace(0, fs/2, num_freq_bins)
T, F = np.meshgrid(time, frequencies)

# 🎨 Figure
fig, (ax_spec, ax_traj) = plt.subplots(1, 2, figsize=(14,6))
plt.subplots_adjust(bottom=0.30)

# ─── Spectrogramme ───────────────────────────────────────────────────────────
spec_plot = ax_spec.pcolormesh(T, F, spec, shading='auto', cmap='magma')
ax_spec.set_title('Spectrogramme')
ax_spec.set_xlabel('Temps (s)')
ax_spec.set_ylabel('Fréquence (Hz)')
ax_spec.set_ylim([0, 5000])
bar_line = ax_spec.axvline(0, color='cyan', linewidth=2)

# état global
current_frame = 0
anim_running = False

# ─── Tracé d’une trame ────────────────────────────────────────────────────────
def draw_frame(frame_index):
    global history
    frame = frame_keys[frame_index]
    partials = data[frame]
    
    omega      = partials[:,0]
    log_amps   = partials[:,1]
    part_ids   = partials[:,4].astype(int)
    freqs_hz   = omega * fs/(2*np.pi)
    
    history.append((freqs_hz, log_amps, part_ids))
    if len(history)>persistence_frames:
        history.pop(0)
    
    # partiels à droite
    ax_traj.clear()
    for i,(f,a,ids) in enumerate(history):
        alpha = 1-(i/persistence_frames)
        for j in range(len(f)):
            c = colors(ids[j]%max_ids)
            ax_traj.scatter(f[j], a[j], c=[c], marker='.', alpha=alpha)
            #ax_traj.plot(f[j], a[j])
    ax_traj.set_xlabel('Fréquences (Hz)')
    ax_traj.set_ylabel('Amplitude (dB)')
    ax_traj.set_title(f'Trame {frame}')
    ax_traj.set_xlim([0,5000])
    ax_traj.set_ylim([-8,0])
    ax_traj.grid(True)
    
    # barre sur spectro
    bar_line.set_xdata([time[frame_index]])
    
    # synchro slider
    frame_slider.eventson = False
    frame_slider.set_val(frame_index)
    frame_slider.eventson = True
    
    fig.canvas.draw_idle()

# ─── Timer pour l’animation ────────────────────────────────────────────────────
timer = fig.canvas.new_timer(interval=interval_ms)
def on_timer():
    global current_frame
    draw_frame(current_frame)
    current_frame += 1
    if current_frame>=num_frames:
        timer.stop()
    # pas de boucle infinie ici, si tu veux, remets current_frame=0

timer.add_callback(on_timer)

# ─── Slider ───────────────────────────────────────────────────────────────────
slider_ax = plt.axes([0.2, 0.15, 0.6, 0.03])
frame_slider = Slider(slider_ax, 'Trame', 0, num_frames-1,
                      valinit=0, valstep=1)
def slider_update(val):
    global anim_running, current_frame, history
    timer.stop()
    anim_running = False
    play_button.label.set_text('Play')
    current_frame = int(val)
    history = []  # si tu veux repartir proprement de la traînée
    draw_frame(current_frame)
frame_slider.on_changed(slider_update)

# ─── Bouton Play/Pause ───────────────────────────────────────────────────────
button_ax = plt.axes([0.82, 0.05, 0.1, 0.05])
play_button = Button(button_ax, 'Play', hovercolor='0.975')
def toggle(event):
    global anim_running
    if anim_running:
        timer.stop()
        play_button.label.set_text('Play')
        anim_running = False
    else:
        timer.start()
        play_button.label.set_text('Pause')
        anim_running = True

play_button.on_clicked(toggle)

# ─── Démarrage initial ───────────────────────────────────────────────────────
draw_frame(0)
plt.show()


