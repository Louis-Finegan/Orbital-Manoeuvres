import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Read from the saved csv (run ThreeBodyTLI beforehand)
df = pd.read_csv("../data/lunar_transfer.csv")
print(f"Length of file: {len(df)}")

# Moon and Spacecraft positions relative to the Earth
x_moon_ec = df["moon_x"] - df["earth_x"]
y_moon_ec = df["moon_y"] - df["earth_y"]
x_sat_ec  = df["sat_x"]  - df["earth_x"]
y_sat_ec  = df["sat_y"]  - df["earth_y"]

# ------------------------------------------------------------
# Plot Trajectory
# ------------------------------------------------------------
plt.figure(figsize=(10, 10))
plt.plot(x_sat_ec, y_sat_ec, linewidth=1, label="Satellite trajectory", color='blue', alpha=0.7)
plt.plot(x_moon_ec, y_moon_ec, linewidth=1, label="Moon trajectory", color='gray', alpha=0.7)
plt.scatter(0, 0, color='red', s=150, label="Earth", zorder=5)
plt.axis("equal")
plt.legend()
plt.title("Trans-Lunar-Injection")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.grid(True, alpha=0.3)
plt.savefig("../figures/trajectory.png", dpi=150)
plt.close()
print("Saved static plot: trajectory.png")

# ------------------------------------------------------------
# Animated Trajectory
# ------------------------------------------------------------
trail_length = 200
step = 1000 # I don't want to kill my computer

fig, ax = plt.subplots(figsize=(10, 10))

sat_trail, = ax.plot([], [], lw=1, alpha=0.7, label="Satellite trail", color='green')
moon_trail, = ax.plot([], [], lw=1, alpha=0.7, label="Moon trail", color='gray')
sat_marker, = ax.plot([], [], 'o', markersize=8, color='green', label="Satellite")
moon_marker, = ax.plot([], [], 'o', markersize=10, color='gray', label="Moon")
ax.scatter(0, 0, color='blue', s=200, label="Earth", zorder=5)

max_dist_ec = max(
    np.max(np.abs(x_sat_ec)), np.max(np.abs(y_sat_ec)),
    np.max(np.abs(x_moon_ec)), np.max(np.abs(y_moon_ec))
) * 1.1
ax.set_xlim(-max_dist_ec, max_dist_ec)
ax.set_ylim(-max_dist_ec, max_dist_ec)
ax.set_aspect('equal')
ax.legend(loc="upper right")
ax.set_title("Trans-Lunar-Injection")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.grid(True, alpha=0.3)

frames = range(0, len(df), step)

def update_ec(frame_idx):
    i = frames[frame_idx]
    sat_trail.set_data(x_sat_ec[:i], y_sat_ec[:i])
    moon_trail.set_data(x_moon_ec[:i], y_moon_ec[:i])
    sat_marker.set_data([x_sat_ec[i]], [y_sat_ec[i]])
    moon_marker.set_data([x_moon_ec[i]], [y_moon_ec[i]])
    if hasattr(ax, 'time_text'):
        ax.time_text.remove()
    time_str = f"Time: {df['time_s'].iloc[i]/86400:.2f} days"
    ax.time_text = ax.text(0.02, 0.98, time_str, transform=ax.transAxes,
                           fontsize=10, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    return sat_trail, moon_trail, sat_marker, moon_marker, ax.time_text

ani_ec = FuncAnimation(fig, update_ec, frames=range(len(frames)),
                       interval=20, blit=False, repeat=True)
ani_ec.save("../figures/lunar_transfer_animation.mp4", writer="ffmpeg", fps=30)
print("Saved Earth-centered animation: lunar_transfer_animation.mp4")
plt.close()