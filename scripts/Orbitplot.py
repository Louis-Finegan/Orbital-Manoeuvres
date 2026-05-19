import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load data
df = pd.read_csv("../data/orbit.csv")

# Spacecraft position relative to the Earth
x_sc_ec = df["sat_x_m"] - df["earth_x_m"]
y_sc_ec = df["sat_y_m"] - df["earth_y_m"]

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
plt.figure(figsize=(8, 8))
plt.plot(x_sc_ec, y_sc_ec, linewidth=1, label="Satellite trajectory")
plt.scatter(0, 0, color='blue', s=100, label="Earth", zorder=5)
plt.axis("equal")
plt.legend()
plt.title("Satellite Motion relative to Earth")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.grid(True, alpha=0.3)
plt.savefig("../figures/trajectory_orbit.png", dpi=150)
plt.close()
print("Saved static plot: trajectory_orbit.png")

# ------------------------------------------------------------
# Animation
# ------------------------------------------------------------
trail_length = 150 
step = 10

# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 8))
trail_line, = ax.plot([], [], lw=1, alpha=0.7, label="Satellite trail")
sat_marker, = ax.plot([], [], 'o', markersize=6, color='green', label="Satellite")

# Earth marker
ax.scatter(0, 0, color='blue', s=150, label="Earth", zorder=5)
max_range = max(np.max(np.abs(x_sc_ec)), np.max(np.abs(y_sc_ec))) * 1.1
ax.set_xlim(-max_range, max_range)
ax.set_ylim(-max_range, max_range)
ax.set_aspect('equal')
ax.legend(loc="upper right")
ax.set_title("Satellite Motion relative to Earth")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.grid(True, alpha=0.3)

def update(frame_idx):
    i = frames[frame_idx] 
    trail_line.set_data(x_sc_ec[:i], y_sc_ec[:i])
    sat_marker.set_data([x_sc_ec[i]], [y_sc_ec[i]])

    return trail_line, sat_marker


frames = range(0, len(df), step) # saves time

ani = FuncAnimation(fig, update, frames=range(len(frames)), interval=20, blit=False, repeat=True)
ani.save("../figures/satellite_orbit.mp4", writer="ffmpeg", fps=30)
print("Saved animation: satellite_orbit.mp4")
plt.close()