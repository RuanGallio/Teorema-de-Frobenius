import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation

# --- 1. Data Generation ---
num_points = 120
x = np.linspace(-2.5, 2.5, num_points)
y = np.linspace(-2.5, 2.5, num_points)
X, Y = np.meshgrid(x, y)
Base_U = 0.5 * (X**2 + Y**2)
C_values = [-2, -1, 0, 1, 2, 3]

# --- 2. Visualization Setup ---
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection="3d")

# Aesthetic
bg_color = "#2b2b2b"
text_color = "white"
fig.patch.set_facecolor(bg_color)
ax.set_facecolor(bg_color)
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
    axis.set_tick_params(colors=text_color)
    axis.label.set_color(text_color)
ax.set_xlabel("x", fontsize=12)
ax.set_ylabel("y", fontsize=12)
ax.set_zlabel("u", fontsize=12)

colormap = cm.get_cmap("plasma", len(C_values))

print("Plotting static surfaces for animation...")
# Plot the static leaves once
for i, C in enumerate(C_values):
    U_leaf = Base_U + C
    color = colormap(i)
    ax.plot_surface(
        X,
        Y,
        U_leaf,
        color=color,
        alpha=0.4,
        edgecolor="none",
        rstride=2,
        cstride=2,
        shade=True,
        antialiased=True,
    )

ax.set_title(
    r"Foliation by Solution Leaves - Rotating View",
    color=text_color,
    fontsize=15,
    pad=20,
)
ax.set_zlim(-3, 6)
ax.set_xlim(-2.5, 2.5)
ax.set_ylim(-2.5, 2.5)

# --- 3. Animation Configuration ---
FPS = 30
DURATION_SEC = 12
TOTAL_FRAMES = FPS * DURATION_SEC
START_ELEV = 25
START_AZIM = -45


def update_view(frame):
    rotation_angle = (frame / TOTAL_FRAMES) * 360
    ax.view_init(elev=START_ELEV, azim=START_AZIM + rotation_angle)
    return (fig,)


print(f"Starting animation rendering ({TOTAL_FRAMES} frames).")
ani = FuncAnimation(
    fig, update_view, frames=TOTAL_FRAMES, interval=1000 / FPS, blit=False
)

# --- 4. Saving the Video ---
output_file = "foliation_rotating_dark.mp4"

try:
    print(f"Saving video to '{output_file}'... please wait.")
    ani.save(output_file, writer="ffmpeg", fps=FPS, bitrate=6000)
    print("Done! Video saved successfully.")
except FileNotFoundError:
    print("\n--- FFmpeg Not Found ---")
    print("Could not save video file. Showing live animation instead.")
    plt.show()
