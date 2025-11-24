import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation

# --- 1. Data ---
num = 120
x = np.linspace(-1.4, 1.4, num)
y = np.linspace(-1.4, 1.4, num)
X, Y = np.meshgrid(x, y)
Exp_term = np.exp(X + Y)
Linear_term = -Y - 1

C_values = [-0.5, 0.2, 1.0, 2.0, 3.5]

# --- 2. Aesthetics ---
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection="3d")
bg_color = "#2b2b2b"
ax.set_facecolor(bg_color)
fig.patch.set_facecolor(bg_color)
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.tick_params(axis="x", colors="white")
ax.tick_params(axis="y", colors="white")
ax.tick_params(axis="z", colors="white")
ax.set_xlabel("x", color="white")
ax.set_ylabel("y", color="white")
ax.set_zlabel("u", color="white")

colormap = cm.get_cmap("plasma", len(C_values))

# --- 3. Plot Leaves ---
for i, C in enumerate(C_values):
    U_leaf = C * Exp_term + Linear_term
    color = colormap(i)
    ax.plot_surface(
        X,
        Y,
        U_leaf,
        color=color,
        alpha=0.5,
        edgecolor="none",
        shade=True,
        rstride=2,
        cstride=2,
    )

ax.set_title(
    r"Foliação: $u = C e^{x+y} - y - 1$",
    color="white",
    fontsize=16,
)
ax.set_zlim(-4, 10)

# --- 4. Animation ---
FPS = 30
frames = 360


def update(frame):
    angle = (frame / frames) * 360
    ax.view_init(elev=15, azim=-45 + angle)
    return (fig,)


ani = FuncAnimation(fig, update, frames=frames, interval=1000 / FPS, blit=False)
print("Saving foliation video...")
output_file = "sys2_foliation.mp4"
try:
    print(f"Saving video to '{output_file}'... please wait.")
    ani.save(output_file, writer="ffmpeg", fps=FPS, bitrate=6000)
    print("Done! Video saved successfully.")
except FileNotFoundError:
    print("\n--- FFmpeg Not Found ---")
    print("Could not save video file. Showing live animation instead.")
    plt.show()
