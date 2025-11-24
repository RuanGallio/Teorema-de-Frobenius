import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Import animation tools
from matplotlib.animation import FuncAnimation

# --- 1. Data Generation (Same as before) ---
# Surface grid
x = np.linspace(-2, 2, 40)
y = np.linspace(-2, 2, 40)
X_mesh, Y_mesh = np.meshgrid(x, y)
U_mesh = 0.5 * (X_mesh**2 + Y_mesh**2)

# Vector grid lying on the surface
x_vec = np.linspace(-1.5, 1.5, 8)
y_vec = np.linspace(-1.5, 1.5, 8)
X_v, Y_v = np.meshgrid(x_vec, y_vec)
U_v = 0.5 * (X_v**2 + Y_v**2)

# Vector field components
# V1 = (1, 0, x)
V1_x, V1_y, V1_u = np.ones_like(X_v), np.zeros_like(X_v), X_v
# V2 = (0, 1, y)
V2_x, V2_y, V2_u = np.zeros_like(Y_v), np.ones_like(Y_v), Y_v

# --- 2. Visualization Setup ---
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

# Plot static elements
surf = ax.plot_surface(
    X_mesh,
    Y_mesh,
    U_mesh,
    cmap="viridis",
    alpha=0.7,
    edgecolor="none",
    rstride=1,
    cstride=1,
)
q1 = ax.quiver(
    X_v,
    Y_v,
    U_v,
    V1_x,
    V1_y,
    V1_u,
    color="red",
    length=0.4,
    normalize=True,
    label=r"$V_1 = \partial_x + x\partial_u$",
)
q2 = ax.quiver(
    X_v,
    Y_v,
    U_v,
    V2_x,
    V2_y,
    V2_u,
    color="blue",
    length=0.4,
    normalize=True,
    label=r"$V_2 = \partial_y + y\partial_u$",
)

# Aesthetics
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u")
ax.set_title(r"Integral Manifold and Tangent Vectors")
ax.set_zlim(0, 4)
ax.legend()

# --- 3. Animation Configuration ---
FPS = 30
DURATION_SEC = 10
TOTAL_FRAMES = FPS * DURATION_SEC
START_ELEV = 30
START_AZIM = -60


def update_view(frame):
    rotation_angle = (frame / TOTAL_FRAMES) * 360
    ax.view_init(elev=START_ELEV, azim=START_AZIM + rotation_angle)
    return (fig,)


print(f"Starting animation rendering ({TOTAL_FRAMES} frames)...")
ani = FuncAnimation(
    fig, update_view, frames=TOTAL_FRAMES, interval=1000 / FPS, blit=False
)

# --- 4. Saving the Video ---
output_file = "manifold_solution_rotating.mp4"

try:
    print(f"Saving video to '{output_file}'... please wait.")
    ani.save(output_file, writer="ffmpeg", fps=FPS, bitrate=6000)
    print("Done! Video saved successfully.")
except FileNotFoundError:
    print("\n--- FFmpeg Not Found ---")
    print("Could not save video file. Showing live animation instead.")
    print("To save videos, please install FFmpeg and ensure it's in your system PATH.")
    plt.show()
