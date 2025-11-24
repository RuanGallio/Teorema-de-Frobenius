import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# --- 1. Data Generation ---
x = np.linspace(-1.3, 1.3, 50)
y = np.linspace(-1.3, 1.3, 50)
X_mesh, Y_mesh = np.meshgrid(x, y)
U_mesh = np.exp(X_mesh + Y_mesh) - Y_mesh - 1

x_vec = np.linspace(-1.0, 1.0, 8)
y_vec = np.linspace(-1.0, 1.0, 8)
X_v, Y_v = np.meshgrid(x_vec, y_vec)
U_v = np.exp(X_v + Y_v) - Y_v - 1

V1_x, V1_y, V1_u = np.ones_like(X_v), np.zeros_like(X_v), 1 + Y_v + U_v
V2_x, V2_y, V2_u = np.zeros_like(X_v), np.ones_like(X_v), Y_v + U_v

# --- 2. Setup Plot ---
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

surf = ax.plot_surface(
    X_mesh,
    Y_mesh,
    U_mesh,
    cmap="viridis",
    alpha=0.8,
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
    length=0.25,
    normalize=True,
    label=r"$X$",
)
q2 = ax.quiver(
    X_v,
    Y_v,
    U_v,
    V2_x,
    V2_y,
    V2_u,
    color="blue",
    length=0.25,
    normalize=True,
    label=r"$Y$",
)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u")
ax.set_title(r"Solution $u = e^{x+y} - y - 1$ and Vector Fields")
ax.set_zlim(-2, 6)  # Taller Z limit
ax.legend()

# --- 3. Animation ---
FPS = 30
frames = 300


def update(frame):
    angle = (frame / frames) * 360
    ax.view_init(elev=25, azim=-55 + angle)
    return (fig,)


ani = FuncAnimation(fig, update, frames=frames, interval=1000 / FPS, blit=False)
print("Saving solution video...")
output_file = "sys2_solution.mp4"
try:
    print(f"Saving video to '{output_file}'... please wait.")
    ani.save(output_file, writer="ffmpeg", fps=FPS, bitrate=6000)
    print("Done! Video saved successfully.")
except FileNotFoundError:
    print("\n--- FFmpeg Not Found ---")
    print("Could not save video file. Showing live animation instead.")
    plt.show()
