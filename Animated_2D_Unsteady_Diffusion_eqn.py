import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parameters
N = 51
L = 1.0
h = L / (N - 1)
dt = 0.0001
alpha = dt / h**2
epsilon = 1e-5
max_iter = 100000
save_interval = 100

# Initialize temperature arrays
T = np.zeros((N, N))
T[N-1, :] = 1.0
T_new = np.zeros((N, N))
T_new[N-1, :] = 1.0

# List to store frames
frames = []

# Simulation loop
iterations = 0
while iterations < max_iter:
    for i in range(1, N-1):
        for j in range(1, N-1):
            T_new[i,j] = T[i,j] + alpha * (T[i+1,j] + T[i-1,j] + T[i,j+1] + T[i,j-1] - 4*T[i,j])
    
    num_error = np.sum(np.abs(T_new[1:N-1, 1:N-1] - T[1:N-1, 1:N-1]))
    T = T_new.copy()
    iterations += 1
    
    if iterations % save_interval == 0:
        frames.append(T.copy())
    
    if num_error < epsilon:
        break

# Create the animation
x_dom = np.linspace(0, L, N)
y_dom = np.linspace(0, L, N)
X, Y = np.meshgrid(x_dom, y_dom)

fig, ax = plt.subplots()
pcm = ax.imshow(frames[0], extent=[0, L, 0, L], origin='lower', vmin=0, vmax=1)
fig.colorbar(pcm, ax=ax, label='Temperature')
ax.set_xlabel('X')
ax.set_ylabel('Y')

def update(frame):
    pcm.set_data(frames[frame])
    ax.set_title(f'Iteration {frame * save_interval}')
    return pcm,

ani = animation.FuncAnimation(fig, update, frames=len(frames), interval=50)

# Save the animation
# To save as MP4, ensure ffmpeg is installed
#ani.save('heat_equation.mp4', writer='ffmpeg')

# Alternatively, save as GIF
#ani.save('heat_equation.gif', writer='pillow')

plt.show()