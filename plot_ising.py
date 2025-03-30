import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd

### Plot Energy & Magnetization vs Beta ###
# Load the data
data = pd.read_csv("ising_results.csv")

# Create figure for graphs
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# Energy vs Beta
axes[0].plot(data["Beta"], data["Average Energy"], marker='o', linestyle='-')
axes[0].set_xlabel("Beta (1/kT)")
axes[0].set_ylabel("Average Energy per Spin")
axes[0].set_title("Energy vs. Beta")
axes[0].grid()

# Magnetization vs Beta
axes[1].plot(data["Beta"], data["Average Magnetization"], marker='s', linestyle='-', color='r')
axes[1].set_xlabel("Beta (1/kT)")
axes[1].set_ylabel("Average Magnetization per Spin")
axes[1].set_title("Magnetization vs. Beta")
axes[1].grid()

plt.tight_layout()
plt.show()  # Show the graphs first


### Animate 2D Ising Model ###

# Read spin configurations from file
def read_spin_configurations(filename):
    configs = []
    with open(filename, "r") as file:
        config = []
        for line in file:
            if "END" in line:
                configs.append(np.array(config))  # Store current config
                config = []  # Reset
            else:
                config.append([int(x) for x in line.split()])
    return configs

# Load spin configurations
spin_configs = read_spin_configurations("ising_spins.txt")

# Create animation figure
fig, ax = plt.subplots()
im = ax.imshow(spin_configs[0], cmap="gray", animated=True)

def update(frame):
    im.set_array(spin_configs[frame])
    return [im]

ani = animation.FuncAnimation(fig, update, frames=len(spin_configs), interval=100)

plt.show()  # Show the animation
