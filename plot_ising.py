import matplotlib.pyplot as plt
import pandas as pd

# Load the data
data = pd.read_csv("ising_results.csv")

# Plot Energy vs Beta
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.plot(data["Beta"], data["Average Energy"], marker='o', linestyle='-')
plt.xlabel("Beta (1/kT)")
plt.ylabel("Average Energy per Spin")
plt.title("Energy vs. Beta")
plt.grid()

# Plot Magnetization vs Beta
plt.subplot(1, 2, 2)
plt.plot(data["Beta"], data["Average Magnetization"], marker='s', linestyle='-', color='r')
plt.xlabel("Beta (1/kT)")
plt.ylabel("Average Magnetization per Spin")
plt.title("Magnetization vs. Beta")
plt.grid()

plt.tight_layout()
plt.show()
