# Imports ====================================================================
import numpy as np
from time import time

from tools.local_emf import local_emf
from tools.plotting import VectorPlot
from tools.GaussmeterLib import (
    EulerRotation,
    read_data,
    data_normal_distribution,
    data_savgol_filter,
    data_rotate,
    data_add_vector,
)

time0 = time()
# Rotation matrix from geographic frame to cage frame
R_G2C = EulerRotation().rz(-23, deg=True)
# Rotation matrix from sensor frame to cage frame
R_S2C = EulerRotation().rx(-180, deg=True)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
local_emf = np.dot(R_G2C, local_emf())

# Data processing ============================================================

# Filenames
filenames = (
    "./data/RackElimination_2023-07-04_12.51.32.dat",
    "./data/RackEliminationControl_2023-07-03_13.00.00.dat"
)

# Pre-allocate storage lists
mean_data_n = [0.]*len(filenames)
sd_data_n = [0.]*len(filenames)

for i, filename in enumerate(filenames):
    # Readout files
    data = read_data(filename, header=True)

    # Apply Savitzky-Golay smoothing
    data = data_savgol_filter(data, 51, 4)

    # Rotate data from the sensor frame S into the cage frame C
    data = data_rotate(data, R_S2C)

    # Normalize the data (in C frame) with respect to the EMF (in C frame):
    data = data_add_vector(data, -local_emf)

    # Distribution properties of normalized disturbance datapoints in C frame
    mean_data_n[i], sd_data_n[i] = data_normal_distribution(data)

# Easy Copy-Paste readout
easy_copy_paste_readout = True
if easy_copy_paste_readout:
    for j in (0, 1, 2):
        for i in range(len(filenames)):
            print(round(mean_data_n[i][j], 3))
        print(" ")

    for j in (0, 1, 2):
        for i in range(len(filenames)):
            print(round(sd_data_n[i][j], 3))
        print(" ")

    for i in range(len(filenames)):
        print(round(np.linalg.norm(mean_data_n[i]), 3))
    print(" ")

# Vector plot
plot_title = """
    Average disturbance vectors with the storage rack present and removed
    Measured on 03-07-2023 and 04-07-2023
    """
vectorplot = VectorPlot()
vectorplot.autoplot(tripod=True, coils=True, walls=True,
                    table=True)
vectorplot.plot_vector(local_emf, linewidth=3, scaling=0.02, alr=0.2, color="orange")
for i in range(len(mean_data_n)):
    vectorplot.plot_vector(mean_data_n[i], color="black",
                           linewidth=3, alr=0.2)

vectorplot.ax.text2D(0.05, 0.9, plot_title,
                     fontsize=14,
                     multialignment="center",
                     transform=vectorplot.ax.transAxes)
vectorplot.show()
vectorplot.fig.savefig("./figures/RackElimination.png", dpi=150)

print("Elapsed time:", time()-time0)
