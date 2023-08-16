"""
Instructions:
1. Run Weekday4_preprocessing.py first
2. Run this script
"""

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
    multipart_filenames
)

time0 = time()
# Rotation matrix from geographic frame to cage frame
R_G2C = EulerRotation().rz(-23, deg=True)
# Rotation matrix from sensor frame to cage frame
R_S2C = EulerRotation().rx(-180, deg=True)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
local_emf = np.dot(R_G2C, local_emf())

# Data processing ============================================================

# Filenames in chronological order
filename_base = "./data/Weekday4_2023-07-13_08.29.31_SPARSE25.dat"

filenames = multipart_filenames(filename_base, 27, breakpoint=9+4)

# For data processing, only use 24 hours of data, so shrink file selection
# In this case, 09:00 on 2023-07-13 to 09:00 on 2023-07-14
filenames = filenames[1:-2]


# Pre-allocate storage lists
mean_data = [np.array([0.]*3)]*len(filenames)
mean_data_n = [np.array([0.]*3)]*len(filenames)
sd_data_n = [np.array([0.]*3)]*len(filenames)


print("Starting data processing...")
for i, filename in enumerate(filenames):

    print("Starting segment {} of {} (elapsed time: {} s)"
          .format(i+1, len(filenames), round(time()-time0, 3)))
    # Readout files
    data = read_data(filename, header=True)

    # Apply Savitzky-Golay smoothing
    data = data_savgol_filter(data, 51, 4)

    # Rotate data from the sensor frame S into the cage frame C
    data = data_rotate(data, R_S2C)

    # Un-normalized data
    mean_data[i], _ = data_normal_distribution(data)

    # Normalize the data (in C frame) with respect to the EMF (in C frame):
    data = data_add_vector(data, -local_emf)

    # Distribution properties of normalized disturbance datapoints in C frame
    mean_data_n[i], sd_data_n[i] = data_normal_distribution(data)


# Easy Copy-Paste readout for Excel file
easy_copy_paste_readout = True

if easy_copy_paste_readout:
    for i in range(len(filenames)):
        print("{} {} {} {} {} {} {}".format(
            round(mean_data_n[i][0], 3),
            round(mean_data_n[i][1], 3),
            round(mean_data_n[i][2], 3),
            round(sd_data_n[i][0], 3),
            round(sd_data_n[i][1], 3),
            round(sd_data_n[i][2], 3),
            round(np.linalg.norm(mean_data_n[i]), 3)
        )
    )


# Vector plot
print("Starting vector plot")

plot_title = """
    Magnetic disturbance vector at the test site during a regular day
    Measured on 13-07-2023, part of a 24+ hour measurement"""

vp = VectorPlot()

# Vector scaling factor: 1 uT corresponds to 2 cm on the 3D plot
vscf = 0.02
vp.autoplot(tripod=True, coils=True, walls=True,
                    table=True)
for i in range(len(mean_data_n)):
    vp.plot_vector(mean_data_n[i], color="black",
                   scaling=vscf, linewidth=3, alr=0.2)
vp.plot_vector(local_emf, linewidth=3, scaling=vscf, alr=0.2, color="orange")

# Mean data vector
mdv = np.mean(mean_data, axis=0)
vp.plot_vector(mdv, scaling=vscf, linewidth=3, color="purple")
# Transparent line from mdv to emf vector
vp.plot_vector((local_emf-mdv)*vscf, origin=mdv*vscf, alpha=0.5, linewidth=1,
               scaling=1, alr=0.001, color="#555")

# Mean normalized data vector (disturbance vector)
mnv = np.mean(mean_data_n, axis=0)
# Transparent line from mdv to disturbance vector
vp.plot_vector((mnv-mdv)*vscf, origin=mdv*vscf, alpha=0.5, linewidth=1,
               scaling=1, alr=0.001, color="#555")

vp.ax.text2D(0.15, 0.9, plot_title,
             fontsize=14,
             multialignment="center",
             transform=vp.ax.transAxes)

vp.show()
vp.fig.savefig("./figures/Weekday4_vectorplot.png", dpi=150)

print("Elapsed time:", time()-time0)
