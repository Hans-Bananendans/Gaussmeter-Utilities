# Imports ====================================================================
import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
from datetime import datetime

from local_emf import local_emf
from TimeplotPyQt import TimeplotPyQt
from GaussmeterAnalysis import (
    EulerRotation,
    read_data,
    data_normal_distribution,
    data_savgol_filter,
    data_rotate,
    data_add_vector,
    time_plot,
)

# Rotation matrix from geographic frame to cage frame
R_G2C = EulerRotation().rz(-23, deg=True)
# Rotation matrix from sensor frame to cage frame
R_S2C = EulerRotation().rx(-180, deg=True)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
local_emf = np.dot(R_G2C, local_emf())

# Data processing ============================================================
filename = "WeekdayTest_2023-07-03_10.41.03_SPARSE50.dat"

# Import data
data = read_data(filename, header=True)

# Apply Savitzky-Golay smoothing
data = data_savgol_filter(data, 51, 4)

# Rotate data from the sensor frame S into the cage frame C
data = data_rotate(data, R_S2C)

# Distribution properties of measurements in C frame
mean_data, sd_data = data_normal_distribution(data)
print("C-frame data: mean, sd:", mean_data.round(3), sd_data.round(3))

# Normalize the data (in C frame) with respect to the EMF (in C frame):
data = data_add_vector(data, -local_emf)

# Distribution properties of normalized disturbance datapoints in C frame
mean_data_n, sd_data_n = data_normal_distribution(data)
print("Normalized data: mean, sd:", mean_data_n.round(3), sd_data_n.round(3))
print("Absolute strength of disturbance:",
      round(np.linalg.norm(mean_data_n), 3),
      "uT (",
      round(100*(np.linalg.norm(mean_data_n)/np.linalg.norm(local_emf)), 2),
      "% of EMF)")
