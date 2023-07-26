# Imports ====================================================================
import numpy as np
from tools.local_emf import local_emf
from tools.plotting import SpectralPlotPyQt
from tools.GaussmeterLib import (
    EulerRotation,
    data_rfft,
    read_data,
    data_normal_distribution,
    data_savgol_filter,
    data_rotate,
    data_add_vector,
)

# Rotation matrix from geographic frame to cage frame
R_G2C = EulerRotation().rz(-23, deg=True)
# Rotation matrix from sensor frame to cage frame
R_S2C = EulerRotation().rx(-180, deg=True)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
local_emf = np.dot(R_G2C, local_emf())

# Data processing ============================================================
filenameON = "./data/2500HzTestLightsOn_2023-07-05_11.04.28.dat"
filenameOFF = "./data/2500HzTestLightsOff_2023-07-05_11.13.36.dat"

# Import data
dataON = read_data(filenameON, header=True)
dataOFF = read_data(filenameOFF, header=True)

# Apply Savitzky-Golay smoothing
dataON = data_savgol_filter(dataON, 51, 4)
dataOFF = data_savgol_filter(dataOFF, 51, 4)

# Rotate data from the sensor frame S into the cage frame C
dataON = data_rotate(dataON, R_S2C)
dataOFF = data_rotate(dataOFF, R_S2C)

# Normalize the data (in C frame) with respect to the EMF (in C frame):
dataON = data_add_vector(dataON, -local_emf)
dataOFF = data_add_vector(dataOFF, -local_emf)

# Value readout
mean_data_ON, sd_data_ON = data_normal_distribution(dataON)
mean_data_OFF, sd_data_OFF = data_normal_distribution(dataOFF)
easy_copy_paste_readout = True
if easy_copy_paste_readout:
    for j in (0, 1, 2):
        for mean in (mean_data_ON, mean_data_OFF):
            print(round(mean[j], 3))
        print(" ")
    for j in (0, 1, 2):
        for sd in (sd_data_ON, sd_data_OFF):
            print(round(sd[j], 3))
        print(" ")
    for mean in (mean_data_ON, mean_data_OFF):
        print(round(np.linalg.norm(mean), 3))
    print(" ")

# Obtain one-sided FFT data
[f, x_t, y_t, z_t] = data_rfft(dataON)
data_f_ON = data_rfft(dataON)

data_f_OFF = data_rfft(dataOFF)

spectralplotON = SpectralPlotPyQt()
# spectralplotON.set_window_size(1920, 1080)
spectralplotON.set_window_size(1200, 800)
spectralplotON.set_pen_alpha(1)
spectralplotON.set_plot_title("Magnetic frequency spectra during day with ceiling lights ON")
spectralplotON.set_xlabel("Frequency [Hz]")
spectralplotON.set_ylabel("Field magnitude [\u03BCT / \u221AHz]")
spectralplotON.save_plot("./figures/Lights On - Spectral Plot.png")
spectralplotON.spectralplot_pyqtgraph(data_f_ON)

spectralplotOFF = SpectralPlotPyQt()
# spectralplotOFF.set_window_size(1920, 1080)
spectralplotOFF.set_window_size(1200, 800)
spectralplotOFF.set_pen_alpha(1)
spectralplotOFF.set_plot_title("Magnetic frequency spectra during day with ceiling lights OFF")
spectralplotOFF.set_xlabel("Frequency [Hz]")
spectralplotOFF.set_ylabel("Field magnitude [\u03BCT / \u221AHz]")
spectralplotOFF.save_plot("./figures/Lights Off - Spectral Plot.png")
spectralplotOFF.spectralplot_pyqtgraph(data_f_OFF)
