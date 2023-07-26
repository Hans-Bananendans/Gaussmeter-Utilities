# Imports ====================================================================
import numpy as np

from tools.local_emf import local_emf
from tools.plotting import TimeplotPyQt, SpectralplotPyQt
from tools.GaussmeterAnalysis import (
    EulerRotation,
    data_rfft,
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

# Obtain one-sided FFT data
# [[x_t1, x_f1], [y_t1, y_f1], [z_t1, z_f1]] = data_rfft(data)
[f, x_t, y_t, z_t] = data_rfft(dataON)
data_f_ON = data_rfft(dataON)

data_f_OFF = data_rfft(dataOFF)


spectralplotON = SpectralplotPyQt()
# spectralplotON.set_window_size(1920, 1080)
spectralplotON.set_window_size(1200, 800)
spectralplotON.set_pen_alpha(1)
spectralplotON.set_plot_title("Magnetic frequency spectra during day with ceiling lights ON")
spectralplotON.set_xlabel("Frequency [Hz]")
spectralplotON.set_ylabel("Field magnitude [\u03BCT / \u221AHz]")
spectralplotON.save_plot("./figures/Lights On - Spectral Plot.png")
spectralplotON.spectralplot_pyqtgraph(data_f_ON)

spectralplotOFF = SpectralplotPyQt()
# spectralplotOFF.set_window_size(1920, 1080)
spectralplotOFF.set_window_size(1200, 800)
spectralplotOFF.set_pen_alpha(1)
spectralplotOFF.set_plot_title("Magnetic frequency spectra during day with ceiling lights OFF")
spectralplotOFF.set_xlabel("Frequency [Hz]")
spectralplotOFF.set_ylabel("Field magnitude [\u03BCT / \u221AHz]")
spectralplotOFF.save_plot("./figures/Lights Off - Spectral Plot.png")
spectralplotOFF.spectralplot_pyqtgraph(data_f_OFF)

# # Plot Spectrogram LightOn XYZ
# fig1 = plt.figure(figsize=(8, 8))
# gs = GridSpec(1,1)
# ax = fig1.add_subplot(gs[0,0])
# # ax.set_ylim(ylim[0], ylim[1])
#
#
# filter_dc: int = 0
#
# ax.plot(f[filter_dc:], np.abs(x_t[filter_dc:]), color='r', label="X", alpha=0.5)
# ax.plot(f[filter_dc:], np.abs(y_t[filter_dc:]), color='g', label="Y", alpha=0.4)
# ax.plot(f[filter_dc:], np.abs(z_t[filter_dc:]), color='b', label="Z", alpha=0.3)
#
# ax.set_title("Spectrogram of magnetic measurements")
# ax.set_xlabel("f [Hz]")
# ax.set_ylabel("Linear magnitude")
# ax.legend()
# ax.grid(True)
# # ax.set_xscale('log')
# ax.set_yscale('log')
# plt.show()


# # Set up plot ================================================================
# plot_title = """
# <h3><b> Magnetic field components of a weekday test performed from {} to {}.
# </b> <br> Data fetched from file:  <samp> {} </samp> </h3>
# """.format(
#     datetime.utcfromtimestamp(data[0][0]).strftime("%d/%m/%y %H:%M:%S"),
#     datetime.utcfromtimestamp(data[0][-1]).strftime("%d/%m/%y %H:%M:%S"),
#     filename
# )
# side_label = "<h3>Magnetic field component B (uT)<h3>"
#
# plot_filename = "WeekdayTest1_timeplot.png"
#
# plot_object = TimeplotPyQt()
# plot_object.set_plot_title(plot_title)
# plot_object.set_side_label(side_label)
# plot_object.set_window_size(1920, 1080)
#
# plot_object.save_plot(plot_filename)
#
# # Generate plot
# plot_object.timeplot_pyqtgraph(data)
