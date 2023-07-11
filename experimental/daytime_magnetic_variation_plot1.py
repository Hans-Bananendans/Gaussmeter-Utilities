# Imports ====================================================================
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore
from datetime import datetime
from scipy.signal import savgol_filter
from copy import deepcopy


# Helper functions =========================================================
def Rx(theta, deg=False):
    if deg:
        theta *= np.pi / 180
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(theta), -np.sin(theta)],
                   [0, np.sin(theta), np.cos(theta)]])
    return Rx


def Ry(theta, deg=False):
    if deg:
        theta *= np.pi / 180
    Ry = np.array([[np.cos(theta), 0, np.sin(theta)],
                   [0, 1, 0],
                   [-np.sin(theta), 0, np.cos(theta)]])
    return Ry


def Rz(theta, deg=False):
    if deg:
        theta *= np.pi / 180
    Rz = np.array([[np.cos(theta), -np.sin(theta), 0],
                   [np.sin(theta), np.cos(theta), 0],
                   [0, 0, 1]])
    return Rz


def read_data(filename, header=False):
    with open(filename, 'r') as f:
        lines = f.readlines()
        if lines[0][0:2] == "!H" or header is True:
            lines = lines[10:]

        t, x, y, z = [0] * len(lines), [0] * len(lines), [0] * len(lines), [0] * len(lines)

        for i in range(len(lines)):
            [tt, xx, yy, zz] = lines[i].split(" ")
            [t[i], x[i], y[i], z[i]] = \
                [np.double(tt), float(xx), float(yy), float(zz)]
        return [t, x, y, z]


def zipAB():
    pass


def zipBA(dataB):
    # Zips data from type B to type A
    # Type A: [
    #   [t_0 ... t_n],
    #   [Bx_0 ... Bx_n],
    #   [By_0 ... By_n],
    #   [Bz_0 ... Bz_n],
    # ]
    #
    # Type B: [
    #   [t_0, np.array(Bx_0, By_0, Bz_0)] ... [t_n, np.array(Bx_n, By_n, Bz_n)]
    # ]
    n = len(dataB)
    dataA = [[0.] * n, [0.] * n, [0.] * n, [0.] * n]
    for i in range(len(dataB)):
        dataA[0][i] = dataB[i][0]
        dataA[1][i] = dataB[i][1][0]
        dataA[2][i] = dataB[i][1][1]
        dataA[3][i] = dataB[i][1][2]
    return dataA


def vector_normal_distribution(vectors: list):
    mean_X = 0
    mean_Y = 0
    mean_Z = 0

    sd_X = 0
    sd_Y = 0
    sd_Z = 0

    for vector in vectors:
        mean_X += vector[0]
        mean_Y += vector[1]
        mean_Z += vector[2]

    mean_X = mean_X / len(vectors)
    mean_Y = mean_Y / len(vectors)
    mean_Z = mean_Z / len(vectors)

    for vector in vectors:
        sd_X += (vector[0] - mean_X) ** 2
        sd_Y += (vector[1] - mean_Y) ** 2
        sd_Z += (vector[2] - mean_Z) ** 2

    sd_X = (sd_X / len(vectors)) ** 0.5
    sd_Y = (sd_Y / len(vectors)) ** 0.5
    sd_Z = (sd_Z / len(vectors)) ** 0.5

    return [np.array([mean_X, mean_Y, mean_Z]), np.array([sd_X, sd_Y, sd_Z])]


def data_savgol_filter(data, windowlength, polyorder=4):
    # Smooths three data channels using a Savitzky-Golay filter
    # Data must be in type A
    x_filtered = savgol_filter(data[1], windowlength, polyorder)
    y_filtered = savgol_filter(data[2], windowlength, polyorder)
    z_filtered = savgol_filter(data[3], windowlength, polyorder)
    return [data[0], x_filtered, y_filtered, z_filtered]


# Set up Qt6 environment =====================================================
app = pg.mkQApp("Plot")
win = pg.GraphicsLayoutWidget(show=True,
                              title="Plotting items examples"
                              )
win.resize(1000, 600)
pg.setConfigOptions(antialias=True)  # Antialiasing for nicer plots


# EMF vector =================================================================
R_G2C = Rz(-23, deg=True)
R_S2C = Rx(-180, deg=True)

# EMF vector at AE building (WMM-2020 model)
B_EMF_WMM2020 = np.array([0.7216, 19.1629, -45.4592])
# EMF vector at AE building (IGRF2020 model)
B_EMF_IGRF2020 = np.array([0.7313, 19.1870, -45.4557])

# Take average of the two as baseline B_EMF:
B_EMF_G = 0.5 * (B_EMF_WMM2020 + B_EMF_IGRF2020)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
B_EMF = np.dot(R_G2C, B_EMF_G)


# Data processing ============================================================
filename = "WeekdayTest_2023-07-03_10.41.03_SPARSE50.dat"

data_r = read_data(filename, header=True)

# print("Dataset contains", len(data_r[0]), "data points.")

data_r = data_savgol_filter(data_r, 51, 4)

data = []
for i in range(len(data_r[0])):
    data.append([
        data_r[0][i],
        np.dot(R_S2C, np.array([data_r[1][i],
                                data_r[2][i],
                                data_r[3][i]]))
    ])

# Put all the vectors in a list
data_vectors = []
for d in data:
    data_vectors.append(d[1])

# mean_data, sd_data = vector_normal_distribution(data_vectors)
# print("C-frame data: mean, sd:", mean_data.round(3), sd_data.round(3))

# Normalize the data (in C frame) with respect to the   EMF (in C frame):
data_n = deepcopy(data)

for i in range(len(data_n)):
    data_n[i][1] = data_n[i][1] - B_EMF

data_n = zipBA(data_n)

# data_n_vectors = []
# for d in data_n:
#     data_n_vectors.append(d[1])
#
# mean_data_n, sd_data_n = vector_normal_distribution(data_n_vectors)


# Plotting ===================================================================

data_t = np.array(data_n[0])
data_t_date = []
for i in range(len(data_t)):
    data_t_date.append(data_t[i])
data_x = np.array(data_n[1])

xmin = data_n[0][0]
xmax = data_n[0][-1]
ymin = min(data_n[1])
ymax = max(data_n[1])
padding = 1.1

# Create plot object:
pen_red = pg.mkPen(color=(255, 0, 0))
pen_green = pg.mkPen(color=(0, 255, 0))
pen_blue = pg.mkPen(color=(0, 0, 255))

p_x = win.addPlot(title="X", x=data_t, y=data_x, pen=pen_red,
                  xmin=xmin, xmax=xmax,
                  axisItems={'bottom': pg.DateAxisItem()})

# p_x.setXRange(xmin*padding, xmax*padding)
# p_x.setYRange(-20, 20)

p_x.showGrid(x=True, y=True)
p_x.setTitle("Title", size="24pt")
# p_x.setLabel("left", "B (uT)", size="18pt")
# p_x.setLabel("bottom", "UNIX time (s)", size="18pt")
# p_x.addLegend()


def generate_time_label(value):
    # TODO: So far no luck
    utc = datetime.utcfromtimestamp(value)
    # return f"{utc.time()}\n{utc.date()}"
    return "{TEST}"


inf1 = pg.InfiniteLine(angle=90,
                       label="{value}",
                       # label="{generate_time_label}",
                       # label=generate_time_label,
                       # label="{{:%Y-%m-%d %H:%M:%S}.format(value)}",
                       # label=datetime.utcfromtimestamp(float("{value}")).time(),
                       # label="Value: {value:.0f}",
                       # label="{:s}".format("test"),  # works
                       pos=xmin,
                       movable=True,
                       bounds=[xmin, xmax],
                       labelOpts={"position": 0.12,
                                  "color": (200, 200, 100),
                                  "fill": (200, 200, 200, 100)}
                       )
p_x.addItem(inf1)

infline_ymin = pg.InfiniteLine(angle=0,
                               label="Min: {value:.2f} uT",
                               pen=pen_red,
                               pos=ymin,
                               movable=False,
                               labelOpts={"position": 0.05,
                                          "color": pen_red.color().getRgb(),
                                          "fill": (0, 0, 0, 0)}
                               )
p_x.addItem(infline_ymin)

infline_ymax = pg.InfiniteLine(angle=0,
                               label="Max: {value:.2f} uT",
                               pen=pen_red,
                               pos=ymax,
                               movable=False,
                               labelOpts={"position": 0.05,
                                          "color": pen_red.color().getRgb(),
                                          "fill": (0, 0, 0, 0)}
                               )
p_x.addItem(infline_ymax)

print(pen_red.color().getRgb())

# Executing PyQtGraph ========================================================
if __name__ == '__main__':
    pg.exec()