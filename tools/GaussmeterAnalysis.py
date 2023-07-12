# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.fft import fft, fftfreq
from datetime import datetime
import pyqtgraph as pg
# from pyqtgraph.Qt import QtCore
from scipy.signal import savgol_filter
from scipy.fft import rfft, rfftfreq
from copy import deepcopy
from time import time

# TODO: Remove (as it has been moved to separate file
class TimeplotPyQt_depricated:
    def __init__(self):
        self.app = pg.mkQApp("Timeplot")
        self.view = pg.GraphicsView()
        self.gl = pg.GraphicsLayout()

        self.view.setCentralItem(self.gl)

        # Default settings
        self.gl.setBorder(100, 100, 100)
        self.view.setWindowTitle("Timeplot")
        self.view.resize(800, 600)

        self.gl.layout.setColumnPreferredWidth(0, 25)
        self.gl.layout.setColumnPreferredWidth(1, 1200)

        pg.setConfigOptions(antialias=True)  # Antialiasing for nicer plots

        self.set_defaults()


    def set_defaults(self):
        self.light_theme = False

        self.padding = [1.1, 1.0]  # Padding factor
        self.force_identical_scale = True
        self.pen_alpha = 0.8

        pen_colors = ((255,   0,   0),
                      (  0, 255,   0),
                      (  0,  75, 255))
        self.generate_pens(pen_colors)

        # Colours of the draggable lines, depending on theme
        if self.light_theme:
            self.label_fill = (255, 255, 255, 180)  # White label backgrounds
            self.label_text_color = (100, 100, 100)  # Darker label text
        else:
            self.label_fill = (0, 0, 0, 180)
            self.label_text_color = (200, 200, 100)

        # Background colour:
        if self.light_theme:
            self.view.setBackground("white")  # White background
        else:
            self.view.setBackground("black")

        self.plot_title = "Plot Title"
        self.side_label = "Side label"

        self.grids = True


    def set_light_theme(self):
        self.light_theme = True

    def set_dark_theme(self):
        self.light_theme = False

    def generate_pens(self, pc):
        self.pen_rgb = pc

        self.pen_rgba = (
            (pc[0][0], pc[0][1], pc[0][2], self.pen_alpha*100),
            (pc[1][0], pc[1][1], pc[1][2], self.pen_alpha*100),
            (pc[2][0], pc[2][1], pc[2][2], self.pen_alpha*100),
        )

    def set_window_size(self, w: int, h: int):
        self.view.resize(w, h)

    def set_window_title(self, window_title: str):
        self.view.setWindowTitle(window_title)

    def set_padding(self, padding_h: float, padding_v: float):
        self.padding = [padding_h, padding_v]

    def set_force_identical_scale(self, force_identical_scale: bool):
        self.force_identical_scale = force_identical_scale

    def set_pen_alpha(self, alpha: float):
        """
        Set alpha of the pen.
        Important: set in a range between 0.0 - 1.0 !
        """
        self.pen_alpha = alpha
        generate_pens(self, self.pen_rgb)

    def set_pen_colors(self, pen_colors):
        """
        Set colors of the pen.
        Must be supplied as a list of listed rgb values, e.g.:
        pen_colors = ((255,   0,   0),
                      (  0, 255,   0),
                      (  0,   0, 255))
        """
        self.generate_pens(pen_colors)

    def set_label_fill(self, rgba):
        """
        Set a label fill RGBA, for example: (255, 255, 255, 180)
        """
        self.label_fill = rgba

    def set_label_text_color(self, rgb):
        """
        Set a label text RGB, for example: (255, 255, 255)
        """
        self.label_text_color = rgb

    def set_background(self, rgb):
        self.view.setBackground(rgb)

    def set_plot_title(self, plot_title: str):
        self.plot_title = plot_title

    def set_side_label(self, side_label: str):
        self.side_label = side_label

    def show_grids(self):
        self.grids = True

    def hide_grids(self):
        self.grids = False

    def timeplot_pyqtgraph(self, data: list):
        # TODO: Plotting stuff

        # # [DEBUG]
        # for data_list in data:
        #     print(type(data_list))
        #     print(type(data_list[0]))

        # Ensure data is of the same length, and of the right type(s)
        for data_list in (data[1], data[2], data[3]):
            assert len(data[0]) == len(data_list)
            assert type(data_list) in (list, np.ndarray)
            assert type(data_list[0]) in (float, int, np.float64)

        # Unpacking data for plotting
        # self.data_t = np.array(data[0])
        # self.data_t_date = []
        # for i in range(len(data[0])):
        #     self.data_t_date.append(self.data[0][i])
        # self.data_x = np.array(data[1])
        # self.data_y = np.array(data[2])
        # self.data_z = np.array(data[3])

        self.trange = [data[0][0], data[0][-1]]
        self.xrange = [min(data[1]), max(data[1])]
        self.yrange = [min(data[2]), max(data[2])]
        self.zrange = [min(data[3]), max(data[3])]

        # Force identical scale if set to True:
        if self.force_identical_scale is True:
            xmid = self.xrange[0] + (self.xrange[1] - self.xrange[0]) / 2
            ymid = self.yrange[0] + (self.yrange[1] - self.yrange[0]) / 2
            zmid = self.zrange[0] + (self.zrange[1] - self.zrange[0]) / 2

            xscale = self.xrange[1] - self.xrange[0]
            yscale = self.yrange[1] - self.yrange[0]
            zscale = self.zrange[1] - self.zrange[0]

            scale = max([xscale, yscale, zscale])

            self.xrange = [xmid - scale / 2, xmid + scale / 2]
            self.yrange = [ymid - scale / 2, ymid + scale / 2]
            self.zrange = [zmid - scale / 2, zmid + scale / 2]

        # Start constructing plot ============================================

        # Title and side label
        self.gl.addLabel(self.plot_title, col=1, colspan=2)
        self.gl.nextRow()
        self.gl.addLabel(self.side_label, angle=-90, rowspan=3)

        # Data plots

        plotlabels = ("X", "Y", "Z")
        p_range = (self.xrange, self.yrange, self.zrange)
        for i, data_array in enumerate((data[1], data[2], data[3])):

            # Main data plot
            p = self.gl.addPlot(
                title=plotlabels[i],
                x=data[0],
                y=np.array(data_array),
                pen=pg.mkPen(color=self.pen_rgba[i]),
                xmin=self.trange[0],
                xmax=self.trange[1],
                axisItems={'bottom': pg.DateAxisItem()}
            )
            p.setYRange(p_range[i][0], p_range[i][1])
            p.showGrid(x=self.grids, y=self.grids)

            # Horizontal nfinite line object (min)
            infline_min = pg.InfiniteLine(
                angle=0,
                label="Min: {value:.2f} uT",
                pen=pg.mkPen(color=self.pen_rgb[i]),
                pos=p_range[i][0],
                movable=False,
                labelOpts={"position": 0.05,
                           "color": self.pen_rgb[i],
                           "fill": self.label_fill}
            )
            p.addItem(infline_min)

            # Horizontal nfinite line object (max)
            infline_max = pg.InfiniteLine(
                angle=0,
                label="Max: {value:.2f} uT",
                pen=pg.mkPen(color=self.pen_rgb[i]),
                pos=p_range[i][1],
                movable=False,
                labelOpts={"position": 0.05,
                           "color": self.pen_rgb[i],
                           "fill": self.label_fill}
            )
            p.addItem(infline_max)

            # # [DEBUG]
            # print(self.label_text_color)
            # print(type(self.label_text_color))

            # Horizontal scanline
            infline_h = pg.InfiniteLine(
                angle=90,
                label="{value:.0f}",
                pos=self.trange[1],
                movable=True,
                bounds=[self.trange[0], self.trange[1]],
                labelOpts={"position": 0.12,
                           "color": self.label_text_color,
                           "fill": self.label_fill}
            )
            p.addItem(infline_h)

            # Vertical scanline
            infline_v = pg.InfiniteLine(
                angle=0,
                label="{value:.02f}",
                pos=p_range[i][1],
                movable=True,
                bounds=[p_range[i][0], p_range[i][1]],
                labelOpts={"position": 0.88,
                           "color": self.label_text_color,
                           "fill": self.label_fill}
            )
            p.addItem(infline_v)

            self.gl.nextRow()

        self.view.show()

        # Executing PyQtGraph
        pg.exec()


class EulerRotation:
    def __init__(self):
        pass

    def rx(self, angle, deg=True):
        if deg:
            angle *= np.pi / 180
        return np.array([[1,             0,              0],
                         [0, np.cos(angle), -np.sin(angle)],
                         [0, np.sin(angle),  np.cos(angle)]])

    def ry(self, angle, deg=True):
        if deg:
            angle *= np.pi / 180
        return np.array([[ np.cos(angle), 0, np.sin(angle)],
                         [             0, 1,             0],
                         [-np.sin(angle), 0, np.cos(angle)]])

    def rz(self, angle, deg=True):
        if deg:
            angle *= np.pi / 180
        return np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle),  np.cos(angle), 0],
                         [            0,              0, 1]])

    def rotateX(self, vector, angle, deg=True):
        return np.dot(self.rx(angle, deg=deg), vector)

    def rotateY(self, vector, angle, deg=True):
        return np.dot(self.ry(angle, deg=deg), vector)

    def rotateZ(self, vector, angle, deg=True):
        return np.dot(self.rz(angle, deg=deg), vector)


def data_rotate(data, rotation_matrix):
    for i in range(len(data[0])):
        vtemp = np.dot(rotation_matrix, np.array([data[1][i],
                                                  data[2][i],
                                                  data[3][i]]))
        data[1][i], data[2][i], data[3][i] = vtemp[0], vtemp[1], vtemp[2]
    return data

def data_add_vector(data, vector):
    assert len(vector) == len(data)-1
    # Loops over entire length of dataset
    for i in range(len(data[0])):
        # Loops over separate data lanes, skipping the first
        for i_data in range(len(data)-1):
            data[i_data+1][i] += vector[i_data]
    return data

def read_data(filename, header=False, len_header=9):
    with open(filename, 'r') as f:
        lines = f.readlines()
        if lines[0][0:2] == "!H" or header is True:
            lines = lines[len_header+1:]

        t, x, y, z = [0] * len(lines), [0] * len(lines), [0] * len(lines), [0] * len(lines)

        for i in range(len(lines)):
            [tt, xx, yy, zz] = lines[i].split(" ")
            [t[i], x[i], y[i], z[i]] = \
                [np.double(tt), float(xx), float(yy), float(zz)]
        return [t, x, y, z]


def data_rfft_legacy(data, sample_rate=-1):
    # Autodetermine sample rate when not given.
    # Assumes constant sample rate.
    if sample_rate == -1:
        sample_rate = round(len(data[1]) / (data[0][-1] - data[0][0]), 3)

    x_t = rfft(data[1])
    x_f = rfftfreq(len(data[1]), 1 / sample_rate)

    y_t = rfft(data[2])
    y_f = rfftfreq(len(data[2]), 1 / sample_rate)

    z_t = rfft(data[3])
    z_f = rfftfreq(len(data[3]), 1 / sample_rate)

    return [[x_t, x_f], [y_t, y_f], [z_t, z_f]]


def data_rfft(data, sample_rate=-1):
    # Autodetermine sample rate when not given.
    # Assumes constant sample rate.
    if sample_rate == -1:
        sample_rate = round(len(data[1]) / (data[0][-1] - data[0][0]), 3)
    n_samples = len(data[0])

    f = rfftfreq(n_samples, 1 / sample_rate)
    # x_t = rfft(data[1] * np.hanning(n_samples))
    # x_t = rfft(data[1])
    x_t = rfft(data[1])
    y_t = rfft(data[2])
    z_t = rfft(data[3])

    return [f, np.abs(x_t), np.abs(y_t), np.abs(z_t)]
    # return [f, x_t, y_t, z_t]


def vector_normal_distribution(vectors: list):
    mean_x = 0
    mean_y = 0
    mean_z = 0

    sd_x = 0
    sd_y = 0
    sd_z = 0

    for vector in vectors:
        mean_x += vector[0]
        mean_y += vector[1]
        mean_z += vector[2]

    mean_x = mean_x / len(vectors)
    mean_y = mean_y / len(vectors)
    mean_z = mean_z / len(vectors)

    for vector in vectors:
        sd_x += (vector[0] - mean_x) ** 2
        sd_y += (vector[1] - mean_y) ** 2
        sd_z += (vector[2] - mean_z) ** 2

    sd_x = (sd_x / len(vectors)) ** 0.5
    sd_y = (sd_y / len(vectors)) ** 0.5
    sd_z = (sd_z / len(vectors)) ** 0.5

    return [np.array([mean_x, mean_y, mean_z]), np.array([sd_x, sd_y, sd_z])]

def data_savgol_filter(data, windowlength, polyorder=4):
    # Smooths three data channels using a Savitzky-Golay filter
    # Data must be in type A
    x_filtered = savgol_filter(data[1], windowlength, polyorder)
    y_filtered = savgol_filter(data[2], windowlength, polyorder)
    z_filtered = savgol_filter(data[3], windowlength, polyorder)
    return [data[0], x_filtered, y_filtered, z_filtered]

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

# from GaussmeterLib import (ReadData,
#                            AnalyzeRate)

# def ReadData(filename, header=False):
#     with open(filename, 'r') as f:
#         lines = f.readlines()
#         if lines[0][0:2] == "!H" or header == True:
#             lines = lines[1:]
        
#         t, x, y, z = [0]*len(lines), [0]*len(lines), [0]*len(lines), [0]*len(lines)
                
#         for i in range(len(lines)):
#             [tt, xx, yy, zz] = lines[i].split(" ")
#             [t[i], x[i], y[i], z[i]] = \
#                 [np.double(tt), float(xx), float(yy), float(zz)]
#         return [t, x, y, z]



def time_plot(data, ylim=(-100, 100)):

    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(6,1)

    ax1 = fig.add_subplot(gs[0,0])

    ax1.set_axis_off()
    pos_x = 0.1
    pos_y = 0.5
    pos_z = 0.9

    ax1.text(pos_x, 0.8, str(round(data[1][len(data)],3)),
             ha="center", fontsize=26)
    ax1.text(pos_x, 0.3, "uT",
             ha="center", fontsize=16)

    ax1.text(pos_y, 0.8, str(round(data[2][len(data)],3)),
             ha="center", fontsize=26)
    ax1.text(pos_y, 0.3, "uT",
             ha="center", fontsize=16)

    ax1.text(pos_z, 0.8, str(round(data[3][len(data)],3)),
             ha="center", fontsize=26)
    ax1.text(pos_z, 0.3, "uT",
             ha="center", fontsize=14)

    ax2 = fig.add_subplot(gs[1:, 0])
    ax2.set_ylim(ylim[0], ylim[1])
    ax2.plot(data[0], data[1], color='red')
    ax2.plot(data[0], data[2], color='green')
    ax2.plot(data[0], data[3], color='blue')

    plt.show()

def analyse_rate(data, n_samples, sr, verbose=0):
    """Analyzes the time between each sample from the UNIX datapoints, and 
        performs basic statistical analysis on them."""
        
    timings = [0]*(n_samples-1)
    
    for i in range(len(data[0])-1):
        timings[i] = data[0][i+1] - data[0][i]
    
    mean = sum(timings)/len(timings)
    
    sd = 0
    for i in range(len(timings)):
        sd += (mean-timings[i])**2
    sd = sd/len(timings)
    
    return timings, mean, sd


def sampling_plot(data, sr):
    rateanalysis = analyse_rate(data, len(data[0]), sr)
    
    sampling_intervals = rateanalysis[0]
    
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(6,1)

    ax2 = fig.add_subplot(gs[1:,0])
    
    ax2.plot(list(range(len(sampling_intervals))), sampling_intervals, color='black')


def create_hhc_elements(coil_sides, coil_spacings):
    # Unpacking coil sides into X/Y/Z
    [coilX_side, coilY_side, coilZ_side] = coil_sides  # [m]
    # Unpacking coil spacings into X/Y/Z
    [coilX_spacing, coilY_spacing, coilZ_spacing] = coil_spacings

    coilXp = Face(Vertex([coilX_spacing / 2, coilX_side / 2, coilX_side / 2]),
                  Vertex([coilX_spacing / 2, -coilX_side / 2, coilX_side / 2]),
                  Vertex([coilX_spacing / 2, -coilX_side / 2, -coilX_side / 2]),
                  Vertex([coilX_spacing / 2, coilX_side / 2, -coilX_side / 2]))
    coilXn = Face(Vertex([-coilX_spacing / 2, coilX_side / 2, coilX_side / 2]),
                  Vertex([-coilX_spacing / 2, -coilX_side / 2, coilX_side / 2]),
                  Vertex([-coilX_spacing / 2, -coilX_side / 2, -coilX_side / 2]),
                  Vertex([-coilX_spacing / 2, coilX_side / 2, -coilX_side / 2]))

    coilYp = Face(Vertex([coilY_side / 2, coilY_spacing / 2, coilY_side / 2]),
                  Vertex([coilY_side / 2, coilY_spacing / 2, -coilY_side / 2]),
                  Vertex([-coilY_side / 2, coilY_spacing / 2, -coilY_side / 2]),
                  Vertex([-coilY_side / 2, coilY_spacing / 2, coilY_side / 2]))
    coilYn = Face(Vertex([coilY_side / 2, -coilY_spacing / 2, coilY_side / 2]),
                  Vertex([coilY_side / 2, -coilY_spacing / 2, -coilY_side / 2]),
                  Vertex([-coilY_side / 2, -coilY_spacing / 2, -coilY_side / 2]),
                  Vertex([-coilY_side / 2, -coilY_spacing / 2, coilY_side / 2]))

    coilZp = Face(Vertex([coilZ_side / 2, coilZ_side / 2, coilZ_spacing / 2]),
                  Vertex([-coilZ_side / 2, coilZ_side / 2, coilZ_spacing / 2]),
                  Vertex([-coilZ_side / 2, -coilZ_side / 2, coilZ_spacing / 2]),
                  Vertex([coilZ_side / 2, -coilZ_side / 2, coilZ_spacing / 2]))
    coilZn = Face(Vertex([coilZ_side / 2, coilZ_side / 2, -coilZ_spacing / 2]),
                  Vertex([-coilZ_side / 2, coilZ_side / 2, -coilZ_spacing / 2]),
                  Vertex([-coilZ_side / 2, -coilZ_side / 2, -coilZ_spacing / 2]),
                  Vertex([coilZ_side / 2, -coilZ_side / 2, -coilZ_spacing / 2]))

    return [[coilXn, coilXp], [coilYn, coilYp], [coilZn, coilZp]]


def create_lab_walls(lx=2, ly=2, lz=2, ox=0.5, oy=0.5, oz=0):
    wall_x = Face(Vertex([lx / 2 + ox, ly / 2 + oy, lz / 2 + oz]),
                  Vertex([lx / 2 + ox, -(ly / 2 + oy), lz / 2 + oz]),
                  Vertex([lx / 2 + ox, -(ly / 2 + oy), -(lz / 2 + oz)]),
                  Vertex([lx / 2 + ox, ly / 2 + oy, -(lz / 2 + oz)]))

    wall_y = Face(Vertex([-(lx / 2 + ox), ly / 2 + oy, lz / 2 + oz]),
                  Vertex([lx / 2 + ox, ly / 2 + oy, lz / 2 + oz]),
                  Vertex([lx / 2 + ox, ly / 2 + oy, -(lz / 2 + oz)]),
                  Vertex([-(lx / 2 + ox), ly / 2 + oy, -(lz / 2 + oz)]))

    return [wall_x, wall_y]


def FourierPlot(data):
    pass

# filename = "output_StressTest1_1200_100_2023-06-19_15.19.02.dat"
# filename2 = "output_StressTest2_50_100_2023-06-19_15.58.24.dat"

# data = ReadData(filename, header=True)
# TimePlot(data)



# sr = 100
# SamplingPlot(data, sr)