# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import mpl_toolkits.mplot3d as mp3d
from scipy.fft import fft, fftfreq
from datetime import datetime
import pyqtgraph as pg
# from pyqtgraph.Qt import QtCore
from scipy.signal import savgol_filter
from scipy.fft import rfft, rfftfreq
from copy import deepcopy
from time import time

from .local_emf import local_emf

from .cp import (
    Vertex,
    Face,
    Geometry,
    Frame,
    plot_face,
    plot_arrow,
    plot_global_tripod,
    plot_frame
)

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

def multipart_filenames(filename_base, n_files):
    """
    Generates filenames for multipart file splits. This function is here so
     that it can both be used for generating these names in other functions,
     and so it can be used by a user to generate lists of filenames with the
     exact same algorithm (to reduce the risk of bugs).
    """
    filenames = []
    for i in range(1, n_files+1):
        filenames.append(
            filename_base[:-4]
            + "_s{}of{}".format(i, n_files)
            + filename_base[-4:]
        )
    return filenames

def read_data(filename, header=False, len_header=9):
    """
    Reads data from data file and formats it to a format that other functions
     can use.
    """
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

def write_data(filename, data, mode="a", verbose=0):
    with open(filename, mode) as fo:
        for i in range(len(data[0])):
            fo.write("{} {} {} {}\n".format(data[0][i],
                                            data[1][i],
                                            data[2][i],
                                            data[3][i]))
    if verbose >= 2:
        print("Appended {} entries to {}".format(len(data[0]), filename))

def data_merge(filenames, filename_merge: str = None, len_header=9, verbose=0):
    """
    Takes a list of filenames of data files to be merged. Doing so will create
     a merge file which copies the header of the first file, and appends the
     data of all files sequentially to this file.

    Use filename_merge to specify a name for the merge file. If unspecified,
     it will be based on the first filename in filenames.

    This function is not optimized for >RAM datasets.
    """

    # If unspecified, base filename of merged file on first file
    if filename_merge is None:
        filename_merge = filenames[0][:-4]+"_MERGED"+filenames[0][-4:]

    # Create empty merged file
    with open(filename_merge, 'x') as fo:
        if verbose >= 2:
            print("Created merge file {}".format(filename_merge))

    # Copy header from first file
    with open(filenames[0], 'r') as f, open(filename_merge, 'a') as fo:
        for i in range(len_header-1):
            fo.write(f.readline())

    # Sequentially append contents of all files in filenames list to the
    # merged file
    for filename in filenames:
        write_data(filename_merge,
                   read_data(filename,
                             header=True,
                             len_header=len_header-2)
                   )
        if verbose >= 2:
            print("Appended contents of {} to merge file".format(filename))

def data_downsample(filenames, downsampling_factor,
                    len_header=9, verbose=0):
    """
    Function that takes a list of filenames, and downsamples them by an integer
     factor.

     This function is not meant for >RAM files, and the whole dataset will be
     read into memory when processing. To downsample a >RAM dataset, first cut
     it into segments using the DataProcessor class before downsampling.
     The downsampled segments can then be merged using the data_merge function.
    """

    # Loop over all filenames in filenames object
    for filename in filenames:

        # Load data
        data = read_data(filename, header=True)

        # Determine the number of chunks in which data can be subdivided:
        chunks, _ = divmod(len(data[0]), downsampling_factor)

        # Pre-allocate object for downsampled data
        data_sparse = [[0.]*chunks, [0.]*chunks, [0.]*chunks, [0.]*chunks]

        # Every <downsampling_factor> sample the data and save to the new set
        i_chunk = 0
        while i_chunk < chunks:
            for v in range(len(data)):
                data_sparse[v][i_chunk] = data[v][0+i_chunk*downsampling_factor]
            i_chunk += 1

        # Generate filename for new file
        flag = "_SPARSE" + str(downsampling_factor)
        filename_sparse = filename[:-4] + flag + filename[-4:]

        # Create empty merged file
        with open(filename_sparse, 'x') as fo:
            if verbose >= 2:
                print("Created merge file {}".format(filename_sparse))

        # Copy header from first file
        with open(filename, 'r') as f, open(filename_sparse, "a") as fo:
            for i in range(len_header):
                fo.write(f.readline())
            if verbose >= 2:
                print("Copied header from {}".format(filename))

        write_data(filename_sparse, data_sparse, mode="a", verbose=verbose)

def data_rfft_legacy(data, sample_rate=None):  # TODO: Remove
    # Autodetermine sample rate when not given.
    # Assumes constant sample rate.
    if sample_rate is None:
        sample_rate = round(len(data[1]) / (data[0][-1] - data[0][0]), 3)

    x_t = rfft(data[1])
    x_f = rfftfreq(len(data[1]), 1 / sample_rate)

    y_t = rfft(data[2])
    y_f = rfftfreq(len(data[2]), 1 / sample_rate)

    z_t = rfft(data[3])
    z_f = rfftfreq(len(data[3]), 1 / sample_rate)

    return [[x_t, x_f], [y_t, y_f], [z_t, z_f]]

def data_rfft(data, sample_rate=None):
    """
    Performs a one-sided (non-complex) Fast Fourier Transform on the data.
     If the sample rate is not specified via sample_rate, it will be
     auto-determined.
    """
    # Autodetermine sample rate when not given.
    # Assumes constant sample rate.
    if sample_rate is None:
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


# def vector_normal_distribution(vectors: list):
#     mean_x = 0
#     mean_y = 0
#     mean_z = 0
#
#     sd_x = 0
#     sd_y = 0
#     sd_z = 0
#
#     for vector in vectors:
#         mean_x += vector[0]
#         mean_y += vector[1]
#         mean_z += vector[2]
#
#     mean_x = mean_x / len(vectors)
#     mean_y = mean_y / len(vectors)
#     mean_z = mean_z / len(vectors)
#
#     for vector in vectors:
#         sd_x += (vector[0] - mean_x) ** 2
#         sd_y += (vector[1] - mean_y) ** 2
#         sd_z += (vector[2] - mean_z) ** 2
#
#     sd_x = (sd_x / len(vectors)) ** 0.5
#     sd_y = (sd_y / len(vectors)) ** 0.5
#     sd_z = (sd_z / len(vectors)) ** 0.5
#
#     return [np.array([mean_x, mean_y, mean_z]), np.array([sd_x, sd_y, sd_z])]

def data_normal_distribution(data: list):
    l_data = len(data[0])

    mean_x = sum(data[1])/l_data
    mean_y = sum(data[2])/l_data
    mean_z = sum(data[3])/l_data

    sd_x = 0
    sd_y = 0
    sd_z = 0

    for i in range(l_data):
        sd_x += (data[1][i] - mean_x) ** 2
        sd_y += (data[2][i] - mean_y) ** 2
        sd_z += (data[3][i] - mean_z) ** 2

    sd_x = (sd_x / l_data) ** 0.5
    sd_y = (sd_y / l_data) ** 0.5
    sd_z = (sd_z / l_data) ** 0.5

    return [np.array([mean_x, mean_y, mean_z]), np.array([sd_x, sd_y, sd_z])]

def data_savgol_filter(data, windowlength, polyorder=4):
    # Smooths three data channels using a Savitzky-Golay filter
    # Data must be in type A
    x_filtered = savgol_filter(data[1], windowlength, polyorder)
    y_filtered = savgol_filter(data[2], windowlength, polyorder)
    z_filtered = savgol_filter(data[3], windowlength, polyorder)
    return [data[0], x_filtered, y_filtered, z_filtered]

# def zipBA(dataB):
#     # Zips data from type B to type A
#     # Type A: [
#     #   [t_0 ... t_n],
#     #   [Bx_0 ... Bx_n],
#     #   [By_0 ... By_n],
#     #   [Bz_0 ... Bz_n],
#     # ]
#     #
#     # Type B: [
#     #   [t_0, np.array(Bx_0, By_0, Bz_0)] ... [t_n, np.array(Bx_n, By_n, Bz_n)]
#     # ]
#     n = len(dataB)
#     dataA = [[0.] * n, [0.] * n, [0.] * n, [0.] * n]
#     for i in range(len(dataB)):
#         dataA[0][i] = dataB[i][0]
#         dataA[1][i] = dataB[i][1][0]
#         dataA[2][i] = dataB[i][1][1]
#         dataA[3][i] = dataB[i][1][2]
#     return dataA

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


# def create_hhc_elements(coil_sides, coil_spacings):
#     # Unpacking coil sides into X/Y/Z
#     [coilX_side, coilY_side, coilZ_side] = coil_sides  # [m]
#     # Unpacking coil spacings into X/Y/Z
#     [coilX_spacing, coilY_spacing, coilZ_spacing] = coil_spacings
#
#     coilXp = Face(Vertex([coilX_spacing / 2, coilX_side / 2, coilX_side / 2]),
#                   Vertex([coilX_spacing / 2, -coilX_side / 2, coilX_side / 2]),
#                   Vertex([coilX_spacing / 2, -coilX_side / 2, -coilX_side / 2]),
#                   Vertex([coilX_spacing / 2, coilX_side / 2, -coilX_side / 2]))
#     coilXn = Face(Vertex([-coilX_spacing / 2, coilX_side / 2, coilX_side / 2]),
#                   Vertex([-coilX_spacing / 2, -coilX_side / 2, coilX_side / 2]),
#                   Vertex([-coilX_spacing / 2, -coilX_side / 2, -coilX_side / 2]),
#                   Vertex([-coilX_spacing / 2, coilX_side / 2, -coilX_side / 2]))
#
#     coilYp = Face(Vertex([coilY_side / 2, coilY_spacing / 2, coilY_side / 2]),
#                   Vertex([coilY_side / 2, coilY_spacing / 2, -coilY_side / 2]),
#                   Vertex([-coilY_side / 2, coilY_spacing / 2, -coilY_side / 2]),
#                   Vertex([-coilY_side / 2, coilY_spacing / 2, coilY_side / 2]))
#     coilYn = Face(Vertex([coilY_side / 2, -coilY_spacing / 2, coilY_side / 2]),
#                   Vertex([coilY_side / 2, -coilY_spacing / 2, -coilY_side / 2]),
#                   Vertex([-coilY_side / 2, -coilY_spacing / 2, -coilY_side / 2]),
#                   Vertex([-coilY_side / 2, -coilY_spacing / 2, coilY_side / 2]))
#
#     coilZp = Face(Vertex([coilZ_side / 2, coilZ_side / 2, coilZ_spacing / 2]),
#                   Vertex([-coilZ_side / 2, coilZ_side / 2, coilZ_spacing / 2]),
#                   Vertex([-coilZ_side / 2, -coilZ_side / 2, coilZ_spacing / 2]),
#                   Vertex([coilZ_side / 2, -coilZ_side / 2, coilZ_spacing / 2]))
#     coilZn = Face(Vertex([coilZ_side / 2, coilZ_side / 2, -coilZ_spacing / 2]),
#                   Vertex([-coilZ_side / 2, coilZ_side / 2, -coilZ_spacing / 2]),
#                   Vertex([-coilZ_side / 2, -coilZ_side / 2, -coilZ_spacing / 2]),
#                   Vertex([coilZ_side / 2, -coilZ_side / 2, -coilZ_spacing / 2]))
#
#     return [[coilXn, coilXp], [coilYn, coilYp], [coilZn, coilZp]]
#
# def create_cuboid(lx, ly, lz, ox=0, oy=0, oz=0):
#     # Define vertices
#     p1 = Vertex([ lx/2,  ly/2,  lz/2])      # .   7___ 6
#     p2 = Vertex([-lx/2,  ly/2,  lz/2])      # .   |\8__\5
#     p3 = Vertex([-lx/2, -ly/2,  lz/2])      # .   | |  |        Y
#     p4 = Vertex([ lx/2, -ly/2,  lz/2])      # .  3| | 2| -------->
#     p5 = Vertex([ lx/2,  ly/2, -lz/2])      # .    \|__|
#     p6 = Vertex([-lx/2,  ly/2, -lz/2])      # .    4  \ 1
#     p7 = Vertex([-lx/2, -ly/2, -lz/2])      # .        \
#     p8 = Vertex([ lx/2, -ly/2, -lz/2])      # .         v X
#
#     # Define faces
#     fA = Face(p1, p4, p3, p2)  # .  E ___
#     fB = Face(p3, p4, p8, p7)  # .   |\ F_\
#     fC = Face(p4, p1, p5, p8)  # . B | |  | D        Y
#     fD = Face(p1, p2, p6, p5)  # .   | | C|  -------->
#     fE = Face(p2, p3, p7, p6)  # .    \|__|
#     fF = Face(p5, p6, p7, p8)  # .      A
#
#     # Assembling geometry
#     cuboid = Geometry([fA, fB, fC, fD, fE, fF])
#
#     # frame = Frame()
#     # # Attaching the cuboid Geometry to 'frame1'
#     # frame.add_geometry(cuboid)
#     # # Translate 'frame1' away from the origin.
#     # frame.translate(ox, oy, oz)
#     return cuboid
#
# def create_table(leg_height=1.23, leg_width=0.04, leg_spacing=0.52,
#                  top_side=0.61, top_height=0.02, ox=0, oy=0, oz=-1.25):
#
#     frame_legs = Frame()
#     frame_top = Frame()
#
#     geometries_legs = []
#
#     # leg1 = create_cuboid(
#     #     leg_width, leg_width, leg_height,
#     #     1*leg_spacing/2, 1*leg_spacing/2, leg_height/2
#     # )
#
#     # leg1.translate(1*leg_spacing/2, 1*leg_spacing/2, oz+leg_height/2)
#     # frame_legs.add_geometry(leg1)
#
#     for i in ((1, 1), (-1, 1), (-1, -1), (1, -1)):
#         leg = create_cuboid(leg_width, leg_width, leg_height)
#         leg.translate(i[0]*leg_spacing/2, i[1]*leg_spacing/2, leg_height/2)
#         frame_legs.add_geometry(leg)
#
#     top = create_cuboid(top_side, top_side, top_height)
#     top.translate(0, 0, leg_height+top_height/2)
#     frame_top.add_geometry(top)
#
#     for frame in (frame_legs, frame_top):
#         frame.translate(ox, oy, oz)
#
#     return frame_legs, frame_top
#
# def create_x_wall(ly=2., lz=2., ox=0., oy=0., oz=0.):
#     wall_x = Face(Vertex([ox,  ly / 2 + oy,  lz / 2 + oz]),
#                   Vertex([ox, -ly / 2 + oy,  lz / 2 + oz]),
#                   Vertex([ox, -ly / 2 + oy, -lz / 2 + oz]),
#                   Vertex([ox,  ly / 2 + oy, -lz / 2 + oz]))
#     return wall_x
#
# def create_y_wall(lx=2., lz=2., ox=0., oy=0., oz=0.):
#     wall_y = Face(Vertex([-lx / 2 + ox, oy,  lz / 2 + oz]),
#                   Vertex([ lx / 2 + ox, oy,  lz / 2 + oz]),
#                   Vertex([ lx / 2 + ox, oy, -lz / 2 + oz]),
#                   Vertex([-lx / 2 + ox, oy, -lz / 2 + oz]))
#     return wall_y
#
# def create_floor(ly=2, lz=2, ox=0, oy=0, oz=0):
#     pass
#     # wall_y = Face(Vertex([-(lx / 2 + ox), ly / 2 + oy, lz / 2 + oz]),
#     #               Vertex([lx / 2 + ox, ly / 2 + oy, lz / 2 + oz]),
#     #               Vertex([lx / 2 + ox, ly / 2 + oy, -(lz / 2 + oz)]),
#     #               Vertex([-(lx / 2 + ox), ly / 2 + oy, -(lz / 2 + oz)]))
#     # return wall_y
#
#
# class VectorPlot:
#     def __init__(self):
#         # Define matplotlib plot structure
#         self.fig = plt.figure(figsize=(15, 10.5))
#         self.ax = mp3d.Axes3D(self.fig, auto_add_to_figure=False)
#         self.fig.add_axes(self.ax)
#         self.ax.clear()
#
#         self.default_plotting_settings()
#
#         self.plotscale = 1.1
#         self.ax.set_xlim(-4 / 3 * self.plotscale, 4 / 3 * self.plotscale)
#         self.ax.set_ylim(-4 / 3 * self.plotscale, 4 / 3 * self.plotscale)
#         self.ax.set_zlim(-self.plotscale, self.plotscale)
#
#         self.ax.set_xlabel('x')
#         self.ax.set_ylabel('y')
#         self.ax.set_zlabel('z')
#
#         # Default camera view
#         self.ax.view_init(elev=20, azim=-50)
#
#     def show(self):
#         plt.show()
#
#     def set_plot_title(self, plot_title: str):
#         self.ax.set_title(plot_title)
#
#     def default_plotting_settings(self):
#         self.plot_settings = {
#             # Tripod properties
#             "show_tripod": True,  # If False, does not plot the tripod
#             "tripod_scale": 1,  # Sets the scale of the tripod
#             "plot_perpendiculars": True,  # Plot perpendiculars
#
#             # Vertex plotting properties:
#             "vertexfill": True,  # If False, vertex will not be plotted
#             "vertexcolour": "#000",  # Specifies the vertex colour
#             "vertexsize": 10,  # Size of the plotted vertex
#             "vertexalpha": 1,  # Opacity of the plotted vertex
#
#             # Face plotting properties:
#             "linefill": True,  # If False, does not plot face lines
#             "linecolour": "#000",  # Colour of face lines
#             "linewidth": 2,  # Thickness of face lines
#             "linealpha": 1,  # Opacity of face lines
#             "facefill": True,  # If False, does not shade the face area
#             "facecolour": "#555",  # Colour of the face area shading
#             "facealpha": 1,  # Opacity of the face area shading
#
#             # Face perpendicular arrow plotting properties:
#             "perpfill": False,  # If True, plot perpendiculars
#             "perpcolour": "#888",  # Specifies the perp. arrow colour
#             "perpscale": 1,  # Size of the plotted perp. arrow
#             "perpalpha": 0.5,  # Opacity of the plotted perp. arrow
#
#             # Illumination:
#             "illumination": False,  # If True, plots illumination intensity
#             "ill_value": 0,  # Used to plot illumination intensity
#             "ill_plane": None,  # If illumination is used, a plane
#
#             # Vector plotting properties:
#             "vectorfill": True,  # If False, does not plot vector arrow
#             "vectorcolour": "#000",  # Colour of vector arrow
#             "vectoralpha": 1,  # Opacity of vector arrow
#             "vectorscale": 1,  # Scale the whole vector by a constant
#             "vectorratio": 0.15  # Vector arrow length ratio
#         }
#
#     def plot_vector(self, tip, origin=(0, 0, 0), linewidth=1.5,
#                     alpha=1, scaling=0.02, alr=0.1, color="purple"):
#         plot_arrow(self.ax, origin, tip, linewidth=linewidth,
#                    alpha=alpha, scaling=scaling, alr=alr, color=color)
#
#     def plot_global_tripod(self, scaling: float = None):
#         if scaling is None:
#             scaling = self.plotscale / 2
#         plot_global_tripod(self.ax, scaling=scaling)
#
#     def plot_hhc_coils(self, coil_sides: list, coil_spacings: list,
#                        coil_thickness: int = 5, coil_alpha: float = 0.5,
#                        coil_colors=("#F00", "#0F0", "#00F")):
#
#         self.coil_sides = coil_sides
#         self.coil_spacings = coil_spacings
#
#         # First make the cage elements
#         cage_objects = create_hhc_elements(coil_sides, coil_spacings)
#
#         for coilX in cage_objects[0]:
#             plot_face(self.ax, coilX, linecolour=coil_colors[0],
#                       linewidth=coil_thickness, linealpha=coil_alpha,
#                       facefill=False, vertexfill=False)
#         for coilY in cage_objects[1]:
#             plot_face(self.ax, coilY, linecolour=coil_colors[1],
#                       linewidth=coil_thickness, linealpha=coil_alpha,
#                       facefill=False, vertexfill=False)
#         for coilZ in cage_objects[2]:
#             plot_face(self.ax, coilZ, linecolour=coil_colors[2],
#                       linewidth=coil_thickness, linealpha=coil_alpha,
#                       facefill=False, vertexfill=False)
#
#     # Plot the lab walls
#     def plot_x_wall(self, ly=2., lz=2., ox=0., oy=0., oz=0.,
#                     linecolour="#333", linewidth=3, linealpha=0.4,
#                     facefill=True, facealpha=0.25):
#         face = create_x_wall(ly=ly, lz=lz, ox=ox, oy=oy, oz=oz)
#         plot_face(self.ax, face,
#                   linecolour=linecolour, linewidth=linewidth,
#                   linealpha=linealpha, facefill=facefill, facealpha=facealpha,
#                   vertexfill=False)
#
#     def plot_y_wall(self, lx=2., lz=2., ox=0., oy=0., oz=0.,
#                     linecolour="#333", linewidth=3, linealpha=0.4,
#                     facefill=True, facealpha=0.25):
#         face = create_y_wall(lx=lx, lz=lz, ox=ox, oy=oy, oz=oz)
#         plot_face(self.ax, face,
#                   linecolour=linecolour, linewidth=linewidth,
#                   linealpha=linealpha, facefill=facefill, facealpha=facealpha,
#                   vertexfill=False)
#
#     # # TODO: Delete
#     # def plot_emf(self, emf=None, linewidth=3,
#     #              alpha=1, scaling=0.02, alr=0.2, color="orange"):
#     #     if emf is None:
#     #         emf = local_emf()
#     #     self.plot_vector(emf, linewidth=linewidth,
#     #                      alpha=alpha, scaling=scaling, alr=alr, color=color)
#
#     def plot_table(self):
#         frame_legs, frame_top = create_table()
#         plot_frame(self.ax, frame_legs,
#                    show_tripod=False,
#                    vertexfill=False,
#                    vertexalpha=0,
#                    linecolour="#CCC",
#                    linealpha=0.3,
#                    facecolour="#CCC",
#                    facealpha=0.2,
#                    )
#         plot_frame(self.ax, frame_top,
#                    show_tripod=False,
#                    vertexfill=False,
#                    vertexalpha=0,
#                    linecolour="#C95",
#                    linealpha=0.3,
#                    facecolour="#C95",
#                    facealpha=0.2,
#                    )
#
#     def autoplot(self,
#                  tripod=True,
#                  coils=False,
#                  walls=False,
#                  table=False,
#                  # emf_vector=False,
#                  delfipq=False,
#                  cubesat3u=False,
#                  cubesat12u=False):
#         if tripod:
#             self.plot_global_tripod()
#         if coils:
#             coil_sides = [1.85, 1.95, 2.05]
#             coil_spacings = [1.0073, 1.0618, 1.1162]
#             self.plot_hhc_coils(coil_sides, coil_spacings)
#         if walls:
#             self.plot_x_wall(3, 2, -1.5, 0, -0.1)
#             self.plot_y_wall(3, 2, 0, 1.5, -0.1)
#         if table:
#             self.plot_table()
#         # if emf_vector:
#         #     self.plot_emf()
#
#
