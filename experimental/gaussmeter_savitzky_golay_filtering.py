#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johan Monster
"""
#%% ============================== Imports =================================

from cp.cp_vertex import Vertex
from cp.cp_face import Face
# from cp_geometry import Geometry
# from cp_frame import Frame
from cp.cp_plotting import plot_face, plot_arrow, plot_global_tripod
#import cp_plotting
# from cp_utilities import d2r
from GaussmeterAnalysis import time_plot

from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d
import numpy as np
from copy import deepcopy
from scipy.signal import savgol_filter
from scipy.fft import rfft, rfftfreq

#%% =================== Create Helmholtz coil elements ======================

def Rx(theta, deg=False):
    if deg:
        theta *= np.pi/180
    Rx = np.array([[1,             0,              0],
                   [0, np.cos(theta), -np.sin(theta)],
                   [0, np.sin(theta),  np.cos(theta)]])
    return Rx

def Ry(theta, deg=False):
    if deg:
        theta *= np.pi/180
    Ry = np.array([[ np.cos(theta), 0, np.sin(theta)],
                   [             0, 1,             0],
                   [-np.sin(theta), 0, np.cos(theta)]])
    return Ry

def Rz(theta, deg=False):
    if deg:
        theta *= np.pi/180
    Rz = np.array([[np.cos(theta), -np.sin(theta), 0],
                   [np.sin(theta),  np.cos(theta), 0],
                   [            0,              0, 1]])
    return Rz


def create_HHC_elements(coil_sides,coil_spacings):
    # Unpacking coil sides into X/Y/Z
    [coilX_side, coilY_side, coilZ_side] = coil_sides # [m]
    # Unpacking coil spacings into X/Y/Z
    [coilX_spacing, coilY_spacing, coilZ_spacing] = coil_spacings
    
    coilXp = Face(Vertex([ coilX_spacing/2, coilX_side/2, coilX_side/2]),
                  Vertex([ coilX_spacing/2,-coilX_side/2, coilX_side/2]),
                  Vertex([ coilX_spacing/2,-coilX_side/2,-coilX_side/2]),
                  Vertex([ coilX_spacing/2, coilX_side/2,-coilX_side/2]))
    coilXn = Face(Vertex([-coilX_spacing/2, coilX_side/2, coilX_side/2]),
                  Vertex([-coilX_spacing/2,-coilX_side/2, coilX_side/2]),
                  Vertex([-coilX_spacing/2,-coilX_side/2,-coilX_side/2]),
                  Vertex([-coilX_spacing/2, coilX_side/2,-coilX_side/2]))
    
    coilYp = Face(Vertex([ coilY_side/2, coilY_spacing/2, coilY_side/2]),
                  Vertex([ coilY_side/2, coilY_spacing/2,-coilY_side/2]),
                  Vertex([-coilY_side/2, coilY_spacing/2,-coilY_side/2]),
                  Vertex([-coilY_side/2, coilY_spacing/2, coilY_side/2]))
    coilYn = Face(Vertex([ coilY_side/2,-coilY_spacing/2, coilY_side/2]),
                  Vertex([ coilY_side/2,-coilY_spacing/2,-coilY_side/2]),
                  Vertex([-coilY_side/2,-coilY_spacing/2,-coilY_side/2]),
                  Vertex([-coilY_side/2,-coilY_spacing/2, coilY_side/2]))
    
    coilZp = Face(Vertex([ coilZ_side/2, coilZ_side/2, coilZ_spacing/2]),
                  Vertex([-coilZ_side/2, coilZ_side/2, coilZ_spacing/2]),
                  Vertex([-coilZ_side/2,-coilZ_side/2, coilZ_spacing/2]),
                  Vertex([ coilZ_side/2,-coilZ_side/2, coilZ_spacing/2]))
    coilZn = Face(Vertex([ coilZ_side/2, coilZ_side/2,-coilZ_spacing/2]),
                  Vertex([-coilZ_side/2, coilZ_side/2,-coilZ_spacing/2]),
                  Vertex([-coilZ_side/2,-coilZ_side/2,-coilZ_spacing/2]),
                  Vertex([ coilZ_side/2,-coilZ_side/2,-coilZ_spacing/2]))
    
    return [[coilXn, coilXp], [coilYn, coilYp], [coilZn, coilZp]]

def create_lab_walls(lx=2, ly=2, lz=2, ox=0.5, oy=0.5, oz=0):

    wall_x = Face(Vertex([ lx/2+ox,   ly/2+oy ,   lz/2+oz ]),
                  Vertex([ lx/2+ox, -(ly/2+oy),   lz/2+oz ]),
                  Vertex([ lx/2+ox, -(ly/2+oy), -(lz/2+oz)]),
                  Vertex([ lx/2+ox,   ly/2+oy , -(lz/2+oz)]))
    
    wall_y = Face(Vertex([ -(lx/2+ox), ly/2+oy,   lz/2+oz ]),
                  Vertex([   lx/2+ox , ly/2+oy,   lz/2+oz ]),
                  Vertex([   lx/2+ox , ly/2+oy, -(lz/2+oz)]),
                  Vertex([ -(lx/2+ox), ly/2+oy, -(lz/2+oz)]))

    return [wall_x, wall_y]

def read_data(filename, header=False):
    with open(filename, 'r') as f:
        lines = f.readlines()
        if lines[0][0:2] == "!H" or header == True:
            lines = lines[10:]
        
        t, x, y, z = [0]*len(lines), [0]*len(lines), [0]*len(lines), [0]*len(lines)
                
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
    dataA = [[0.]*n, [0.]*n, [0.]*n, [0.]*n]
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
    
    mean_X = mean_X/len(vectors)
    mean_Y = mean_Y/len(vectors)
    mean_Z = mean_Z/len(vectors)
        
    for vector in vectors:
        sd_X += (vector[0]-mean_X)**2
        sd_Y += (vector[1]-mean_Y)**2
        sd_Z += (vector[2]-mean_Z)**2
    
    sd_X = (sd_X/len(vectors))**0.5
    sd_Y = (sd_Y/len(vectors))**0.5
    sd_Z = (sd_Z/len(vectors))**0.5
    
    return [np.array([mean_X, mean_Y, mean_Z]), np.array([sd_X, sd_Y, sd_Z])]


#%% ========================= Plotting settings ============================
       
# # Tripod properties
# show_tripod=True,     # If False, does not plot the tripod
# tripod_scale=1,       # Sets the scale of the tripod
# plot_perpendiculars=True, # Plot perpendiculars
       
# # Vertex plotting properties:
# vertexfill = True,    # If False, vertex will not be plotted
# vertexcolour="#000",  # Specifies the vertex colour
# vertexsize=10,        # Size of the plotted vertex
# vertexalpha=1,        # Opacity of the plotted vertex

# # Face plotting properties:
# linefill=True,        # If False, does not plot face lines
# linecolour="#000",    # Colour of face lines
# linewidth=2,          # Thickness of face lines
# linealpha=1,          # Opacity of face lines
# facefill=True,        # If False, does not shade the face area
# facecolour="#555",    # Colour of the face area shading
# facealpha=1,          # Opacity of the face area shading
    
# # Face perpendicular arrow plotting properties:
# perpfill = False,     # If True, plot perpendiculars
# perpcolour="#888",    # Specifies the perp. arrow colour
# perpscale=1,          # Size of the plotted perp. arrow
# perpalpha=0.5,        # Opacity of the plotted perp. arrow             
    
# # Illumination:
# illumination=False,   # If True, plots illumination intensity
# ill_value=0,          # Used to plot illumination intensity
# ill_plane=None,       # If illumination is used, a plane

       
# # Vector plotting properties:
# vectorfill=True,      # If False, does not plot vector arrow
# vectorcolour="#000",  # Colour of vector arrow
# vectoralpha=1,        # Opacity of vector arrow
# vectorscale=1,        # Scale the whole vector by a constant
# vectorratio=0.15      # Vector arrow length ratio


#%% ======================== Setting up the plot ===========================

# fig = plt.figure(figsize=(15, 10.5))
# ax = mp3d.Axes3D(fig, auto_add_to_figure=False)
# fig.add_axes(ax)
# ax.clear()

# plotscale = 1.1
# ax.set_xlim(-4/3*plotscale, 4/3*plotscale)
# ax.set_ylim(-4/3*plotscale, 4/3*plotscale)
# ax.set_zlim(-plotscale, plotscale)

# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')

# # Default camera view
# ax.view_init(elev=20, azim=-50)


#%% ======================== Plotting lab objects ========================== 


# # Plotting the global XYZ tripod
# plot_global_tripod(ax, scaling=plotscale/2)


# # Create Helmholtz Cage coil objects.
# # Returns: Face objects, ordered like: [ [X- X+] ,[Y- Y+], [Z- Z+] ]
# coil_sides = [1.85, 1.95, 2.05]                 # [m]
# coil_spacings = [1.0073, 1.0618, 1.1162]        # [m]

# cage_objects = create_HHC_elements(coil_sides, coil_spacings)

# # Plot the Helmholtz cage coil pairs:
# for coilX in cage_objects[0]:
#     plot_face(ax, coilX, linecolour="#F00", linewidth=5, linealpha=0.5, 
#               facefill=False, vertexfill=False)
# for coilY in cage_objects[1]:
#     plot_face(ax, coilY, linecolour="#0F0", linewidth=5, linealpha=0.5,
#               facefill=False, vertexfill=False)
# for coilZ in cage_objects[2]:
#     plot_face(ax, coilZ, linecolour="#00F", linewidth=5, linealpha=0.5,
#               facefill=False, vertexfill=False)
    
# # Plot the lab walls
# walls = create_lab_walls(lx=2, ly=2, lz=2, ox=0.5, oy=0.5, oz=0.1)
# for wall in walls:
#     plot_face(ax, wall, linecolour="#333", linewidth=3, linealpha=0.5,
#               facefill=True, facealpha=0.3, vertexfill=False)

    
#%% =========================== Add EMF vector =============================
R_G2C = Rz(-(90+23), deg=True)
# R_G2C = Rz(-(90+23), deg=True)
R_S2C = np.dot(Rz(-90, deg=True), Rx(-180, deg=True))
# R_S2C = np.dot(Rz(-90, deg=True), Rx(-180, deg=True))

# EMF vector at AE building (WMM-2020 model)
B_EMF_WMM2020  = np.array([0.7216, 19.1629, -45.4592])
# EMF vector at AE building (IGRF2020 model)
B_EMF_IGRF2020 = np.array([0.7313, 19.1870, -45.4557])

# Take average of the two as baseline B_EMF:
B_EMF_G = 0.5 * (B_EMF_WMM2020 + B_EMF_IGRF2020)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
B_EMF = np.dot(R_G2C, B_EMF_G)
print("C-frame EMF vector:", B_EMF.round(3))

# plot_arrow(ax, [0,0,0], B_EMF, alpha=1, scaling=0.02, alr=0.1, color="orange")


#%% =========================== Add Measured =============================

# filename = "WeekendTest_2023-06-30_15.07.51_SPARSE50.dat"
# filename = "WeekdayTest_2023-07-03_10.41.03.dat"
# filename = "300HzTest_2023-07-05_10.24.48.dat"
filename1 = "2500HzTestLightOn_2023-07-05_11.04.28.dat"
filename2 = "2500HzTestLightOff_2023-07-05_11.13.36.dat"

data1 = read_data(filename1, header=True)
data2 = read_data(filename2, header=True)

print("Dataset contains",len(data1[0]),"data points.")

# mG to uT conversion:
# for i in range(1, len(data_r)):
#     for j in range(len(data_r[0])):
#         data_r[i][j] *= 0.1

# time_plot(np.array(data_r), ylim=[-50, 50], 
          # title="Weekday test - Measurements in cage frame")

def data_savgol_filter(data, windowlength, polyorder=4):
    # Smooths three data channels using a Savitzky-Golay filter
    # Data must be in type A
    x_filtered = savgol_filter(data[1], windowlength, polyorder)
    y_filtered = savgol_filter(data[2], windowlength, polyorder)
    z_filtered = savgol_filter(data[3], windowlength, polyorder)
    return [data[0], x_filtered, y_filtered, z_filtered]


def data_rfft(data, sample_rate=-1):
    # Autodetermine sample rate when not given.
    # Assumes constant sample rate.
    if sample_rate == -1:
        sample_rate = round(len(data[1])/(data[0][-1]-data[0][0]),3)
    
    x_t = rfft(data[1])
    x_f = rfftfreq(len(data[1]), 1/sample_rate)

    y_t = rfft(data[2])
    y_f = rfftfreq(len(data[2]), 1/sample_rate)

    z_t = rfft(data[3])
    z_f = rfftfreq(len(data[3]), 1/sample_rate)
    
    return [[x_t, x_f], [y_t, y_f], [z_t, z_f]]


def data_rotate(dataA, R):
    dataB = []
    for i in range(len(dataA[0])):
        dataB.append([
            dataA[0][i], 
            np.dot(R, np.array([dataA[1][i], dataA[2][i], dataA[3][i]]))
            ])
    return zipBA(dataB)

#%% 

# Rotate data to cage frame (C):
data1 = data_rotate(data1, R_S2C)
data2 = data_rotate(data2, R_S2C)

# Smooth data with a Savitzky-Golay filter
data1 = data_savgol_filter(data1, 51, 4)
data2 = data_savgol_filter(data2, 51, 4)

[[x_t1, x_f1], [y_t1, y_f1], [z_t1, z_f1]] = data_rfft(data1)
[[x_t2, x_f2], [y_t2, y_f2], [z_t2, z_f2]] = data_rfft(data2)


#%% Plot time domain:
time_plot(np.array(data1), ylim=[-50, 50], 
           title="Weekday test - Measurements in cage frame")

#%% Plot Spectrogram LightOn XYZ

fig1 = plt.figure(figsize=(8, 8))
gs = GridSpec(1,1)
ax = fig1.add_subplot(gs[0,0])
# ax.set_ylim(ylim[0], ylim[1])

filter_dc: int = 3

ax.plot(x_f1[filter_dc:], np.abs(x_t1[filter_dc:]), color='r', label="X")
ax.plot(y_f1[filter_dc:], np.abs(y_t1[filter_dc:]), color='g', label="Y")
ax.plot(z_f1[filter_dc:], np.abs(z_t1[filter_dc:]), color='b', label="Z")

ax.set_title("Spectrogram of magnetic measurements")
ax.set_xlabel("f [Hz]")
ax.set_ylabel("Linear magnitude")
ax.legend()
ax.grid(True)


#%% Plot Spectrograms LightOn/LightOff


def plot_spectrogram_comparison(x_on, x_off, c1="red", c2="black", 
                                label1 = "on", label2 = "off",
                                title="Spectrogram", filter_dc=3):
    [x_on_t, x_on_f] = x_on
    [x_off_t, x_off_f] = x_off   
    
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(1,1)
    ax = fig.add_subplot(gs[0,0])
    # ax.set_ylim(ylim[0], ylim[1])

    ax.plot(x_on_f[filter_dc:], np.abs(x_on_t[filter_dc:]), 
            color=c1, label=label1, linewidth=2)
    ax.plot(x_off_f[filter_dc:], np.abs(x_off_t[filter_dc:]), 
            color=c2, label=label2, linewidth=2)

    ax.set_title(title)
    ax.set_xlabel("f [Hz]")
    ax.set_ylabel("Linear magnitude")
    ax.legend()
    ax.grid(True)

filter_dc = 3

plot_spectrogram_comparison([x_t1, x_f1], [x_t2, x_f2], c1="r", filter_dc=filter_dc)
plot_spectrogram_comparison([y_t1, y_f1], [y_t2, y_f2], c1="g", filter_dc=filter_dc)
plot_spectrogram_comparison([z_t1, z_f1], [z_t2, z_f2], c1="b", filter_dc=filter_dc)

# #%% Vector plot

# data1B = []
# for i in range(len(data1[0])):
#     data1B.append([
#         data1[0][i], 
#         np.array([data1[1][i], data1[2][i], data1[3][i]])
#         ])
# data1_vectors = []
# for d in data1B:
#     data1_vectors.append(d[1])

# mean_data1, sd_data1 = vector_normal_distribution(data1_vectors)
# # print("C-frame data: mean, sd:", mean_data1.round(3), sd_data1.round(3))

# # Plot only average measured vector:
# plot_arrow(ax, [0,0,0], list(mean_data1), alpha=1, scaling=0.02, alr=0.1, 
#            color="purple")

# # time_plot(zipBA(data), ylim=[-50, 50])

# # Normalize the data (in C frame) with respect to the EMF (in C frame):
# data1_n = deepcopy(data1)
    
# for i in range(len(data1_n)):
#     data1_n[i][1] = data1_n[i][1] - B_EMF

# time_plot(zipBA(data1_n), ylim=[-15, 20], 
#           title="Rack Elimination Test - Normalized disturbance vector")

# data1_n_vectors = []
# for d in data1_n:
#     data1_n_vectors.append(d[1])
    
# mean_data1_n, sd_data1_n = vector_normal_distribution(data1_n_vectors)
# # print("Normalized data: mean, sd:", mean_data1_n.round(3), sd_data1_n.round(3))
# # print("Absolute strength of disturbance:", 
# #       round(np.linalg.norm(mean_data1_n),3),
# #       "uT (", 
# #       round(100*(np.linalg.norm(mean_data1_n)/np.linalg.norm(B_EMF)),2),
# #       "% of EMF)")

# plot_arrow(ax, [0,0,0], mean_data1_n, linewidth=3, alpha=1, scaling=0.02, 
#            alr=0.3, color="black")
