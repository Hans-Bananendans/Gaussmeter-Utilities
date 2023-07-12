#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johan Monster
"""
#%% ============================== Imports =================================

from cp import (
    Vertex,
    Face,
    plot_face,
    plot_arrow,
    plot_global_tripod,
)
from GaussmeterAnalysis import time_plot

from local_emf import local_emf
from TimeplotPyQt import TimeplotPyQt
from SpectralplotPyQt import SpectralplotPyQt
from GaussmeterAnalysis import (
    EulerRotation,
    data_rfft,
    read_data,
    vector_normal_distribution,
    data_savgol_filter,
    data_rotate,
    data_add_vector,
    time_plot,
    create_hhc_elements,
    create_lab_walls,
)

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d
import numpy as np
from copy import deepcopy


#%% =================== Create Helmholtz coil elements ======================


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


#%% ========================= Plotting settings ============================
       
# Tripod properties
show_tripod=True,     # If False, does not plot the tripod
tripod_scale=1,       # Sets the scale of the tripod
plot_perpendiculars=True, # Plot perpendiculars
       
# Vertex plotting properties:
vertexfill = True,    # If False, vertex will not be plotted
vertexcolour="#000",  # Specifies the vertex colour
vertexsize=10,        # Size of the plotted vertex
vertexalpha=1,        # Opacity of the plotted vertex

# Face plotting properties:
linefill=True,        # If False, does not plot face lines
linecolour="#000",    # Colour of face lines
linewidth=2,          # Thickness of face lines
linealpha=1,          # Opacity of face lines
facefill=True,        # If False, does not shade the face area
facecolour="#555",    # Colour of the face area shading
facealpha=1,          # Opacity of the face area shading
    
# Face perpendicular arrow plotting properties:
perpfill = False,     # If True, plot perpendiculars
perpcolour="#888",    # Specifies the perp. arrow colour
perpscale=1,          # Size of the plotted perp. arrow
perpalpha=0.5,        # Opacity of the plotted perp. arrow             
    
# Illumination:
illumination=False,   # If True, plots illumination intensity
ill_value=0,          # Used to plot illumination intensity
ill_plane=None,       # If illumination is used, a plane

       
# Vector plotting properties:
vectorfill=True,      # If False, does not plot vector arrow
vectorcolour="#000",  # Colour of vector arrow
vectoralpha=1,        # Opacity of vector arrow
vectorscale=1,        # Scale the whole vector by a constant
vectorratio=0.15      # Vector arrow length ratio


#%% ======================== Setting up the plot ===========================

fig = plt.figure(figsize=(15, 10.5))
ax = mp3d.Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)
ax.clear()

plotscale = 1.1
ax.set_xlim(-4/3*plotscale, 4/3*plotscale)
ax.set_ylim(-4/3*plotscale, 4/3*plotscale)
ax.set_zlim(-plotscale, plotscale)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# Default camera view
ax.view_init(elev=20, azim=-50)


#%% ======================== Plotting lab objects ========================== 


# Plotting the global XYZ tripod
plot_global_tripod(ax, scaling=plotscale/2)


# Create Helmholtz Cage coil objects.
# Returns: Face objects, ordered like: [ [X- X+] ,[Y- Y+], [Z- Z+] ]
coil_sides = [1.85, 1.95, 2.05]                 # [m]
coil_spacings = [1.0073, 1.0618, 1.1162]        # [m]

cage_objects = create_hhc_elements(coil_sides, coil_spacings)

# Plot the Helmholtz cage coil pairs:
for coilX in cage_objects[0]:
    plot_face(ax, coilX, linecolour="#F00", linewidth=5, linealpha=0.5, 
              facefill=False, vertexfill=False)
for coilY in cage_objects[1]:
    plot_face(ax, coilY, linecolour="#0F0", linewidth=5, linealpha=0.5,
              facefill=False, vertexfill=False)
for coilZ in cage_objects[2]:
    plot_face(ax, coilZ, linecolour="#00F", linewidth=5, linealpha=0.5,
              facefill=False, vertexfill=False)
    
# Plot the lab walls
walls = create_lab_walls(lx=2, ly=2, lz=2, ox=0.5, oy=0.5, oz=0.1)
for wall in walls:
    plot_face(ax, wall, linecolour="#333", linewidth=3, linealpha=0.5,
              facefill=True, facealpha=0.3, vertexfill=False)

    
#%% =========================== Add EMF vector =============================

# Rotation matrix from geographic frame to cage frame
R_G2C = EulerRotation().rz(-23, deg=True)
# Rotation matrix from sensor frame to cage frame
R_S2C = EulerRotation().rx(-180, deg=True)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
local_emf = np.dot(R_G2C, local_emf())

print("C-frame EMF vector:", local_emf.round(3))

plot_arrow(ax, [0,0,0], local_emf, alpha=1, scaling=0.02, alr=0.1, color="orange")


#%% =========================== Add Measured =============================

filename = "WeekdayTest_2023-07-03_10.41.03_SPARSE50.dat"

data_r = read_data(filename, header=True)

print("Dataset contains",len(data_r[0]), "data points.")

# mG to uT conversion:
# for i in range(1, len(data_r)):
#     for j in range(len(data_r[0])):
#         data_r[i][j] *= 0.1

# time_plot(data_r, ylim=[-50, 50], 
#           title="Rack Elimination Test - Measurements in cage frame")

# Alternative data format with field data as vectors
# These vectors are converted from the sensor frame (S) straight to the 
# cage frame (C).


# def data_savgol_filter(data, windowlength, polyorder=4):
#     # Smooths three data channels using a Savitzky-Golay filter
#     # Data must be in type A
#     x_filtered = savgol_filter(data[1], windowlength, polyorder)
#     y_filtered = savgol_filter(data[2], windowlength, polyorder)
#     z_filtered = savgol_filter(data[3], windowlength, polyorder)
#     return [data[0], x_filtered, y_filtered, z_filtered]

data_r = data_savgol_filter(data_r, 51, 4)


# time_plot(data_r, ylim=[-50, 50],
#           title="Rack Elimination Test - Measurements in cage frame")

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

mean_data, sd_data = vector_normal_distribution(data_vectors)
print("C-frame data: mean, sd:", mean_data.round(3), sd_data.round(3))

# Plot all points (resource intensive!)
# for data_vector in data_vectors:
#     plot_arrow(ax, [0,0,0], data_vector, alpha=1, scaling=0.02, alr=0.1, color="purple")

# Plot only average measured vector:
plot_arrow(ax, [0,0,0], mean_data, alpha=1, scaling=0.02, alr=0.1, 
           color="purple")

# time_plot(zipBA(data), ylim=[-50, 50])


# Normalize the data (in C frame) with respect to the EMF (in C frame):
data_n = deepcopy(data)
    
for i in range(len(data_n)):
    data_n[i][1] = data_n[i][1] - local_emf

#
# time_plot(zipBA(data_n), ylim=[-50, 50],
#           title="Weekday Test - Normalized disturbance vector")

data_n_vectors = []
for d in data_n:
    data_n_vectors.append(d[1])
    
mean_data_n, sd_data_n = vector_normal_distribution(data_n_vectors)
print("Normalized data: mean, sd:", mean_data_n.round(3), sd_data_n.round(3))
print("Absolute strength of disturbance:", 
      round(np.linalg.norm(mean_data_n),3),
      "uT (", 
      round(100*(np.linalg.norm(mean_data_n)/np.linalg.norm(local_emf)),2),
      "% of EMF)")

plot_arrow(ax, [0,0,0], mean_data_n, linewidth=3, alpha=1, scaling=0.02, 
           alr=0.3, color="black")

plt.show()