# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 11:53:51 2023

@author: Main
"""

"""
from GaussmeterLib import ReadData, AnalyzeRate
from matplotlib.gridspec import GridSpec


Coordinates of AE building:
    51.99002134132983, 4.375289736921873
    altitude: 30 m (estimated altitude of 8 floors)

 - Mean values of Bx, By, Bz
 - Variance/SD of Bx, By, Bz
 - Find local magnetic field vector
 - Subtract local EMF from measurements by datapoint
 - Describe mean magnitudes and SD of normalized data
 - Describe rate of change in terms of angle and magnitude

"""
from GaussmeterAnalysis import time_plot
import matplotlib.pyplot as plt
import numpy as np
from array import array
from copy import deepcopy


#%% FUNCTIONS

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


#%% # Rotation frames

R_G2C = Rz(-(90+23), deg=True)
R_C2G = -R_G2C

R_S2C = np.dot(Rz(-90, deg=True), Rx(-180, deg=True))
R_C2S = -R_S2C


#%% Local EMF (calculated 03-07-2023)

# EMF vector at AE building (WMM-2020 model)
B_EMF_WMM2020  = np.array([0.7216, 19.1629, -45.4592])
# EMF vector at AE building (IGRF2020 model)
B_EMF_IGRF2020 = np.array([0.7313, 19.1870, -45.4557])

# Take average of the two as baseline B_EMF:
B_EMF_G = 0.5 * (B_EMF_WMM2020 + B_EMF_IGRF2020)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
B_EMF = np.dot(R_G2C, B_EMF_G)

#%% 
filename = "ChamberTest3_2023-06-30_15.04.31_SPARSE50.dat"

##%% DEBUG

data_r = read_data(filename, header=True)


# mG to uT conversion:
for i in range(1, len(data_r)):
    for j in range(len(data_r[0])):
        data_r[i][j] *= 0.1

# time_plot(data_r, ylim=[-50, 50])


# Alternative data format with field data as vectors
# These vectors are converted from the sensor frame (S) straight to the 
# cage frame (C).

data = []
for i in range(len(data_r[0])):
    data.append([
        data_r[0][i], 
        np.dot(R_S2C, np.array([data_r[1][i], 
                                data_r[2][i], 
                                data_r[3][i]]))
        ])

time_plot(zipBA(data), ylim=[-50, 50])

# Normalize the data (in C frame) with respect to the EMF (in C frame):
data_n = deepcopy(data)
    
for i in range(len(data_n)):
    data_n[i][1] = data_n[i][1] - B_EMF


time_plot(zipBA(data_n), ylim=[-50, 50])



# for i in range(len(data_n)):
#     # data_n_r.append([data_n[i][0], *list(data_n[i][1])])
#     data_n_r[0].append(data_n[i][0])
#     data_n_r.append([
#         data_n[i][0],
#         data_n[i][1][0],
#         data_n[i][1][1],
#         data_n[i][1][2]])

#%% Mean values and SDs of B_ambient


# meanX = sum(data[1])/len(data[1])
# meanY = sum(data[2])/len(data[2])
# meanZ = sum(data[3])/len(data[3])

# sigX = 0
# sigY = 0
# sigZ = 0

# for i in range(len(data[0])):
#     sigX += (data[1][i] - meanX)**2
#     sigY += (data[2][i] - meanY)**2
#     sigZ += (data[3][i] - meanZ)**2
    
# sigX = sigX/len(data[0])
# sigY = sigY/len(data[0])
# sigZ = sigZ/len(data[0])

# print(round(meanX,3), round(sigX**0.5,3))
# print(round(meanY,3), round(sigY**0.5,3))
# print(round(meanZ,3), round(sigZ**0.5,3))

# B_ambient = np.array([meanX, meanY, meanZ])
