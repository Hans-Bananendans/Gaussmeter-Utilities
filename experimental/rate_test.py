# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 11:53:51 2023

@author: Main
"""

"""
from GaussmeterLib import ReadData, AnalyzeRate
from matplotlib.gridspec import GridSpec
"""
from GaussmeterAnalysis import sampling_plot, time_plot
import matplotlib.pyplot as plt
import numpy as np
from array import array

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

def analyse_rate(data, n_samples, sr, verbose=0):
    """Analyses the time between each sample from the UNIX datapoints, and 
        performs basic statistical analysis on them."""
        
    timings = array("d", [0]*(n_samples-1))
    
    for i in range(len(data[0])-1):
        timings[i] = data[0][i+1] - data[0][i]
    
    mean = sum(timings)/len(timings)
    
    sd = 0
    for i in range(len(timings)):
        sd += (mean-timings[i])**2
    sd = sd/len(timings)
    
    return timings, mean, sd

sr = 2

filename = "ChamberTest3_2023-06-30_15.07.51_SPARSE50.dat"

data = read_data(filename, header=True)

# mG to uT conversion:
for i in range(len(data)):
    for j in range(len(data[0])):
        data[i][j] *= 0.1

sampling_plot(data, sr)
time_plot(data, ylim=[-50, 50])


meanX = sum(data[1])/len(data[1])
meanY = sum(data[2])/len(data[2])
meanZ = sum(data[3])/len(data[3])

sigX = 0
sigY = 0
sigZ = 0

for i in range(len(data[0])):
    sigX += (data[1][i] - meanX)**2
    sigY += (data[2][i] - meanY)**2
    sigZ += (data[3][i] - meanZ)**2
    
sigX = sigX/len(data[0])
sigY = sigY/len(data[0])
sigZ = sigZ/len(data[0])

print(round(meanX,3), round(sigX**0.5,3))
print(round(meanY,3), round(sigY**0.5,3))
print(round(meanZ,3), round(sigZ**0.5,3))

B_ambient = np.array([meanX, meanY, meanZ])
