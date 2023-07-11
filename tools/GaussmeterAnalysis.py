# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.fft import fft, fftfreq
from datetime import datetime

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

def time_plot(data, ylim=[-100, 100], 
              title="Title",
              labelx = "X", labely = "Y", labelz = "Z",
              ):
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(1,1)
    
    data_normalized = []
    for i in range(len(data[0])):
        data_normalized.append(data[0][i]-data[0][0])
    
    ax2 = fig.add_subplot(gs[0,0])
    ax2.set_ylim(ylim[0], ylim[1])
    # ax2.plot(data[0], data[1], color='red', label = labelx)
    # ax2.plot(data[0], data[2], color='green', label = labely)
    # ax2.plot(data[0], data[3], color='blue', label = labelz)  
    
    ax2.plot(data_normalized, data[1], color='red', label = labelx)
    ax2.plot(data_normalized, data[2], color='green', label = labely)
    ax2.plot(data_normalized, data[3], color='blue', label = labelz)  
    
    plt.legend()
    plt.title(title)
    plt.xlabel("Time [s]")
    plt.ylabel("B [uT]")
    
    start = datetime.utcfromtimestamp(data[0][0]).strftime('%Y-%m-%d %H:%M:%S')
    end = datetime.utcfromtimestamp(data[0][-1]).strftime('%Y-%m-%d %H:%M:%S')
    o_x = round(data_normalized[-1]/8)
    o_y = round((ylim[1]-ylim[0])/12)
    o_x = 0
    # o_y = 0
    plt.annotate(start, 
                 (data_normalized[0], ylim[0]), 
                 (data_normalized[0]-o_x, ylim[0]-o_y))
    plt.annotate(end, 
                 (data_normalized[-1], ylim[0]),
                 (data_normalized[-1]-o_x, ylim[0]-o_y))


def time_plot2(data, ylim=[-100, 100]):
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

    ax2 = fig.add_subplot(gs[1:,0])
    ax2.set_ylim(ylim[0], ylim[1])
    ax2.plot(data[0], data[1], color='red')
    ax2.plot(data[0], data[2], color='green')
    ax2.plot(data[0], data[3], color='blue')

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
    
    # print(len(list(range(len(sampling_intervals)))))
    # print(type(list(range(len(sampling_intervals)))))
    
    # print(len(list(sampling_intervals)))
    # print(type(list(sampling_intervals)))
    
    ax2.plot(list(range(len(sampling_intervals))), sampling_intervals, color='black')
    
    
def FourierPlot(data):
    pass

# filename = "output_StressTest1_1200_100_2023-06-19_15.19.02.dat"
# filename2 = "output_StressTest2_50_100_2023-06-19_15.58.24.dat"

# data = ReadData(filename, header=True)
# TimePlot(data)



# sr = 100
# SamplingPlot(data, sr)