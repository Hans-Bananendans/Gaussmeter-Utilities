# Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from scipy.signal import savgol_filter, filtfilt, iirnotch
from scipy.fft import rfft, rfftfreq


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

def multipart_filenames(filename_base, n_files, breakpoint=4):
    """
    Generates filenames for multipart file splits. This function is here so
     that it can both be used for generating these names in other functions,
     and so it can be used by a user to generate lists of filenames with the
     exact same algorithm (to reduce the risk of bugs).
    """
    filenames = []
    for i in range(1, n_files+1):
        filenames.append(
            filename_base[:-breakpoint]
            + "_s{}of{}".format(i, n_files)
            + filename_base[-breakpoint:]
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


def data_notch_filter(data, fn, fs, Q=None, nw=1):
    # If no notch quality factor Q is given, determine an appropriate one
    # using equation 15 from Nguyen2018 (DOI:
    # It will choose a Q with significant attenuation only in a range of
    # 0.5 Hz around the notch. This can be changed by setting "nw"
    if Q is None:
        Q = np.tan(np.pi*fn/fs) / (np.tan(np.pi*(fn+nw/2)/fs) - np.tan(np.pi*(fn-nw/2)/fs))

    if Q <= 0:
        raise ValueError(f"Q = {Q} but Q must be larger than 0")

    # Create narrow-band, second-order notch filter at 50 Hz
    b, a = iirnotch(fn, Q, fs)

    # Apply filter to data signals, using Gustafsson method for better
    # performance at the edges
    for i in (1, 2, 3):
        data[i] = filtfilt(b, a, data[i], method="gust")

    return data

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
    x_filtered = savgol_filter(data[1], windowlength, polyorder)
    y_filtered = savgol_filter(data[2], windowlength, polyorder)
    z_filtered = savgol_filter(data[3], windowlength, polyorder)
    return [data[0], x_filtered, y_filtered, z_filtered]

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
