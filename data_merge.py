# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 11:53:51 2023

@author: Main
"""

import numpy as np
from array import array

from tools.GaussmeterLib import read_data

# def read_data(filename, header=False):
#     with open(filename, 'r') as f:
#         lines = f.readlines()
#         if lines[0][0:2] == "!H" or header == True:
#             lines = lines[10:]
#
#         t, x, y, z = [0]*len(lines), [0]*len(lines), [0]*len(lines), [0]*len(lines)
#
#         for i in range(len(lines)):
#             [tt, xx, yy, zz] = lines[i].split(" ")
#             [t[i], x[i], y[i], z[i]] = \
#                 [np.double(tt), float(xx), float(yy), float(zz)]
#         return [t, x, y, z]

filenames = (
    "./data/test_data_00.dat",
    "./data/test_data_01.dat",
    "./data/test_data_02.dat"
)

def write_data(filename, data, mode="a", len_header=10, verbose=0):
    with open(filename, mode) as fo:
        for i in range(len(data[0])):
            fo.write("{} {} {} {}\n".format(data[0][i],
                                            data[1][i],
                                            data[2][i],
                                            data[3][i]))


def data_merge(filenames, filename_merge: str = None, len_header=10, verbose=0):

    # If unspecified, base filename of merged file on first file
    if filename_merge is None:
        filename_merge = filenames[0][:-4]+"_MERGED"+filenames[0][-4:]

    # Create empty merged file
    with open(filename_merge, 'x') as fo:
        if verbose >= 2:
            print("Created merge file {}".format(filename_merge))

    # Copy header from first file
    with open(filenames[0], 'r') as f, open(filename_merge, 'a') as fo:
        for i in range(len_header):
            fo.write(f.readline())

    # Sequentially append contents of all files in filenames list to the
    # merged file
    for filename in filenames:
        write_data(filename_merge,
                   read_data(filename,
                             header=True,
                             len_header=len_header-1),
                   verbose=verbose)
        if verbose >= 2:
            print("Appended contents of {} to merge file".format(filename))

data_merge(filenames, verbose=3)

    # with open(filename_merge, 'x') as fs:
    #     pass
    #
    # with open(filename_sparse, 'a') as fs:
    #     for i_line in range(len_header - 1):
    #         fs.write(header_list[i_line])
    #     fs.write("SPARSE" + str(downsampling_factor) + "\n")
    #     fs.write(header_list[-1])

# def generate_filenames(filename_base, n_files):
#     filenames = []
#     for i in range(1, n_files+1):
#         filenames.append(
#             filename_base[:-4]
#             + "_s{}of{}".format(i, n_files)
#             + filename_base[-4:]
#         )
#     return filenames


# filenames = generate_filenames(filename_base, 25)

# filenames = [
#     "WeekdayTest_2023-07-03_10.41.03.dat"
#     ]


# for filename in filenames:
#     data = read_data(filename, header=True)
#
#     #%%
#
#     downsampling_factor = 25
#     chunks, _ = divmod(len(data[0]), downsampling_factor)
#
#     data_sparse = [[0.]*chunks, [0.]*chunks, [0.]*chunks, [0.]*chunks]
#
#     #%%
#     i_chunk = 0
#
#     while i_chunk < chunks:
#         for v in range(len(data)):
#             data_sparse[v][i_chunk] = data[v][0+i_chunk*downsampling_factor]
#         i_chunk += 1
#
#
#     #%% Save data
#     filename_sparse = filename[:-4]+"_SPARSE"+str(downsampling_factor)+filename[-4:]
#
#
#     len_header = 9
#
#     with open(filename, 'r') as f:
#         header_list = f.readlines(1000)[0:9]
#
#         with open(filename_sparse, 'x') as fs:
#             pass
#
#         with open(filename_sparse, 'a') as fs:
#             for i_line in range(len_header-1):
#                 fs.write(header_list[i_line])
#             fs.write("SPARSE"+str(downsampling_factor)+"\n")
#             fs.write(header_list[-1])
#
#             for i in range(len(data_sparse[0])):
#                 fs.write("{} {} {} {}\n".format(data_sparse[0][i],
#                                                 data_sparse[1][i],
#                                                 data_sparse[2][i],
#                                                 data_sparse[3][i]))
