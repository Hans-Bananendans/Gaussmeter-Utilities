from datetime import datetime
from time import time
import psutil
import os
import sys
from itertools import islice
from bitstring import BitArray

time0 = time()

# First open the file and read the number of lines

# filename = ".\data\LightOff_2023-07-04_10.57.37.dat" #
# filename = ".\data\WeekdayTest_2023-07-03_10.41.03.dat" # 8.6M lines
filename = ".\data\Weekday2500Hz_LightsOff_2023-07-06_19.12.02.dat" # 216M lines
# filename = ".\data\\test_data.dat"  # 500 lines

len_header = 9


# class DataProcessor:
#     def __init__(self, filename, n_lines_header=9):
#         # self.f = open(filename)
#         self.filename = filename
#         self.n_lines_header = n_lines_header
#
#
#         # self.header = self.get_header()
#
#     def refresh(self):
#         """Must be called EVERY time iterator f is used."""
#         self.f = open(filename)
#
#     def print_header(self):
#         for line in islice(self.f, 0, self.n_lines_header):
#             print(line[:-1])
#
#     def get_header(self):
#         header = []
#         for line in islice(self.f, 0, self.n_lines_header):
#             header.append(line[:-1])
#         return header
#
#     def get_n_lines_header(self):
#         return self.n_lines_header
#
#     def set_n_lines_header(self, n_lines_header: int):
#         self.n_lines_header = n_lines_header
#
#     def get_line(self, i, skip_header=True):
#         with open(self.filename, 'r') as f:
#             for line in islice(f,
#                                self.n_lines_header+i-1,
#                                self.n_lines_header+i):
#                 print(line[:-1])
#
#     def get_line2(self, n, skipheader=True):
#         for i, line in enumerate(self.f):
#         # check if the line number is specified in the lines to read array
#             if i == n:
#                 # print the required line number
#                 print(line)
#             break
#
#     def get_line3(self, n, skipheader=True):
#         with open(self.filename, 'r') as f:
#             for i, line in enumerate(f):
#                 # check if the line number is specified in the lines to read array
#                 if i == n:
#                     # print the required line number
#                     print(line[:-1])
#                     break
#
#     def find_first_occurrence(self, value, sample_rate):
#         # get_line(self,self.)
#         pass
#
# def find_timebins(filename, binlen=3600, sr=100, len_header=9,
#                   round_timeunit=True, verbose=0):
#
#     print(islice)
#     """Opens a file
#         - skips the header
#         - determines number of lines in file
#         - read first line, and establish datetime of first datapoint
#         - determines the line indices of bins of a given binlen
#         - returns a list of line indices of the start and end of each bin
#         - if round_timeunit = True, will round to minutes, hours, or days, but
#             only if binlen is 60, 3600, or 86400, or this option will be
#             ignored.
#         - Console feedback on verbose >= 1: Stages, memory usage, time elapsed
#     """
#
#     # if False:
#     #     for i in
#
#     n_lines = 0
#
#     with open(filename, 'r') as f:
#
#         for line in islice(f, 0, 10):
#             print(line[:-1])
#         print(islice)
#         # f.seek(int(step))
#         # row = parseRow(f.readline())
#         # print(row[0])
#
#         # for i, line in enumerate(f):
#         #     if i < 20:
#         #         print(line[:-1])
#         #         n_lines += 1
#         #         f.seek(1,1)
#
#
#         #     n_lines += 1 # Number of lines in file
#
#         # int(row.split(None, 1)[0]), row
#         #
#         # f.seek(0,0)
#         # print(f.readline()[:-1])
#         # f.seek(8,1)
#         # print(f.readline()[:-1])
#         # f.seek(-8,1)
#         # print(f.readline()[:-1])
#
#             # if i == len_header - 2:  # Should be first line!
#             #     print(line[:-1])
#             #
#             # if i < 20:
#             #     print(line[:-1])
#             #
#             # if i == 20:
#             #     test = line
#             #     print(test)
#             #     print(type(test))
#             #     print(len(test))


# fp = DataProcessor(filename)
# # fp.print_header()
#
#
# time1 = time()
# fp.get_line(10000000)
# fp.get_line(10000000)
# print("Elapsed time:", round(1000*(time()-time1), 3), 'ms')
#
# time2 = time()
# fp.get_line3(10000000)
# fp.get_line3(10000000)
# # fp.find_first(1688468300, 100)
#
# # print("RAM:",
# #       round(psutil.Process().memory_info().rss / (1024 * 1024), 3),
# #       "MB")
# print("Elapsed time:", round(1000*(time()-time2), 3), 'ms')





class DataProcessor:
    def __init__(self, filename, sample_rate, n_lines_header=10):
        # self.f = open(filename)
        self.filename = filename
        self.sample_rate = sample_rate
        self.n_lines_header = n_lines_header
        self.maxlinelength = 150

        self.filesize = os.path.getsize(filename)

        self.header_size = self.get_header_size()

    def get_header_size(self):
        with open(self.filename, 'rb') as f:
            header_size = 0
            newline_counts = 0
            while newline_counts < self.n_lines_header:
                if f.read(1) == b"\n":
                    newline_counts += 1
                    header_size += 1
                else:
                    header_size += 1
        return header_size

    def print_header(self):
        with open(self.filename, 'rb') as f:
            print(f.read(self.header_size))

    def set_maxlinelength(self, maxlinelength):
        self.maxlinelength = maxlinelength

    def _seek_to_startline(self, f):
        imax = 0
        while True:
            f.seek(-2, 1)
            if f.read(1) == b"\n":
                break
            elif imax > self.maxlinelength:
                if verbose >= 1:
                    print("Couldn't find a newline after",
                          self.maxlinelength, "characters. Please use",
                          "DataProcessor.set_maxlinelength() if your",
                          "lines are long (>150 characters).")
                break
            else:
                imax += 1

    def _seek_to_nextline(self, f):
        imax = 0
        while True:
            f.seek(0, 1)
            if f.read(1) == b"\n":
                break
            elif imax > self.maxlinelength:
                if verbose >= 1:
                    print("Couldn't find a newline after",
                          self.maxlinelength, "characters. Please use",
                          "DataProcessor.set_maxlinelength() if your",
                          "lines are long (>150 characters).")
                break
            else:
                imax += 1

    def return_by_line_byte(self, line_byte: int):
        if line_byte <= self.header_size:
            return -1
        else:
            with open(self.filename, 'rb') as f:
                f.seek(line_byte)
                self._seek_to_startline(f)
                return self._return_number(f)

    @staticmethod
    def _return_number(f):
        row = f.readline()
        number_b = row.split(b" ", 1)[0]
        # print(number_b, type(number_b))
        number = float(str(number_b, "utf-8"))
        return number

    def return_line_byte(self, n: float, max_iterations=50, verbose=0):
        """
        Finds the byte number of the start of the line containing this
        number. It does this with a binary search, resulting in very low
        computation time and memory overhead.
        See also: https://en.wikipedia.org/wiki/Binary_search_algorithm
        """
        if verbose >= 1:
            time_start = time()
        found = False

        # def seek_to_startline(f):
        #     imax = 0
        #     while True:
        #         f.seek(-2, 1)
        #         if f.read(1) == b"\n":
        #             break
        #         elif imax > self.maxlinelength:
        #             if verbose >= 1:
        #                 print("Couldn't find a newline after",
        #                       self.maxlinelength, "characters. Please use",
        #                       "DataProcessor.set_maxlinelength() if your",
        #                       "lines are long (>150 characters).")
        #             break
        #         else:
        #             imax += 1
        #
        # def seek_to_nextline(f):
        #     imax = 0
        #     while True:
        #         f.seek(0, 1)
        #         if f.read(1) == b"\n":
        #             break
        #         elif imax > self.maxlinelength:
        #             print("stopped, found none")
        #             break
        #         else:
        #             imax += 1

        with open(self.filename, 'rb') as f:

            # Starting step in middle of file

            lowerbound = self.header_size+1
            upperbound = self.filesize

            iterations = 0
            while iterations <= max_iterations:  # Loop: After 50 tries, give up

                # if step <= self.header_size:
                #     if verbose >= 1:
                #         print("Number cannot be found in set (too small")
                #     break
                # if step > self.filesize:
                #     if verbose >= 1:
                #         print("Number cannot be found in set (too large)")
                #     break

                # New failure conditions
                if lowerbound > upperbound:
                    if verbose >= 1:
                        print("Number cannot be found in set.")
                    break

                if iterations == 50:
                    if verbose >= 1:
                        print("Number could not be found after", iterations,
                              "iterations")
                    break

                step = int(lowerbound + (upperbound - lowerbound) / 2.)
                if verbose >= 2:
                    print("Focus on byte", step)

                f.seek(step)                # Focus on this byte of the file

                self._seek_to_startline(f)        # Find the start of the line

                number = self._return_number(f)   # Get the number that this line starts with

                # Now make a decision:
                if (n - number) > (1/self.sample_rate*1.1):          # If you found a number significantly smaller than the desired one (note the 10% margin)
                    iterations += 1
                    lowerbound = step+1
                    # step = int(step + step / (2**iterations-1))     # Try a higher number
                    if verbose >= 2:
                        print("Found number", number, "-> Trying a LARGER number.")

                elif (n - number) < 0:       # If instead you found a number larger than the desired one
                    iterations += 1
                    upperbound = step-1
                    # step = int(step - step / (2**iterations-1))     # Look at the middle of the lower half of the division
                    if verbose >= 2:
                        print("Found number", number, "-> Trying a SMALLER number.")

                # if (n - number) > 1.0:          # If you found a number significantly smaller than the desired one
                #     iterations += 1
                #     step = int(step + step / (2**iterations-1))     # Try a higher number
                #     if verbose >= 2:
                #         print("Found number", number, "-> Trying a LARGER number.")
                #
                # elif (n - number) < 0:       # If instead you found a number larger than the desired one
                #     iterations += 1
                #     step = int(step - step / (2**iterations-1))     # Look at the middle of the lower half of the division
                #     if verbose >= 2:
                #         print("Found number", number, "-> Trying a SMALLER number.")

                else:
                    found = True
                    if verbose >= 1:
                        print("I found number", number,
                              "at position", f.tell(),
                              "after", iterations, "iterations in",
                              round(1000*(time()-time_start), 3), "ms")
                    # TODO:                         # Sneak upward towards the number
                    break
            if found:
                return f.tell()
            else:
                return -1

    def previous_line_byte(self, byte_number):
        with open(self.filename, 'rb') as f:
            f.seek(byte_number-1)
            self._seek_to_startline(f)
            if f.tell() <= self.header_size:
                return -1
            else:
                return f.tell()

    def next_line_byte(self, byte_number):
        with open(self.filename, 'rb') as f:
            f.seek(byte_number+1)
            self._seek_to_nextline(f)
            if f.tell() > self.filesize:
                return -1
            else:
                return f.tell()

    def get_start_time(self, unix=True):
        start_time = self.return_by_line_byte(fp.get_header_size() + 1)
        if unix:
            return start_time
        else:
            return datetime.utcfromtimestamp(start_time)


    def find_line_bytes_hourly(self):
        """
        This function will:
         - Look at the dataset, and find the Unix time of the first entry
         - It will also remember the line_byte of this entry
         - Then it will predict the Unix time corresponding to the next full hour
         - It will then fetch the line_byte at this time.
        :return:
        """
        unix0 = self.get_start_time(unix=True)





fp = DataProcessor(filename, 5)
# fp.print_header()
print(fp.return_line_byte(1688670722+10000, verbose=0))
print(fp.get_start_time())

print("Elapsed time:", round(1000*(time()-time0), 3), 'ms')
    # print("RAM:",
    #       round(psutil.Process().memory_info().rss / (1024 * 1024), 3),
    #       "MB")
    # print("Elapsed time:", round(1000*(time()-time0), 3), 'ms')
    # print("File has", n_lines, "lines.")


# Given a sample rate and a dataset containing Unix timestamps:
#  - Cut the datafile into
#       - Unix time interval bins of a given length OR
#       - bins of n_samples
#  - Round the bins to whole hours, if possible (optional)
#  - Read out the data of the chunk, remember the cutoff line
#  - Do some processing, save some data
#  - Go on to the next chunk



# For inspiration, check:
# https://stackoverflow.com/questions/12859488/using-python-to-search-extremely-large-text-file