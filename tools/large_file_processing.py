import os
from time import time
from datetime import datetime, timezone

from .GaussmeterLib import multipart_filenames

class DataProcessor:
    def __init__(self, filename, sample_rate, n_lines_header=10):
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

    def _seek_to_startline(self, f, verbose=1):
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

    def _seek_to_nextline(self, f, verbose=1):
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

        with open(self.filename, 'rb') as f:

            # Starting step in middle of file

            lowerbound = self.header_size+1
            upperbound = self.filesize

            iterations = 0
            while iterations <= max_iterations:  # Loop: After 50 tries, give up

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
        start_time = self.return_by_line_byte(self.get_header_size() + 1)
        if unix:
            return start_time
        else:
            return datetime.utcfromtimestamp(start_time)

    def get_end_time(self, unix=True):
        end_time = self.return_by_line_byte(self.filesize)
        if unix:
            return end_time
        else:
            return datetime.utcfromtimestamp(end_time)

    @staticmethod
    def unix2datetime(unix):
        return datetime.utcfromtimestamp(unix)

    @staticmethod
    def datetime2unix(dt):
        return dt.timestamp()

    @staticmethod
    def unix_hour(unix):
        dt = datetime.utcfromtimestamp(unix)
        return "{}:{}:{}".format(str(dt.hour).rjust(2, "0"),
                                 str(dt.minute).rjust(2, "0"),
                                 str(dt.second).rjust(2, "0"))

    def find_line_bytes_hourly(self, verbose=0):
        """
        This function will:
         - Look at the dataset, and find the Unix time of the first entry
         - It will also remember the line_byte of this entry
         - Then it will predict the Unix time corresponding to the next full hour
         - It will then fetch the line_byte at this time.
         - It will place the line_byte of the start and end of the interval into
            a list for storage.
         - It will do this for as long as
        """
        unix_start = self.get_start_time(unix=True)
        unix_end = self.get_end_time(unix=True)

        dt_start = self.unix2datetime(unix_start)
        unix_mid_start = self.datetime2unix(
            datetime(dt_start.year, dt_start.month, dt_start.day,
                     dt_start.hour + 1, 0, 0, tzinfo=timezone.utc)
        )

        duration_head = unix_mid_start-unix_start
        n_whole_hours = int((unix_end - unix_mid_start) / 3600)
        duration_mid = n_whole_hours * 3600
        duration_tail = (unix_end - unix_mid_start) - duration_mid

        duration_total = unix_end-unix_start

        if verbose >= 2:
            print("Found {} whole 1-hour segments, and a header and tail."
                  .format(n_whole_hours))

        if n_whole_hours == 0:
            if verbose >= 1:
                print("There is no need to split this file in hourly segments!")
            return -1
        else:
            # Pre-allocate array
            line_bytes_hourly = [[0, 0] for i in range(n_whole_hours + 2)]

            line_bytes_hourly[0][0] = self.return_line_byte(
                self.get_start_time(unix=True)
            )

            line_bytes_hourly[-1][1] = self.return_line_byte(
                self.get_end_time(unix=True)
            )

            for hour in range(n_whole_hours+1):
                line_bytes_hourly[hour+1][0] = \
                    self.return_line_byte(unix_mid_start + 3600 * hour)
                line_bytes_hourly[hour][1] = line_bytes_hourly[hour+1][0] - 1

            if verbose >= 1:
                for i, pair in enumerate(line_bytes_hourly):
                    if i == 0:
                        print("Header size :",
                              round((pair[1] - pair[0])/1024**2), "MB",
                              "       (start:",
                              self.unix_hour(unix_start),
                              ")")
                    elif i == len(line_bytes_hourly)-1:
                        print("Tail size :",
                              round((pair[1] - pair[0])/1024**2), "MB",
                              "       (end:",
                              self.unix_hour(unix_end),
                              ")")
                    else:
                        print("Interval", i+1, "size:",
                              round((pair[1] - pair[0])/1024**2), "MB",
                              "       (",
                              self.unix_hour(unix_mid_start+3600*(i-1)),
                              ")")

            return line_bytes_hourly

    # def segment_split(self, segments: list[list[int]], verbose=0):
    def segment_split(self, verbose=0):

        # First find the line_byte boundaries of the intervals
        intervals = self.find_line_bytes_hourly(verbose=verbose)

        # Loop over every interval
        for n in range(len(intervals)):

            # Generate output file name
            filename_n = multipart_filenames(self.filename, len(intervals))[n]

            # Open data file
            with open(self.filename, 'rb') as f:

                # In the data to the byte that belongs to the start of the interval
                f.seek(intervals[n][0], 1)
                len_interval = intervals[n][1] - intervals[n][0]

                with open(filename_n, 'w') as o:
                    o.write(
                        str(f.read(len_interval).replace(b'\r', b''), "utf-8")
                    )
            if verbose >= 2:
                print("Saved {} MB to file: {}"
                      .format(len_interval/1024**2, filename_n))

        return len(intervals)

# For inspiration, check:
# https://stackoverflow.com/questions/12859488/using-python-to-search-extremely-large-text-file
