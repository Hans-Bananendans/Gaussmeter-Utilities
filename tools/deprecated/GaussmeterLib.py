# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 13:29:48 2023

@author: Johan Monster

Library for Gaussmeter utilities.
"""
NAME = "Gaussmeter Recorder"
VERSION = "1.0.0"
VERSION_DATE = "28-06-2023"


#%% IMPORT DEPENDENCIES

import os
import sys
import ctypes
import threading
import PyDAQmx
import numpy as np
from tqdm import tqdm
from array import array
from time import time, sleep
from datetime import datetime
from argparse import ArgumentParser
from colorama import init, Fore, Style; init(autoreset=True)

import psutil

#%% DataRecorder Class

class DataRecorder():
    # Important distinction:
    # BUFFER_SIZE refers to the number of samples for each channel per buffer.
    # LEN_BUFFER refers to the actual length of the buffer, which is BUFFER_SIZE
    #   multiplied by the number of channels.
    # Example: For a buffer with a size of 100 samples, reading from 3 channels:
    #   BUFFER_SIZE = 100
    #   LEN_BUFFER = 300
    
    # Default parameters
    __n_channels: int       = 1
    __n_samples: int        = 3600
    __sr: int               = 60        # S/s
    __recording_time: float = 60.       # s
    __n_buffers: int        = 2
    __buffer_size: int      = 120       # S/channel
    __allow_n_samples_rounding: bool = True
    __v_min: float          = -1.25     # V
    __v_max: float          = 1.25      # V
    __rounding_t: int       = 9
    __rounding_v: int       = 3
    __scale_v: float        = 1.
    __name: str             = ""
    __filetype: str()       = ".dat"
    __use_header: bool      = True,         
    __header: str           = "!H",
    __filename: str         = "default.dat"
    __verbose: int          = 0
    __dt: float             = 0.        # s
    
    
    def __init__(self,
                 recording_time: float = 60.,    # sec
                 sr: int = 60,                  # S/s
                 name: str = "",                # Name of the recording session
                 filetype: str = ".dat",        # File extension of the data file
                 allow_n_samples_rounding: bool = True,
                 use_header: bool = True,         
                 header: str = "!H",              
                 verbose: int = 0               # Verbosity level
                 ):
        """
        Initialize DataRecorder object
        
        Parameters
        ----------
        recording_time : float, optional
            DURATION of recording in [s]. The default is 60.0.
        sr : int, optional
            SAMPLING RATE of the recording in [S/s]. The default is 60.
        name : str, optional
            NAME of the recording session, will be added to the filename of 
            the data file. The default is "".
        filetype : str, optional
            The FILE EXTENSION of the data file. The default is ".dat".
        allow_n_samples_rounding : bool, optional
            If allowed, this will increase the total recorded slightly such 
            that n_samples is an integer multiple of the buffer_size. 
            Turning this off may result in a number of zeroes at the end of 
            the data file. The default is True.
        use_header : bool, optional
            Whether to make the first line of the data file a header.
            The default is True.
        header : str, optional
            The CONTENTS of the header. The default is "!H".
        verbose : bool, optional
            The VERBOSITY LEVEL of the program. The default is 0.
        """
        
        self.__recording_time = recording_time
        self.__sr = sr
        self.__name = name
        self.__filetype = filetype
        self.__allow_n_samples_rounding = allow_n_samples_rounding
        self.__use_header = use_header
        self.__header = header
        self.__verbose=verbose     
        
        self.__dt = 1/self.__sr
        # Recompute n_samples based on given recording_time and sr
        self.generate_n_samples()
        
        # Generate filename
        self.generate_filename()
        
    
    #%% SETTINGS METHODS
    
    def set_sampling_settings(self, 
                              n_channels: int = __n_channels, 
                              n_samples: int = __n_samples,
                              sr: int = __sr,
                              n_buffers: int = __n_buffers,
                              buffer_size: int = __buffer_size,
                              v_min: int = __v_min,
                              v_max: int = __v_max,
                              allow_n_samples_rounding: bool = 
                                  __allow_n_samples_rounding
                              ):
        
        """
        Specifies additional sampling parameters of the recording object.

        Parameters
        ----------
        n_channels : int, optional
            NUMBER OF VOLTAGE CHANNELS to record. The default is 1.
        n_samples : int, optional
            TOTAL NUMBER OF SAMPLES to gather during the recording session.
        sr : int, optional
            SAMPLING RATE of the recording in [S/s].
        n_buffers : int, optional
            The NUMBER OF BUFFERS to use. Use at least 2. The default is 2.
        buffer_size : int, optional
            Number of datapoints in the buffer per channel.
        v_min : int, optional
            The MINIMUM VOLTAGE of the recording range. The default is -1.25V.
        v_max : int, optional
            The MAXIMUM VOLTAGE of the recording range. The default is +1.25V.
        allow_n_samples_rounding : bool, optional
            If allowed, this will increase the total recorded slightly such 
            that n_samples is an integer multiple of the buffer_size. 
            Turning this off may result in a number of zeroes at the end of 
            the data file.
        """
        
        self.__n_channels = n_channels
        self.__n_samples = n_samples
        self.__sr = sr
        self.__n_buffers = n_buffers
        self.__buffer_size = buffer_size
        self.__v_min = v_min
        self.__v_max = v_max
        
        self.__dt = 1/self.__sr
    
    def set_recording_time(self, recording_time: float):
        """
        Set the recording time in [s].

        Parameters
        ----------
        recording_time : float
            The desired recording time in seconds.
        """
        
        self.__recording_time = recording_time
        self.generate_n_samples()
        
    
    def set_data_settings(self,
                          name: str = __name,
                          filetype: str = __filetype,
                          use_header: bool = __use_header,         
                          header: str = __header,
                          rounding_t: int = __rounding_t,
                          rounding_v: int = __rounding_v,
                          scale_v: float = __scale_v
                          ):
        """
        Specifies data file parameters of the recording object.

        Parameters
        ----------
        name : str, optional
            NAME of the recording session, will be added to the filename of 
            the data file. The default is "".
        filetype : str, optional
            The FILE EXTENSION of the data file. The default is ".dat".
        use_header : bool, optional
            Whether to make the first line of the data file a header.
            The default is True.
        header : str, optional
            The CONTENTS of the header. The default is "!H".
        rounding_t : int, optional
            The number of rounding decimals of the UNIX TIME measurements.
            The default is 9.
        rounding_v : int, optional
            The number of rounding decimals of the VOLTAGE measurements. 
            The default is 3.
        scale_v : float, optional
            SCALING FACTOR by which all voltage measurements will be 
            multiplied. The default is 1.
        """
        
        self.__name = name
        self.__filetype = filetype
        self.__use_header = use_header       
        self.__header = header
        self.__rounding_t = rounding_t
        self.__rounding_v = rounding_v
        self.__scale_v = scale_v
        
        self.__filename = self.generate_filename()


    #%% GENERATION METHODS
    
    def generate_n_samples(self):
        self.__n_samples = self.__recording_time*self.__sr
        
        
    def generate_filename(self):     
        # Remove space characters from name to avoid complications
        name = self.__name.replace(" ", "_")
        
        timestamp = str(datetime.utcfromtimestamp(time()).strftime('%Y-%m-%d_%H.%M.%S'))
        filetype = "." + self.__filetype.lstrip(".")
        
        if name == "output" or name == "":
            filename = "output"+"_"
        else:
            filename = name+"_"
            
        # filename += str(self.__recording_time)+"_"+str(self.__sr)+"_"
        filename += timestamp+filetype 
        
        return filename
  
    def generate_header(self, timestamp: float):
        """ 
        Writes a static header to a specified output file. 
        Will create the file, and throw an error if it already exists.
        """
        with open(self.__filename, 'x') as output_file:
            # Write header
            # header = f"{'!H time UNIX [s]'.rjust(14, ' ')} {'Bx [uT]'.rjust(10, ' ')} {'Bx [uT]'.rjust(10, ' ')} {'Bx [uT]'.rjust(10, ' ')} \n"
            header  = "!H________RECORDING PROPERTIES__________\n"
            header += "Date: "+ \
                str(datetime.utcfromtimestamp(timestamp).strftime('%Y-%m-%d'))
            header += "   Start time: "+ \
                str(datetime.utcfromtimestamp(timestamp).strftime('%H:%M:%S'))
            header += "\n\n"
            header += "Duration:          "+str(self.__recording_time)
            header += "s ("+str(round(self.__recording_time/3600,2))+" h)\n"
            header += "Total samples:     "+str(self.__n_samples)+"\n"
            header += "Sampling rate:     "+str(self.__sr)+" S/s \n"
            header += "Channels:          "+str(self.__n_channels)
            header += " ["+str(self.__v_min)+", "+str(self.__v_max)+"] V\n"
            header += "Buffer size:       "+str(self.__buffer_size)
            header += " x "+str(self.__n_buffers)+"\n"
            header += "________________________________________\n"
                
                
            output_file.write(header)
        output_file.close()

    # def generate_header_old(self):
    #     """ 
    #     Writes a static header to a specified output file. 
    #     Will create the file if it does not exist.
    #     """
    #     # TODO: Expand function to display recording properties card
    #     with open(self.__filename, 'x') as output_file:
    #         # Write header
    #         header = f"{'!H time UNIX [s]'.rjust(14, ' ')} {'Bx [uT]'.rjust(10, ' ')} {'Bx [uT]'.rjust(10, ' ')} {'Bx [uT]'.rjust(10, ' ')} \n"
    #         output_file.write(header)
    #     output_file.close()

    def progress_update(self, t0, n_loops, i_loop):
        t_elapsed = time()-t0
        t_eta = t_elapsed*(n_loops)/(i_loop+1)
        
        print("Progress:", 
              str(Fore.GREEN+str(i_loop+1).rjust(3," ") +" / "+str(n_loops)),
              "  |  Time elapsed:", Fore.CYAN+str(round(t_elapsed)),
              "s - remaining:", Fore.CYAN+str(round(t_eta-t_elapsed)),"s")

    
    def generate_simulation_summary(self):
        if self.__verbose >= 4:
            print("generate_simulation_summary() called.")
            
        # now = time()
        # n_samples = duration*sr
        # name = GenerateFilename(name, duration, sr, folder=folder, filetype=filetype, verbose=verbose)
       
        # estimated_filesize = EstimateFilesize(n_samples)
        # estimated_memory = EstimateMemory(n_samples, sr)
        
        # print(Fore.CYAN + " ==== Simulation summary ==== ")
        # print("Selected duration:       ", duration, "s") 
        # print("Sample rate:             ", sr, "S/s")     
        # print("Total samples:           ", n_samples, "samples")
        # print("Estimated memory usage:  ", estimated_memory, "MB")  # TODO add memory estimation
        # print("Estimated file size:     ", estimated_filesize, "MB")
        # print(" ")
        # print("Current time:            ", datetime.utcfromtimestamp(now).strftime('%Y-%m-%d_%H.%M.%S'))
        # print("Expected time at finish: ", datetime.utcfromtimestamp(now+duration).strftime('%Y-%m-%d_%H.%M.%S'))
        # print(" ")
        # print("Saving to:", str(Fore.CYAN + name))
        # print(Fore.CYAN + " ============================ ")   

    #%% RECORDING METHODS

    def setup_channels(self, taskHandle: ctypes.c_void_p):
        """
        Sets up a number (up to four) of analog voltage channels within the 
        context of a DAQmx Task object. 

        Parameters
        ----------
        taskHandle : ctypes.c_void_p
            Handle of the DAQmx Task object.

        Returns
        -------
        channels : list
            List with all analog voltage channel objects.
        """
            
        hw_channel_names = ("Dev1/ai0", "Dev1/ai1", "Dev1/ai2", "Dev1/ai3")
        
        channels = []
        
        for channel in range(self.__n_channels):
            # Create a channel to measure voltage, and add the channel to the 
            # task specified by the taskhandle.
            
            chan = PyDAQmx.DAQmxCreateAIVoltageChan(
                taskHandle,                 # TaskHandle object to link task
                hw_channel_names[channel],  # Name of the physical channel
                "",                         # Name assigned to virtual channel
                PyDAQmx.DAQmx_Val_Diff,     # Input terminal config -> look at https://documentation.help/NI-DAQmx-C-Functions/DAQmxCreateAIVoltageChan.html
                self.__v_min,               # minVal of voltage in <units>
                self.__v_max,               # maxVal of voltage in <units>
                PyDAQmx.DAQmx_Val_Volts,    # Meaning of <units>
                None                        # Name of custom scale to apply to channel
                )
            channels.append(chan)
        
        return channels
    
    def setup_buffers(self):
        buffers = []
        tstart_buffers = []
        
        buffer_length = self.__buffer_size*self.__n_channels
        
        for buffer in range(self.__n_buffers):
            buffers.append(np.zeros((buffer_length,), dtype=np.float64))
            tstart_buffers.append(0.0)
        
        if self.__verbose >= 3:
            print("Successfully set up", len(buffers),
                  "buffers of length", self.__buffer_size,".")
        
        return buffers, tstart_buffers
        

    def idle(self, start_event: threading.Event, predelay: float):
        """
        Idling thread implementing a predelay pause.

        Parameters
        ----------
        start_event : threading.Event
            Threading event that marks the moment the recording may start.
        predelay : float
            Predelay in seconds.

        Returns
        -------
        None.

        """
        sleep(predelay)
        start_event.set()
        if self.__verbose >= 3:
            print("Starting simulation after a",predelay,"s predelay.")
    
    
    def record(self, predelay: float = 0):
        
        t0 = time() # Mark start of recording function
        
        # Verbosity warning:
        self.verbosity_warning()
        
        # If buffer_size is unspecified (-1), then auto-set it to:
        #   - if n_samples < 10000: buffer_size = int(sr)
        #   - if n_samples >= 10000: buffer_size = 1% of n_samples
        if self.__buffer_size == -1:
            if self.__n_samples < 10_000:
                self.__buffer_size = int(self.__sr)
            elif self.__n_samples >= 10_000:
                self.__buffer_size = int(0.01*self.__n_samples)
                
        # Round up sample number to nearest specified buffer size, if this 
        # is allowed and required, respectively.
        if self.__allow_n_samples_rounding:
            dmq, dmm = divmod(self.__n_samples, self.__buffer_size)
            if dmm != 0:
                if self.__verbose >= 1:
                    print("Note: Number of samples was rounded up from",
                          self.__n_samples, "to", (dmq+1)*self.__n_samples)
                self.__n_samples = (dmq+1)*self.__buffer_size
        
        # Number of loops to perform based on n_samples and buffer_size
        n_loops = int(self.__n_samples/self.__buffer_size)
        
        # Pre-fill data file with header
        if self.__use_header:
            self.generate_header(t0)
        
        # Initializing global buffer objects:
        global buffers, tstart_buffers
        # buffers = []
        # tstart_buffers = []
        buffers, tstart_buffers = self.setup_buffers()
        len_buffers = len(buffers)
        
        # global buffer1, buffer2, buffers
        # buffer1 = np.zeros((BUFFER_SIZE*N_CHANNELS,), dtype=np.float64)
        # buffer2 = np.zeros((BUFFER_SIZE*N_CHANNELS,), dtype=np.float64)
        # buffers = [buffer1, buffer2]
        # len_buffers = len(buffers)
        # global tstart_buffer1, tstart_buffer2, tstart_buffers
        # tstart_buffer1 = 0.0
        # tstart_buffer2 = 0.0
        # tstart_buffers = [tstart_buffer1, tstart_buffer2]
        
        # Declaration of Task handle object
        taskHandle = PyDAQmx.TaskHandle()
        # TODO: Clarify this step
        read = PyDAQmx.DAQmxTypes.int32()
        
        try:
            # DAQmx Configure Code
            # Create the task object, give it a name, and pass the taskhandle object
            #   by reference (using byref() from ctypes)
            PyDAQmx.DAQmxCreateTask("Record",ctypes.byref(taskHandle))
        
            # Configure the channels and route them to the task object
            self.setup_channels(taskHandle)
            
            # Configure clock settings on datalogger
            PyDAQmx.DAQmxCfgSampClkTiming(
                taskHandle,                     # TaskHandle object to link task
                "",                             # The source terminal of the sample clock (if nothing, it uses internal device clock)
                self.__sr,                      # The sampling rate in samples per second per channel
                PyDAQmx.DAQmx_Val_Rising,       # Active edge, essentially whether samples are gathered at the rising or falling edge of the clock signal
                PyDAQmx.DAQmx_Val_ContSamps,    # Specify whether finite number of samples are gathered, whether continuous samples are gathered until the task is stopped, or whether to use hardware timed single point sample mode.
                self.__buffer_size              # Number of samples to acquire or generate for EACH channel in the task (only if FiniteSamps mode is on). If ContSamps is on, this value determines buffer size.
                )
            
            # Start Task (overhead ~1000 ns)
            PyDAQmx.DAQmxStartTask(taskHandle)  # Transitions the task from the committed state to the running state, which begins measurement or generation.
    
            if self.__verbose >= 3:
                print("Setup complete.")
                
            # If predelay is specified, idle in thread until it is time to start
            leftover_delay = time()-t0 - predelay
            if leftover_delay > 0:
                start_event = threading.Event()
                th_idle = threading.Thread(
                    target=self.idle, 
                    args=(start_event, leftover_delay)
                    )
                th_idle.daemon = True
                th_idle.start()
                th_idle.join()      
                
            if self.__verbose >= 3:
                print("Starting loop...")
            t1 = time()
            
            try:
                for i_loop in range(n_loops):
                    
                    th_daq  = threading.Thread(
                        target=self.daq, 
                        args=(taskHandle, i_loop%len_buffers, read)
                        )
                    th_write = threading.Thread(
                        target=self.write, 
                        args=(self.__filename,(i_loop-1)%len_buffers)
                        )
                    th_daq.daemon = True
                    th_write.daemon = True
                    
                    th_daq.start()
                    if i_loop != 0:
                        th_write.start()
                        th_write.join() # New DAQ thread may not start before writing is done.
                    th_daq.join() # Wait for DAQ to finish before starting new loop
            
                    if self.__verbose >= 2:
                        self.progress_update(t1, n_loops, i_loop)
            except (KeyboardInterrupt, SystemExit):
                print("Threads terminated by KeyboardInterrupt.")
                # This implementation may lead to import lock, see also:
                # https://stackoverflow.com/questions/46290045/import-silently-kills-thread/46354248#46354248
                # TODO: See if this requires addressing and if so, how exactly.
                
            if self.__verbose >= 4:
                print("Loop finished.")  
        except PyDAQmx.DAQError as err:
            print("DAQmx Error:",err)
        finally:
            if taskHandle:
                # Terminate DAQmx task object
                PyDAQmx.DAQmxStopTask(taskHandle)
                PyDAQmx.DAQmxClearTask(taskHandle)
                if self.__verbose >= 4:
                    print("Task terminated succesfully.")  
            if self.__verbose >= 4:
                print("Program done.")              
        # ========================

    def daq(self, taskHandle, i_buffer, read):
        """ 
        Threaded data acquisition function.
        """
        if self.__verbose >= 4:
            print("daq() called. Reading to buffer", i_buffer,"...")
            if self.__verbose >= 5:
                t_start = time()
        
        tstart_buffers[i_buffer] = time()
        
        # DAQmx Read Code
        PyDAQmx.DAQmxReadAnalogF64(         # Reads multiple floating-point samples from a task that contains one or more analog input channels.
            taskHandle,                     # TaskHandle object to link task
            self.__buffer_size,             # The number of samples, per channel, to read. If readArray does not contain enough space, this function returns as many samples as fit in readArray.
            1.5*(self.__buffer_size/self.__sr), # The timeout in seconds. (here 1.5x of expected duration)
            PyDAQmx.DAQmx_Val_GroupByChannel, # fillMode: interleaving is OFF with DAQmx_Val_GroupByChannel and ON with DAQmx_Val_GroupByScanNumber
            buffers[i_buffer],              # readArray
            self.__buffer_size*self.__n_channels,# length of the readArray
            ctypes.byref(read),             # data type of readArray passed by reference?
            None                            # Reserved field?
            )
        
        if self.__verbose >= 4:
            print("daq() finished.")
            if self.__verbose >= 5:
                t_end = time()
                t_eta = t_start+self.__buffer_size/self.__sr
                print("Expected duration:",round(t_eta-t_start,6),"s")
                print("Actual duration:  ",round(t_end-t_start,6),"s")
                print("Time shift",
                      round(1_000_000*(t_end-t_eta),1),
                      "ns")
                print("Delay factor:     ",
                      round(100*((t_end-t_start)/(t_eta-t_start)-1),1),
                      "%")
    

    def write(self, filename, i_buffer, empty_buffer=True):
        """ 
        Threaded data writing function.
        """
        if self.__verbose >= 4:
            print("write() called. Writing from buffer", i_buffer,"...")
            if self.__verbose >= 5:
                t0 = time()
        
        with open(filename, 'a') as output_file:
            for sample in range(self.__buffer_size):
                # TODO: Make universal for n_channels:
                # Can construct "{} {} ... {}\n" string and then
                # array = round(buffers[i_buffer][0][buffer_line],rounding_v)
                # then use "{} {} ... {}\n".format(t_term, *array)
    
                line = "{} {} {} {}\n".format(
                    round(tstart_buffers[i_buffer]
                          +sample*self.__dt,self.__rounding_t),
                    round(self.__scale_v 
                          * buffers[i_buffer][sample],
                          self.__rounding_v),
                    round(self.__scale_v 
                          * buffers[i_buffer][sample+1*self.__buffer_size], 
                          self.__rounding_v),
                    round(self.__scale_v 
                          * buffers[i_buffer][sample+2*self.__buffer_size], 
                          self.__rounding_v),
                    )
    
                output_file.write(line)
            
        output_file.close()
        
        if self.__verbose >= 5:
            t1 = time()
        
        # Empty the buffer, for small overhead gain increased "debuggability".
        if empty_buffer:
            buffers[i_buffer] = np.zeros((self.__buffer_size*self.__n_channels,), 
                    dtype=np.float64)
            tstart_buffers[i_buffer] = 0.0
    
        if self.__verbose >= 5:
            t2 = time()
    
        if self.__verbose >= 4:
            print("write() finished.")
            if self.__verbose >= 5:
                print("Total time taken:          ",
                      round(1000*(time()-t0),3), "ms")
                print("Time spent writing:        ",
                      round(1000*(t1-t0),3), "ms")
                print("Time spent emptying buffer:",
                      round(1000_000*(t2-t1),3), "ns")
                
    def verbosity_warning(self):
        if self.__verbose > 1:
            print("Note: Higher levels of verbosity may add additional overhead.")
    
#%% OTHER FUNCTIONS    
    
    
# def TestProgressBar(n_samples, sr):
#     progress_bar_disable = False
#     progress_items = 100
#     progress_bar = tqdm(range(progress_items), ncols=80, 
#                     disable = progress_bar_disable,
#                     colour="CYAN",  position=0, leave=True,
#                     desc="TEST BAR", postfix=" ".ljust(16," "),
#                     bar_format="{l_bar}{bar}| {elapsed} - ETA:{remaining}{postfix}")
#     update_interval = n_samples/100
#     for i in range(n_samples):
#         if i%update_interval == 0:
#             progress_bar.update(1)
#         sleep(1/sr)
        


# def SetupProgressBar(n_samples, progress_bar_disable=False):
#     progress_bar_disable = False
#     progress_items = 100
#     progress_bar = tqdm(range(progress_items), ncols=80, 
#                     disable = progress_bar_disable,
#                     colour="CYAN",  position=0, leave=True,
#                     desc="Recording", postfix=" ".ljust(24," "),
#                     bar_format="{l_bar}{bar}| {elapsed} - ETA:{remaining}{postfix}")
#     update_interval = n_samples/100
#     return progress_bar, update_interval

# def SimulationSummary(duration, sr, name="", folder="", filetype=".dat", verbose=0):

#     now = time()
#     n_samples = duration*sr
#     name = GenerateFilename(name, duration, sr, folder=folder, filetype=filetype, verbose=verbose)
   
#     estimated_filesize = EstimateFilesize(n_samples)
#     estimated_memory = EstimateMemory(n_samples, sr)
    
#     print(Fore.CYAN + " ==== Simulation summary ==== ")
#     print("Selected duration:       ", duration, "s") 
#     print("Sample rate:             ", sr, "S/s")     
#     print("Total samples:           ", n_samples, "samples")
#     print("Estimated memory usage:  ", estimated_memory, "MB")  # TODO add memory estimation
#     print("Estimated file size:     ", estimated_filesize, "MB")
#     print(" ")
#     print("Current time:            ", datetime.utcfromtimestamp(now).strftime('%Y-%m-%d_%H.%M.%S'))
#     print("Expected time at finish: ", datetime.utcfromtimestamp(now+duration).strftime('%Y-%m-%d_%H.%M.%S'))
#     print(" ")
#     print("Saving to:", str(Fore.CYAN + name))
#     print(Fore.CYAN + " ============================ ")   


# def ProgressReport(S, n_samples, sr, t_start):
#     progress = round(S/n_samples*100,0)
#     t_actual = time()-t_start
#     t_expected = t_start+1/sr*S
#     t_behind = t_expected-t_actual
#     delay_factor = t_actual/t_expected
#     t_eta = delay_factor*(t_start+n_samples/sr)
    
#     print("Progress:", 
#           progress,
#           "%. | Running for:", 
#           datetime.utcfromtimestamp(t_actual).strftime('%H:%M:%S'),
#           "(",
#           round(t_behind,1),
#           "s behind schedule) | ETA:",
#           datetime.utcfromtimestamp(t_eta).strftime('%H:%M:%S')
#           )

def print_version(version, date, name, verbose=0):
    """Prints version number if verbosity >= 0, and additionally prints
        program name and version date if verbosity is >= 1."""
    if verbose == 0:
        print(version)    
    elif verbose >= 1:
        print(Fore.CYAN + name)
        print(version)
        print(Fore.CYAN + date)
    
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