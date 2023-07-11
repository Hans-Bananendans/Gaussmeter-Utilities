import PyDAQmx
import ctypes
import threading
import numpy as np
# from array import array
from time import time, sleep


# Important distinction:
# BUFFER_SIZE refers to the number of samples for each channel per buffer.
# LEN_BUFFER refers to the actual length of the buffer, which is BUFFER_SIZE
#   multiplied by the number of channels.
# Example: For a buffer with a size of 100 samples, reading from 3 channels:
#   BUFFER_SIZE = 100
#   LEN_BUFFER = 300

def setup_channels(taskHandle, n_channels, v_min, v_max):
    """ 
    Sets up up to four voltage channels on a NI-DAQmx compatible device
    Specify:
        taskHandle:             handle of the contextual task object
        n_channels <int>:       number of channels to set up.
        v_min <float>:          lower voltage bound in V
        v_max <float>:          upper voltage bound in V
        
    Estimated overhead: ~2000 ns per channel
    """
        
    hw_channel_names = ("Dev1/ai0", "Dev1/ai1", "Dev1/ai2", "Dev1/ai3")
    
    channels = []
    
    for channel in range(n_channels):
        # Create a channel to measure voltage, and add the channel to the task
        #   specified by the taskhandle.
        
        chan = PyDAQmx.DAQmxCreateAIVoltageChan(
            taskHandle,                 # TaskHandle object to link task
            hw_channel_names[channel],  # Name of the physical channel
            "",                         # Name assigned to virtual channel
            PyDAQmx.DAQmx_Val_Diff,     # Input terminal config -> look at https://documentation.help/NI-DAQmx-C-Functions/DAQmxCreateAIVoltageChan.html
            v_min,                      # minVal of voltage in <units>
            v_max,                       # maxVal of voltage in <units>
            PyDAQmx.DAQmx_Val_Volts,    # Meaning of units
            None                        # Name of custom scale to apply to channel
            )
        channels.append(chan)
    
    return channels

def verbosity_warning(verbose):
    if verbose > 1:
        print("Note: Higher levels of verbosity may add additional overhead.")

def write_header(filename):
    """ 
    Writes a static header to a specified output file. 
    Will create the file if it does not exist.
    """
    with open(filename, 'a') as output_file:
        # Write header
        header = f"{'!H time UNIX [s]'.rjust(14, ' ')} {'Bx [uT]'.rjust(10, ' ')} {'Bx [uT]'.rjust(10, ' ')} {'Bx [uT]'.rjust(10, ' ')} \n"
        output_file.write(header)
    output_file.close()

def progress_update(t0, n_loops, i_loop):
    t_elapsed = time()-t0
    t_eta = t_elapsed*(n_loops)/(i_loop+1)
    
    print("Progress:",i_loop+1, "/", n_loops,
          "  |  Time elapsed:",round(t_elapsed),
          "s - remaining:", round(t_eta-t_elapsed),"s")
    

def daq(taskHandle, n_channels, sr, buffer_size, i_buffer, read, verbose):
    """ 
    Threaded data acquisition function.
    """
    if verbose >= 3:
        print("daq() called. Reading to buffer", i_buffer,"...")
        if verbose >= 4:
            t_start = time()
    
    tstart_buffers[i_buffer] = time()
    
    # DAQmx Read Code
    PyDAQmx.DAQmxReadAnalogF64(         # Reads multiple floating-point samples from a task that contains one or more analog input channels.
        taskHandle,                     # TaskHandle object to link task
        buffer_size,                      # The number of samples, per channel, to read. If readArray does not contain enough space, this function returns as many samples as fit in readArray.
        1.5*(buffer_size/sr),             # The timeout in seconds. (here 1.5x of expected duration)
        PyDAQmx.DAQmx_Val_GroupByChannel, # fillMode: interleaving is OFF with DAQmx_Val_GroupByChannel and ON with DAQmx_Val_GroupByScanNumber
        buffers[i_buffer],              # readArray
        buffer_size*n_channels,         # length of the readArray
        ctypes.byref(read),             # data type of readArray passed by reference?
        None                            # Reserved field?
        )
    
    if verbose >= 3:
        print("daq() finished.")
        if verbose >= 4:
            t_end = time()
            t_eta = t_start+buffer_size/sr
            print("Expected duration:",round(t_eta-t_start,6),"s")
            print("Actual duration:  ",round(t_end-t_start,6),"s")
            print("Time shift",
                  round(1_000_000*(t_end-t_eta),1),
                  "ns")
            print("Delay factor:     ",
                  round(100*((t_end-t_start)/(t_eta-t_start)-1),1),
                  "%")
            
    
def write(filename, n_channels, sr, buffer_size, i_buffer, rounding_t, rounding_v, scale=1, empty_buffer=True, verbose=0):
    """ 
    Threaded data writing function.
    """
    if verbose >= 3:
        print("write() called. Writing to buffer", i_buffer,"...")
        if verbose >= 4:
            t0 = time()
    
    dt = 1/sr
    with open(filename, 'a') as output_file:
        for sample in range(buffer_size):
            # TODO: Make universal for n_channels:
            # Can construct "{} {} ... {}\n" string and then
            # array = round(buffers[i_buffer][0][buffer_line],rounding_v)
            # then use "{} {} ... {}\n".format(t_term, *array)

            line = "{} {} {} {}\n".format(
                round(tstart_buffers[i_buffer]+sample*dt,rounding_t),
                round(scale * buffers[i_buffer][sample], rounding_v),
                round(scale * buffers[i_buffer][sample+1*buffer_size], rounding_v),
                round(scale * buffers[i_buffer][sample+2*buffer_size], rounding_v),
                )

            output_file.write(line)
        
    output_file.close()
    
    if verbose >= 4:
        t1 = time()
    
    # Empty the buffer, for small overhead gain increased "debuggability".
    if empty_buffer:
        buffers[i_buffer] = np.zeros((buffer_size*n_channels,), dtype=np.float64)
        tstart_buffers[i_buffer] = 0.0

    if verbose >= 4:
        t2 = time()

    if verbose >= 3:
        print("write() finished.")
        if verbose >= 4:
            print("Total time taken:          ",
                  round(1000*(time()-t0),3), "ms")
            print("Time spent writing:        ",
                  round(1000*(t1-t0),3), "ms")
            print("Time spent emptying buffer:",
                  round(1000*(t2-t1),3), "ms")
    
def main():
    N_CHANNELS = 3
    N_SAMPLES = 12000
    SR = 300
    BUFFER_SIZE = 3000
    SCALE = 1000 # To scale from V to mV
    # readout_period=0.0
    V_MIN=-1.25
    V_MAX=1.25
    verbose=2
    rounding_t = 9
    rounding_v = 3
    # cf = 0.0009
    filetype=".dat"
    filename="test_data_"+str(round(time()))+filetype
    header = True
    
    t0 = time()
    
    # Round up sample number to nearest specified buffer size:
    dmq, dmm = divmod(N_SAMPLES, BUFFER_SIZE)
    if dmm != 0:
        if verbose >= 1:
            print("Note: Number of samples was rounded up from",N_SAMPLES,
                  "to",(dmq+1)*BUFFER_SIZE)
        N_SAMPLES = (dmq+1)*BUFFER_SIZE
    n_loops = int(N_SAMPLES/BUFFER_SIZE)
        
    if SR == 0:
        dt = 0
    else:
        dt = 1/SR
    
    if header:
        write_header(filename)
    
    # Initializing global buffer objects:
    global buffer1, buffer2, buffers
    buffer1 = np.zeros((BUFFER_SIZE*N_CHANNELS,), dtype=np.float64)
    buffer2 = np.zeros((BUFFER_SIZE*N_CHANNELS,), dtype=np.float64)
    buffers = [buffer1, buffer2]
    len_buffers = len(buffers)
    
    global tstart_buffer1, tstart_buffer2, tstart_buffers
    tstart_buffer1 = 0.0
    tstart_buffer2 = 0.0
    tstart_buffers = [tstart_buffer1, tstart_buffer2]
    
    
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
        setup_channels(taskHandle, N_CHANNELS, V_MIN, V_MAX)
        
        # Configure clock settings on datalogger
        PyDAQmx.DAQmxCfgSampClkTiming(
            taskHandle,                     # TaskHandle object to link task
            "",                             # The source terminal of the sample clock (if nothing, it uses internal device clock)
            SR,                             # The sampling rate in samples per second per channel
            PyDAQmx.DAQmx_Val_Rising,       # Active edge, essentially whether samples are gathered at the rising or falling edge of the clock signal
            PyDAQmx.DAQmx_Val_ContSamps,    # Specify whether finite number of samples are gathered, whether continuous samples are gathered until the task is stopped, or whether to use hardware timed single point sample mode.
            BUFFER_SIZE                     # Number of samples to acquire or generate for EACH channel in the task (only if FiniteSamps mode is on). If ContSamps is on, this value determines buffer size.
            )
        
        # Start Task (overhead ~1000 ns)
        PyDAQmx.DAQmxStartTask(taskHandle)  # Transitions the task from the committed state to the running state, which begins measurement or generation.

        if verbose >= 4:
            print("Setup complete. Starting loop...")
        try:
            for i_loop in range(n_loops):
                
                # daq(taskHandle, n_channels, sr, buffer_size, i_buffer, read, verbose)
                th_daq  = threading.Thread(
                    target=daq, 
                    args=(
                        taskHandle,
                        N_CHANNELS,
                        SR,
                        BUFFER_SIZE,
                        i_loop%len_buffers,
                        read,
                        verbose
                        )
                    )
                # write(filename, n_channels, sr, buffer_size, i_buffer, rounding_t, rounding_v, empty_buffer=True, verbose=0)
                th_write = threading.Thread(
                    target=write, 
                    args=(
                        filename,
                        N_CHANNELS,
                        SR,
                        BUFFER_SIZE,
                        (i_loop+1)%len_buffers,
                        rounding_t,
                        rounding_v,
                        SCALE,
                        True,
                        verbose
                        )
                    )
                th_daq.daemon = True
                th_write.daemon = True
                
                th_daq.start()
                if i_loop != 0:
                    th_write.start()
                    th_write.join() # New DAQ thread may not start before writing is done.
                th_daq.join() # Wait for DAQ to finish before starting new loop
        
                if verbose >= 2:
                    progress_update(t0, n_loops, i_loop)
        except (KeyboardInterrupt, SystemExit):
            print("Threads terminated by KeyboardInterrupt.")
            # This implementation may lead to import lock, see also:
            # https://stackoverflow.com/questions/46290045/import-silently-kills-thread/46354248#46354248
            # TODO: See if this requires addressing and if so, how exactly.
            
        if verbose >= 4:
            print("Loop finished.")  
    except PyDAQmx.DAQError as err:
        print("DAQmx Error:",err)
    finally:
        if taskHandle:
            # DAQmx Stop Code
            t7 = time() # OH: 0 ns
            PyDAQmx.DAQmxStopTask(taskHandle)
            t8 = time() # OH: 0 ns
            PyDAQmx.DAQmxClearTask(taskHandle)
            t9 = time() # OH: 0 ns
            if verbose >= 4:
                print("Task terminated succesfully.")  
        if verbose >= 4:
            print("Program done.")              
    # ========================
          

#%% MAIN PROGRAM
if __name__ == "__main__":
    main()
    
    