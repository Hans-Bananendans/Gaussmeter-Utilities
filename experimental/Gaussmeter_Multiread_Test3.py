import nidaqmx
import threading
import numpy as np
from array import array
from time import time, sleep



def setup_channels(task, n_channels, v_min, v_max):
    """ 
    Sets up up to four voltage channels on a NI-DAQmx compatible device
    Specify:
        task <nidaqmx.Task>:    contextual nidaqmx task object
        n_channels <int>:       number of channels to set up.
        v_min <float>:          lower voltage bound in V
        v_max <float>:          upper voltage bound in V
    """
        
    ai0 = task.ai_channels.add_ai_voltage_chan(
        "Dev1/ai0",
        terminal_config=nidaqmx.constants.TerminalConfiguration.DIFF,
        min_val=v_min,
        max_val=v_max,
        units=nidaqmx.constants.VoltageUnits.VOLTS,
        )
    ai1 = task.ai_channels.add_ai_voltage_chan(
        "Dev1/ai1",
        terminal_config=nidaqmx.constants.TerminalConfiguration.DIFF,
        min_val=v_min,
        max_val=v_max,
        units=nidaqmx.constants.VoltageUnits.VOLTS,
        )
    ai2 = task.ai_channels.add_ai_voltage_chan(
        "Dev1/ai2",
        terminal_config=nidaqmx.constants.TerminalConfiguration.DIFF,
        min_val=v_min,
        max_val=v_max,
        units=nidaqmx.constants.VoltageUnits.VOLTS,
        )
    ai3 = task.ai_channels.add_ai_voltage_chan(
        "Dev1/ai3",
        terminal_config=nidaqmx.constants.TerminalConfiguration.DIFF,
        min_val=v_min,
        max_val=v_max,
        units=nidaqmx.constants.VoltageUnits.VOLTS,
        )
    if n_channels == 1:
        return (ai0)
    elif n_channels == 2:
        return (ai0, ai1)
    elif n_channels == 3:
        return (ai0, ai1, ai2)
    elif n_channels == 4:
        return (ai0, ai1, ai2, ai3)
    else:
        return 0 #TODO: Handle exception


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


def daq_bs(writing_done, dt, verbose):
    
    print("daq_bs() called.")
    print("daq_bs: event signalled.")
    # buffer1_full.set()
    # buffer1_full = True
    print("daq() finished.")


def daq(task, channels, buffer_size, i_buffer, dt, verbose):
    """ 
    Threaded data acquisition function.
    """
    if verbose >= 3:
        print("daq() called. Reading to buffer", i_buffer,"...")
        if verbose >= 4:
            stopwatch1 = time()
            
    cf = 0.0
    S=0
    t_start = time()
    
    if dt == 0:
        # Sample with unbounded sample rate
        while S < buffer_size:
            t_begin = time()
            
            samps = task.read(number_of_samples_per_channel=1)
            
            # t1 = time()
            buffers[i_buffer][0][S] = time()
            # t2 = time()
            buffers[i_buffer][1][S] = samps[0][0]
            # t3 = time()
            buffers[i_buffer][2][S] = samps[1][0]
            # t4 = time()
            buffers[i_buffer][3][S] = samps[2][0]
            # t5 = time()

            
            S += 1
            # t6 = time()
            # Idle until it is time for next sampling
            # t_left = time() - t_begin + cf
            # if t_left < dt:
            #     sleep(dt-t_left)
            # print("t1:",100_000*(t1-t_begin))
            # print("t2:",100_000*(t2-t1))
            # print("t3:",100_000*(t3-t2))
            # print("t4:",100_000*(t4-t3))
            # print("t5:",100_000*(t5-t4))
            # print("t6:",100_000*(t6-t5))
            
    
    elif dt > 0:
        # Sample with bounded sample rate
        while S < buffer_size:
            t_begin = time()
            
            samps = task.read(number_of_samples_per_channel=1)
            buffers[i_buffer][0][S] = time()
            buffers[i_buffer][1][S] = samps[0][0]
            buffers[i_buffer][2][S] = samps[1][0]
            buffers[i_buffer][3][S] = samps[2][0]
                
            S += 1
            
            # Idle until it is time for next sampling
            t_left = time() - t_begin + cf
            if t_left < dt:
                sleep(dt-t_left)
    else:
        pass #TODO: Handle exception

    if verbose >= 3:
        print("daq() finished.")
        if verbose >= 4:
            if dt == 0:
                print("Sampled",buffer_size,
                      "samples to buffer",i_buffer,
                      "at unbound rate in",1000*(time()-stopwatch1),"ms"
                      )
            else:
                print("Sampled",buffer_size,
                      "samples to buffer",i_buffer,
                      "at",round(1/dt),
                      "S/s in",1000*(time()-stopwatch1),"ms"
                      )
                
            t_avg = 1000*(time()-stopwatch1)/buffer_size
            if dt != 0:
                print("t_exp per sample:   ",dt*buffer_size,"ms")
            print("t_avg per sample:   ",t_avg,"ms")
            if dt != 0:
                print("Average slowdown:",100*(t_avg/(dt*buffer_size)-1),"%")
    
def write(buffer_size, i_buffer, filename, rounding_t, rounding_v, empty=False, verbose=0):
    """ 
    Threaded data writing function.
    """
    if verbose >= 3:
        print("write() called. Writing to buffer", i_buffer,"...")
        if verbose >= 4:
            stopwatch1 = time()
    
    with open(filename, 'a') as output_file:
        for buffer_line in range(buffer_size):

            # line = f"{round(data[0][i],tr)}"" {round(data[1][i],dr)} {round(data[2][i],dr)} {round(data[3][i],dr)}\n"
            # line = "%d %f %f %f".format(
            #     buffers[i_buffer][0][buffer_line],
            #     buffers[i_buffer][1][buffer_line],
            #     buffers[i_buffer][2][buffer_line],
            #     buffers[i_buffer][3][buffer_line])
            output_file.write("{} {} {} {}\n".format(
                round(buffers[i_buffer][0][buffer_line],rounding_t),
                round(buffers[i_buffer][1][buffer_line],rounding_v),
                round(buffers[i_buffer][2][buffer_line],rounding_v),
                round(buffers[i_buffer][3][buffer_line],rounding_v))
                )
        
    output_file.close()
    
    if verbose >= 4:
        stopwatch2 = time()
        
    if empty:
        buffers[i_buffer] = [array("d", [0]*buffer_size),
                             array("f", [0]*buffer_size), 
                             array("f", [0]*buffer_size), 
                             array("f", [0]*buffer_size)]  
    if verbose >= 4:
        stopwatch3 = time()

    if verbose >= 3:
        print("write() finished.")
        if verbose >= 4:
            print("Total time taken:          ",
                  round(1000*(time()-stopwatch1),3), "ms")
            print("Time spent writing:        ",
                  round(1000*(stopwatch2-stopwatch1),3), "ms")
            print("Time spent emptying buffer:",
                  round(1000*(stopwatch3-stopwatch2),3), "ms")
    
def main():
    N_SAMPLES = 10
    SR = 0
    BUFFER_SIZE = 5
    readout_period=0.0
    V_MIN=-1.25
    V_MAX=1.25
    verbose=4
    rounding_t = 6
    rounding_v = 6
    cf = 0.0009
    filetype=".dat"
    filename="test_data_"+str(round(time()))+filetype
    header = True
    
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
    
    # Initializing buffer objects:
    global buffer1, buffer2
    buffer1 = [array("d", [0]*BUFFER_SIZE), array("f", [0]*BUFFER_SIZE), 
               array("f", [0]*BUFFER_SIZE), array("f", [0]*BUFFER_SIZE)]
    buffer2 = [array("d", [0]*BUFFER_SIZE), array("f", [1]*BUFFER_SIZE), 
               array("f", [0]*BUFFER_SIZE), array("f", [0]*BUFFER_SIZE)]
    global buffers 
    buffers = [buffer1, buffer2]
    len_buffers = len(buffers)
    # Initializing event flags:
    # global buffer1_full, buffer2_full
    # buffer1_full = threading.Event()
    # buffer2_full = threading.Event()
    # buffer_full = (buffer1_full, buffer2_full)
    # global write_done
    # write_done = threading.Event()

    # Initializing event flags:
    # write_lock = threading.Lock()
    
    # th_daq  = threading.Thread(target=daq)
    # th_write = threading.Thread(target=write)
    
    # Initialize nidaqmx Task and run
    with nidaqmx.Task() as task:   
        
        (ai_x, ai_y, ai_z) = setup_channels(task, 3, V_MIN, V_MAX)
        # Specifying channels to analog input channels of USB-6008

        if verbose >= 4:
            print("Setup complete. Starting loop...")
        try:
            for i in range(n_loops):
                th_daq  = threading.Thread(
                    target=daq, 
                    args=(
                        task,
                        (ai_x, ai_y, ai_z),
                        BUFFER_SIZE,
                        i%len_buffers,
                        dt,
                        verbose
                        )
                    )
                th_write = threading.Thread(
                    target=write, 
                    args=(
                        BUFFER_SIZE,
                        (i+1)%len_buffers,
                        filename,
                        rounding_t,
                        rounding_v,
                        True,
                        verbose
                        )
                    )
                th_daq.daemon = True
                th_write.daemon = True
                
                th_daq.start()
                if i != 0:
                    th_write.start()
                    th_write.join() # New DAQ thread may not start before writing is done.
                th_daq.join() # Wait for DAQ to finish before starting new loop
        except (KeyboardInterrupt, SystemExit):
            print("Threads terminated by KeyboardInterrupt.")
            # This implementation may lead to import lock, see also:
            # https://stackoverflow.com/questions/46290045/import-silently-kills-thread/46354248#46354248
            # TODO: See if this requires addressing and if so, how exactly.
            
        if verbose >= 4:
            print("Program done.")              

#%% MAIN PROGRAM
if __name__ == "__main__":
    main()
    
    