# -*- coding: utf-8 -*-

# Source:
# 

# from PyDAQmx import *
from time import time
import numpy as np
import PyDAQmx
import ctypes

t0 = time()
# Declaration of Task handle object
taskHandle = PyDAQmx.TaskHandle()

# TODO: Clarify this step
read = PyDAQmx.DAQmxTypes.int32()
data = np.zeros((300,), dtype=np.float64)

t1 = time() #OH: 0 ns

try:
    # DAQmx Configure Code
    # Create the task object, give it a name, and pass the taskhandle object
    #   by reference (using byref() from ctypes)
    PyDAQmx.DAQmxCreateTask("Measure_ai0",ctypes.byref(taskHandle))
    
    t2 = time() #OH: 4000 ns
    
    # Create a channel to measure voltage, and add the channel to the task
    #   specified by the taskhandle.
    PyDAQmx.DAQmxCreateAIVoltageChan(
        taskHandle,                 # TaskHandle object to link task
        "Dev1/ai0",                 # Name of the physical channel
        "",                         # Name assigned to virtual channel
        PyDAQmx.DAQmx_Val_Diff,     # Input terminal config -> look at https://documentation.help/NI-DAQmx-C-Functions/DAQmxCreateAIVoltageChan.html
        -1.25,                      # minVal of voltage in <units>
        1.25,                       # maxVal of voltage in <units>
        PyDAQmx.DAQmx_Val_Volts,    # Meaning of units
        None                        # Name of custom scale to apply to channel
        )
    
    PyDAQmx.DAQmxCreateAIVoltageChan(
        taskHandle,                 # TaskHandle object to link task
        "Dev1/ai1",                 # Name of the physical channel
        "",                         # Name assigned to virtual channel
        PyDAQmx.DAQmx_Val_Diff,     # Input terminal config -> look at https://documentation.help/NI-DAQmx-C-Functions/DAQmxCreateAIVoltageChan.html
        -1.25,                      # minVal of voltage in <units>
        1.25,                       # maxVal of voltage in <units>
        PyDAQmx.DAQmx_Val_Volts,    # Meaning of units
        None                        # Name of custom scale to apply to channel
        )
    
    PyDAQmx.DAQmxCreateAIVoltageChan(
        taskHandle,                 # TaskHandle object to link task
        "Dev1/ai2",                 # Name of the physical channel
        "",                         # Name assigned to virtual channel
        PyDAQmx.DAQmx_Val_Diff,     # Input terminal config -> look at https://documentation.help/NI-DAQmx-C-Functions/DAQmxCreateAIVoltageChan.html
        -1.25,                      # minVal of voltage in <units>
        1.25,                       # maxVal of voltage in <units>
        PyDAQmx.DAQmx_Val_Volts,    # Meaning of units
        None                        # Name of custom scale to apply to channel
        )
    
    
    t3 = time() #OH: 2000 ns
    
    PyDAQmx.DAQmxCfgSampClkTiming(
        taskHandle,                     # TaskHandle object to link task
        "",                             # The source terminal of the sample clock (if nothing, it uses internal device clock)
        10.0,                           # The sampling rate in samples per second per channel
        PyDAQmx.DAQmx_Val_Rising,       # Active edge, essentially whether samples are gathered at the rising or falling edge of the clock signal
        PyDAQmx.DAQmx_Val_ContSamps,    # Specify whether finite number of samples are gathered, whether continuous samples are gathered until the task is stopped, or whether to use hardware timed single point sample mode.
        100                             # Number of samples to acquire or generate for EACH channel in the task (only if FiniteSamps mode is on). If ContSamps is on, this value determines buffer size.
        )


    t4 = time() # OH: 0 ns
    
    # DAQmx Start Code
    PyDAQmx.DAQmxStartTask(taskHandle)  # Transitions the task from the committed state to the running state, which begins measurement or generation.

    t5 = time() #OH: 1000 ns
    
    # DAQmx Read Code
    PyDAQmx.DAQmxReadAnalogF64(         # Reads multiple floating-point samples from a task that contains one or more analog input channels.
        taskHandle,                     # TaskHandle object to link task
        100,                            # The number of samples, per channel, to read. If readArray does not contain enough space, this function returns as many samples as fit in readArray.
        15.0,                           # The timeout in seconds.
        PyDAQmx.DAQmx_Val_GroupByChannel, # fillMode: interleaving is OFF with DAQmx_Val_GroupByChannel and ON with DAQmx_Val_GroupByScanNumber
        data,                           # readArray?
        300,                            # uInt32 size of the readArray (why is it here??)
        ctypes.byref(read),             # data type of readArray passed by reference?
        None                            # Reserved field?
        )
    
    t6 = time() # OH: ~sr*n_samples_channel
    print(f"Acquired {read.value} points\n")
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

t = [t0, t1, t2, t3, t4, t5, t6, t7, t8, t9]

for i in range(len(t)-1):
    print(str(i+1), ":", str(round(100_000*(t[i+1]-t[i]))))

for i in range(int(len(data)/3)):
    print("{} {} {}".format(data[i], data[i+100], data[i+200]))
