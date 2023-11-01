"""
author: Anahita AS
date: Oct 31, 2023
"""
import sys
import os
path = os.path.abspath("./")
sys.path.append(path)
from tools.ContrastTools import ReadResults, Slice
import argparse
import matplotlib.pyplot as plt
import numpy as np

def main(Args):
    
    """
    ` the goal of this function is to plot the inflow contrast function assigned to the input of the pipe
    """

    #* Read xdmf files from results_AdvectionDiffusion Folder
    NFiles = 100 # >>> Number of temporal samples
    [MeshDict, N] = ReadResults(Args.InputFolder, 'xdmf', NFiles) # >>> Reading the Mesh Files and storing them into a dictionary
    print(N)
    #* Getting Bounds of the mesh along the x-axis
    min_val = MeshDict[0].GetBounds()[0] # >>> x_min
    
    #* Defining the origin and normal of the slicer
    normal = (1., 0., 0.) # >>> i
    origin = (min_val, 0., 0.) # >>> Inflow
    
    #* Getting time-attenuation curve at the inlet using cross-sectional average
    #* looping over the input meshs
    TempAttAve = [] # >>>Temporal Attenuation Average Value across the cross-section
    
    for mesh in MeshDict.values():
        
        #* taking the average value in the inlet slice
        [SectionAverage, CenterVal] = Slice(mesh, normal, origin) # >>> taking the inflow slice and reading the average and center value of the concentration array
        TempAttAve.append(SectionAverage) # >>> Updating Temporal Attenuation Average list
        
    CycleLength = Args.CycleDuration/Args.Increment # >>> Number of Time Steps in a Cycle
    t = np.arange(0,N,int(N/NFiles))/CycleLength # >>>Time (s)

    #* Plotting inflow
    plt.figure()
    plt.plot(t.reshape(-1,1), np.array(TempAttAve).reshape(-1,1))
    plt.xlabel('time (s)')
    plt.ylabel('Inflow Contrast')
    plt.title('Contrast Agent Inflow Function')
    plt.show()


if __name__=="__main__":
	#* Define argparse arguments
    parser = argparse.ArgumentParser(description="This script implement the contrast dispersion pipeline")
    
    #* Take Input Folder path
    parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest='InputFolder')
    #* Define the duration of each cycle
    parser.add_argument('-CycleDuration', '--CycleDuration', type=int, required=False, dest='CycleDuration', default=1000)
    #* Define the Increment
    parser.add_argument('-Increment', '--Increment', type=int, required=False, dest='Increment', default=20)
    #* Parse arguments
    args = parser.parse_args()
    
    main(args)