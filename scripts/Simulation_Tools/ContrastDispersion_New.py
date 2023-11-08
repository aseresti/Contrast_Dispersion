"""
author: Anahita A. Seresti

Date: Oct24, 2023
"""

import numpy as np
import sys
import os
path = os.path.abspath("./")
sys.path.append(path)
from tools.ContrastTools import ReadResults, Slice, Model1D
import argparse
import matplotlib.pyplot as plt
import math

class ContrastDispersion():
    
    """
    `This class takes the results_AdvectionDiffusion Folder and returns the 
    
    `averaged velocity using Contrast Dispersion Method provided by Eslami 2022.
    """

    #* takes argparse args as input
    def __init__(self,args):
        self.args = args
    
    def main(self):
        
        #* Read xdmf files from results_AdvectionDiffusion Folder
        print('--- Reading the input meshes within the input folder')
        NFiles = 20 # >>> Number of temporal samples
        [MeshDict, N] = ReadResults(self.args.InputFolder, 'xdmf', NFiles) # >>> Reading the Mesh Files and storing them into a dictionary
        print('--- The number of read files within the input folder is: ',N)
        print('-'*25)

        #* Getting Bounds of the mesh along the x-axis
        min_val = MeshDict[0].GetBounds()[0] # >>> x_min
        max_val = MeshDict[0].GetBounds()[1] # >>> x_max

        #* Defining the origin and normal of the slicer
        normal = (1., 0., 0.) # >>> i
        origin = (min_val, 0., 0.) # >>> Inflow
        
        #* Getting time-attenuation curve at the inlet using cross-sectional average
        #* looping over the input meshs
        TempAttAve = [] # >>>Temporal Attenuation Average Value across the cross-section
        TempAttCent = [] # >>>Temporal Attenuation in the center of the cross-section
        
        print('--- Finding the temporal-attenuation curve by taking the inflow of the mesh')
        for mesh in MeshDict.values():
            
            #* taking the average value in the inlet slice
            [SectionAverage, CenterVal] = Slice(mesh, normal, origin) # >>> taking the inflow slice and reading the average and center value of the concentration array
            TempAttAve.append(SectionAverage) # >>> Updating Temporal Attenuation Average list
            TempAttCent.append(CenterVal) # >>> Updating Temporal Attenuation Center List
        
        #* Get space-attenuation curve
        NSlice = 50 # >>>Number of slices along the vessel
        interval = (max_val - min_val)/NSlice # >>> The interval between slices along the lumen
        SpatialAttAve = [] # >>>Spatial Attenuation Average Value across the cross-section
        SpatialAttCent = [] # >>>Spatial Attenuation in the center of the cross-section
        
        print('--- Finding the spatial-attenuation curve in the peak by slicing along the pipe')
        #* Moving along the pipe
        peak = math.ceil(NFiles/20*16)
        #peak = NFiles-1

        for i in range(NSlice):
            
            #* Define the cross-sectional origin along the pipe
            origin = (min_val + i*interval, 0, 0) # >>> Updating the origin of the slice: Moving along the lumen
            [SectionAverage, CenterVal] = Slice(MeshDict[peak], normal, origin) # >>> taking the slices and reading the average and center value of the concentration array
            SpatialAttAve.append(SectionAverage) # >>> Updating Spatial Attenuation Average List
            SpatialAttCent.append(CenterVal) # >>> Updating Spatial Attenuation Center List
        
        print('-'*25)

        #* fit linear model on data
        print('--- Fitting a linear model on temporal and spatial data')
        lag = 6 # >>> Time lag of the contrast diffusion
        CycleLength = self.args.CycleDuration/self.args.Increment # >>> Number of Time Steps in a Cycle
        t = np.arange(0,N,int(N/NFiles))/CycleLength # >>>Time (s)
        t = t[lag:] # >>> Applying time lag on the time array
        TempAttAve = np.array(TempAttAve[lag:])
        TempAttCent = np.array(TempAttCent[lag:])
        x = np.arange(min_val, max_val, interval) # >>> Space(cm)
        
        #`Linear Model (Lumen Mean)`
        [dCdtAve, TempPredAve] = Model1D(t, TempAttAve)
        [dCdxAve, SpacePredAve] = Model1D(x, SpatialAttAve)
        
        #`Linear Model (Lumen Center)`
        [dCdtCent, TempPredCent] = Model1D(t, TempAttCent)
        [dCdxCent, SpacePredCent] = Model1D(x, SpatialAttCent)
        
        #* Plotting Curves
        print('--- Plotting temporal and spatial attenuation curves')
        plt.figure(figsize=(13,5))
        
        #`temp`
        plt.subplot(121)
        #plt.scatter(t.reshape(-1,1), TempAttAve.reshape(-1,1), label = "Inflow Mean")
        #plt.plot(t.reshape(-1,1), TempPredAve.reshape(-1,1), color = 'red', label = "Linear (Inflow Mean)")
        plt.scatter(t.reshape(-1,1), TempAttCent.reshape(-1,1), label = "Inflow Center")
        plt.plot(t.reshape(-1,1), TempPredCent.reshape(-1,1), color = 'red', label = "Linear (Inflow Center)")
        plt.title("Temporal Attenuation Curve")
        plt.xlabel("time (s)")
        #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        
        #`space`
        plt.subplot(122)
        #plt.scatter(x.reshape(-1,1), np.array(SpatialAttAve).reshape(-1,1), label = "Lumen Mean")
        #plt.plot(x.reshape(-1,1), np.array(SpacePredAve).reshape(-1,1), color = 'red', label = "Linear (Lumen Mean)")
        plt.scatter(x.reshape(-1,1), np.array(SpatialAttCent).reshape(-1,1), label = "Lumen Center")
        plt.plot(x.reshape(-1,1), np.array(SpacePredCent).reshape(-1,1), color = 'red', label = "Linear (Linear Center)")
        plt.title("Spatial Attenuation Curve")
        plt.xlabel("length (cm)")
        #plt.legend(bbox_to_anchor=(1.05, 1), loc='best', borderaxespad=0.)
        plt.show()

        #* return velocity
        VelocityAve = abs(dCdtAve/dCdxAve)
        VelocityCent = abs(dCdtCent/dCdxCent)
        return VelocityAve, VelocityCent



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
    
    #* call the class
    [VelocityAve, VelocityCent] = ContrastDispersion(args).main()
    print('The sectional average velocity is: ', VelocityAve)
    print('The centerline velocity is:', VelocityCent)
