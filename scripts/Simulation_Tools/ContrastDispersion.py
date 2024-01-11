"""
author: Anahita A. Seresti

Date: Oct24, 2023
"""

import sys
import os
import argparse
import matplotlib.pyplot as plt
import numpy as np


path = os.path.abspath("./")
sys.path.append(path)

from tools.ContrastTools import ReadResults, Slice, Model1D

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
        NFiles = 1000
        [MeshDict, N] = ReadResults(self.args.InputFolder, 'xdmf', NFiles) # >>> Reading the Mesh Files and storing them into a dictionary
        print('--- The number of files within the input folder is: ',N)
        print('-'*25)

        #* Getting Bounds of the mesh along the x-axis
        min_val = MeshDict[0].GetBounds()[0] # >>> x_min
        max_val = MeshDict[0].GetBounds()[1] # >>> x_max

        #* Defining the origin and normal of the slicer
        normal = (1., 0., 0.) # >>> i
        origin = (min_val, 0., 0.) # >>> Inflow
        origin_outlet = (max_val, 0., 0.) # >>> Outflow
        
        #* Getting time-attenuation curve at the inlet using cross-sectional average
        #* looping over the input meshs
        TempAttAve = [] # >>>Temporal Attenuation Average Value across the cross-section
        TempAttCent = [] # >>>Temporal Attenuation in the center of the cross-section
        TempAttCent_outlet = []
        
        print('--- Finding the temporal-attenuation curve by taking the inflow of the mesh')
        for mesh in MeshDict.values():
            
            #* taking the average value in the inlet slice
            [SectionAverage, CenterVal] = Slice(mesh, normal, origin) # >>> taking the inflow slice and reading the average and center value of the concentration array
            TempAttAve.append(SectionAverage) # >>> Updating Temporal Attenuation Average list
            TempAttCent.append(CenterVal) # >>> Updating Temporal Attenuation Center List
            [SectionAverage, CenterVal] = Slice(mesh, normal, origin_outlet)
            TempAttCent_outlet.append(CenterVal)
        

        desired_outlet_concentration = 0.22 # >>> 25% of the max input concentration
        
        for time, concentration in zip(np.arange(0,N,int(N/NFiles)), TempAttCent_outlet):
            if concentration >= 0:
                delay = time
            if concentration >= desired_outlet_concentration:
                peak = time
                break
            #elif concentration >= desired_outlet_concentration:
            #    print("WARNING: The running-time of the Advection-Diffusion simulation was not Enough for the outlet concentration to reach 25% of the max inlet concentration")
            #    peak = N-1
        CycleLength = self.args.CycleDuration/self.args.Increment # >>> Number of Time Steps in a Cycle
        print("---")
        print("The transition delay from inlet to outlet is: ", delay/CycleLength)
        print("Peak mesh (22% of the input maximum) happens at: ", peak/CycleLength)
        print("concentration at the reported time: ", concentration)
        print("---")
        PeakMesh = MeshDict[peak]

        

        
        #* Get space-attenuation curve
        NSlice = 50 # >>>Number of slices along the vessel
        interval = (max_val - min_val)/NSlice # >>> The interval between slices along the lumen
        SpatialAttAve = [] # >>>Spatial Attenuation Average Value across the cross-section
        SpatialAttCent = [] # >>>Spatial Attenuation in the center of the cross-section
        
        print('--- Finding the spatial-attenuation curve in the peak by slicing along the pipe')
        #* Moving along the peak pipe

        for i in range(NSlice):
            
            #* Define the cross-sectional origin along the pipe
            origin = (min_val + i*interval, 0, 0) # >>> Updating the origin of the slice: Moving along the lumen
            [SectionAverage, CenterVal] = Slice(PeakMesh, normal, origin) # >>> taking the slices and reading the average and center value of the concentration array
            SpatialAttAve.append(SectionAverage) # >>> Updating Spatial Attenuation Average List
            SpatialAttCent.append(CenterVal) # >>> Updating Spatial Attenuation Center List
        
        print('-'*25)

        #* fit linear model on data
        print('--- Fitting a linear model on temporal and spatial data')
        lag = 200 # >>> Time lag of the contrast diffusion
        lag_x = 20
        t = np.arange(0,N, int(N/NFiles))/CycleLength # >>>Time (s)
        t_new = t[lag:-lag] # >>> Applying time lag on the time array
        TempAttAve_new = np.array(TempAttAve[lag:-lag])
        TempAttCent_new = np.array(TempAttCent[lag:-lag])
        TempAttAve =  np.array(TempAttAve)
        TempAttCent =  np.array(TempAttCent)
        x = np.arange(min_val, max_val, interval) # >>> Space(cm)
        x_new = x[lag_x:]
        SpatialAttCent_new = SpatialAttCent[lag_x:]
        SpatialAttAve_new = SpatialAttAve[lag_x:]
        
        #`Linear Model (Lumen Mean)`
        [dCdtAve, TempPredAve] = Model1D(t_new, TempAttAve_new)
        [dCdxAve, SpacePredAve] = Model1D(x_new, SpatialAttAve_new)
        
        #`Linear Model (Lumen Center)`
        [dCdtCent, TempPredCent] = Model1D(t_new, TempAttCent_new)
        [dCdxCent, SpacePredCent] = Model1D(x_new, SpatialAttCent_new)
        
        #* Plotting Curves
        print('--- Plotting temporal and spatial attenuation curves')
        plt.figure(figsize=(13,5))
        
        #`temp`
        plt.subplot(121)
        #plt.scatter(t.reshape(-1,1), TempAttAve.reshape(-1,1), label = "Inflow Mean")
        #plt.plot(t_new.reshape(-1,1), TempPredAve.reshape(-1,1), color = 'green', label = "Linear (Inflow Mean)")
        plt.scatter(t.reshape(-1,1), TempAttCent.reshape(-1,1), label = "Inflow Center")
        plt.plot(t_new.reshape(-1,1), TempPredCent.reshape(-1,1), color = 'red', label = "Linear (Inflow Center)")
        plt.title("Temporal Attenuation Curve")
        plt.xlabel("time (s)")
        #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        
        #`space`
        plt.subplot(122)
        #plt.scatter(x.reshape(-1,1), np.array(SpatialAttAve).reshape(-1,1), label = "Lumen Mean")
        #plt.plot(x_new.reshape(-1,1), np.array(SpacePredAve).reshape(-1,1), color = 'green', label = "Linear (Lumen Mean)")
        plt.scatter(x.reshape(-1,1), np.array(SpatialAttCent).reshape(-1,1), label = "Lumen Center")
        plt.plot(x_new.reshape(-1,1), np.array(SpacePredCent).reshape(-1,1), color = 'red', label = "Linear (Linear Center)")
        plt.title("Spatial Attenuation Curve")
        plt.xlabel("length (cm)")
        #plt.legend(bbox_to_anchor=(1.05, 1), loc='best', borderaxespad=0.)
        
        print('--- Plotting Inlet Outlet Concentration')
        #`Inlet & Outlet vs time`
        plt.figure(figsize=(5,7))
        plt.plot(t.reshape(-1,1), TempAttCent.reshape(-1,1),color="black", label = "Inflow")
        plt.plot(t.reshape(-1,1), np.array(TempAttCent_outlet).reshape(-1,1),color="red", label = "Outflow")
        plt.title("Inlet and Outlet Concentration")
        plt.xlabel("time(s)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='best', borderaxespad=0.)

        plt.show()


        #* Writing the output files
        print('-'*25)
        print("--- generating the output time-concentration text files")

        parent_dir = os.path.dirname(self.args.InputFolder)
        with open(f'{parent_dir}/time_concentration.txt', 'w') as writefile:
            writefile.writelines('t, Inflow Contrast, Outflow Contrast\n')
            writefile.writelines(f'{t[i]}, {TempAttCent[i]}, {TempAttCent_outlet[i]}\n' for i in range(np.size(t)))

        print("--- generating the output distance-concentration text files")
        with open(f'{parent_dir}/distance_concentration.txt', 'w') as writefile:
            writefile.writelines('x, Contrast\n')
            writefile.writelines(f'{x[i]}, {SpatialAttCent[i]}\n' for i in range(np.size(x)))

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
