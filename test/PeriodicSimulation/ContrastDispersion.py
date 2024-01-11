"""
author: Anahita A. Seresti

Date: Dec14, 2023

Adapted for periodic simulations
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
    
    def uprise(self):
        """
        `Selects the uprising part of the contrast inflow
        """
        total_time = 120
        rise_time = 8
        start_time_1 = 1
        start_time_2 = 40
        start_time_3 = 80
        
        #* Read xdmf files from results_AdvectionDiffusion Folder
        print('--- Reading the input meshes within the input folder')
        self.NFiles = 1200 # >>> Number of temporal samples
        [MeshDict, self.N, PeakMesh] = ReadResults(self.args.InputFolder, 'xdmf', self.NFiles) # >>> Reading the Mesh Files and storing them into a dictionary
        print('--- The number of files within the input folder is: ',self.N)
        print('-'*25)

        mapping = int(self.NFiles/total_time)
        batch_1 = batch_2 = batch_3 = dict()
        for i in range(rise_time*mapping):
            batch_1[i] = MeshDict[mapping*start_time_1 + i]
            batch_2[i] = MeshDict[mapping*start_time_2 + i]
            batch_3[i] = MeshDict[mapping*start_time_3 + i]
        
        
        self.Peak1 = MeshDict[mapping*(rise_time+start_time_1)]
        self.Peak2 = MeshDict[mapping*(rise_time+start_time_2)]
        self.Peak3 = MeshDict[mapping*(rise_time+start_time_3)]
        #print(mapping*(rise_time+start_time_1), mapping*(rise_time+start_time_2), mapping*(rise_time+start_time_3))

        self.CycleLength = self.args.CycleDuration/self.args.Increment # >>> Number of Time Steps in a Cycle
        self.t1 = np.arange(mapping*start_time_1,mapping*(rise_time+start_time_1))/self.CycleLength*10
        self.t2 = np.arange(mapping*start_time_2,mapping*(rise_time+start_time_2))/self.CycleLength*10
        self.t3 = np.arange(mapping*start_time_3,mapping*(rise_time+start_time_3))/self.CycleLength*10
        

        return MeshDict, batch_1, batch_2, batch_3


    def slicer(self, MeshDict, PeakMesh):
        """
        `slice the mesh in inlet and along the centerline
        """

        #* Getting Bounds of the mesh along the x-axis
        self.min_val = MeshDict[0].GetBounds()[0] # >>> x_min
        self.max_val = MeshDict[0].GetBounds()[1] # >>> x_max

        #* Defining the origin and normal of the slicer
        normal = (1., 0., 0.) # >>> i
        origin = (self.min_val, 0., 0.) # >>> Inflow
        origin_outlet = (self.max_val, 0., 0.) # >>> Outflow

        #* Getting time-attenuation curve at the inlet using cross-sectional average
        #* looping over the input meshs
        TempAttCent = [] # >>>Temporal Attenuation in the center of the cross-section
        TempAttCent_outlet = []
        
        print('-'*25)
        print('--- Finding the temporal-attenuation curve by taking the inflow of the mesh')
        for mesh in MeshDict.values():
            
            #* taking the average value in the inlet slice
            [SectionAverage, CenterVal] = Slice(mesh, normal, origin)
            TempAttCent.append(CenterVal) # >>> Updating Temporal Attenuation Center List
            [SectionAverage, CenterVal] = Slice(mesh, normal, origin_outlet)
            TempAttCent_outlet.append(CenterVal)
        
        #* Get space-attenuation curve
        NSlice = 50 # >>>Number of slices along the vessel
        self.interval = (self.max_val - self.min_val)/NSlice # >>> The interval between slices along the lumen
        SpatialAttCent = [] # >>>Spatial Attenuation in the center of the cross-section
        
        print('--- Finding the spatial-attenuation curve in the peak by slicing along the pipe')
        #* Moving along the peak pipe

        for i in range(NSlice):
            
            #* Define the cross-sectional origin along the pipe
            origin = (self.min_val + i*self.interval, 0, 0) # >>> Updating the origin of the slice: Moving along the lumen
            [SectionAverage, CenterVal] = Slice(PeakMesh, normal, origin) # >>> taking the slices and reading the average and center value of the concentration array
            SpatialAttCent.append(CenterVal) # >>> Updating Spatial Attenuation Center List
        
        print('-'*25)

        return TempAttCent, TempAttCent_outlet, SpatialAttCent
    

    def fit(self, TempAttCent, SpatialAttCent, t):
        """
        `fit linear model on temporal and spatial data
        """
        
        TempAttCent =  np.array(TempAttCent)
        SpatialAttCent =  np.array(SpatialAttCent)
        self.x = np.arange(self.min_val, self.max_val, self.interval) # >>> Space(cm)
        
        #* Linear Model (Lumen Center)
        print('--- Fitting a linear model on temporal data')
        [dCdtCent, TempPredCent]  = Model1D(t, TempAttCent)

        print('--- Fitting a linear model on spatial data')
        self.xlag = 20
        [dCdxCent, SpacePredCent] = Model1D(self.x[self.xlag:], SpatialAttCent[self.xlag:])

        return dCdtCent, TempPredCent, dCdxCent, SpacePredCent

    def main(self):
        
        #* Read xdmf files from results_AdvectionDiffusion Folder
        [MeshDict, batch_1, batch_2, batch_3] = self.uprise()
        
        
        [TempAttCent, TempAttCent_outlet, SpatialAttCent] = self.slicer(MeshDict,self.Peak3)
        
        t = np.arange(0,self.N,int(self.N/self.NFiles))/self.CycleLength # >>>Time (s)
        
        #* Plotting Curves
        print('--- Plotting Inlet Outlet Concentration')
        #`Inlet & Outlet vs time`
        plt.figure(figsize=(5,7))
        plt.plot(t.reshape(-1,1), np.array(TempAttCent).reshape(-1,1),color="black", label = "Inflow")
        plt.plot(t.reshape(-1,1), np.array(TempAttCent_outlet).reshape(-1,1),color="red", label = "Outflow")
        plt.title("Inlet and Outlet Concentration")
        plt.xlabel("time(s)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='best', borderaxespad=0.)

        plt.show()

        print('--- Plotting temporal and spatial attenuation curves')
        
        plt.figure(figsize=(13,5))
        
        #`temp`
        plt.subplot(121)
        plt.scatter(t.reshape(-1,1), np.array(TempAttCent).reshape(-1,1), label = "Inflow Center")
        plt.title("Temporal Attenuation Curve")
        plt.xlabel("time (s)")

        #>>> Period = 1
        [TempAttCent, TempAttCent_outlet, SpatialAttCent] = self.slicer(batch_1,self.Peak1)
        [dCdtCent, TempPredCent, dCdxCent, SpacePredCent] = self.fit(TempAttCent, SpatialAttCent, t=self.t1)
        plt.plot(self.t1.reshape(-1,1), np.array(TempPredCent).reshape(-1,1), color = 'red', label = "Linear (Inflow Center)")
        plt.scatter(self.t1[-1], TempAttCent[-1], color = 'orange')

        #`space`
        plt.subplot(122)
        plt.scatter(self.x.reshape(-1,1), np.array(SpatialAttCent).reshape(-1,1), color = "blue", label = "Lumen Center")
        plt.plot(self.x[self.xlag:].reshape(-1,1), np.array(SpacePredCent).reshape(-1,1),  color = "blue", label = "Linear (Linear Center)")
        plt.title("Spatial Attenuation Curve")
        plt.xlabel("length (cm)")

        #* Writing the output files
        print('-'*25)
        print("--- generating the output time-concentration text files")

        parent_dir = os.path.dirname(self.args.InputFolder)
        with open(f'{parent_dir}/time_concentration_1.txt', 'w') as writefile:
            writefile.writelines('t, Inflow Contrast, Outflow Contrast\n')
            writefile.writelines(f'{self.t1[i]}, {TempAttCent[i]} \n' for i in range(np.size(self.t1)))

        print("--- generating the output distance-concentration text files")
        with open(f'{parent_dir}/distance_concentration_1.txt', 'w') as writefile:
            writefile.writelines('x, Contrast\n')
            writefile.writelines(f'{self.x[i]}, {SpatialAttCent[i]}\n' for i in range(np.size(self.x)))

        velocity_p1 =  abs(dCdtCent/dCdxCent)
        
        #>>> Period = 2
        [TempAttCent, TempAttCent_outlet, SpatialAttCent] = self.slicer(batch_2,self.Peak2)
        [dCdtCent, TempPredCent, dCdxCent, SpacePredCent] = self.fit(TempAttCent, SpatialAttCent, t=self.t2)
        plt.subplot(121)
        plt.plot(self.t2.reshape(-1,1), np.array(TempPredCent).reshape(-1,1), color = 'red', label = "Linear (Inflow Center)")
        plt.scatter(self.t2[-1], TempAttCent[-1], color = 'orange')
        
        #`space`
        plt.subplot(122)
        plt.scatter(self.x.reshape(-1,1), np.array(SpatialAttCent).reshape(-1,1), color = "red", label = "Lumen Center")
        plt.plot(self.x[self.xlag:].reshape(-1,1), np.array(SpacePredCent).reshape(-1,1), color = "red", label = "Linear (Linear Center)")
        plt.title("Spatial Attenuation Curve")
        plt.xlabel("length (cm)")

        #* Writing the output files
        print('-'*25)
        print("--- generating the output time-concentration text files")

        parent_dir = os.path.dirname(self.args.InputFolder)
        with open(f'{parent_dir}/time_concentration_2.txt', 'w') as writefile:
            writefile.writelines('t, Inflow Contrast, Outflow Contrast\n')
            writefile.writelines(f'{self.t2[i]}, {TempAttCent[i]} \n' for i in range(np.size(self.t2)))

        print("--- generating the output distance-concentration text files")
        with open(f'{parent_dir}/distance_concentration_2.txt', 'w') as writefile:
            writefile.writelines('x, Contrast\n')
            writefile.writelines(f'{self.x[i]}, {SpatialAttCent[i]}\n' for i in range(np.size(self.x)))

        velocity_p2 =  abs(dCdtCent/dCdxCent)

        #>>> Period = 3
        [TempAttCent, TempAttCent_outlet, SpatialAttCent] = self.slicer(batch_3,self.Peak3)
        [dCdtCent, TempPredCent, dCdxCent, SpacePredCent] = self.fit(TempAttCent, SpatialAttCent, t=self.t3)
        plt.subplot(121)
        plt.plot(self.t3.reshape(-1,1), np.array(TempPredCent).reshape(-1,1), color = 'red', label = "Linear (Inflow Center)")
        plt.scatter(self.t3[-1], TempAttCent[-1], color = 'orange')

        #`space`
        plt.subplot(122)
        plt.scatter(self.x.reshape(-1,1), np.array(SpatialAttCent).reshape(-1,1), color = "green", label = "Lumen Center")
        plt.plot(self.x[self.xlag:].reshape(-1,1), np.array(SpacePredCent).reshape(-1,1), color = "green", label = "Linear (Linear Center)")
        plt.title("Spatial Attenuation Curve")
        plt.xlabel("length (cm)")
        plt.show()


        #* Writing the output files
        print('-'*25)
        print("--- generating the output time-concentration text files")

        parent_dir = os.path.dirname(self.args.InputFolder)
        with open(f'{parent_dir}/time_concentration_3.txt', 'w') as writefile:
            writefile.writelines('t, Inflow Contrast, Outflow Contrast\n')
            writefile.writelines(f'{self.t3[i]}, {TempAttCent[i]} \n' for i in range(np.size(self.t3)))

        print("--- generating the output distance-concentration text files")
        with open(f'{parent_dir}/distance_concentration_3.txt', 'w') as writefile:
            writefile.writelines('x, Contrast\n')
            writefile.writelines(f'{self.x[i]}, {SpatialAttCent[i]}\n' for i in range(np.size(self.x)))

        velocity_p3 =  abs(dCdtCent/dCdxCent)
        
        #* return velocity
        
        return velocity_p1, velocity_p2, velocity_p3



if __name__=="__main__":
	#* Define argparse arguments
    parser = argparse.ArgumentParser(description="This script implement the contrast dispersion pipeline")
    
    #* Take Input Folder path
    parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest='InputFolder')
    #* Define the duration of each cycle
    parser.add_argument('-CycleDuration', '--CycleDuration', type=int, required=False, dest='CycleDuration', default=1000)
    #* Define the Increment
    parser.add_argument('-Increment', '--Increment', type=int, required=False, dest='Increment', default=10)
    #* Parse arguments
    args = parser.parse_args()
    
    #* call the class
    [velocity_p1, velocity_p2, velocity_p3] = ContrastDispersion(args).main()
    print(velocity_p1, velocity_p2, velocity_p3)
    
