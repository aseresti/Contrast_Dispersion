"""
author: Anahita A. Seresti

Date: Oct24, 2023
"""

import numpy as np
import vtk
import sys
import os
path = os.path.abspath("./")
print(path)
sys.path.append(path)
from tools.utilities import vtk_to_numpy, ReadXDMFFile
import os
from tools.ContrastTools import get_order, ReadResults, AverageSlice, Model1D
import argparse
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

class ContrastDispersion():
    
    """
    This class takes the results_AdvectionDiffusion Folder and returns the 
    
    averaged velocity using Contrast Dispersion Method provided by Eslami 2022.
    """

    # takes argparse args as input
    def __init__(self,args):
        self.args = args
    
    def main(self):
        # Read xdmf files from results_AdvectionDiffusion Folder
        NFiles = 10
        [MeshDict, N] = ReadResults(self.args.InputFolder, 'xdmf', NFiles)

        # Getting Bounds of the mesh 
        min_val = MeshDict[0].GetBounds()[0]
        max_val = MeshDict[0].GetBounds()[1]
        
        # Defining the origin and normal of the slicer
        normal = (1., 0., 0.)
        origin = (min_val, 0., 0.)
        
        # Getting time-attenuation curve at the inlet using cross-sectional average
        # looping over the input meshs
        TemporalAttenuation = []
        for mesh in MeshDict.values():
            
            # taking the average value in the inlet slice
            [average, slice] = AverageSlice(mesh, normal, origin)
            TemporalAttenuation.append(average)
            # ! exploring slice attributes 
            print(dir(slice.GetOutputPort()))
            print(slice.GetOutput().GetPointData().GetArray(0).GetValue(1))
            # TODO: Find  the pointid of the origin coordinate and find the array value at that pointid
            # print(slice.GetPoint(i for i in origin))
            exit(1)
        
        # fit linear model on the temporal curve
        plt.figure(1)
        t = np.arange(0,N,int(N/NFiles))/1000
        [plt, dc_dt] = Model1D(t, TemporalAttenuation)
        
        # Get space-attenuation curve
        NSlice = 50
        interval = (max_val - min_val)/NSlice
        SpatialAttenuation = []

        # Moving along the pipe
        for i in range(NSlice):

            # Define the cross-sectional origin along the pipe
            origin = (min_val + i*interval)
            [average, slice] = AverageSlice(MeshDict[NFiles-1], normal,origin)
            SpatialAttenuation.append(average)

            # TODO: Get the centerline value

        # fit linear model on the spatial data
        x = np.arange(min_val, max_val, interval)
        plt.figure(2)
        [plt, dc_dx] = Model1D(x, SpatialAttenuation)

        # TODO: Output four figures of time and space for centerline and cross-sectional average
        
        # return velocity
        velocity_sectional = abs(dc_dt/dc_dx)
        return velocity_sectional



if __name__=="__main__":
	# Define argparse arguments
    parser = argparse.ArgumentParser(description="This script implement the contrast dispersion pipeline")
    
    # Take Input Folder path
    parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest='InputFolder')
    # Define the duration of each cycle
    parser.add_argument('-CycleDuration', '--CycleDuration', type=int, required=False, dest='CycleDuration', default=1000)
    # Parse arguments
    args = parser.parse_args()
    # call the class
    velocity_sectional = ContrastDispersion(args).main()
    print(velocity_sectional)
