#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 13:46:57 2023

@author: ana
"""


import numpy as np
import glob
import argparse
from math import sqrt, atan, pi
from utilities import ReadVTPFile
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, ifft

class ImageAnalysisAdvectionDiffusionAlongCL():
    def __init__(self,Args):
        self.Args = Args
        
    def MAFilter(self,input_):
        
        window_length = 8
        FilteredSignal = np.zeros(np.size(input_)-window_length)
        for point in range(np.size(input_)-window_length):
            FilteredSignal[point] = np.average(input_[point:point+window_length])
        return FilteredSignal
        
    def Main(self):
        
        #Store the cl file names inside the input folder
        FolderName = self.Args.InputFolder
        FileNames  = glob.glob(f'{FolderName}/*.vtp')
        FileNames  = sorted(FileNames)
        NFiles = FileNames.__len__()
        Npoints = ReadVTPFile(FileNames[0]).GetNumberOfPoints()
        
        #Read the averaged pixel values of each cl file and storing them into an excel sheet
        PixelValArray = np.zeros((Npoints, NFiles))
        for i in range(NFiles):
            CLFile = ReadVTPFile(FileNames[i])
            PixelValArray[:,i] = CLFile.GetPointData().GetArray("AvgPixelValue")
        
        #RadiusArray = CLFile.GetPointData().GetArray("MaximumInscribedSphereRadius")

        # Taking the linear part of the lumen and the upslope samples
        upslope_delay = 3
        upslope_peak = 7
        upslope = PixelValArray[0:-3, upslope_delay:upslope_peak] 
        
        #Calculating Velocity
        #D = 0.01
        point_1 = CLFile.GetPoint(1)
        point_2 = CLFile.GetPoint(2)
        x_step = sqrt((point_1[0] - point_2[0])**2+(point_1[1] - point_2[1])**2+(point_1[2] - point_2[2])**2) #Calculate the distance between the CL points
        
        #print(x_step)
        time_step = 3.2#4 #assuming the heart-rate to be 60-per-min and the pictures are taken every 4-cycle
        (nx, nt) = upslope.shape
        (fx, ft) = (fftfreq(nx, x_step), fftfreq(nt, time_step))
            
        dC_dt = 1j*2*pi*ft*fft(upslope)
                
        dC_dx = 1j*2*pi*fx*(fft(upslope.transpose()))
                
        v = abs(ifft(dC_dt))/abs(ifft(dC_dx).transpose())
        print(v[1:,:])
        print('maximum velocity is', np.amax(v[1:,:])/10, 'cm/s')
        print('average velocity is', np.average(v[1:,:])/10, 'cm/s')
        print('median velocity', np.median(v[1:,:])/10, 'cm/s')        
        
if __name__ == "__main__":
    #Description
    parser = argparse.ArgumentParser(description = "This script gets the centerline files with the averaged pixel value array and calculates the advection diffusion velocity")
    
    #InputFiles
    parser.add_argument('-InputFolder', '--InputFolder', type = str, required = True, dest = "InputFolder", help = "The Folder containing the centerline files with averaged pixel values")
    
    args = parser.parse_args()
    ImageAnalysisAdvectionDiffusionAlongCL(args).Main()
    