#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 13:02:53 2023

@author: ana
"""



import numpy as np
import glob
import argparse
from math import sqrt, atan, pi
from utilities import ReadVTPFile
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

class ImageAnalysisAdvectionDiffusionAlongCL():
    def __init__(self,Args):
        self.Args = Args
        
    def MAFilter(self,input_):
        
        window_length = 1
        FilteredSignal = np.zeros(np.size(input_)-window_length)
        for point in range(np.size(input_)-window_length):
            FilteredSignal[point] = np.average(input_[point:point+window_length])
        return FilteredSignal
        
    def Main(self):
        
        #Store the cl file names inside the input folder
        FileName = self.Args.InputFile
                
        #Read the averaged pixel values of each cl file and storing them into an excel sheet
        CLFile = ReadVTPFile(FileName)
        PixelValArray = CLFile.GetPointData().GetArray("AvgPixelValue")
        
        #RadiusArray = CLFile.GetPointData().GetArray("MaximumInscribedSphereRadius")

        #D = 0.01
        point_1 = CLFile.GetPoint(1)
        point_2 = CLFile.GetPoint(2)
        x_step = sqrt((point_1[0] - point_2[0])**2+(point_1[1] - point_2[1])**2+(point_1[2] - point_2[2])**2) #Calculate the distance between the CL points
        
        print(x_step)
        
        
        x_MAFiltered = ImageAnalysisAdvectionDiffusionAlongCL(args).MAFilter(np.transpose(PixelValArray))
        nx = np.size(x_MAFiltered)
        x_arr = np.arange(0, nx*x_step, x_step)
        model = LinearRegression()
        model.fit(x_arr.reshape(-1, 1),x_MAFiltered.reshape(-1, 1))
        dc_dx = model.coef_[0][0]
        plt.scatter(x_arr, x_MAFiltered, color = 'black')
        pred = model.predict(x_arr.reshape(-1, 1))
        plt.plot(x_arr.reshape(-1, 1), pred.reshape(-1, 1), color = 'red')
        plt.text(11, 1240, f"slope = {int(model.coef_[0][0]*100)/100}", rotation = -15)
        plt.title("Averaged Contrast over Space")
        plt.xlabel("distance from inflow cap (mm)")
        plt.ylabel("Contrast (HU)")
        #plt.axis([0, 30, 150, 200])
        plt.show()
        #velocity = dc_dt/dc_dx
        #print(abs(velocity), 'mm/sec')
        print(dc_dx)
        
        
if __name__ == "__main__":
    #Description
    parser = argparse.ArgumentParser(description = "This script gets the centerline files with the averaged pixel value array and calculates the advection diffusion velocity")
    
    #InputFiles
    parser.add_argument('-InputFile', '--InputFile', type = str, required = True, dest = "InputFile", help = "The Centerline Files with averaged pixel values")
    
    args = parser.parse_args()
    ImageAnalysisAdvectionDiffusionAlongCL(args).Main()
    