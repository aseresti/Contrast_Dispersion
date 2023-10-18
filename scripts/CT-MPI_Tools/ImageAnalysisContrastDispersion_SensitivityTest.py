#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 08:11:52 2023

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
    '''      
    def MAFilter(self,input_):
        
        window_length = 8
        FilteredSignal = np.zeros(np.size(input_)-window_length)
        for point in range(np.size(input_)-window_length):
            FilteredSignal[point] = np.average(input_[point:point+window_length])
        return FilteredSignal
    '''        
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
        v_length = PixelValArray.shape[0]
        # Taking the linear part of the lumen and the upslope samples
        upslope_delay = 3
        upslope_peak = 7
        
        upslope = PixelValArray[0:int(v_length/5), upslope_delay:upslope_peak] 
        fullarray = PixelValArray[0:int(v_length/5), :]
        #print(x_step)
        time_step = 3.2#4 #assuming the heart-rate to be 60-per-min and the pictures are taken every 4-cycle
        (nx, nt) = upslope.shape
        t_arr = np.arange(time_step*upslope_delay, (nt+upslope_delay)*time_step, time_step)
        t_denoised = upslope.sum(axis = 0)/nx
        model = LinearRegression()
        model.fit(t_arr.reshape(-1, 1), t_denoised.reshape(-1, 1))
        dc_dt = model.coef_[0][0]
        t_plot = fullarray.sum(axis = 0)/nx
        nt1 = PixelValArray.shape[1]
        t_arr2 = np.arange(0, nt1*time_step, time_step)
        plt.scatter(t_arr2, t_plot, color = 'black')
        pred = model.predict(t_arr.reshape(-1, 1))
        plt.plot(t_arr.reshape(-1, 1), pred.reshape(-1, 1), color = 'red')
        #plt.text(6, 250, f"slope = {int(model.coef_[0][0]*100)/100}", rotation = 70)
        plt.title("Contrast over time")
        plt.xlabel("time (sec)")
        plt.ylabel("Contrast (HU)")
        print(f"The slope in the 1st fifth of the vessel length = {dc_dt}")
        
        upslope = PixelValArray[int(v_length/5):int(2*v_length/5), upslope_delay:upslope_peak] 
        fullarray = PixelValArray[int(v_length/5):int(2*v_length/5), :]
        #print(x_step)
        time_step = 3.2#4 #assuming the heart-rate to be 60-per-min and the pictures are taken every 4-cycle
        (nx, nt) = upslope.shape
        t_arr = np.arange(time_step*upslope_delay, (nt+upslope_delay)*time_step, time_step)
        t_denoised = upslope.sum(axis = 0)/nx
        model = LinearRegression()
        model.fit(t_arr.reshape(-1, 1), t_denoised.reshape(-1, 1))
        dc_dt = model.coef_[0][0]
        t_plot = fullarray.sum(axis = 0)/nx
        nt1 = PixelValArray.shape[1]
        t_arr2 = np.arange(0, nt1*time_step, time_step)
        plt.scatter(t_arr2, t_plot, color = 'blue')
        pred = model.predict(t_arr.reshape(-1, 1))
        plt.plot(t_arr.reshape(-1, 1), pred.reshape(-1, 1), color = 'red')
        #plt.text(6, 250, f"slope = {int(model.coef_[0][0]*100)/100}", rotation = 70)
        print(f"The slope in the 2nd fifth of the vessel length = {dc_dt}")
        
        upslope = PixelValArray[int(2*v_length/5):int(3*v_length/5), upslope_delay:upslope_peak] 
        fullarray = PixelValArray[int(2*v_length/5):int(3*v_length/5), :]
        #print(x_step)
        time_step = 3.2#4 #assuming the heart-rate to be 60-per-min and the pictures are taken every 4-cycle
        (nx, nt) = upslope.shape
        t_arr = np.arange(time_step*upslope_delay, (nt+upslope_delay)*time_step, time_step)
        t_denoised = upslope.sum(axis = 0)/nx
        model = LinearRegression()
        model.fit(t_arr.reshape(-1, 1), t_denoised.reshape(-1, 1))
        dc_dt = model.coef_[0][0]
        t_plot = fullarray.sum(axis = 0)/nx
        nt1 = PixelValArray.shape[1]
        t_arr2 = np.arange(0, nt1*time_step, time_step)
        plt.scatter(t_arr2, t_plot, color = 'green')
        pred = model.predict(t_arr.reshape(-1, 1))
        plt.plot(t_arr.reshape(-1, 1), pred.reshape(-1, 1), color = 'red')
        #plt.text(6, 250, f"slope = {int(model.coef_[0][0]*100)/100}", rotation = 70)
        print(f"The slope in the 3rd fifth of the vessel length = {dc_dt}")
        
        upslope = PixelValArray[int(3*v_length/5):int(4*v_length/5), upslope_delay:upslope_peak] 
        fullarray = PixelValArray[int(3*v_length/5):int(4*v_length/5), :]
        #print(x_step)
        time_step = 3.2#4 #assuming the heart-rate to be 60-per-min and the pictures are taken every 4-cycle
        (nx, nt) = upslope.shape
        t_arr = np.arange(time_step*upslope_delay, (nt+upslope_delay)*time_step, time_step)
        t_denoised = upslope.sum(axis = 0)/nx
        model = LinearRegression()
        model.fit(t_arr.reshape(-1, 1), t_denoised.reshape(-1, 1))
        dc_dt = model.coef_[0][0]
        t_plot = fullarray.sum(axis = 0)/nx
        nt1 = PixelValArray.shape[1]
        t_arr2 = np.arange(0, nt1*time_step, time_step)
        plt.scatter(t_arr2, t_plot, color = 'purple')
        pred = model.predict(t_arr.reshape(-1, 1))
        plt.plot(t_arr.reshape(-1, 1), pred.reshape(-1, 1), color = 'red')
        #plt.text(6, 250, f"slope = {int(model.coef_[0][0]*100)/100}", rotation = 70)
        print(f"The slope in the 4th fifth of the vessel length = {dc_dt}")
        
        upslope = PixelValArray[int(4*v_length/5):v_length-3, upslope_delay:upslope_peak] 
        fullarray = PixelValArray[int(4*v_length/5):v_length-3, :]
        #print(x_step)
        time_step = 3.2#4 #assuming the heart-rate to be 60-per-min and the pictures are taken every 4-cycle
        (nx, nt) = upslope.shape
        t_arr = np.arange(time_step*upslope_delay, (nt+upslope_delay)*time_step, time_step)
        t_denoised = upslope.sum(axis = 0)/nx
        model = LinearRegression()
        model.fit(t_arr.reshape(-1, 1), t_denoised.reshape(-1, 1))
        dc_dt = model.coef_[0][0]
        t_plot = fullarray.sum(axis = 0)/nx
        nt1 = PixelValArray.shape[1]
        t_arr2 = np.arange(0, nt1*time_step, time_step)
        plt.scatter(t_arr2, t_plot, color = 'cyan')
        pred = model.predict(t_arr.reshape(-1, 1))
        plt.plot(t_arr.reshape(-1, 1), pred.reshape(-1, 1), color = 'red')
        #plt.text(6, 250, f"slope = {int(model.coef_[0][0]*100)/100}", rotation = 70)
        print(f"The slope in the last fifth of the vessel length = {dc_dt}")
        
        plt.show()
        ##########
        
        
if __name__ == "__main__":
    #Description
    parser = argparse.ArgumentParser(description = "This script gets the centerline files with the averaged pixel value array and calculates the advection diffusion velocity")
    
    #InputFiles
    parser.add_argument('-InputFolder', '--InputFolder', type = str, required = True, dest = "InputFolder", help = "The Folder containing the centerline files with averaged pixel values")
    
    args = parser.parse_args()
    ImageAnalysisAdvectionDiffusionAlongCL(args).Main()
    