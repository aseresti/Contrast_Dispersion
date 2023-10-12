#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:04:22 2023

@author: ana
"""

import glob
import os
import argparse


class ConvertDicomtoVTI():
    def __init__(self,Args):
        self.Args = Args
    def main(self):
        filenames = glob.glob(f'{self.Args.InputFolderName}/*.dcm')
        filenames = sorted(filenames)
        self.N_file_per_cycle = int(len(filenames)/self.Args.NumberOfCycles)
        

        #self.N_img_per_cycle = int(len(filenames)/self.Args.NumberOfCycles)
        for i in range(0,self.Args.NumberOfCycles):
            directoryname = f'perfusion_image_cycle_{i+1}'
            pathDicom = f'{self.Args.InputFolderName}/{directoryname}'
            os.system(f"mkdir {pathDicom}")
            for j in range((i)*self.N_file_per_cycle,(i+1)*self.N_file_per_cycle-1):
                os.system(f'cp {filenames[j]} {pathDicom}')
            
            print(f'--- Looping over cycle: {i}')
            filenames_ = glob.glob(f'{pathDicom}/*.dcm')
            filenames_ = sorted(filenames_)
            os.system(f'vmtkimagereader -ifile {filenames_[0]} --pipe vmtkimagewriter -ofile {self.Args.InputFolderName}/Cycle_Image_{i+1}.vti')

if __name__ == '__main__':
    #descreption
    parser = argparse.ArgumentParser(description='Thsi script takes a dicom folder with N cycles and outputs an averaged vti image')
    #Input
    parser.add_argument('-InputFolder', '--InputFolder', type = str, required = True, dest = 'InputFolderName', help = 'The name of the folder with all of the dicom files')
    #NumberOfCycles
    parser.add_argument('-NofCycle', '--NumberOfCycles', type = int, required = True, dest = 'NumberOfCycles', help = 'The number of perfusion images that are in the dicom folder')
    args = parser.parse_args()
    ConvertDicomtoVTI(args).main()
