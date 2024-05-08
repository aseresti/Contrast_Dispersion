"""This script is developed to validate the interpolation function implemented
on CT-MPI data.
"""

import os
import argparse
import sys

import numpy as np
import vtk
import matplotlib.pyplot as plt

#from utilities import ReadVTUFile, ReadVTPFile, WriteVTUFile, GetCentroid, vtk_to_numpy, ThresholdInBetween
path = os.path.abspath("/Users/ana/Documents/AnahitaSeresti/02_ContrastStudy/Contrast_Dispersion/scripts/CT_MPI_Tools")
sys.path.append(path)

path = os.path.abspath("/Users/ana/Documents/AnahitaSeresti/02_ContrastStudy/Contrast_Dispersion/scripts")
sys.path.append(path)

from CT_MPI_Tools.ContrastDispersionAlongVessel import ContrastDispersionAlongVessel

def ImplementPipeline(data):

    data.interpolation_factor = 2
    data.interpolation_peak = data.interpolation_factor * data.Args.peak - 2 # the index of the peak point of the time attenuation curve after interpolation
    data.interpolation_pre_peak = data.interpolation_peak - int(data.interpolation_factor/2) # the index of the pre peak point of the time attenuation curve after interpolation

    print(data.interpolation_peak)

    InputVolumesDict = data.ReadInputVolumes()

    CenterLineContrastDict = {volume[0]: None for volume in InputVolumesDict[:]}
    print("--- Extracting the average contrast along the centerline of the vessel")
    count = 1
    for volume in InputVolumesDict:
        print(f"--- Reading the centerline of the volume: {volume[0]}")
        Clip1: vtk.vtkUnstructuredGrid = None # Clear Clip1 before reassignment
        CenterLineContrast = np.empty([data.NPoints,1])
        [Clip1, CenterLineContrast[0]] = data.SphereClip(data.CLFile.GetPoint(0), volume[1])
        for point in range(1,data.NPoints):
            [Clip2, CenterLineContrast[point]] = data.SphereClip(data.CLFile.GetPoint(point), volume[1])
            Clip1 = data.AppendVolumes(Clip1, Clip2)
            
        count+=1
            
        CenterLineContrastDict[volume[0]] = CenterLineContrast
        
    #CenterLineContrastDict = dict(sorted(CenterLineContrastDict.items()), key=lambda x: x[0])
    ClipOutputFile = "ClippedVolume.vtu"
    print(f"--- Storing the centerline clipper at peak in a vtu file at {data.Args.InputFolderName}/{data.OutputFolderName}/{ClipOutputFile}")
    #WriteVTUFile(f"{data.Args.InputFolderName}/{data.OutputFolderName}/{ClipOutputFile}", PeakCylinderClip)


    return CenterLineContrastDict

def ExtractTemporalData(t,TheDict):
    TemporalData = np.empty(len(t))
    for i, (key, value) in enumerate(TheDict.items()):
        TemporalData[i] = value[0]
    return TemporalData



if __name__ == "__main__":

    parser1 = argparse.ArgumentParser()
    parser1.add_argument("-InputFolderName", "--InputFolderName", type=str, required=False, default="/Users/ana/Documents/AnahitaSeresti/02_ContrastStudy/Image_Based_Volumes/Post-CABG/CABG10B/step6_InterpolationValidation",  dest="InputFolderName")
    parser1.add_argument("-HeartBeat", "--HeartBeat", type=float, required=False, dest="HeartBeat", default=37.5)
    parser1.add_argument("-delay", "--delay", type=int, required=False, dest="delay", default=1, help="bolus time delay before upslope rises")
    parser1.add_argument("-peak", "--peak", type=int, required=False, dest="peak",default=3, help="the number of image on witch the upslope reaches its peak")

    args1 = parser1.parse_args()

    parser2 = argparse.ArgumentParser()
    parser2.add_argument("-InputFolderName", "--InputFolderName", type=str, required=False, default="/Users/ana/Documents/AnahitaSeresti/02_ContrastStudy/Image_Based_Volumes/Post-CABG/CABG10B/Step5_Projected_Volume_descending_aorta",  dest="InputFolderName")
    parser2.add_argument("-HeartBeat", "--HeartBeat", type=float, required=False, dest="HeartBeat", default=75)
    parser2.add_argument("-delay", "--delay", type=int, required=False, dest="delay", default=3, help="bolus time delay before upslope rises")
    parser2.add_argument("-peak", "--peak", type=int, required=False, dest="peak",default=7, help="the number of image on witch the upslope reaches its peak")

    args2 = parser2.parse_args()

    data_complete = ContrastDispersionAlongVessel(args2)
    data_sparse = ContrastDispersionAlongVessel(args1)
    data_complete.ExtractCeneterLine()
    data_sparse.CLFile = data_complete.CLFile
    data_sparse.NPoints = data_complete.NPoints

    CenterLineContrastDict_complete = ImplementPipeline(data_complete)
    CenterLineContrastDict_sparse = ImplementPipeline(data_sparse)

    (x1,t1,up_slope1) = data_complete.CreateCoords()
    (x2,t2,up_slope2) = data_sparse.CreateCoords()
    
    interpolate_CenterlineDict = data_sparse.TemporalInterpolation(CenterLineContrastDict_sparse, t2)

    indeces: np.array = np.arange(0,len(t2),1)
    indeces_times2: np.array = np.arange(0,len(t2), 1./data_sparse.interpolation_factor)

    new_t: np.array = np.interp(indeces_times2, indeces, t2)

    SparseSignal = ExtractTemporalData(t2, CenterLineContrastDict_sparse)
    interpolateSignal = ExtractTemporalData(new_t,interpolate_CenterlineDict)
    ActualSignal = ExtractTemporalData(t1, CenterLineContrastDict_complete)
    
    plt.figure(figsize=(13,7))
    plt.scatter(t1.reshape(-1,1), ActualSignal.reshape(-1,1), marker='o', label='Actual Signal')
    plt.scatter(t2.reshape(-1,1), SparseSignal.reshape(-1,1), marker = '*', label = 'Sparse Signal Before Interpolation')
    plt.scatter(new_t.reshape(-1,1), interpolateSignal.reshape(-1,1), marker='+', label= 'Sparse Signal After Interpolation')
    plt.title('Temporal Attenuation Curve')
    plt.xlabel('time (s)')
    plt.legend()
    plt.show()


