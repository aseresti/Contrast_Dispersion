""" This script implements the Contrast Dispersion in the Descending Aorta, which is segmented out of the CT-MPI images of the patient. 
The script takes a folder containing the volume files onto which the pixel values are projected, and the surface file of the aorta. 
It extracts the centerline along the aorta and calculates the average pixel value inside a sphere clip that moves along the centerline. 
It interpolates points along the centerline in time to compensate for the dynamic shuttle mode in 64-row CT scanners. 
It calculates dC/dt using the average value of the first clip over time and dC/dx using the average pixel value along the centerline 
(one half taken from the peak, and the other half taken from the pre-peak). The user needs to specify the patient's heartbeats, 
the delay (or the timestep at which the upslope begins), and the peak (or the timestep at which the upslope reaches the peak). 
For future use, the script might require improvements.

Author: ana @ github.com/aseresti
Date: April 2024
"""
import os
import argparse

import numpy as np
import vtk
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple
#from itertools import zip_longest

from utilities import ReadVTUFile, ReadVTPFile, WriteVTUFile, GetCentroid, vtk_to_numpy, ThresholdInBetween
from ContrastTools import Model1D, MAFilter


class ContrastDispersionAlongVessel():
    def __init__(self, args: argparse) -> None:
        """the initial operations  on input data

        Args:
            args (argparse): Must contain the following instances: an InputFolderName, HeartBeat, delay, peak
        """
        self.Args = args
        
        # takes the input folder name to read the Input Volumes, Input Surface
        filenames: str = os.listdir(self.Args.InputFolderName)
        Volumefilenames: List[str] = [filename for filename in filenames if "vtu" in filename]
        self.VolumeFileNames: List[str] = sorted(Volumefilenames)

        SurfaceFileName: List[str] = [filename for filename in filenames if "surface" in filename]
        self.SurfaceFileName: str = SurfaceFileName[0]

        slicenames: List[str] = [filename for filename in filenames if "slice" in filename]
        self.SliceNames: List[str] = sorted(slicenames)

        # creating the output folder to store the script output files
        self.OutputFolderName: str = "ClippingOutput"
        if self.OutputFolderName not in filenames:
            print(f"*** Making A New Directory: {self.OutputFolderName}")
            os.system(f"mkdir {self.Args.InputFolderName}/{self.OutputFolderName}")
        else:
            print(f"*** {self.OutputFolderName} Already Exists!")

        self.TemporalInterval = 4 * 60 / self.Args.HeartBeat # the temporal interval between each timestep
        
        # Adjustable Parameters (depending on the case, the user might need to change these parameters)
        self.MAFilter_Length = 0 # should you prefer to smooth the data along the centerline, you can use the MAFilter or Moving average filter
        self.radius = 3 # the radius of the sphere clip 
        # to-do: set the radius of the sphere clip as a percentage of the redius of the descending aorta
        self.CL_Point_Length = 4 # the distance between the two CL points in mm (Try to avoid overlapping sphere clips)
        self.interpolation_factor = 2  # To compensate the dynamic shuttle mode in CT scanners, at least you need to interpolate it with a factor of 2. 4 is recommended.
        self.interpolation_peak = self.interpolation_factor * self.Args.peak - 2 # the index of the peak point of the time attenuation curve after interpolation
        self.interpolation_pre_peak = self.interpolation_peak - int(self.interpolation_factor/2) # the index of the pre peak point of the time attenuation curve after interpolation

    def ExtractCeneterLine(self) -> None:
        """Extracts the centerline outof the surface of the vessel. Reads and stores the centerline and the number of points along it as features of the instance of the class
        """

        CL_File_Name = "aorta_cl.vtp"
        # to first extract the centerline, then smooth it and finally resample it using the provided length between each CL point
        os.system(f"vmtkcenterlines -ifile {self.Args.InputFolderName}/{self.SurfaceFileName} -ofile {self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name} -endpoints 1 -resampling 1 -resamplingstep 4")
        os.system(f"vmtkcenterlinesmoothing -ifile {self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name} -ofile {self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name} -factor 0.1 -iterations 200")
        os.system(f"vmtkcenterlineresampling -ifile {self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name} -ofile {self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name} -length {self.CL_Point_Length}")
        
        # save the Centerline for further visualizations
        self.CLFile = ReadVTPFile(f"{self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name}")
        self.NPoints = self.CLFile.GetNumberOfPoints() # the number of points along the centerline; specifies the number of samples taken along the vessel for C(x)

    def ExtractPieceWiseCeneterLine(self) -> np.array:
        center = []
        count = 0
        for SliceName in self.SliceNames:
            slice = ReadVTPFile(f"{self.Args.InputFolderName}/{SliceName}")
            print(GetCentroid(slice))
            center.append(GetCentroid(slice))
            count += 1

        CLPoints = []
        for i in range(1,count-1):
            x0 = center[i-1][0]; y0 = center[i-1][1]; z0 = center[i-1][2]
            x1 = center[i][0]; y1 = center[i][1]; z1 = center[i][2]
            for X in np.arange(x0,x1,0.5):
                Y = (y1-y0)/(x1-x0)*(X-x0) + y0
                Z = (z1-z0)/(x1-x0)*(X-x0) + z0
                CLPoints.append([X,Y,Z])
        
        return CLPoints
    
    def CenterLineMeshSection(self,volume: vtk.vtkUnstructuredGrid) -> np.array:
        """ Extracts the surface mesh sections along the centerline and reads the average pixel values.

        Args:
            volume (vtk.vtkUnstructuredGrid): the vtu file onto which the pixel values are projected
        Returns:
            Contrast_Along_CL (np.array): An array including the average pixel values along the centerline of the vessel.
        """

        CL_File_Name = "aorta_cl.vtp"
        MS_File_Name = "aorta_mesh_section.vtp"
        os.system(f"vmtkcenterlinemeshsections -centerlinesfile {self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name} -ifile {self.Args.InputFolderName}/{self.OutputFolderName}/{volume} -ofile {self.Args.InputFolderName}/{self.OutputFolderName}/{MS_File_Name}")
        
        SurfaceSections = ReadVTPFile(f"{self.Args.InputFolderName}/{self.OutputFolderName}/{MS_File_Name}")
        Nstart,Nend=SurfaceSections.GetPointData().GetArray("SectionIds").GetRange()
        Nstart=int(Nstart)
        Nend=int(Nend)
        Contrast_Along_CL = np.empty([Nend-Nstart,1])
        
        for i in range(Nstart,Nend):
			#Get the Section of the Slice
            section_ = ThresholdInBetween(SurfaceSections,"SectionIds",i,i)
            ContrastArray = vtk_to_numpy(section_.GetPointData().GetArray("scalars"))
            Contrast_Along_CL[i] = np.average(ContrastArray)
        
        return Contrast_Along_CL
    


    def SphereClip(self, center: vtk.vtkPoints, volume: vtk.vtkUnstructuredGrid) -> Tuple[vtk.vtkUnstructuredGrid, float]:
        """ clips the sphere out of the given volume on the given center point 

        Args:
            center (vtk.vtkPoints): (x, y, z) of the point from points along the centerline
            volume (vtk.vtkUnstructuredGrid): the vtu file onto which the pixel values are projected

        Returns:
            clipper.GetOutput() (vtk.vtkUnstructuredGrid): The clipped sphere containing the pixel value data and other specifics of the given volume
            AveragePixelValue (float): The average pixel value inside the clipped sphere
        """

        #define the Sphere
        Sphere = vtk.vtkSphere()
        Sphere.SetCenter(center)
        Sphere.SetRadius(self.radius)

        #Implement vtkclipping filter with the defined "sphere" as the clipping function
        clipper = vtk.vtkClipDataSet()
        clipper.SetInputData(volume)
        clipper.SetClipFunction(Sphere)
        clipper.InsideOutOn()
        clipper.GetOutputInformation(1)
        clipper.Update()

        #WriteVTUFile(f"{self.Args.InputFolderName}/{self.OutputFileName}/clipper.vtu", clipper.GetOutput())
        
        # calculating the average pixel value iside each Sphere
        ArrayName = clipper.GetOutput().GetPointData().GetArrayName(0)
        SphereOutput = clipper.GetOutput().GetPointData().GetArray(ArrayName)
        SphereOutput = vtk_to_numpy(SphereOutput)
        AveragePixelValue = np.average(SphereOutput)

        return clipper.GetOutput(), AveragePixelValue

    def AppendVolumes(self,volume1: vtk.vtkUnstructuredGrid, volume2: vtk.vtkUnstructuredGrid) -> vtk.vtkUnstructuredGrid:
        """Appends two given volumes, it is used to have all of the sphere clips stored in one single file

        Args:
            volume1 (vtk.vtkUnstructuredGrid): The clipped volume
            volume2 (vtk.vtkUnstructuredGrid): The clipped volume

        Returns:
            vtk.vtkUnstructuredGrid: The appended volume
        """
        append_filter = vtk.vtkAppendFilter()
        append_filter.AddInputData(volume1) 
        append_filter.AddInputData(volume2)
        append_filter.Update()

        return append_filter.GetOutput()
    
    def ContrastDisperssion(self, Inflow_Contrast_Temporal: np.array, CenterLineContrastDict: Dict[str,np.array]) -> None:
        """ Implements the contrast Dispersion method taking the temporal and spatial data

        Args:
            Inflow_Contrast_Temporal (np.array): The average pixel values along in the center of the Inflow of the vessel
            CenterLineContrastDict (Dict[np.array]): A dictionary containing the average pixel value along the centerline of the vessel in different time points
        """

        [x, t, UpSlope] = self.CreateCoords()
        Inflow_Upslope: np.array[float] = Inflow_Contrast_Temporal[self.Args.delay: self.Args.peak] # the images that are taken during the upslope time

        # find the velocity before interpolation (just in case, if the vessel doesn't have a jump in pixel values need interpolation)
        CenterLine_Contrast: np.array = CenterLineContrastDict[self.VolumeFileNames[self.Args.peak]] # taking the centerline extracted from the image at the peak of the upslope
        [xslope, xpred] = Model1D(x,CenterLine_Contrast) # calculating dC/dx
        [tslope, tpred] = Model1D(UpSlope, Inflow_Upslope) # calculating dC/dt
        VelocityPredicted: float = tslope/xslope # calculate the velocity of Contrast Dispersion
        print("the predicted velocity is: ", abs(VelocityPredicted), " mm/s")

        #Plotting Curves C(x) and C(t)
        plt.figure(figsize=(13,5))

        # plot the contrast along the time
        plt.subplot(121)
        plt.scatter(t.reshape(-1,1), Inflow_Contrast_Temporal.reshape(-1,1), label = "Inflow Center")
        plt.plot(UpSlope.reshape(-1,1), tpred.reshape(-1,1), color = 'red', label = 'Linear (Inflow Cenetr)')
        plt.title("Temporal Attenuation Curve")
        plt.xlabel("time (s)")
        plt.legend()

        # plot the contrast along the centerline which was wxtracted out of the image at the peak up-slope
        plt.subplot(122)
        plt.scatter(x.reshape(-1,1), CenterLine_Contrast.reshape(-1,1), label = "Lumen Centerline")
        plt.plot(x.reshape(-1,1), xpred.reshape(-1,1), color = 'red', label='Linear (Linear Center)')
        plt.title("Spatial Attenuation Curve")
        plt.xlabel("length (mm)")
        plt.legend()
        
        plt.show()

        # to-do: add after interpolation spatial data
        # Write the text output file to store the contrast in time and along the centerline
        TextFileName = "Output.txt"
        TextFile_data = [
            ["x_array (mm)"] + list(item[0] for item in x),
            ["Contrast_Along_CL (HU)"] + list(CenterLine_Contrast),
            ["t_array_UpSlope (s)"] + list(UpSlope),
            ["Inflow Contrast (HU)"] + list(item[0] for item in Inflow_Upslope)
        ]

        # Write the transposed data to the output file
        with open(f"{self.Args.InputFolderName}/{self.OutputFolderName}/{TextFileName}", "w") as writefile:
            for row in TextFile_data:
                writefile.write(", ".join(str(item) for item in row) + "\n")

        # extracting the centerline after the interpolation; half of the points are extracted from the peak and the other half from the pre-peak and being concatenated
        interp_CenterLineContrastDict: Dict[str,np.array] = self.TemporalInterpolation(CenterLineContrastDict,t)
        dict_item_interp = list(interp_CenterLineContrastDict.items())
        CenterLine_Contrast0: np.array = dict_item_interp[self.interpolation_pre_peak][1]
        CenterLine_Contrast1: np.array = dict_item_interp[self.interpolation_peak][1]
        CenterLine_Contrast2 = np.concatenate((CenterLine_Contrast1[:int(self.NPoints/2)], CenterLine_Contrast0[int(self.NPoints/2):]))
        
        [xslope, xpred] = Model1D(x,CenterLine_Contrast2) # calculating dC/dx from the concatenated centerline points
        VelocityPredicted: float = tslope/xslope # calculating the velocity after interpolation

        # plot the peak and pre-peak centerlines and the concatenated centerline
        plt.figure(figsize=(10,5))
        plt.plot(x.reshape(-1,1), CenterLine_Contrast2.reshape(-1,1), color = 'orange', label = 'concatenated signal')
        plt.scatter(x.reshape(-1,1), CenterLine_Contrast1.reshape(-1,1), color= 'green', label = 'peak centerline', marker='^')
        plt.scatter(x.reshape(-1,1), CenterLine_Contrast0.reshape(-1,1), color= 'blue', label = 'pre-peak centerline', marker='^')
        plt.plot(x.reshape(-1,1), xpred.reshape(-1,1), color = 'red', label='Linear (Linear Center)')
        plt.title('The contrast along the centerline of the vessel')
        plt.xlabel('length (mm)')
        plt.legend()
        plt.show()
        
        print("the predicted velocity after interpolation is: ", abs(VelocityPredicted), " mm/s")


    def TemporalInterpolation(self,CenterLineContrastDict: Dict[str,np.array], t: np.array) -> Dict[str,np.array]:
        """Takes the average pixel values along the centerline in different time points and interpolates them in time 

        Args:
            CenterLineContrastDict (Dict[np.array]): A dictionary containing the average pixel value along the centerline of the vessel in several time steps
            t (np.array): the time domain coordinate

        Returns:
            Dict[np.array: A dictionary containing the average pixel value along the centerline of the vessel in interpolated time steps
        """

        dict_item: list = sorted(CenterLineContrastDict.items())

        indeces: np.array = np.arange(0,len(t),1)
        indeces_times2: np.array = np.arange(0,len(t), 1./self.interpolation_factor)

        new_t: np.array = np.interp(indeces_times2, indeces, t)
        New_CenterLineContrastDict: Dict[str,np.array] = {f'interp_{i}': np.empty([self.NPoints,1]) for i in range(0,len(indeces_times2))} #interpolated data

        for point in range(self.NPoints):
            Points = np.empty([len(t)])
            count = 0

            for CenterLineContrast in dict_item:
                Points[count] = CenterLineContrast[1][point]
                count +=1
            

            interpolateSignal = np.interp(new_t, t, Points)
            
            for i, (key, value) in enumerate(New_CenterLineContrastDict.items()):
                value[point] = interpolateSignal[i]
        
        # plotting the time attenuation curve of the initial signal and the interpolated signal; specifying the peak and pre-peak time-points; The plot is based on the last point on the centerline
        plt.figure(figsize=(13,7))
        plt.scatter(t.reshape(-1,1), Points.reshape(-1,1), marker='o', label='Before Interpolation')
        plt.scatter(new_t.reshape(-1,1), interpolateSignal.reshape(-1,1), marker='*', label= 'After Interpolation')
        plt.scatter(new_t[self.interpolation_peak], interpolateSignal[self.interpolation_peak], marker='+', label='peak point')
        plt.scatter(new_t[self.interpolation_pre_peak], interpolateSignal[self.interpolation_pre_peak], marker='+', label='pre peak point')
        plt.title('Temporal Attenuation Curve')
        plt.xlabel('time (s)')
        plt.legend()
        plt.show()
        
        return New_CenterLineContrastDict

    def CreateCoords(self) -> Tuple[np.array, np.array, np.array]:
        """ Generates the spatial and temporal coordinates

        Returns:
            np.array: spatial and temporal coordinates
        """
        # Create the spatial coordinate (x) along the centerline of the vessel
        point0 = self.CLFile.GetPoint(0)
        point = self.CLFile.GetPoint(self.NPoints-1)
        Length = (
                (point[0]-point0[0])**2 +
                (point[1]-point0[1])**2 +
                (point[2]-point0[2])**2)**0.5
        print(f"***** The length of the centerline is: {Length} mm")
        CL_Coord: np.array = np.empty([self.NPoints - self.MAFilter_Length,1])
        CL_Coord[0] = 0
        for n in range(1,self.NPoints - self.MAFilter_Length):
            #point0 = self.CLFile.GetPoint(n-1)
            point = self.CLFile.GetPoint(n)
            CL_Coord[n] = (
                (point[0]-point0[0])**2 +
                (point[1]-point0[1])**2 +
                (point[2]-point0[2])**2)**0.5
        
        # Create the time coordinate based on the heartbeat; each images in CT-MPI are by 4 heartbeats apart
        TimePoints: int = self.VolumeFileNames.__len__()
        Time_Coord: np.array = np.array([i*self.TemporalInterval for i in range(TimePoints)])

        # seperating the upslope part of the coordinate
        UpSlope: np.array = Time_Coord[self.Args.delay:self.Args.peak]

        return CL_Coord, Time_Coord, UpSlope

    def ReadInputVolumes(self) -> Dict[str,vtk.vtkUnstructuredGrid]:
        """ Takes the input volume names to read it as a vtu file and returns the vtk.vtkUnstructuredGrid

        Returns:
            Dict[vtk.vtkUnstructuredGrid]: A dictionary containing the input volumes
        """
        InputVolumesDict = dict()

        for volume_ in self.VolumeFileNames:
            volume = ReadVTUFile(f"{self.Args.InputFolderName}/{volume_}")
            if volume.GetPointData().GetArray("scalars"):
                InputVolumesDict[f"{volume_}"] = volume
            else:
                ArrayName = volume.GetPointData().GetArrayName(0)
                volume.GetPointData().GetArray(ArrayName).SetName("scalars")
                InputVolumesDict[f"{volume_}"] = volume
        
        return sorted(InputVolumesDict.items())

    def main(self) -> None:
        """ The core of the class to handle the process on sequence
        """
        print("-"*10, "***" "-"*10)
        self.ExtractCeneterLine()
        print("-"*10, "***" "-"*10)
        print("--- number of point along the centerline = ", self.NPoints)

        print("-"*10, "***" "-"*10)
        print("--- Extracting the temporal attenuation curve at the inflow of the vessel")
        InputVolumesDict: Dict[str,vtk.vtkUnstructuredGrid] = self.ReadInputVolumes()
        Inflow_Contrast_Temporal: np.array[float] = np.empty([self.VolumeFileNames.__len__(),1])

        count = 0
        for volume in InputVolumesDict:
            print(f"--- Reading the center of the Inflow of {volume[0]}")
            Inflow_Clip = self.SphereClip(self.CLFile.GetPoint(0), volume[1])
            Inflow_Contrast_Temporal[count] = Inflow_Clip[1]
            count += 1
        
        InflowClipFileName: str = 'InflowClip.vtu'
        print(f"--- Storing the inflow clipper in a vtu file at {self.Args.InputFolderName}/{self.OutputFolderName}/{InflowClipFileName}")
        WriteVTUFile(f"{self.Args.InputFolderName}/{self.OutputFolderName}/{InflowClipFileName}", Inflow_Clip[0])
        
        # Reading and Storing the centerline of the every files in the upslope
        CenterLineContrastDict: Dict[str,np.array] = {volume[0]: None for volume in InputVolumesDict[:]}
        print("--- Extracting the average contrast along the centerline of the vessel")
        count = 1
        for volume in InputVolumesDict:
            print(f"--- Reading the centerline of the volume: {volume[0]}")
            Clip1: vtk.vtkUnstructuredGrid = None # Clear Clip1 before reassignment
            CenterLineContrast = np.empty([self.NPoints,1])
            [Clip1, CenterLineContrast[0]] = self.SphereClip(self.CLFile.GetPoint(0), volume[1])
            for point in range(1,self.NPoints):
                [Clip2, CenterLineContrast[point]] = self.SphereClip(self.CLFile.GetPoint(point), volume[1])
                Clip1 = self.AppendVolumes(Clip1, Clip2)
            
            count+=1
            if count == self.Args.peak:
                PeakCylinderClip = Clip1
            
            CenterLineContrastDict[volume[0]] = CenterLineContrast
        
        #CenterLineContrastDict = dict(sorted(CenterLineContrastDict.items()), key=lambda x: x[0])
        ClipOutputFile = "ClippedVolume.vtu"
        print(f"--- Storing the centerline clipper at peak in a vtu file at {self.Args.InputFolderName}/{self.OutputFolderName}/{ClipOutputFile}")
        WriteVTUFile(f"{self.Args.InputFolderName}/{self.OutputFolderName}/{ClipOutputFile}", PeakCylinderClip)

        print("-"*10, "***" "-"*10)
        print("--- Implementing the contrast dispersion and the predicted velocity")
        self.ContrastDisperssion(Inflow_Contrast_Temporal, CenterLineContrastDict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputFolderName", "--InputFolderName", type=str, required=True, dest="InputFolderName")
    parser.add_argument("-HeartBeat", "--HeartBeat", type=float, required=True, dest="HeartBeat")
    parser.add_argument("-delay", "--delay", type=int, required=True, dest="delay", help="bolus time delay before upslope rises")
    parser.add_argument("-peak", "--peak", type=int, required=True, dest="peak", help="the number of image on witch the upslope reaches its peak")

    args = parser.parse_args()

    ContrastDispersionAlongVessel(args).main()


