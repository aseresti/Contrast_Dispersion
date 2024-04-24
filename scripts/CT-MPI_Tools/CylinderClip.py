
import os
import argparse

import numpy as np
import vtk
import matplotlib.pyplot as plt
#from itertools import zip_longest

from utilities import ReadVTUFile, ReadVTPFile, WriteVTUFile, GetCentroid, vtk_to_numpy, ThresholdInBetween
from ContrastTools import Model1D, MAFilter


class CylinderClipAlongCL():
    def __init__(self, args) -> None:
        self.Args = args
        
        filenames = os.listdir(self.Args.InputFolderName)
        Volumefilenames = [filename for filename in filenames if "vtu" in filename]
        self.VolumeFileNames = sorted(Volumefilenames)

        SurfaceFileName = [filename for filename in filenames if "surface" in filename]
        self.SurfaceFileName = SurfaceFileName[0]

        slicenames = [filename for filename in filenames if "slice" in filename]
        self.SliceNames = sorted(slicenames)

        self.OutputFolderName = "ClippingOutput"
        if self.OutputFolderName not in filenames:
            print(f"*** Making A New Directory: {self.OutputFolderName}")
            os.system(f"mkdir {self.Args.InputFolderName}/{self.OutputFolderName}")
        else:
            print(f"*** {self.OutputFolderName} Already Exists!")

        self.TemporalInterval = 60 / self.Args.HeartBeat
        
        # Adjustable Parameters
        self.MAFilter_Length = 0
        self.radius = 7 # the radius of the sphere clip
        self.CL_Point_Length = 4 # the distance between the two CL points
        self.interpolation_factor = 4 
        self.interpolation_peak = self.interpolation_factor * self.Args.peak - 2 # the index of the peak point of the time attenuation curve after interpolation
        self.interpolation_pre_peak = self.interpolation_peak - int(self.interpolation_factor/2) # the index of the pre peak point of the time attenuation curve after interpolation

    def ExtractCeneterLine(self):
        CL_File_Name = "aorta_cl.vtp"
        os.system(f"vmtkcenterlines -ifile {self.Args.InputFolderName}/{self.SurfaceFileName} -ofile {self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name} -endpoints 1 -resampling 1 -resamplingstep 4")
        os.system(f"vmtkcenterlinesmoothing -ifile {self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name} -ofile {self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name} -factor 0.1 -iterations 200")
        os.system(f"vmtkcenterlineresampling -ifile {self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name} -ofile {self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name} -length {self.CL_Point_Length}")
        
        self.CLFile = ReadVTPFile(f"{self.Args.InputFolderName}/{self.OutputFolderName}/{CL_File_Name}")
        self.NPoints = self.CLFile.GetNumberOfPoints()

    def ExtractPieceWiseCeneterLine(self):
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
    
    def CenterLineMeshSection(self,volume):
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
    


    def SphereClip(self, center, volume):
        #define the Sphere
        Sphere = vtk.vtkSphere()
        Sphere.SetCenter(center)
        Sphere.SetRadius(self.radius)

        #Implement vtkclipping filter "sphere"
        clipper = vtk.vtkClipDataSet()
        clipper.SetInputData(volume)
        clipper.SetClipFunction(Sphere)
        clipper.InsideOutOn()
        clipper.GetOutputInformation(1)
        clipper.Update()

        #WriteVTUFile(f"{self.Args.InputFolderName}/{self.OutputFileName}/clipper.vtu", clipper.GetOutput())
        
        #AverageFilter on each Sphere
        ArrayName = clipper.GetOutput().GetPointData().GetArrayName(0)
        SphereOutput = clipper.GetOutput().GetPointData().GetArray(ArrayName)
        SphereOutput = vtk_to_numpy(SphereOutput)
        averagePixelValue = np.average(SphereOutput)

        return clipper.GetOutput(), averagePixelValue

    def AppendVolumes(self,volume1,volume2):
        append_filter = vtk.vtkAppendFilter()
        append_filter.AddInputData(volume1)
        append_filter.AddInputData(volume2)
        append_filter.Update()
        return append_filter.GetOutput()
    
    def ContrastDisperssion(self, Inflow_Contrast_Temporal, CenterLineContrastDict):
        [x, t, UpSlope] = self.CreateCoords()
        Inflow_Upslope = Inflow_Contrast_Temporal[self.Args.delay: self.Args.peak]

        CenterLine_Contrast = CenterLineContrastDict[self.VolumeFileNames[self.Args.peak]]
        [xslope, xpred] = Model1D(x,CenterLine_Contrast)
        [tslope, tpred] = Model1D(UpSlope, Inflow_Upslope)
        VelocityPredicted = tslope/xslope
        print("the predicted velocity is: ", abs(VelocityPredicted), " mm/s")

        #Plotting Curves
        plt.figure(figsize=(13,5))

        plt.subplot(121)
        plt.scatter(t.reshape(-1,1), Inflow_Contrast_Temporal.reshape(-1,1), label = "Inflow Center")
        plt.plot(UpSlope.reshape(-1,1), tpred.reshape(-1,1), color = 'red', label = 'Linear (Inflow Cenetr)')
        plt.title("Temporal Attenuation Curve")
        plt.xlabel("time (s)")
        plt.legend()

        plt.subplot(122)
        plt.scatter(x.reshape(-1,1), CenterLine_Contrast.reshape(-1,1), label = "Lumen Centerline")
        plt.plot(x.reshape(-1,1), xpred.reshape(-1,1), color = 'red', label='Linear (Linear Center)')
        plt.title("Spatial Attenuation Curve")
        plt.xlabel("length (mm)")
        plt.legend()
        
        plt.show()

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

        interp_CenterLineContrastDict = self.UpslopeInterpolation(CenterLineContrastDict,t)
        dict_item_interp = list(interp_CenterLineContrastDict.items())
        CenterLine_Contrast0 = dict_item_interp[self.interpolation_pre_peak][1]
        CenterLine_Contrast1 = dict_item_interp[self.interpolation_peak][1]
        CenterLine_Contrast2 = np.concatenate((CenterLine_Contrast1[:int(self.NPoints/2)], CenterLine_Contrast0[int(self.NPoints/2):]))
        
        [xslope, xpred] = Model1D(x,CenterLine_Contrast2)
        VelocityPredicted = tslope/xslope

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


    def UpslopeInterpolation(self,CenterLineContrastDict, t):
        dict_item = sorted(CenterLineContrastDict.items())

        new_length = len(t) * self.interpolation_factor #interpolate the centerline in time domain because of the dynamic shuttle mode
        new_t = np.linspace(t[0], t[-1], new_length)

        New_CenterLineContrastDict = {f'interp_{i}': np.empty([self.NPoints,1]) for i in range(0,new_length)}

        for point in range(self.NPoints):
            Points = np.empty([len(t)])
            count = 0

            for CenterLineContrast in dict_item:
                Points[count] = CenterLineContrast[1][point]
                count +=1
            

            interpolateSignal = np.interp(new_t, t, Points)
            
            for i, (key, value) in enumerate(New_CenterLineContrastDict.items()):
                value[point] = interpolateSignal[i]
        
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


    def CreateCoords(self):
        
        point0 = self.CLFile.GetPoint(0)
        point = self.CLFile.GetPoint(self.NPoints-1)
        Length = (
                (point[0]-point0[0])**2 +
                (point[1]-point0[1])**2 +
                (point[2]-point0[2])**2)**0.5
        print(f"***** The length of the centerline is: {Length} mm")
        CL_Coord = np.empty([self.NPoints - self.MAFilter_Length,1])
        CL_Coord[0] = 0
        for n in range(1,self.NPoints - self.MAFilter_Length):
            #point0 = self.CLFile.GetPoint(n-1)
            point = self.CLFile.GetPoint(n)
            CL_Coord[n] = (
                (point[0]-point0[0])**2 +
                (point[1]-point0[1])**2 +
                (point[2]-point0[2])**2)**0.5
            
        TimePoints = self.VolumeFileNames.__len__()
        Time_Coord = np.array([i*self.TemporalInterval for i in range(TimePoints)])
        UpSlope = Time_Coord[self.Args.delay:self.Args.peak]

        return CL_Coord, Time_Coord, UpSlope

    def ReadInputVolumes(self):
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

    def main(self):
        self.ExtractCeneterLine()
        print("number of point along the centerline = ", self.NPoints)

        InputVolumesDict = self.ReadInputVolumes()
        Inflow_Contrast_Temporal = np.empty([self.VolumeFileNames.__len__(),1])

        count = 0
        for volume in InputVolumesDict:
            print(f"--- Reading the Inflow center of {volume[0]}")
            Inflow_Clip = self.SphereClip(self.CLFile.GetPoint(0), volume[1])
            Inflow_Contrast_Temporal[count] = Inflow_Clip[1]
            count += 1
        InflowClipFileName = 'InflowClip.vtu'
        WriteVTUFile(f"{self.Args.InputFolderName}/{self.OutputFolderName}/{InflowClipFileName}", Inflow_Clip[0])
        # Read the centerline of the every files in the upslope
        CenterLineContrastDict = {volume[0]: None for volume in InputVolumesDict[self.Args.delay:self.Args.peak+1]}

        count = 1
        for volume in InputVolumesDict:
            print(f"--- Reading the centerline of the volume: {volume[0]}")
            Clip1 = None #Clear Clip1 before reassignment
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
        WriteVTUFile(f"{self.Args.InputFolderName}/{self.OutputFolderName}/{ClipOutputFile}", PeakCylinderClip)

        #CenterLine_Contrast_ = self.CenterLineMeshSection(ClipOutputFile)
        
        self.ContrastDisperssion(Inflow_Contrast_Temporal, CenterLineContrastDict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputFolderName", "--InputFolderName", type=str, required=True, dest="InputFolderName")
    parser.add_argument("-HeartBeat", "--HeartBeat", type=float, required=True, dest="HeartBeat")
    parser.add_argument("-delay", "--delay", type=int, required=True, dest="delay", help="bolus time delay before upslope rises")
    parser.add_argument("-peak", "--peak", type=int, required=True, dest="peak", help="the number of image on witch the upslope reaches its peak")

    args = parser.parse_args()

    #CylinderClipAlongCL(args).ExtractPieceWiseCeneterLine()
    #CylinderClipAlongCL(args).ReadInputVolumes()
    CylinderClipAlongCL(args).main()


