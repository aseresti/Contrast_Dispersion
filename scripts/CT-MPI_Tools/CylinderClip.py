import vtk
import os
import argparse
import numpy as np
from utilities import ReadVTUFile, ReadVTPFile, vtk_to_numpy, WriteVTUFile

class CylinderClipAlongCL():
    def __init__(self, args) -> None:
        self.Args = args
        
        filenames = os.listdir(self.Args.InputFolderName)
        Volumefilenames = [filename for filename in filenames if "vtu" in filename]
        self.VolumeFileNames = sorted(Volumefilenames)

        SurfaceFileName = [filename for filename in filenames if "vtp" in filename]
        self.SurfaceFileName = SurfaceFileName[0]
        self.OutputFileName = "ClippingOutput"
        os.system(f"mkdir {self.Args.InputFolderName}/{self.OutputFileName}")

    def ExtractCeneterLine(self):
        CL_File_Name = "aorta_cl.vtp"
        os.system(f"vmtkcenterlines -ifile {self.Args.InputFolderName}/{self.SurfaceFileName} -ofile {self.Args.InputFolderName}/{self.OutputFileName}/{CL_File_Name} -endpoints 1 -resampling 1 -resamplingstep 15")
        os.system(f"vmtkcenterlinesmoothing -ifile {self.Args.InputFolderName}/{self.OutputFileName}/{CL_File_Name} -ofile {self.Args.InputFolderName}/{self.OutputFileName}/{CL_File_Name} -factor 0.1 -iterations 200")
        os.system(f"vmtkcenterlineresampling -ifile {self.Args.InputFolderName}/{self.OutputFileName}/{CL_File_Name} -ofile {self.Args.InputFolderName}/{self.OutputFileName}/{CL_File_Name} -length 0.05")
        
        CLFile = ReadVTPFile(f"{self.Args.InputFolderName}/{self.OutputFileName}/{CL_File_Name}")
        return CLFile


    def SphereClip(self, center, volume):
        #define the Sphere
        Sphere = vtk.vtkSphere()
        Sphere.SetCenter(center)
        Sphere.SetRadius(2)

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
        print(averagePixelValue)
        return clipper.GetOutput()

    def AppendVolumes(self,volume1,volume2):
        append_filter = vtk.vtkAppendFilter()
        append_filter.AddInputData(volume1)
        append_filter.AddInputData(volume2)
        append_filter.Update()
        return append_filter.GetOutput()
    
    def ContrastDisperssion(self):
        pass

    def ReadInputVolumes(self):
        InputVolumesDict = dict()

        for volume_ in self.VolumeFileNames:
            volume = ReadVTUFile(f"{self.Args.InputFolderName}/{volume_}")
            if volume.GetOutput().GetPointData().GetArray("scalars"):
                InputVolumesDict[f"volume_"] = volume
            else:
                ArrayName = volume.GetPointData().GetArrayName(0)
                volume.GetPointData().GetArray(ArrayName).SetName("scalars")
        
        return sorted(InputVolumesDict.items())

    def main(self):
        CLFile = self.ExtractCeneterLine()
        NPoints = CLFile.GetNumberOfPoints()
        print("number of point along the centerline = ", NPoints)

        InputVolumesDict = self.ReadInputVolumes()
        Inflow_Temporal = np.empty([self.VolumeFileNames.__len__(),1])

        for volume in InputVolumesDict.items():
            Inflow_Clip = self.SphereClip(CLFile.GetPoint(0),volume)

        Clip1 = self.SphereClip(CLFile.GetPoint(0))
        for point in range(1,NPoints-2):
            Clip2 = self.SphereClip(CLFile.GetPoint(point))
            Clip1 = self.AppendVolumes(Clip1, Clip2)
        
        ClipOutputFile = "ClippedVolume.vtu"
        WriteVTUFile(f"{self.Args.InputFolderName}/{self.OutputFileName}/{ClipOutputFile}", Clip1)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputFolderName", "--InputFolderName", type=str, required=True, dest="InputFolderName")
    args = parser.parse_args()

    CylinderClipAlongCL(args).ReadInputVolumes()
    #CylinderClipAlongCL(args).main()


