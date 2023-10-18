import numpy as np
import argparse
#import os
import re
from utilities import GetCentroid, ReadVTIFile, ReadVTPFile, ReadVTUFile, ThresholdInBetween, PrintProgress, numpy_to_vtk, WriteVTUFile

class PixelPickerAlongCenterline():
    def __init__(self,Args):
        self.Args = Args
        
    def main(self):
        
        #Load vti image
        print(f"---Loading vti Image ---: {self.Args.InputImageFile}")
        Image = ReadVTIFile(self.Args.InputImageFile)
        #Load vtp surface model
        #print(f"---Loading volume File ---: {self.Args.InputSurfaceFile}")
        #Surface = ReadVTUFile(self.Args.InputSurfaceFile)
        num = re.findall(r'\d+', self.Args.InputImageFile)
        Output_txt_File = f'Perfusion_{num[-1]}_output.txt'
        outfile = open(Output_txt_File,'w')
        outfile.write("Length, Pixel_Value_Average\n")
        #extract centerline of the lumen
        '''
        CL_FileName = self.Args.InputSurfaceFile.replace(".vtp","_cl.vtp")
        os.system(f"vmtkcenterlines -ifile {self.Args.InputSurfaceFile} -ofile {CL_FileName} -endpoints 1 -resampling 1 -resamplingstep 0.09")
        #extract cross sections along with the centerline
        print("---Extracting the cross sections along the centerline")
        os.system(f"vmtkcenterlinemeshsections -centerlinesfile {CL_FileName} -ifile {self.Args.InputVolumeFile} -ofile ./Section.vtp")
        '''
        #Reading sections
        SurfaceSection = ReadVTPFile("Section.vtp")
        Nstart,Nend = SurfaceSection.GetPointData().GetArray("SectionIds").GetRange()
        Nstart=int(Nstart)
        Nend=int(Nend)
        
        #Find the image points enclosed with the cross section
        progress_ = 0
        for i in range(Nstart,Nend):
            section_ = ThresholdInBetween(SurfaceSection,"SectionIds",i,i)  
            PointScalar_ = np.array([])
            progress_ = PrintProgress(i, Nend, progress_)
            for j in range(section_.GetNumberOfPoints()):
                 point_ = section_.GetPoint(j)
                 pointId_ = Image.FindPoint(point_)
                 intensity_ = Image.GetPointData().GetScalars().GetValue(pointId_)
                 PointScalar_ = np.append(PointScalar_, intensity_)
            
            Section_Averaged_ = np.average(PointScalar_)

            #Calculating the distance of each section from the point of Origin
            Centroid_=GetCentroid(section_)
            if i==Nstart: 
                Dist_=0
                Centroid_old_=Centroid_
            else:
                Dist_+= np.sqrt( (Centroid_[0]-Centroid_old_[0])**2 + (Centroid_[1]-Centroid_old_[1])**2 + (Centroid_[2]-Centroid_old_[2])**2 )
                #updating the output text file with the (Distance, AveragedPixelIntensity)    
                outfile.write(f'{Dist_}, {Section_Averaged_}\n')
                Centroid_old_=Centroid_
        #Projecting the pixel intensity into the volume file
        File1 = ReadVTUFile(self.Args.InputVolumeFile)
        Npoints = File1.GetNumberOfPoints()
        PixelIntensity = np.zeros(Npoints)
        for i in range(0,Npoints):
            PointId = Image.FindPoint(File1.GetPoint(i))
            PixelIntensity[i] = Image.GetPointData().GetScalars().GetValue(PointId)
        PixelIntensityVTK = numpy_to_vtk(PixelIntensity)
        PixelIntensityVTK.SetName("Pixel_Intensity_Value")
        File1.GetPointData().AddArray(PixelIntensityVTK)
        OutputVolumeFile = f'./OutputVolumeFile/OutputVolume_{num[-1]}.vtu'
        WriteVTUFile(OutputVolumeFile, File1)
        #Remove unnecasary files
        #os.system("rm -f Section.vtp RCA_CL.vtp")
        #os.system(f"cp {self.Args.InputVolumeFile} /Users/ana/Documents/AnahitaSeresti/CoronaryPerfusion_Stanford/CABG11B/SVB/OutputVolumeFile")
        
if __name__ == "__main__":
    #Description
    parser = argparse.ArgumentParser( description = "This script will find the averaged pixel intensity of a cross section along the centerline of a lumen")
    
    #InputImageName
    parser.add_argument('-InputSurfaceFile', '--InputSurfaceFile', type = str, required = True, dest = "InputSurfaceFile", help = "File name of the vtp Surface")
    
    #InputSurfaceName
    parser.add_argument('-InputVolumeFile', '--InputVolumeFile', type = str, required = True, dest = "InputVolumeFile", help = "Name of the vtu file containing the mesh")
    
    #InputImageName
    parser.add_argument('-InputImageFile', '--InputImageFile', type = str, required = True, dest = "InputImageFile", help = "Name of the vtu file containing the mesh")
    
    #OutputFileName
    #parser.add_argument('-OutputFileName', '--OutputFileName', type = str, required = True, dest = "OutputFileName", help = "Name of the output text file")
    
    args = parser.parse_args()
    PixelPickerAlongCenterline(args).main()