"""This script takes the registered cycle images and
outputs the averaged imaged"""

import numpy as np
import argparse
from vmtk import vmtkscripts

class CalculateAveragedImage():
    def __init__(self, Args):
        self.Args = Args
        #filenames = glob.glob(f'{self.Args.InputFolderName}/Register_*')
        #filenames = sorted(filenames)
        #self.N_file_per_cycle = int(len(filenames)/self.Args.NumberOfCycles)
        

    def main(self):
        #image_name = f'Cycle_Image_0.vtk'
        #image_path = f'{self.Args.InputFolderName}/{image_name}'
            
            
        #print(f'------------------------------------Cycle_Image_0')
                       
        #reader = vmtkscripts.vmtkImageReader()
        #reader.InputFileName = image_path
        #reader.Execute()
            
        #imageNumpyAdaptor = vmtkscripts.vmtkImageToNumpy()
        #imageNumpyAdaptor.Image = reader.Image
        #imageNumpyAdaptor.Execute()

        #numpyImage = imageNumpyAdaptor.ArrayDict

        #print(numpyImage)
        #Averaged_data = np.array(numpyImage['PointData']['ImageScalars']) 
            
        
        Averaged_data = 0

        for i in range(1,self.Args.NumberOfCycles):
            
            image_name = f'Cycle_Image_{i}_Registered.vtk'
            image_path = f'{self.Args.InputFolderName}/{image_name}'
            
            
            print(f'------------------------------------Cycle_Image_{i}')
                       
            reader = vmtkscripts.vmtkImageReader()
            reader.InputFileName = image_path
            reader.Execute()
            
            imageNumpyAdaptor = vmtkscripts.vmtkImageToNumpy()
            imageNumpyAdaptor.Image = reader.Image
            imageNumpyAdaptor.Execute()

            numpyImage = imageNumpyAdaptor.ArrayDict

            print(numpyImage)
            Averaged_data = Averaged_data + np.array(numpyImage['PointData']['scalars']) 
            
        Averaged_data = Averaged_data/self.Args.NumberOfCycles
        numpyImage['PointData']['scalars'] = Averaged_data
        print(numpyImage)

        output_image = vmtkscripts.vmtkNumpyToImage()
        output_image.ArrayDict = numpyImage # ArrayDict_
        output_image.Execute()

        output_vti = vmtkscripts.vmtkImageWriter()
        output_vti.Image = output_image.Image
        output_vti.OutputFileName = f'{self.Args.OutputFileName}.vti'
        output_vti.Execute()

if __name__ == '__main__':
    #descreption
    parser = argparse.ArgumentParser(description='Thsi script takes a dicom folder with N cycles and outputs an averaged vti image')
    #Input
    parser.add_argument('-InputFolder', '--InputFolder', type = str, required = True, dest = 'InputFolderName', help = 'The name of the folder with all of the dicom files')
    #NumberOfCycles
    parser.add_argument('-NofCycle', '--NumberOfCycles', type = int, required = True, dest = 'NumberOfCycles', help = 'The number of perfusion images that are in the dicom folder')
    #OutputFileName
    parser.add_argument('-OutputFile', '-OutputFileName', type = str, required = True, default = 'Perfusionaveraged', dest = 'OutputFileName', help = 'The name of the output averaged file in vti format')
    args = parser.parse_args()
    CalculateAveragedImage(args).main()

    