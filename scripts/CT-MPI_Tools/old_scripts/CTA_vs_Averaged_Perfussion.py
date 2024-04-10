import numpy as np
from vmtk import vmtkscripts

CTA_path = './CABG4A/CTA_SimVascular_4A/Images/dicom_CTA_diastole.vti'
Perfusion_path = './CABG4A/Perfusion_Averaged_SimVascular_4A/Images/CABG4A.vti'
output_path = './CABG4A/'

reader = vmtkscripts.vmtkImageReader()
reader.InputFileName = CTA_path
reader.Execute()
            
imageNumpyAdaptor = vmtkscripts.vmtkImageToNumpy()
imageNumpyAdaptor.Image = reader.Image
imageNumpyAdaptor.Execute()

numpyImage = imageNumpyAdaptor.ArrayDict

print(numpyImage)
numpyImage['Origin'] = [0, 0, 0]

output_image = vmtkscripts.vmtkNumpyToImage()
output_image.ArrayDict = numpyImage # ArrayDict_
output_image.Execute()

output_vti = vmtkscripts.vmtkImageWriter()
output_vti.Image = output_image.Image
output_vti.OutputFileName = f'{output_path}/CTA_CABG4A.vtk'
output_vti.Execute()



reader = vmtkscripts.vmtkImageReader()
reader.InputFileName = Perfusion_path
reader.Execute()
            
imageNumpyAdaptor = vmtkscripts.vmtkImageToNumpy()
imageNumpyAdaptor.Image = reader.Image
imageNumpyAdaptor.Execute()

numpyImage_2 = imageNumpyAdaptor.ArrayDict

numpyImage_2['Origin'] = [0, 0, 0]
print(numpyImage_2)

output_image = vmtkscripts.vmtkNumpyToImage()
output_image.ArrayDict = numpyImage_2 # ArrayDict_
output_image.Execute()

output_vti = vmtkscripts.vmtkImageWriter()
output_vti.Image = output_image.Image
output_vti.OutputFileName = f'{output_path}/Averaged_CABG4A.vtk'
output_vti.Execute()