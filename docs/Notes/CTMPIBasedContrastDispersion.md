# Contrast Dispersion Pipeline Based on CT-MPI Images

In order to implement the Contrast Dispersion pipeline on CTMPI dataset of the patient, the following steps are crucial:

 1. Convert Dicom Images to VTI: The CT-MPI dataset contains a stack of Dicom images that should be divided into groups based on the number of the cycle images. To learn that number in your dataset, Load the Dicom files on Slicer and count the number of cycle images. Use the [ConvertDicomToVTI.py](../../scripts/CT_MPI_Tools/ConvertDicomToVTI.py) script using the following command to convert the dicom files into several vti image files:

 ```bash
 python scripts/CT_MPI_Tools/ConvertDicomToVTI.py -InputFolder path/to/the/folder/containing/the/dicom/images -NofCycle number_of_cycle_images_in_the_stack
 ```

 2. Register Images: Take the image with at the peak of the time attenuation curve, when contrast reaches the peak, and register all of the other images based on that image on 3DSlicer. Note: to load the images on 3DSlicer as Volume, they have to be in vtk format. After saving them in a seperate folder convert them to vti.

 3. Make a SimVascular project based the peak image you selected in the previous step. Use VMTK toolkit to segement the descending aorta. Load the segmentation into the project and create a mesh with filled caps out of that segmentation.

 4. Projection: Use [ProjectImagetoMesh](../../scripts/CT_MPI_Tools/ImageAnalysisProjectImageToMesh.py) script to project every single one of the vti images in the stack to project them into the mesh you created in SimVascular. Use the script as what follows:

 ```bash
 python scripts/CT_MPI_Tools/ImageAnalysisProjectImageToMesh.py -InputFileName1 /path/to/Image.vti -InputFileName2 /path/to/SimVascularMesh.vtu -OutputFileName descending_aorta_#.vtu
 ```

 5. Extract Surface: extract the surface of the `descending_aorta_#.vtu` using ParaView and save it as a `.vtp` file. Make sure that it includes the word `surface` in its name.

 6. Contrast Dispersion: run the [ContrastDispersion](../../scripts/CT_MPI_Tools/ContrastDispersionAlongVessel.py) script by providing it with the folder containing of the projected volumes, `descending_aorta_#.vtu`, and its surface file, `descending_aorta_surface.vtp` as the following:

 ```bash
 python scripts/CT_MPI_Tools/ContrastDispersionAlongVessel.py -InputFolderName /path/to/the/input/folder -HeartBeat patiant's_heartbeat -delay # -peak #
 ```

 To find delay and peak integer values wich actually indicate the start and point of the upslope, wun the script one time and count the number of points on TAC form 0 to get to the start and end point of the upslope. Note: Upslope must be linear. Do not include the peak on the curve in your counting.