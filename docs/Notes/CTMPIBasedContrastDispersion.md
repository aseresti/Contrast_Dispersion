# Contrast Dispersion Pipeline Based on CT-MPI Images

In order to implement the Contrast Dispersion pipeline on CTMPI dataset of the patient, the following steps are crucial:

1. The CT-MPI dataset contains a stack of Dicom images that should be divided into groups based on the number of the cycle images. To learn that number in your dataset, Load the Dicom files on Slicer and count the number of cycle images. Use the [ConvertDicomToVTI.py](../../scripts/CT_MPI_Tools/ConvertDicomToVTI.py) script using the following command to convert the dicom files into several vti image files:

```bash
python scripts/CT_MPI_Tools/ConvertDicomToVTI.py -InputFolder path/to/the/folder/containing/the/dicom/images -NofCycle number_of_cycle_images_in_the_stack
```
