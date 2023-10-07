import vtk
import sys
import os
path = os.path.abspath("./")
sys.path.append(path)
from tools.utilities import WriteVTPFile
import fenics
from tools.ContrastTools import slicer
from tools.utilities import ReadXDMFFile
import numpy as np
path = "/Volumes/CIMBL_HardDrive1/AnahitaAbbasnejad/StenoticPipeSimulations/Re100_P2P1/results_AdvectionDiffusion/Concentration_1000.xdmf"
Mesh = ReadXDMFFile(path)
Xmin_val = Mesh.GetBounds()[0]
Xmax_val = Mesh.GetBounds()[1]
Ymin_val = Mesh.GetBounds()[2]
Ymax_val = Mesh.GetBounds()[3]
Zmin_val = Mesh.GetBounds()[4]
Zmax_val = Mesh.GetBounds()[5]
XArray = vtk.vtkFloatArray()
YArray = vtk.vtkFloatArray()
ZArray = vtk.vtkFloatArray()
for x in np.arange(0,1):
    XArray.InsertNextValue(x)
for y in np.arange(Ymin_val,Ymax_val,100):
    YArray.InsertNextValue(y)
for z in np.arange(Zmin_val,Zmax_val,100):
    ZArray.InsertNextValue(z)
normal = (1., 0., 0.)
origin = (Xmin_val+(Xmax_val-Xmin_val)/2,0.,0.)
slice = slicer(Mesh,normal,origin)
square = vtk.vtkRectilinearGrid()
square.SetDimensions(1, 100, 100)
square.SetXCoordinates(XArray)
square.SetYCoordinates(YArray)
square.SetZCoordinates(ZArray)
#square = fenics.UnitSquareMesh(2,2)
PSet = vtk.vtkRectilinearGridToPointSet()
PSet.SetInputData(square)
Filter = vtk.vtkAppendFilter()
Filter.AddInputConnection(PSet.GetOutputPort())
Filter.Update()
Prob = vtk.vtkProbeFilter()
Prob.SetInputData(Filter.GetOutput())
Prob.SetSourceData(slice)
Prob.Update()
PolyData = vtk.vtkAppendPolyData()
PolyData.AddInputConnection(PSet.GetOutputPort())
WriteVTPFile('./square.vtp',PolyData.GetOutput())
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('./square.vtu')
writer.SetInputData(Prob.GetOutput())
writer.Update()