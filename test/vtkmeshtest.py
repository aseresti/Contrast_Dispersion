import vtk
import sys
import os
path = os.path.abspath("./")
sys.path.append(path)
from tools.utilities import WriteVTUFile, vtk_to_numpy, numpy_to_vtk
import fenics
from tools.ContrastTools import slicer, lowpass
from tools.utilities import ReadXDMFFile
import numpy as np
from scipy.fft import fftn, ifftn, fftfreq
import matplotlib.pyplot as plt
from scipy import special
from skimage import filters

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
Xrange = np.linspace(Xmin_val+(Xmax_val-Xmin_val)/20,Xmin_val+(Xmax_val-Xmin_val)/20+1)
Yrange = np.linspace(Ymin_val,Ymax_val,100)
Zrange = np.linspace(Zmin_val,Zmax_val,100)
for x in Xrange:
    XArray.InsertNextValue(x)
for y in Yrange:
    YArray.InsertNextValue(y)
for z in Zrange:
    ZArray.InsertNextValue(z)
normal = (1., 0., 0.)
origin = (Xmin_val+(Xmax_val-Xmin_val)/20,0.,0.)
slice = slicer(Mesh,normal,origin)
square = vtk.vtkRectilinearGrid()
square.SetDimensions(1, 100, 100)
square.SetXCoordinates(XArray)
square.SetYCoordinates(YArray)
square.SetZCoordinates(ZArray)
PSet = vtk.vtkRectilinearGridToPointSet()
PSet.SetInputData(square)
#square = fenics.UnitSquareMesh(2,2)
Filter = vtk.vtkAppendFilter()
#Filter.SetInputData(PSet.GetOutput())
Filter.AddInputConnection(PSet.GetOutputPort())
Filter.Update()
Prob = vtk.vtkProbeFilter()
Prob.SetInputData(Filter.GetOutput())
Prob.SetSourceData(slice)
Prob.Update()
#PolyData = vtk.vtkAppendPolyData()
#PolyData.AddInputConnection(Prob.GetOutputPort())
#WriteVTPFile('./square.vtp',PolyData.GetOutput())
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('./square.vtu')
writer.SetInputData(Prob.GetOutput())
writer.Update()
ProbOutput = Prob.GetOutput()
Array = ProbOutput.GetPointData().GetArray('f_22')
Array  = vtk_to_numpy(Array)
Arrayfft = abs(fftn(Array))
#Arrayfft = numpy_to_vtk(abs(Arrayfft))
fy = fftfreq(100, (Ymax_val-Ymin_val)/100)
fz = fftfreq(100, (Ymax_val-Ymin_val)/100)

[Fy, Fz] = np.meshgrid(fy,fz)
ax = plt.axes(projection = "3d")
#ax.scatter(Fy, Fz, Array.reshape(100,100))
#plt.show()

def circularLowpassKernel(omega_c, N):  # omega = cutoff frequency in radians (pi is max), N = horizontal size of the kernel, also its vertical size.
  with np.errstate(divide='ignore',invalid='ignore'):
    kernel = np.fromfunction(lambda x, y: omega_c*special.j1(omega_c*np.sqrt((x - (N - 1)/2)**2 + (y - (N - 1)/2)**2))/(2*np.pi*np.sqrt((x - (N - 1)/2)**2 + (y - (N - 1)/2)**2)), [N, N])
  if N % 2:
    kernel[(N - 1)//2, (N - 1)//2] = omega_c**2/(4*np.pi)
  return kernel

kernelN = 6  # Horizontal size of the kernel, also its vertical size.
omega_c = np.pi/4  # Cutoff frequency in radians <= pi
kernel = circularLowpassKernel(omega_c, kernelN)

#ArrayFiltered = fftconvolve(Array.reshape(100,100),kernel)
[Y, Z] = np.meshgrid(YArray,ZArray)
ax2 = plt.axes(projection = "3d")
#ax2.scatter(Y, Z, ArrayFiltered)
#plt.show()

fs = 10000
nyq = 0.5*fs

ArrayFiltered = filters.butterworth(Array,0.02,False,2)
ArrayFiltered = numpy_to_vtk(ArrayFiltered)
Prob.GetOutput().GetPointData().AddArray(ArrayFiltered)
WriteVTUFile("./squarefft.vtu", Prob.GetOutput())
