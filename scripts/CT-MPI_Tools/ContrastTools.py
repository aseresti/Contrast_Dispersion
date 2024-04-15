import vtk
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fftfreq, fft2, ifft2
from math import pi
from scipy.signal import butter,filtfilt
import glob
import re
import os
import sys
from utilities import ReadXDMFFile, vtk_to_numpy

def MAFilter(input, window_length):
    FilteredSignal = np.zeros(np.size(input)-window_length)
    for point in range(np.size(input)-window_length):
        FilteredSignal[point] = np.average(input[point:point+window_length])
    return FilteredSignal

def Slice(Mesh,normal,origin):
    plane = vtk.vtkPlane()
    plane.SetNormal(normal)
    plane.SetOrigin(origin)
    slicer = vtk.vtkExtractGeometry()
    slicer.SetExtractInside(1)
    slicer.ExtractBoundaryCellsOn()
    slicer.SetImplicitFunction(plane)
    slicer.SetInputData(Mesh)
    slicer.Update()
    array_ = slicer.GetOutput().GetPointData().GetArray(0)
    #* Defining a Point Locator to Find the PointID of the Point Close to the Origin of the Slice
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(slicer.GetOutput())
    PointId = locator.FindClosestPoint(origin)
    #* Getting the Average Contrast and CenterLine Value of the Slice
    CneterlineValue = array_.GetValue(PointId)
    array_ = vtk_to_numpy(array_)
    AverageValue = np.average(array_)
    return AverageValue, CneterlineValue

def Model1D(x,array):
    array = np.array(array)
    model = LinearRegression()
    model.fit(x.reshape(-1,1),array.reshape(-1,1))
    pred = model.predict(x.reshape(-1,1))
    slope = model.coef_[0][0]
    return slope, pred

def ModelPoly(x,array):
    array = np.array(array)
    poly = PolynomialFeatures(degree=2, include_bias=False)
    poly_features = poly.fit_transform(x)
    model = LinearRegression()
    model.fit(poly_features, array)
    coef2 = model.coef_[0][1]
    coef1 = model.coef_[0][0]
    d2c_dx2 = 2*coef2
    dc_dx = np.array(coef2*x + coef1)
    dc_dx = dc_dx[0][0]
    pred = model.predict(poly_features)
    plt.scatter(x, array)
    plt.plot(x, pred, color = 'red')
    return plt, dc_dx, d2c_dx2

def Modelfft(x,x_step,array):
    fx = fftfreq(x.shape[0], d = x_step)
    array_fft = np.array(fft(array))
    darray = 1j*2*pi*fx*array_fft
    d2array = 1j*2*pi*fx*darray
    return fx, array_fft, darray, d2array

def lowpass(data, cutoff, fs, order):
    nyq = 0.5*fs
    normal_cutoff = cutoff/nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

def extract_numbers(file):
    return int(re.findall(r'\d+', file)[0])

def get_order(path, ext):
    filenames = os.listdir(path)
    filenames = [f for f in filenames if f"{ext}" in f]
    filenames = sorted(filenames, key = extract_numbers)
    filenames = [f"{path}/{f}" for f in filenames]
    return filenames


def ReadResults(foldername,ext,NFiles):
    filenames = get_order(foldername, ext)
    N = len(filenames)
    counter = 0
    MeshDict = dict()
    for i in range(0,N,int(N/NFiles)):
        MeshDict[counter] = ReadXDMFFile(filenames[i])
        counter += int(N/NFiles)
    #PeakMesh = ReadXDMFFile(filenames[N-1])
    return MeshDict, N

def ReadResults2(foldername,ext):
    filenames = get_order(foldername, ext)
    N = len(filenames)
    counter = 0
    MeshDict = dict()
    for i in range(0,N):
        MeshDict[counter] = ReadXDMFFile(filenames[i])
        counter += 1
    return MeshDict, N
