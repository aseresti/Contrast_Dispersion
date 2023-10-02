import vtk
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fftfreq
from math import pi
from scipy.signal import butter,filtfilt
import glob
from utilities import ReadXDMFFile

def MAFilter(input, window_length):
    FilteredSignal = np.zeros(np.size(input)-window_length)
    for point in range(np.size(input)-window_length):
        FilteredSignal[point] = np.average(input[point:point+window_length])
    return FilteredSignal

def slicer(Mesh,normal,origin):
    plane = vtk.vtkPlane()
    plane.SetNormal(normal)
    plane.SetOrigin(origin)
    slicer = vtk.vtkExtractGeometry()
    slicer.SetExtractInside(1)
    slicer.ExtractBoundaryCellsOn()
    slicer.SetImplicitFunction(plane)
    slicer.SetInputData(Mesh)
    slicer.Update()
    return slicer.GetOutput()

def Model1D(x,array):
    model = LinearRegression()
    model.fit(x.reshape(-1,1),array.reshape(-1,1))
    pred = model.predict(x.reshape(-1,1))
    plt.scatter(x.reshape(-1,1), array.reshape(-1,1))
    plt.plot(x.reshape(-1,1), pred.reshape(-1,1), color = 'red')
    slope = model.coef_[0][0]
    return plt, slope

def ModelPoly(x,array):
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
    return array_fft, darray, d2array

def lowpass(data, cutoff, fs, order):
    nyq = 0.5*fs
    normal_cutoff = cutoff/nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

def ReadResults(foldername,NFiles):
    filenames = glob.glob(f'{foldername}/*.xdmf')
    filenames = sorted(filenames)
    N = len(filenames)
    counter = 0
    for i in range(0,N,int(N/NFiles)):
        counter += 1
        MeshDict = {counter:ReadXDMFFile(filenames[i])}
    return MeshDict