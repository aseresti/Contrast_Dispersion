import numpy as np
import argparse
from scipy import interpolate
import glob
import pandas as pd
from scipy.fft import fft2, ifft2, fftfreq
from math import pi

class VelocityAlongCLPerfusion():
    def __init__(self,Args):
        self.Args = Args
        
    def Main(self):
        
        #read input text files
        FolderName = self.Args.InputFolderName
        FileNames  = glob.glob(f'{FolderName}/*.txt')
        FileNames  = sorted(FileNames)

        max_lengths = []
        min_lengths = []
        Length = []
        PixelValue = []
        for filename in FileNames:
            with open(filename) as f:
                data = [line.strip() for line in f]
                length_ = []
                value_ = []
                for line in data:
                    x = line.split(", ")
                    length_.append(x[0])
                    value_.append(x[-1])
                length_.remove(length_[0])
                length_ = [float(y) for y in length_]
                max_lengths.append(max(item for item in length_))
                min_lengths.append(min(item for item in length_))
                Length.append(length_)
                value_.remove(value_[0])
                value_ = [float(z) for z in value_]
                PixelValue.append(value_)
        '''
        #interpolate data along time
        x_step = 0.5
        x_new = np.arange(max(item for item in min_lengths), min(item for item in max_lengths), x_step)
        count = 0
        #Store data into an array
        PixelValue_array = np.zeros((x_new.shape[0],len(PixelValue)))
        #PixelValue_array[:,0] = x_new
        for value in PixelValue:
            f = interpolate.interp1d(Length[count], value)
            Pix = f(x_new)
            PixelValue_array[:,count] = Pix
            count += 1
        '''
        df = pd.DataFrame(PixelValue)#_array)
        df.to_excel(excel_writer = "./PerfusionAlongCL.xlsx")
        '''
        #converting time-space domain to the 2D Fourier domain
        time_step = 0.1
        time_new = np.arange(0,PixelValue_array.shape[1],time_step)
        time_old = np.arange(0,PixelValue_array.shape[1],1)
        PixelValue_array_intr = np.zeros((PixelValue_array.shape[0], time_new.shape[0]))
        #for index in range(0,PixelValue_array.shape[0]):
        f1 = interpolate.interp2d(time_old, x_new, PixelValue_array)
        PixelValue_array_intr = f1(time_new, x_new)
        '''
            
        '''    
        df1 = pd.DataFrame(PixelValue_array_intr)
        df1.to_excel(excel_writer = "./PerfusionAlongCLintr.xlsx")
        '''
        '''    
        #Advection Diffusion in Fourier domain
        D = 0.01
        Contrast_fft = np.array(fft2(PixelValue_array_intr))
        (nx, nt) = PixelValue_array_intr.shape
        (fx, ft) = (fftfreq(nx, d=x_step), fftfreq(nt, d=time_step))
        dC_dt_fft = 1j*2*pi*ft*Contrast_fft
        dC_dx_fft = 1j*2*pi*fx*Contrast_fft.transpose()
        d2C_dx2_fft = 1j*2*pi*fx*dC_dx_fft
        dC_dx_fft = dC_dx_fft.transpose()
        d2C_dx2_fft = d2C_dx2_fft.transpose()
        
        #velocity_fft = np.dot((D*d2C_dx2_fft - dC_dt_fft),dC_dx_fft)
        velocity_fft = (D*d2C_dx2_fft - dC_dt_fft)/dC_dx_fft

        #Back to the time-space domain: Velocity Along the Centerline
        velocity = ifft2(velocity_fft[1:,:])

        f2 = interpolate.interp2d(time_new, x_new[1:], np.absolute(velocity))
        velocity_intr = f2(time_old, x_new[1:])
        
        df1 = pd.DataFrame(velocity_intr)
        df1.to_excel(excel_writer = "./VelocityAlongCL.xlsx")
        #print(np.absolute(velocity))
        '''   
if __name__ == "__main__":
    #Description
    parser = argparse.ArgumentParser(description = "This script gets averaged pixel values along the centerline and calculates velocity")
    
    #InputFolderName
    parser.add_argument('-InputFolderName', '--InputFolderName', type = str, required = True, dest = "InputFolderName", help = "The name of the folder containing the averaged pixel values text file")
    
    args = parser.parse_args()
    VelocityAlongCLPerfusion(args).Main()