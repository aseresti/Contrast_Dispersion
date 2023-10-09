'''
The goal of this script is to take the path to the Simulations
and calculates the velocity based on different models.
Author: Anahita A. Seresti
Date: Oct 2nd, 2023
'''
import pandas as pd
import glob
import sys
import os
path = os.path.abspath("./")
sys.path.append(path)
from tools.utilities import ReadXDMFFile
from tools.ContrastTools import get_order
from math import pi
import re
import argparse
from scripts.ContrastDispersion import ContrastDispersion
from scripts.ContrastDispersionFiltered import ContrastDispersionFilter
class ContrastMain():
    def __init__(self,args):
        self.Args = args
    def main(self):
        foldernames = get_order(f'{self.Args.InputPath}', 'Re')
        #foldernames = sorted(glob.glob(f'{self.Args.InputPath}/Re*'))
        MeshPath = f'{foldernames[0]}/results_AdvectionDiffusion/*.xdmf'
        Mesh = ReadXDMFFile(glob.glob(MeshPath)[0])
        Diameter = Mesh.GetBounds()[3]-Mesh.GetBounds()[2]
        Area = pi*(Diameter/2)**2
        #print(Diameter)
        viscosity = 0.04
        velocity_ContrastDispersion = []
        velocity_Filtered = []
        RealVelocity = []
        Re = []
        for folder in foldernames:
            print(folder)
            self.Args.InputFolder = f'{folder}/results_AdvectionDiffusion'
            [dc_dt,dc_dx,velocity_] = ContrastDispersion(self.Args).main()
            velocity_ContrastDispersion.append(velocity_)
            [dc_dt,dc_dx,velocity_] = ContrastDispersionFilter(self.Args).main()
            velocity_Filtered.append(velocity_)
            Re_ = re.findall(r'\d+',folder)
            Re_ = Re_[1]                           
            Re.append(Re_)
            RealVelocity.append(viscosity*int(Re_)/Diameter)
        AssignedFlow = [v*Area for v in RealVelocity]
        df = pd.DataFrame({'Simulation Re #': Re, 'Assigned Flow': AssignedFlow, 'Assigned Velocity': RealVelocity, 'Velocity': velocity_ContrastDispersion, 'Velocity (Filtered Contrast Method)': velocity_Filtered})
        df.to_csv(f'{self.Args.OutputPath}')
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="This script takes the ")
    parser.add_argument('-InputPath', '--InputPath', type=str, required=True, dest="InputPath", help="path to the simulations")
    parser.add_argument('-OutputPath', '--OutputPath', type=str, required=False, default="./Results.csv", dest="OutputPath", help="output file path")
    args = parser.parse_args()
    ContrastMain(args).main()