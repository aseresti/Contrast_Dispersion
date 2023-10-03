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
import re
import argparse
from scripts.ContrastDispersion import ContrastDispersion
class ContrastMain():
    def __init__(self,args):
        self.Args = args
    def main(self):
        foldernames = sorted(glob.glob(f'{self.Args.InputPath}/Re*'))
        MeshPath = f'{foldernames[0]}/results_AdvectionDiffusion/*.xdmf'
        Mesh = ReadXDMFFile(glob.glob(MeshPath)[0])
        Diameter = Mesh.GetBounds()[3]-Mesh.GetBounds()[2]
        viscosity = 0.04
        dcdx = []
        dcdt = []
        velocity = []
        RealVelocity = []
        Re = []
        for folder in foldernames:
            self.Args.InputFolder = f'{folder}/Results_AdvectionDiffusion/*.xdmf'
            [dc_dt,dc_dx,velocity_] = ContrastDispersion(self.Args).main()
            dcdt.append(dc_dt)
            dcdx.append(dc_dx)
            velocity.append(velocity_)
            Re_ = re.findall('r\d+', folder)
            Re.append(Re_)
            RealVelocity.append(viscosity*Re/Diameter)
        df = pd.DataFrame({'Simulation Re #': Re, 'Assigned Velocity': RealVelocity, 'Velocity': velocity})
        df.to_csv(f'{self.Args.OutputPath}')
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="This script takes the ")
    parser.add_argument('-InputPath', '--InputPath', type=str, required=True, dest="InputPath", help="path to the simulations")
    parser.add_argument('-OutputPath', '--OutputPath', type=str, required=True, dest="OutputPath", help="output file path")
    args = parser.parse_args()
    ContrastMain(args).main()