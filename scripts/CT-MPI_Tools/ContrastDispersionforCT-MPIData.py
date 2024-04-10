"""The purpose of this script is to get the folder containing
the CT-MPI dataset of a patient. The images have to be registered beforehand.
The script asks for the users input by opening a window to slice the model
in the begining and in the end of where the process should be done.
It takes the ceneterline of the segment and drives the contrast dispersion analysis
on that segment.

author: ana
date: Apr 3, 2024
"""

import os
import glob
import vtk
import numpy as np

