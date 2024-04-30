"""This script is developed to validate the interpolation function implemented
on CT-MPI data.
"""

import os
import argparse

import numpy as np
import vtk
import matplotlib.pyplot as plt

from ..CT_MPI_Tools.CylinderClip.py