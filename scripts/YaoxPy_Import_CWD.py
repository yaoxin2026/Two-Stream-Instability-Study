import os

import numpy
import scipy
import h5py


import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.serif": ["Times"],
    "axes.axisbelow":True
})


import YaoxPy_Import_Funs as yaoxpy



############################################################
############################################################

CGS={'me': 9.1094e-28,
     'mp': 1.6726231e-24,
     'mpi0': 2.406175207e-25,
     'e': 4.8032068e-10,
     'c': 29979245800.0,
     'eps0': 1,
     'mu0': 1,
     'thomcs': 6.65245e-25,
     'h': 6.6260755e-27,
     'hbar': 1.054572669e-27,
     'kb': 1.380658e-16,
     'r0': 2.8179e-13,
     'G': 6.67428e-8,
    }

############################################################
############################################################

path_ws      = os.path.dirname(os.getcwd())
path_fig     = yaoxpy.dir_mk(os.path.join(path_ws,"fig"))
path_data    = yaoxpy.dir_mk(os.path.join(path_ws,"data"))

print("path_data =",path_data)
print("path_fig  =",path_fig)

path_data_pic    = path_data



fig_fmts = [[".png",350],[".pdf",350]]




############################################################
############################################################
print("*"*65)
print("*"*65)