import os,sys

import YaoxPy_PIC as yaoxpy_vis

############################################################
############################################################

import numpy

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.serif": ["Times"],
    "axes.axisbelow":True
})


import h5py

import numpy



############################################################

CGS={'me': 9.1094e-28,
     'mp': 1.6726231e-24,
     'mpi0': 2.406175207e-25,
     'e': 4.8032068e-10,
     'c': 29979245800.0,
    }


############################################################

YXColorBlue = "#003171" #Prussian blue color (Konjo-iro)
YXColorRed  = "#FF0000"


symbol_axis = ["x","y","z"]


fig_fmts = [[".png",350],[".pdf",550]]
#fig_fmts = [[".png",350]]



############################################################
############################################################

path_ws      = os.path.dirname(os.getcwd())
path_fig     = yaoxpy_vis.dir_mk(os.path.join(path_ws,"fig"))
path_data    = yaoxpy_vis.dir_mk(os.path.join(path_ws,"data"))


print("path_data =",path_data)
print("path_fig  =",path_fig)

path_data_pic    = os.path.join(path_data,"yaoxpic_v25_counter_simu")


list_dirname = ["yaoxpic_v25_counter_1","yaoxpic_v25_counter_2","yaoxpic_v25_counter_3"]


#src_dirname  = list_dirname[0]
#path_pic     = path_data+"/%s/yaoxpic_py/data"%(src_dirname)
#print("path_pic  =",path_pic)


Timestep_ALL = 10000


############################################################
############################################################
#from shutil import copyfile

#if not os.path.exists(os.path.join(path_data,"parameters.h5")):
#   copyfile(os.path.join(path_pic,"parameters.h5"),os.path.join(path_data,"parameters.h5"))



############################################################
# FFT
############################################################

Timestep_PIC = 10500
Timestep_FFT = 50

NS_DT   = 10

NS_FFT  = 1024

#Dx     = dx
#Dt     = NS_DT*dt

list_Timestep_FFT=numpy.arange(0,Timestep_PIC-NS_FFT*NS_DT,Timestep_FFT)


list_Timestep_FFT=[0,50]


list_Timestep_Particle=numpy.arange(0,Timestep_PIC+500,500)





############################################################
############################################################

from YaoxPy_Wave_Equations import (
    dispersion_relation_langmuir,
    growth_rate_langmuir,
    dispersion_relation_ion_acoustic,
    growth_rate_ion_acoustic,
    dispersion_relation_MHD_transverse,
)



vthe0 = 0.03*CGS["c"]
vthe1 = 0.03*CGS["c"]

mu    = 1836

wpe   = 5e9

wpi   = wpe/numpy.sqrt(mu)

vthi  = vthe0/numpy.sqrt(mu)


wpeL  = 0.0875*wpe 

############################################################

norm_w = wpe
norm_k = wpe/CGS["c"]


kanchor = 5.0

kmax    = 20.0
dk      = 0.002
N       = int(numpy.ceil(kmax/dk))
k       = numpy.arange(-N,N+1, dtype=float)*dk*norm_k



#################### Run1

alpha0_run1 = 0.5
alpha1_run1 = 0.5

vd0_run1    = -0.3*CGS["c"]
vd1_run1    =  0.3*CGS["c"]

wpe0_run1   = wpe*numpy.sqrt(alpha0_run1)
wpe1_run1   = wpe*numpy.sqrt(alpha1_run1)


filename = "data_wave_equation_reduced_Run1.h5"
H5FILE_R=h5py.File(os.path.join(path_data,filename),"r")
k_reduced_run1           = H5FILE_R["/roots/k"][()]
kb_reduced_run1          = H5FILE_R["/roots/kb"][()]
wroots_real_reduced_run1 = H5FILE_R["/roots/real"][()]
wroots_imag_reduced_run1 = H5FILE_R["/roots/imag"][()]
H5FILE_R.close()


filename = "data_wave_equation_full_maxwellian_Run1.h5"
H5FILE_R=h5py.File(os.path.join(path_data,filename),"r")
k_run1           = H5FILE_R["/roots/k"][()]
kb_run1          = H5FILE_R["/roots/kb"][()]
wroots_real_run1 = H5FILE_R["/roots/real"][()]
wroots_imag_run1 = H5FILE_R["/roots/imag"][()]
H5FILE_R.close()


cs,wIA_run1= dispersion_relation_ion_acoustic(k,wpe0_run1,vthe0,wpi)

wIA_linear_run1 = cs*numpy.abs(k)


wT_run1  = dispersion_relation_MHD_transverse(k,wpe)



#################### Run2

alpha0_run2 = 4.0/5.0
alpha1_run2 = 1.0/5.0

vd0_run2    = -0.075*CGS["c"]
vd1_run2    =  0.3*CGS["c"]

wpe0_run2   = wpe*numpy.sqrt(alpha0_run2)
wpe1_run2   = wpe*numpy.sqrt(alpha1_run2)


filename = "data_wave_equation_reduced_Run2.h5"
H5FILE_R=h5py.File(os.path.join(path_data,filename),"r")
k_reduced_run2           = H5FILE_R["/roots/k"][()]
kb_reduced_run2          = H5FILE_R["/roots/kb"][()]
wroots_real_reduced_run2 = H5FILE_R["/roots/real"][()]
wroots_imag_reduced_run2 = H5FILE_R["/roots/imag"][()]
H5FILE_R.close()


filename = "data_wave_equation_full_maxwellian_Run2.h5"
H5FILE_R=h5py.File(os.path.join(path_data,filename),"r")
k_run2           = H5FILE_R["/roots/k"][()]
kb_run2          = H5FILE_R["/roots/kb"][()]
wroots_real_run2 = H5FILE_R["/roots/real"][()]
wroots_imag_run2 = H5FILE_R["/roots/imag"][()]
H5FILE_R.close()



cs,wIA_run2 = dispersion_relation_ion_acoustic(k,wpe0_run2,vthe0,wpi)

wIA_linear_run2 = cs*numpy.abs(k)


wT_run2  = dispersion_relation_MHD_transverse(k,wpe)




#################### Run3

alpha0_run3 = 20.0/21.0
alpha1_run3 = 1.0/21.0

vd0_run3    = -0.015*CGS["c"]
vd1_run3    =  0.3*CGS["c"]

wpe0_run3   = wpe*numpy.sqrt(alpha0_run3)
wpe1_run3   = wpe*numpy.sqrt(alpha1_run3)


filename = "data_wave_equation_reduced_Run3.h5"
H5FILE_R=h5py.File(os.path.join(path_data,filename),"r")
k_reduced_run3           = H5FILE_R["/roots/k"][()]
kb_reduced_run3          = H5FILE_R["/roots/kb"][()]
wroots_real_reduced_run3 = H5FILE_R["/roots/real"][()]
wroots_imag_reduced_run3 = H5FILE_R["/roots/imag"][()]
H5FILE_R.close()


filename = "data_wave_equation_full_maxwellian_Run3.h5"
H5FILE_R=h5py.File(os.path.join(path_data,filename),"r")
k_run3           = H5FILE_R["/roots/k"][()]
kb_run3          = H5FILE_R["/roots/kb"][()]
wroots_real_run3 = H5FILE_R["/roots/real"][()]
wroots_imag_run3 = H5FILE_R["/roots/imag"][()]
H5FILE_R.close()



cs,wIA_run3 = dispersion_relation_ion_acoustic(k,wpe0_run3,vthe0,wpi)

wIA_linear_run3 = cs*numpy.abs(k)


wT_run3  = dispersion_relation_MHD_transverse(k,wpe)

wL_run3  = dispersion_relation_langmuir(k,wpe0_run3,vd0_run3,vthe1)


wBeam_run3 = vd1_run3*k



############################################################
print("*"*65)
print("*"*65)