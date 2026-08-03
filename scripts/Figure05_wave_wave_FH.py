from YaoxPy_Import_CWD import *

import h5py


import scipy

############################################################
############################################################

Timestep_ALL           = 10500

Timestep_PIC           = 10500


Timestep_per_Field     = 5
Timestep_per_Particle  = 500

############################################################

FFT_Sample_Num         = 2048

FFT_Sample_Timestep    = 100

FFT_Sample_Num_dt      = Timestep_per_Field

list_Timestep_FFT      = numpy.arange(0,Timestep_PIC-FFT_Sample_Num*FFT_Sample_Num_dt,FFT_Sample_Timestep)

list_Timestep_Particle = numpy.arange(0,Timestep_PIC+Timestep_per_Particle,Timestep_per_Particle)


list_Timestep_FFT      = [100]

print("list_Timestep_FFT      =",list_Timestep_FFT)
print("list_Timestep_Particle =",list_Timestep_Particle)


############################################################


list_dirname = ["yaoxpic_v25_counter_3"]



############################################################


dirid=0
dir_tmp=list_dirname[dirid]


path_pic_tmp = os.path.join(path_data_pic,dir_tmp)
path_pic_tmp = os.path.join(path_pic_tmp,"data")


print("path_pic_tmp =",path_pic_tmp)


list_parameters = yaoxpy_vis.pic_parameter_read(path_pic_tmp,path_pic_tmp)
    
dt  = list_parameters["dt"]
dx  = list_parameters["dx"]
nx  = list_parameters["nx"]
ny  = list_parameters["ny"]
wpe = list_parameters["wpe"]
wce = list_parameters["wce"]
ncpu     = list_parameters["ncpu"]
nspecies = list_parameters["nspecies"]


print("dt        = %14.4e"%(dt))
print("dx        = %10.4f"%(dx))
print("nx        = %5d"%(nx))
print("ny        = %5d"%(ny))
print("ncpu      = %5d"%(ncpu))
print("nspecies  = %5d"%(nspecies))
print("wpe       = %14.4e"%(wpe))
print("wce       = %14.4e"%(wce))
    

m0     = list_parameters["m0"]
macro0 = list_parameters["macro0"]
vd0    = list_parameters["vd0"]
vthe0  = list_parameters["vth0"]

m2     = list_parameters["m2"]
macro2 = list_parameters["macro2"]
vd2   = list_parameters["vd2"]
vthe2 = list_parameters["vth2"]

Dx    = dx
Dt    = dt*FFT_Sample_Num_dt
    
cs    = CGS["c"]
    
de    = CGS["c"]/wpe
rhon0 = (wpe/CGS["e"])**2*CGS["me"]/4.0/numpy.pi
J0    = CGS["e"]*rhon0*vthe0*CGS["c"]
B0    = CGS["me"]*CGS["c"]/CGS["e"]*wce
    
print("*"*20)

############################################################

wpe0  = wpe0_run3

wpeL  = 0.875*wpe



print("wpe0 = %.4f wpe"%(wpe0/wpe))
print("wpeL = %.4f wpe"%(wpeL/wpe))


############################################################
############################################################

Nx = nx
Ny = ny

xmin,xmax =-21,21
ymin,ymax =-0.2,4.2

xticks=numpy.arange(-20,20+5,5)
xticklabel=[r"$%d$"%(tmp) for tmp in xticks]

yticks=numpy.arange(0,4+1,1)
yticklabel=[r"$%d$"%(tmp) for tmp in yticks]




############################################################
############################################################

vmin1,vmax1=-7.5,1.5


cmap = mpl.cm.jet
#cmap = mpl.cm.seismic
#cmap = mpl.cm.rainbow

#cmap = mpl.cm.turbo


Timestep=0

h5name_R="fft2d_timestep_{:d}_period_{:d}dt_samples_{:d}_wkxwky_eletbz.h5".format(Timestep,FFT_Sample_Num_dt,FFT_Sample_Num)



norm_k = wpe/CGS["c"]
norm_w = wpeL

#norm_w = wpe



unit_kxy_symbol=r"\cdot d_e"
unit_kw_symbol=r"/\omega_{pe1}"

kx=numpy.fft.fftshift(numpy.fft.fftfreq(Ny,Dx))*2.0*numpy.pi/norm_k
ky=numpy.fft.fftshift(numpy.fft.fftfreq(Nx,Dx))*2.0*numpy.pi/norm_k
kw=numpy.fft.fftshift(numpy.fft.fftfreq(FFT_Sample_Num,Dt))*2.0*numpy.pi/norm_w

print("kx = ",numpy.min(kx),numpy.max(kx))
print("ky = ",numpy.min(ky),numpy.max(ky))
print("kw = ",numpy.min(kw),numpy.max(kw))



############################################################
############################################################

quiver_y2x = 9.2



hfig = plt.figure(figsize=(16,11))

margin=[0.06,0.08,0.04,0.08,0.06,0.08]
barbox=[0.03,0.015,0.65]
windows=[3,2]
size=[1,1]


lwid = 1.5


############################################################
############################################################
rid,cid = 1,0
axes_pos = [ 0.0600,  0.0800,  0.4000,  0.5600]

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[2,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))

dir_tmp = list_dirname[0]
path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
path_data_tmp = os.path.join(path_pic_tmp,"data")

H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
path_dataset_r="/psd/El_wkx"
data2d = H5FILE_R[path_dataset_r][()]
#data2d = data2d/B0/B0
data2d = yaoxpy_vis.data_zero_replace(data2d) 
data2d = numpy.log10(data2d)
data2d = numpy.fliplr(data2d)

H5FILE_R.close()

print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin1,vmax=vmax1)


######################################### harmonic

wharmonic = numpy.zeros(len(k))

haxe.plot(k/norm_k,wharmonic+1.0,linestyle="--",linewidth=1.0,color="w")

#haxe.plot(k/norm_k,wharmonic+2.0,linestyle="--",linewidth=1.0,color="w")



######################################### Bohm-Gross

haxe.plot(k/norm_k,wL_run3/norm_w,linestyle="-.",linewidth=lwid,color="r",label=r"$Langmuir$")

haxe.annotate(r"$Langmuir$",xy=(-8,1.25),xytext=(-15,1.0),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")



########################################

m = 0
wm_real = numpy.copy(wroots_real_run3[:,m])

index = wm_real>=0

haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="--",linewidth=lwid,color="r",label=r"$beam\textendash beam$")


haxe.annotate(r"$Langmuir\textendash beam$",xy=(2.0,1.5),xytext=(-8.5,1.7),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r",fontsize=16)




m = 2
wm_real = numpy.copy(wroots_real_run3[:,m])

index = wm_real>=0

haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="--",linewidth=lwid,color="m",label=r"$hybrid$")


haxe.annotate(r"$beam\textendash modified$",xy=(12,1.05),xytext=(13,0.8),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m",fontsize=16)


######################################## Beam

#wBeam0 = vd0_run3*k
#haxe.plot(k/norm_k, wBeam0/norm_w,linestyle="--",linewidth=lwid,color="g")

wBeam1 = vd1_run3*k
haxe.plot(k/norm_k, wBeam1/norm_w,linestyle="--",linewidth=lwid,color="g")


haxe.annotate(r"$beam$",xy=(11.0,3.5),xytext=(12.5,3.3),arrowprops=dict(facecolor="g",edgecolor="g",width=0.2,headwidth=3.0,headlength=3.0),color="g",fontsize=16)


######################################## IA

haxe.plot(k/norm_k,wIA_run3/norm_w,linestyle="--",linewidth=lwid,color="b",label=r"$Ion-Acoustic$")



######################################## coalescence

lx0_L = 4.0
haxe.quiver(0.0,0.0,lx0_L,quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="r")

haxe.text(3.2,0.25,r"$L$",color="r",fontsize=28)



#####

lx0_S = 9.5
haxe.quiver(0.0,0.0,lx0_S,0.0,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="b")

print("S = ",lx0_S*cs/wpe0)


haxe.text(11.0,-0.1,r"$S$",color="b",fontsize=28)


haxe.quiver(lx0_S,0.0,(lx0_L-lx0_S),quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="grey")

haxe.quiver(0.0,0.0,(lx0_L-lx0_S),quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="w")

haxe.text(-6.0,0.4,r"$F$",color="w",fontsize=28)



########################################

haxe.set_xlim(xmin,xmax)
haxe.set_ylim(ymin,ymax)
#haxe.set_zlim(zmin,zmax)

haxe.set_xticks(xticks)
haxe.set_xticklabels(xticklabel)
haxe.set_yticks(yticks)
haxe.set_yticklabels(yticklabel)

haxe.tick_params(direction='in', length=5, width=1, colors='k')

haxe.tick_params(axis="x",labelsize=20)
haxe.tick_params(axis="y",labelsize=20)
        

        

# label
xlim_tmp=haxe.get_xlim()
ylim_tmp=haxe.get_ylim()
haxe.text(xlim_tmp[0]+0.05*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.9*(ylim_tmp[1]-ylim_tmp[0]),r"$(b1)Decay\ Process\ L\to S+F$",color="w",fontsize=20)


haxe.set_xlabel(r"$k_{\parallel}\cdot d_e$",fontsize=24)


haxe.set_ylabel(r"$\omega/\omega_{pe}^L$",fontsize=24)






















############################################################
############################################################
rid,cid = 1,1
axes_pos = [ 0.5200,  0.0800,  0.4000,  0.5600]

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[2,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))



dir_tmp = list_dirname[0]
path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
path_data_tmp = os.path.join(path_pic_tmp,"data")

H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
path_dataset_r="/psd/El_wkx"
data2d = H5FILE_R[path_dataset_r][()]
#data2d = data2d/B0/B0
data2d = yaoxpy_vis.data_zero_replace(data2d) 
data2d = numpy.log10(data2d)
data2d = numpy.fliplr(data2d)

H5FILE_R.close()

print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin1,vmax=vmax1)


axes_pos_bar = [ 0.9500,  0.1780,  0.0150,  0.3640]
#axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[2,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
haxe_bar = hfig.add_axes(axes_pos_bar)
hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)

print("rid,cid      = (%d, %d)"%(rid,cid))
print("axes_pos_bar = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos_bar[0],axes_pos_bar[1],axes_pos_bar[2],axes_pos_bar[3]))





######################################### harmonic

wharmonic = numpy.zeros(len(k))

haxe.plot(k/norm_k,wharmonic+1.0,linestyle="--",linewidth=1.0,color="w")

haxe.plot(k/norm_k,wharmonic+2.0,linestyle="--",linewidth=1.0,color="w")



######################################### Bohm-Gross

haxe.plot(k/norm_k,wL_run3/norm_w,linestyle="-.",linewidth=lwid,color="r",label=r"$Langmuir$")

haxe.annotate(r"$Langmuir$",xy=(-8,1.25),xytext=(-15,1.0),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")



########################################


m = 2
wm_real = numpy.copy(wroots_real_run3[:,m])

index = wm_real>=0

haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="--",linewidth=lwid,color="m")


haxe.annotate(r"$beam\textendash modified$",xy=(12,1.05),xytext=(13,0.8),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m",fontsize=16)





######################################## coalescence


lx0_L=4.0
haxe.quiver(0.0,0.0,lx0_L,quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="r",zorder=4)

haxe.text(3.0,0.3,r"$L$",color="r",fontsize=28)

print("L =",1.0*wpe/wpe0)


lx0_F=-4.5
ly0_F=1.0

haxe.quiver(0.0,0.0,lx0_F,ly0_F*quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="w")

#haxe.text(-6.5,1.3,r"$L-S\to F$",color="w",fontsize=24)



haxe.quiver(lx0_F,ly0_F,lx0_L,quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="grey")


haxe.quiver(0.0,0.0,(lx0_L+lx0_F),(ly0_F+1.0)*quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="w")




haxe.text(-6,0.4,r"$F$",color="w",fontsize=28)


haxe.text(-2.3,2.1,r"$H$",color="w",fontsize=28)








#haxe.legend(loc="upper center",frameon=True,fontsize=10)


#haxe.grid(linestyle="--",linewidth=0.2,color="grey")


        
haxe.set_xlim(xmin,xmax)
haxe.set_ylim(ymin,ymax)
#haxe.set_zlim(zmin,zmax)

haxe.set_xticks(xticks)
haxe.set_xticklabels(xticklabel)
haxe.set_yticks(yticks)
haxe.set_yticklabels(yticklabel)

haxe.tick_params(direction='in', length=5, width=1, colors='k')

haxe.tick_params(axis="x",labelsize=20)
haxe.tick_params(axis="y",labelsize=20)


# label
xlim_tmp=haxe.get_xlim()
ylim_tmp=haxe.get_ylim()
haxe.text(xlim_tmp[0]+0.05*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.9*(ylim_tmp[1]-ylim_tmp[0]),r"$(b2)Coalescence\ Process\ F+L\to H$",color="w",fontsize=20)



haxe.set_xlabel(r"$k_{\parallel}\cdot d_e$",fontsize=24)


#haxe.set_ylabel(r"$\omega/\omega_{pe}^L$",fontsize=24)





############################################################
############################################################


dir_tmp = list_dirname[0]
path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
path_data_tmp = os.path.join(path_pic_tmp,"data")



filename = "timestep_timeseries.h5"

H5FILE_R=h5py.File(os.path.join(path_data_tmp,filename),"r")
dt = H5FILE_R["/dt"][()]
EE = H5FILE_R["/energy/E"][()]
EB = H5FILE_R["/energy/B"][()]
EK = H5FILE_R["/energy/kinetic"][()]

MOMENT0 = H5FILE_R["/species_0/moment"][()]
MOMENT2 = H5FILE_R["/species_2/moment"][()]

timestep_state =H5FILE_R["/timestep_state"][()]
H5FILE_R.close()


print(EE.shape)

EE=numpy.array(EE)
EE=EE.T

EB=numpy.array(EB)
EB=EB.T

EK=numpy.array(EK)
EK=EK.T

time0 = numpy.arange(EE.shape[0])*wpe*dt


index=timestep_state>0
time0=time0[index]
EE=EE[index,:]
EB=EB[index,:]
EK=EK[index,:]


#norm_E = EK[0,2]
norm_E = 1.0

print("norm_E =",norm_E)


##### E+B
DEB = EE[:,0]+EE[:,1]+EE[:,2]+EB[:,0]+EB[:,1]+EB[:,2]
DEB = DEB-DEB[0]

##### Ek-all
DEK = EK[:,0]+EK[:,1]+EK[:,2]+EK[:,3]
DEK = DEK-DEK[0]





h5name_R="energy_extract_wavemode_components_4.h5"


H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")


field_name = 'EL'
ELangmuir = H5FILE_R[f"/energy/{field_name}/Langmuir"][()]
EBeam     = H5FILE_R[f"/energy/{field_name}/Beam"][()]
EIA       = H5FILE_R[f"/energy/{field_name}/IA"][()]

field_name = 'ET'
EH1       = H5FILE_R[f"/energy/{field_name}/H1"][()]
EH2       = H5FILE_R[f"/energy/{field_name}/H2"][()]


for field_name in ["Ez","Bx","By","Bz"]:
    EH1  += H5FILE_R[f"/energy/{field_name}/H1"][()]
    EH2  += H5FILE_R[f"/energy/{field_name}/H2"][()]

H5FILE_R.close()


time = numpy.arange(len(ELangmuir))*FFT_Sample_Num_dt*wpe*dt



############################################################
############################################################
rid,cid = 0,0
axes_pos = [ 0.0600,  0.7200,  0.4000,  0.2400]

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))


##### Ek
xtmp=time0
ytmp=-1.0*DEK/norm_E


fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),1000)
ytmp=fun_intp(xtmp)


index=numpy.logical_and(xtmp>5,xtmp<155)

xtmp=xtmp[index]
ytmp=ytmp[index]

ytmp=numpy.log10(ytmp)

haxe.plot(xtmp,ytmp,linestyle="--",linewidth=2.0,color="r",label=r"$-\Delta \mathcal{E}_{k}$")



##### ELangmuir
xtmp=time
ytmp=ELangmuir/norm_E

fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),1000)
ytmp=fun_intp(xtmp)


xtmp_min=xtmp[ytmp==numpy.min(ytmp)]

index=numpy.logical_and(xtmp>xtmp_min,xtmp<155)

xtmp=xtmp[index]
ytmp=ytmp[index]

ytmp=numpy.log10(ytmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="cyan",label=r"$Langmuir$")




##### Ebeam
xtmp=time
ytmp=EBeam/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),1000)
ytmp=fun_intp(xtmp)


xtmp_min=xtmp[ytmp==numpy.min(ytmp)]

index=numpy.logical_and(xtmp>xtmp_min,xtmp<155)


xtmp=xtmp[index]
ytmp=ytmp[index]

ytmp=numpy.log10(ytmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="g",label=r"$Beam$")



index=ytmp==numpy.max(ytmp)
xtmp_max=xtmp[index][0]
haxe.axvline(x=xtmp_max,linestyle="--",linewidth=0.8,color="g")
haxe.text(xtmp_max-6,2.0,r"$t=%.2f$"%(xtmp_max),rotation=90,color="g")

print("beam = ",xtmp_max)


##### EIA

xtmp=time
ytmp=EIA/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),1000)
ytmp=fun_intp(xtmp)


xtmp_min=xtmp[ytmp==numpy.min(ytmp)]

index=numpy.logical_and(xtmp>xtmp_min,xtmp<155)


xtmp=xtmp[index]
ytmp=ytmp[index]

ytmp=numpy.log10(ytmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="b",label=r"$Ion\textendash acoustic$")



haxe.legend(loc="center right",frameon=True,fontsize=16)



xmin,xmax =0,160
xticks=numpy.arange(0,160+40,40)
xticklabel=[r"$%d$"%(tmp) for tmp in xticks]

ymin,ymax =-0.1,6.1
yticks=numpy.arange(0,7,1)
yticklabel=[r"$10^{%d}$"%(tmp) for tmp in yticks]

haxe.set_xlim(xmin,xmax)
#haxe.set_xticks(xticks)
#haxe.set_xticklabels(xticklabel)


haxe.set_ylim(ymin,ymax)
haxe.set_yticks(yticks)
haxe.set_yticklabels(yticklabel)


haxe.tick_params(direction='in', length=5, width=1, colors='k')

haxe.tick_params(axis="x",labelsize=18)
haxe.tick_params(axis="y",labelsize=18)

# label
xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.04*(xlim[1]-xlim[0]),ylim[0]+0.85*(ylim[1]-ylim[0]),r"$(a1)$",fontsize=20)



haxe.set_xlabel(r"$t\cdot \omega_{pe}$",fontsize=24)

haxe.set_ylabel(r"$\mathcal{E}\ [erg]$",fontsize=24)








############################################################
############################################################
rid,cid = 0,1
axes_pos = [ 0.5200,  0.7200,  0.4000,  0.2400]

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))


xtmp=time
ytmp=EH1/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),1000)
ytmp=fun_intp(xtmp)

xtmp_min=xtmp[ytmp==numpy.min(ytmp)]

index=numpy.logical_and(xtmp>xtmp_min,xtmp<155)



xtmp=xtmp[index]
ytmp=ytmp[index]

ytmp=numpy.log10(ytmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="coral",label=r"$F$")


index=[]
for k in numpy.arange(1,len(ytmp)-1):
   if (ytmp[k]>ytmp[k-1] and ytmp[k]>ytmp[k+1]) and xtmp[k]>30:
      index.append(k)
index=numpy.array(index)


xtmp_max=xtmp[index[0]]
haxe.axvline(x=xtmp_max,linestyle="--",linewidth=0.8,color="coral")
haxe.text(xtmp_max-6,1.0,r"$t=%.2f$"%(xtmp_max),rotation=90,color="coral")






xtmp=time
ytmp=EH2/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),1000)
ytmp=fun_intp(xtmp)

xtmp_min=xtmp[ytmp==numpy.min(ytmp)]

index=numpy.logical_and(xtmp>xtmp_min,xtmp<155)



xtmp=xtmp[index]
ytmp=ytmp[index]

ytmp=numpy.log10(ytmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="m",label=r"$H$")


index=[]
for k in numpy.arange(1,len(ytmp)-1):
   if (ytmp[k]>ytmp[k-1] and ytmp[k]>ytmp[k+1]) and xtmp[k]>30:
      index.append(k)
index=numpy.array(index)



xtmp_max=xtmp[index[0]]
haxe.axvline(x=xtmp_max,linestyle="--",linewidth=0.8,color="m")
haxe.text(xtmp_max+2,-0.2,r"$t=%.2f$"%(xtmp_max),rotation=90,color="m")
print("H2   = ",xtmp_max)






'''
xtmp=time
ytmp=EIA

index=xtmp<=35

print(ytmp[index][-1])

for i in range(len(index)):
    ytmp[i]=i*ytmp[index][-1]/35.0


fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),1000)
ytmp=fun_intp(xtmp)

index=numpy.logical_and(xtmp>5,xtmp<155)
haxe.plot(xtmp[index],ytmp[index],linestyle="-",linewidth=2.0,color="b",label=r"$Ion-Acoustic$")
'''




haxe.legend(loc="center right",frameon=True,fontsize=16)



xmin,xmax =0,160
xticks=numpy.arange(0,160+40,40)
xticklabel=[r"$%d$"%(tmp) for tmp in xticks]


ymin,ymax =-1.1,3.1
yticks=numpy.arange(-1,4,1)
yticklabel=[r"$10^{%d}$"%(tmp) for tmp in yticks]

haxe.set_xlim(xmin,xmax)
#haxe.set_xticks(xticks)
#haxe.set_xticklabels(xticklabel)


haxe.set_ylim(ymin,ymax)
haxe.set_yticks(yticks)
haxe.set_yticklabels(yticklabel)


haxe.tick_params(direction='in', length=5, width=1, colors='k')

haxe.tick_params(axis="x",labelsize=18)
haxe.tick_params(axis="y",labelsize=18)

# label
xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.04*(xlim[1]-xlim[0]),ylim[0]+0.85*(ylim[1]-ylim[0]),r"$(a2)$",fontsize=20)



haxe.set_xlabel(r"$t\cdot \omega_{pe}$",fontsize=24)

#haxe.set_ylabel(r"$\mathcal{E}/\mathcal{E}_{k0}$",fontsize=24)








############################################################
#plt.show()
    
fig_path = path_fig
fig_name = os.path.splitext(os.path.basename(__file__))[0]
    
yaoxpy_vis.fig_save(plt,figpath=fig_path,figname=fig_name,extension=fig_fmts)