from YaoxPy_Import_CWD import *

import h5py

from scipy.ndimage import gaussian_filter

sigmax,sigmay=0.5,1.0

############################################################
############################################################

Timestep_ALL            = 10500

Timestep_PIC            = 10500


Timestep_per_Field      = 5
Timestep_per_Particle   = 500

############################################################

FFT_Sample_Num          = 256

#FFT_Sample_Timestep    = 100

FFT_Sample_Num_dt       = 5

#list_Timestep_FFT      = numpy.arange(0,Timestep_PIC-FFT_Sample_Num*FFT_Sample_Num_dt,FFT_Sample_Timestep)

#list_Timestep_Particle = numpy.arange(0,Timestep_PIC+Timestep_per_Particle,Timestep_per_Particle)


list_Timestep_FFT      = [0,1250,4000]

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

wpeL   = 0.875*wpe

print("wpe0 = %.4f wpe"%(wpeL/wpe))



############################################################

Nx = nx
Ny = ny

xmin,xmax =-21,21
ymin,ymax =-0.2,4.2

xticks=numpy.arange(-20,20+10,10)
xticklabel=[r"$%d$"%(tmp) for tmp in xticks]

yticks=numpy.arange(0,4+1,1)
yticklabel=[r"$%d$"%(tmp) for tmp in yticks]




############################################################
############################################################



cmap = mpl.cm.jet
#cmap = mpl.cm.seismic
#cmap = mpl.cm.rainbow

#cmap = mpl.cm.turbo





norm_k = wpe/CGS["c"]
norm_w = wpeL

#norm_w = wpe



unit_kxy_symbol=r"\cdot d_e"
unit_kw_symbol=r"/\omega_{pe}^L"

kx=numpy.fft.fftshift(numpy.fft.fftfreq(Ny,Dx))*2.0*numpy.pi/norm_k
ky=numpy.fft.fftshift(numpy.fft.fftfreq(Nx,Dx))*2.0*numpy.pi/norm_k
kw=numpy.fft.fftshift(numpy.fft.fftfreq(FFT_Sample_Num,Dt))*2.0*numpy.pi/norm_w

print("kx = ",numpy.min(kx),numpy.max(kx))
print("ky = ",numpy.min(ky),numpy.max(ky))
print("kw = ",numpy.min(kw),numpy.max(kw))


############################################################
############################################################

hfig = plt.figure(figsize=(16,13))

margin=[0.05,0.06,0.04,0.06,0.03,0.06]
barbox=[0.016,0.008,0.65]
windows=[3,3]
size=[1,1]




lwid = 1.2

############################################################
############################################################
rid,cid = 0,0

vmin1,vmax1=-7.5,1.5

dir_tmp = list_dirname[0]
path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
path_data_tmp = os.path.join(path_pic_tmp,"data")


for cid in [0,1,2]:

    if cid==0:
       axes_pos = [ 0.0500,  0.7000,  0.2767,  0.2600]
    elif cid==1:
       axes_pos = [ 0.3567,  0.7000,  0.2767,  0.2600]
    elif cid==2:
       axes_pos     = [ 0.6633,  0.7000,  0.2767,  0.2600]
       axes_pos_bar = [ 0.9560,  0.7455,  0.0080,  0.1690]

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))


    Timestep = list_Timestep_FFT[cid]

    h5name_R = "fft2d_timestep_{:d}_period_{:d}dt_samples_{:d}_wkxwky_eletbz.h5".format(Timestep,FFT_Sample_Num_dt,FFT_Sample_Num)


    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/El_wkx"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy_vis.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    data2d_smoothed = gaussian_filter(data2d, sigma=(sigmax,sigmay))

    him=haxe.imshow(data2d_smoothed,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin1,vmax=vmax1)

    if cid==2:
       #axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       haxe_bar = hfig.add_axes(axes_pos_bar)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)

       print("rid,cid      = (%d, %d)"%(rid,cid))
       print("axes_pos_bar = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos_bar[0],axes_pos_bar[1],axes_pos_bar[2],axes_pos_bar[3]))


    ######################################### Bohm-Gross
    haxe.plot(k/norm_k,wL_run3/norm_w,linestyle="-.",linewidth=lwid,color="r",label=r"$Langmuir$")

    haxe.annotate(r"$Langmuir$",xy=(-8,1.25),xytext=(-13,0.75),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")


    ######################################## beam-beam

    m = 0
    wm_real = numpy.copy(wroots_real_run3[:,m])

    index = wm_real>=0

    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="--",linewidth=lwid,color="r",label=r"$beam\textendash beam$")


    haxe.annotate(r"$Langmuir\textendash beam$",xy=(2.0,1.5),xytext=(-13,1.8),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")




    m = 2
    wm_real = numpy.copy(wroots_real_run3[:,m])

    index = wm_real>=0

    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="--",linewidth=lwid,color="m",label=r"$hybrid$")


    haxe.annotate(r"$beam\textendash modified$",xy=(6.0,0.7),xytext=(8.0,0.3),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")




    ######################################## beam
    k2    = numpy.arange(0,30+0.01,0.01)*norm_k
    wbeam = vd1_run3*k2
    haxe.plot(k2/norm_k,wbeam/norm_w,linestyle="--",linewidth=lwid,color="g",label=r"$beam\ \omega=v_{b}\cdot k$")

    haxe.annotate(r"$beam$",xy=(7.0,2.5),xytext=(0.0,2.7),arrowprops=dict(facecolor="g",edgecolor="g",width=0.4,headwidth=4.0,headlength=4.0),color="g")




    #########################################
    if cid==2:
       kL = 4.0530
       wL = 0.8317*wpe/wpeL
       for n in [1,2,3,4]:
           haxe.text(n*kL,n*wL,r"$L_{%d}$"%(n),color="w",fontsize=24)






    #haxe.axhline(y=0.9,linestyle='--',linewidth=1.0,color="w")
    #haxe.axhline(y=1.8,linestyle='--',linewidth=1.0,color="w")


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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.9*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{l}$"%(chr(ord('a')+cid),rid+1),color="w",fontsize=20)



    haxe.set_xlabel(r"$k_{\parallel}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)


    haxe.set_title(r"$t\cdot\omega_{pe}=%.2f-%.2f$"%(Timestep*dt*wpe,(Timestep+FFT_Sample_Num*FFT_Sample_Num_dt)*dt*wpe),fontsize=24)











############################################################
############################################################
rid,cid = 1,0

vmin1,vmax1=-9,-1




dir_tmp = list_dirname[0]
path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
path_data_tmp = os.path.join(path_pic_tmp,"data")


for cid in [0,1,2]:

    if cid==0:
       axes_pos = [ 0.0500,  0.3800,  0.2767,  0.2600]
    elif cid==1: 
       axes_pos = [ 0.3567,  0.3800,  0.2767,  0.2600]
    elif cid==2: 
       axes_pos     = [ 0.6633,  0.3800,  0.2767,  0.2600]
       axes_pos_bar = [ 0.9560,  0.4255,  0.0080,  0.1690]

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")
    
    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))


    Timestep = list_Timestep_FFT[cid]

    h5name_R = "fft2d_timestep_{:d}_period_{:d}dt_samples_{:d}_wkxwky_eletbz.h5".format(Timestep,FFT_Sample_Num_dt,FFT_Sample_Num)


    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Et_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy_vis.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    data2d_smoothed = gaussian_filter(data2d, sigma=(sigmax,sigmay))

    him=haxe.imshow(data2d_smoothed,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin1,vmax=vmax1)

    if cid==2:
       #axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       haxe_bar = hfig.add_axes(axes_pos_bar)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)

       print("rid,cid      = (%d, %d)"%(rid,cid))
       print("axes_pos_bar = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos_bar[0],axes_pos_bar[1],axes_pos_bar[2],axes_pos_bar[3]))


    if cid>0:
       haxe.axhline(y=1.0,linestyle='--',linewidth=0.8,color="w")
       haxe.axhline(y=2.0,linestyle='--',linewidth=0.8,color="w")


    if cid>0:
       haxe.annotate(r"$F$",xy=(-2.0,1.0),xytext=(-5,0.55),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r",fontsize=20)

       haxe.annotate(r"$H$",xy=(-2.5,2.0),xytext=(-6,1.55),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r",fontsize=20)





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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.9*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{t}$"%(chr(ord('a')+cid),rid+1),color="w",fontsize=20)


    haxe.set_xlabel(r"$k_{\perp}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)






############################################################
############################################################
rid,cid = 2,0

vmin1,vmax1=-8,0




dir_tmp = list_dirname[0]
path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
path_data_tmp = os.path.join(path_pic_tmp,"data")


for cid in [0,1,2]:
    
    if cid==0:
       axes_pos = [ 0.0500,  0.0600,  0.2767,  0.2600]
    elif cid==1:
       axes_pos = [ 0.3567,  0.0600,  0.2767,  0.2600]
    elif cid==2:
       axes_pos     = [ 0.6633,  0.0600,  0.2767,  0.2600]
       axes_pos_bar = [ 0.9560,  0.1055,  0.0080,  0.1690]

    

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")
    
    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))

    Timestep = list_Timestep_FFT[cid]

    h5name_R = "fft2d_timestep_{:d}_period_{:d}dt_samples_{:d}_wkxwky_eletbz.h5".format(Timestep,FFT_Sample_Num_dt,FFT_Sample_Num)


    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Bz_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy_vis.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    data2d_smoothed = gaussian_filter(data2d, sigma=(sigmax,sigmay))
    
    him=haxe.imshow(data2d_smoothed,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin1,vmax=vmax1)

    if cid==2:
       #axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       haxe_bar = hfig.add_axes(axes_pos_bar)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)

       print("rid,cid      = (%d, %d)"%(rid,cid))
       print("axes_pos_bar = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos_bar[0],axes_pos_bar[1],axes_pos_bar[2],axes_pos_bar[3]))


    if cid>0:
       haxe.axhline(y=1.0,linestyle='--',linewidth=0.8,color="w")
       haxe.axhline(y=2.0,linestyle='--',linewidth=0.8,color="w")


    if cid>0:
       haxe.annotate(r"$F$",xy=(-2.0,1.0),xytext=(-5,0.55),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r",fontsize=20)

       haxe.annotate(r"$H$",xy=(-2.5,2.0),xytext=(-6,1.55),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r",fontsize=20)





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

    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.9*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)B_{z}$"%(chr(ord('a')+cid),rid+1),color="w",fontsize=20)

    haxe.set_xlabel(r"$k_{\perp}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)





############################################################
#plt.show()
    
fig_path = path_fig
fig_name = os.path.splitext(os.path.basename(__file__))[0]
    
yaoxpy_vis.fig_save(plt,figpath=fig_path,figname=fig_name,extension=fig_fmts)