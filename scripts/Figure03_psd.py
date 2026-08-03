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

list_dirname = ["yaoxpic_v25_counter_1","yaoxpic_v25_counter_2","yaoxpic_v25_counter_3"]


list_run_id = [1,2,3]


############################################################


dirid=1
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


de    = CGS["c"]/wpe
rhon0 = (wpe/CGS["e"])**2*CGS["me"]/4.0/numpy.pi
J0    = CGS["e"]*rhon0*vthe0*CGS["c"]
B0    = CGS["me"]*CGS["c"]/CGS["e"]*wce
    
print("*"*20)


############################################################
############################################################

Nx = nx
Ny = ny


xmin,xmax =-22,22
ymin,ymax =-0.5,4.5

xticks=numpy.arange(-20,20+10,10)
xticklabel=[r"$%d$"%(tmp) for tmp in xticks]

yticks=numpy.arange(0,4+1,1)
yticklabel=[r"$%d$"%(tmp) for tmp in yticks]


xmin2,xmax2 =0,22
ymin2,ymax2 =-25,-10

xticks2=numpy.arange(0,20+5,5)
xticklabel2=[r"$%d$"%(tmp) for tmp in xticks2]

#yticks2=numpy.arange(-25,-10+5,5)
#yticklabel2=[r"$%d$"%(tmp) for tmp in yticks2]



vmin1,vmax1=-10,2

vmin2,vmax2=-11,1


lwid = 1.0


cmap = mpl.cm.jet
#cmap = mpl.cm.seismic
#cmap = mpl.cm.rainbow

#cmap = mpl.cm.turbo


for tid in [0]:

    Timestep=list_Timestep_FFT[tid]
    
    print("*"*65)
    print("Timestep = %d : %.4f - %.4f"%(Timestep,Timestep*dt*wpe,FFT_Sample_Num*Dt*wpe))

    h5name_R="fft2d_timestep_{:d}_period_{:d}dt_samples_{:d}_wkxwky_eletbz.h5".format(Timestep,FFT_Sample_Num_dt,FFT_Sample_Num)

    norm_k=wpe/CGS["c"]
    norm_w=wpe

    unit_kxy_symbol=r"\cdot d_e"
    unit_kw_symbol=r"/\omega_{pe}"

    kx=numpy.fft.fftshift(numpy.fft.fftfreq(Ny,Dx))*2.0*numpy.pi/norm_k
    ky=numpy.fft.fftshift(numpy.fft.fftfreq(Nx,Dx))*2.0*numpy.pi/norm_k
    kw=numpy.fft.fftshift(numpy.fft.fftfreq(FFT_Sample_Num,Dt))*2.0*numpy.pi/norm_w

    print("kx = ",numpy.min(kx),numpy.max(kx))
    print("ky = ",numpy.min(ky),numpy.max(ky))
    print("kw = ",numpy.min(kw),numpy.max(kw))


    K_TMP=numpy.arange(xmin,xmax+0.01,0.01)*norm_k


    ############################################################
    ############################################################


    hfig = plt.figure(figsize=(16,13))

    margin=[0.05,0.06,0.04,0.06,0.03,0.06]
    barbox=[0.016,0.008,0.65]
    windows=[4,3]
    size=[1,1]



    ############################################################
    ############################################################
    rid,cid = 0,0

    axes_pos = [ 0.0500,  0.7800,  0.2767,  0.1800]
    
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos)

    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))



    dir_tmp = list_dirname[cid]
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


    if cid==2:
       axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       haxe_bar = hfig.add_axes(axes_pos_bar)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)

        
    ##### wave mode

    haxe.annotate(r"$Beam\textendash Beam$",xy=(-4.5,1.8),xytext=(-18.0,1.2),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")

    haxe.annotate(r"$ion\textendash inertia$",xy=(11,0.1),xytext=(12.5,0.55),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")

    haxe.annotate(r"$ion\textendash inertia$",xy=(-12,0.1),xytext=(-18,0.55),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")


    #haxe.annotate(r"$IA$",xy=(15,0.1),xytext=(10.5,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    #haxe.annotate(r"$IA$",xy=(-15,0.1),xytext=(-13,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")


    #haxe.annotate(r"$hybrid$",xy=(7,0.2),xytext=(7.5,0.65),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")


    #haxe.annotate(r"$EA$",xy=(14.5,0.1),xytext=(9.0,0.45),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")

    #haxe.annotate(r"$beam$",xy=(10,2.8),xytext=(12.5,2.35),arrowprops=dict(facecolor="g",edgecolor="g",width=0.2,headwidth=3.0,headlength=3.0),color="g")


        
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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{l}$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    if rid==0:
       haxe.set_title(r"$Run%d$"%(list_run_id[cid]),fontsize=24)
        


    haxe.set_xlabel(r"$k_{\parallel}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)




    ############################################################
    ############################################################
    rid,cid = 1,0

    axes_pos = [ 0.0500,  0.5400,  0.2767,  0.1800]

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos)

    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))



    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Et_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy_vis.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin2,vmax=vmax2)


    if cid==2:
       axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       haxe_bar = hfig.add_axes(axes_pos_bar)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)


    ##### wave mode

    haxe.annotate(r"$T$",xy=(-3.5,3.3),xytext=(-7.0,2.8),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)


        
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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{t}$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    haxe.set_xlabel(r"$k_{\perp}%s$"%(unit_kxy_symbol),fontsize=24)
    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)



    ############################################################
    ############################################################
    rid,cid = 2,0

    axes_pos = [ 0.0500,  0.3000,  0.2767,  0.1800]

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos)

    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))



    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Bz_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy_vis.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin2,vmax=vmax2)


    if cid==2:
       axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       haxe_bar = hfig.add_axes(axes_pos_bar)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)


    ##### wave mode

    haxe.annotate(r"$T$",xy=(-3.5,3.3),xytext=(-7.0,2.8),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)


        
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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)B_z$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    haxe.set_xlabel(r"$k_{\perp}%s$"%(unit_kxy_symbol),fontsize=24)
    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)





    
    ############################################################
    ############################################################
    rid,cid = 3,0

    axes_pos = [ 0.0500,  0.0600,  0.2767,  0.1800]

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))


    ######################################## beam-beam

    m = 0
    wm_real = numpy.copy(wroots_real_run1[:,m])

    index = wm_real>=0

    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="r",label=r"$beam\textendash beam$")




    m = 2
    wm_real = numpy.copy(wroots_real_run1[:,m])

    index = wm_real<=0
    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="--",linewidth=lwid,color="m")

    index = wm_real>=0
    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="m")


    '''
    m = 3
    wm_real = numpy.copy(wroots_real_run1[:,m])

    index = wm_real<=0
    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="--",linewidth=lwid,color="g")
    
    index = wm_real>=0
    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="g")

    m = 4
    wm_real = numpy.copy(wroots_real_run1[:,m])

    index = wm_real<=0
    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="--",linewidth=lwid,color="g")
    
    index = wm_real>=0
    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="g")
    '''


    m = 5
    wm_real = numpy.copy(wroots_real_run1[:,m])

    index = wm_real<=0
    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="--",linewidth=lwid,color="m")
    
    index = wm_real>=0
    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="m")



    ######################################## Beam

    wBeam0 = vd0_run1*k
    haxe.plot(k/norm_k, wBeam0/norm_w,linestyle="--",linewidth=lwid,color="g")

    wBeam1 = vd1_run1*k
    haxe.plot(k/norm_k, wBeam1/norm_w,linestyle="--",linewidth=lwid,color="g")


    ######################################## IA

    #haxe.plot(k/norm_k,wIA_run3/norm_w,linestyle="-",linewidth=lwid,color="b",label=r"$Ion\textendash Acoustic$")


    ######################################## MHD waves

    haxe.plot(k/norm_k,wT_run1/norm_w,color=YXColorBlue,linestyle="-",linewidth=1.0,label=r"$T$")


    ##### wave mode

    haxe.annotate(r"$Beam\textendash Beam$",xy=(-4.5,1.8),xytext=(-18.0,1.2),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")

    haxe.annotate(r"$ion\textendash inertia$",xy=(11,0.1),xytext=(12.5,0.55),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")

    #haxe.annotate(r"$IA$",xy=(15,0.1),xytext=(10.5,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    haxe.annotate(r"$ion\textendash inertia$",xy=(-12,0.1),xytext=(-18,0.55),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")

    haxe.annotate(r"$beam$",xy=(11,3.0),xytext=(13.5,2.4),arrowprops=dict(facecolor="g",edgecolor="g",width=0.2,headwidth=3.0,headlength=3.0),color="g")


    haxe.annotate(r"$T$",xy=(-3.5,3.3),xytext=(-7.0,2.8),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)


    ######################################## harmonic

    #k2 = numpy.arange(-30,30+0.01,0.01)*norm_k

    #w  = numpy.ones(len(k2))
    #haxe.plot(k2/norm_k,w,linestyle="--",linewidth=1.0,color="r")

    #w  = 2.0*numpy.ones(len(k2))
    #haxe.plot(k2/norm_k,w,linestyle="--",linewidth=1.0,color="r")

    #w  = 3.0*numpy.ones(len(k2))
    #haxe.plot(k2/norm_k,w,linestyle="--",linewidth=1.0,color="r")


    #if cid>0:
    #   haxe.annotate(r"$F$",xy=(3.8,1.1),xytext=(6.0,1.4),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")
    #
    #   haxe.annotate(r"$H$",xy=(3.8,2.1),xytext=(6.0,2.4),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")
    #
    #   haxe.annotate(r"$H_3$",xy=(3.8,3.1),xytext=(6.0,3.4),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.grid(linestyle="--",linewidth=0.2,color="w")


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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)

    haxe.set_xlabel(r"$k%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)





    



    
    ############################################################
    ############################################################
    rid,cid = 0,1

    axes_pos = [ 0.3567,  0.7800,  0.2767,  0.1800]
    
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos)


    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))



    dir_tmp = list_dirname[cid]
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


    if cid==2:
       axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       haxe_bar = hfig.add_axes(axes_pos_bar)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)

        
    ##### wave mode

    #haxe.annotate(r"$Langmuir$",xy=(1.5,1.2),xytext=(-8,1.65),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.annotate(r"$Beam\textendash Beam$",xy=(-5.0,1.5),xytext=(-10.5,2.15),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")

    haxe.annotate(r"$ion\textendash inertia$",xy=(-12,0.1),xytext=(-18,0.55),arrowprops=dict(facecolor="lime",edgecolor="lime",width=0.2,headwidth=3.0,headlength=3.0),color="lime")
    
    #haxe.annotate(r"$IA$",xy=(15,0.1),xytext=(10.5,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    #haxe.annotate(r"$IA$",xy=(-15,0.1),xytext=(-13,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")


    haxe.annotate(r"$hybrid$",xy=(8,0.55),xytext=(10,0.9),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")


        
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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{l}$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    if rid==0:
       haxe.set_title(r"$Run%d$"%(list_run_id[cid]),fontsize=24)
        


    haxe.set_xlabel(r"$k_{\parallel}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)




    ############################################################
    ############################################################
    rid,cid = 1,1

    axes_pos = [ 0.3567,  0.5400,  0.2767,  0.1800]

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos)


    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))



    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    print(os.path.join(path_data_tmp,h5name_R))

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Et_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy_vis.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin2,vmax=vmax2)


    if cid==2:
       axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       haxe_bar = hfig.add_axes(axes_pos_bar)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)


    ##### wave mode

    haxe.annotate(r"$T$",xy=(-3.5,3.3),xytext=(-7.0,2.8),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)

        
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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{t}$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    haxe.set_xlabel(r"$k_{\perp}%s$"%(unit_kxy_symbol),fontsize=24)
    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)



    
    ############################################################
    ############################################################
    rid,cid = 2,1
    axes_pos = [ 0.3567,  0.3000,  0.2767,  0.1800]

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos)

    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))



    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    print(os.path.join(path_data_tmp,h5name_R))

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Bz_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy_vis.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin2,vmax=vmax2)


    if cid==2:
       axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       haxe_bar = hfig.add_axes(axes_pos_bar)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)


    ##### wave mode

    haxe.annotate(r"$T$",xy=(-3.5,3.3),xytext=(-7.0,2.8),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)

        
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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)B_z$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    haxe.set_xlabel(r"$k_{\perp}%s$"%(unit_kxy_symbol),fontsize=24)
    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)



    ############################################################
    ############################################################
    rid,cid = 3,1
    axes_pos = [ 0.3567,  0.0600,  0.2767,  0.1800]

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    #axes_pos[0]+=0.05
    haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")


    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))



    ######################################## beam-beam

    m = 0
    wm_real = numpy.copy(wroots_real_run2[:,m])

    index = wm_real>=0

    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="r",label=r"$beam\textendash beam$")


    m = 2
    wm_real = numpy.copy(wroots_real_run2[:,m])

    index = wm_real<0
    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="--",linewidth=lwid,color="m")

    index = wm_real>=0
    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="m",label=r"$beam\ \omega=v_{d2}k$")

    xtmp=k[index]/norm_k
    ytmp=wm_real[index]/norm_w

    spline = scipy.interpolate.UnivariateSpline(xtmp,ytmp,s=5.0)
    xtmp = numpy.linspace(0.0,numpy.max(xtmp),400)
    ytmp = spline(xtmp)

    haxe.plot(xtmp,ytmp,linestyle="--",linewidth=lwid,color="m")


    '''
    m = 3
    wm_real = numpy.copy(wroots_real_run2[:,m])

    index = wm_real>=0

    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="g")
    '''


    m = 4
    wm_real = numpy.copy(wroots_real_run2[:,m])

    index = wm_real>=0

    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="lime")

    m = 5
    wm_real = numpy.copy(wroots_real_run2[:,m])

    index = wm_real>=0

    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="coral")



    ######################################## Beam

    wBeam0 = vd0_run2*k
    haxe.plot(k/norm_k, wBeam0/norm_w,linestyle="--",linewidth=lwid,color="g")

    wBeam1 = vd1_run2*k
    haxe.plot(k/norm_k, wBeam1/norm_w,linestyle="--",linewidth=lwid,color="g")


    ######################################## IA

    #haxe.plot(k/norm_k,wIA_run3/norm_w,linestyle="-",linewidth=lwid,color="b",label=r"$Ion\textendash Acoustic$")



    ######################################## MHD waves

    haxe.plot(k/norm_k,wT_run2/norm_w,color=YXColorBlue,linestyle="-",linewidth=1.0,label=r"$T$")


    ##### wave mode

    #haxe.annotate(r"$Langmuir$",xy=(1.5,1.2),xytext=(-8,1.65),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.annotate(r"$Beam\textendash Beam$",xy=(-5.0,1.5),xytext=(-10.5,2.15),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")

    haxe.annotate(r"$ion\textendash inertia$",xy=(-12,0.1),xytext=(-18,0.55),arrowprops=dict(facecolor="lime",edgecolor="lime",width=0.2,headwidth=3.0,headlength=3.0),color="lime")

    #haxe.annotate(r"$IA$",xy=(15,0.1),xytext=(10.5,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    #haxe.annotate(r"$IA$",xy=(-15,0.1),xytext=(-13,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")


    haxe.annotate(r"$hybrid$",xy=(8,0.55),xytext=(10,0.9),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")




    haxe.annotate(r"$beam$",xy=(10,2.8),xytext=(12.0,2.1),arrowprops=dict(facecolor="g",edgecolor="g",width=0.2,headwidth=3.0,headlength=3.0),color="g")


    haxe.annotate(r"$T$",xy=(-3.5,3.3),xytext=(-7.0,2.8),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)


    ######################################## harmonic

    k2 = numpy.arange(-30,30+0.01,0.01)*norm_k

    #w  = numpy.ones(len(k2))
    #haxe.plot(k2/norm_k,w,linestyle="--",linewidth=1.0,color="r")

    #w  = 2.0*numpy.ones(len(k2))
    #haxe.plot(k2/norm_k,w,linestyle="--",linewidth=1.0,color="r")

    #w  = 3.0*numpy.ones(len(k2))
    #haxe.plot(k2/norm_k,w,linestyle="--",linewidth=1.0,color="r")


    #if cid>0:
    #   haxe.annotate(r"$F$",xy=(3.8,1.1),xytext=(6.0,1.4),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")
    #
    #   haxe.annotate(r"$H$",xy=(3.8,2.1),xytext=(6.0,2.4),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")
    #
    #   haxe.annotate(r"$H_3$",xy=(3.8,3.1),xytext=(6.0,3.4),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.grid(linestyle="--",linewidth=0.2,color="w")


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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)

    haxe.set_xlabel(r"$k%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)


    


    
    ############################################################
    ############################################################
    rid,cid = 0,2

    axes_pos = [ 0.6633,  0.7800,  0.2767,  0.1800]
    
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos)

    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))



    dir_tmp = list_dirname[cid]
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


    if cid==2:
       axes_pos_bar = [ 0.9560,  0.8115,  0.0080,  0.1170]

       #axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       haxe_bar = hfig.add_axes(axes_pos_bar)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)

       print("rid,cid      = (%d, %d)"%(rid,cid))
       print("axes_pos_bar = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos_bar[0],axes_pos_bar[1],axes_pos_bar[2],axes_pos_bar[3]))





        
    ##### wave mode

    haxe.annotate(r"$Langmuir\textendash Beam$",xy=(1.5,1.3),xytext=(-14,2.0),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")




    haxe.annotate(r"$IA$",xy=(16,-0.08),xytext=(12,-0.6),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    haxe.annotate(r"$IA$",xy=(-15,0.1),xytext=(-13,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")


    haxe.annotate(r"$beam\textendash modified$",xy=(6.5,0.7),xytext=(8,0.1),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")




    haxe.annotate(r"$F$",xy=(-2.0,0.8),xytext=(-5,0.3),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.annotate(r"$H$",xy=(-2.5,1.7),xytext=(-6,1.2),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")



    #haxe.text(6.5,1.1,r"$L_1$",color=YXColorBlue,fontsize=20)

    #haxe.text(10,2.0,r"$L_2$",color=YXColorBlue,fontsize=20)

    #haxe.text(13.5,2.9,r"$L_3$",color=YXColorBlue,fontsize=20)

    #haxe.text(17,3.8,r"$L_4$",color=YXColorBlue,fontsize=20)


    kF=4.4
    wF=0.96
    haxe.text(1.0*kF,1.0*wF,r"$L_1$",color=YXColorBlue,fontsize=20)
    haxe.text(2.0*kF,2.0*wF,r"$L_2$",color=YXColorBlue,fontsize=20)
    haxe.text(3.0*kF,3.0*wF,r"$L_3$",color=YXColorBlue,fontsize=20)
    haxe.text(4.0*kF,4.0*wF,r"$L_4$",color=YXColorBlue,fontsize=20)



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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{l}$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    if rid==0:
       haxe.set_title(r"$Run%d$"%(list_run_id[cid]),fontsize=24)
        


    haxe.set_xlabel(r"$k_{\parallel}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)




    ############################################################
    ############################################################
    rid,cid = 1,2
    axes_pos = [ 0.6633,  0.5400,  0.2767,  0.1800]

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos)

    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))




    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Et_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy_vis.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin2,vmax=vmax2)


    if cid==2:
       axes_pos_bar = [ 0.9560,  0.5715,  0.0080,  0.1170]

       #axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       haxe_bar = hfig.add_axes(axes_pos_bar)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)

       print("rid,cid      = (%d, %d)"%(rid,cid))
       print("axes_pos_bar = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos_bar[0],axes_pos_bar[1],axes_pos_bar[2],axes_pos_bar[3]))



    ##### wave mode

    haxe.annotate(r"$T$",xy=(-3.5,3.3),xytext=(-7.0,2.8),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)

    haxe.annotate(r"$F$",xy=(-2.0,0.8),xytext=(-5,0.3),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.annotate(r"$H$",xy=(-2.5,1.7),xytext=(-6,1.2),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")


        
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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{t}$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    haxe.set_xlabel(r"$k_{\perp}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)




    ############################################################
    ############################################################
    rid,cid = 2,2
    axes_pos = [ 0.6633,  0.3000,  0.2767,  0.1800]

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos)

    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))



    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Bz_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy_vis.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin2,vmax=vmax2)


    if cid==2:
       axes_pos_bar = [ 0.9560,  0.3315,  0.0080,  0.1170]

       #axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       haxe_bar = hfig.add_axes(axes_pos_bar)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)

       print("rid,cid      = (%d, %d)"%(rid,cid))
       print("axes_pos_bar = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos_bar[0],axes_pos_bar[1],axes_pos_bar[2],axes_pos_bar[3]))


    
    ##### wave mode

    haxe.annotate(r"$T$",xy=(-3.5,3.3),xytext=(-7.0,2.8),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)

    haxe.annotate(r"$F$",xy=(-2.0,0.8),xytext=(-5,0.3),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.annotate(r"$H$",xy=(-2.5,1.7),xytext=(-6,1.2),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")


        
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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)B_z$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    haxe.set_xlabel(r"$k_{\perp}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)


    
    ############################################################
    ############################################################
    rid,cid = 3,2
    axes_pos = [ 0.6633,  0.0600,  0.2767,  0.1800]

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

    print("rid,cid  = (%d, %d)"%(rid,cid))
    print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))




    ######################################## Bohm-Gross

    haxe.plot(k/norm_k,wL_run3/norm_w,linestyle="-.",linewidth=1.5,color="r",label=r"$Langmuir$")


    ######################################## beam-beam

    m = 0
    wm_real = numpy.copy(wroots_real_run3[:,m])

    index = wm_real>=0

    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="r",label=r"$beam\textendash beam$")


    m = 2
    wm_real = numpy.copy(wroots_real_run3[:,m])

    index = wm_real<0
    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="--",linewidth=lwid,color="m")

    index = wm_real>=0
    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="m",label=r"$beam\ \omega=v_{d2}k$")

    xtmp=k[index]/norm_k
    ytmp=wm_real[index]/norm_w

    spline = scipy.interpolate.UnivariateSpline(xtmp,ytmp,s=5.0)
    xtmp = numpy.linspace(0.0,numpy.max(xtmp),400)
    ytmp = spline(xtmp)

    haxe.plot(xtmp,ytmp,linestyle="--",linewidth=lwid,color="m")


    '''
    m = 3
    wm_real = numpy.copy(wroots_real_run3[:,m])

    index = wm_real>=0

    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="g")
    '''


    m = 4
    wm_real = numpy.copy(wroots_real_run3[:,m])

    index = wm_real>=0

    haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="-",linewidth=lwid,color="lime")

    ######################################## Beam

    wBeam0 = vd0_run3*k
    haxe.plot(k/norm_k, wBeam0/norm_w,linestyle="--",linewidth=lwid,color="g")

    wBeam1 = vd1_run3*k
    haxe.plot(k/norm_k, wBeam1/norm_w,linestyle="--",linewidth=lwid,color="g")




    ######################################## IA

    haxe.plot(k/norm_k,wIA_run3/norm_w,linestyle="-",linewidth=lwid,color="b",label=r"$Ion\textendash Acoustic$")



    ######################################## MHD waves

    haxe.plot(k/norm_k,wT_run3/norm_w,color=YXColorBlue,linestyle="-",linewidth=1.0,label=r"$T$")


    ######################################## thermal

    k2 = numpy.arange(-1000,1000+1,1)*0.01*norm_k


    wthermal = k2*vthe1

    haxe.plot(k2/norm_k,wthermal/norm_w,color="coral",linestyle="--",linewidth=1.0,label=r"$\omega=v_{the2}k$")




    ##### wave mode

    haxe.annotate(r"$Langmuir$",xy=(13,1.1),xytext=(8,1.7),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.annotate(r"$Langmuir\textendash Beam$",xy=(1.5,1.3),xytext=(-14,2.0),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")


    haxe.annotate(r"$IA$",xy=(16,-0.08),xytext=(12,-0.6),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    haxe.annotate(r"$IA$",xy=(-15,0.1),xytext=(-13,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    haxe.annotate(r"$beam\textendash modified$",xy=(8.5,0.8),xytext=(10.5,0.15),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")


    haxe.annotate(r"$beam$",xy=(10,2.8),xytext=(12.0,2.1),arrowprops=dict(facecolor="g",edgecolor="g",width=0.2,headwidth=3.0,headlength=3.0),color="g")

    haxe.annotate(r"$T$",xy=(-3.5,3.3),xytext=(-7.0,2.8),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)

    
    haxe.text(0.5,-0.4,r"$\omega=v_{the}\cdot k$",color="coral",rotation=8)





    ######################################## harmonic

    #k2 = numpy.arange(-20,20+0.01,0.01)*norm_k

    #w  = numpy.ones(len(k2))
    #haxe.plot(k2/norm_k,w,linestyle="--",linewidth=1.0,color="r")

    #w  = 2.0*numpy.ones(len(k2))
    #haxe.plot(k2/norm_k,w,linestyle="--",linewidth=1.0,color="r")


    #haxe.annotate(r"$F$",xy=(-2.0,0.8),xytext=(-5,0.3),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    #haxe.annotate(r"$H$",xy=(-2.5,1.7),xytext=(-6,1.2),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")



    haxe.grid(linestyle="--",linewidth=0.2,color="w")


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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)

    haxe.set_xlabel(r"$k%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)

    

    
    ############################################################
    #plt.show()
    
    fig_path = path_fig
    fig_name = os.path.splitext(os.path.basename(__file__))[0]+"_Timestep_%04d"%(Timestep)


    yaoxpy_vis.fig_save(plt,figpath=fig_path,figname=fig_name,extension=fig_fmts)