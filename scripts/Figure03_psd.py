from YaoxPy_Import_CWD import *

from YaoxPy_Wave_Equations_Two_Electrons import *

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


list_parameters = yaoxpy.pic_parameter_read(path_pic_tmp,path_pic_tmp)
    
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




mu    = 1836
vthe  = vthe0

ud_para = 0.0
ud_perp = 0.2

gamma_v = 1.0/numpy.sqrt(1.0-(ud_para**2+ud_perp**2))

print("gamma = %.4f, 1/gamma = %.4f, 1 - 1/gamma = %.4f"%(gamma_v,1.0/gamma_v,1.0-1.0/gamma_v))


print("*"*65)

print("w_UH = %.8f wce"%(numpy.sqrt(wpe**2+wce**2)/wpe))

print("w_X  = %.8f wce"%(0.5*(wce+numpy.sqrt(4.0*wpe**2+wce**2))/wpe))



############################################################


#alpha0 = 0.5
#alpha1 = 0.5

#vd0    = -0.2*CGS["c"]
#vd1    = 0.2*CGS["c"]

vthe0  = 0.03*CGS["c"]
vthe1  = 0.03*CGS["c"]

mu     = 1836

wpe    = 5e9

wpi    = wpe/numpy.sqrt(mu)


list_color04 = ["r","coral","m"]


list_color12 = ["r","g","g","coral","m","grey","grey","m","coral","g","g","r"]

YXColorBlue = "#003171"


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


line_wid = 1.0


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

    margin=[0.04,0.06,0.04,0.1,0.03,0.06]
    barbox=[0.016,0.008,0.65]
    windows=[4,3]
    size=[1,1]



    ############################################################
    ############################################################
    rid,cid = 0,0
    
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos = [0.04,0.79,0.28,0.17]
    haxe=hfig.add_axes(axes_pos)


    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/El_wkx"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin1,vmax=vmax1)


    #if cid==2:
    #   axes_pos = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
    #   haxe_bar = hfig.add_axes(axes_pos)
    #   hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)

        
    ##### wave mode

    #haxe.annotate(r"$Langmuir$",xy=(-2.5,1.2),xytext=(-14,1.3),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")


    haxe.annotate(r"$Beam\textendash Beam$",xy=(-3.5,1.4),xytext=(-18.5,1.0),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")

    haxe.annotate(r"$IA$",xy=(15,0.1),xytext=(10.5,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    haxe.annotate(r"$IA$",xy=(-15,0.1),xytext=(-13,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")


    haxe.annotate(r"$hybrid$",xy=(7,0.2),xytext=(7.5,0.65),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")


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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{\iota}$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    if rid==0:
       haxe.set_title(r"$Run%d$"%(list_run_id[cid]),fontsize=24)
        


    haxe.set_xlabel(r"$k_{\parallel}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)




    ############################################################
    ############################################################
    rid,cid = 1,0
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos = [0.04,0.56,0.28,0.17]
    haxe=hfig.add_axes(axes_pos)


    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Et_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin2,vmax=vmax2)


    #if cid==2:
    #   axes_pos = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
    #   haxe_bar = hfig.add_axes(axes_pos)
    #   hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)


    ##### wave mode

    haxe.annotate(r"$T$",xy=(-3.5,2.9),xytext=(-7.0,2.5),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)


        
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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{\tau}$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    haxe.set_xlabel(r"$k_{\perp}%s$"%(unit_kxy_symbol),fontsize=24)
    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)



    ############################################################
    ############################################################
    rid,cid = 2,0
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos = [0.04,0.33,0.28,0.17]
    haxe=hfig.add_axes(axes_pos)


    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Bz_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin2,vmax=vmax2)


    #if cid==2:
    #   axes_pos = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
    #   haxe_bar = hfig.add_axes(axes_pos)
    #   hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)


    ##### wave mode

    haxe.annotate(r"$T$",xy=(-3.5,2.9),xytext=(-7.0,2.5),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)


    r'''
    ######################################## cone

    vd0    = 0.3*CGS["c"]
    #vd1    = 0.3*CGS["c"]

    k = numpy.arange(-30,30+0.01,0.01)*norm_k

    w = k*vd0
    haxe.plot(k/norm_k,w/norm_w,linestyle="-.",linewidth=1.5,color="w",label=r"$Langmuir(BG)$")

    haxe.plot(k/norm_k,-1.0*w/norm_w,linestyle="-.",linewidth=1.5,color="w",label=r"$Langmuir(BG)$")
    '''





        
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
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos = [0.04,0.1,0.28,0.17]
    haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")


    ######################################## Bohm-Gross

    #k = numpy.arange(-30,30+0.01,0.01)*norm_k

    #w = numpy.sqrt(wpe*wpe+3.0*numpy.power(k*vthe0,2.0))

    #haxe.plot(k/norm_k,w/norm_w,linestyle="-.",linewidth=1.5,color="r",label=r"$Langmuir(BG)$")


    ########################################
    #list_color04 = ["b","coral","g"]

    list_color04=["r","coral","m"]

    alpha0_run1 = 0.5
    alpha1_run1 = 0.5
    vd0_run1    = -0.3*CGS["c"]
    vd1_run1    =  0.3*CGS["c"]

    k = numpy.arange(-3000,3000+1,1)*0.01*norm_k

    w_run1 = wave_equation_two_electrons_twelveth_solve(wpe, mu, alpha0_run1, alpha1_run1, vd0_run1, vd1_run1, vthe0, vthe1, k,
    a=1.0, b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8)

    list_wavemode = ["Langmuir","beam","beam","beam"+"$-$"+"like","EA"]

    for i in [0]:
        wtmp=numpy.real(w_run1[i,:])
        haxe.plot(k/norm_k,wtmp/norm_w,linestyle="-",linewidth=1.0,color=list_color12[i])

    for i in [1,10]:
        wtmp=numpy.real(w_run1[i,:])

        index=wtmp/norm_w>=0.0
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="-",linewidth=1.0,color=list_color12[i])

        index=numpy.logical_and(wtmp/norm_w<=0.0,wtmp/norm_w>=-0.5)
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="--",linewidth=1.0,color=list_color12[i])

    for i in [4,7]:
        wtmp=numpy.real(w_run1[i,:])

        index=wtmp/norm_w>=0.0
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="-",linewidth=1.0,color=list_color12[i])

        index=numpy.logical_and(wtmp/norm_w<=0.0,wtmp/norm_w>=-0.5)
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="--",linewidth=1.0,color=list_color12[i])


    for i in [4,7]:
        wtmp=numpy.real(w_run1[i,:])

        index=wtmp/norm_w>=-100
        #haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="--",linewidth=1.0,color=list_color12[i])

        xtmp=k[index]/norm_k
        ytmp=wtmp[index]/norm_w

        spline = scipy.interpolate.UnivariateSpline(xtmp,ytmp,s=5.0)
        xtmp = numpy.linspace(0.0,numpy.max(xtmp),400)
        ytmp = spline(xtmp)

        haxe.plot(xtmp,ytmp,linestyle="-.",linewidth=1.0,color=list_color12[i])





    ######################################## IA
    w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0_run1,alpha1_run1,vthe0,vthe1,k)
    haxe.plot(k/norm_k,w/norm_w,linestyle="-",linewidth=1.0,color="b",label=r"$IA$")



    ######################################## MHD waves
    theta=0

    W_TMP=yaoxpy.plasma_waves_MHD(K_TMP,wpe,wce,mu,CGS["c"],theta)


    haxe.plot(K_TMP/norm_k,W_TMP[:,1]/norm_w,color=YXColorBlue,linestyle="-",linewidth=1.0,label=r"$T$")



    ######################################## thermal mode

    #k2 = numpy.arange(-30,30+0.01,0.01)*norm_k

    #w  = 3.0*vthe0*k2
    #haxe.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=1.0,color="grey")
    #haxe.plot(k2/norm_k,-1.0*w/norm_w,linestyle="--",linewidth=1.0,color="grey",label=r"$\omega=\pm 3v_{the0}\cdot k$")

    #w = numpy.power(3.0*numpy.power(vthe0*k2*wpe,2.0),0.25)
    #haxe.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=1.0,color="grey")
    #haxe.plot(k2/norm_k,-1.0*w/norm_w,linestyle="--",linewidth=1.0,color="grey",label=r"$\omega=\pm (\sqrt{3}\omega_{pe}v_{the0}\cdot k)^{1/2}$")



    #haxe.legend(loc="upper center",frameon=True,fontsize=10)


    ##### wave mode

    #haxe.annotate(r"$Langmuir$",xy=(-2.5,1.2),xytext=(-14,1.3),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")


    haxe.annotate(r"$Beam\textendash Beam$",xy=(-3.5,1.4),xytext=(-18.5,1.0),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")


    haxe.annotate(r"$IA$",xy=(15,0.1),xytext=(10.5,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    haxe.annotate(r"$IA$",xy=(-15,0.1),xytext=(-13,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")


    haxe.annotate(r"$hybrid$",xy=(7,0.2),xytext=(7.5,0.65),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")


    #haxe.annotate(r"$EA$",xy=(14.5,0.1),xytext=(9.0,0.45),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")

    haxe.annotate(r"$beam$",xy=(10,2.8),xytext=(12.5,2.35),arrowprops=dict(facecolor="g",edgecolor="g",width=0.2,headwidth=3.0,headlength=3.0),color="g")


    haxe.annotate(r"$T$",xy=(-3.5,2.9),xytext=(-7.0,2.5),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)

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
    
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos = [0.35,0.79,0.28,0.17]
    haxe=hfig.add_axes(axes_pos)


    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/El_wkx"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin1,vmax=vmax1)


    #if cid==2:
    #   axes_pos = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
    #   haxe_bar = hfig.add_axes(axes_pos)
    #   hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)

        
    ##### wave mode

    #haxe.annotate(r"$Langmuir$",xy=(1.5,1.2),xytext=(-8,1.65),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.annotate(r"$Beam\textendash Beam$",xy=(-5.0,1.5),xytext=(-10.5,2.0),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")



    haxe.annotate(r"$IA$",xy=(15,0.1),xytext=(10.5,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    haxe.annotate(r"$IA$",xy=(-15,0.1),xytext=(-13,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")


    haxe.annotate(r"$hybrid$",xy=(7,0.45),xytext=(8.5,0.85),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")

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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{\iota}$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    if rid==0:
       haxe.set_title(r"$Run%d$"%(list_run_id[cid]),fontsize=24)
        


    haxe.set_xlabel(r"$k_{\parallel}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)




    ############################################################
    ############################################################
    rid,cid = 1,1
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos = [0.35,0.56,0.28,0.17]
    haxe=hfig.add_axes(axes_pos)


    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    print(os.path.join(path_data_tmp,h5name_R))

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Et_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin2,vmax=vmax2)


    #if cid==2:
    #   axes_pos = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
    #   haxe_bar = hfig.add_axes(axes_pos)
    #   hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)


    ##### wave mode

    haxe.annotate(r"$T$",xy=(-3.5,2.9),xytext=(-7.0,2.5),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)

        
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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{\tau}$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    haxe.set_xlabel(r"$k_{\perp}%s$"%(unit_kxy_symbol),fontsize=24)
    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)



    
    ############################################################
    ############################################################
    rid,cid = 2,1
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos = [0.35,0.33,0.28,0.17]
    haxe=hfig.add_axes(axes_pos)


    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    print(os.path.join(path_data_tmp,h5name_R))

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Bz_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin2,vmax=vmax2)


    #if cid==2:
    #   axes_pos = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
    #   haxe_bar = hfig.add_axes(axes_pos)
    #   hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)


    ##### wave mode

    haxe.annotate(r"$T$",xy=(-3.5,2.9),xytext=(-7.0,2.5),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)

        
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
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    #axes_pos[0]+=0.05
    axes_pos = [0.35,0.1,0.28,0.17]
    haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")


    ######################################## Bohm-Gross

    #k = numpy.arange(-30,30+0.01,0.01)*norm_k

    #w = numpy.sqrt(wpe*wpe+3.0*numpy.power(k*vthe0,2.0))

    #haxe.plot(k/norm_k,w/norm_w,linestyle="-.",linewidth=1.5,color="r",label=r"$Langmuir(BG)$")


    ########################################
    #list_color04 = ["b","coral","g"]

    list_color04=["r","coral","m"]

    nop0 = 400
    nop1 = 100

    alpha0_run2 = nop0/(nop0+nop1)
    alpha1_run2 = nop1/(nop0+nop1)

    vd0_run2    = -0.3*CGS["c"]*nop1/nop0
    vd1_run2    = 0.3*CGS["c"]


    k = numpy.arange(-3000,3000+1,1)*0.01*norm_k

    w_run2 = wave_equation_two_electrons_twelveth_solve(wpe, mu, alpha0_run2, alpha1_run2, vd0_run2, vd1_run2, vthe0, vthe1, k,
    a=1.0, b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8)

    list_wavemode = ["Langmuir","beam","beam","beam"+"$-$"+"like","EA"]

    for i in [0]:
        wtmp=numpy.real(w_run2[i,:])
        haxe.plot(k/norm_k,wtmp/norm_w,linestyle="-",linewidth=1.0,color=list_color12[i])

    for i in [1,10]:
        wtmp=numpy.real(w_run2[i,:])

        index=wtmp/norm_w>=0.0
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="-",linewidth=1.0,color=list_color12[i])

        index=numpy.logical_and(wtmp/norm_w<=0.0,wtmp/norm_w>=-0.5)
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="--",linewidth=1.0,color=list_color12[i])

    for i in [4,7]:
        wtmp=numpy.real(w_run2[i,:])

        index=wtmp/norm_w>=0.0
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="-",linewidth=1.0,color=list_color12[i])

        index=numpy.logical_and(wtmp/norm_w<=0.0,wtmp/norm_w>=-0.5)
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="--",linewidth=1.0,color=list_color12[i])
    
    
    for i in [4,7]:
        wtmp=numpy.real(w_run2[i,:])

        index=wtmp/norm_w>=-100
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="-.",linewidth=1.0,color=list_color12[i])

        xtmp=k[index]/norm_k
        ytmp=wtmp[index]/norm_w

        spline = scipy.interpolate.UnivariateSpline(xtmp,ytmp,s=5.0)
        xtmp = numpy.linspace(0.0,numpy.max(xtmp),400)
        ytmp = spline(xtmp)

        haxe.plot(xtmp,ytmp,linestyle="-.",linewidth=1.0,color=list_color12[i])
     




    ######################################## IA
    w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0_run2,alpha1_run2,vthe0,vthe1,k)
    haxe.plot(k/norm_k,w/norm_w,linestyle="-",linewidth=1.0,color="b",label=r"$IA$")




    ######################################## beam mode

    #k2   = numpy.arange(0,30+0.01,0.01)*norm_k

    #w = vb*k2
    #haxe.plot(k2/norm_k,w/norm_w,linestyle="-",linewidth=2.0,color="g",label=r"$beam\ \omega=v_{b}\cdot k$")




    ######################################## MHD waves
    theta=0

    W_TMP=yaoxpy.plasma_waves_MHD(K_TMP,wpe,wce,mu,CGS["c"],theta)


    haxe.plot(K_TMP/norm_k,W_TMP[:,1]/norm_w,color=YXColorBlue,linestyle="-",linewidth=1.0,label=r"$T$")



    ######################################## thermal mode

    #k2 = numpy.arange(-30,30+0.01,0.01)*norm_k

    #w  = 3.0*vthe0*k2
    #haxe.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=1.0,color="grey")
    #haxe.plot(k2/norm_k,-1.0*w/norm_w,linestyle="--",linewidth=1.0,color="grey",label=r"$\omega=\pm 3v_{the0}\cdot k$")

    #w = numpy.power(3.0*numpy.power(vthe0*k2*wpe,2.0),0.25)
    #haxe.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=1.0,color="grey")
    #haxe.plot(k2/norm_k,-1.0*w/norm_w,linestyle="--",linewidth=1.0,color="grey",label=r"$\omega=\pm (\sqrt{3}\omega_{pe}v_{the0}\cdot k)^{1/2}$")



    #haxe.legend(loc="upper center",frameon=True,fontsize=10)






    ##### wave mode

    #haxe.annotate(r"$Langmuir$",xy=(1.5,1.2),xytext=(-8,1.65),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.annotate(r"$Beam\textendash Beam$",xy=(-5.0,1.5),xytext=(-10.5,2.0),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")



    haxe.annotate(r"$IA$",xy=(15,0.1),xytext=(10.5,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    haxe.annotate(r"$IA$",xy=(-15,0.1),xytext=(-13,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")


    haxe.annotate(r"$hybrid$",xy=(7,0.55),xytext=(8.5,0.9),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")




    #haxe.annotate(r"$EA$",xy=(14.5,0.1),xytext=(9.0,0.45),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")

    haxe.annotate(r"$beam$",xy=(10,2.8),xytext=(12.5,2.35),arrowprops=dict(facecolor="g",edgecolor="g",width=0.2,headwidth=3.0,headlength=3.0),color="g")


    haxe.annotate(r"$T$",xy=(-3.5,2.9),xytext=(-7.0,2.5),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)


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
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos = [0.66,0.79,0.28,0.17]
    haxe=hfig.add_axes(axes_pos)


    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/El_wkx"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin1,vmax=vmax1)


    if cid==2:
       #axes_pos = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       axes_pos = [0.956,0.81975,0.008,0.1105 ]
       haxe_bar = hfig.add_axes(axes_pos)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)




        
    ##### wave mode

    #haxe.annotate(r"$Langmuir$",xy=(2,1.2),xytext=(-9,2),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.annotate(r"$Beam\textendash Beam$",xy=(1.5,1.3),xytext=(-12.5,2.0),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")




    haxe.annotate(r"$IA$",xy=(15,0.1),xytext=(10.5,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    haxe.annotate(r"$IA$",xy=(-15,0.1),xytext=(-13,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")


    haxe.annotate(r"$hybrid$",xy=(13,0.85),xytext=(15,0.4),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")

    #haxe.annotate(r"$beam$",xy=(10,2.8),xytext=(12.5,2.35),arrowprops=dict(facecolor="g",edgecolor="g",width=0.2,headwidth=3.0,headlength=3.0),color="g")


    haxe.annotate(r"$F$",xy=(-2.0,0.8),xytext=(-5,0.3),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.annotate(r"$H$",xy=(-2.5,1.7),xytext=(-6,1.2),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")



    #haxe.text(6.5,1.1,r"$L_1$",color=YXColorBlue,fontsize=20)

    #haxe.text(10,2.0,r"$L_2$",color=YXColorBlue,fontsize=20)

    #haxe.text(13.5,2.9,r"$L_3$",color=YXColorBlue,fontsize=20)

    #haxe.text(17,3.8,r"$L_4$",color=YXColorBlue,fontsize=20)


    kl=4.2
    wl=0.95
    haxe.text(1.0*kl,1.0*wl,r"$L_1$",color=YXColorBlue,fontsize=20)
    haxe.text(2.0*kl,2.0*wl,r"$L_2$",color=YXColorBlue,fontsize=20)
    haxe.text(3.0*kl,3.0*wl,r"$L_3$",color=YXColorBlue,fontsize=20)
    haxe.text(4.0*kl,4.0*wl,r"$L_4$",color=YXColorBlue,fontsize=20)



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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{\iota}$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    if rid==0:
       haxe.set_title(r"$Run%d$"%(list_run_id[cid]),fontsize=24)
        


    haxe.set_xlabel(r"$k_{\parallel}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)




    ############################################################
    ############################################################
    rid,cid = 1,2
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos = [0.66,0.56,0.28,0.17]
    haxe=hfig.add_axes(axes_pos)


    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Et_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin2,vmax=vmax2)


    if cid==2:
       #axes_pos = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       axes_pos = [0.956,0.58975,0.008,0.1105 ]
       haxe_bar = hfig.add_axes(axes_pos)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)


    ##### wave mode

    haxe.annotate(r"$T$",xy=(-3.5,2.9),xytext=(-7.0,2.5),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)

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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{\tau}$"%(chr(ord('a')+cid),rid+1),color="k",fontsize=20)


    haxe.set_xlabel(r"$k_{\perp}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)




    ############################################################
    ############################################################
    rid,cid = 2,2
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos = [0.66,0.33,0.28,0.17]
    haxe=hfig.add_axes(axes_pos)


    dir_tmp = list_dirname[cid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Bz_wky"
    data2d = H5FILE_R[path_dataset_r][()]
    #data2d = data2d/B0/B0
    data2d = yaoxpy.data_zero_replace(data2d) 
    data2d = numpy.log10(data2d)
    data2d = numpy.fliplr(data2d)

    H5FILE_R.close()

    print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

    him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin2,vmax=vmax2)


    if cid==2:
       #axes_pos = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[1,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
       axes_pos = [0.956,0.35975,0.008,0.1105 ]
       haxe_bar = hfig.add_axes(axes_pos)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)




    r'''
    ######################################## cone
    
    vd1_run3    = 0.3*CGS["c"]

    k = numpy.arange(-30,30+0.01,0.01)*norm_k

    w = k*vd1_run3


    haxe.plot(k/norm_k,w/norm_w,linestyle="-.",linewidth=1.5,color="w",label=r"$Langmuir(BG)$")
    '''



    ##### wave mode

    haxe.annotate(r"$T$",xy=(-3.5,2.9),xytext=(-7.0,2.5),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)

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
    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos = [0.66,0.1,0.28,0.17]
    haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")


    ######################################## Bohm-Gross

    k = numpy.arange(-30,30+0.01,0.01)*norm_k

    wpeL=0.875*wpe

    #w = numpy.sqrt(wpe*wpe+3.0*numpy.power(k*vthe0,2.0))
    w = numpy.sqrt(wpeL*wpeL+3.0*numpy.power(k*vthe0,2.0))

    haxe.plot(k/norm_k,w/norm_w,linestyle="-.",linewidth=1.5,color="r",label=r"$Langmuir(BG)$")


    ########################################
    #list_color04 = ["b","coral","g"]

    list_color04=["r","coral","m"]

    nop0 = 400
    nop1 = 20

    alpha0_run3 = nop0/(nop0+nop1)
    alpha1_run3 = nop1/(nop0+nop1)

    norm_w = wpe*numpy.sqrt(alpha0_run3)

    vd0_run3    = -0.015*CGS["c"]
    vd1_run3    = 0.3*CGS["c"]



    k = numpy.arange(-3000,3000+1,1)*0.01*norm_k

    w_run3 = wave_equation_two_electrons_twelveth_solve(wpe, mu, alpha0_run3, alpha1_run3, vd0_run3, vd1_run3, vthe0, vthe1, k,
    a=1.0, b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8)

    list_wavemode = ["Langmuir","beam","beam","beam"+"$-$"+"like","EA"]

    for i in [0]:
        wtmp=numpy.real(w_run3[i,:])
        haxe.plot(k/norm_k,wtmp/norm_w,linestyle="-",linewidth=1.0,color=list_color12[i])

    for i in [1,10]:
        wtmp=numpy.real(w_run3[i,:])

        index=wtmp/norm_w>=0.0
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="-",linewidth=1.0,color=list_color12[i])

        index=numpy.logical_and(wtmp/norm_w<=0.0,wtmp/norm_w>=-0.5)
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="--",linewidth=1.0,color=list_color12[i])

    for i in [4,7]:
        wtmp=numpy.real(w_run3[i,:])

        index=wtmp/norm_w>=0.0
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="-",linewidth=1.0,color=list_color12[i])

        index=numpy.logical_and(wtmp/norm_w<=0.0,wtmp/norm_w>=-0.5)
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="--",linewidth=1.0,color=list_color12[i])


    for i in [4,7]:
        wtmp=numpy.real(w_run3[i,:])

        index=wtmp/norm_w>=-100
        haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="-.",linewidth=1.0,color=list_color12[i])

        xtmp=k[index]/norm_k
        ytmp=wtmp[index]/norm_w

        spline = scipy.interpolate.UnivariateSpline(xtmp,ytmp,s=5.0)
        xtmp = numpy.linspace(0.0,numpy.max(xtmp),400)
        ytmp = spline(xtmp)

        haxe.plot(xtmp,ytmp,linestyle="-.",linewidth=1.0,color=list_color12[i])


    ######################################## IA
    w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0_run3,alpha1_run3,vthe0,vthe1,k)
    haxe.plot(k/norm_k,w/norm_w,linestyle="-",linewidth=1.0,color="b",label=r"$IA$")



    ######################################## beam mode

    #k2   = numpy.arange(0,30+0.01,0.01)*norm_k

    #w = vb*k2
    #haxe.plot(k2/norm_k,w/norm_w,linestyle="-",linewidth=2.0,color="g",label=r"$beam\ \omega=v_{b}\cdot k$")



    ######################################## MHD waves
    theta=0

    W_TMP=yaoxpy.plasma_waves_MHD(K_TMP,wpe,wce,mu,CGS["c"],theta)


    haxe.plot(K_TMP/norm_k,W_TMP[:,1]/norm_w,color=YXColorBlue,linestyle="-",linewidth=1.0,label=r"$T$")



    ######################################## thermal mode

    #k2 = numpy.arange(-30,30+0.01,0.01)*norm_k

    #w  = 3.0*vthe0*k2
    #haxe.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=1.0,color="grey")
    #haxe.plot(k2/norm_k,-1.0*w/norm_w,linestyle="--",linewidth=1.0,color="grey",label=r"$\omega=\pm 3v_{the0}\cdot k$")

    #w = numpy.power(3.0*numpy.power(vthe0*k2*wpe,2.0),0.25)
    #haxe.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=1.0,color="grey")
    #haxe.plot(k2/norm_k,-1.0*w/norm_w,linestyle="--",linewidth=1.0,color="grey",label=r"$\omega=\pm (\sqrt{3}\omega_{pe}v_{the0}\cdot k)^{1/2}$")



    #haxe.legend(loc="upper center",frameon=True,fontsize=10)






    ##### wave mode

    haxe.annotate(r"$Langmuir$",xy=(13,1.3),xytext=(8,1.7),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")

    haxe.annotate(r"$Beam\textendash Beam$",xy=(1.5,1.2),xytext=(-14.5,1.7),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")


    haxe.annotate(r"$IA$",xy=(15,0.1),xytext=(10.5,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    haxe.annotate(r"$IA$",xy=(-15,0.1),xytext=(-13,0.45),arrowprops=dict(facecolor="b",edgecolor="b",width=0.2,headwidth=3.0,headlength=3.0),color="b")

    haxe.annotate(r"$hybrid$",xy=(13,0.85),xytext=(15,0.4),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m")

    haxe.annotate(r"$beam$",xy=(10,2.8),xytext=(12.5,2.35),arrowprops=dict(facecolor="g",edgecolor="g",width=0.2,headwidth=3.0,headlength=3.0),color="g")


    haxe.annotate(r"$T$",xy=(-3.5,2.9),xytext=(-7.0,2.5),arrowprops=dict(facecolor=YXColorBlue,edgecolor=YXColorBlue,width=0.2,headwidth=3.0,headlength=3.0),color=YXColorBlue)

    
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
    #plt.show()
    
    fig_path = path_fig
    fig_name = os.path.splitext(os.path.basename(__file__))[0]+"_Timestep_%04d"%(Timestep)


    yaoxpy.fig_save(plt,figpath=fig_path,figname=fig_name,extension=fig_fmts)