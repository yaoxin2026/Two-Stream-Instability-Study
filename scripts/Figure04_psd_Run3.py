from YaoxPy_Import_CWD import *

from YaoxPy_Wave_Equations_Two_Electrons import *

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


list_Timestep_FFT      = [0,1500,4000]

print("list_Timestep_FFT      =",list_Timestep_FFT)

############################################################


list_dirname = ["yaoxpic_v25_counter_3"]


############################################################


dirid=0
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



r'''
mu    = 1836
vthe  = vthe0

ud_para = 0.0
ud_perp = 0.2

gamma_v = 1.0/numpy.sqrt(1.0-(ud_para**2+ud_perp**2))

print("gamma = %.4f, 1/gamma = %.4f, 1 - 1/gamma = %.4f"%(gamma_v,1.0/gamma_v,1.0-1.0/gamma_v))


print("*"*65)

print("w_UH = %.8f wce"%(numpy.sqrt(wpe**2+wce**2)/wce))

print("w_X  = %.8f wce"%(0.5*(wce+numpy.sqrt(4.0*wpe**2+wce**2))/wce))
'''


############################################################

nop0 = 400
nop1 = 20

alpha0 = nop0/(nop0+nop1)
alpha1 = nop1/(nop0+nop1)



vd0    = -0.015*CGS["c"]
vd1    = 0.3*CGS["c"]

vthe0  = 0.03*CGS["c"]
vthe1  = 0.03*CGS["c"]

mu     = 1836

wpe    = 5e9

wpi    = wpe/numpy.sqrt(mu)


wpe0   = wpe*numpy.sqrt(alpha0)

wpeL   = 0.875*wpe


print("wpe0 = %.4f wpe"%(wpeL/wpe))
print("wpe0 = %.4f wpe"%(numpy.sqrt(alpha0)*0.9))



list_color04=["r","coral","m"]
list_color08=["r","g","g","coral","m"]


list_color12=["r","g","g","coral","m","grey","grey","b","coral","m","g","r"]



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


dw_nyquist = 2.0*numpy.pi/(FFT_Sample_Num*Dt)

print("dw_nyquist =",dw_nyquist/wpe)

K_TMP=numpy.arange(xmin,xmax+0.01,0.01)*norm_k



############################################################

k = numpy.arange(-3000,3000+1,1)*0.01*norm_k

w_run3 = wave_equation_two_electrons_twelveth_solve(wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, k,
    a=1.0, b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8)




############################################################
############################################################

quiver_y2x = 10.5



hfig = plt.figure(figsize=(16,13))

margin=[0.06,0.06,0.06,0.08,0.03,0.08]
barbox=[0.016,0.008,0.65]
windows=[3,3]
size=[1,1]




############################################################
############################################################
rid,cid = 0,0


vmin1,vmax1=-7.5,1.5



dir_tmp = list_dirname[0]
path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
path_data_tmp = os.path.join(path_pic_tmp,"data")


list_axes_pos=[[0.06,0.70666667,0.27333333,0.23333333],
               [0.36333333,0.70666667,0.27333333,0.23333333],
               [0.66666667,0.70666667,0.27333333,0.23333333]]


for cid in [0,1,2]:

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos=list_axes_pos[cid]
    haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

    Timestep = list_Timestep_FFT[cid]

    h5name_R = "fft2d_timestep_{:d}_period_{:d}dt_samples_{:d}_wkxwky_eletbz.h5".format(Timestep,FFT_Sample_Num_dt,FFT_Sample_Num)


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
       axes_pos = [0.956,0.7475,0.008,0.15166667]
       haxe_bar = hfig.add_axes(axes_pos)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)




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
        haxe.plot(k/norm_k,wtmp/norm_w,linestyle="--",linewidth=1.0,color=list_color12[i])

    haxe.annotate(r"$Beam\textendash Beam$",xy=(-10.5,1.4),xytext=(-14.5,1.75),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")




    ######################################## beam mode
    k2 = numpy.arange(0,30+0.01,0.01)*norm_k
    w  = vd1*k2
    haxe.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=1.0,color="g",label=r"$beam\ \omega=v_{b}\cdot k$")

    haxe.annotate(r"$beam$",xy=(7.5,2.5),xytext=(1.0,2.8),arrowprops=dict(facecolor="g",edgecolor="g",width=0.4,headwidth=4.0,headlength=4.0),color="g")


    ######################################### Bohm-Gross
    k = numpy.arange(-30,30+0.01,0.01)*norm_k
    #w = numpy.sqrt(wpe0*wpe0+3.0*numpy.power(k*vthe0,2.0))

    w = numpy.sqrt(wpeL*wpeL+3.0*numpy.power(k*vthe0,2.0))

    haxe.plot(k/norm_k,w/norm_w,linestyle="-.",linewidth=1.2,color="r")

    haxe.annotate(r"$Langmuir$",xy=(-7.0,0.9),xytext=(-16,0.4),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")


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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.9*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{\iota}$"%(chr(ord('a')),cid+1),color="w",fontsize=20)



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


list_axes_pos=[[0.06,0.39333333,0.27333333,0.23333333],
               [0.36333333,0.39333333,0.27333333,0.23333333],
               [0.66666667,0.39333333,0.27333333,0.23333333]]




for cid in [0,1,2]:

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos=list_axes_pos[cid]
    haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")
    
    Timestep = list_Timestep_FFT[cid]

    h5name_R = "fft2d_timestep_{:d}_period_{:d}dt_samples_{:d}_wkxwky_eletbz.h5".format(Timestep,FFT_Sample_Num_dt,FFT_Sample_Num)


    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Et_wky"
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
       axes_pos = [0.956,0.43416667,0.008,0.15166667]
       haxe_bar = hfig.add_axes(axes_pos)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)


    if cid>0:
       haxe.axhline(y=1.0,linestyle='--',linewidth=0.8,color="w")
       haxe.axhline(y=2.0,linestyle='--',linewidth=0.8,color="w")


    if cid>0:
       haxe.annotate(r"$F$",xy=(-2.0,1.0),xytext=(-5,0.65),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r",fontsize=20)

       haxe.annotate(r"$H$",xy=(-2.5,2.0),xytext=(-6,1.65),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r",fontsize=20)





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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.9*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)E_{\tau}$"%(chr(ord('a')+rid),cid+1),color="w",fontsize=20)

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


list_axes_pos=[[0.06,0.08,0.27333333,0.23333333],
               [0.36333333,0.08,0.27333333,0.23333333],
               [0.66666667,0.08,0.27333333,0.23333333]]




for cid in [0,1,2]:

    #axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
    axes_pos=list_axes_pos[cid]
    haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")
    
    Timestep = list_Timestep_FFT[cid]

    h5name_R = "fft2d_timestep_{:d}_period_{:d}dt_samples_{:d}_wkxwky_eletbz.h5".format(Timestep,FFT_Sample_Num_dt,FFT_Sample_Num)


    H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
    path_dataset_r="/psd/Bz_wky"
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
       
       axes_pos = [0.956,0.12083333,0.008,0.15166667]
       haxe_bar = hfig.add_axes(axes_pos)
       hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)


    if cid>0:
       haxe.axhline(y=1.0,linestyle='--',linewidth=0.8,color="w")
       haxe.axhline(y=2.0,linestyle='--',linewidth=0.8,color="w")


    if cid>0:
       haxe.annotate(r"$F$",xy=(-2.0,1.0),xytext=(-5,0.65),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r",fontsize=20)

       haxe.annotate(r"$H$",xy=(-2.5,2.0),xytext=(-6,1.65),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r",fontsize=20)





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
    haxe.text(xlim_tmp[0]+0.06*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.9*(ylim_tmp[1]-ylim_tmp[0]),r"$(%s%s)B_z$"%(chr(ord('a')+rid),cid+1),color="w",fontsize=20)

    haxe.set_xlabel(r"$k_{\perp}%s$"%(unit_kxy_symbol),fontsize=24)

    if cid==0:
       haxe.set_ylabel(r"$\omega%s$"%(unit_kw_symbol),fontsize=24)





############################################################
#plt.show()
    
fig_path = path_fig
fig_name = os.path.splitext(os.path.basename(__file__))[0]
    
yaoxpy.fig_save(plt,figpath=fig_path,figname=fig_name,extension=fig_fmts)