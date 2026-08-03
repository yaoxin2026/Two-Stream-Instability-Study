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
############################################################

Nx = nx
Ny = ny


vd1 = vd1_run3


############################################################
############################################################

vmin1,vmax1=-7.5,1.5


cmap  = mpl.cm.jet
#cmap = mpl.cm.seismic
#cmap = mpl.cm.rainbow

#cmap = mpl.cm.turbo


Timestep = 0

h5name_R = "fft2d_timestep_{:d}_period_{:d}dt_samples_{:d}_wkxwky_eletbz.h5".format(Timestep,FFT_Sample_Num_dt,FFT_Sample_Num)



norm_k = wpe/CGS["c"]
norm_w = wpe

#norm_w = wpe



unit_kxy_symbol=r"\cdot d_e"
unit_kw_symbol=r"/\omega_{pe}"

kx=numpy.fft.fftshift(numpy.fft.fftfreq(Ny,Dx))*2.0*numpy.pi/norm_k
ky=numpy.fft.fftshift(numpy.fft.fftfreq(Nx,Dx))*2.0*numpy.pi/norm_k
kw=numpy.fft.fftshift(numpy.fft.fftfreq(FFT_Sample_Num,Dt))*2.0*numpy.pi/norm_w


dkx_nyquist = 2.0*numpy.pi/Dx/Ny
dkw_nyquist = 2.0*numpy.pi/Dt/FFT_Sample_Num

print("kx = ",numpy.min(kx),numpy.max(kx))
print("ky = ",numpy.min(ky),numpy.max(ky))
print("kw = ",numpy.min(kw),numpy.max(kw))

print("dkx_nyquist =",dkx_nyquist/norm_k,kx[1]-kx[0])
print("dkw_nyquist =",dkw_nyquist/norm_w,kw[1]-kw[0])




KKW,KKX=numpy.meshgrid(kw,kx,indexing="ij")


#K_TMP=numpy.arange(xmin,xmax+0.01,0.01)*norm_k


sigma_w  = 0.01*wpe/norm_w




############################################################
############################################################

quiver_y2x = 4.4


hfig = plt.figure(figsize=(19,7))

margin=[0.06,0.02,0.06,0.16,0.06,0.12]
barbox=[0.1,0.015,0.65]
windows=[2,3]
size=[1,1]


############################################################
############################################################
rid,cid = 0,2
axes_pos = [ 0.7133,  0.6100,  0.2667,  0.3300]

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))





dir_tmp = list_dirname[0]
path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
path_data_tmp = os.path.join(path_pic_tmp,"data")

H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
path_dataset_r="/psd/El_wkx"
data2d  = H5FILE_R[path_dataset_r][()]
data2d  = numpy.fliplr(data2d)
#data2d = data2d/B0/B0
#data2d = yaoxpy_vis.data_zero_replace(data2d) 
#data2d = numpy.log10(data2d)
H5FILE_R.close()



print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))


print("data2d   =",data2d.shape)
print("kx,ky,kw =",len(kx),len(ky),len(kw))




list_color=["darkblue","m","g","r"]


list_kf=[]

no=0
for rho in [0.5,0.6,0.7,0.8]:

    WTMP = KKX*norm_k*vd1*rho/norm_w
    DIS  = numpy.abs(KKW-WTMP)/sigma_w
    MASK = numpy.exp(-0.5*DIS**2)

    data2d_mask=MASK*data2d
    data2d_mask=numpy.sum(data2d_mask**2,axis=0)*dkw_nyquist
    data2d_mask=numpy.log10(data2d_mask)

    index=kx>=0

    ktmp=kx[index]
    ptmp=data2d_mask[index]
    #haxe.plot(ktmp,ptmp,color=list_color[no])

    p_smooth  = scipy.ndimage.gaussian_filter1d(ptmp,sigma=5.0,mode='nearest')

    #haxe.plot(ktmp,p_smooth,color=list_color[no])

    pf=numpy.max(p_smooth)
    kf=ktmp[p_smooth==pf][0]
    print("rho, kf,pf = %.2f, %.4f, %.4f"%(rho,kf,pf))

    list_kf.append(kf)
    
    
    #for k in [1,2,3,4]:
    #    haxe.axvline(x=k*kfa,linestyle="--",color="r")

    haxe.plot(ktmp/kf,ptmp,color=list_color[no],label=r"$v_{ph}=%.1fv_{d2},k_F=%.2fd_e^{-1}$"%(rho,kf))

    #haxe.plot(ktmp/kf,p_smooth,color=list_color[no],label=r"$\rho=%.1f:k_F=%.2fd_e^{-1}$"%(rho,kf))

    no+=1







haxe.legend(loc="upper right",frameon=False,fontsize=12)


for i in range(1,5):
    haxe.axvline(x=i,linestyle="--",linewidth=0.5,color="grey")



xmin,xmax=0,5
xticks=numpy.arange(1,5)
haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)


ymin,ymax=-5,20
yticks=numpy.arange(-5,20+5,5)
yticklabel=[r"$10^{%d}$"%(tmp) for tmp in yticks]
haxe.set_ylim(ymin,ymax)
haxe.set_yticks(yticks)
haxe.set_yticklabels(yticklabel)


haxe.tick_params(direction='in', length=5, width=1, colors='k')

haxe.tick_params(axis="x",labelsize=20)
haxe.tick_params(axis="y",labelsize=20)



# label
xlim_tmp=haxe.get_xlim()
ylim_tmp=haxe.get_ylim()
haxe.text(xlim_tmp[0]+0.05*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(c1)$",color="k",fontsize=18)




haxe.set_xlabel(r"$k_{\parallel}/k_F$",fontsize=20)

haxe.set_ylabel(r"$Power\ \mathcal{P}(k_{\parallel})$",fontsize=20)













############################################################
############################################################
rid,cid = 1,2
axes_pos = [ 0.7133,  0.1600,  0.2667,  0.3300]

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))




dir_tmp = list_dirname[0]
path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
path_data_tmp = os.path.join(path_pic_tmp,"data")

H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
path_dataset_r="/psd/El_wkx"
data2d  = H5FILE_R[path_dataset_r][()]
data2d  = numpy.fliplr(data2d)
#data2d = data2d/B0/B0
#data2d = yaoxpy_vis.data_zero_replace(data2d) 
#data2d = numpy.log10(data2d)
H5FILE_R.close()



print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))


print("data2d   =",data2d.shape)
print("kx,ky,kw =",len(kx),len(ky),len(kw))




list_color=["darkblue","m","g","r"]


list_wf=[]

no=0
for rho in [0.5,0.6,0.7,0.8]:

    WTMP = KKX*norm_k*vd1*rho/norm_w
    DIS  = numpy.abs(KKW-WTMP)/sigma_w
    MASK = numpy.exp(-0.5*DIS**2)

    data2d_mask=MASK*data2d
    data2d_mask=numpy.sum(data2d_mask**2,axis=1)*dkx_nyquist
    data2d_mask=numpy.log10(data2d_mask)
    
    index=numpy.logical_and(kw>=0,kw<=6)

    ktmp=kw[index]
    ptmp=data2d_mask[index]
    #haxe.plot(ktmp,ptmp,color=list_color[no])

    
    p_smooth  = scipy.ndimage.gaussian_filter1d(ptmp,sigma=5.0,mode='nearest')

    #haxe.plot(ktmp,p_smooth,linestyle="-",color=list_color[no])
    
    
    pf=numpy.max(p_smooth)
    kf=ktmp[p_smooth==pf][0]
    print("rho, wf,pf = %.2f, %.4f, %.4f"%(rho,kf,pf))


    list_wf.append(kf)
    
    
    haxe.plot(ktmp/kf,ptmp,color=list_color[no],label=r"$v_{ph}=%.1fv_{d2},\omega_F=%.2fd_e^{-1}$"%(rho,kf))
    
    #haxe.plot(ktmp/kf,p_smooth,color=list_color[no],label=r"$\rho=%.1f:\omega_F=%.2fd_e^{-1}$"%(rho,kf))


    no+=1






haxe.legend(loc="upper right",frameon=False,fontsize=12)



for i in range(1,5):
    haxe.axvline(x=i,linestyle="--",linewidth=0.5,color="grey")




xmin,xmax=0,5
xticks=numpy.arange(1,5)
haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)


ymin,ymax=-15,10
yticks=numpy.arange(-15,10+5,5)
yticklabel=[r"$10^{%d}$"%(tmp) for tmp in yticks]
haxe.set_ylim(ymin,ymax)
haxe.set_yticks(yticks)
haxe.set_yticklabels(yticklabel)


haxe.tick_params(direction='in', length=5, width=1, colors='k')

haxe.tick_params(axis="x",labelsize=20)
haxe.tick_params(axis="y",labelsize=20)



# label
xlim_tmp=haxe.get_xlim()
ylim_tmp=haxe.get_ylim()
haxe.text(xlim_tmp[0]+0.05*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.85*(ylim_tmp[1]-ylim_tmp[0]),r"$(c2)$",color="k",fontsize=18)



haxe.set_xlabel(r"$\omega/\omega_F$",fontsize=20)

haxe.set_ylabel(r"$Power\ \mathcal{P}(\omega)$",fontsize=20)






lwid = 1.2




############################################################
############################################################
############################################################
############################################################
rid,cid = 0,1
axes_pos = [ 0.3767,  0.1600,  0.2667,  0.7800]

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[2,1],windows=windows,margin=margin)
axes_pos[0] -= 0.01

haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")


print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))



dir_tmp = list_dirname[0]
path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
path_data_tmp = os.path.join(path_pic_tmp,"data")

H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R),"r")
path_dataset_r="/psd/El_wkx"
data2d = H5FILE_R[path_dataset_r][()]
data2d = numpy.fliplr(data2d)

'''
print("data2d =",data2d.shape)
print("KKX    =",KKX.shape)

#sigma_w  = 0.05*wpe/norm_w
WTMP = KKX*norm_k*vd1*0.6/norm_w
DIS  = numpy.abs(KKW-WTMP)/sigma_w
MASK = numpy.exp(-0.5*DIS**2)

data2d*=MASK
'''

data2d = yaoxpy_vis.data_zero_replace(data2d) 
data2d = numpy.log10(data2d)


H5FILE_R.close()

print("data = %.4f - %.4f"%(numpy.min(data2d),numpy.max(data2d)))

him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto",vmin=vmin1,vmax=vmax1)

#him=haxe.imshow(data2d,origin="lower",extent=[numpy.min(ky),numpy.max(ky),numpy.min(kw),numpy.max(kw)],cmap=cmap,aspect="auto")


axes_pos_bar = [ 0.4333,  0.0450,  0.1733,  0.0150]
#axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[2,1],windows=windows,margin=margin,orient='horizontal',barbox=barbox)
haxe_bar = hfig.add_axes(axes_pos_bar)
hfig.colorbar(him,cax=haxe_bar,orientation='horizontal',shrink=0.5)

print("rid,cid       = (%d, %d)"%(rid,cid))
print("axes_pos_bar = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos_bar[0],axes_pos_bar[1],axes_pos_bar[2],axes_pos_bar[3]))






######################################### Bohm-Gross

haxe.plot(k/norm_k,wL_run3/norm_w,linestyle="-.",linewidth=lwid,color="r",label=r"$Langmuir$")

haxe.annotate(r"$Langmuir$",xy=(-7.5,1.1),xytext=(-18,0.8),arrowprops=dict(facecolor="r",edgecolor="r",width=0.2,headwidth=3.0,headlength=3.0),color="r")




######################################## Beam

#wBeam0 = vd0_run3*k
#haxe.plot(k/norm_k, wBeam0/norm_w,linestyle="--",linewidth=lwid,color="g")

wBeam1 = vd1_run3*k
haxe.plot(k/norm_k, wBeam1/norm_w,linestyle="--",linewidth=lwid,color="g")


haxe.annotate(r"$beam$",xy=(11.0,3.1),xytext=(13.5,2.9),arrowprops=dict(facecolor="g",edgecolor="g",width=0.2,headwidth=3.0,headlength=3.0),color="g",fontsize=16)


########################################

m = 2
wm_real = numpy.copy(wroots_real_run3[:,m])

index = wm_real>=0

haxe.plot(k[index]/norm_k,wm_real[index]/norm_w,linestyle="--",linewidth=lwid,color="m",label=r"$hybrid$")


haxe.annotate(r"$beam\textendash modified$",xy=(12,0.9),xytext=(13,0.6),arrowprops=dict(facecolor="m",edgecolor="m",width=0.2,headwidth=3.0,headlength=3.0),color="m",fontsize=16)






######################################## coalescence
k    = numpy.arange(0,30+0.01,0.01)*norm_k

lwid = 1.6




##########
w = k*vd1*0.5
haxe.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="w")

kL = list_kf[0]
wL = list_wf[0]

for n in [1,2,3,4]:
    haxe.scatter(n*kL,n*wL,marker="o",s=35,color="w")


##########
w = k*vd1*0.6
haxe.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="w")

kL = list_kf[1]
wL = list_wf[1]

for n in [1,2,3,4]:
    haxe.scatter(n*kL,n*wL,marker="o",s=35,color="w")


##########
w = k*vd1*0.7
haxe.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="w")

kL = list_kf[2]
wL = list_wf[2]

for n in [1,2,3,4]:
    haxe.scatter(n*kL,n*wL,marker="o",s=35,color="w")



##########
w = k*vd1*0.8
haxe.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="w")

kL = list_kf[3]
wL = list_wf[3]

for n in [1,2,3,4]:
    haxe.scatter(n*kL,n*wL,marker="o",s=35,color="w")



for n in [1,2,3,4]:
    haxe.text(n*kL-4,n*wL,r"$L_{%d}$"%(n),color="w",fontsize=24)






#haxe.legend(loc="upper center",frameon=True,fontsize=10)

#haxe.grid(linestyle="--",linewidth=0.2,color="grey")


xmin,xmax =-26,26
ymin,ymax =-0.2,4.2

xticks=numpy.arange(-25,25+5,5)
xticklabel=[r"$%d$"%(tmp) for tmp in xticks]

yticks=numpy.arange(0,4+1,1)
yticklabel=[r"$%d$"%(tmp) for tmp in yticks]




        
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
haxe.text(xlim_tmp[0]+0.05*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.9*(ylim_tmp[1]-ylim_tmp[0]),r"$(b)Wave\ Steepening\ nL\to L_n$",color="w",fontsize=18)




haxe.set_xlabel(r"$k_{\parallel}\cdot d_e$",fontsize=20)

haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=20)






############################################################
############################################################
############################################################
############################################################
rid,cid = 0,0
axes_pos = [ 0.0600,  0.1600,  0.2667,  0.7800]

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[2,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))



h5name_R1 = "fft2d_timestep_{:d}_period_{:d}dt_samples_{:d}_txy.h5".format(Timestep,FFT_Sample_Num_dt,FFT_Sample_Num)

dir_tmp = list_dirname[0]
path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
path_data_tmp = os.path.join(path_pic_tmp,"data")




H5FILE_R=h5py.File(os.path.join(path_data_tmp,h5name_R1),"r")
path_dataset_r="/data/Ex_tx"
data2d = H5FILE_R[path_dataset_r][()]
H5FILE_R.close()

data2d=data2d.T

t=numpy.arange(data2d.shape[1])*Dt
x=numpy.arange(data2d.shape[0])*Dx


norm_t = 1.0/wpe
norm_x = Dx

tmin0,tmax0 = 0,numpy.max(t)/norm_t
xmin0,xmax0 = 0,numpy.max(x)/norm_x

him=haxe.imshow(data2d,origin="lower",extent=[tmin0,tmax0,xmin0,xmax0],cmap=mpl.cm.coolwarm,aspect="auto")

axes_pos_bar = [ 0.1067,  0.0450,  0.1733,  0.0150]
#axes_pos_bar = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[2,1],windows=windows,margin=margin,orient='horizontal',barbox=barbox)
haxe_bar = hfig.add_axes(axes_pos_bar)
hfig.colorbar(him,cax=haxe_bar,orientation='horizontal',shrink=0.5)

print("rid,cid      = (%d, %d)"%(rid,cid))
print("axes_pos_bar = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos_bar[0],axes_pos_bar[1],axes_pos_bar[2],axes_pos_bar[3]))






#x_vd1   = t*vd1_run3
#x_vthe1 = t*vthe1

#haxe.plot(t/norm_t,x_vd1/norm_x,linestyle="-",linewidth=1.0,color="g")
#haxe.plot(t/norm_t,x_vthe1/norm_x,linestyle="-",linewidth=1.0,color="r")

vd1   = vd1_run3


lwid = 1.6


#####
x_vd1 = 0.5*t*vd1-330*norm_x

index = numpy.logical_and(x_vd1>=0, x_vd1<=1024*norm_x)

haxe.plot(t[index]/norm_t,x_vd1[index]/norm_x,linestyle="--",linewidth=lwid,color="w")

haxe.text(130,25,r"$v_{ph}=0.5v_{d2}$",color="w",fontsize=16,rotation=25)



#####
x_vd1 = 0.6*t*vd1-214*norm_x


index = numpy.logical_and(x_vd1>=0, x_vd1<=1024*norm_x)

haxe.plot(t[index]/norm_t,x_vd1[index]/norm_x,linestyle="--",linewidth=lwid,color="w")

haxe.text(100,130,r"$v_{ph}=0.6v_{d2}$",color="w",fontsize=16,rotation=25)




#####
x_vd1 = 0.7*t*vd1+635*norm_x

index = numpy.logical_and(x_vd1>=0, x_vd1<=1024*norm_x)

haxe.plot(t[index]/norm_t,x_vd1[index]/norm_x,linestyle="--",linewidth=lwid,color="w")

haxe.text(45,750,r"$v_{ph}=0.7v_{d2}$",color="w",fontsize=16,rotation=30)




#####
x_vd1 = 0.80*t*vd1+822*norm_x

index = numpy.logical_and(x_vd1>=0, x_vd1<=1024*norm_x)

haxe.plot(t[index]/norm_t,x_vd1[index]/norm_x,linestyle="--",linewidth=lwid,color="w")

haxe.text(10,815,r"$v_{ph}=0.8v_{d2}$",color="w",fontsize=16,rotation=35)












#haxe.legend(loc="upper center",frameon=True,fontsize=10)

#haxe.grid(linestyle="--",linewidth=0.2,color="grey")


xmin,xmax = tmin0,165
ymin,ymax = xmin0,xmax0

xticks=numpy.arange(0,160+40,40)
xticklabel=[r"$%d$"%(tmp) for tmp in xticks]

yticks=numpy.arange(0,1000+200,200)
yticklabel=[r"$%d$"%(tmp) for tmp in yticks]

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
haxe.text(xlim_tmp[0]+0.05*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.9*(ylim_tmp[1]-ylim_tmp[0]),r"$(a)\overline{E_{l}(x,t)}$",color="k",fontsize=18)





haxe.set_xlabel(r"$t\cdot\omega_{pe}$",fontsize=20)

haxe.set_ylabel(r"$x/\Delta x$",fontsize=20)







############################################################
#plt.show()
    
fig_path = path_fig
fig_name = os.path.splitext(os.path.basename(__file__))[0]
    
yaoxpy_vis.fig_save(plt,figpath=fig_path,figname=fig_name,extension=fig_fmts)