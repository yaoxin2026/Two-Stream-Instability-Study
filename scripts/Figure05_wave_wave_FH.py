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

xmin,xmax =-26,26
ymin,ymax =-0.2,4.2

xticks=numpy.arange(-25,25+5,5)
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

quiver_y2x = 11.8



hfig = plt.figure(figsize=(16,11))

margin=[0.08,0.08,0.04,0.08,0.06,0.08]
barbox=[0.03,0.015,0.65]
windows=[3,2]
size=[1,1]


############################################################
############################################################
rid,cid = 1,0
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[2,1],windows=windows,margin=margin)
axes_pos = [0.08,0.08,0.39,0.56]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")



dir_tmp = list_dirname[0]
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

#haxe.axhline(y=0.9,linestyle='--',linewidth=1.0,color="w")
#haxe.axhline(y=1.8,linestyle='--',linewidth=1.0,color="w")





######################################### Bohm-Gross
k = numpy.arange(-30,30+0.01,0.01)*norm_k
#w = numpy.sqrt(wpe0*wpe0+3.0*numpy.power(k*vthe0,2.0))

w = numpy.sqrt(wpeL*wpeL+3.0*numpy.power(k*vthe0,2.0))

haxe.plot(k/norm_k,w/norm_w,linestyle="-.",linewidth=1.2,color="r")

haxe.annotate(r"$Langmuir$",xy=(-7.0,1.0),xytext=(-15.5,0.7),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")

haxe.axhline(y=1.0,linestyle='--',linewidth=0.8,color="w")
haxe.axhline(y=2.0,linestyle='--',linewidth=0.8,color="w")


########################################

for i in [0]:
    wtmp=numpy.real(w_run3[i,:])
    haxe.plot(k/norm_k,wtmp/norm_w,linestyle="--",linewidth=1.0,color=list_color12[i])



for i in [4]:
    wtmp=numpy.real(w_run3[i,:])

    index=wtmp/norm_w>=-100
    haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="--",linewidth=1.0,color="b")

    xtmp=k[index]/norm_k
    ytmp=wtmp[index]/norm_w

    spline = scipy.interpolate.UnivariateSpline(xtmp,ytmp,s=5.0)
    xtmp = numpy.linspace(0.0,numpy.max(xtmp),400)
    ytmp = spline(xtmp)

    haxe.plot(xtmp,ytmp,linestyle="-.",linewidth=1.0,color="b")




haxe.annotate(r"$Beam-Beam$",xy=(-10.0,1.5),xytext=(-8.5,1.73),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")


haxe.annotate(r"$hybrid$",xy=(8.0,1.0),xytext=(9.0,0.65),arrowprops=dict(facecolor="b",edgecolor="b",width=0.4,headwidth=4.0,headlength=4.0),color="b")




######################################## IA
w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0,alpha1,vthe0,vthe1,k)
haxe.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=1.0,color="b")



index=k/norm_k>5.0
print("IA : wl = %.4fwpe"%(w[index][0]/norm_w))

#haxe.axhline(y=1.0,linestyle="--",linewidth=1.0,color="w")
#haxe.axhline(y=0.9,linestyle="--",linewidth=1.0,color="w")
#haxe.axhline(y=0.9*2.0,linestyle="--",linewidth=1.0,color="w")







######################################## beam mode
#k2   = numpy.arange(0,30+0.01,0.01)*norm_k
#w = vb*k2
#haxe.plot(k2/norm_k,w/norm_w,linestyle="-",linewidth=2.0,color="g",label=r"$beam\ \omega=v_{b}\cdot k$")








######################################## coalescence

lx0_L=6.5
haxe.quiver(0.0,0.0,lx0_L,quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="r")

haxe.text(3.5,0.25,r"$L$",color="r",fontsize=28)



#####

lx0_S = 9.5
haxe.quiver(0.0,0.0,lx0_S,0.0,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="b")

print("S = ",lx0_S*cs/wpe0)


haxe.text(11.0,-0.1,r"$S$",color="b",fontsize=28)


haxe.quiver(lx0_S,0.0,(lx0_L-lx0_S),quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="grey")

haxe.quiver(0.0,0.0,(lx0_L-lx0_S),quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="w")

haxe.text(-5.0,0.5,r"$F$",color="w",fontsize=28)






########################################
w=-5*numpy.ones(len(k))

haxe.plot(k/norm_k,w,linestyle="--",linewidth=1.0,label=r"$Beam\textendash Beam$",color="r")

haxe.plot(k/norm_k,w,linestyle="-.",linewidth=1.0,label=r"$Langmuir$",color="r")

#haxe.plot(k/norm_k,w,linestyle="--",linewidth=1.0,label=r"$Beam$",color="g")




haxe.plot(k/norm_k,w,linestyle="--",linewidth=1.0,label=r"$IA$",color="b")


haxe.plot(k/norm_k,w,linestyle="--",linewidth=1.0,label=r"$Beam\textendash modified\ Langmuir\ (hybrid)$",color="cyan")

haxe.plot(k/norm_k,w,linestyle="--",linewidth=1.0,label=r"$EA$",color="k")



r'''
haxe.plot(k/norm_k,w,linestyle="--",linewidth=2.0,label=r"$ \omega=k\cdot v_b$",color="g")

haxe.plot(k/norm_k,w,linestyle="-",linewidth=2.0,label=r"$Beam$ - $like$",color="coral")

haxe.plot(k/norm_k,w,linestyle="--",linewidth=2.0,label=r"$IA$",color="grey")

#haxe.plot(k/norm_k,w,linestyle="--",linewidth=2.0,label=r"$\omega=k\cdot c_s$",color="b")
'''


#haxe.legend(loc="upper left",bbox_to_anchor=(-0.04,1.35),ncol=3,frameon=False,prop={"size":17})

#haxe.legend(loc="upper left",bbox_to_anchor=(1.11,1.05),ncol=1,frameon=False,prop={"size":18})









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
haxe.text(xlim_tmp[0]+0.05*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.9*(ylim_tmp[1]-ylim_tmp[0]),r"$(b1)Decay\ Process\ L\to S+F$",color="w",fontsize=24)


haxe.set_xlabel(r"$k_{\parallel}\cdot d_e$",fontsize=24)


haxe.set_ylabel(r"$\omega/\omega_{pe}^L$",fontsize=24)






















############################################################
############################################################
rid,cid = 1,1
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[2,1],windows=windows,margin=margin)
axes_pos = [0.53,0.08,0.39,0.56]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")




dir_tmp = list_dirname[0]
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


#axes_pos = yaoxpy_vis.fig_axes_position_colorbar(index=[rid,cid],size=[2,1],windows=windows,margin=margin,orient='vertical',barbox=barbox)
axes_pos = [0.95,0.178,0.015,0.364]
haxe_bar = hfig.add_axes(axes_pos)
hfig.colorbar(him,cax=haxe_bar,orientation='vertical',shrink=0.5)









######################################### Bohm-Gross
k = numpy.arange(-30,30+0.01,0.01)*norm_k
#w = numpy.sqrt(wpe0*wpe0+3.0*numpy.power(k*vthe0,2.0))

w = numpy.sqrt(wpeL*wpeL+3.0*numpy.power(k*vthe0,2.0))

haxe.plot(k/norm_k,w/norm_w,linestyle="-.",linewidth=1.2,color="r")


haxe.annotate(r"$Langmuir$",xy=(-7.0,1.0),xytext=(-15.5,0.7),arrowprops=dict(facecolor="r",edgecolor="r",width=0.4,headwidth=4.0,headlength=4.0),color="r")



haxe.axhline(y=1.0,linestyle='--',linewidth=0.8,color="w")
haxe.axhline(y=2.0,linestyle='--',linewidth=0.8,color="w")




########################################

'''
for i in [0]:
    wtmp=numpy.real(w_run3[i,:])
    haxe.plot(k/norm_k,wtmp/norm_w,linestyle="--",linewidth=1.0,color=list_color12[i])


for i in [1,10]:
    wtmp=numpy.real(w_run3[i,:])

    index=wtmp/norm_w>=0.0
    haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="-",linewidth=1.0,color=list_color12[i])

    index=numpy.logical_and(wtmp/norm_w<=0.0,wtmp/norm_w>=-0.5)
    haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="--",linewidth=1.0,color=list_color12[i])
'''

'''
for i in [4]:
    wtmp=numpy.real(w_run3[i,:])

    index=wtmp/norm_w>=0.0
    haxe.plot(k[index]/norm_k,wtmp[index]/norm_w,linestyle="--",linewidth=1.0,color=list_color12[i])
'''


######################################## IA
#w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0,alpha1,vthe0,vthe1,k)
#haxe.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=1.0,color="b",label=r"$IA$")





######################################## beam mode
#k2   = numpy.arange(0,30+0.01,0.01)*norm_k
#w = vb*k2
#haxe.plot(k2/norm_k,w/norm_w,linestyle="-",linewidth=2.0,color="g",label=r"$beam\ \omega=v_{b}\cdot k$")




######################################## coalescence


lx0_L=0.0
haxe.quiver(0.0,0.0,lx0_L,quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="r")

haxe.text(2.5,0.3,r"$L$",color="r",fontsize=28)

print("L =",1.0*wpe/wpe0)


lx0_F=-3.0
ly0_F=1.0

haxe.quiver(0.0,0.0,lx0_F,ly0_F*quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="w")

#haxe.text(-6.5,1.3,r"$L-S\to F$",color="w",fontsize=24)



haxe.quiver(lx0_F,ly0_F,lx0_L,quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="grey")


haxe.quiver(0.0,0.0,lx0_F,(ly0_F+1.0)*quiver_y2x,scale=1.0,scale_units="x",width=0.01,linewidth=1.2,color="w")



#haxe.text(-6.5,2.3,r"$F+L\to H$",color="w",fontsize=24)



haxe.text(-5,0.5,r"$F$",color="w",fontsize=28)


haxe.text(-3.5,2.1,r"$H$",color="w",fontsize=28)








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
haxe.text(xlim_tmp[0]+0.05*(xlim_tmp[1]-xlim_tmp[0]),ylim_tmp[0]+0.9*(ylim_tmp[1]-ylim_tmp[0]),r"$(b2)Coalescence\ Process\ F+L\to H$",color="w",fontsize=24)



haxe.set_xlabel(r"$k_{\parallel}\cdot d_e$",fontsize=24)


#haxe.set_ylabel(r"$\omega/\omega_{pe}^L$",fontsize=24)





############################################################
############################################################
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
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
axes_pos = [0.08,0.72,0.39,0.24]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")



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




'''
index=ytmp==numpy.max(ytmp)
xtmp_max=xtmp[index][0]
haxe.axvline(x=xtmp_max,linestyle="--",linewidth=0.8,color="g")
haxe.text(xtmp_max-6,3.0,r"$t=%.2f$"%(xtmp_max),rotation=90,color="g")

print("beam = ",xtmp_max)

'''






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
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
axes_pos = [0.53,0.72,0.39,0.24]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")



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
    
yaoxpy.fig_save(plt,figpath=fig_path,figname=fig_name,extension=fig_fmts)