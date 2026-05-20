from YaoxPy_Import_CWD import *

############################################################
############################################################

list_dirname = ["yaoxpic_v25_counter_1","yaoxpic_v25_counter_2","yaoxpic_v25_counter_3"]


list_Timestep_Particle = [0,600,800,1000,1200,1600,2000,3000,4000,6000,8000,10000]

list_Timestep_Particle_ion = [10000]


list_color=yaoxpy.colors_generate(len(list_Timestep_Particle)+1,mpl.cm.rainbow)

list_color = list_color[::-1]

list_color2= ["grey"]

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
vd2    = list_parameters["vd2"]
vthe2  = list_parameters["vth2"]
    
Dx    = dx
#Dt   = dt*FFT_Sample_Num_dt
    
cs    = CGS["c"]
    
de    = CGS["c"]/wpe
rhon0 = (wpe/CGS["e"])**2*CGS["me"]/4.0/numpy.pi
J0    = CGS["e"]*rhon0*vthe0*CGS["c"]
B0    = CGS["me"]*CGS["c"]/CGS["e"]*wce
    
print("*"*20)


############################################################
############################################################

hfig = plt.figure(figsize=(15,11))

margin=[0.06,0.04,0.03,0.08,0.1,0.05]
barbox=[0.05,0.01,0.8]
windows=[3,2]
size=[1,1]


############################################################
############################################################
rid,cid=0,0
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
#print("axes_pos =",axes_pos)

axes_pos = [0.06,0.70666667,0.4,0.26333333]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = %d,%d"%(rid,cid))
print("axes_pos =",axes_pos)


for tid in range(len(list_Timestep_Particle)):
    
    Timestep=list_Timestep_Particle[tid]

    print("Timestep = %d"%(Timestep))

    dir_tmp = list_dirname[rid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    filename = "timestep_%d_particle_vdf.h5"%(Timestep)
    
    H5FILE_R=h5py.File(os.path.join(path_data_tmp,filename),"r")


    sid=0
    bin_vpara = H5FILE_R["/bin/vpara"][()]
    bin_vperp = H5FILE_R["/bin/vperp"][()]


    vpara0 = H5FILE_R["/species_%d/vpara/vx"%(sid)][()]

    sid=2
    vpara2 = H5FILE_R["/species_%d/vpara/vx"%(sid)][()]


    H5FILE_R.close()


    nop1=400
    nop2=400

    alpha1=nop1/(nop1+nop2)
    alpha2=nop2/(nop1+nop2)

    vpara0 = vpara0*alpha1
    vpara2 = vpara2*alpha2

    #haxe.plot(bin_vpara,vpara0+vpara2,linestyle="-",linewidth=2.0,color="r",label=r"$f_{e}(v_{\parallel})=f_{e0}(v_{\parallel})+f_{eb}(v_{\parallel})$")

    #haxe.plot(bin_vpara,vpara2,linestyle="-",linewidth=2.0,color=list_color[tid],label=r"$t=%.2f\omega_{pe}^{-1}$"%(Timestep*dt*wpe))

    #haxe.plot(bin_vpara,vpara0+vpara2,linestyle="-",linewidth=2.0,color=list_color[tid],label=r"$timestep=%d$"%(Timestep))

    haxe.plot(bin_vpara,vpara0+vpara2,linestyle="-",linewidth=2.5,color=list_color[tid],label=r"$t\cdot \omega_{pe}=%.0f$"%(Timestep*dt*wpe))



r'''
for tid in range(len(list_Timestep_Particle_ion)):
    
    Timestep=list_Timestep_Particle_ion[tid]

    #Timestep=0

    dir_tmp = list_dirname[0]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    filename = "timestep_%d_particle_vdf_ion.h5"%(Timestep)
    
    H5FILE_R=h5py.File(os.path.join(path_data_tmp,filename),"r")



    bin_vpara = H5FILE_R["/bin/vpara"][()]
    bin_vperp = H5FILE_R["/bin/vperp"][()]

    sid=1
    vpara1 = H5FILE_R["/species_%d/vpara/vx"%(sid)][()]

    sid=3
    vpara3 = H5FILE_R["/species_%d/vpara/vx"%(sid)][()]


    H5FILE_R.close()

    print("vpara =",numpy.max(vpara1),numpy.max(vpara3))


    vpara1 = vpara1*alpha1
    vpara3 = vpara3*alpha2


    haxe.plot(bin_vpara,(vpara1+vpara3)*0.025,linestyle="-",linewidth=2.5,color=list_color2[tid],label=r"$t\cdot \omega_{pe}=%.0f\ (ion)$"%(Timestep*dt*wpe))
'''


haxe.grid(True,linestyle="--",linewidth=1,color="w")

vb1=-0.3
vb2=0.3

#haxe.axvline(x=0.0,linestyle="--",linewidth=1.0,color=YXColorBlue)
haxe.axvline(x=vb1,linestyle="--",linewidth=1.0,color="r")
haxe.axvline(x=vb2,linestyle="--",linewidth=1.0,color="r")


haxe.legend(loc="upper right",frameon=False,ncol=2,fontsize=11)

#haxe.legend(loc="upper right",bbox_to_anchor=(1.45,1.1),ncol=1,frameon=False,prop={"size":16})


xmin,xmax=-0.43,0.43
haxe.set_xlim(xmin,xmax)


ymin,ymax=-1,8
haxe.set_ylim(ymin,ymax)
haxe.set_yticks(numpy.arange(0,8+2,2))


xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.04*(xlim[1]-xlim[0]),ylim[0]+0.85*(ylim[1]-ylim[0]),r"$(a1)Run1$",fontsize=20)


#haxe.set_xlabel(r"$v_{\parallel}/c$",fontsize=24)
haxe.set_ylabel(r"$f(v_{\parallel})$",fontsize=24)




############################################################
############################################################
rid,cid=1,0
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
#print("axes_pos =",axes_pos)

axes_pos = [0.06,0.39333333,0.4,0.26333333]

haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = %d,%d"%(rid,cid))
print("axes_pos =",axes_pos)

for tid in range(len(list_Timestep_Particle)):
    
    Timestep=list_Timestep_Particle[tid]

    #Timestep=0

    dir_tmp = list_dirname[rid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    filename = "timestep_%d_particle_vdf.h5"%(Timestep)
    
    H5FILE_R=h5py.File(os.path.join(path_data_tmp,filename),"r")


    sid=0
    bin_vpara = H5FILE_R["/bin/vpara"][()]
    bin_vperp = H5FILE_R["/bin/vperp"][()]


    vpara0 = H5FILE_R["/species_%d/vpara/vx"%(sid)][()]

    sid=2
    vpara2 = H5FILE_R["/species_%d/vpara/vx"%(sid)][()]


    H5FILE_R.close()


    nop1=400
    nop2=100

    alpha1=nop1/(nop1+nop2)
    alpha2=nop2/(nop1+nop2)

    vpara0 = vpara0*alpha1
    vpara2 = vpara2*alpha2

    #haxe.plot(bin_vpara,vpara0+vpara2,linestyle="-",linewidth=2.0,color="r",label=r"$f_{e}(v_{\parallel})=f_{e0}(v_{\parallel})+f_{eb}(v_{\parallel})$")

    #haxe.plot(bin_vpara,vpara2,linestyle="-",linewidth=2.0,color=list_color[tid],label=r"$t=%.2f\omega_{pe}^{-1}$"%(Timestep*dt*wpe))

    #haxe.plot(bin_vpara,vpara0+vpara2,linestyle="-",linewidth=2.0,color=list_color[tid],label=r"$timestep=%d$"%(Timestep))

    haxe.plot(bin_vpara,vpara0+vpara2,linestyle="-",linewidth=2.5,color=list_color[tid],label=r"$t\cdot \omega_{pe}=%.0f$"%(Timestep*dt*wpe))



r'''
for tid in range(len(list_Timestep_Particle_ion)):
    
    Timestep=list_Timestep_Particle_ion[tid]

    #Timestep=0

    dir_tmp = list_dirname[0]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    filename = "timestep_%d_particle_vdf_ion.h5"%(Timestep)
    
    H5FILE_R=h5py.File(os.path.join(path_data_tmp,filename),"r")



    bin_vpara = H5FILE_R["/bin/vpara"][()]
    bin_vperp = H5FILE_R["/bin/vperp"][()]

    sid=1
    vpara1 = H5FILE_R["/species_%d/vpara/vx"%(sid)][()]

    sid=3
    vpara3 = H5FILE_R["/species_%d/vpara/vx"%(sid)][()]


    H5FILE_R.close()

    print("vpara =",numpy.max(vpara1),numpy.max(vpara3))



    vpara1 = vpara1*alpha1
    vpara3 = vpara3*alpha2

    haxe.plot(bin_vpara,(vpara1+vpara3)*0.025,linestyle="-",linewidth=2.5,color=list_color2[tid],label=r"$t\cdot \omega_{pe}=%.0f\ (ion)$"%(Timestep*dt*wpe))
'''


haxe.grid(True,linestyle="--",linewidth=1,color="w")

vb1=-0.3/4.0
vb2=0.3

#haxe.axvline(x=0.0,linestyle="--",linewidth=1.0,color=YXColorBlue)
haxe.axvline(x=vb1,linestyle="--",linewidth=1.0,color="r")
haxe.axvline(x=vb2,linestyle="--",linewidth=1.0,color="r")

#haxe.legend(loc="upper right",frameon=False)


xmin,xmax=-0.43,0.43
haxe.set_xlim(xmin,xmax)


ymin,ymax=-1,12
haxe.set_ylim(ymin,ymax)
haxe.set_yticks(numpy.arange(0,12+3,3))


xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.04*(xlim[1]-xlim[0]),ylim[0]+0.85*(ylim[1]-ylim[0]),r"$(b1)Run2$",fontsize=20)


#haxe.set_xlabel(r"$v_{\parallel}/c$",fontsize=24)
haxe.set_ylabel(r"$f(v_{\parallel})$",fontsize=24)




############################################################
############################################################
rid,cid=2,0
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
#print("axes_pos =",axes_pos)

axes_pos = [0.06,0.08,0.4,0.26333333]

haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = %d,%d"%(rid,cid))
print("axes_pos =",axes_pos)

for tid in range(len(list_Timestep_Particle)):
    
    Timestep=list_Timestep_Particle[tid]

    #Timestep=0

    dir_tmp = list_dirname[rid]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    filename = "timestep_%d_particle_vdf.h5"%(Timestep)
    
    H5FILE_R=h5py.File(os.path.join(path_data_tmp,filename),"r")


    sid=0
    bin_vpara = H5FILE_R["/bin/vpara"][()]
    bin_vperp = H5FILE_R["/bin/vperp"][()]


    vpara0 = H5FILE_R["/species_%d/vpara/vx"%(sid)][()]

    sid=2
    vpara2 = H5FILE_R["/species_%d/vpara/vx"%(sid)][()]


    H5FILE_R.close()


    nop1=400
    nop2=20

    alpha1=nop1/(nop1+nop2)
    alpha2=nop2/(nop1+nop2)

    vpara0 = vpara0*alpha1
    vpara2 = vpara2*alpha2*2.0

    #haxe.plot(bin_vpara,vpara0+vpara2,linestyle="-",linewidth=2.0,color="r",label=r"$f_{e}(v_{\parallel})=f_{e0}(v_{\parallel})+f_{eb}(v_{\parallel})$")

    #haxe.plot(bin_vpara,vpara2,linestyle="-",linewidth=2.0,color=list_color[tid],label=r"$t=%.2f\omega_{pe}^{-1}$"%(Timestep*dt*wpe))

    #haxe.plot(bin_vpara,vpara0+vpara2,linestyle="-",linewidth=2.0,color=list_color[tid],label=r"$timestep=%d$"%(Timestep))

    haxe.plot(bin_vpara,vpara0+vpara2,linestyle="-",linewidth=2.5,color=list_color[tid],label=r"$t\cdot \omega_{pe}=%.0f$"%(Timestep*dt*wpe))



r'''
for tid in range(len(list_Timestep_Particle_ion)):
    
    Timestep=list_Timestep_Particle_ion[tid]

    #Timestep=0

    dir_tmp = list_dirname[0]
    path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
    path_data_tmp = os.path.join(path_pic_tmp,"data")

    filename = "timestep_%d_particle_vdf_ion.h5"%(Timestep)
    
    H5FILE_R=h5py.File(os.path.join(path_data_tmp,filename),"r")



    bin_vpara = H5FILE_R["/bin/vpara"][()]
    bin_vperp = H5FILE_R["/bin/vperp"][()]

    sid=1
    vpara1 = H5FILE_R["/species_%d/vpara/vx"%(sid)][()]

    sid=3
    vpara3 = H5FILE_R["/species_%d/vpara/vx"%(sid)][()]


    H5FILE_R.close()

    print("vpara =",numpy.max(vpara1),numpy.max(vpara3))


    vpara1 = vpara1*alpha1
    vpara3 = vpara3*alpha2

    haxe.plot(bin_vpara,(vpara1+vpara3)*0.025,linestyle="-",linewidth=2.5,color=list_color2[tid],label=r"$t\cdot \omega_{pe}=%.0f\ (ion)$"%(Timestep*dt*wpe))
'''


haxe.grid(True,linestyle="--",linewidth=1,color="w")

vb1=-0.3/20.0
vb2=0.3

#haxe.axvline(x=0.0,linestyle="--",linewidth=1.0,color=YXColorBlue)
haxe.axvline(x=vb1,linestyle="--",linewidth=1.0,color="r")
haxe.axvline(x=vb2,linestyle="--",linewidth=1.0,color="r")


#haxe.legend(loc="upper right",frameon=False)


xmin,xmax=-0.43,0.43
haxe.set_xlim(xmin,xmax)


ymin,ymax=-1,15
haxe.set_ylim(ymin,ymax)
haxe.set_yticks(numpy.arange(0,15+3,3))


xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.04*(xlim[1]-xlim[0]),ylim[0]+0.85*(ylim[1]-ylim[0]),r"$(c1)Run3$",fontsize=20)


haxe.set_xlabel(r"$v_{\parallel}/c$",fontsize=24)
haxe.set_ylabel(r"$f(v_{\parallel})$",fontsize=24)



############################################################
############################################################
############################################################
############################################################
rid,cid=0,1
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)

axes_pos = [0.56,0.70666667,0.4,0.26333333]

haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = %d,%d"%(rid,cid))
print("axes_pos =",axes_pos)

dir_tmp = list_dirname[rid]
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

time = numpy.arange(EE.shape[0])*wpe*dt


index=timestep_state>0
time=time[index]
EE=EE[index,:]
EB=EB[index,:]
EK=EK[index,:]

print("EE =",EE.shape)
print("EB =",EB.shape)
print("EK =",EK.shape)



norm_E = EK[0,2]

print("norm_E =",norm_E)


print("*"*40+" Run1")


####################

NUM_INTP = 800


##### E+B
DEB = EE[:,0]+EE[:,1]+EE[:,2]+EB[:,0]+EB[:,1]+EB[:,2]
DEB = DEB-DEB[0]

xtmp=time
ytmp=DEB/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

#haxe.plot(xtmp,ytmp,linestyle="--",linewidth=2.0,color="darkblue")

### vertical t0
t0=[]
for i in range(1,len(ytmp)-1):
    if ytmp[i]>=ytmp[i-1] and ytmp[i]>=ytmp[i+1]:
       print("t0 = %.4f"%(xtmp[i]))
       t0.append(xtmp[i])
    
    if len(t0)>0:
       break

t0=numpy.array(t0)

for i in range(1):
    haxe.axvline(x=t0[i],linestyle="--",linewidth=1.0,color="grey")
    haxe.text(t0[i]-10,0.05,r"$t_%d=%.2f$"%(i,t0[i]),color="grey",rotation=90)

for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, E+B = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))


### scatter
xtmp=numpy.linspace(0,numpy.max(xtmp),50)
ytmp=fun_intp(xtmp)

index=numpy.logical_or(xtmp<20,xtmp>35)
#haxe.scatter(xtmp[index],ytmp[index],marker="o",s=25,color="darkblue",label=r"$\Delta \mathcal{E}_{E}+\Delta \mathcal{E}_{B}$")
haxe.scatter(xtmp[index],ytmp[index],marker="o",s=25,color="darkblue")

xtmp=numpy.array([21,21.8,22.4,23,23.7,25,27.5,28.6,30,32,35])
ytmp=fun_intp(xtmp)
haxe.scatter(xtmp,ytmp,marker="o",s=25,color="darkblue")



##### Ex
DEE = EE[:,0]
DEE = DEE-DEE[0]

xtmp=time
ytmp=DEE/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="m",label=r"$\Delta \mathcal{E}_{E_{x}}$")

for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, Ex  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))


##### Ey
DEE = EE[:,1]
DEE = DEE-DEE[0]

xtmp=time
ytmp=DEE/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="b",label=r"$\Delta \mathcal{E}_{E_{y}}$")


for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, Ey  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))




##### Bz
DEB = EB[:,2]
DEB = DEB-DEB[0]

xtmp=time
ytmp=DEB/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="grey",label=r"$\Delta \mathcal{E}_{B_{z}}$")

for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, Bz  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))



haxe.scatter(numpy.arange(10),numpy.zeros(10)-10,marker="o",s=25,color="darkblue",label=r"$\Delta \mathcal{E}_{E}+\Delta \mathcal{E}_{B}$")

#haxe.scatter(numpy.arange(10),numpy.zeros(10)-10,marker="o",s=25,color="None",label=r"$\ $")



####################


##### Ek0
DEK = EK[:,0]
DEK = DEK-DEK[0]

xtmp=time
ytmp=DEK/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="g",label=r"$\Delta \mathcal{E}_{k1}^{e}$")

for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, K0  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))




##### Ek2
DEK = EK[:,2]
DEK = DEK-DEK[0]

xtmp=time
ytmp=-1.0*DEK/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="r",label=r"$-\Delta \mathcal{E}_{k2}^{e}$")


for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, K2  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))



##### Ek-all
DEK = EK[:,0]+EK[:,1]+EK[:,2]+EK[:,3]
DEK = DEK-DEK[0]

xtmp=time
ytmp=-1.0*DEK/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="--",linewidth=2.0,color="r",label=r"$-\Delta \mathcal{E}_{k}$")


for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, K-all = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))



haxe.legend(loc="upper right",frameon=False,ncol=2,fontsize=13)

haxe.grid(linestyle="--",linewidth=0.5,color="w")



xmin,xmax =-5,165
xticks=numpy.arange(0,160+20,20)
xticklabel=[r"$%d$"%(tmp) for tmp in xticks]

ymin,ymax =-0.1,0.2
yticks=numpy.arange(-0.1,0.2+0.05,0.05)
print(yticks)
yticklabel=[r"$%.2f$"%(tmp) for tmp in yticks]

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)
haxe.set_xticklabels(xticklabel)


haxe.set_ylim(ymin,ymax)
haxe.set_yticks(yticks)
#haxe.set_yticklabels(yticklabel)



haxe.tick_params(direction='in', length=5, width=1, colors='k')

haxe.tick_params(axis="x",labelsize=18)
haxe.tick_params(axis="y",labelsize=18)

# label
xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.04*(xlim[1]-xlim[0]),ylim[0]+0.85*(ylim[1]-ylim[0]),r"$(a2)$",fontsize=20)



#haxe.set_xlabel(r"$t\cdot \omega_{pe}$",fontsize=24)

haxe.set_ylabel(r"$\Delta \mathcal{E}/\mathcal{E}_{k0}$",fontsize=24)










############################################################
############################################################
rid,cid=1,1
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)

axes_pos = [0.56,0.39333333,0.4,0.26333333]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = %d,%d"%(rid,cid))
print("axes_pos =",axes_pos)


dir_tmp = list_dirname[rid]
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

time = numpy.arange(EE.shape[0])*wpe*dt


index=timestep_state>0
time=time[index]
EE=EE[index,:]
EB=EB[index,:]
EK=EK[index,:]

print("EE =",EE.shape)
print("EB =",EB.shape)
print("EK =",EK.shape)



norm_E = EK[0,2]

print("norm_E =",norm_E)


print("*"*40+" Run2")


####################

NUM_INTP = 800


##### E+B
DEB = EE[:,0]+EE[:,1]+EE[:,2]+EB[:,0]+EB[:,1]+EB[:,2]
DEB = DEB-DEB[0]

xtmp=time
ytmp=DEB/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

#haxe.plot(xtmp,ytmp,linestyle="--",linewidth=2.0,color="darkblue")



### vertical t0
t0=[]
for i in range(1,len(ytmp)-1):
    if (ytmp[i]>ytmp[i-1] and ytmp[i]>ytmp[i+1]) and xtmp[i]>25:
       print("t0 = %.4f"%(xtmp[i]))
       t0.append(xtmp[i])
    
    if len(t0)>0:
       break

t0=numpy.array(t0)

for i in range(1):
    haxe.axvline(x=t0[i],linestyle="--",linewidth=1.0,color="grey")
    haxe.text(t0[i]-10,0.2,r"$t_%d=%.2f$"%(i,t0[i]),color="grey",rotation=90)

for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, E+B = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))


### scatter
xtmp=numpy.linspace(0,numpy.max(xtmp),50)
ytmp=fun_intp(xtmp)

index=numpy.logical_or(xtmp<20,xtmp>35)
#haxe.scatter(xtmp[index],ytmp[index],marker="o",s=25,color="darkblue",label=r"$\Delta \mathcal{E}_{E}+\Delta \mathcal{E}_{B}$")
haxe.scatter(xtmp[index],ytmp[index],marker="o",s=25,color="darkblue")

xtmp=numpy.array([21,21.8,22.4,23,23.7,25,27.5,28.6,30,32,35])
ytmp=fun_intp(xtmp)
haxe.scatter(xtmp,ytmp,marker="o",s=25,color="darkblue")



##### Ex
DEE = EE[:,0]
DEE = DEE-DEE[0]

xtmp=time
ytmp=DEE/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="m",label=r"$\Delta \mathcal{E}_{E_{x}}$")

for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, Ex  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))


##### Ey
DEE = EE[:,1]
DEE = DEE-DEE[0]

xtmp=time
ytmp=DEE/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="b",label=r"$\Delta \mathcal{E}_{E_{y}}$")


for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, Ey  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))




##### Bz
DEB = EB[:,2]
DEB = DEB-DEB[0]

xtmp=time
ytmp=DEB/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="grey",label=r"$\Delta \mathcal{E}_{B_{z}}$")

for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, Bz  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))



haxe.scatter(numpy.arange(10),numpy.zeros(10)-10,marker="o",s=25,color="darkblue",label=r"$\Delta \mathcal{E}_{E}+\Delta \mathcal{E}_{B}$")

#haxe.scatter(numpy.arange(10),numpy.zeros(10)-10,marker="o",s=25,color="None",label=r"$\ $")



####################


##### Ek0
DEK = EK[:,0]
DEK = DEK-DEK[0]

xtmp=time
ytmp=DEK/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="g",label=r"$\Delta \mathcal{E}_{k1}^{e}$")

for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, K0  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))



### vertical t1
t1=[]
for i in range(1,len(ytmp)-1):
    if (ytmp[i]>ytmp[i-1] and ytmp[i]>ytmp[i+1]) and xtmp[i]>25:
       #print("t1 = %.4f"%(xtmp[i]))
       t1.append(xtmp[i])
    
    if len(t1)>0:
       break

t1=numpy.array(t1)

for t1_tmp in t1:
    index=xtmp>=t1_tmp
    print("t1 = %6.2f, K0 = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))







##### Ek2
DEK = EK[:,2]
DEK = DEK-DEK[0]

xtmp=time
ytmp=-1.0*DEK/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="r",label=r"$-\Delta \mathcal{E}_{k2}^{e}$")


for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, K2  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))



### vertical t1
t1=[]
for i in range(1,len(ytmp)-1):
    if (ytmp[i]>ytmp[i-1] and ytmp[i]>ytmp[i+1]) and xtmp[i]>25:
       #print("t1 = %.4f"%(xtmp[i]))
       t1.append(xtmp[i])
    
    if len(t1)>0:
       break

t1=numpy.array(t1)

for t1_tmp in t1:
    index=xtmp>=t1_tmp
    print("t2 = %6.2f, K2 = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))





##### Ek-all
DEK = EK[:,0]+EK[:,1]+EK[:,2]+EK[:,3]
DEK = DEK-DEK[0]

xtmp=time
ytmp=-1.0*DEK/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="--",linewidth=2.0,color="r",label=r"$-\Delta \mathcal{E}_{k}$")


for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, K-all = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))



#haxe.legend(loc="upper right",frameon=False,ncol=2,fontsize=12)

haxe.grid(linestyle="--",linewidth=0.5,color="w")



xmin,xmax =-5,165
xticks=numpy.arange(0,160+20,20)
xticklabel=[r"$%d$"%(tmp) for tmp in xticks]

ymin,ymax =-0.03,0.3
yticks=numpy.arange(-0.0,0.6+0.1,0.1)
print(yticks)
yticklabel=[r"$%.2f$"%(tmp) for tmp in yticks]

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)
haxe.set_xticklabels(xticklabel)


haxe.set_ylim(ymin,ymax)
haxe.set_yticks(yticks)
#haxe.set_yticklabels(yticklabel)



haxe.tick_params(direction='in', length=5, width=1, colors='k')

haxe.tick_params(axis="x",labelsize=18)
haxe.tick_params(axis="y",labelsize=18)

# label
xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.04*(xlim[1]-xlim[0]),ylim[0]+0.85*(ylim[1]-ylim[0]),r"$(b2)$",fontsize=20)



#haxe.set_xlabel(r"$t\cdot \omega_{pe}$",fontsize=24)

haxe.set_ylabel(r"$\Delta \mathcal{E}/\mathcal{E}_{k0}$",fontsize=24)







############################################################
############################################################
rid,cid=2,1
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=size,windows=windows,margin=margin)
axes_pos = [0.56,0.08,0.4,0.26333333]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = %d,%d"%(rid,cid))
print("axes_pos =",axes_pos)


dir_tmp = list_dirname[rid]
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

time = numpy.arange(EE.shape[0])*wpe*dt


index=timestep_state>0
time=time[index]
EE=EE[index,:]
EB=EB[index,:]
EK=EK[index,:]

print("EE =",EE.shape)
print("EB =",EB.shape)
print("EK =",EK.shape)



norm_E = EK[0,2]

print("norm_E =",norm_E)


print("*"*40+" Run3")


####################

NUM_INTP = 800


##### E+B
DEB = EE[:,0]+EE[:,1]+EE[:,2]+EB[:,0]+EB[:,1]+EB[:,2]
DEB = DEB-DEB[0]

xtmp=time
ytmp=DEB/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

#haxe.plot(xtmp,ytmp,linestyle="--",linewidth=2.0,color="darkblue")

### vertical t0
t0=[]
for i in range(1,len(ytmp)-1):
    if (ytmp[i]>ytmp[i-1] and ytmp[i]>ytmp[i+1]) and xtmp[i]>25:
       print("t0 = %.4f"%(xtmp[i]))
       t0.append(xtmp[i])
    
    if len(t0)>1:
       break

t0=numpy.array(t0)

for i in range(2):
    haxe.axvline(x=t0[i],linestyle="--",linewidth=1.0,color="grey")
    haxe.text(t0[i]-10,0.2,r"$t_%d=%.2f$"%(i,t0[i]),color="grey",rotation=90)

for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, E+B = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))


### scatter
xtmp=numpy.linspace(0,numpy.max(xtmp),50)
ytmp=fun_intp(xtmp)

index=numpy.logical_or(xtmp<20,xtmp>35)
#haxe.scatter(xtmp[index],ytmp[index],marker="o",s=25,color="darkblue",label=r"$\Delta \mathcal{E}_{E}+\Delta \mathcal{E}_{B}$")
haxe.scatter(xtmp[index],ytmp[index],marker="o",s=25,color="darkblue")

xtmp=numpy.array([28,32,35,37,40.5,42.5,48])
ytmp=fun_intp(xtmp)
haxe.scatter(xtmp,ytmp,marker="o",s=25,color="darkblue")



##### Ex
DEE = EE[:,0]
DEE = DEE-DEE[0]

xtmp=time
ytmp=DEE/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="m",label=r"$\Delta \mathcal{E}_{E_{x}}$")

for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, Ex  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))


##### Ey
DEE = EE[:,1]
DEE = DEE-DEE[0]

xtmp=time
ytmp=DEE/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="b",label=r"$\Delta \mathcal{E}_{E_{y}}$")


for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, Ey  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))




##### Bz
DEB = EB[:,2]
DEB = DEB-DEB[0]

xtmp=time
ytmp=DEB/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="grey",label=r"$\Delta \mathcal{E}_{B_{z}}$")

for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, Bz  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))



haxe.scatter(numpy.arange(10),numpy.zeros(10)-10,marker="o",s=25,color="darkblue",label=r"$\Delta \mathcal{E}_{E}+\Delta \mathcal{E}_{B}$")

#haxe.scatter(numpy.arange(10),numpy.zeros(10)-10,marker="o",s=25,color="None",label=r"$\ $")



####################


##### Ek0
DEK = EK[:,0]
DEK = DEK-DEK[0]

xtmp=time
ytmp=DEK/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="g",label=r"$\Delta \mathcal{E}_{k1}^{e}$")

for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, K0  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))




##### Ek2
DEK = EK[:,2]
DEK = DEK-DEK[0]

xtmp=time
ytmp=-1.0*DEK/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="-",linewidth=2.0,color="r",label=r"$-\Delta \mathcal{E}_{k2}^{e}$")


for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, K2  = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))



##### Ek-all
DEK = EK[:,0]+EK[:,1]+EK[:,2]+EK[:,3]
DEK = DEK-DEK[0]

xtmp=time
ytmp=-1.0*DEK/norm_E
fun_intp=scipy.interpolate.interp1d(xtmp,ytmp,kind = 'quadratic')
xtmp=numpy.linspace(0,numpy.max(xtmp),NUM_INTP)
ytmp=fun_intp(xtmp)

haxe.plot(xtmp,ytmp,linestyle="--",linewidth=2.0,color="r",label=r"$-\Delta \mathcal{E}_{k}$")


for t0_tmp in t0:
    index=xtmp>=t0_tmp
    print("t = %6.2f, K-all = %8.4f"%(xtmp[index][0],ytmp[index][0]*100))



#haxe.legend(loc="upper right",frameon=False,ncol=2,fontsize=12)

haxe.grid(linestyle="--",linewidth=0.5,color="w")



xmin,xmax =-5,165
xticks=numpy.arange(0,160+20,20)
xticklabel=[r"$%d$"%(tmp) for tmp in xticks]

ymin,ymax =-0.02,0.42
yticks=numpy.arange(-0.0,0.4+0.1,0.1)
print(yticks)
yticklabel=[r"$%.2f$"%(tmp) for tmp in yticks]

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)
haxe.set_xticklabels(xticklabel)


haxe.set_ylim(ymin,ymax)
haxe.set_yticks(yticks)
#haxe.set_yticklabels(yticklabel)



haxe.tick_params(direction='in', length=5, width=1, colors='k')

haxe.tick_params(axis="x",labelsize=18)
haxe.tick_params(axis="y",labelsize=18)

# label
xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.04*(xlim[1]-xlim[0]),ylim[0]+0.85*(ylim[1]-ylim[0]),r"$(c2)$",fontsize=20)



haxe.set_xlabel(r"$t\cdot \omega_{pe}$",fontsize=24)

haxe.set_ylabel(r"$\Delta \mathcal{E}/\mathcal{E}_{k0}$",fontsize=24)




############################################################

#plt.show()
    
fig_path = path_fig
fig_name = os.path.splitext(os.path.basename(__file__))[0]


yaoxpy.fig_save(plt,figpath=fig_path,figname=fig_name,extension=fig_fmts)