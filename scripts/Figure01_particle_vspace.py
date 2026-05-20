from YaoxPy_Import_CWD import *


############################################################
############################################################

list_dirname = ["yaoxpic_v25_counter_1","yaoxpic_v25_counter_2","yaoxpic_v25_counter_3"]


color_beam1 = "#003171"
color_beam2 = "coral"
color_vdf   = "r"

color_axis  = "gray"



############################################################
############################################################

hfig = plt.figure(figsize=(15,9))

margin=[0.0,0.05,0.03,0.1,0.04,0.05]
barbox=[0.05,0.01,0.8]
windows=[3,5]
size=[1,1]


############################################################
rid,cid=0,0
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[3,3],windows=windows,margin=margin)
#axes_pos[0]-=0.04

axes_pos = [-0.04,0.1,0.554,0.87]

haxe=hfig.add_axes(axes_pos,projection="3d")

print("rid,cid  = %d,%d"%(rid,cid))
print("axes_pos =",axes_pos)



Timestep=0

dir_tmp = list_dirname[0]
path_pic_tmp  = os.path.join(path_data_pic,dir_tmp)
path_data_tmp = os.path.join(path_pic_tmp,"data")

filename = "timestep_%d_particle.h5"%(Timestep)
    
H5FILE_R=h5py.File(os.path.join(path_data_tmp,filename),"r")

sid=0
vx0 = H5FILE_R["/species_%d/vx"%(sid)][()]
vy0 = H5FILE_R["/species_%d/vy"%(sid)][()]
vz0 = H5FILE_R["/species_%d/vz"%(sid)][()]
    
sid=2
vx2 = H5FILE_R["/species_%d/vx"%(sid)][()]
vy2 = H5FILE_R["/species_%d/vy"%(sid)][()]
vz2 = H5FILE_R["/species_%d/vz"%(sid)][()]

H5FILE_R.close()

'''
dis=numpy.sqrt((vx0-numpy.mean(vx0))**2+vy0**2+vz0**2)
index=numpy.arange(0,len(vx0),200)
him0=haxe.scatter(vx0[index],vy0[index],vz0[index],marker="o",s=5,c=dis[index],cmap=mpl.cm.jet,rasterized=True)
'''


index=numpy.arange(0,len(vx0),200)
him2=haxe.scatter(vx0[index],vy0[index],vz0[index],marker="o",s=5,rasterized=True,color=color_beam1)


index=numpy.arange(0,len(vx2),200)
him2=haxe.scatter(vx2[index],vy2[index],vz2[index],marker="o",s=5,rasterized=True,color=color_beam2)




########## axis
xmin,xmax=0.35,0.45
yaoxpy.arrow3d(haxe,[-0.50,-0.08,0.05],[-0.38,-0.08,0.05],col=color_axis,ra=0.0)
yaoxpy.arrow3d(haxe,[-0.15,-0.08,0.05],[0.195,-0.08,0.05],col=color_axis,ra=0.0)
yaoxpy.arrow3d(haxe,[0.45,-0.08,0.05],[0.55,-0.08,0.05],col=color_axis,ra=0.35)
haxe.text(0.55,-0.05,0.0,r"$v_{\parallel}$",color=color_axis,fontsize=25)


ymin,ymax=-0.02,-0.25
yaoxpy.arrow3d(haxe,[0.011,ymin,0.0],[0.011,ymax,0.0],col=color_axis,ra=0.15)
haxe.text(0.0,-0.35,0.0,r"$v_{\perp2}$",color=color_axis,fontsize=25)


zmin,zmax=-0.02,0.25
yaoxpy.arrow3d(haxe,[0.0,0.0,zmin],[0.0,0.0,zmax],col=color_axis,ra=0.13)
haxe.text(-0.02,0.0,0.28,r"$v_{\perp1}$",color=color_axis,fontsize=25)



haxe.text(0.0,0.0,-0.09,r"$O$",color=color_axis,fontsize=20)




#haxe.text(-0.22,0.0,0.1,r'$background$',color=YXColorBlue,fontsize=25)


haxe.text(-0.48,-0.05,-0.27,r'$Population 1$',color=color_beam1,fontsize=25)

haxe.text(0.18,-0.05,-0.27,r'$Population 2$',color=color_beam2,fontsize=25)



xmin,xmax=-0.45,0.45
haxe.set_xlim(xmin,xmax)

ymin,ymax=-0.40,0.40
haxe.set_ylim(ymin,ymax)

zmin,zmax=-0.40,0.40
haxe.set_zlim(zmin,zmax)


haxe.text(-0.25,-0.1,0.35,r"$(a)$",color="k",fontsize=20)


haxe.set_axis_off()

#haxe.text(-0.07,0.0,0.45,r"$Timestep = %4d : t=%3.1f\omega_{pe}^{-1}$"%(Timestep,Timestep*dt*wpe),color="k")



#haxe.set_title(r"$Cloud$")

#haxe.set_xlabel(r"$v_x$",fontsize=20)
#haxe.set_ylabel(r"$v_y$",fontsize=20)
#haxe.set_zlabel(r"$v_z$",fontsize=20)


haxe.view_init(elev=30,azim=-65,roll=0)








############################################################
rid,cid=0,3
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,2],windows=windows,margin=margin)
#axes_pos[0]-=0.05
#axes_pos[2]+=0.05

axes_pos = [0.544,0.71333333,0.406,0.25666667]

haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = %d,%d"%(rid,cid))
print("axes_pos =",axes_pos)

Timestep=0

dir_tmp = list_dirname[0]
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




haxe.plot(bin_vpara,vpara0+vpara2,linestyle="-",linewidth=1.5,color=color_vdf,label=r"$f_{e}(v_{\parallel})=f_{e1}(v_{\parallel})+f_{e2}(v_{\parallel})$")


#haxe.plot(bin_vpara,vpara0,linestyle="--",linewidth=3.0,color=YXColorBlue)

#haxe.plot(bin_vpara,vpara2,linestyle="--",linewidth=3.0,color="coral")



vb1=-0.3
vb2=0.3

index=bin_vpara<=(vb1+vb2)*0.5
x=bin_vpara[index]
y=vpara0[index]
index=numpy.arange(0,len(x),10)
haxe.scatter(x[index],y[index],marker="o",s=55,color=color_beam1,label=r"$f_{e1}(v_{\parallel})$")

index=bin_vpara>=(vb1+vb2)*0.5
x=bin_vpara[index]
y=vpara2[index]
index=numpy.arange(0,len(x),10)
haxe.scatter(x[index],y[index],marker="o",s=55,color=color_beam2,label=r"$f_{e2}(v_{\parallel})$")


#haxe.annotate(r"$\frac{\partial f_e}{\partial v_{\parallel}}>0$",xy=(0.15,0.8),xytext=(0.06,1.1),arrowprops=dict(facecolor="r",edgecolor="r",width=0.8,headwidth=4.0,headlength=4.0),color="r",fontsize=20)

haxe.quiver(-0.15,4,-0.8,-0.3,scale_units="width",width=0.005,linewidth=0.8,color=color_vdf)

haxe.quiver(0.15,4,0.8,-0.3,scale_units="width",width=0.005,linewidth=0.8,color=color_vdf)

haxe.text(-0.08,5,r"$v_{\parallel}\cdot\frac{\partial f_e}{\partial v_{\parallel}}>0$",color=color_vdf,fontsize=20)




haxe.grid(True,linestyle="--",linewidth=1,color="w")

#haxe.axvline(x=0.0,linestyle="--",linewidth=1.0,color=YXColorBlue)
haxe.axvline(x=vb1,linestyle="--",linewidth=1.0,color=color_beam1)
haxe.axvline(x=vb2,linestyle="--",linewidth=1.0,color=color_beam2)

haxe.legend(loc="upper right",frameon=False,fontsize=12)



xmin,xmax=-0.42,0.42
haxe.set_xlim(xmin,xmax)


ymin,ymax=-1,16
haxe.set_ylim(ymin,ymax)
haxe.set_yticks(numpy.arange(0,15+3,3))


xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(b1)Run1$",fontsize=20)






#haxe.set_xlabel(r"$v_{\parallel}/c$",fontsize=24)
haxe.set_ylabel(r"$f(v_{\parallel})$",fontsize=24)








############################################################
rid,cid=1,3
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,2],windows=windows,margin=margin)
#axes_pos[0]-=0.05
#axes_pos[2]+=0.05

axes_pos=[0.544,0.40666667,0.406,0.25666667]


haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = %d,%d"%(rid,cid))
print("axes_pos =",axes_pos)

Timestep=0

dir_tmp = list_dirname[1]
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

haxe.plot(bin_vpara,vpara0+vpara2,linestyle="-",linewidth=1.5,color=color_vdf,label=r"$f_{e}(v_{\parallel})=f_{e1}(v_{\parallel})+f_{e2}(v_{\parallel})$")


#haxe.plot(bin_vpara,vpara0,linestyle="--",linewidth=3.0,color=YXColorBlue)

#haxe.plot(bin_vpara,vpara2,linestyle="--",linewidth=3.0,color="coral")


vb1=-0.3/4.0
vb2=0.3

index=bin_vpara<=(vb1+vb2)*0.5
x=bin_vpara[index]
y=vpara0[index]
index=numpy.arange(0,len(x),10)
haxe.scatter(x[index],y[index],marker="o",s=55,color=color_beam1,label=r"$f_{e1}(v_{\parallel})$")

index=bin_vpara>=(vb1+vb2)*0.5
x=bin_vpara[index]
y=vpara2[index]
index=numpy.arange(0,len(x),10)
haxe.scatter(x[index],y[index],marker="o",s=55,color=color_beam2,label=r"$f_{e2}(v_{\parallel})$")


#haxe.annotate(r"$\frac{\partial f_e}{\partial v_{\parallel}}>0$",xy=(0.15,0.8),xytext=(0.06,1.1),arrowprops=dict(facecolor="r",edgecolor="r",width=0.8,headwidth=4.0,headlength=4.0),color="r",fontsize=20)

haxe.quiver(0.05,3.5,-0.8,-0.6,scale_units="width",width=0.005,linewidth=0.8,color=color_vdf)

haxe.quiver(0.2,3.5,0.8,-0.6,scale_units="width",width=0.005,linewidth=0.8,color=color_vdf)

haxe.text(0.06,5,r"$v_{\parallel}\cdot\frac{\partial f_e}{\partial v_{\parallel}}>0$",color=color_vdf,fontsize=20)




haxe.grid(True,linestyle="--",linewidth=1,color="w")

#haxe.axvline(x=0.0,linestyle="--",linewidth=1.0,color=YXColorBlue)
haxe.axvline(x=vb1,linestyle="--",linewidth=1.0,color=color_beam1)
haxe.axvline(x=vb2,linestyle="--",linewidth=1.0,color=color_beam2)

haxe.legend(loc="upper right",frameon=False,fontsize=12)



xmin,xmax=-0.42,0.42
haxe.set_xlim(xmin,xmax)


ymin,ymax=-1,16
haxe.set_ylim(ymin,ymax)
haxe.set_yticks(numpy.arange(0,15+3,3))


xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(b2)Run2$",fontsize=20)






#haxe.set_xlabel(r"$v_{\parallel}/c$",fontsize=24)
haxe.set_ylabel(r"$f(v_{\parallel})$",fontsize=24)









############################################################
rid,cid=2,3
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,2],windows=windows,margin=margin)
#axes_pos[0]-=0.05
#axes_pos[2]+=0.05

axes_pos = [0.544,0.1,0.406,0.25666667]

haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

print("rid,cid  = %d,%d"%(rid,cid))
print("axes_pos =",axes_pos)

Timestep=0

dir_tmp = list_dirname[2]
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

haxe.plot(bin_vpara,vpara0+vpara2,linestyle="-",linewidth=1.5,color=color_vdf,label=r"$f_{e}(v_{\parallel})=f_{e1}(v_{\parallel})+f_{e2}(v_{\parallel})$")


#haxe.plot(bin_vpara,vpara0,linestyle="--",linewidth=3.0,color=YXColorBlue)

#haxe.plot(bin_vpara,vpara2,linestyle="--",linewidth=3.0,color="coral")



vb1=-0.3/20.0
vb2=0.3

index=bin_vpara<=(vb1+vb2)*0.5
x=bin_vpara[index]
y=vpara0[index]
index=numpy.arange(0,len(x),10)
haxe.scatter(x[index],y[index],marker="o",s=55,color=color_beam1,label=r"$f_{e1}(v_{\parallel})$")

index=bin_vpara>=(vb1+vb2)*0.5
x=bin_vpara[index]
y=vpara2[index]
index=numpy.arange(0,len(x),10)
haxe.scatter(x[index],y[index],marker="o",s=55,color=color_beam2,label=r"$f_{e2}(v_{\parallel})$")


#haxe.annotate(r"$\frac{\partial f_e}{\partial v_{\parallel}}>0$",xy=(0.15,0.8),xytext=(0.06,1.1),arrowprops=dict(facecolor="r",edgecolor="r",width=0.8,headwidth=4.0,headlength=4.0),color="r",fontsize=20)

#haxe.quiver(-0.08,3,-0.8,-0.5,scale_units="width",width=0.005,linewidth=0.8,color=color_vdf)

haxe.quiver(0.2,3.5,0.8,-0.6,scale_units="width",width=0.005,linewidth=0.8,color=color_vdf)

haxe.text(0.08,5,r"$v_{\parallel}\cdot\frac{\partial f_e}{\partial v_{\parallel}}>0$",color=color_vdf,fontsize=20)





haxe.grid(True,linestyle="--",linewidth=1,color="w")


#haxe.axvline(x=0.0,linestyle="--",linewidth=1.0,color=YXColorBlue)
haxe.axvline(x=vb1,linestyle="--",linewidth=1.0,color=color_beam1)
haxe.axvline(x=vb2,linestyle="--",linewidth=1.0,color=color_beam2)

haxe.legend(loc="upper right",frameon=False,fontsize=12)


xmin,xmax=-0.42,0.42
haxe.set_xlim(xmin,xmax)


ymin,ymax=-1,16
haxe.set_ylim(ymin,ymax)
haxe.set_yticks(numpy.arange(0,15+3,3))


xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(b3)Run3$",fontsize=20)






haxe.set_xlabel(r"$v_{\parallel}/c$",fontsize=24)
haxe.set_ylabel(r"$f(v_{\parallel})$",fontsize=24)









############################################################

#plt.show()
    
fig_path = path_fig
fig_name = os.path.splitext(os.path.basename(__file__))[0]


yaoxpy.fig_save(plt,figpath=fig_path,figname=fig_name,extension=fig_fmts)