from YaoxPy_Import_CWD import *

from YaoxPy_Wave_Equations_Two_Electrons import *

from matplotlib.patches import Rectangle,ConnectionPatch

import scipy

##################################################################

alpha0 = 400/420
alpha1 = 20/420

vd0    = -0.015*CGS["c"]
vd1    =  0.3*CGS["c"]


vthe0  = 0.03*CGS["c"]
vthe1  = 0.03*CGS["c"]

mu     = 1836

wpe    = 5e9

wpi    = wpe/numpy.sqrt(mu)


######################################################

norm_w = wpe

norm_k = wpe/CGS["c"]

dk     = 0.01*norm_k


kmin,kmax = -35,35

ymin,ymax = -3,3

symmetry_mode = 2


list_color06=["r","coral","m","m","coral","r"]

list_color12=["r","g","m","coral","m","grey","grey","m","coral","m","g","r"]



lwid = 1.3


######################################################


xmin1,xmax1=-22,22
xticks1=numpy.arange(-20,20+5,5)

ymin1,ymax1=-0.5,4.5
yticks1=numpy.arange(0,4+1,1)


xmin2,xmax2=-20,20
xticks2=numpy.arange(-20,20+5,5)

ymin2,ymax2=-0.01,0.02
yticks2=numpy.arange(-0.01,0.02+0.01,0.01)




######################################################
######################################################

hfig = plt.figure(figsize=(15,9))


margin=[0.05,0.04,0.14,0.1,0.06,0.04]
barbox=[0.05,0.01,0.8]
windows=[8,2]
size=[1,1]




######################################################
######################################################


rid,cid=0,0

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[4,1],windows=windows,margin=margin)
axes_pos = [0.05,0.5,0.425,0.36 ]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")




########################################
k = numpy.arange(-3000,3000+1,1)*0.01*norm_k

w = wave_equation_two_electrons_twelveth_solve(
    wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, k,
    a=1.0, b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8
)

wr=w.real


for i in [0]:
    haxe.plot(k/norm_k,wr[i,:]/norm_w,linestyle="-",linewidth=lwid,color=list_color12[i])

for i in [1,3,4]:
    index=k>=0
    haxe.plot(k[index]/norm_k,wr[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color12[i])

for i in [8,7]:
    index=k<=0
    haxe.plot(k[index]/norm_k,wr[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color12[i])


for i in [10]:
    index=k<=0
    haxe.plot(k[index]/norm_k,wr[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color12[i])


######################################## Bohm-Gross

k = numpy.arange(-3000,3000+1,1)*0.01*norm_k

wpe0=wpe*alpha0
w = numpy.sqrt(wpe0*wpe0+3.0*numpy.power(k*vthe0,2.0))

haxe.plot(k/norm_k,w/norm_w,linestyle="-.",linewidth=lwid,color="r")



######################################## EA

k = numpy.arange(0,3000+1,1)*0.01*norm_k

#w = wave_equation_EA(wpe,mu,alpha0,alpha1,vd1,vthe0,vthe1,k)
#haxe.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="k")

w = wave_equation_EA_gary(wpe,mu,alpha0,alpha1,vd1,vthe0,vthe1,k)
haxe.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="k")


######################################## IA

k = numpy.arange(-3000,3000+1,1)*0.01*norm_k

#w = wave_equation_two_electrons_IA_solve(wpe,mu,alpha0,alpha1,vthe0,vthe1,k)
#haxe.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="grey")

w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0,alpha1,vthe0,vthe1,k)
haxe.plot(k/norm_k,w/norm_w,linestyle="-",linewidth=lwid,color="b")


k2 = numpy.arange(500,3000+1,1)*0.01*norm_k
cs = wave_equation_two_electrons_IA_cs(wpe,mu,alpha0,alpha1,vthe0,vthe1)
w  = wpe+numpy.zeros(len(k2))
haxe.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="m")

#haxe.text(10,0.6,r"$\omega=k\cdot c_s+\omega_{pe1}$",color="m")



########################################
w=-5*numpy.ones(len(k))

haxe.plot(k/norm_k,w,linestyle="-",linewidth=lwid,label=r"$Beam\textendash Beam$",color="r")

haxe.plot(k/norm_k,w,linestyle="-.",linewidth=lwid,label=r"$Langmuir$",color="r")


haxe.plot(k/norm_k,w,linestyle="-",linewidth=lwid,label=r"$Beam$",color="g")

haxe.plot(k/norm_k,w,linestyle="--",linewidth=lwid,label=r"$ \omega=v_b\cdot k$",color="g")

haxe.plot(k/norm_k,w,linestyle="-",linewidth=lwid,label=r"$Beam\textendash like$",color="coral")

haxe.plot(k/norm_k,w,linestyle="-",linewidth=lwid,label=r"$\ $",color="w")


haxe.plot(k/norm_k,w,linestyle="-",linewidth=lwid,label=r"$hybrid$",color="m")

haxe.plot(k/norm_k,w,linestyle="--",linewidth=lwid,label=r"$\omega=\omega_{pe}$",color="m")

haxe.plot(k/norm_k,w,linestyle="--",linewidth=lwid,label=r"$EA$",color="k")

haxe.plot(k/norm_k,w,linestyle="-",linewidth=lwid,label=r"$\ $",color="w")


haxe.plot(k/norm_k,w,linestyle="-",linewidth=lwid,label=r"$IA$",color="b")

haxe.plot(k/norm_k,w,linestyle="--",linewidth=lwid,label=r"$\omega=c_s\cdot k$",color="b")


#haxe.plot(k/norm_k,w,linestyle="--",linewidth=lwid,label=r"$IA$",color="grey")

#haxe.plot(k/norm_k,w,linestyle="--",linewidth=lwid,label=r"$\omega=k\cdot c_s$",color="b")



haxe.legend(loc="upper left",bbox_to_anchor=(-0.03,1.35),ncol=6,frameon=False,prop={"size":17})

#haxe.legend(loc="upper left",bbox_to_anchor=(1.11,1.05),ncol=1,frameon=False,prop={"size":17})






########################################
haxe.set_xlim(xmin1,xmax1)
haxe.set_ylim(ymin1,ymax1)

haxe.set_xticks(xticks1)
haxe.set_yticks(yticks1)

haxe.grid(linestyle="--",linewidth=0.8,color="w")


xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.9*(ylim[1]-ylim[0]),r"$(a1)$",fontsize=20)


#haxe.set_xlabel(r"$k\cdot d_e$",fontsize=24)
haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=24)




######################################################
######################################################

rid,cid=0,1

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[2,1],windows=windows,margin=margin)
axes_pos = [0.535,0.7,0.425,0.16 ]
haxe2=hfig.add_axes(axes_pos,facecolor="whitesmoke")



########################################
k = numpy.arange(-3000,3000+1,1)*0.01*norm_k


w = wave_equation_two_electrons_twelveth_solve(
    wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, k,
    a=1.0, b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8
)

wr=numpy.real(w)



#for i in [0]:
#    haxe2.plot(k/norm_k,wr[i,:]/norm_w,linestyle="-",linewidth=lwid,color=list_color12[i])

for i in [1,3,4]:
    index=k>=0
    haxe2.plot(k[index]/norm_k,wr[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color12[i])

for i in [8,7]:
    index=k<=0
    haxe2.plot(k[index]/norm_k,wr[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color12[i])


for i in [10]:
    index=k<=0
    haxe2.plot(k[index]/norm_k,wr[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color12[i])





######################################## IA

k = numpy.arange(-3000,3000+1,1)*0.01*norm_k

#w = wave_equation_two_electrons_IA_solve(wpe,mu,alpha0,alpha1,vthe0,vthe1,k)
#haxe2.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="grey")

w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0,alpha1,vthe0,vthe1,k)
haxe2.plot(k/norm_k,w/norm_w,linestyle="-",linewidth=lwid,color="b")


#haxe2.text(6,-0.003,r"$\omega=c_sk/\sqrt{1+\lambda_D^2k^2}$",color="b",rotation=9,fontsize=16)





cs = wave_equation_two_electrons_IA_cs(wpe,mu,alpha0,alpha1,vthe0,vthe1)
w  = numpy.abs(k*cs)
haxe2.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="b")

#haxe2.text(10,0.012,r"$\omega=k\cdot c_s$",rotation=10,color="b")


print("ion-acoustic speed v = %.8fc"%(cs/CGS["c"]))
print("ion-acoustic speed w = %.8fwpe"%(cs*norm_k/norm_w*20))








haxe2.axhline(y=0.0,linestyle="--",linewidth=0.5,color="k")


haxe2.grid(linestyle="--",linewidth=0.8,color="w")

haxe2.set_xlim(xmin2,xmax2)
haxe2.set_ylim(ymin2,ymax2)

haxe2.set_xticks(xticks2)
haxe2.set_yticks(yticks2)



xlim=haxe2.get_xlim()
ylim=haxe2.get_ylim()
haxe2.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(a2)$",fontsize=20)


haxe2.set_xlabel(r"$k\cdot d_e$",fontsize=24)
#haxe2.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=24)





##################################################

color_box="grey"

haxe.add_patch(Rectangle((xmin2,ymin2),xmax2-xmin2,ymax2-ymin2,facecolor="none",edgecolor=color_box,linestyle="--"))

haxe2.spines["bottom"].set_color(color_box)
haxe2.spines["top"].set_color(color_box)
haxe2.spines["left"].set_color(color_box)
haxe2.spines["right"].set_color(color_box)


con1 = ConnectionPatch(xyA=(xmax2,ymin2),
                       xyB=(xmin2,ymin2),
                       coordsA="data",
                       coordsB="data",
                       axesA=haxe,
                       axesB=haxe2,
                       linestyle="--",
                       color=color_box)

con2 = ConnectionPatch(xyA=(xmax2,ymax2),
                       xyB=(xmin2,ymax2),
                       coordsA="data",
                       coordsB="data",
                       axesA=haxe,
                       axesB=haxe2,
                       linestyle="--",
                       color=color_box)


haxe2.add_artist(con1)
haxe2.add_artist(con2)







######################################################
######################################################


rid,cid=4,0

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[4,1],windows=windows,margin=margin)
axes_pos = [0.05,0.1,0.425,0.36 ]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")




########################################
k = numpy.arange(-3000,3000+1,1)*0.01*norm_k


w = two_stream_wave_equation_electron_sixth_solve(
    wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, k,
    a=1.0, b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8
)

wr=numpy.real(w)

for i in [0]:
    haxe.plot(k/norm_k,wr[i,:]/norm_w,linestyle="-",linewidth=lwid,color=list_color06[i])

for i in [1,2]:
    index=k>=0
    haxe.plot(k[index]/norm_k,wr[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color06[i])


for i in [4,3]:
    index=k<=0
    haxe.plot(k[index]/norm_k,wr[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color06[i])





######################################## Bohm-Gross
k = numpy.arange(-3000,3000+1,1)*0.01*norm_k

wpe0=wpe*alpha0
w = numpy.sqrt(wpe0*wpe0+3.0*numpy.power(k*vthe0,2.0))

haxe.plot(k/norm_k,w/norm_w,linestyle="-.",linewidth=lwid,color="r")




######################################## EA
k = numpy.arange(0,3000+1,1)*0.01*norm_k

#w = wave_equation_EA(wpe,mu,alpha0,alpha1,vd1,vthe0,vthe1,k)
#haxe.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="k")

w = wave_equation_EA_gary(wpe,mu,alpha0,alpha1,vd1,vthe0,vthe1,k)
haxe.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="k")




######################################## beam mode
k2 = numpy.arange(-30+0.01,0.0,0.01)*norm_k

w  = vd0*k2
haxe.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="g")


k2 = numpy.arange(0,30+0.01,0.01)*norm_k
w  = vd1*k2
haxe.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="g")






######################################## IA

k = numpy.arange(-3000,3000+1,1)*0.01*norm_k

#w = wave_equation_two_electrons_IA_solve(wpe,mu,alpha0,alpha1,vthe0,vthe1,k)
#haxe.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="grey")

w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0,alpha1,vthe0,vthe1,k)
haxe.plot(k/norm_k,w/norm_w,linestyle="-",linewidth=lwid,color="b")



########################################
haxe.set_xlim(xmin1,xmax1)
haxe.set_ylim(ymin1,ymax1)

haxe.set_xticks(xticks1)
haxe.set_yticks(yticks1)

haxe.grid(linestyle="--",linewidth=0.8,color="w")





xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.9*(ylim[1]-ylim[0]),r"$(b1)$",fontsize=20)


haxe.set_xlabel(r"$k\cdot d_e$",fontsize=24)
haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=24)



######################################################
######################################################

rid,cid=4,1

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[2,1],windows=windows,margin=margin)
axes_pos = [0.535,0.3,0.425,0.16 ]
haxe2=hfig.add_axes(axes_pos,facecolor="whitesmoke")



########################################
k = numpy.arange(-3000,3000+1,1)*0.01*norm_k


w = two_stream_wave_equation_electron_sixth_solve(
    wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, k,
    a=1.0, b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8
)

wr=numpy.real(w)



#for i in [0]:
#    haxe.plot(k/norm_k,wr[i,:]/norm_w,linestyle="-",linewidth=lwid,color=list_color06[i])


for i in [1,2]:
    index=k>=0
    haxe2.plot(k[index]/norm_k,wr[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color06[i])


for i in [4,3]:
    index=k<=0
    haxe2.plot(k[index]/norm_k,wr[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color06[i])




######################################## Bohm-Gross

#k = numpy.arange(-3000,3000+1,1)*0.01*norm_k

#w = numpy.sqrt(wpe*wpe+3.0*numpy.power(k*vthe0,2.0))

#haxe2.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="r")



######################################## beam mode

k1 = numpy.arange(-3000,0+1,1)*0.01*norm_k
w  = vd0*k1
haxe2.plot(k1/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="g")


k = numpy.arange(0,3000+1,1)*0.01*norm_k
w  = vd1*k2
haxe2.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="g")




######################################## IA

k = numpy.arange(-3000,3000+1,1)*0.01*norm_k

#w = wave_equation_two_electrons_IA_solve(wpe,mu,alpha0,alpha1,vthe0,vthe1,k)
#haxe2.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="grey")

w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0,alpha1,vthe0,vthe1,k)
haxe2.plot(k/norm_k,w/norm_w,linestyle="-",linewidth=lwid,color="b")



cs = wave_equation_two_electrons_IA_cs(wpe,mu,alpha0,alpha1,vthe0,vthe1)
w  = numpy.abs(k*cs)
haxe2.plot(k/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="b")

#haxe2.text(10,0.012,r"$\omega=k\cdot c_s$",rotation=10,color="b")


haxe2.axhline(y=0.0,linestyle="--",linewidth=0.5,color="k")

haxe2.grid(linestyle="--",linewidth=0.8,color="w")

haxe2.set_xlim(xmin2,xmax2)
haxe2.set_ylim(ymin2,ymax2)

haxe2.set_xticks(xticks2)
haxe2.set_yticks(yticks2)



xlim=haxe2.get_xlim()
ylim=haxe2.get_ylim()
haxe2.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(b2)$",fontsize=20)


haxe2.set_xlabel(r"$k\cdot d_e$",fontsize=24)
#haxe2.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=24)





##################################################

color_box="grey"

haxe.add_patch(Rectangle((xmin2,ymin2),xmax2-xmin2,ymax2-ymin2,facecolor="none",edgecolor=color_box,linestyle="--"))

haxe2.spines["bottom"].set_color(color_box)
haxe2.spines["top"].set_color(color_box)
haxe2.spines["left"].set_color(color_box)
haxe2.spines["right"].set_color(color_box)


con1 = ConnectionPatch(xyA=(xmax2,ymin2),
                       xyB=(xmin2,ymin2),
                       coordsA="data",
                       coordsB="data",
                       axesA=haxe,
                       axesB=haxe2,
                       linestyle="--",
                       color=color_box)

con2 = ConnectionPatch(xyA=(xmax2,ymax2),
                       xyB=(xmin2,ymax2),
                       coordsA="data",
                       coordsB="data",
                       axesA=haxe,
                       axesB=haxe2,
                       linestyle="--",
                       color=color_box)


haxe2.add_artist(con1)
haxe2.add_artist(con2)








############################################################
#plt.show()

fig_path = path_fig
fig_name = os.path.splitext(os.path.basename(__file__))[0]


yaoxpy.fig_save(plt,figpath=fig_path,figname=fig_name,extension=fig_fmts)