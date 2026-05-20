from YaoxPy_Import_CWD import *

from YaoxPy_Wave_Equations_Two_Electrons import *


##################################################################
##################################################################


vthe0  = 0.03*CGS["c"]
vthe1  = 0.03*CGS["c"]

mu     = 1836

wpe    = 5e9

wpi    = wpe/numpy.sqrt(mu)


######################################################

norm_w = wpe

norm_k = wpe/CGS["c"]

dk     = 0.01*norm_k

list_color03 = ["r","g","m"]

######################################################

k = numpy.arange(-3000,3000+1,1)*0.01*norm_k


##### Run1
alpha0_run1 = 0.5
alpha1_run1 = 0.5
vd0_run1    = -0.3*CGS["c"]
vd1_run1    =  0.3*CGS["c"]

w_run1 = wave_equation_two_electrons_twelveth_solve(
    wpe, mu, alpha0_run1, alpha1_run1, vd0_run1, vd1_run1, vthe0, vthe1, k,
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

wr_run1=w_run1.real



##### Run2
alpha0_run2 = 4.0/5.0
alpha1_run2 = 1.0/5.0
vd0_run2    = -0.075*CGS["c"]
vd1_run2    =  0.3*CGS["c"]

w_run2 = wave_equation_two_electrons_twelveth_solve(
    wpe, mu, alpha0_run2, alpha1_run2, vd0_run2, vd1_run2, vthe0, vthe1, k,
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

wr_run2=w_run2.real


##### Run3
alpha0_run3 = 20.0/21.0
alpha1_run3 = 1.0/21.0
vd0_run3    = -0.015*CGS["c"]
vd1_run3    =  0.3*CGS["c"]

w_run3 = wave_equation_two_electrons_twelveth_solve(
    wpe, mu, alpha0_run3, alpha1_run3, vd0_run3, vd1_run3, vthe0, vthe1, k,
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

wr_run3=w_run3.real



######################################################
######################################################

xmin,xmax=-32,32
xticks=numpy.arange(-30,30+5,5)



lwid=1.3


######################################################
######################################################

hfig = plt.figure(figsize=(13,9))


margin=[0.08,0.03,0.06,0.1,0.1,0.04]
barbox=[0.05,0.01,0.8]
windows=[5,2]
size=[1,1]



######################################################
######################################################
rid,cid=0,0

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
axes_pos = [0.08,0.804,0.395,0.136]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")



for i in [0]:
    
    index=k>=-100*norm_k

    haxe.plot(k[index]/norm_k,wr_run2[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    haxe.plot(k[index]/norm_k,wr_run3[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

    haxe.plot(k[index]/norm_k,wr_run1[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])



##### legend

y=numpy.zeros(len(k))-10

haxe.plot(k/norm_k,y,linestyle="-",linewidth=lwid,color=list_color03[0],label=r"$Run1$")

haxe.plot(k/norm_k,y,linestyle="-",linewidth=lwid,color=list_color03[1],label=r"$Run2$")

haxe.plot(k/norm_k,y,linestyle="-",linewidth=lwid,color=list_color03[2],label=r"$Run3$")


haxe.legend(loc="lower right",fontsize=10)




haxe.grid(linestyle="--",linewidth=0.8,color="w")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.5,10.5)
haxe.set_yticks(numpy.arange(0,10+2,2))



xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)Beam\textendash Beam\ mode$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)


haxe.set_title(r"$Dispersion\ Relation$")



if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\gamma/\omega_{pe}$",fontsize=20)



######################################################
######################################################
rid,cid=0,1

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
axes_pos = [0.575,0.804,0.395,0.136]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

for i in [0]:
    index=k>=-100*norm_k
    
    wtmp=w_run2[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])


    wtmp=w_run3[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

    wtmp=w_run1[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])



haxe.axhline(y=0.0,linestyle="--",linewidth=0.5,color="k")

haxe.grid(linestyle="--",linewidth=0.8,color="w")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.27,0.02)
haxe.set_yticks(numpy.arange(-0.25,0.0+0.05,0.05))




xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)



haxe.set_title(r"$Growth\ Rate$")




if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\gamma/\omega_{pe}$",fontsize=20)




######################################################
######################################################
rid,cid=1,0

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
axes_pos = [0.08,0.628,0.395,0.136]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")


for i in [1]:
    
    index=k>=0

    haxe.plot(k[index]/norm_k,wr_run2[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    haxe.plot(k[index]/norm_k,wr_run3[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

    haxe.plot(k[index]/norm_k,wr_run1[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])



for i in [10]:
    
    index=k<=0

    haxe.plot(k[index]/norm_k,wr_run2[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    haxe.plot(k[index]/norm_k,wr_run3[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

    haxe.plot(k[index]/norm_k,wr_run1[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])



haxe.grid(linestyle="--",linewidth=0.8,color="w")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)


haxe.set_ylim(-0.5,10.5)
haxe.set_yticks(numpy.arange(0,10+2,2))




xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)Beam\ mode$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)




if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\gamma/\omega_{pe}$",fontsize=20)



######################################################
######################################################
rid,cid=1,1

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
axes_pos = [0.575,0.628,0.395,0.136]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

for i in [1]:
    index=k>=0

    wtmp=w_run2[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    wtmp=w_run3[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

    wtmp=w_run1[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])



for i in [10]:
    index=k<=0
   

    wtmp=w_run2[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    wtmp=w_run3[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])


    wtmp=w_run1[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])



haxe.axhline(y=0.0,linestyle="--",linewidth=0.5,color="k")

haxe.grid(linestyle="--",linewidth=0.8,color="w")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)


haxe.set_ylim(-7e-5,1e-5)
haxe.set_yticks(numpy.arange(-6,0+1,1)*1e-5)




xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)


if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\gamma/\omega_{pe}$",fontsize=20)




######################################################
######################################################
rid,cid=2,0

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
axes_pos = [0.08,0.452,0.395,0.136]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")


for i in [3]:
    
    index=k>=0

    haxe.plot(k[index]/norm_k,wr_run2[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    haxe.plot(k[index]/norm_k,wr_run3[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

    haxe.plot(k[index]/norm_k,wr_run1[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])



for i in [8]:
    
    index=k<=0

    haxe.plot(k[index]/norm_k,wr_run2[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    haxe.plot(k[index]/norm_k,wr_run3[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

    haxe.plot(k[index]/norm_k,wr_run1[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])





haxe.grid(linestyle="--",linewidth=0.8,color="w")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.5,10.5)
haxe.set_yticks(numpy.arange(0,10+2,2))


xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)Beam\textendash like\ mode$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)




if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\gamma/\omega_{pe}$",fontsize=20)



######################################################
######################################################
rid,cid=2,1

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
axes_pos = [0.575,0.452,0.395,0.136]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

for i in [3]:
    index=k>=0

    wtmp=w_run2[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    wtmp=w_run3[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

    wtmp=w_run1[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])


for i in [8]:
    index=k<=0
    

    wtmp=w_run2[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    wtmp=w_run3[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

    wtmp=w_run1[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])


haxe.axhline(y=0.0,linestyle="--",linewidth=0.5,color="k")

haxe.grid(linestyle="--",linewidth=0.8,color="w")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.27,0.02)
haxe.set_yticks(numpy.arange(-0.25,0.0+0.05,0.05))



xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)



if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\gamma/\omega_{pe}$",fontsize=20)










######################################################
######################################################
rid,cid=3,0


#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
axes_pos = [0.08,0.276,0.395,0.136]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")


for i in [4]:
    
    index=k>=0

    haxe.plot(k[index]/norm_k,wr_run2[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    haxe.plot(k[index]/norm_k,wr_run3[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

    haxe.plot(k[index]/norm_k,wr_run1[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])


for i in [4]:
    
    index=k/norm_k>=-100

    xtmp=k[index]/norm_k
    ytmp=wr_run1[i,index]/norm_w

    spline = scipy.interpolate.UnivariateSpline(xtmp,ytmp,s=3.0)
    xtmp = numpy.linspace(0,numpy.max(xtmp),400)
    ytmp = spline(xtmp)

    haxe.plot(xtmp,ytmp,linestyle="-.",linewidth=1.5,color="r")



    xtmp=k[index]/norm_k
    ytmp=wr_run2[i,index]/norm_w

    spline = scipy.interpolate.UnivariateSpline(xtmp,ytmp,s=3.0)
    xtmp = numpy.linspace(0.0,numpy.max(xtmp),400)
    ytmp = spline(xtmp)

    haxe.plot(xtmp,ytmp,linestyle="-.",linewidth=1.5,color="g")


    xtmp=k[index]/norm_k
    ytmp=wr_run3[i,index]/norm_w

    spline = scipy.interpolate.UnivariateSpline(xtmp,ytmp,s=3.0)
    xtmp = numpy.linspace(0.0,numpy.max(xtmp),400)
    ytmp = spline(xtmp)

    haxe.plot(xtmp,ytmp,linestyle="-.",linewidth=1.5,color="m")








for i in [7]:
    
    index=k<=0

    haxe.plot(k[index]/norm_k,wr_run2[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    haxe.plot(k[index]/norm_k,wr_run3[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

    haxe.plot(k[index]/norm_k,wr_run1[i,index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])



######################################## IA

k2 = numpy.arange(500,3000+1,1)*0.01*norm_k
cs = wave_equation_two_electrons_IA_cs(wpe,mu,alpha0_run3,alpha1_run3,vthe0,vthe1)
w  = k2*cs+wpe
haxe.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="k")

print("cs = %.4fc"%(cs/CGS["c"]))

haxe.text(10,0.7,r"$\omega=\omega_{pe}$")



haxe.grid(linestyle="--",linewidth=0.8,color="w")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)


haxe.set_ylim(-0.1,1.3)
haxe.set_yticks(numpy.arange(0,1.2+0.2,0.2))



xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)Hybrid\ mode$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)




if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\gamma/\omega_{pe}$",fontsize=20)



######################################################
######################################################
rid,cid=3,1

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
axes_pos = [0.575,0.276,0.395,0.136]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")

for i in [4]:
    index=k>=0
    

    wtmp=w_run2[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    wtmp=w_run3[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

    wtmp=w_run1[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])



for i in [7]:
    index=k<=0

    wtmp=w_run2[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

    wtmp=w_run3[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])


    wtmp=w_run1[i,:]
    gamma=two_stream_growth_rate_electron(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,wtmp)
    haxe.plot(k[index]/norm_k,gamma[index]/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])



haxe.axhline(y=0.0,linestyle="--",linewidth=0.5,color="k")

haxe.grid(linestyle="--",linewidth=0.8,color="w")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.085,0.005)
haxe.set_yticks(numpy.arange(-0.08,0.0+0.02,0.02))





xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)



if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\gamma/\omega_{pe}$",fontsize=20)











######################################################
######################################################
rid,cid=4,0

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
axes_pos = [0.08,0.1,0.395,0.136]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")


lwid1=0.5

for i in [4]:
    
    index=k>=0

    haxe.plot(k[index]/norm_k,wr_run2[i,index]/norm_w,linestyle="--",linewidth=lwid1,color=list_color03[1])

    haxe.plot(k[index]/norm_k,wr_run3[i,index]/norm_w,linestyle="--",linewidth=lwid1,color=list_color03[2])

    haxe.plot(k[index]/norm_k,wr_run1[i,index]/norm_w,linestyle="--",linewidth=lwid1,color=list_color03[0])








for i in [7]:
    
    index=k<=0

    haxe.plot(k[index]/norm_k,wr_run2[i,index]/norm_w,linestyle="--",linewidth=lwid1,color=list_color03[1])

    haxe.plot(k[index]/norm_k,wr_run3[i,index]/norm_w,linestyle="--",linewidth=lwid1,color=list_color03[2])

    haxe.plot(k[index]/norm_k,wr_run1[i,index]/norm_w,linestyle="--",linewidth=lwid1,color=list_color03[0])





######################################## IA

k2 = numpy.arange(-3000,3000+1,1)*0.01*norm_k
cs = wave_equation_two_electrons_IA_cs(wpe,mu,alpha0_run3,alpha1_run3,vthe0,vthe1)
w  = numpy.abs(k2)*cs
haxe.plot(k2/norm_k,w/norm_w,linestyle="--",linewidth=lwid,color="k")

print("cs = %.4fc"%(cs/CGS["c"]))

haxe.text(10,0.012,r"$\omega=k\cdot c_s$",rotation=10)


#k = numpy.arange(-3100,3100+1,1)*0.01*norm_k

w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0_run2,alpha1_run2,vthe0,vthe1,k)
haxe.plot(k/norm_k,w/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])

w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0_run3,alpha1_run3,vthe0,vthe1,k)
haxe.plot(k/norm_k,w/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])

w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0_run1,alpha1_run1,vthe0,vthe1,k)
haxe.plot(k/norm_k,w/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])








haxe.grid(linestyle="--",linewidth=0.8,color="w")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.005,0.045)
haxe.set_yticks(numpy.arange(0,0.04+0.01,0.01))





xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)Ion\textendash Acoustic\ mode$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)




if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\gamma/\omega_{pe}$",fontsize=20)



######################################################
######################################################
rid,cid=4,1

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
axes_pos = [0.575,0.1,0.395,0.136]
haxe=hfig.add_axes(axes_pos,facecolor="whitesmoke")



w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0_run2,alpha1_run2,vthe0,vthe1,k)

gamma=wave_equation_two_electrons_IA_growth_rate_solve_2(wpe,mu,alpha0_run2,alpha1_run2,vd0_run2,vd1_run2,vthe0,vthe1,k,w)

haxe.plot(k/norm_k,gamma/norm_w,linestyle="-",linewidth=lwid,color=list_color03[1])




w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0_run3,alpha1_run3,vthe0,vthe1,k)

gamma=wave_equation_two_electrons_IA_growth_rate_solve_2(wpe,mu,alpha0_run3,alpha1_run3,vd0_run3,vd1_run3,vthe0,vthe1,k,w)

haxe.plot(k/norm_k,gamma/norm_w,linestyle="-",linewidth=lwid,color=list_color03[2])


w = wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0_run1,alpha1_run1,vthe0,vthe1,k)

gamma=wave_equation_two_electrons_IA_growth_rate_solve_2(wpe,mu,alpha0_run1,alpha1_run1,vd0_run1,vd1_run1,vthe0,vthe1,k,w)

haxe.plot(k/norm_k,gamma/norm_w,linestyle="-",linewidth=lwid,color=list_color03[0])





haxe.axhline(y=0.0,linestyle="--",linewidth=0.5,color="k")

haxe.grid(linestyle="--",linewidth=0.8,color="w")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)


haxe.set_ylim(-0.0065,0.0005)
haxe.set_yticks(numpy.arange(-0.006,0.0+0.002,0.002))



xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)



if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\gamma/\omega_{pe}$",fontsize=20)




############################################################
#plt.show()

fig_path = path_fig
fig_name = os.path.splitext(os.path.basename(__file__))[0]


yaoxpy.fig_save(plt,figpath=fig_path,figname=fig_name,extension=fig_fmts)