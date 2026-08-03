from YaoxPy_Import_CWD import *

############################################################
############################################################



gIA_run1 = growth_rate_ion_acoustic(k,wpe0_run1,vd0_run1,vthe0,wpi,vthi)

gIA_run2 = growth_rate_ion_acoustic(k,wpe0_run2,vd0_run2,vthe0,wpi,vthi)

gIA_run3 = growth_rate_ion_acoustic(k,wpe0_run3,vd0_run3,vthe0,wpi,vthi)

glangmuir_run3 = growth_rate_langmuir(k,wpe0_run3,vd0_run3,vthe1)



############################################################
############################################################

xmin,xmax=-21,21
xticks=numpy.arange(-20,20+5,5)




#list_colors=yaoxpy_vis.colors_generate(8,cmap=mpl.cm.rainbow)

list_colors  = ["r","r","m","green","lime","coral"]

list_colors1 = ["r","r","m","green","green","m"]



lwid_0 = 1.0
lwid_1 = 0.9

############################################################
############################################################

hfig = plt.figure(figsize=(13,11))


margin=[0.09,0.03,0.05,0.07,0.12,0.04]
barbox=[0.05,0.01,0.8]
windows=[5,2]
size=[1,1]


############################################################
############################################################
rid,cid = 0,0

axes_pos = [ 0.0900,  0.8060,  0.3800,  0.1440]

#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="white")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))


# reduced cold limit

lsty = "--"


km = k_reduced_run1
kb = kb_reduced_run1


for m in [0,3,4,2,5]:

    wm_real = numpy.copy(wroots_real_reduced_run1[:,m])
    wm_imag = numpy.copy(wroots_imag_reduced_run1[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_real)
    index       = wtmp<=0
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid_1,color=list_colors1[m])



# full maxwellian

list_modes   = [f"beam\\textendash beam",
                f"beam\\textendash beam",
                f"ion\\textendash inertia",
                f"beam\\textendash like",
                f"beam\\textendash like",
                f"ion\\textendash inertia"]



lsty = "-"


km = k_run1
kb = kb_run1


for m in [0,3,5]:

    wm_real = numpy.copy(wroots_real_run1[:,m])
    wm_imag = numpy.copy(wroots_imag_run1[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_real)
    index       = wtmp<=0
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid_1,color=list_colors1[m],label=r"$%s$"%(list_modes[m]))



for m in [2,4]:

    wm_real = numpy.copy(wroots_real_run1[:,m])
    wm_imag = numpy.copy(wroots_imag_run1[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_real)
    index       = wtmp<=0
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid_1,color=list_colors1[m])

haxe.axvline(kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")
haxe.axvline(-kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")



haxe.legend(loc="upper center",frameon=False,fontsize=12)


haxe.grid(linestyle="--",linewidth=0.05,color="gainsboro")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.2,8.2)
haxe.set_yticks(numpy.arange(0,8+1,1))



xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)Run1$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)


haxe.set_title(r"$Dispersion\ Relation$")



if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega_r/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\omega_i/\omega_{pe}$",fontsize=20)




############################################################
############################################################
rid,cid = 0,1
axes_pos = [ 0.5900,  0.8060,  0.3800,  0.1440]
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="white")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))


# full maxwellian

lsty = "-"


km = k_run1
kb = kb_run1


print("kb = ",kb)

for m in [0]:

    wm_real = numpy.copy(wroots_real_run1[:,m])
    wm_imag = numpy.copy(wroots_imag_run1[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_imag)
    index       = wm_real<=0
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan
    
    if m==0:
       lwid=lwid_0
    else:
       lwid=lwid_1

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid,color=list_colors1[m])



for m in [3,4,2,5]:

    wm_real = numpy.copy(wroots_real_run1[:,m])
    wm_imag = numpy.copy(wroots_imag_run1[:,m])

    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_imag)
    index       = numpy.abs(ktmp)>=kb[0]
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle="--",linewidth=lwid_1,color="grey")
    

    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_imag)
    index       = numpy.abs(ktmp)<=kb[0]
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan

    if m==0:
       lwid=lwid_0
    else:
       lwid=lwid_1

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid,color=list_colors1[m])






haxe.axvline(kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey",label=r"$|k_b|=%.2fd_e^{-1}$"%(kb[0]/norm_k))
haxe.axvline(-kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")



haxe.legend(loc="upper right",frameon=False,fontsize=12)


haxe.grid(linestyle="--",linewidth=0.05,color="gainsboro")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.8,0.4)
haxe.set_yticks(numpy.arange(-0.8,0.4+0.2,0.2))



xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)



haxe.set_title(r"$Growth\ Rate$")




if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega_r/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\omega_i/\omega_{pe}$",fontsize=20)






############################################################
############################################################
rid,cid = 1,0
axes_pos = [ 0.0900,  0.6220,  0.3800,  0.1440]
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="white")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))


# reduced cold limit

lsty = "--"


km = k_reduced_run2
kb = kb_reduced_run2


for m in [0,3,4,2,5]:

    wm_real = numpy.copy(wroots_real_reduced_run2[:,m])
    wm_imag = numpy.copy(wroots_imag_reduced_run2[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_real)
    index       = wtmp<=0
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid_1,color=list_colors[m])



# full maxwellian

list_modes   = [f"beam\\textendash beam",
                f"beam\\textendash beam",
                f"hybrid",
                f"beam\\textendash like",
                f"beam\\textendash like",
                f"ion\\textendash inertia"]




lsty = "-"


km = k_run2
kb = kb_run2


for m in [0,3,4,2,5]:

    wm_real = numpy.copy(wroots_real_run2[:,m])
    wm_imag = numpy.copy(wroots_imag_run2[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_real)
    index       = wm_real<=0
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid_1,color=list_colors[m],label=r"$%s$"%(list_modes[m]))



haxe.axvline(kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")
#haxe.axvline(-kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")



haxe.legend(loc="upper center",frameon=False,fontsize=12)


haxe.grid(linestyle="--",linewidth=0.05,color="gainsboro")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.2,8.2)
haxe.set_yticks(numpy.arange(0,8+1,1))



xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)Run2$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)


#haxe.set_title(r"$Dispersion\ Relation$")



if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega_r/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\omega_i/\omega_{pe}$",fontsize=20)




############################################################
############################################################
rid,cid = 1,1
axes_pos = [ 0.5900,  0.6220,  0.3800,  0.1440]
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="white")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))





# full maxwellian

lsty = "-"


km = k_run2
kb = kb_run2


print("kb = ",kb)


for m in [0,3,4,2,5]:

    wm_real = numpy.copy(wroots_real_run2[:,m])
    wm_imag = numpy.copy(wroots_imag_run2[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_imag)
    index       = wm_real<=0
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan

    if m==0:
       lwid=lwid_0
    else:
       lwid=lwid_1

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid,color=list_colors[m])




haxe.axvline(kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey",label=r"$k_b=%.2fd_e^{-1}$"%(kb[0]/norm_k))
#haxe.axvline(-kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")



haxe.legend(loc="upper right",frameon=False,fontsize=12)


haxe.grid(linestyle="--",linewidth=0.05,color="gainsboro")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.8,0.4)
haxe.set_yticks(numpy.arange(-0.8,0.4+0.2,0.2))



xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)



#haxe.set_title(r"$Growth\ Rate$")




if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega_r/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\omega_i/\omega_{pe}$",fontsize=20)










############################################################
############################################################
rid,cid = 2,0
axes_pos = [ 0.0900,  0.4380,  0.3800,  0.1440]
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="white")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))


# reduced cold limit

lsty = "--"


km = k_reduced_run3
kb = kb_reduced_run3


for m in [0,3,4,2,5]:

    wm_real = numpy.copy(wroots_real_reduced_run3[:,m])
    wm_imag = numpy.copy(wroots_imag_reduced_run3[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_real)
    index       = wm_real<=0
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid_1,color=list_colors[m])



# full maxwellian

list_modes   = [f"beam\\textendash beam",
                f"beam\\textendash beam",
                f"hybrid",
                f"beam\\textendash like",
                f"ion\\textendash inertia",
                f"ion\\textendash inertia"]




lsty = "-"


km = k_run3
kb = kb_run3


for m in [0,3,2,4,5]:

    wm_real = numpy.copy(wroots_real_run3[:,m])
    wm_imag = numpy.copy(wroots_imag_run3[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_real)
    index       = wm_real<=0
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid_1,color=list_colors[m],label=r"$%s$"%(list_modes[m]))



haxe.axvline(kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")
#haxe.axvline(-kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")



haxe.legend(loc="upper center",frameon=False,fontsize=12)


haxe.grid(linestyle="--",linewidth=0.05,color="gainsboro")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.2,8.2)
haxe.set_yticks(numpy.arange(0,8+1,1))



xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)Run3$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)


#haxe.set_title(r"$Dispersion\ Relation$")



if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega_r/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\omega_i/\omega_{pe}$",fontsize=20)




############################################################
############################################################
rid,cid = 2,1
axes_pos = [ 0.5900,  0.4380,  0.3800,  0.1440]
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="white")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))





# full maxwellian

lsty = "-"


km = k_run3
kb = kb_run3


print("kb = ",kb)


for m in [0,3,4,2,5]:

    wm_real = numpy.copy(wroots_real_run3[:,m])
    wm_imag = numpy.copy(wroots_imag_run3[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_imag)
    index       = wm_real<=0
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan

    if m==0:
       lwid=lwid_0
    else:
       lwid=lwid_1

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid,color=list_colors[m])




haxe.axvline(kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey",label=r"$k_b=%.2fd_e^{-1}$"%(kb[0]/norm_k))
#haxe.axvline(-kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")



haxe.legend(loc="upper right",frameon=False,fontsize=12)


haxe.grid(linestyle="--",linewidth=0.05,color="gainsboro")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.8,0.4)
haxe.set_yticks(numpy.arange(-0.8,0.4+0.2,0.2))



xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)



#haxe.set_title(r"$Growth\ Rate$")




if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega_r/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\omega_i/\omega_{pe}$",fontsize=20)







############################################################
############################################################
rid,cid = 3,0
axes_pos = [ 0.0900,  0.2540,  0.3800,  0.1440]
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="white")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))





# full maxwellian

list_modes   = [f"beam\\textendash beam",
                f"beam\\textendash beam",
                f"hybrid",
                f"beam\\textendash like",
                f"ion\\textendash inertia",
                f"ion\\textendash inertia"]



lwid = lwid_1
lsty = "-"


km = k_run3
kb = kb_run3


for m in [0,2]:

    wm_real = numpy.copy(wroots_real_run3[:,m])
    wm_imag = numpy.copy(wroots_imag_run3[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_real)
    #index       = wm_real<0
    #ktmp[index] = numpy.nan
    #wtmp[index] = numpy.nan

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid,color=list_colors[m],label=r"$%s$"%(list_modes[m]))



haxe.axvline(kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")
#haxe.axvline(-kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")



haxe.legend(loc="upper center",frameon=False,fontsize=12)






# Langmuir


haxe.plot(k/norm_k,wL_run3/norm_w,linestyle="--",linewidth=lwid_0,color="r",label=r"$Langmuir$")



haxe.plot(k/norm_k,wBeam_run3/norm_w,linestyle="--",linewidth=lwid_1,color="grey",label=r"$\omega=v_{d2}k$")





haxe.legend(loc="upper right",frameon=False,fontsize=12)



haxe.grid(linestyle="--",linewidth=0.05,color="gainsboro")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.1,3.1)
haxe.set_yticks(numpy.arange(0,3+1,1))



xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)Run3$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)



#haxe.set_title(r"$Dispersion\ Relation$")



if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega_r/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\omega_i/\omega_{pe}$",fontsize=20)







############################################################
############################################################
rid,cid = 3,1
axes_pos = [ 0.5900,  0.2540,  0.3800,  0.1440]
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="white")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))




# full maxwellian

lwid = lwid_1
lsty = "-"


km = k_run3
kb = kb_run3


print("kb = ",kb)


for m in [0,2]:

    wm_real = numpy.copy(wroots_real_run3[:,m])
    wm_imag = numpy.copy(wroots_imag_run3[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_imag)
    index       = wm_real<0
    ktmp[index] = numpy.nan
    wtmp[index] = numpy.nan

    if m==0:
       lwid=lwid_0
    else:
       lwid=lwid_1

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid,color=list_colors[m])




haxe.axvline(kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey",label=r"$k_b=%.2fd_e^{-1}$"%(kb[0]/norm_k))
#haxe.axvline(-kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")



haxe.legend(loc="upper right",frameon=False,fontsize=12)


#haxe.plot(k/norm_k,glangmuir_run3/norm_w,linestyle="-",linewidth=lwid,color="r")





haxe.grid(linestyle="--",linewidth=0.05,color="gainsboro")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.8,0.4)
haxe.set_yticks(numpy.arange(-0.8,0.4+0.2,0.2))


xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)


#haxe.set_title(r"$Growth\ Rate$")




if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega_r/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\omega_i/\omega_{pe}$",fontsize=20)














############################################################
############################################################
rid,cid = 4,0
axes_pos = [ 0.0900,  0.0700,  0.3800,  0.1440]
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="white")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))




# full maxwellian

list_modes   = [f"beam\\textendash beam",
                f"beam\\textendash beam",
                f"hybrid",
                f"beam\\textendash like",
                f"ion\\textendash inertia",
                f"ion\\textendash inertia"]



lwid = lwid_0
lsty = "-"


km = k_run3
kb = kb_run3


for m in [4,5]:

    wm_real = numpy.copy(wroots_real_run3[:,m])
    wm_imag = numpy.copy(wroots_imag_run3[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_real)
    #index       = wm_real<0
    #ktmp[index] = numpy.nan
    #wtmp[index] = numpy.nan

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid,color=list_colors[m],label=r"$%s$"%(list_modes[m]))



#haxe.axvline(kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")
#haxe.axvline(-kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")



haxe.legend(loc="upper center",frameon=False,fontsize=12)






# ion-acoustic



lwid = lwid_1


haxe.plot(k/norm_k,wIA_run3/norm_w,linestyle="--",linewidth=lwid,color="b",label=r"$ion\textendash acoustic$")


haxe.plot(k/norm_k,-1.0*wIA_run3/norm_w,linestyle="--",linewidth=lwid,color="b")


#haxe.plot(k/norm_k,wIA_linear_run3/norm_w,linestyle="--",linewidth=lwid,color="b",label=r"$\omega=c_sk$")





haxe.legend(loc="upper right",frameon=False,fontsize=12)



haxe.grid(linestyle="--",linewidth=0.05,color="gainsboro")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)

haxe.set_ylim(-0.015,0.015)
haxe.set_yticks(numpy.arange(-0.015,0.014+0.005,0.005))



xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)Run3$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)



#haxe.set_title(r"$Dispersion\ Relation$")



if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega_r/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\omega_i/\omega_{pe}$",fontsize=20)




############################################################
############################################################
rid,cid = 4,1
axes_pos = [ 0.5900,  0.0700,  0.3800,  0.1440]
#axes_pos=yaoxpy_vis.fig_axes_position(index=[rid,cid],size=[1,1],windows=windows,margin=margin)
haxe=hfig.add_axes(axes_pos,facecolor="white")

print("rid,cid  = (%d, %d)"%(rid,cid))
print("axes_pos = [%7.4f, %7.4f, %7.4f, %7.4f]"%(axes_pos[0],axes_pos[1],axes_pos[2],axes_pos[3]))





# full maxwellian

lsty = "-"


km = k_run3
kb = kb_run3


print("kb = ",kb)


for m in [4,5]:

    wm_real = numpy.copy(wroots_real_run3[:,m])
    wm_imag = numpy.copy(wroots_imag_run3[:,m])


    ktmp        = numpy.copy(km)
    wtmp        = numpy.copy(wm_imag)
    #index       = wm_real<0
    #ktmp[index] = numpy.nan
    #wtmp[index] = numpy.nan

    if m==0:
       lwid=lwid_0
    else:
       lwid=lwid_1

    haxe.plot(ktmp/norm_k,wtmp/norm_w,linestyle=lsty,linewidth=lwid,color=list_colors[m])


#haxe.axvline(kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey",label=r"$|k_b|=%.2fd_e^{-1}$"%(kb[0]/norm_k))
#haxe.axvline(-kb[0]/norm_k,linestyle="--",linewidth=0.5,color="grey")




#haxe.plot(k/norm_k,gIA_run3/norm_w,linestyle="-",linewidth=lwid,color="b")


#haxe.legend(loc="upper right",frameon=False,fontsize=12)



haxe.grid(linestyle="--",linewidth=0.05,color="gainsboro")

haxe.set_xlim(xmin,xmax)
haxe.set_xticks(xticks)


haxe.set_ylim(-0.0035,0.0035)
haxe.set_yticks(numpy.arange(-0.003,0.003+0.001,0.001))




xlim=haxe.get_xlim()
ylim=haxe.get_ylim()
haxe.text(xlim[0]+0.05*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),r"$(%s%d)$"%(chr(ord("a")+rid),cid+1),color="k",fontsize=16)


#haxe.set_title(r"$Growth\ Rate$")




if rid==(windows[0]-1):
   haxe.set_xlabel(r"$k\cdot d_e$",fontsize=20)

if cid==0:
   haxe.set_ylabel(r"$\omega_r/\omega_{pe}$",fontsize=20)
else:
   haxe.set_ylabel(r"$\omega_i/\omega_{pe}$",fontsize=20)











############################################################
############################################################
#plt.show()

fig_path = path_fig
fig_name = os.path.splitext(os.path.basename(__file__))[0]


yaoxpy_vis.fig_save(plt,figpath=fig_path,figname=fig_name,extension=fig_fmts)