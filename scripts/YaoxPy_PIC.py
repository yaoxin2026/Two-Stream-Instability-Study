import os

import h5py

import numpy

############################################################
############################################################
import shutil


def pic_parameter_read(path_simu_data,path_data):
   
    print(os.path.join(path_data,"parameters.h5"))
    
    if not os.path.exists(os.path.join(path_data,"parameters.h5")):
       #shutil.copyfile(os.path.join(os.path.dirname(path_simu_data),"pic2d"),os.path.dirname(path_data))
       shutil.copyfile(os.path.join(path_simu_data,"parameters.h5"),os.path.join(path_data,"parameters.h5"))


    H5FILE_R=h5py.File(os.path.join(path_data,"parameters.h5"),"r")

    list_parameters = {}

    list_var=["dx","dt","wpe","wce"]
    for var_name in list_var:
        var_value=H5FILE_R["/parameter/%s"%(var_name)][()]
        list_parameters.update({"%s"%(var_name):var_value})
    
    list_var=["ndimension","nx","ny","ncpux","ncpuy","nspecies"]
    for var_name in list_var:
        var_value=H5FILE_R["/parameter/%s"%(var_name)][()]
        list_parameters.update({"%s"%(var_name):int(var_value)})

    
    if int(list_parameters["ndimension"])==2:
       ncpu=int(list_parameters["ncpux"])*int(list_parameters["ncpuy"])
    elif int(list_parameters["ndimension"])==3:
       ncpu=int(list_parameters["ncpux"])*int(list_parameters["ncpuy"])*int(list_parameters["ncpuz"])
       
    list_parameters.update({"ncpu":ncpu})
    
    nspecies=int(list_parameters["nspecies"])
    
    list_var=["vd","vth","m","macro","ppc","q"]
    for sid in range(nspecies):
        for var_name in list_var:
            var_value=H5FILE_R["/parameter/species_%d/%s"%(sid,var_name)][()]
            list_parameters.update({"%s%d"%(var_name,sid):var_value})
    
    H5FILE_R.close()
    
    #print("list_parameters =",list_parameters["vd0"])
    
    return list_parameters








############################################################
############################################################


def dir_mk(path0,path_sub=None):
    if path_sub is None:
       filepath=path0
    else:
       filepath=os.path.join(path0,path_sub)
    
    try:
        if not os.path.exists(filepath):
           os.makedirs(filepath)
    except OSError:
        pass
       
    return filepath


############################################################
############################################################

def arrow3d(haxe,coor1,coor2,col='r',lines='-',linew=2.0,lent=1.0,ra=0.08):
    coor1=numpy.array(coor1)
    coor2=numpy.array(coor2)
    haxe.quiver(coor1[0],coor1[1],coor1[2],coor2[0]-coor1[0],coor2[1]-coor1[1],coor2[2]-coor1[2],length=lent,arrow_length_ratio=ra,color=col,linewidth=linew,linestyle=lines)


############################################################
# figsave
############################################################
# this_file=__file__
# figname : without extension


list_fig_extension = [".png",".pdf",".svg",".eps",".tiff"]


def fig_save(plt,figpath=None,figname=None,this_file=None,extension=".png",dpi=300,if_info=True):
    
    if figpath is None:
       print("Error : argument \'figpath\' is default !")
       exit(1)
    
    if figname is None and this_file is None:
       print("Error : one of arguments \'figname\' and \'this_file\' must be not default !")
       exit(1)

    print("*"*65)
    print("plotting ...")

    if not os.path.exists(figpath):
       os.makedirs(figpath)
    
    if figname is None and not this_file is None:
       figname=os.path.splitext(os.path.basename(this_file))[0]
    
    

    if fig_extension(extension)==0:
       
       if figname[-4:] in list_fig_extension:
          figname=figname[:-4]
       
       figname_tmp=figname+extension

       plt.savefig(os.path.join(figpath,figname_tmp),dpi=dpi)

       if if_info:
          #print("*"*65)
          print("save to  : path = ",figpath)
          print("           file = ",figname_tmp)
          #print("*"*65)

    else:
       if if_info:
          #print("*"*65)
          print("save to  : path = ",figpath)

       for exttmp in extension:
           if fig_extension(extension)==1:
              figname_tmp=figname+exttmp
           elif fig_extension(extension)==2:
              figname_tmp=figname+exttmp[0]
              dpi=exttmp[1]
           
           plt.savefig(os.path.join(figpath,figname_tmp),dpi=dpi)
    
           if if_info:
              print("           file = ",figname_tmp)


    print("*"*65)
       
    #plt.clf()
    plt.close()




def fig_extension(extension):
    
    if isinstance(extension,str):
       if extension in list_fig_extension:
          return 0
       else:
          print("Error : \"extension\" must be one in",list_fig_extension)
          exit(1)
    elif isinstance(extension,list):
       if isinstance(extension[0],str):
          for exttmp in extension:
              if not (exttmp in list_fig_extension):
                 print("Error : \"extension\" must be one in",list_fig_extension)
                 exit(1)
          
          return 1
       elif isinstance(extension[0],list):

          for exttmp in extension:
              if not ((isinstance(exttmp,list) and len(exttmp)==2) and (isinstance(exttmp[0],str) and isinstance(exttmp[1],int))):
                 print("Error : \"extension\" must be in formats as [[\".png\",300],[\".pdf\",300],...] with 300 as dpi !")
                 exit(1)
          
          return 2




import matplotlib as mpl

def colors_generate(N,cmap=mpl.cm.rainbow,mode=None):
    # cmap = plt.colormaps["rainbow"]
    # cmap = mpl.cm.jet,mpl.cm.rainbow,mpl.cm.nipy_spectral
    colors_tmp = cmap(numpy.linspace(0.0,1.0,N))

    if mode=="rgb":
       colors_tmp = [mpl.colors.to_rgb(tmp) for tmp in colors_tmp]
    elif mode=="rgba":
       colors_tmp = [mpl.colors.to_rgba(tmp) for tmp in colors_tmp]
    
    return colors_tmp



def data_zero_replace(data):
    index=numpy.where(data==0)
    if numpy.sum(index):
       data[index]=data.min()+1e-32
    return data


############################################################
############################################################
print("*"*65)
print("*"*65)



