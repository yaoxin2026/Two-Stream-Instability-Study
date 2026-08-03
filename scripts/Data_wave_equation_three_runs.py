import os
import numpy

from YaoxPy_Wave_Equations import (
    CGS,
    save_solution_hdf5,
    solve_six_full_maxwellian_branches,
    solve_six_reduced_branches,
)


from YaoxPy_Import_CWD import (
    path_data
)

############################################################
############################################################


mu  = 1836.0

wpe = 5.0e9
wpi = wpe / numpy.sqrt(mu)

vthe0 = 0.03 * CGS["c"]
vthe1 = 0.03 * CGS["c"]
vthi  = vthe0 / numpy.sqrt(mu)


norm_k = wpe / CGS["c"]
norm_w = wpe


kanchor = 5.0

kmax    = 20.0
dk      = 0.002
N       = int(numpy.ceil(kmax/dk))
k       = numpy.arange(-N,N+1, dtype=float)*dk*norm_k


#path_data = "."



############################################################
############################################################

########## Run1
alpha0 = 0.5
alpha1 = 0.5

vd1  =  0.3*CGS["c"]
vd0  =  -vd1*alpha1/alpha0

wpe0 = wpe*numpy.sqrt(alpha0)
wpe1 = wpe*numpy.sqrt(alpha1)


roots_reduced = solve_six_reduced_branches(k, wpe0, vd0, wpe1, vd1, wpi,)

roots_full    = solve_six_full_maxwellian_branches(k,
                                                   roots_reduced,
                                                   wpe0, vd0, vthe0,
                                                   wpe1, vd1, vthe1,
                                                   wpi, vthi,
                                                   anchor_abs_K=kanchor)


save_solution_hdf5(os.path.join(path_data,f"data_wave_equation_reduced_Run1.h5"),roots_reduced)
save_solution_hdf5(os.path.join(path_data,f"data_wave_equation_full_maxwellian_Run1.h5"),roots_full)


########## Run2
alpha0 = 4.0/5.0
alpha1 = 1.0/5.0

vd1  =  0.3*CGS["c"]
vd0  =  -vd1*alpha1/alpha0

wpe0 = wpe*numpy.sqrt(alpha0)
wpe1 = wpe*numpy.sqrt(alpha1)

roots_reduced = solve_six_reduced_branches(k, wpe0, vd0, wpe1, vd1, wpi,)

roots_full    = solve_six_full_maxwellian_branches(k,
                                                   roots_reduced,
                                                   wpe0, vd0, vthe0,
                                                   wpe1, vd1, vthe1,
                                                   wpi, vthi,
                                                   anchor_abs_K=kanchor)

save_solution_hdf5(os.path.join(path_data,f"data_wave_equation_reduced_Run2.h5"),roots_reduced)
save_solution_hdf5(os.path.join(path_data,f"data_wave_equation_full_maxwellian_Run2.h5"),roots_full)



########## Run3
alpha0 = 20.0/21.0
alpha1 = 1.0/21.0

vd1  =  0.3*CGS["c"]
vd0  =  -vd1*alpha1/alpha0

wpe0 = wpe*numpy.sqrt(alpha0)
wpe1 = wpe*numpy.sqrt(alpha1)

roots_reduced = solve_six_reduced_branches(k, wpe0, vd0, wpe1, vd1, wpi,)

roots_full    = solve_six_full_maxwellian_branches(k,
                                                   roots_reduced,
                                                   wpe0, vd0, vthe0,
                                                   wpe1, vd1, vthe1,
                                                   wpi, vthi,
                                                   anchor_abs_K=kanchor)

save_solution_hdf5(os.path.join(path_data,f"data_wave_equation_reduced_Run3.h5"),roots_reduced)
save_solution_hdf5(os.path.join(path_data,f"data_wave_equation_full_maxwellian_Run3.h5"),roots_full)
