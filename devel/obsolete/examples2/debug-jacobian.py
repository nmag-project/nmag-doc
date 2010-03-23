# (C) 2006 Dr. Thomas Fischbacher
#
# Debugging our Jacobian matrix

import nmag2 as nmag, nmesh as nm, sys, math, time
import nfem
# import nfem.visual

import ocaml #need this as long as we use ocaml.probe_field

from nmag2.si_units import SI

#nmag.set_global_log_level('debug')
#nmag.set_log_level('debug')

intensive_param_by_name={"H_x":0.0,"H_y":0.0,"H_z":0.0}
# Not applying any external field.

mat_Py = nmag.MagMaterial("Py",
                          Ms=SI(1e6,"A/m"),
                          J=SI(13.0e-12,"J/m"),
                          )

# We afterwards hack the parameters of that material definition
# that actually enter the LLG:

print "LLG length control: ", mat_Py.su_llg_gamma 
# THIS RETURNS: LLG length control:  1e-24
# That is completely bogus!
# BUGFIX REQUIRED!

mat_Py.su_llg_coeff1=  -0.5*1.0 # M x H
mat_Py.su_llg_coeff2= -0.02*0.0 # M x (M x H) -- NOTE that the coefficient MUST be negative here.
# mat_Py.su_llg_gamma=2.0 # length control
mat_Py.su_llg_gamma=0.2 # length control

sim=nmag.SimulationContext("debug-jacobian",do_demag=False)

# sim.timestepper_tuning_params=[1e-6,1e-6,2,300] # These are the defaults...
sim.timestepper_tuning_params=[1e-4,1e-4,2,300] # These are the defaults...

sim.load_mesh("./debug-jacobian.mesh", [("Py", mat_Py)])

print "DDD loaded mesh!"

def initial_magnetization(coords,mag_type):
    # Note: length is not one here, but the system will internally
    # normalize M in set_magnetization!
    print "Setting magnetization for coords=",coords," mag_type=",mag_type
    return [0.0,
            0.8*math.cos(coords[0]/6.0),
            0.6*math.sin(coords[0]/6.0),
            ]

sim.set_magnetization(initial_magnetization)
#sim.set_magnetization([0.0,0.0,1.0])

# XXX Note: this actually should not be necessary! get_M should make
# sure that a timestepper is being created!

sim._ensure_have_timestepper()

print "Field at x=4.5:",ocaml.probe_field(sim.fields["m"],"",[4.5])
print "M0: ",sim.get_m([0.5]) # XXX Note: returns None - this clearly is wrong!

time_start=time.time()

# nfem.visual.fields2vtkfile(sim.fields.values(),'jacobi-0.vtk',mesh=sim.mesh)


for i in range(1,10): # range(1,1000):
    target_time=sim.advance_time(intensive_param_by_name,10*0.05*i)
    print "t=%f Field at x=4.5:" % target_time,ocaml.probe_field(sim.fields["m"],"",[4.5])

    #write SI data into default file (results.ndt)
    sim.save_data_table()

    #write SU data into special file 
    sim.save_data_table_su(filename='results_su.ndt')

    # nfem.visual.fields2vtkfile(sim.fields.values(),'jacobi-%d.vtk' % i,mesh=sim.mesh)

# DDD ...now, at the end: "activate the omega-13"

sim.advance_time(intensive_param_by_name,-100.0)

# PROBLEMS:
#
# t=0.950000 Field at x=4.5: [('m_Py', [0.0, -0.070478999969441875, -0.023382641308906715])]
#
# Time-stepping makes magnetization slowly vanish. Actually, the last term in our extended LLG
# should prevent that!

""" Just documenting some results:

mat_Py.su_llg_coeff1=  -0.5 # M x H
mat_Py.su_llg_coeff2= -0.02 # M x (M x H) -- NOTE that the coefficient MUST be negative here.
mat_Py.su_llg_gamma=2.0 # length control

This gives a Jacobi re-computation pattern of the following form:

tf@alpha:~/ocaml/examples$ grep Computed /tmp/run.log 
JACOBI T=0.000000 JT=-10000000000000000159028911097599180468360808563945281389781327557747838772170381060813469985856815104.000000: Computed!
JACOBI T=8.565703 JT=0.000000: Computed!
JACOBI T=10.856515 JT=8.565703: Computed!
JACOBI T=12.655843 JT=10.856515: Computed!
JACOBI T=15.313185 JT=12.655843: Computed!
JACOBI T=18.499731 JT=15.313185: Computed!
JACOBI T=18.550025 JT=18.499731: Computed!
JACOBI T=18.540595 JT=18.550025: Computed!
JACOBI T=21.514439 JT=18.540595: Computed!
JACOBI T=26.087441 JT=21.514439: Computed!
JACOBI T=26.147476 JT=26.087441: Computed!
JACOBI T=26.135904 JT=26.147476: Computed!
JACOBI T=27.140596 JT=26.135904: Computed!
JACOBI T=28.800988 JT=27.140596: Computed!
JACOBI T=30.099500 JT=28.800988: Computed!
JACOBI T=31.632170 JT=30.099500: Computed!
JACOBI T=33.186126 JT=31.632170: Computed!
JACOBI T=35.191287 JT=33.186126: Computed!
JACOBI T=39.131542 JT=35.191287: Computed!
JACOBI T=39.213749 JT=39.131542: Computed!
JACOBI T=39.203960 JT=39.213749: Computed!
JACOBI T=39.638673 JT=39.203960: Computed!
JACOBI T=39.629990 JT=39.638673: Computed!

mat_Py.su_llg_coeff1=  -0.0 # M x H
mat_Py.su_llg_coeff2= -0.02 # M x (M x H) -- NOTE that the coefficient MUST be negative here.
mat_Py.su_llg_gamma=2.0 # length control

JACOBI T=0.000000 JT=-10000000000000000159028911097599180468360808563945281389781327557747838772170381060813469985856815104.000000: Computed!
JACOBI T=25.068991 JT=0.000000: Computed!
JACOBI T=26.572215 JT=25.068991: Computed!
JACOBI T=26.339279 JT=26.572215: Computed!
JACOBI T=27.137482 JT=26.339279: Computed!
JACOBI T=28.739430 JT=27.137482: Computed!
JACOBI T=28.538401 JT=28.739430: Computed!
JACOBI T=29.358774 JT=28.538401: Computed!
JACOBI T=31.340416 JT=29.358774: Computed!
JACOBI T=31.136907 JT=31.340416: Computed!
JACOBI T=31.940376 JT=31.136907: Computed!
JACOBI T=33.555610 JT=31.940376: Computed!
JACOBI T=33.355238 JT=33.555610: Computed!
JACOBI T=34.156580 JT=33.355238: Computed!
JACOBI T=35.704616 JT=34.156580: Computed!
JACOBI T=35.501275 JT=35.704616: Computed!
JACOBI T=36.281759 JT=35.501275: Computed!
JACOBI T=37.687769 JT=36.281759: Computed!
JACOBI T=37.476843 JT=37.687769: Computed!
JACOBI T=38.291323 JT=37.476843: Computed!
JACOBI T=39.907470 JT=38.291323: Computed!
JACOBI T=39.707128 JT=39.907470: Computed!
JACOBI T=40.499010 JT=39.707128: Computed!
JACOBI T=41.634776 JT=40.499010: Computed!
JACOBI T=41.421761 JT=41.634776: Computed!
JACOBI T=42.997701 JT=41.421761: Computed!
JACOBI T=42.682996 JT=42.997701: Computed!
JACOBI T=43.356324 JT=42.682996: Computed!
JACOBI T=44.482758 JT=43.356324: Computed!
JACOBI T=44.895495 JT=44.482758: Computed!
JACOBI T=45.320437 JT=44.895495: Computed!
JACOBI T=45.771658 JT=45.320437: Computed!
JACOBI T=46.215794 JT=45.771658: Computed!
JACOBI T=46.653974 JT=46.215794: Computed!
JACOBI T=47.081220 JT=46.653974: Computed!
JACOBI T=47.529153 JT=47.081220: Computed!
JACOBI T=47.968532 JT=47.529153: Computed!
JACOBI T=48.407637 JT=47.968532: Computed!
JACOBI T=48.839384 JT=48.407637: Computed!
JACOBI T=49.285394 JT=48.839384: Computed!
JACOBI T=49.722679 JT=49.285394: Computed!


mat_Py.su_llg_coeff1=  -0.5 # M x H
mat_Py.su_llg_coeff2= -0.00 # M x (M x H) -- NOTE that the coefficient MUST be negative here.
mat_Py.su_llg_gamma=2.0 # length control

--> ocaml weird low level exception

"""
