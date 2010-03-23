#
# Copied from sphere3, and zeroed anisotropy and made initial magnetisation uniform.
#

#
# (C) 2006 Dr. Thomas Fischbacher
# Relaxation of the homogeneously magnetized sphere

import nmag2 as nmag, nmesh as nm, sys, math, time
import nfem
import nfem.visual

import ocaml #need this as long as we use ocaml.probe_field

from nmag2.si_units import SI

import nmag2.si_units as si_units

#produce more output for now
nmag.set_global_log_level('debug')
nmag.set_log_level('debug')

intensive_param_by_name={"H_x":0.0,"H_y":0.0,"H_z":0.0}
# very slightly pulling in x-direction

def test_anisotropy_energy(v): # direction is unit-normalized!
    print "DDD anisotropy dir=",v
    result=(1.0*v[0]*v[0]+2.0*v[1]*v[1]+3.5*v[2]*v[2])
    return result*0.0

mat_Py = nmag.MagMaterial("Py",
                          Ms=SI(1e6,"A/m"),
                          exchange_coupling=SI(13.0e-12,"J/m"),
                          #Ms=1.,      
                          #J=13.,      
                          anisotropy_order=2,
                          anisotropy=test_anisotropy_energy
                          )

sim=nmag.SimulationContext("sphere")

# sim.timestepper_tuning_params=[1e-6,1e-6,2,300] # These are the defaults...
sim.timestepper_tuning_params=[1e-6,1e-6,4,300]

sim.set_magnetization([0.0,1.0,0.0])
# ^ just to show we can place set_magnetization() where we want...

sim.defregion("Py", nm.ellipsoid([3.0,3.0,3.0]), mag_mat=mat_Py)

def initial_magnetization(coords,mag_type):
    return [1,0,0]
    return [0.8*math.sin(coords[0]/3.0)*1e6,
            0.8*math.cos(coords[0]/3.0)*1e6,
            0.6*1e6
            ]

sim.generate_mesh(([-5.0,-5.0,-5.0],[5.0,5.0,5.0]), # bounding box
                  a0=1.0,
                  max_steps=1200,
                  cache_name="sphere3"
                  )

#
#sim.set_magnetization(initial_magnetization)
#
#import ocaml
#print "Field at origin (initial):",ocaml.probe_field(sim.fields["m"],"",[0.0,0.0,0.0])
#
#for i in range(1,3):
#    target_time=sim.advance_time(intensive_param_by_name,0.05*i)
#    print "Field at origin (T=%f):"%target_time,ocaml.probe_field(sim.fields["m"],"",[0.0,0.0,0.0])

# Re-set magnetization:

sim.set_magnetization(initial_magnetization)

nfem.visual.fields2vtkfile(sim.fields.values(),'sphere-0.vtk',mesh=sim.mesh)

time_start=time.time()


for i in range(1,5):
    target_time=sim.advance_time(intensive_param_by_name,0.05*i)

    #write SI data into default file (results.ndt)
    sim.save_data_table()

    #write SU data into special file 
    sim.save_data_table_su(filename='results_su.ndt')

    nfem.visual.fields2vtkfile(sim.fields.values(),'sphere-%d.vtk' % i,mesh=sim.mesh)
    print "Field at origin (T=%f):"%target_time,sim.get_m([0.0,0.0,0.0])
time_end=time.time()

print "SIM RUN TIME: ",time_end-time_start," sec."


if True:
    f=sim._get_nfem_fields()
    print "Fields are",f



    print "probe magnetisation along line:"

    nm = 1e-9
    data = []
    for i in range(-10,11):
        x = float(i)/10*3*nm
        data.append(sim.get_M([x,0,0]))
    print data

    print "Demag field at [0,0,0] =", sim.get_H_demag([0,0,0])
    print "Exchange field at [0,0,0] =", sim.get_H_exchange([0,0,0])
    print "anisotropy field at [0,0,0] =", sim.get_H_anisotropy([0,0,0])

import sys
sys.exit(0)



su=nmag.simulation_units

from math import pi

a=SI(3e-9,"m") 
print "a=",a
print "a**1",a**1
print "a**2",a**2
print "a**3",a**3
print "a**4",a**4
print "a**0",a**0
vol = a**3*pi*4./3. #m^3

print "Volume=",vol," -> ",su.of(vol),"su"

demagfactor = -1./3. #theoretical value
demagfactor = -0.32815 #from simulation

Ms=SI(1e6,"A/m")

#Ms = 795774.72 #Example from Magpar, expect energydensity=132629.12 J/m^3

Hdem = Ms*demagfactor #A/m
print "Hdem",Hdem

mu0 = si_units.mu0 #correct value: 4*pi*1e-7  # kg m / A^2 s^2
print "correct mu0=",mu0
#mu0 = mu0/(4*pi/10)
#print "using broken mu0 (factor 4pi/10 missing): broken mu0=",mu0

B = mu0*Hdem     # kg m / (A^2 s^2) * A/m = kg / A s^2

print "B=",B

energy = -1*B*Ms*vol # kg/A s^2 * A/m * m^3 = kg m^2 / s^2

print "energy=",energy

print "Energy:",energy, "J"
print "energy density", energy / vol

print "Simulation units: energy:",su.of(energy), "su"
print "Simulation units: energy density", su.of(energy/vol) ,"su"

#print "to see energy density in si units, run this on results.ndt: gawk '{print $1 "\t" $20}' results.ndt"
