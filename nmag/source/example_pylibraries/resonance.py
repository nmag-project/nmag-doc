import nmag, math, scipy.optimize
from nmag import SI, every

A_mag = SI(0.05e6, "A/m") # amplitude

sim = nmag.Simulation()

# NOTE: this simulation makes little sense if we do not specify the material-dependent
# damping constant! (But neither does it without H_demag!)

Py = nmag.MagMaterial(name = 'Py',
                      Ms = SI(1e6, 'A/m'),
                      llg_damping=0.02,
                      exchange_coupling = SI(13.0e-12, 'J/m'))

sim.load_mesh('rod.nmesh.h5',[('rod', Py)],unit_length=SI(1e-9,'m'))

r=open("resonance.log","w")
r2=open("resonance-plot.log","w")
r2.write("import pylab\n\ndata=[\n")

def simulate_resonance(freq,
                       amplitude=A_mag,
                       dt=SI(1e-12,"s"), # 1 ps time step is reasonable here.
                       time_steps=6000
                       # We simulate six full oscillations at 1 GHz
                       ):
    sim.set_m([0,0,1])
    angle_step=(2*math.pi*freq*dt).value
    result=[]
    for i in range(time_steps):
        h_ext=[amplitude.value*math.sin(i*angle_step),0,0]
        r.write("T=%s h_ext=%s" %(str(i*dt),str(h_ext)))
        sim.set_H_ext(h_ext,SI(1,"A/m"))
        
        sim.advance_time(i*dt)
        avg=sim.get_subfield_average_siv("M","Py") # XXX Why is this not in the manual?
        # XXX Note: counter-intuitive: get_subfield_average_siv("M","M_Py") does not work!
        r.write(" AVG=%s\n"% repr(avg))
        r2.write(" [%f,%f,%f],\n" %(i,avg[0],avg[1]))
        r2.flush()
        result.append((i*dt,h_ext,avg))
    return result
                       
def amplitude(freq_GHz):
    a_data=simulate_resonance(SI(freq_GHz*1e9,"1/s"))
    o_data=[[a[0].value,a[2][0],a[2][1],a[2][2]] for a in a_data]
    (fit_freq,fit_params)=nmag.fit_oscillation(o_data[5000:]) # skip lead-in
    a=math.sqrt(sum([p[1]*p[1] for p in fit_params]))
    print "Freq: %.2f GHz -- Amplitude: %f"%(freq_GHz,a)
    r.write("Freq: %.2f GHz -- Amplitude: %f\n"%(freq_GHz,a)) # DDD
    r.flush() # DDD
    r2.write("  ]\n\n\n")
    r2.flush()
    return -a # we minimize the *negative* amplitude here!

resonance=scipy.optimize.fmin(amplitude,[2.0],xtol=0.05,ftol=0.05)

print "Resonance frequency: %.1f GHz" % resonance

