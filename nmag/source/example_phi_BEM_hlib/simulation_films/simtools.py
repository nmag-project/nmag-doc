import nmag
from nmag import SI, every

import commands, os, time 


def get_memory_of_process(pid=None):
    """Retrieves the memory requirements of the process (in Kilobyte) 
    specified by pid by extracting them using the command \"ps ux\". 
    If no pid is given, the current process will be used.

    """

    if pid==None:
        pid=os.getpid()

    output=commands.getoutput("ps ux")
    output=output.split('\n')
    #some checks whether 'ps ux' gives the expected output
    headline = output[0].split()
    assert 'RSS' in headline
    assert headline.index('RSS')==5
    for i in range(1,len(output)):
        line=output[i].split()
        if int(line[1])==pid:
            memory_rss = float(line[5])
            break

    return memory_rss


def run_simulation(simname,meshfile,hmatrix=None,tol=0.000001):
    
    T_start = time.time()
    
    sim = nmag.Simulation(name=simname,phi_BEM=hmatrix)
    Ni = nmag.MagMaterial(name="Ni",
                          Ms=SI(493380, "A/m"),
                          exchange_coupling=SI(7.2e-12, "J/m"),
                          llg_damping=1.0)
    ps = SI(1e-12,"s")
    sim.load_mesh(meshfile, [("thinfilm", Ni)], unit_length=SI(1e-9, "m"))
    sim.set_m([0.05,0.02,1.00])
    sim.set_params(ts_abs_tol=tol,ts_rel_tol=tol)
    sim.set_H_ext([0,0,0], SI('A/m'))
    dt = SI(5e-12, "s")

    T_end = time.time()
    T_setup = T_end - T_start

    T_start = time.time()
    sim.relax(save=[('averages',every('time', 10*ps))])
    T_end = time.time()
    T_sim = T_end - T_start
    sim.save_data(['m'])

    number_nodes = len(sim.mesh.points)
    mem_rss =  get_memory_of_process() / 1024.0
    

    return [number_nodes,mem_rss,T_setup,T_sim]

