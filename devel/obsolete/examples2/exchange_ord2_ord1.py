
import nfem, nmesh, math, sys

# One problem with the computation of the exchange field is that we do
# not sample the initial magnetization configuration properly.
# (Actually, we should do a least-squares fit, rather than just
# sampling at the support points.)
#
# One thing we can do in order to check for the importance of this
# effect is to sample M at second order, and translate that information
# back to first order by means of a second-order-M -> first-order-M
# projector.


nfem.set_default_dimension(3)
nfem.set_default_order(1)


bar = nmesh.box([-5.0,-2.0,-2.0],[5.0,2.0,2.0])


the_mesh = nmesh.mesh(objects = [bar],
                      cache_name="bar-mesh3",
                      a0=0.7,
                      bounding_box=[[-6.0,-3.0,-3.0],[6.0,3.0,3.0]],
                      neigh_force_scale = 1.,
                      # density = density,
                      initial_settling_steps = 50,
                      max_relaxation = 4,
                      # callback=(my_function, N),
                      # max_steps=677
                      max_steps=500
                      )

nfem.set_default_mesh(the_mesh)

element_M = nfem.make_element("M",[3])
element_M2 = nfem.make_element("M2",[3],ord=2)
element_H = nfem.make_element("H",[3])

mwe_M      = nfem.make_mwe("mwe_M", [(1,element_M)])
mwe_M2      = nfem.make_mwe("mwe_M2", [(1,element_M2)])
mwe_H      = nfem.make_mwe("mwe_H", [(1,element_H)])

diffop_v_laplace = nfem.diffop("-<d/dxj H(k) || d/dxj M(k)>, j:3, k:3")

print "Made MWEs"
sys.stdout.flush()

# Note that this magnetization is "mostly zero":
def fun_M0(dof_name_indices,position):
    if dof_name_indices[1][0]==1: # y-direction
        x=position[0]
        return 1.0/(1.0+(x*x/4.0))
    else:
        return 0.0

# MMM

field_M2=nfem.make_field(mwe_M2,fun_M0)

print "Made M2"
sys.stdout.flush()


diffop_m_m2 = nfem.diffop("<M(j) || M2(j)>, j:3")
prematrix_m_m2=nfem.prematrix(diffop_m_m2,mwe_M,mwe_M2)

print "Made M2->M prematrix"
sys.stdout.flush()


m_m2=nfem.prematrix_applicator(prematrix_m_m2)

field_M=nfem.cofield_to_field(m_m2(field_M2))

print "Have field_M"
sys.stdout.flush()

# MMM


prematrix_v_laplace=nfem.prematrix(diffop_v_laplace,mwe_H,mwe_M)

v_laplace=nfem.prematrix_applicator(prematrix_v_laplace)

field_H=nfem.cofield_to_field(v_laplace(field_M))

range_i=100

# nfem.field_print_contents(field_H)

for i in range(0,range_i):
    pos=[-4.5+9.0*i/(0.0+range_i),0.0,0.0]
    field_val=nfem.probe_field(field_H,"H",pos)
    print pos[0], " ", field_val[1]


#for i in range(0,range_i):
#    pos=[-4.5+9.0*i/(0.0+range_i),0.0,0.0]
#    print pos, "=> ", nfem.probe_field(field_H,"H",pos)
