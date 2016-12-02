import nmag

m=nmag.get_subfield_from_h5file('bar_dat.h5','m_py',row=0)

pos=nmag.get_subfield_positions_from_h5file('bar_dat.h5','m_Py')

site=nmag.get_subfield_sites_from_h5file('bar_dat.h5','m_Py')

assert m.shape == pos.shape

assert len(m) == len(site)

for i in range(len(m)):
    print site[i], pos[i], m[i]
