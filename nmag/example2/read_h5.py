import nmag

#read data, positions, and sites from h5 file
m=nmag.get_subfield_from_h5file('bar_dat.h5','m_Py',id=0)
pos=nmag.get_subfield_positions_from_h5file('bar_dat.h5','m_Py')
site=nmag.get_subfield_sites_from_h5file('bar_dat.h5','m_Py')

#can carry out some sanity checks (but is not necessary)
assert m.shape == pos.shape
assert len(m) == len(site)

#print the data
for i in range(len(m)):
    print site[i], pos[i], m[i]
