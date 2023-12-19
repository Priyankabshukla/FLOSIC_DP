import numpy as np
import dpdata
import ase
from ase.io import read,write
import os
import matplotlib.pyplot as plt
a=read(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config/Allconfig/POSCAR-137') ###Atoms object
print(a.get_positions())
b=a.get_positions()*0.529177
a.set_positions(b)
write(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config/test/frmorb-137',a,format='xyz')
print(a.get_positions())

fod=read(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config/test/FOD')
fod_b=fod.get_positions()*0.529177
fod.set_positions(fod_b)
write(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config/test/fodAng',fod,format='xyz')

s=np.loadtxt(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config/test/frmorb-137',skiprows=2,usecols=(1,2,3))
write=np.savetxt(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config/test/FRMORB-{i}',k,fmt='%1.10f',header='5 0',comments='',newline='\n')

####POSCAR to xyz using ase###
angtoBohr=[]
for i in range(10):
    a=read(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O/perturb_FOD_no_atoms/files/POSCAR-{i}') ###Atoms object
    print('a:',a)
    print(a.get_positions())
    write(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O/perturb_FOD_no_atoms/files/frmorb-{i}',a,format='xyz')
