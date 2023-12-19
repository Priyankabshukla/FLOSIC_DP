import numpy as np
import dpdata
import ase
from ase.io import read,write
import os
import matplotlib.pyplot as plt
from ase import Atoms
import shutil

#Dataset generation
#Use dpdata to create different configurations of water perturbing only the FODs and not H2O atoms. This is because for now, we want to include forces in the loss function.
#Pert_FOD = 200 (0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8) = 20 config each

#Orginal position of H2O atoms that should go in the CLUSTER file.
#Unit cell of 10 Å
#O 4.9999571664500007 5.0000000000000000 5.2207225380500004
#H 4.9999571664500007 6.4299727700000000 4.1138589180500009
#H 4.9999571664500007 3.5700272300000000 4.1138589180500009
#Original positions of FODs that will be perturbed
#He 4.9999470468500009 5.0026075396999996 5.2162526237500009
#He 5.0001625276500006 5.9905207688999997 4.4153197348500006
#He 5.0001692131500004 4.0072238908999998 4.4479858855500005
#He 6.1701970542500009 4.9935068753999996 5.8861382761500005
#He 3.8298029457500009 4.9935210320000003 5.886141081950000

def perturb(infile,outfile,pert_num,atom_per):
    perturbed_system = dpdata.System(f'{infile}').perturb(pert_num=pert_num,
    cell_pert_fraction=0.0,
    atom_pert_distance=atom_per,
    atom_pert_style='normal')
    perturbed_system.to('vasp/poscar',f'{outfile}')
    
pert_dist=[0.1,0.15,0.2,0.25,0.3,0.4,0.45,0.5,0.6,0.7]

for i in range(0,20):
    perturb('/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/PertFOD/POSCAR',
            f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/PertFOD/POSCAR-{i}',
           20,0.1)
           
for i in range(20):
    a=read(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/PertFOD/POSCAR-{i}') ###Atoms object
    write(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/PertFOD/frmorb-{i}',a,format='xyz')

## ASE xyz to FRMORB
## Create FRMORB file from perturb FOD
for i in range(0,20):
    s=np.loadtxt(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/PertFOD/frmorb-{i}',skiprows=2,usecols=(1,2,3))
    write=np.savetxt(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/FRMORB-files/FRMORB-{i}',s,fmt='%1.10f',header='5 0',comments='',newline='\n')

## Create a single CLUSTER file for pert FOD systems
target_PATH='/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/CLUSTER-file'

infile=np.loadtxt(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/CLUSTER-org',skiprows=2,usecols=(1,2,3))
#     print(len(infile))
if len(infile) == 8:
    X_Atoms=np.array(infile[-3:,0])
    Y_Atoms=np.array(infile[-3:,1])
    Z_Atoms=np.array(infile[-3:,2])
elif len(infile) == 3:
    X_Atoms=np.array(infile[-3:,0])
    Y_Atoms=np.array(infile[-3:,1])
    Z_Atoms=np.array(infile[-3:,2])

atom_no=np.array(['8','1','1'])
keyword_ALL=np.array(['ALL','ALL','ALL'])
AB=np.zeros(keyword_ALL.size, dtype=[('var1',float),('var2',float),('var3', float),('var4', 'U6'),('var5','U6')]) ##structured array
AB['var1'] = X_Atoms
AB['var2'] = Y_Atoms
AB['var3'] = Z_Atoms
AB['var4'] = atom_no
AB['var5'] = keyword_ALL


np.savetxt(f'{target_PATH}/CLUSTER',AB,fmt='%10.10f %10.10f %10.10f %10s %10s ',
                   header='LDA-PW91*LDA-PW91\nGRP\n3  Number of inequivalent atoms',
                  footer='0.00000000      0.00000000      0.00000000   ALL\n0.0     0.000     Charge and Moment\n ',comments='')
                  
#Visualize

for i in range(20):
    atom_fod=Atoms('He5',positions=np.loadtxt(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/FRMORB-files/FRMORB-{i}',skiprows=1,usecols=(0,1,2)))
    atom_cluster=Atoms('OHH',positions=np.loadtxt(open(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/CLUSTER-file/CLUSTER','rt').readlines()[:6:], skiprows=3,usecols=(0,1,2)))
    atom_fod += atom_cluster
    write(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/Visualize/{i}.xyz',atom_fod,format='xyz')
    
### Perform LDA DFT calculation
for i in range(0,20):
#     print(i)
    os.chdir('/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03')
    os.chdir('RUN_FLOSIC')
    os.mkdir(f'{i}')
    os.chdir(f'{i}')
    os.system(f'cp /bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/FRMORB-files/FRMORB-{i} .')
    os.system(f'cp /bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/CLUSTER-file/CLUSTER .')
    os.system('cp /bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/nrlmolDft.opa.batch .')
    os.system('cp /bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/NRLMOL_INPUT.DAT .')
    f = open("nrlmolDft.opa.batch", "rt")
    data = f.read()
    data = data.replace('runs',str(i))
    f.close()

    f = open("nrlmolDft.opa.batch","wt")
    f.write(data)
    f.close()
    
    os.system('sbatch nrlmolDft.opa.batch')
    
    
fail=[]
success=[]
path = '/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/RUN_FLOSIC/'
for d in os.listdir(path):
    print(d)
    if len(os.listdir(path+d)) > 4:
        slurm_file = [x for x in os.listdir(path+d) if x.split('.')[-1] == 'out']
        slurm_path = path+d+'/'+slurm_file[0]
        if os.stat(slurm_path).st_size == 0:
            print('Success: ', d)
            success.append(d)
        else:
            s = open(slurm_path,'r').readlines()[-2]
            print('Failed :::',d, s.split(' ')[-6])
            fail.append(d)
        
#### Get FLOSIC forces and energies from records file
path="/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/RUN_FLOSIC/"
energy=[]
force=[]

for p in range(200):
    force_counter=0
    no_atoms=0
    small_e=[]
    small_f=[]
    with open(f'{path}/{p}/records') as f:
        lines=f.readlines()
        for i,line in enumerate(lines):
            split_lines=line.split()
            if len(split_lines)==1:
                small_e.append(float(split_lines[0]))
                energy.append(small_e)  ## Correct this
            if len(split_lines)==2:
                n_up=int(split_lines[0]) ##Up FODs
                n_down=int(split_lines[1]) ##Down FODs
                no_atoms= n_up+n_down
            if i>=int(no_atoms)+3:
                small_f.append(split_lines)
        force.append(small_f)
        
    
print('energy',np.array(energy,dtype=float).T)
print('forces',np.array(force,dtype=float))
print('energy',np.array(energy).shape)
print('forces',np.array(force).shape)

plt.figure(figsize=(8,8))
plt.plot(np.sort(energy, axis=0),'ro')
plt.xlabel('Configuration',size=20)
plt.ylabel('FLOSIC SPE (Ha)', size=20)
plt.tick_params(direction = 'in', right = True, top = True)
plt.tick_params(axis='y',labelsize = 20,length=10, width=2)
plt.tick_params(axis='x',labelsize = 20,length=10, width=2)
#plt.savefig('/ihome/kjohnson/pbs13/PhD/FLOSIC/Python/DP_FLOSIC/H2O_perturbFOD.03/EnergyProfile.png',dpi=400,bbox_inches='tight')


## Create SOAP descriptors
import dscribe
from dscribe.descriptors import SOAP
from ase.visualize import view

# Setting up the SOAP descriptor
species=["He"]
soap = SOAP(
    species=species,
    periodic=False,
    rcut=10.0,
    sigma=0.25,
    nmax=6,
    lmax=6,
)
print(soap.get_number_of_features())

# Generate dataset of energies and forces
from ase import Atoms
n_samples=200
traj=[]
n_atoms=5
energies=np.array([])
forces=np.array([])
r=np.arange(n_samples)

for i in range(20):
    a=Atoms('He5',
            positions=Atoms('He5',positions=np.loadtxt(f'/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/RUN_FLOSIC/{i}/FRMORB-{i}'
                                                      ,skiprows=1,usecols=(0,1,2))).get_positions())
    traj.append(a)

## Get center of mass for optimized FOD coordinates

h2o_fod_pos=read('/bgfs/kjohnson/pbs13/FLOSIC/Projects/DP_FLOSIC/H2O_200config.03/fodPlusCluster.xyz',format='xyz')
h2o_fod_pos.get_center_of_mass()

derivatives, descriptors = soap.derivatives(
    traj,
    positions=[[h2o_fod_pos.get_center_of_mass()]]*len(r),
    method="analytical"
)
print(derivatives.shape)
print(descriptors.shape)

plt.figure(figsize=(8,8))
for i in range(len(r)):
    plt.plot(descriptors[i].T)
plt.xlabel('Features',size=20)
plt.ylabel('Descriptors D(r)', size=20)
plt.tick_params(direction = 'in', right = True, top = True)
plt.tick_params(axis='y',labelsize = 20,length=10, width=2)
plt.tick_params(axis='x',labelsize = 20,length=10, width=2)

# Save to disk for later training
np.save("/ihome/kjohnson/pbs13/PhD/FLOSIC/Python/DP_FLOSIC/H2O_perturbFOD.03/r.npy", r)
np.save("/ihome/kjohnson/pbs13/PhD/FLOSIC/Python/DP_FLOSIC/H2O_perturbFOD.03/E.npy", energies)
np.save("/ihome/kjohnson/pbs13/PhD/FLOSIC/Python/DP_FLOSIC/H2O_perturbFOD.03/D.npy", descriptors)
np.save("/ihome/kjohnson/pbs13/PhD/FLOSIC/Python/DP_FLOSIC/H2O_perturbFOD.03/dD_dr.npy", derivatives)
np.save("/ihome/kjohnson/pbs13/PhD/FLOSIC/Python/DP_FLOSIC/H2O_perturbFOD.03/F.npy", forces)


