import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#########################################################################################
# Plot the normal plot
file_name = 'result.csv'
dir = 'plot'

dat = pd.read_csv(file_name)
 
#Bad structure

# bad = dat.loc[dat.iloc[:,3] > 1.0, :]
lst_do = dat.iloc[:,0].unique()
# lst_bad = bad.iloc[:,0].unique().tolist()

for molecule in lst_do:
    # if molecule in lst_bad:
        # continue
    tmp = dat.loc[dat.iloc[:,0] == molecule,:]
    plt.plot(tmp.iloc[:,1], tmp.iloc[:,2])
    plt.xlim(-180,180)
    plt.title('The torsional profile of molecule {} - unscaled'.format(molecule))
    plt.xlabel('Torsion angle (degree)')
    plt.ylabel('Energy (kJ/mol)')
    plt.savefig(dir + '/' + str(molecule) + 'a.png',dpi=300)
    plt.clf()
    tmp.iloc[:,2] = tmp.iloc[:,2] / 4.184
    tmp.iloc[:,2] = tmp.iloc[:,2] - tmp.iloc[:,2].min()
    plt.plot(tmp.iloc[:,1], tmp.iloc[:,2])
    plt.xlim(-180,180)
    plt.title('The torsional profile of molecule {} - scaled'.format(molecule))
    plt.xlabel('Torsion angle (degree)')
    plt.ylabel('Energy (kcal/mol)')
    plt.savefig(dir + '/' + str(molecule) + 'b.png',dpi=300)
    plt.clf()
    


#########################################################################################
# Plot with ref

file_name = 'result_86_fix.csv'
ref_name = 'ccsd.csv'
dir = 'plot'

dat = pd.read_csv(file_name)
lst_do = dat.iloc[:,0].unique()
ref = pd.read_csv(ref_name, index_col = 0)

f = open('map.txt','r')
map = {}
lines = f.readlines()
f.close()
for line in lines:
    line = line.split()
    if len(line) == 2:
        map[int(line[0])] = int(line[1])
    

for molecule in lst_do:
    compute = dat.loc[dat.iloc[:,0] == molecule, :]
    compute.iloc[:,2] = compute.iloc[:,2] / 4.184
    compute.iloc[:,2] = compute.iloc[:,2] - compute.iloc[:,2].min()
    reference = ref.iloc[:,molecule]
    plt.plot(compute.iloc[:,1], compute.iloc[:,2], c = 'b', label = 'OpenMM')
    plt.plot(compute.iloc[:,1], reference, c = 'r', label = 'CCSD')
    plt.title('OpenMM and CCSD energy of molecule {} - {}'.format(molecule, map[molecule]))
    plt.xlabel('Torsion angle (degree)')
    plt.ylabel('Energy (kcal/mol)')
    plt.legend(loc = 'upper right')
    plt.savefig(dir + '/' + str(molecule) + '.png', dpi = 300)
    plt.clf()
    


