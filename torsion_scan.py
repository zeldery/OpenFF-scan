import os
import numpy as np
import pandas as pd


from simtk.openmm.app import PDBFile
from openforcefield.topology import Topology
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule
from simtk import openmm, unit
from simtk.openmm import CustomTorsionForce

from openmoltools.forcefield_generators import generateForceFieldFromMolecules, generateOEMolFromTopologyResidue

pdb_format = 'scan_{}.pdb'
sdf_format = 'scan_{}.sdf'

def read_data(file_name):
    f = open(file_name,'r')
    lines = f.readlines()
    f.close()
    list_name = []
    list_list_atom = []
    for line in lines:
        tmp = line.split()
        if len(tmp) != 5:
            continue
        list_name.append(tmp[0])
        list_atom = []
        for i in range(1,5):
            list_atom.append(int(tmp[i]))
        list_list_atom.append(list_atom)
    return list_name, list_list_atom

def write_result(file_name, lst_name, lst_torsion, results):
    f = open(file_name, 'w')
    f.write('Name,Angle,Energy,Restrain_Energy\n')
    for i, name in enumerate(lst_name):
        for j, torsion in enumerate(lst_torsion):
            f.write('{},{},{},{}\n'.format(name,torsion, result[i][j][0], result[i][j][1]))
    f.close()

#def get_pdb_list(dir):
#    lst_dir = os.listdir(dir)
#    raw_names = []
#    for name in lst_dir:
#        raw_names.append(name[:-4])
#    return raw_names

def extract_energy(simulation, float_conversion = True):
    forcegroups = {}
    n = simulation.system.getNumForces()
    for i in range(n):
        force = simulation.system.getForce(i)
        force.setForceGroup(i)
        forcegroups[force] = i
    energies = {}
    energies_float = []
    for f,i in forcegroups.items():
        data = simulation.context.getState(getEnergy = True, groups = 2**i).getPotentialEnergy()
        energies[f] = data
        energies_float.append(float(str(data).split()[0]))
    if float_conversion:
        return energies_float 
    else:
        return energies

def extract_torsion_angle(simulation, atom_list):
    # Need to change
    return 0.0

def float_quantity(item):
    return 10*float(str(item).split()[0])

def get_atom(topo):
    lst_atom = []
    for atom in topo.atoms():
        lst_atom.append(atom.element.symbol)
    return lst_atom

def write_xyz(coordinate, topology, file_name):
    f = open(file_name, 'w')
    lst_atom = get_atom(topology)
    f.write(str(len(lst_atom)) + '\n\n')
    for i, atom in enumerate(coordinate):
        f.write(lst_atom[i] + ' ' )
        for coor in atom:
            f.write(str(float_quantity(coor)) + ' ')
        f.write('\n')
    f.close()

def write_xml(topology, file_name):
    residues = [ residue for residue in topology.residues()]
    residue = residues[0]
    molOE = generateOEMolFromTopologyResidue(residue, geometry=False, tripos_atom_names = True)
    molOE.SetTitle('MOL')
    ffxml = generateForceFieldFromMolecules([molOE])
    f = open(file_name, 'w')
    f.write(ffxml)
    f.close()

def make_restrain_torsion(list_atom, angle, magnitude = 1e6):
    '''Create the external harmonic torsion force'''
    restrain_force = CustomTorsionForce('k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); pi = 3.1415926')
    restrain_force.addPerTorsionParameter('k')
    restrain_force.addPerTorsionParameter('theta0')
    restrain_force.addTorsion(list_atom[0]-1, list_atom[1]-1, list_atom[2]-1, list_atom[3]-1,
                              [magnitude*unit.kilojoule/unit.mole, angle/180*3.1415926])
    return restrain_force


# def scan(name, lst_angle, pdb_dir, sdf_dir):
    # time_step = 2*unit.femtoseconds  # simulation timestep
    # temperature = 300*unit.kelvin  # simulation temperature
    # friction = 1/unit.picosecond  # collision rate
    # minimize_tolerance = 1e-5 * unit.kilojoule/unit.mole
    # minimize_iteration_step = 1000000
    # forcefield = ForceField('openff-1.1.0.offxml')
    
    # pdbfile = PDBFile(pdb_dir + '/' + pdb_format.format(name))
    # uni_mol = Molecule.from_file(sdf_dir + '/' + sdf_format.format(name))
    
    # list_energy = []
    # for angle in lst_angle:
        # # Load the structure

        # topo = pdbfile.topology
        # topo_ff = Topology.from_openmm(topo, [uni_mol])
        # system = forcefield.create_openmm_system(topo_ff)
        # restrain_force = make_restrain_torsion(list_atoms[i], float(angle))
        # system.addForce(restrain_force)
        
        # integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
        # simulation = openmm.app.Simulation(topo, system, integrator)
        # positions = pdbfile.getPositions()
        # simulation.context.setPositions(positions)
        # simulation.context.setVelocitiesToTemperature(temperature)
        # simulation.minimizeEnergy(tolerance = minimize_tolerance,
                                  # maxIterations = minimize_iteration_step)
        # energy_list = extract_energy(simulation)
        # sum_energy = 0.0
        # for j in range(len(energy_list)-1):
            # sum_energy += energy_list[j]
        # list_energy.append([sum_energy, energy_list[-1]])
    # return list_energy

def minimize(dat_file, lst_angle, pdb_dir, sdf_dir, coor_dir = None, xml_dir = None):
    # The simulation configuration
    time_step = 2*unit.femtoseconds  # simulation timestep
    temperature = 300*unit.kelvin  # simulation temperature
    friction = 1/unit.picosecond  # collision rate
    minimize_tolerance = 1e-5 * unit.kilojoule/unit.mole
    minimize_iteration_step = 1000000
    forcefield = ForceField('openff-1.1.1.offxml')
    
    list_name, list_atoms = read_data(dat_file)
    list_energies = []
    for i, name in enumerate(list_name):
        pdbfile = PDBFile(pdb_dir + '/' + pdb_format.format(name))
        uni_mol = Molecule.from_file(sdf_dir + '/' + sdf_format.format(name))
        
        list_energy = []
        
        previous_structure = pdbfile.getPositions()
        for angle in lst_angle:
            # Load the structure

            topo = pdbfile.topology
            topo_ff = Topology.from_openmm(topo, [uni_mol])
            system = forcefield.create_openmm_system(topo_ff)
            restrain_force = make_restrain_torsion(list_atoms[i], float(angle), 1e6)
            system.addForce(restrain_force)
            
            integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
            simulation = openmm.app.Simulation(topo, system, integrator)
            #positions = pdbfile.getPositions()
            #simulation.context.setPositions(positions)
            simulation.context.setPositions(previous_structure)
            #simulation.context.setVelocitiesToTemperature(temperature)
            simulation.minimizeEnergy(tolerance = minimize_tolerance,
                                      maxIterations = minimize_iteration_step)
            energy_list = extract_energy(simulation)
            sum_energy = 0.0
            for j in range(len(energy_list)-1):
                sum_energy += energy_list[j]
            list_energy.append([sum_energy, energy_list[-1]])
            previous_structure = simulation.context.getState(getPositions = True).getPositions()
            if coor_dir != None:
                write_xyz(previous_structure, topo, 
                          coor_dir + '/' + str(name) + '_' + str(angle) + '.xyz')
            if xml_dir != None:
                write_xml(topo,
                          xml_dir + '/' + str(name) + '_' + str(angle) + '.xml')
        list_energies.append(list_energy)
        
        # Check
        f = open('check.log','a')
        f.write(str(i) + '\n')
        f.close()
        
    return list_energies
            
        
if __name__ == '__main__':
    lst_name, _ = read_data('indicies.txt')
    result = minimize('indicies.txt',range(-180, 181, 10), 'pdb', 'sdf', 'coor','xml')
    write_result('result.csv',lst_name, range(-180, 181, 10), result)
    















