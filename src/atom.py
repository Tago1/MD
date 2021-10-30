#!/usr/bin/env python3
# function: define a Atom class for MD
# editted by PenghuDu, 20211013

import numpy as np
import math

Boltz = 1.38064852 * 10**(-23) # unit in kg m^2 s^-2 K^-1
Avogadro = 6.02214086 * 10**23    # unit in mol^-1
eV2Joule = 1.60217662 * 10**(-19)
kinetic_metal2iso = 10 / Avogadro

class Atom:
    'atomic information in primary cell for MD simulation'

    # define pubic attributes
    dict = {'He': 4.0026,
            }   # element and its corresponding mass, unit in g/mol
    num_atoms = 0
    symb_atom = []
    pos_atom = []
    vel_atom = []
    init_vel = []

    ## defing a function to read atomic symble and postion from input configuration file.
    def read_position(self, config):

        f_in = open(config, 'r')
        lines = f_in.readlines()
        f_in.close()
        Num_lines = len(lines)

        for i in range(Num_lines):
            if '%ATOMIC_POSTION' in lines[i]:
#                print('The atomic positions in input configuration:')

                Atom.num_atoms = 0   
                Atom.symb_atom = []             
                Atom.pos_atom = []
                for j in range(Num_lines):
                    temp = lines[i+j+1].strip()
                    if len(temp) == 0:
                        break
                    else:
                        Atom.num_atoms += 1
                        Atom.symb_atom.append(str(temp.split()[0]))
                        Atom.pos_atom.append([float(temp.split()[1:][i]) for i in range(3)])
#                        print(Atom.num_atoms-1, Atom.symb_atom[num_atoms-1], *Atom.pos_atom[num_atoms-1])
                            
                print('The number of atoms in input configuration:', Atom.num_atoms)

        return Atom.symb_atom, Atom.pos_atom

## define a function to read atomic symbol and velocity from input configuration file.
    def read_velocity(self, config):

        f_in = open(config, 'r')
        lines = f_in.readlines()
        f_in.close()
        Num_lines = len(lines)

        for i in range(Num_lines):
            if '%ATOMIC_VELOCITY' in lines[i]:
#                print('The initial velocity of atoms in input configuration:')

                Atom.num_atoms = 0
                Atom.symb_atom = []
                Atom.vel_atom = []
                for j in range(Num_lines):
                    temp = lines[i+j+1].strip()
                    if len(temp) == 0:
                        break
                    else:
                        Atom.num_atoms += 1
                        Atom.vel_atom.append([float(temp.split()[1:][i]) for i in range(3)])
#                        print(Atom.num_atoms-1, Atom.symb_atom[num_atoms-1], *Atom.vel_atom[num_atoms-1])

#                print('The number of atoms in input configuration:', Atom.num_atoms)

        return Atom.vel_atom

    ## defing a function to initiate velocity
    def initiate_velocity(self, param):
        input_config = param['geo_dir']
        style_initiate = param['read_vel']
        
        f_in = open(input_config, 'r')
        lines = f_in.readlines()
        f_in.close()
        Num_lines = len(lines)

        for i in range(Num_lines):
            if '%ATOMIC_VELOCITY' in lines[i]:
#                print('The initial velocity of atoms in input configuration:')

                Atom.num_atoms = 0
                Atom.symb_atom = []
                Atom.vel_atom = []
                for j in range(Num_lines):
                    temp = lines[i+j+1].strip()
                    if len(temp) == 0:
                        break
                    else:
                        Atom.num_atoms += 1
                        Atom.vel_atom.append([float(temp.split()[1:][i]) for i in range(3)])
#                        print(Atom.num_atoms-1, Atom.symb_atom[num_atoms-1], *Atom.vel_atom[num_atoms-1])

#                print('The number of atoms in input configuration:', Atom.num_atoms)

        Atom.init_vel = [ [0] * 3 ] * Atom.num_atoms
        vel_random = 0.5 * np.random.rand(Atom.num_atoms, 3) # A uniform distribution random of (-0.5, 0.5) to initialize the velocities of each atom in the x, y, and z directions

        if style_initiate == 1:
            Atom.init_vel = Atom.vel_atom
            
        if style_initiate == 0:
            temp = param['temp']
            kinetic_energy = 0

            mass_sum = 0
            for i in range(Atom.num_atoms): 
                kinetic_energy_i = 0
                for k in range(3):
                    momentum_k = 0
                    for j in range(Atom.num_atoms):
                        mass_j = Atom.dict[Atom.symb_atom[j]]
                        mass_sum += mass_j

                        momentum_j_k = mass_j * vel_random[j][k]
                        momentum_k += momentum_j_k

                    vel_center_k = momentum_k / mass_sum
                    vel_random[i][k] = vel_random[i][k] - vel_center_k

                    kinetic_energy_i_k  = 0.5 * mass_j * (vel_random[j][k]) ** 2
                    kinetic_energy_i += kinetic_energy_i_k

                kinetic_energy += kinetic_energy_i

            factor = (3 * Atom.num_atoms / 2) * Boltz * temp / (kinetic_energy * kinetic_metal2iso)
            Atom.init_vel = math.sqrt(factor) * vel_random

        return Atom.init_vel
