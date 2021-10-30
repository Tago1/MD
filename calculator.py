#!/usr/bin/env python
# function: define a function to calculate evolved physical quantities
# editted by PenghuDu, 20211012

import math
from cell import Cell
from atom import Atom

## some physical constant and unit conversion 
Boltz = 1.38064852 * 10**(-23) # unit in kg m^2 s^-2 K^-1
Avogadro = 6.02214086 * 10**23    # unit in mol^-1
eV2Joule = 1.60217662 * 10**(-19)
kinetic_metal2iso = 10 / Avogadro


## define a funtion to abtain shortest distance from atom pair index, considering PBC
def shortest_dist(pos_i, pos_j, cell):
    delta_ij = [0] * 3
    inbox = True

    for k in range(3):  # x, y, z three directions
        delta_ij[k] = pos_i[k] - pos_j[k]
        while delta_ij[k] <= (-1 * cell[k][k] / 2): # 'while' circulation equals to many times of 'if';
            delta_ij[k] += cell[k][k]
            inbox = False
        while delta_ij[k] > (cell[k][k] / 2):
            delta_ij[k] -= cell[k][k]
            index = False
    
    dist_ij = math.sqrt(delta_ij[0]**2 + delta_ij[1]**2 + delta_ij[2]**2)

    return dist_ij, inbox

## define a function to restrict atom in main box
def restrict_in_box(pos, cell):
    for k in range(3):  # x, y, z three directions
        while pos[k] <= (-1 * cell[k][k]): 
            pos[k] += cell[k][k]
        while pos[k] > cell[k][k]:
            pos[k] -= cell[k][k]    
    
    return pos


## define a function to calculate atoimc forces and total potential energy of one step
def calculate_force_potential(param, cell, position, nlist, neighbor_list):
    natom = param['natom']
    input_config = param['geo_dir']
    rcut = param['rcut']
    neighbor_n = param['neighbor_n']
    epsilon = param['epsilon']
    sigma = param['sigma']

    tot_pot = 0
    force = [ [0] * 3 ] * natom

    for i in range(natom):
        pot_i = 0
        force_i = [0] * 3

        for j in neighbor_list[i]:

            dist_ij, inbox = shortest_dist(position[i], position[j], cell)
            if dist_ij > rcut:
                H = 0
            else:
                H = 1 

            force_i_j = [0] * 3
            for k in range(3): 
                force_i_j[k] = 4*epsilon * (12*(sigma/dist_ij)**12 - 6*(sigma/dist_ij)**6) * (float(position[i][k]) - float(position[j][k])) / (dist_ij ** 2) * H
                force_i[k] += force_i_j[k]

            if inbox == True: # only atoms in box and in rcut contribute to total petential energy which is defined as all energy of one box;
                pot_ij = 4*epsilon * ((sigma/dist_ij)**12 - (sigma/dist_ij)**6)
                pot_cut = 4*epsilon * ((sigma/rcut)**12 - (sigma/rcut)**6)
                pot_trunc_ij = (pot_ij - pot_cut) * H   # truncated potential by rcut, minus pot_cut for continuity in rcut. 
                pot_i += (0.5 * pot_trunc_ij)

        force[i] = force_i
        tot_pot += pot_i

    return tot_pot, force

## define a function to calculate atomic velocity and position, kinetic energy and temporary temperature of one step
def calculate_velocity_position(param, force, old_position, old_velocity):
    natom = param['natom']
    input_config = param['geo_dir'] 
    step = param['step']
    symb, _ = Atom().read_position(input_config)

    mass= [0] * natom
    accelerator = [ [0] * 3 ] * natom
    new_position = [ [0] * 3 ] * natom
    new_velocity = [ [0] * 3 ] * natom

    kinetic = 0
    for i in range(natom):
        mass[i] = Atom().dict[symb[i]]
        kinetic_i = 0
        for j in range(3):     
            accelerator[i][j] = force[i][j] / mass[i] * eV2Joule * Avogadro * 10**12 * 10**(-14)
            new_velocity[i][j] = old_velocity[i][j] + accelerator[i][j] * step
            new_position[i][j] = old_position[i][j] + (old_velocity[i][j] * step + 0.5 * accelerator[i][j] * step**2)
            kinetic_i_j = 0.5 * mass[i] * old_velocity[i][j]**2
            kinetic_i += kinetic_i_j
        kinetic += kinetic_i    # unit in g/mol A^2/ps^2 
    kinetic = kinetic * kinetic_metal2iso   # unit in J
    temp = kinetic / (0.5 * natom * Boltz) 
    kinetic = kinetic / eV2Joule    # unit in eV

    return new_position, new_velocity, kinetic, temp