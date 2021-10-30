#!/usr/bin/env python3
# funtion: define a function to creat neighbor table for input configuration and output 'neighbor_list.dat';
# eitted by PenghuDu, 20211012;

import math
from cell import Cell
from atom import Atom


## define a funtion to calculation interatomic distance
def interatomic_dist(pos_i, pos_j): # pos_i and pos_j are both a list of size 3
    dist = math.sqrt( ( pos_i[0] - pos_j[0] )**2
                    + ( pos_i[1] - pos_j[1] )**2 
                    + ( pos_i[2] - pos_j[2] )**2)
    return dist
    

## define a creat neighbor list and write to text file 'neighbor_list'
def create_neighbor(param, cell, position):

    natom = param['natom']
    neighbor_n = param['neighbor_n']
    input_config = param['geo_dir']
    rcut = param['rcut']
    neighbor_r = param['neighbor_r']

    symb, _ = Atom().read_position(input_config)    # _ stands for non-sense variable name

    nlist = [0] * natom
    list = [ [0] * neighbor_n for _ in range(natom)]

    f_out = open('neigbor_list2.dat', 'w')

    for i in range(natom):
        for j in range(i+1, natom):
            delta_ij = [0] * 3

            for k in range(3):
                delta_ij[k] = position[i][k] - position[j][k]
                if delta_ij[k] <= (-1 * cell[k][k] / 2):
                    delta_ij[k] += cell[k][k]
                if delta_ij[k] > (cell[k][k] / 2):
                    delta_ij[k] -= cell[k][k]

            dist_ij = math.sqrt(delta_ij[0]**2 + delta_ij[1]**2 + delta_ij[2]**2)
            if dist_ij <= (rcut + neighbor_r):
                nlist[i] += 1 
                nlist[j] += 1
                list[i][nlist[i]-1] = j
                list[j][nlist[j]-1] = i

#        print(i, symble[i], *position[i], nlist[i], list[i])
        print(i, symb[i], *position[i], nlist[i], list[i], file=f_out)

    f_out.close()

    return nlist, list