#!/usr/bin/env python3
# funtion: creat an complete neighbor list for input configuration
# eitted by PenghuDu, 20211013

import math
import time

import numpy

import calculator
import input_output
import neighbor
from atom import Atom
from cell import Cell

start = time.clock()    # start time

param = input_output.read_paramter('md.in')    # read MD parameters from 'md.in'
natom = param['natom']
input_config = param['geo_dir']

cell = Cell().read_cell(input_config) 
symb, position = Atom().read_position(input_config)
nlist, list = neighbor.create_neighbor(param, cell, position)   # calculate atomic neighbor list for input configuration;

end = time.clock()  # stop time
print('Running time: %s seconds'%(end-start))   # print running time

