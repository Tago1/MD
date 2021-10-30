#!/usr/bin/env python3
# function: define a Cell class for MD
# editted by PenghuDu, 20211013

class Cell:
    'primary cell for MD simulation'
    
    ## defing a function to read cell parameters from input configuration file.
    def read_cell(self, config):
        # read all lines from input configuration
        f_in = open(config, 'r+')
        lines = f_in.readlines()
        f_in.close()
        Num_lines = len(lines)
    
        # create a 3*3 matrix to save cell vectors from '%CELL_PARAMETER' partition of input configuration and return it.
        for i in range(Num_lines):
            if '%CELL_PARAMETER' in lines[i]:
                cell = []
                for j in range(Num_lines):
                    temp = lines[i+j+1].strip()
                    if '%ATOMIC_POSTION' in temp:
                        break
                    else:
                        if len(temp) != 0:
                            cell.append([float(temp.split()[i]) for i in range(3)])
                print('The cell lattice in input configuration:\n', cell)

        return cell