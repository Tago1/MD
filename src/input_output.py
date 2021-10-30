#!/usr/bin/env python3
# funtion: define some functions to process input and output;
# eitted by PenghuDu, 20211016;


# funtion: define a function to read MD parameters from input file and print;
def read_paramter(input):

    # read each line of input file
    f_in = open(input, 'r+')
    lines = f_in.readlines()
    Num_lines = len(lines)
    f_in.close()

    # construct a dictionary to save MD input parameters
    param = {}
    
    for i in range(Num_lines):
        if 'id' in lines[i]:
            param['id'] = str(lines[i].strip().split()[1])
        if 'natom' in lines[i]:
            param['natom'] = int(lines[i].strip().split()[1])
        if 'geo_dir' in lines[i]:
            param['geo_dir'] = str(lines[i].strip().split()[1])
        if 'rcut' in lines[i]:
            param['rcut'] = float(lines[i].strip().split()[1])
        if 'neighbor_r' in lines[i]:
            param['neighbor_r'] = float(lines[i].strip().split()[1])
        if 'neighbor_n' in lines[i]:
            param['neighbor_n'] = int(lines[i].strip().split()[1])
        if 'epsilon' in lines[i]:
            param['epsilon'] = float(lines[i].strip().split()[1])
        if 'sigma' in lines[i]:
            param['sigma'] = float(lines[i].strip().split()[1])
        if 'mass' in lines[i]:
            param['mass'] = float(lines[i].strip().split()[1])
        if 'read_vel' in lines[i]:
            param['read_vel'] = int(lines[i].strip().split()[1])
        if 'temp' in lines[i]:
            param['temp'] = int(lines[i].strip().split()[1])
        if 'step' in lines[i]:
            param['step'] = int(lines[i].strip().split()[1])
        if 'nstep' in lines[i]:
            param['nstep'] = int(lines[i].strip().split()[1])

                   
    print(param)    # print all MD paramter on screen;
    
    return param

## define a funtion to print quantities of each step
def output(istep, symb, position, velocity, force, tot_pot, tot_kin, temp):
    f_pos = open('position.dat', 'a+')
    f_vel = open('velocity.dat', 'a+')
    f_force = open('force.dat', 'a+')
    f_log = open('run.log', 'w+')
    print('\n#MD step:', istep, file=f_pos)
    print('\n#MD step:', istep, file=f_vel)
    print('\n#MD step:', istep, file=f_force)
    print('\n#MD step:', istep, file=f_log)
    print('%ATOMIC_POSTION', file=f_pos)
    print('%ATOMIC_VELOCITY', file=f_vel)
    print('%ATOMIC_FORCE', file=f_force)
    for i in range(len(symb)):
        print(symb[i], *position[i], file=f_pos)
        print(symb[i], *velocity[i], file=f_vel)
        print(symb[i], *force[i], file=f_force)

    print((tot_pot+tot_kin), tot_pot, tot_kin, temp, file=f_log)
    f_pos.close()
    f_vel.close()
    f_force.close()
    f_log.close()

    return 1
    
    
    
    