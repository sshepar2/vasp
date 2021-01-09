import sys
import os
import io

'''This script is for:

    1) Removing selective dynamics from a POSCAR or CONTCAR file '-r'
    2) Changing which atoms are free to move and which are fixed '-F' or '-T'
    
    Ex1:
    >> python select.py -T POSCAR
    >> [1,23]
    The above adds 'T  T  T' next to the first 23 atoms in your POSCAR file.
    If the Selective Dynamics line is not already added, this will also add 'Selective Dynamics'
    in the appropriate place.
    
    Ex2:
    >> python select.py -F CONTCAR
    >> []
    Th above adds 'F  F  F' next to all atoms in the CONTCAR file. Will also add 'Selective Dynamics'
    in the appropriate place if not already there.
'''

file = sys.argv[2]

infile = open(sys.argv[2],'r')
data = infile.readlines()
infile.close()

atoms = data[6].split()
num_atoms = 0
for i in range(len(atoms)):
    num_atoms = num_atoms + int(atoms[i])


switch = data[7].split()
letter = switch[0][0]
#turning on/changing selective dynamics
# ex: python thisscript.py -T CONTCAR

if sys.argv[1] == '-t' or sys.argv[1] == '-T' or sys.argv[1] == '-f' or sys.argv[1] == '-F':
    which = input('To which atoms? (ex: [], for all, or [1,6] for atoms 1-6, or ([1,6],[9,15]), etc.: ')
    if sys.argv[1] == '-t' or sys.argv[1] == '-T':
        mode = 'T'
        antimode = 'F'
    elif sys.argv[1] == '-f' or sys.argv[1] == '-F':
        mode = 'F'
        antimode = 'T'
    if which == []:  # changing all
        if letter == 'c' or letter == 'C' or letter == 'D' or letter == 'd':
            bu = open('old','w+')
            bu.writelines(data)
            bu.close()
            os.remove(sys.argv[2])
            line = 'Selective Dynamics\n'
            f = open(sys.argv[2],'a+')
            f.write(data[0])
            f.write(data[1])
            f.write(data[2])
            f.write(data[3])
            f.write(data[4])
            f.write(data[5])
            f.write(data[6])
            f.write(line)
            f.write(data[7])
            for i in range(num_atoms):
                line = data[i+8].split()
                f.write("%20s %20s %20s %4s %4s %4s\n" %(line[0],line[1],line[2],mode,mode,mode))
            f.close()
        else:
            for i in range(num_atoms):
                line = data[i+9].split()
                data[i+9] = str("%20s %20s %20s %4s %4s %4s\n" % (line[0], line[1], line[2], mode, mode, mode))
            f = open(sys.argv[2],'w')
            f.writelines(data)
            f.close()

    if which != [] and type(which) == list:
        if len(which) == 1: #changing one
            atom = which[0]
            if letter == 'c' or letter == 'C' or letter == 'D' or letter == 'd':
                line = 'Selective Dynamics\n'
                f = open(sys.argv[2],'a')
                f.write(data[0])
                f.write(data[1])
                f.write(data[2])
                f.write(data[3])
                f.write(data[4])
                f.write(data[5])
                f.write(data[6])
                f.write(line)
                f.write(data[7])
                for i in range(atom-1):
                    line = data[i+8].split()
                    f.write("%20s %20s %20s %4s %4s %4s\n" %(line[0], line[1], line[2], antimode, antimode, antimode))
                line = data[atom-1+8].split()
                f.write("%20s %20s %20s %4s %4s %4s\n" %(line[0], line[1], line[2], mode, mode, mode))
                for i in range(atom,num_atoms):
                    line = data[i+8].split()
                    f.write("%20s %20s %20s %4s %4s %4s\n" %(line[0], line[1], line[2], antimode, antimode, antimode))
                f.close()
            else:
                line = data[atom-1+9].split()
                data[atom-1+9] = str("%20s %20s %20s %4s %4s %4s\n" % (line[0], line[1], line[2], mode, mode, mode))
                f = open(sys.argv[2],'w')
                f.writelines(data)
                f.close()

        if len(which) == 2: #changing one range
            start = which[0]
            end = which[1]
            if letter == 'c' or letter == 'C' or letter == 'D' or letter == 'd':
                line = 'Selective Dynamics\n'
                f = open(sys.argv[2],'a')
                f.write(data[0])
                f.write(data[1])
                f.write(data[2])
                f.write(data[3])
                f.write(data[4])
                f.write(data[5])
                f.write(data[6])
                f.write(line)
                f.write(data[7])
                for i in range(start-1):
                    line = data[i+8].split()
                    f.write("%20s %20s %20s %4s %4s %4s\n" %(line[0], line[1], line[2], antimode, antimode, antimode))
                for i in range(start-1,end):
                    line = data[i+8].split()
                    f.write("%20s %20s %20s %4s %4s %4s\n" %(line[0], line[1], line[2], mode, mode, mode))
                for i in range(end,num_atoms):
                    line = data[i+8].split()
                    f.write("%20s %20s %20s %4s %4s %4s\n" %(line[0], line[1], line[2], antimode, antimode, antimode))
                f.close()
            else:
                for i in range(start-1,end):
                    line = data[i+9].split()
                    data[i+9] = str("%20s %20s %20s %4s %4s %4s\n" % (line[0], line[1], line[2], mode, mode, mode))
                f= open(sys.argv[2],'w+')
                f.writelines(data)
                f.close()
               
