# -*- coding: utf-8 -*-

'''
   Only works well with Python 2.7 at the moment. This script takes a POSCAR or CONTCAR
   file (VASP input and output position files) and manipulates the atoms.
   
   Options:
       0-translate a set of atoms
       1-rotate a set of atoms
       2-stretch the unit cell
       3-shear the unit cell

   The updated coordinates are output in a file called POSCAR or CONTCAR
   depending on which you input.  The original POSCAR or CONTCAR file is
   saved in a file called POSCAR_old or CONTCAR_old.  Another file called
   z.py is output as well.  This is for use in another script and can be
   ignored or deleted.  The file contains the z component of the 3rd
   lattice vector.

   Acronyms:
       DC - Direct Coordinates
       CC - Cartesian Coordinates
       SD - Selective Dynamics
       cc - counter-clockwise

'''

import io
import os
import sys
import ast
from math import sqrt
from math import acos
from math import sin
from math import cos
from math import radians
from math import degrees

#####################################################################################################################
# DEFINING FUNCTIONS BEGIN ##########################################################################################
#####################################################################################################################

def mag(vec):
    
    '''
    Calculates the norm squared of a vector in CC
    '''

    return sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

def dot(vec1,vec2):
    
    '''
    Calculates dot product of two vectors in CC
    '''
    
    vec = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]
    return vec


def cross(vec1,vec2):

    '''
    Takes the cross product of two vectors in CC
    '''
    vec = [vec1[1]*vec2[2]-vec2[1]*vec1[2],vec2[0]*vec1[2]-vec1[0]*vec2[2],vec1[0]*vec2[1]-vec2[0]*vec1[1]]
    return vec


def build_c_d(basis):
    
    '''
    The purpose of this function is to build the matrix needed to convert from
    CC to DC using the current basis. The matrix should be saved to the 
    variable 'cdmat' i.e. use this function as done below.
    
    cdmat = build_c_d(basis)
    
    cdmat is used in the c_d function which performs the conversion from CC to DC.
    Input basis form:        [ [ax,ay,az], [bx,by,bz], [cx,cy,cz] ]
    Output matrix form: inv( [ [ax,bx,cx], [ay,by,cy], [az,bz,cz] ] )
    '''

    det = basis[0][0]*(basis[1][1]*basis[2][2]-basis[2][1]*basis[1][2]) \
            - basis[1][0]*(basis[0][1]*basis[2][2]-basis[2][1]*basis[0][2])\
            + basis[2][0]*(basis[0][1]*basis[1][2]-basis[1][1]*basis[0][2])
    cdmat = [[(basis[1][1]*basis[2][2]-basis[2][1]*basis[1][2])/det,\
            (basis[2][0]*basis[1][2]-basis[1][0]*basis[2][2])/det,\
            (basis[1][0]*basis[2][1]-basis[2][0]*basis[1][1])/det],\
            [(basis[2][1]*basis[0][2]-basis[0][1]*basis[2][2])/det,\
            (basis[0][0]*basis[2][2]-basis[2][0]*basis[0][2])/det,\
            (basis[2][0]*basis[0][1]-basis[0][0]*basis[2][1])/det],\
            [(basis[0][1]*basis[1][2]-basis[1][1]*basis[0][2])/det,\
            (basis[1][0]*basis[0][2]-basis[0][0]*basis[1][2])/det,\
            (basis[0][0]*basis[1][1]-basis[1][0]*basis[0][1])/det]]
    return cdmat


def c_d(points):
    
    '''
    Converts a list of points [ [p1x,p1y,p1z], [p2x,p2y,p2z],... ] from CC to DC.
    Be sure the variable 'cdmat' has been set for the current basis using the 
    build_c_d function.
    '''

    direct_points = []
    for i in range(len(points)):
        vec = points[i]
        new_vec_a = cdmat[0][0]*vec[0] + cdmat[0][1]*vec[1] + cdmat[0][2]*vec[2]
        new_vec_b = cdmat[1][0]*vec[0] + cdmat[1][1]*vec[1] + cdmat[1][2]*vec[2]
        new_vec_c = cdmat[2][0]*vec[0] + cdmat[2][1]*vec[1] + cdmat[2][2]*vec[2]
        new_vec = [new_vec_a,new_vec_b,new_vec_c]
        direct_points.append(new_vec)
    points = direct_points
    return points


def d_c(points):
    
    '''
    Does the reverse of the c_d function. Make sure the current 'basis' variable
    has been built. basis = [ [a],[b],[c] ]
    '''
    
    cart_points = []
    for i in range(len(points)):
        vec = points[i]
        new_vec_a = basis[0][0]*vec[0] + basis[1][0]*vec[1] + basis[2][0]*vec[2]
        new_vec_b = basis[0][1]*vec[0] + basis[1][1]*vec[1] + basis[2][1]*vec[2] 
        new_vec_c = basis[0][2]*vec[0] + basis[1][2]*vec[1] + basis[2][2]*vec[2] 
        new_vec = [new_vec_a,new_vec_b,new_vec_c]
        cart_points.append(new_vec)
    points = cart_points
    return points


def rot(vec,t,p):
    
    '''
    Rotates the vector 'vec' by an angle 't' IN DEGREES cc about the vector 'p'
    All vectors should be in CC.  cc - counterclockwise.
    '''
    
    t = radians(t)
    p_mag = sqrt(p[0]**2+p[1]**2+p[2]**2)
    p_unit = [p[0]/p_mag,p[1]/p_mag,p[2]/p_mag]
    px = p_unit[0]
    py = p_unit[1]
    pz = p_unit[2]
    d = sqrt(py**2 + pz**2)
    if py == 0 and pz == 0:
        R = [[1,0,0],[0,cos(t),-sin(t)],[0,sin(t),cos(t)]]
    else:
        R = [[d**2*cos(t)+px**2,px*py*(1-cos(t))-pz*sin(t),px*pz*(1-cos(t))+py*sin(t)],\
                [px*py*(1-cos(t))+pz*sin(t),(px**2*py**2+pz**2)*cos(t)/(d**2)+py**2,\
                (px**2-1)*py*pz*cos(t)/(d**2)-(py**2+pz**2)*px*sin(t)/(d**2)+py*pz],\
                [px*pz*(1-cos(t))-py*sin(t),(px**2-1)*py*pz*cos(t)/(d**2)+(py**2+pz**2)*px*sin(t)/(d**2)+py*pz,\
                (px**2*pz**2+py**2)*cos(t)/(d**2)+pz**2]]
    rvec = []
    for j in range(3):
        temp = vec[0]*R[j][0]+vec[1]*R[j][1]+vec[2]*R[j][2]
        rvec.append(temp)
    return rvec

#####################################################################################################################
# DEFINING FUNCTIONS END ############################################################################################
#####################################################################################################################


# initial prompt
dowhat = int(input('translate[0], rotate[1], stretch[2], or shear[3]?: '))
if len(sys.argv) == 1:
    towhat = int(input('CONTCAR[0] or POSCAR[1]?: '))
else:
    towhat = 2

#####################################################################################################################
# COLLECTING PRELIMINARY DATA  ######################################################################################
#####################################################################################################################

# opening file
if towhat == 0: # (file name CONTCAR)
    f = open('CONTCAR','r')
elif towhat == 1: # (file name POSCAR)
    f = open('POSCAR','r')
elif towhat == 2: # (file name specified by arg[1], manual input)
    f = open(sys.argv[1],'r')

# reading all data from file    
data = f.readlines()

# closing file
f.close()

# saving important data

# list of # of each atom type
atoms = data[6].split()
# calculating total atoms in system
num_atoms = 0
for i in range(len(atoms)):
    num_atoms = num_atoms + int(atoms[i])

# determining type of coordinates used, Cartesian (CC) of Direct (DC)
# and checking whether Selective Dynamics (SD) is activated
switch = data[7].split()
letter = switch[0][0]
if letter == 'D' or letter == 'd' or letter == 'C' or letter == 'c':
    points = []
    for i in range(num_atoms):
        row = []
        point = data[i+8].split()
        for j in range(3):
            row.append(float(point[j]))
        points.append(row)
elif letter == 'S' or letter == 's':
    points = []
    for i in range(num_atoms):
        row = []
        point = data[i+9].split()
        for j in range(3):
            row.append(float(point[j]))
        points.append(row)
    switch = data[8].split()
    letter = switch[0][0]

# saving basis vector data
avec = data[2].split(); bvec = data[3].split(); cvec = data[4].split()

# converting strings to floats and re-storing basis vectors
a = []; b = []; c = []

for i in range(3):
    a.append(float(avec[i]))
    b.append(float(bvec[i]))
    c.append(float(cvec[i]))

# building matrix of basis vectors [ [ax,ay,az],[bx,by,bz],[cx,cy,cz] ]
basis = []
basis.append(a); basis.append(b) ;basis.append(c)

# calculating magnitude of each basis vector
mag_a = mag(a); mag_b = mag(b); mag_c = mag(c)

# calculating alpha, beta, gamma angles between basis vectors
alpha = acos((b[0]*c[0]+b[1]*c[1]+b[2]*c[2])/(mag_b*mag_c))
beta = acos((a[0]*c[0]+a[1]*c[1]+a[2]*c[2])/(mag_a*mag_c))
gamma = acos((b[0]*a[0]+b[1]*a[1]+b[2]*a[2])/(mag_b*mag_a))

# building the matrix needed to convert from CC to DC
cdmat = build_c_d(basis)
    

# this file is used for a script doing a sequence of calls on transform.py (it exports the magnitude of the original 
# c basis vector to a file called 'z.py')
f = open('z.py','w+')
f.write(str(mag_c))
f.close


#####################################################################################################################
#####################################################################################################################
# BEGIN CODING FOR TRANSFORMATIONS ##################################################################################
#####################################################################################################################
#####################################################################################################################

#####################################################################################################################
# [2] - STRETCH BEGIN ###############################################################################################
#####################################################################################################################

'''
This comes in handy for slab systems, or systems with a vacuum region which you wish to change. I personally use this
to change my periodic unit cell with vacuum in the VASP program, to a unit cell without vacuum for a program which
will use fixed boundary conditions. Studying transport properties, for example.

The 'stretch' code works by first prompting for new basis vector magnitudes as [[new_mag_a],[new_mag_b],[new_mag_c]]. 
Notice the strange unecessary use of a list of lists of length 1. This is done so the option of [[],[],[]] will leave 
all basis vectors at the current length, likewise for a and b if [[],[],[30]] is fed to the prompt.

Since atoms in DC scale the same way (multply by the same factor) as the basis vectors, DC are preferred. Also, there 
is an option (which would need to be toggled off, align = 'T' -> 'F') which aligns the basis vectors as best as 
possible with the x,y,z coorinate system. 'align' is one of the first varables defined in this section. The aligning
is done in the following two steps:

1. rotating all basis vectors by the angle between basis vecter a and the x axis.
2. rotating all basis vectors by the angle between basis vector b and the xy-plane.

The coordates will now be close to: [ [ax,0,0], [bx,by,0], [cx,cy,cz] ]

Atoms in DC are also invariant to rotations.

Therefore, the stretch transformation is written for DC, and the first step is to convert them to DC if needed.

Each basis vectors is scaled by a factor given by the new_mag/old_mag. Then each individual atom coordinate in DC 
will be scaled by the inverse of that factor, so as to not strectch apart with the system. Before that it is checked 
that all atomic coordinates in DC are between 0-1, so they are constrained to modulo 1.

Since rotating the basis vectors is optional, here we proceed as if the transformation has been completed. If the
atomic coordinated were originally given in CC, then they are converted back, but not before rebuilding the basis
variable. 

If the option to rotate the basis to align nicely with x,y,z is activated (align = 'T'), then the coordinate 
transform matrix (cdmat) must be recalculated since the coordinates will then be transormed back to direct for 
rotations otherwise the DC will be incorrect since they will be converted using the wrong basis.

The rotations are performed as described, then again the basis variable and cdmat variable are rebuilt in case
needed. Atom coordinates are returned to original format and the new file is written with the new basis vectors
and new atom coordinates in the stretched cell system.

Avoid stretching such that atoms, once inside the unit cell, would end up outside the unit cell, say by shrinking
the basis vector by a large amount. The DC of the atom is scaled up and will be larger than one, since the modulo
was performed prior to the scaling. Now the atom will wrap around the cell an appear in the wrong place.

'''

# prompting for input
if dowhat == 2:
    new_mags = ast.literal_eval(input('Enter new cell dimensions to stretch to (ex: [[],[],[20]], means stretch c to 20Å): '))
    
    # align option, set to 'T' to align basis along x,y,z after stretching, and anything else, e.g. 'F' to leave them alone
    align = 'F'

    # converting atom coordinates to DC if provided in CC
    if letter == 'C' or letter == 'c':
        points = c_d(points)

    # scaling each basis vector based on new magnitudes from raw input
    scalers = [1,1,1]

    # scaling basis vector a
    if new_mags[0] != []:
        scalers[0] = new_mags[0][0]/mag_a
        new_a = [float('%f' % (scalers[0]*a[0])), float('%f' % (scalers[0]*a[1])), float('%f' % (scalers[0]*a[2]))]
        a = new_a
    
    # scaling basis vector b
    if new_mags[1] != []:
        scalers[1] = new_mags[1][0]/mag_b
        new_b = [float('%f' % (scalers[1]*b[0])), float('%f' % (scalers[1]*b[1])), float('%f' % (scalers[1]*b[2]))]
        b = new_b

    # scaling basis vector c
    if new_mags[2] != []:
        scalers[2] = new_mags[2][0]/mag_c
        new_c = [float('%f' % (scalers[2]*c[0])), float('%f' % (scalers[2]*c[1])), float('%f' % (scalers[2]*c[2]))]
        c = new_c

    # for this scaling method to work properly all coordinates need to be between 0 and 1
    for i in range(num_atoms):
        for j in range(3):
            points[i][j] = points[i][j] % 1

    # scaling coordinates due to length change of basis vectors
    for i in range(num_atoms):
        for j in range(3):
            points[i][j] = points[i][j]/scalers[j]

    # building new basis prior to the DC to CC conversion
    basis = []
    basis.append(a)
    basis.append(b)
    basis.append(c)

    # not technically needed for DC to CC, but if another CC to DC is performed it will be.
    cdmat = build_c_d(basis)

    # converting atom coordinates back to CC if provided in CC
    if letter == 'C' or letter == 'c':
        points = d_c(points)
    
    ###############################################################################################
    # BEGIN ROTATIONS (OPTIONAL)
    #
    # The following is optional and align axes so that 
    # a = (ax,0,0), b = (bx,by,0), c = (cx,cy,cz)
    ###############################################################################################
    
    # converting atom coordinates to DC if in provided in CC
    if align == 'T':
        if letter == 'C' or letter == 'c':
            points = c_d(points)
        
    
        # rotation 1) 

        # defining angle and vector to rotate basis vector a about to get along x axis
        v_ax = [0,a[2],-1.0*a[1]]
        if mag(v_ax) != 0:
            t_ax = degrees(acos(a[0]/mag(a)))
        
            # rotating basis vectors by rotation 1)
            a_ax = rot(a,t_ax,v_ax); a = a_ax
            b_ax = rot(b,t_ax,v_ax); b = b_ax
            c_ax = rot(c,t_ax,v_ax); c = c_ax
         
    
        # rotation 2)

        # defining angle and vector to rotate basis vector b about to get b xy-plane, while keeping a along x. (rotate about x)
        v_bxy = [1,0,0]
        t_bxy = degrees(acos(sqrt(b[1]**2+b[2]**2)/mag(b)))
        
        # rotating basis vectors by rotation 2)
        a_bxy = rot(a,t_bxy,v_bxy); a = a_bxy
        b_bxy = rot(b,t_bxy,v_bxy); b = b_bxy
        c_bxy = rot(c,t_bxy,v_bxy); c = c_bxy
        

        # rebuilding basis for any conversion
        basis = []
        basis.append(a)
        basis.append(b)
        basis.append(c)

        # rebuilding cdmat incase of CC to DC conversion occurs again
        cdmat = build_c_d(basis)

        
        # converting atom coordinates back to CC if provided in CC
        if letter == 'C' or letter == 'c':
            points = d_c(points)

    #################################################################################################
    # END ROTATIONS (OPTIONAL)
    #################################################################################################

    #################################################################################################
    # BEGIN WRITING OUTPUT
    #################################################################################################

    # creating backup files and writing original data to them
    if towhat == 0:
        f = open('CONTCAR_old','w+')
        f.writelines(data)
        f.close()
    elif towhat == 1:
        f = open('POSCAR_old','w+')
        f.writelines(data)
        f.close()
    
    # recalling if SD was activated
    select = data[7].split()
    sel = select[0][0]

    # replacing old basis vector data with new basis vectors
    data[2] = str('%22.16f %22.16f %22.16f\n' % (a[0], a[1], a[2]))
    data[3] = str('%22.16f %22.16f %22.16f\n' % (b[0], b[1], b[2]))
    data[4] = str('%22.16f %22.16f %22.16f\n' % (c[0], c[1], c[2]))

    # replacing old atom coordinates with new coordinates, including SD flags if necessary
    if sel == 'S' or sel == 's':
        for i in range(num_atoms):
            line = data[i+9].split()
            data[i+9] = str('%20.16f %20.16f %20.16f %4s %4s %4s\n' % (points[i][0], points[i][1], points[i][2], line[3], line[4], line[5]))
    else:
        for i in range(num_atoms):
            data[i+8] = str('%20.16f %20.16f %20.16f\n' % (points[i][0], points[i][1], points[i][2]))
    
    # writing new data to file with name of original input file
    if towhat == 0 or towhat == 2:
        f = open('CONTCAR','w+')
        f.writelines(data)
        f.close()
        #os.system("open -a 'vesta' CONTCAR") 
    elif towhat == 1:
        f = open('POSCAR','w+')
        f.writelines(data)
        f.close()
        #os.system("open -a 'vesta' POSCAR")

#####################################################################################################################
# [2] - STRETCH END #################################################################################################
#####################################################################################################################


#####################################################################################################################
# [0] - TRANSLATE BEGIN #############################################################################################
#####################################################################################################################

'''

    Translates a set of atoms in some direction and by some magnitude specified by an input vector

'''

if dowhat == 0:
    transvec = ast.literal_eval(input('Enter translation vector (ex: ["D",0.5,0.2,0.3]) where "D" means the translation vector you\'ve provided is in direct units: '))
    atom_set = ast.literal_eval(input('Which atoms? (ex: [1,5] for atoms 1 to 5, or [] for all atoms): '))


    if towhat == 0:
        f = open('CONTCAR_old','w+')
        f.writelines(data)
        f.close()
    elif towhat == 1:
        f = open('POSCAR_old','w+')
        f.writelines(data)
        f.close()


    if atom_set == []:
        atom_num = num_atoms
        atom_start = 1
    else:
        atom_num = atom_set[1] - atom_set[0] + 1
        atom_start = atom_set[0]

    sub_points = points[atom_start-1:atom_start-1+atom_num]

    #if input file in cartesian, this converts to direct coordinates
    if letter == 'C' or letter == 'c':
        sub_points = c_d(sub_points)
        
    if transvec[0] == "C" or transvec[0] == "c":
        transvec[1] = transvec[1]/mag_a
        transvec[2] = transvec[2]/mag_b
        transvec[3] = transvec[3]/mag_c

    for i in range(atom_num):
        for j in range(3):
            sub_points[i][j] = sub_points[i][j] + transvec[j+1]
            
    #converting back to cartesian if input file was cartesian
    if letter == 'C' or letter == 'c':
        sub_points = d_c(sub_points)
        
    select = data[7].split()
    sel = select[0][0]
            
    if sel == 'S' or sel == 's':
        for i in range(atom_num):
            line = data[i+9+atom_start-1].split()
            data[i+9+atom_start-1] = str('%20.16f %20.16f %20.16f %4s %4s %4s\n' % (sub_points[i][0], sub_points[i][1], sub_points[i][2], line[3], line[4], line[5]))
    else:
        for i in range(atom_num):
            data[i+8+atom_start-1] = str('%20.16f %20.16f %20.16f\n' % (sub_points[i][0], sub_points[i][1], sub_points[i][2]))
        
    if towhat == 0 or towhat == 2:
        f = open('CONTCAR','w+')
        f.writelines(data)
        f.close()
        #os.system("open -a 'vesta' CONTCAR")
    elif towhat == 1:
        f = open('POSCAR','w+')
        f.writelines(data)
        f.close()
        #os.system("open -a 'vesta' POSCAR")


#####################################################################################################################
# [1] - ROTATE BEGIN ################################################################################################
#####################################################################################################################

'''
    Rotate a set of atoms by some angle about a vector defined by inputing two points. Providing two points is useful
    since generally the rotation axis tends to come from the position of two atoms in the system, whose coordinates
    are available in the CONTCAR/POSCAR file. It is also difficult to intuitively know the vector components when
    the basis vectors are not orthogonal.

    The rotation is performed on CC always. After converting if necessary the rotation is performed in 3 major steps.
    1 - Perform a translation on all atoms, given by the distance between the origin and the tail of the rotation
        vector.
    2 - Perform the rotation using a rotation matrix build using the function 'rot' defined at the start of the file,
        calculated from the input rotation angle.
    3 - Translate all atoms back (reversing step 1)
 
'''



if dowhat == 1:
    rotvec = ast.literal_eval(input('Enter two points to define vector and one angle (in degrees) to rotate by (ex: [45,"D",[0,0,0],[1,0,0]] would rotate atoms 45 degrees counterclockwise about the x axis): '))
    atom_set = ast.literal_eval(input('Which atoms? (ex: [1,5] for atoms 1 to 5, or [] for all atoms): '))

    if towhat == 0:
        f = open('CONTCAR_old','w+')
        f.writelines(data)
        f.close()
    elif towhat == 1:
        f = open('POSCAR_old','w+')
        f.writelines(data)
        f.close()

    if atom_set == []:
        atom_num = num_atoms
        atom_start = 1
    else:
        atom_num = atom_set[1] - atom_set[0] + 1
        atom_start = atom_set[0]

    sub_points = points[atom_start-1:atom_start-1+atom_num]


    #converting to cartesian coordinates if given in Direct
    if letter == 'D' or letter == 'd':
        sub_points = d_c(sub_points)
    
    p1 = rotvec[2]
    p2 = rotvec[3]
    t = rotvec[0] # keep in degrees, the rot function converts it to degrees


    #converting input points to cartesian
    if rotvec[1] == 'D' or rotvec[1] == 'd':
        p1 = [rotvec[2][0]*a[0]+rotvec[2][1]*b[0]+rotvec[2][2]*c[0],\
                rotvec[2][0]*a[1]+rotvec[2][1]*b[1]+rotvec[2][2]*c[1],\
                rotvec[2][0]*a[2]+rotvec[2][1]*b[2]+rotvec[2][2]*c[2]]
        
        p2 = [rotvec[3][0]*a[0]+rotvec[3][1]*b[0]+rotvec[3][2]*c[0],\
                rotvec[3][0]*a[1]+rotvec[3][1]*b[1]+rotvec[3][2]*c[1],\
                rotvec[3][0]*a[2]+rotvec[3][1]*b[2]+rotvec[3][2]*c[2]]

    p = [(p2[0]-p1[0]),(p2[1]-p1[1]),(p2[2]-p1[2])]
    
    #p_mag = sqrt(p[0]**2+p[1]**2+p[2]**2)
    #p_unit = [p[0]/p_mag,p[1]/p_mag,p[2]/p_mag]

    #px = p_unit[0]
    #py = p_unit[1]
    #pz = p_unit[2]
    #d = sqrt(py**2 + pz**2)
        
    #Rotation Matrix
    #if py == 0 and pz == 0:
        #R = [[1,0,0],[0,cos(t),-sin(t)],[0,sin(t),cos(t)]]
    #else:
        #R = [[d**2*cos(t)+px**2,px*py*(1-cos(t))-pz*sin(t),px*pz*(1-cos(t))+py*sin(t)],
                #[px*py*(1-cos(t))+pz*sin(t),(px**2*py**2+pz**2)*cos(t)/(d**2)+py**2,(px**2-1)*py*pz*cos(t)/(d**2)-(py**2+pz**2)*px*sin(t)/(d**2)+py*pz],
                #[px*pz*(1-cos(t))-py*sin(t),(px**2-1)*py*pz*cos(t)/(d**2)+(py**2+pz**2)*px*sin(t)/(d**2)+py*pz,(px**2*pz**2+py**2)*cos(t)/(d**2)+pz**2]]
        
    #translate to origin T
    for i in range(atom_num):
        for j in range(3):
            sub_points[i][j] = sub_points[i][j] - p1[j]

    #rotating R
    rot_points = []
    for i in range(atom_num):
        rot_point = rot(sub_points[i],t,p)
        rot_points.append(rot_point)
    sub_points = rot_points

    #rot_points = []
    #for i in range(atom_num):
        #rot_point = []
        #for j in range(3):
            #temp = sub_points[i][0]*R[j][0]+sub_points[i][1]*R[j][1]+sub_points[i][2]*R[j][2]
            #rot_point.append(temp)
        #rot_points.append(rot_point)
    #sub_points = rot_points

                
    #translate back T^-1
    for i in range(atom_num):
        for j in range(3):
            sub_points[i][j] = sub_points[i][j] + p1[j]
        

    #converting from Cartesian back to direct if necessary if originally in direct
    if letter == 'D' or letter == 'd':
        sub_points = c_d(sub_points)

    select = data[7].split()
    sel = select[0][0]

    #################################################################################################
    # BEGIN WRITING OUTPUT
    #################################################################################################

    if sel == 'S' or sel == 's':
        for i in range(atom_num):
            line = data[i+9+atom_start-1].split()
            data[i+9+atom_start-1] = str('%20.16f %20.16f %20.16f %4s %4s %4s\n' % (sub_points[i][0], sub_points[i][1], sub_points[i][2], line[3], line[4], line[5]))
    else:
        for i in range(atom_num):
            data[i+8+atom_start-1] = str('%20.16f %20.16f %20.16f\n' % (sub_points[i][0], sub_points[i][1], sub_points[i][2]))
        
    if towhat == 0 or towhat == 2:
        f = open('CONTCAR','w+')
        f.writelines(data)
        f.close()
        #os.system("open -a 'vesta' CONTCAR")
    elif towhat == 1:
        f = open('POSCAR','w+')
        f.writelines(data)
        f.close()
        #os.system("open -a 'vesta' POSCAR")

#####################################################################################################################
# [1] - ROTATE END ##################################################################################################
#####################################################################################################################

#####################################################################################################################
# [3] - SHEAR BEGIN #################################################################################################
#####################################################################################################################

'''
    This is still being tested, but can be useful for shearing bulk systems by transforming the atomic coordinates as
    an angle between basis vectors is changed.

'''

if dowhat == 3:
    new_angles = ast.literal_eval(input('Enter new cell angles in degrees to shear to (ex: [[],[],[90]], means shear gamma to 90 degrees): '))
    
    # we want to preserve atom positions with respect to each other
    if letter == 'D' or letter == 'd':
        points = d_c(points)
    
    # setting new alpha, beta, and gamma angles
    new_g = gamma; new_b = beta; new_a = alpha

    if new_angles[0] != []:
        new_a = radians(new_angles[0][0])

    if new_angles[1] != []:
        new_b = radians(new_angles[1][0])

    if new_angles[2] != []:
        new_g = radians(new_angles[2][0])

    dg = degrees(gamma - new_g); db = degrees(new_b - beta); da = degrees(new_a - alpha)
    

    
    # rotation 1) rotate a to x

    # defining angle and vector to rotate basis vector a about to get along x axis
    v_ax = [0,a[2],-1.0*a[1]]
    if mag(v_ax) != 0:
        t_ax = degrees(acos(a[0]/mag_a))
        # rotating basis vectors by rotation 1)
        a_ax = rot(a,t_ax,v_ax); a = a_ax
        b_ax = rot(b,t_ax,v_ax); b = b_ax
        c_ax = rot(c,t_ax,v_ax); c = c_ax

    if dg != 0:

        # rotation 2) rotate b,c by dg about bxa
        v_ab = cross(b,a)
        # rotating b,c basis vectors
        b_ab = rot(b,dg,v_ab); b = b_ab
        c_ab = rot(c,dg,v_ab); c = c_ab
    
        # the desired gamma angle should now be locked in
    
    if da != 0:

        # rotation 3) rotate b about a (x) by theta_yz(da)
        mag_byz = sqrt(b[1]**2+b[2]**2)
        mag_cyz = sqrt(c[1]**2+c[2]**2)
        theta_yz = acos((mag_b*mag_c*cos(da) - dot(b,c))/(mag_byz*mag_cyz))
        # rotating b
        b_a = rot(b,theta_yz,a); b = b_a
    
        # the desired alpha angle should now be locked in
    
    if db != 0:
        
        # rotation 4) rotate all basis vectors about b x y such that it aligns along y
        v_by = cross(b,[0,1,0])
        if mag(v_by) != 0:
            t_by = acos(b[1]/mag_b)
            # rotating basis by rotation 4)
            a_by = rot(a,t_by,v_by); a = a_by
            b_by = rot(b,t_by,v_by); b = b_by
            c_by = rot(c,t_by,v_by); c = c_by

        # rotation 5) rotating c about b (y) by theta_xz(db)
        mag_axz = sqrt(a[0]**2+a[2]**2)
        mag_cxz = sqrt(c[0]**2+c[2]**2)
        theta_xz = acos((mag_a*mag_c*cos(db) - dot(a,c))/(mag_axz*mag_cxz))
        # rotating c by 5)
        c_b = rot(c,theta_xz,b); c = c_b
    
        # the desired beta angle should now be locked in

    # rotation 6), 7) realigning basis with x,y,z
    v_ax = [0,a[2],-1.0*a[1]]
    if mag(v_ax) != 0:
        t_ax = degrees(acos(a[0]/mag_a))
        # rotating basis vectors by rotation 6)
        a_ax = rot(a,t_ax,v_ax); a = a_ax
        b_ax = rot(b,t_ax,v_ax); b = b_ax
        c_ax = rot(c,t_ax,v_ax); c = c_ax

    # defining angle and vector to rotate basis vector b about to get b xy-plane, while keeping a along x. (rotate about x)
    v_bxy = [1,0,0]
    t_bxy = degrees(acos(sqrt(b[1]**2+b[2]**2)/mag(b)))
    # rotating basis vectors by rotation 7)
    a_bxy = rot(a,t_bxy,v_bxy); a = a_bxy
    b_bxy = rot(b,t_bxy,v_bxy); b = b_bxy
    c_bxy = rot(c,t_bxy,v_bxy); c = c_bxy

    # Rebuilding basis vector matrix and cdmat matrix
    basis = []
    basis.append(a); basis.append(b); basis.append(c)

    cdmat = build_c_d(basis)

    # conveting back to DC if provided in DC
    if letter == 'D' or letter == 'd':
        points = c_d(points)

    #################################################################################################
    # BEGIN WRITING OUTPUT
    #################################################################################################

    # creating backup files and writing original data to them
    if towhat == 0:
        f = open('CONTCAR_old','w+')
        f.writelines(data)
        f.close()
    elif towhat == 1:
        f = open('POSCAR_old','w+')
        f.writelines(data)
        f.close()

    # recalling if SD was activated
    select = data[7].split()
    sel = select[0][0]

    # replacing old basis vector data with new basis vectors
    data[2] = str('%22.16f %22.16f %22.16f\n' % (a[0], a[1], a[2]))
    data[3] = str('%22.16f %22.16f %22.16f\n' % (b[0], b[1], b[2]))
    data[4] = str('%22.16f %22.16f %22.16f\n' % (c[0], c[1], c[2]))

    # replacing old atom coordinates with new coordinates, including SD flags if necessary
    if sel == 'S' or sel == 's':
        for i in range(num_atoms):
            line = data[i+9].split()
            data[i+9] = str('%20.16f %20.16f %20.16f %4s %4s %4s\n' % (points[i][0], points[i][1], points[i][2], line[3], line[4], line[5]))
    else:
        for i in range(num_atoms):
            data[i+8] = str('%20.16f %20.16f %20.16f\n' % (points[i][0], points[i][1], points[i][2]))

    # writing new data to file with name of original input file
    if towhat == 0 or towhat == 2:
        f = open('CONTCAR','w+')
        f.writelines(data)
        f.close()
        #s.system("open -a 'vesta' CONTCAR")
    elif towhat == 1:
        f = open('POSCAR','w+')
        f.writelines(data)
        f.close()
        #os.system("open -a 'vesta' POSCAR")

#####################################################################################################################
# [3] - SHEAR END ###################################################################################################
#####################################################################################################################



