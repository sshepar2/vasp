
'''
    >> python /path/to/vasprun.py

    Use the command above with the vasprun.xml file in the current directory. This
    will parse the vasprun.xml file and create a series of .mat objects which can 
    then be loaded using Matlab for easy plotting.
    
    vasprun.xml is chosen over other VASP output files which contain band structure data
    because it includes information not available in the others. The vasprun.xml is also 
    already formatted for parsing. The drawback is the file tends to be very large.

    Output files in order as they are output:
        matbands.mat   - dim(bands,kpts), ordered eigenvalues
        matkpts.mat    - dim(kpts), properly scaled, unique k-points.
        fermien.mat    - contains the fermi energy.
        ispin.mat      - contains information on whether spin polarization is included.
        symlines.mat   - contains the values of kpts at which high-symmetry points occur. 
                        Used for plotting vertical bars at these points.    
        matpbands.mat  - contains the atomic / orbital angular momentum projected
                        eigenvalue data.
        matbands1.mat  - same as matbands.mat except for the second spin component.
        matpbands1.mat - same as matpbands except for the second spin component.

'''

from xml.dom import minidom
import xml.etree.ElementTree as ET
from math import sqrt
import numpy as np
import scipy.io

tree = ET.parse("vasprun.xml")
root = tree.getroot()

vasp = minidom.parse("vasprun.xml")

gen = vasp.getElementsByTagName("generator")[0]

gen_i = gen.getElementsByTagName("i")

#for i in range(len(gen_i)):
#    info = gen_i[i].firstChild.data
#    print(info)

structure_1 = vasp.getElementsByTagName("structure")[0]
crystal_1 = structure_1.getElementsByTagName("crystal")[0]
rec_basis = crystal_1.getElementsByTagName("varray")[1]
rec_vec = rec_basis.getElementsByTagName("v")
rec_1 = rec_basis.getElementsByTagName("v")[0]
rec_2 = rec_basis.getElementsByTagName("v")[1]
rec_3 = rec_basis.getElementsByTagName("v")[2]

b1 = rec_1.firstChild.data.split()
for i in range(len(b1)):
    b1[i] = float(b1[i])

b2 = rec_2.firstChild.data.split()
for i in range(len(b2)):
    b2[i] = float(b2[i])

b3 = rec_3.firstChild.data.split()
for i in range(len(b3)):
    b3[i] = float(b3[i])

#getting number of divisions along each path
div = root.findall(".//*[@param='listgenerated']/i")
divs = int(div[0].text)

#getting number of paths
sym = root.findall(".//*[@param='listgenerated']/v")
paths = len(sym)-1

#getting number of paths
efermi = root.findall(".//*[@name='efermi']")
fermien = float(efermi[0].text)

#getting number of bands
num_band = root.findall(".//*[@name='NBANDS']")
num_bands = int(num_band[0].text)

#getting ISPIN tag
ispin = root.findall(".//*[@name='ISPIN']")
spin = int(ispin[0].text)

#getting number of ions
ion = root.findall("./atominfo/atoms")
ions = int(ion[0].text)

#getting ion info
ioninfo = root.findall("./atominfo/array/set/rc/c")
ionsinfo = []
for i in range(len(ioninfo)):
    ionsinfo.append(ioninfo[i].text)

#getting all kpoints list for which eigenvalues were generated
kpts = root.findall(".//*[@name='kpointlist']/v")
kmat = []
for i in range(len(kpts)):
    k = kpts[i].text.split()
    kmat.append(k)

kall = len(kmat)
repeated = paths-1
ktot = kall-repeated

#converting all entries from str to float
for i in range(len(kmat)):
    for j in range(3):
        kmat[i][j] = float(kmat[i][j])

#removing repeated kpoint entries

for i in range(len(kmat)-1):
    if kmat[i] == kmat[i+1]:
        kmat[i+1] = []

for i in range(len(kmat)):
    kmat = [x for x in kmat if x != []]

#getting eigenvalues
eigen = root.findall("./calculation/eigenvalues/array/set/set/set/r")
eigenlist = []
for i in range(len(eigen)):
    E = eigen[i].text.split()
    eigenlist.append(E[0])

#converting eigenvalues to floats
for i in range(len(eigenlist)):
    eigenlist[i] = float(eigenlist[i])

#getting projected bands (s py pz px dxy dyz dz2 dxz x2-y2) for 9 orbitals
peigen = root.findall("./calculation/projected/array/set/set/set/set/r")
peigenlist = []
for i in range(len(peigen)):
    P = peigen[i].text.split()
    peigenlist.append(P)

#converting str's in peigenlist to floats
for i in range(len(peigenlist)):
    for j in range(len(peigenlist[i])):
        peigenlist[i][j] = float(peigenlist[i][j])

#incase spin is included
eigenlist1 = []
peigenlist1 = []
if spin == 1:
    print('non-spin-polarized')
elif spin == 2:
    eigenlist1 = eigenlist[len(eigenlist)//2:]
    eigenlist = eigenlist[:len(eigenlist)//2]
    peigenlist1 = peigenlist[len(peigenlist)//2:]
    peigenlist = peigenlist[:len(peigenlist)//2]
    print('spin-polarized')
    
#removing duplicates in peigenlist
for i in range(repeated):
    for j in range(ions*num_bands):
        peigenlist[(i+1)*ions*num_bands*divs+j] = []
peigenlist = [x for x in peigenlist if x !=[]]


#reordering peigenlist ([ion][orbital][band][kpoint])
pbandmat = []
for i in range(ions):
    orbital = []
    for o in range(len(peigenlist[i])):
        bands = []
        for b in range(num_bands):
            values = []
            for k in range(ktot):
                values.append(peigenlist[i+ions*num_bands*k+ions*b][o])
            bands.append(values)
        orbital.append(bands)
    pbandmat.append(orbital)

#reordering eigenvalues into bands
bandmat = []

for j in range(num_bands):
    band_j = []
    for i in range(kall):
        val = eigenlist[i*num_bands+j]
        band_j.append(val)
    bandmat.append(band_j)

#removing duplicates in bandmat
for j in range(len(bandmat)):
    for i in range(repeated):
        bandmat[j][i*divs+divs] = ''
        
for i in range(len(bandmat)):
    bandmat[i] = [x for x in bandmat[i] if x !='']


'''Now we have a list of lists, with each list corresponding to a band (k vs E).
The lists are as long as there are unique (removed duplicates) k-points.  The
list is as long as there are bands.'''

'''Example:  3 paths specified in KPOINTS, with 30 divisions along each path.
This is a total of 90 kpoints.  Since there are 3 connected paths there will be
 2 repeated kpoints next to eachother where 2 paths meet.  That means there are
 88 unique kpoints.  If vasp calculated 16 bands then the dimensions of bandmat
 would be a list of 16 lists, with 88 eigenvalues in each of those lists.'''

#Making list to plot with eigenvalues for correct incrementation
k_plot = [0]
for i in range(len(kmat)-1):
    db1 = (kmat[i][0]-kmat[i+1][0])
    db2 = (kmat[i][1]-kmat[i+1][1])
    db3 = (kmat[i][2]-kmat[i+1][2])
    b1db1 = []
    b2db2 = []
    b3db3 = []
    for j in range(3):
        b1db1.append(db1*b1[j])
        b2db2.append(db2*b2[j])
        b3db3.append(db3*b3[j])
    val = sqrt((b1db1[0]+b2db2[0]+b3db3[0])**2+(b1db1[1]+b2db2[1]+b3db3[1])**2+(b1db1[2]+b2db2[2]+b3db3[2])**2)
    k_plot.append(val+k_plot[i])


highsym = []
for i in range(len(sym)):
    highsym.append(k_plot[i*divs-i])
    print('High symmetry point at: ' + str((k_plot[i*divs-i])))

#making matlab object of dim(bands,kpoints)
matbands = np.asarray(bandmat)
scipy.io.savemat('/mnt/c/Users/Stuart Shepard/Documents/MATLAB/matbands.mat',mdict={'matbands' : matbands})

#making matlab object of kpoint spacing
matkpts = np.asarray(k_plot)
scipy.io.savemat('/mnt/c/Users/Stuart Shepard/Documents/MATLAB/matkpts.mat',mdict={'matkpts' :matkpts})

#making matlab object of fermi energy
fermi = np.asarray(fermien)
scipy.io.savemat('/mnt/c/Users/Stuart Shepard/Documents/MATLAB/fermi.mat',mdict={'fermi' :fermi})

#making matlab object of spin
ispin = np.asarray(spin)
scipy.io.savemat('/mnt/c/Users/Stuart Shepard/Documents/MATLAB/ispin.mat',mdict={'ispin' :ispin})

#making matlab object of vertical lines as high symmetry points
symlines = np.asarray(highsym)
scipy.io.savemat('/mnt/c/Users/Stuart Shepard/Documents/MATLAB/symlines.mat',mdict={'symlines' :symlines})

#making matlab object of projected band structure
matpbands = np.asarray(pbandmat)
scipy.io.savemat('/mnt/c/Users/Stuart Shepard/Documents/MATLAB/matpbands.mat',mdict={'matpbands' :matpbands})


if spin == 2:
    #reordering eigenvalues of second spin component into lists of bands
    bandmat1 = []
    for j in range(num_bands):
        band_j = []
        for i in range(kall):
            val = eigenlist1[i*num_bands+j]
            band_j.append(val)
        bandmat1.append(band_j)
    #changing duplicates entries to ''
    for j in range(len(bandmat1)):
        for i in range(repeated):
            bandmat1[j][i*divs+divs] = ''
    #removing '' entries
    for i in range(len(bandmat1)):
        bandmat1[i] = [x for x in bandmat1[i] if x !='']
    #making matlab object of dim(bands,kpoints)
    matbands1 = np.asarray(bandmat1)
    scipy.io.savemat('/mnt/c/Users/Stuart Shepard/Documents/MATLAB/matbands1.mat',mdict={'matbands1' :matbands1})
    #spin component 2 of projected bands
    for i in range(repeated):
        for j in range(ions*num_bands):
            peigenlist1[(i+1)*ions*num_bands*divs+j] = []
    peigenlist1 = [x for x in peigenlist1 if x !=[]]
    pbandmat1 = []
    for i in range(ions):
        orbital = []
        for o in range(len(peigenlist[i])):
            bands = []
            for b in range(num_bands):
                values = []
                for k in range(ktot):
                    values.append(peigenlist1[i+ions*num_bands*k+ions*b][o])
                bands.append(values)
            orbital.append(bands)
        pbandmat1.append(orbital)
    matpbands1 = np.asarray(pbandmat1)
    scipy.io.savemat('/mnt/c/Users/Stuart Shepard/Documents/MATLAB/matpbands1.mat',mdict={'matpbands1' :matpbands1})
    




    

