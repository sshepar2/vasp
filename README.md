# vasp
Various scripts to manipulate input/output data for/from the VASP program. 

Includes scripts for plotting band structure (no spin, spin), density of states (no spin, spin, spin-orbit coupling),




## transform.py 


General purpose code to translate, rotate, stretch, and shear unit cells / atoms in the POSCAR format used by VASP.
See transform.py for more details on functionality.

Most of these transformations can be easily implemented in an excel spreadsheet, but if needs to be done repetetively 
is very useful especially if performing a series of transformations on the cell / atoms.

A regular bash script can be used for input and sucessively call transform.py multiple times.

Can be implemented into an automated workflow easily.


#############
# #
#############
