# vasp
Various scripts to manipulate input/output data for/from the VASP program. 

Includes scripts for plotting band structure (no spin, collinear spin), density of states (no spin, collinear spin, 
non-collinear spin), visualizing Kohn Sham Bloch wave functions.




## transform.py 

General purpose code to translate, rotate, stretch, and shear unit cells / atoms in the POSCAR format used by VASP.
See transform.py for more details on functionality.

Most of these transformations can be easily implemented in an excel spreadsheet, but if needs to be done repetetively 
is very useful especially if performing a series of transformations on the cell / atoms.

A regular bash script can be used for input and sucessively call transform.py multiple times.

Can be easily implemented into an automated workflow.



## vasprun.py

Code which parses the vasprun.xml file to collect band structure data and output easy to manipulate matlab objects for plotting.

