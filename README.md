# vasp
Various scripts to manipulate input/output data for/from the VASP program. 

Includes scripts for plotting band structure (no spin, collinear spin), density of states (no spin, collinear spin, 
non-collinear spin), visualizing Kohn Sham Bloch wave functions.


## Atomic Positions and Unit Cell

### transform.py 

General purpose code to translate, rotate, stretch, and shear unit cells / atoms in the POSCAR format used by VASP.
See transform.py for more details on functionality.

Most of these transformations can be easily implemented in an excel spreadsheet, but if needs to be done repetetively 
is very useful especially if performing a series of transformations on the cell / atoms.

A regular bash script can be used for input and sucessively call transform.py multiple times.

Can be easily implemented into an automated workflow.


## Spectral Properties

### vasprun.py

Code which parses the vasprun.xml file to collect band structure data and output easy to manipulate matlab objects for plotting.



### bands_xml.m

Matlab script which imports .mat files generated by vasprun.py for plotting band structure. Requires editing by user to get
plots with the desired data.


### plotDOS.m and plotDOS.mlx

Matlab script and live script version used to plot density of states from the DOSCAR file.


## Real Space

### psinks.m

Matlab script used to plot real space wave functions from a text file created from the WAVECAR file.

### locpot.m

Simple Matlab script which is used to plot the real space local potential from the LOCPOT file.








