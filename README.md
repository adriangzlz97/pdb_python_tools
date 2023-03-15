# pdb_python_tools
Small tools for analyzing pdb/cif files.  
The objective of this project is to set up a python script that reads in two pdb/cif files and outputs the atoms/residues sorted by the Å they have been displaced. It will serve as a tool to diagnose structure modeling and make sure unexpected changes are not happening after a global refinement.
If time allows, I might add some other things that would be useful.

I am aware there might be available tools that do similar things, as well as libraries that do some of the things done here. This is a project for learning and some decisions are made with that aim in mind!  

# Requirements
Python3 and the modules sys, math, numpy (optional), mpi4py (optional)

## track_xyz.py usage
This script will find the maximum xyz coordinate change between two pdb/cif files. The output will be a tabulated table sorted by the residues with the largest coordinate change.
For now, the files should be aligned first with some other program (e.g. ChimeraX)  
To include hetatm (ignored by default) add as a argument after the inputs:  
-HETATM  
To include hydrogens (ignored by default) add as a argument after the inputs:  
-ignore-hydrogens-false  
### Non-mpi
python track_xyz.py *filename1 filename2* > output.txt
### Mpi
mpiexec -n *number of mpi processes* python track_xyz.py *filename1 filename2* > output.txt  
If you encounter an error about host, try:  
mpiexec -n *number of mpi processes* -host localhost python track_xyz_mpi.py *filename1 filename2* > output.txt  

## find_contacts.py usage
This script will find all atoms between chains within a certain distance (for a given chain). The output will be a tabulated table following the original pdb atom order.  
To include hetatm (ignored by default) add as a argument after the inputs:  
-HETATM  
To include hydrogens (ignored by default) add as a argument after the inputs:  
-ignore-hydrogens-false  
### Non-mpi
python find_contacts.py *filename distance chainid* > output.txt
### Mpi
mpiexec -n *number of mpi processes* python find_contacts_mpi.py *filename1 distance chainid* > output.txt  
If you encounter an error about host, try:  
mpiexec -n *number of mpi processes* -host localhost python find_contacts_mpi.py *filename1 distance chainid* > output.txt  
