# pdb_python_tools
Small tools for analyzing pdb/cif files.  
The objective of this project is to set up a python script that reads in two pdb/cif files and outputs the atoms/residues sorted by the Å they have been displaced. It will serve as a tool to diagnose structure modeling and make sure unexpected changes are not happening after a global refinement.
If time allows, I might add some other things that would be useful.

I am aware there might be available tools that do similar things, as well as libraries that do some of the things done here. This is a project for learning and some decisions are made with that aim in mind!  
  
# For usage of any of these tools, run them with -h  

# Requirements
Python3 and the modules argparse, math, numpy (optional), mpi4py (optional)

## track_xyz.py usage
This script will find the maximum xyz coordinate change between two pdb/cif files. The output will be a tabulated table sorted by the residues with the largest coordinate change.
For now, the files should be aligned first with some other program (e.g. ChimeraX)  
To include hetatm (ignored by default) add as an argument after the inputs:  
--HETATM or -HET  
To include hydrogens (ignored by default) add as an argument after the inputs:  
--hydrogens or -hy  
### Non-mpi
python track_xyz.py *filename1 filename2* > output.txt
### Mpi
mpiexec -n *number of mpi processes* python track_xyz.py *filename1 filename2* > output.txt  
If you encounter an error about host, try:  
mpiexec -n *number of mpi processes* -host localhost python track_xyz_mpi.py *filename1 filename2* > output.txt  

## find_contacts.py usage
This script will find all atoms between chains within a certain distance (for a given chain). The output will be a tabulated table following the original pdb atom order.  
To include hetatm (ignored by default) add as an argument after the inputs:  
--HETATM or -HET  
To include hydrogens (ignored by default) add as an argument after the inputs:  
--hydrogens or -hy  
To only include possible polar contacts add as an argument after the inputs:  
--polar_only or -p  
### Non-mpi
python find_contacts.py *filename distance chainid* > output.txt
### Mpi
mpiexec -n *number of mpi processes* python find_contacts_mpi.py *filename1 distance chainid* > output.txt  
If you encounter an error about host, try:  
mpiexec -n *number of mpi processes* -host localhost python find_contacts_mpi.py *filename1 distance chainid* > output.txt  

## analyze_differences.py usage  
his script will find the maximum xyz coordinate change between two pdb/cif files. The output will be a tabulated table sorted by the residues with the largest coordinate change.
For now, the files should be aligned first with some other program (e.g. ChimeraX)  
To include hetatm (ignored by default) add as an argument after the inputs:  
--HETATM or -HET  
To include hydrogens (ignored by default) add as an argument after the inputs:  
--hydrogens or -hy  
It can be run with mpi or not, but mpi4py is required.  

# Test files  
In test_files folder, you can try the test files. They are ribosomes structures, so they are quite big and the program track_xyz will take quite a while (approx 25-30 min) unless you use mpi with a few cores. There are also outputs of what you would get with the command indicated in the .txt files.  
This test is comparing two cryoEM ribosome structures frozen at different time-points. You can track the changes between both structures with the track_xyz program.  
The other test finds all the contacts of the mRNA (chain 4) with other chains with find_contacts. 