# pdb_python_tools
Small tools for analyzing pdb files
The objective of this project is to set up a python script that reads in two pdb files and outputs the atoms/residues sorted by the Å they have been displaced. It will serve as a tool to diagnose structure modeling and make sure unexpected changes are not happening after a global refinement.
If time allows, I might add some other things that would be useful.

I am aware there might be available tools that do similar things, as well as libraries that do some of the things done here. This is a project for learning and some decisions are made with that aim in mind!

## pdb_track_xyz_changes.py Usage
This script will find the maximum xyz coordinate change between two pdb files. The output will be a tabulated table sorted by the residues with the largest coordinate change.
For now, the files should be in pdb format.
### Non-mpi
python pdb_track_xyz_changes.py <filename1> <filename2> > output.txt
### Mpi
mpiexec -n <number of mpi processes> python pdb_track_xyz_changes.py <filename1> <filename2> > output.txt
If you encounter an error about host, try:
mpiexec -n <number of mpi processes> -host localhost python pdb_track_xyz_changes.py <filename1> <filename2> > output.txt
