#!/bin/env/python
import sys
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
from pdb_python_tools import get_atoms_from_cif
from pdb_python_tools import compare_pdb_mpi
from pdb_python_tools import find_max_res
from pdb_python_tools import get_resi_from_cif
from pdb_python_tools import compare_resi_pdb_mpi
from pdb_python_tools import get_resi_from_pdb
from pdb_python_tools import Residue
import numpy as np
from mpi4py import MPI
# Set up MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
hetatm = 0
hydrogens = 0
# Check for flags
for i in sys.argv:
    if i == "-HETATM":
        hetatm = i
    if i == "-ignore-hydrogens-false":
        hydrogens = i


# Check format and parse with appropriate function
if ".pdb" in sys.argv[2]:
    pdb2 = get_resi_from_pdb(sys.argv[2], hetatm, hydrogens)
elif ".cif" in sys.argv[2]:
    pdb2 = get_resi_from_cif(sys.argv[2], hetatm, hydrogens)
# Check format of the other pdb, parse and split in equal chunks for the mpi processes
if rank == 0:
    if ".pdb" in sys.argv[1]:
        pdb1 = get_resi_from_pdb(sys.argv[1], hetatm, hydrogens)
    elif ".cif" in sys.argv[1]:
        pdb1 = get_resi_from_cif(sys.argv[1], hetatm, hydrogens)
    df_pdb1 = np.array_split(pdb1,size)
    
else:
    df_pdb1 = None

# Compare both pdbs through mpi    
pdb1 = compare_resi_pdb_mpi(df_pdb1,pdb2)
# Root process does the rest
if rank == 0:
    # Find maximum xyz change per residue and produce a list of those atoms
    flat_pdb1 = [item for sublist in pdb1 for item in sublist]
    for resi in flat_pdb1:
        resi.max_xyz = max(resi.atom_list, key=lambda x: x.xyz_change)
        resi.average_xyz = sum(atom.xyz_change for atom in resi.atom_list) / len(resi.atom_list)

    flat_pdb1.sort(key=lambda x: x.max_xyz.xyz_change, reverse=True)
    # Print table
    print("Chain\tResidue\tResidue name\tMax_Distance\tMax_atom\tAverage_distance\tCA/C1'_distance")
    for resi in flat_pdb1:
        if resi.max_xyz.xyz_change > 0.01:
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (resi.chainid, resi.seqid, resi.restyp, resi.max_xyz.xyz_change, resi.max_xyz.altid, resi.average_xyz, resi.CA_xyz))