#!/bin/env/python
import sys
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
from pdb_python_tools import get_atoms_from_cif
import math
import numpy as np
from mpi4py import MPI
from pdb_python_tools import find_contacts_mpi

# Set up MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
hetatm = 0
hydrogens = 0
polar = False

# Check for flags
for i in sys.argv:
    if i == "-HETATM":
        hetatm = i
    if i == "-ignore-hydrogens-false":
        hydrogens = i
    if i == "-polar_only":
        polar = True

# Check format and parse with appropriate function
if ".pdb" in sys.argv[1]:
    pdb = get_atoms_from_pdb(sys.argv[1], hetatm, hydrogens)
elif ".cif" in sys.argv[1]:
    pdb = get_atoms_from_cif(sys.argv[1], hetatm, hydrogens)

# Gather chainid
chain = sys.argv[3]

# Distribute over mpi processes
if rank == 0:
    df_pdb = np.array_split(pdb,size)  
else:
    df_pdb = None

# Gather max distance to consider
distance = float(sys.argv[2])

# Find the contacts within that distance through mpi
atom_pairs = find_contacts_mpi(pdb, df_pdb, distance, chain, polar)

# Print table on root
if rank == 0:
    print("Chain1\tResidue1\tResidue1 number\tAtom1\tChain2\tResidue2\tResidue2 number\tAtom2\tDistance")
    for i in atom_pairs:
        if len(i) == 3:
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (i[0].chainid, i[0].restyp, i[0].seqid, i[0].altid, i[1].chainid, i[1].restyp, i[1].seqid, i[1].altid, i[2]))