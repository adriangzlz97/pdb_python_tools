#!/bin/env/python
import sys
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
import math
import numpy as np
from mpi4py import MPI
from pdb_python_tools import find_contacts_mpi

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
hetatm = 0
if len(sys.argv) >= 4:
    hetatm = sys.argv[3]
pdb = get_atoms_from_pdb(sys.argv[1],hetatm)
if rank == 0:
    df_pdb = np.array_split(pdb,size)
    
else:
    df_pdb = None
distance = float(sys.argv[2])
atom_pairs = find_contacts_mpi(pdb, df_pdb, distance)
if rank == 0:
    print("Chain1\tResidue1\tResidue1 number\tAtom1\tChain2\tResidue2\tResidue2 number\tAtom2\tDistance")
    for i in atom_pairs:
        if len(i) == 3:
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (i[0].chainid, i[0].restyp, i[0].seqid, i[0].altid, i[1].chainid, i[1].restyp, i[1].seqid, i[1].altid, i[2]))