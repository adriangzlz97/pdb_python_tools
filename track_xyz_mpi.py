#!/bin/env/python
import sys
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
from pdb_python_tools import get_atoms_from_cif
from pdb_python_tools import compare_pdb_mpi
from pdb_python_tools import find_max_res
import numpy as np
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
hetatm = 0
hydrogens = 0
for i in sys.argv:
    if i == "-HETATM":
        hetatm = i
    if i == "-ignore-hydrogens-false":
        hydrogens = i
if ".pdb" in sys.argv[2]:
    pdb2 = get_atoms_from_pdb(sys.argv[2], hetatm, hydrogens)
elif ".cif" in sys.argv[2]:
    pdb2 = get_atoms_from_cif(sys.argv[2], hetatm, hydrogens)
if rank == 0:
    if ".pdb" in sys.argv[1]:
        pdb1 = get_atoms_from_pdb(sys.argv[1], hetatm, hydrogens)
    elif ".cif" in sys.argv[1]:
        pdb1 = get_atoms_from_cif(sys.argv[1], hetatm, hydrogens)
    df_pdb1 = np.array_split(pdb1,size)
    
else:
    df_pdb1 = None
    
pdb1 = compare_pdb_mpi(df_pdb1,pdb2)
if rank == 0:
    flat_pdb1 = [item for sublist in pdb1 for item in sublist]
    resi_list = find_max_res(flat_pdb1)
    resi_list.sort(key=lambda x: x.xyz_change, reverse=True)
    print("Chain\tResidue\tResidue name\tDistance")
    for i in resi_list:
        if i.xyz_change > 0.0009:
            print("%s\t%s\t%s\t%s" % (i.chainid, i.seqid, i.restyp, i.xyz_change))