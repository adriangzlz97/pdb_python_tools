#!/bin/env/python
import sys
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
from pdb_python_tools import compare_pdb_mpi
from pdb_python_tools import find_max_res
import numpy as np
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
hetatm = 0
if len(sys.argv) >= 4:
    hetatm = sys.argv[3]
pdb2 = get_atoms_from_pdb(sys.argv[2], hetatm)
if rank == 0:
    pdb1 = get_atoms_from_pdb(sys.argv[1], hetatm)
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
        if i.xyz_change > 0:
            print("%s\t%s\t%s\t%s" % (i.chainid, i.seqid, i.restyp, i.xyz_change))