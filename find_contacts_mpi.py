#!/bin/env/python
import sys
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
import math
import numpy as np
def find_contacts_mpi(pdb, df_pdb, distance):
    """
    Finds all atoms from different chains within a specific distance and returns a list of pairs.
    """
    atom_pairs = []
    # MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # Scatter the lists
    sc_pdb = comm.scatter(df_pdb, root=0)
    # Iterate through the part of the lists assigned to each rank
    for i in zip(sc_pdb):
        for atom1 in i:
            for atom2 in pdb:
                if atom1.chainid != atom2.chainid:
                    #Get coordinates from each atom
                    x1, y1, z1 = atom1.x, atom1.y, atom1.z
                    x2, y2, z2 = atom2.x, atom2.y, atom2.z
                    #Calculate vector distance
                    xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                    xyz = math.sqrt(xyz)
                    if xyz <= distance:
                       atom_pairs += [[atom1,atom2,xyz]]
    # Gather results on rank 0
    atom_pairs = comm.gather(atom_pairs, root=0)
    # Remove empty and duplicates
    if rank == 0:
        while [] in atom_pairs:
                atom_pairs.remove([])
        for i in atom_pairs:
            for j in atom_pairs:
                if len(i[0]) == 3:
                    if [[i[0][1].atomid,i[0][0].atomid,i[0][2]]] == [[j[0][1].atomid,j[0][0].atomid,j[0][2]]]:
                        atom_pairs.remove(i)
        #Clean list
        n = 0
        for i in atom_pairs:
            atom_pairs[n] = i[0]
            n += 1
        print(atom_pairs)
        return(atom_pairs)
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
pdb = get_atoms_from_pdb(sys.argv[1])
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