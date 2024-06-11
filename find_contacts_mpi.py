#!/usr/bin/env python3
from pdb_python_tools import Atom
from pdb_python_tools import Residue
from pdb_python_tools import get_resi_from_pdb
from pdb_python_tools import get_resi_from_cif
import argparse
import numpy as np
from mpi4py import MPI
from pdb_python_tools_mpi import find_contacts_resi_mpi

# Set up MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
hetatm = 0
hydrogens = 0
polar = False
simple = False

# Check for flags
parser = argparse.ArgumentParser(
                    prog='find_contacts_mpi.py',
                    description='Find possible contacts between chains for a given chain and within a given distance',
                    epilog='Usage: pdb1/cif1 -arguments')
parser.add_argument('pdb', help='coordinate file (pdb/cif)')
parser.add_argument('--chain', help='chain id to analyze', required=True)
parser.add_argument('--distance', help='distance to check',type=float, required=True)
parser.add_argument('-HET','--HETATM', action='store_true', dest='hetatm', help='include hetatms')
parser.add_argument('-hy','--hydrogens', action='store_true', dest='hydrogens', help='include hydrogens')
parser.add_argument('-p','--polar_only', action='store_true', dest='polar', help='check only polar')
parser.add_argument('-s','--simple', action='store_true', dest='simple', help='simple output: just display the residues involved and closest distance')
args = parser.parse_args()
pdb = args.pdb
chain = args.chain
distance = args.distance
hetatm = args.hetatm
hydrogens = args.hydrogens
polar = args.polar
simple = args.simple

# Check format and parse with appropriate function
if ".pdb" in pdb:
    pdb = get_resi_from_pdb(pdb, hetatm, hydrogens)
elif ".cif" or ".mmcif" in pdb:
    pdb = get_resi_from_cif(pdb, hetatm, hydrogens)

# Distribute over mpi processes
if rank == 0:
    df_pdb = np.array_split(pdb,size)  
else:
    df_pdb = None

# Find the contacts within that distance through mpi
atom_pairs = find_contacts_resi_mpi(pdb, df_pdb, distance, chain, polar)

# Print table on root
if rank == 0:
    if simple:
        n=0
        for i in atom_pairs:
            if n == 0:
                atom_pairs_simple = [i]
                n+=1
            else:
                count = -1
                for j in atom_pairs_simple:
                    check = False
                    count += 1
                    if i[0].seqid == j[0].seqid:
                        if i[1].seqid == j[1].seqid:
                            check = True
                            if i[2] < j[2]:
                                atom_pairs_simple[count][2] = i[2]
                                break
                if not check:
                    atom_pairs_simple.append(i)
        print("Residue1\tResidue1 number\tChain2\tResidue2\tResidue2 number\tDistance")
        for i in atom_pairs_simple:
            print("%s\t%s\t%s\t%s\t%s\t%s" % (i[0].restyp,i[0].seqid,i[1].chainid,i[1].restyp,i[1].seqid,i[2]))
    else:
        # Print table
        print("Chain1\tResidue1\tResidue1 number\tAtom1\tChain2\tResidue2\tResidue2 number\tAtom2\tDistance")
        for i in atom_pairs:
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (i[0].chainid, i[0].restyp, i[0].seqid, i[0].altid, i[1].chainid, i[1].restyp, i[1].seqid, i[1].altid, i[2]))