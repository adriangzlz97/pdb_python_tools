#!/bin/env/python
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
from pdb_python_tools import get_atoms_from_cif
from pdb_python_tools import find_contacts
import argparse

# Check for flags
parser = argparse.ArgumentParser(
                    prog='find_contacts.py',
                    description='Find possible contacts between chains for a given chain and within a given distance',
                    epilog='Usage: pdb1/cif1 -arguments')
parser.add_argument('pdb', help='coordinate file (pdb/cif)')
parser.add_argument('--chain', help='chain id to analyze', required=True)
parser.add_argument('--distance', help='distance to check',type=float, required=True)
parser.add_argument('-HET','--HETATM', action='store_true', dest='hetatm', help='include hetatms')
parser.add_argument('-hy','--hydrogens', action='store_true', dest='hydrogens', help='include hydrogens')
parser.add_argument('-p','--polar_only', action='store_true', dest='polar', help='check only polar')
args = parser.parse_args()
pdb = args.pdb
chain = args.chain
distance = args.distance
hetatm = args.hetatm
hydrogens = args.hydrogens
polar = args.polar

# Check format and parse with appropriate function
if ".pdb" in pdb:
    pdb = get_atoms_from_pdb(pdb, hetatm, hydrogens)
elif ".cif" in pdb:
    pdb = get_atoms_from_cif(pdb, hetatm, hydrogens)

# Find the contacts within that distance
atom_pairs = find_contacts(pdb, distance, chain, polar)

# Print table
print("Chain1\tResidue1\tResidue1 number\tAtom1\tChain2\tResidue2\tResidue2 number\tAtom2\tDistance")
for i in atom_pairs:
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (i[0].chainid, i[0].restyp, i[0].seqid, i[0].altid, i[1].chainid, i[1].restyp, i[1].seqid, i[1].altid, i[2]))