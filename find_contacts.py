#!/bin/env/python
import sys
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
from pdb_python_tools import get_atoms_from_cif
import math
from pdb_python_tools import find_contacts

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

# Gather chainid
chain = sys.argv[3]

# Check format and parse with appropriate function
if ".pdb" in sys.argv[1]:
    pdb = get_atoms_from_pdb(sys.argv[1], hetatm, hydrogens)
elif ".cif" in sys.argv[1]:
    pdb = get_atoms_from_cif(sys.argv[1], hetatm, hydrogens)

# Gather max distance to find
distance = float(sys.argv[2])

# Find the contacts within that distance
atom_pairs = find_contacts(pdb, distance, chain, polar)

# Print table
print("Chain1\tResidue1\tResidue1 number\tAtom1\tChain2\tResidue2\tResidue2 number\tAtom2\tDistance")
for i in atom_pairs:
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (i[0].chainid, i[0].restyp, i[0].seqid, i[0].altid, i[1].chainid, i[1].restyp, i[1].seqid, i[1].altid, i[2]))