#!/bin/env/python
import sys
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
import math
from pdb_python_tools import find_contacts

hetatm = 0
hydrogens = 0
for i in sys.argv:
    if i == "-HETATM":
        hetatm = i
    if i == "-ignore-hydrogens-false":
        hydrogens = i
chain = sys.argv[3]
pdb = get_atoms_from_pdb(sys.argv[1], hetatm, hydrogens)
distance = float(sys.argv[2])
atom_pairs = find_contacts(pdb, distance, chain)
print("Chain1\tResidue1\tResidue1 number\tAtom1\tChain2\tResidue2\tResidue2 number\tAtom2\tDistance")
for i in atom_pairs:
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (i[0].chainid, i[0].restyp, i[0].seqid, i[0].altid, i[1].chainid, i[1].restyp, i[1].seqid, i[1].altid, i[2]))