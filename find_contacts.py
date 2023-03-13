#!/bin/env/python
import sys
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
import math

def find_contacts(pdb, distance):
    """
    Finds all atoms from different chains within a specific distance and returns a list of pairs.
    """
    atom_pairs = []
    for atom1 in pdb:
        for atom2 in pdb:
            if atom1.chainid != atom2.chainid:
                #Get coordinates from each atom
                x1, y1, z1 = atom1.x, atom1.y, atom1.z
                x2, y2, z2 = atom2.x, atom2.y, atom2.z
                #Calculate vector distance
                xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                xyz = math.sqrt(xyz)
                if xyz <= distance:
                    #Ignore duplicate
                    if [atom2,atom1,xyz] not in atom_pairs:
                        atom_pairs += [[atom1,atom2,xyz]]
    return(atom_pairs)


pdb = get_atoms_from_pdb(sys.argv[1])
distance = float(sys.argv[2])
atom_pairs = find_contacts(pdb, distance)
print("Chain1\tResidue1\tResidue1 number\tAtom1\tChain2\tResidue2\tResidue2 number\tAtom2\tDistance")
for i in atom_pairs:
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (i[0].chainid, i[0].restyp, i[0].seqid, i[0].altid, i[1].chainid, i[1].restyp, i[1].seqid, i[1].altid, i[2]))