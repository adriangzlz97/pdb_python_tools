#!/bin/env/python

"""
#Define atom class:
    Attributes: atomid (id for the atom), element (C, N, O ...), altid (id within residue), restyp (residue type), 
    chainid (chain id), seqid (residue number), x, y, z (location), occ (occupancy), biso (b factor), xyz_change (movement compared to another pdb)
    
"""    
class Atom:
    def __init__(self, atomid, element, altid, restyp, chainid, seqid, x, y, z, occ, biso, xyz_change):
        self.atomid = atomid
        self.element = element
        self.altid = altid
        self.restyp = restyp
        self.chainid = chainid
        self.seqid = seqid
        self.x = x
        self.y = y
        self.z = z
        self.occ = occ
        self.biso = biso
        self.xyz_change = xyz_change
#Function to parse through the pdb file and obtain a list of atoms
def get_atoms_from_pdb(file):
    lines = file.readlines()
    pdb = []
    count = 0
    for line in lines:
        if "ATOM" in line:
            line = line.split()
            pdb += [Atom(line[1],line[2], line[3], line[5], line[6], line[7], float(line[10]), float(line[11]), float(line[12]), line[13], line[14],0 )]
            count += 1
    return(pdb)
def displace_xyz(pdb, distance):
    for atom in pdb:
        atom.x += distance
        atom.y += distance
        atom.z += distance


file1 = open("test.pdb", 'r')
pdb1 = get_atoms_from_pdb(file1)
print(pdb1[0].x)
displace_xyz(pdb1,1)
print(pdb1[0].x)