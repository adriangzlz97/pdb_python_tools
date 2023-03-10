#!/bin/env/python
import math
class Atom:
    """
    Define atom class.
    
    Attributes
    ----------
    atomid : id for the atom
    element : C, N, O ...
    altid : id within residue
    restyp : residue type 
    chainid : chain id
    seqid : residue number 
    x, y, z : location
    occ : occupancy
    biso : b factor
    xyz_change : movement compared to another pdb
    
    """   
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
    def print_info(self):
        """
        Prints the attributes of each atom tab separated. For testing purposes.
        
        """
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.atomid, self.element, self.altid, self.restyp, self.chainid, self.seqid, self.x, self.y, self.z, self.occ, self.biso, self.xyz_change))
#Function to parse through the pdb file and obtain a list of atoms
def get_atoms_from_pdb(file):
    """
    Parses through a pdb file and generates a list of atoms.

    Inputs
    ------
    String indicating the pdb file to parse.

    Returns
    -------
    List of atoms as Atom class.
    """
    file = open(file, 'r')
    lines = file.readlines()
    pdb = []
    count = 0
    for line in lines:
        if "ATOM" in line:
            line = line.split()
            pdb += [Atom(line[1],line[2], line[3], line[5], line[6], int(line[8]), float(line[10]), float(line[11]), float(line[12]), line[13], line[14],0 )]
            count += 1
    return(pdb)
import random
def displace_xyz(pdb):
    """
    Displaces x, y and z by a random number between 1 and 10. For testing purposes.

    Input
    -----
    Parsed pdb file : list of Atoms (class)
    
    """
    for atom in pdb:
        atom.x += random.randint(1,10)
        atom.y += random.randint(1,10)
        atom.z += random.randint(1,10)
#Define function to compare atom distances and write into the atom attribute
def compare_pdb_xyz(pdb1, pdb2):
    """
    Compares two lists of Atoms (class)

    Inputs
    ------
    pdb1, pdb2 : List of Atoms (class)

    Returns
    -------
    Modifies self.xyz_change from pdb1 atoms (class) based on the
    x, y, z change between pdb1 and pdb2.

    """
    for atom1 in pdb1:
        for atom2 in pdb2:
            #Make sure it is the same atom being compared
            if atom1.atomid == atom2.atomid:
                #Get coordinates from each atom
                x1, y1, z1 = atom1.x, atom1.y, atom1.z
                x2, y2, z2 = atom2.x, atom2.y, atom2.z
                #Calculate vector distance
                xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                xyz = math.sqrt(xyz)
                #Write distance to attribute
                atom1.xyz_change = xyz

def find_max_res(pdb):
    """
    Finds the atom with most xyz_change and returns a list of these atoms sorted by the xyz_change.

    Inputs
    ------
    pdb : List of Atoms (class)

    Returns
    -------
    resi_list_max : list of Atoms from pdb with the largest xyz_change per residue.

    """
    atom_p = Atom(0,0,0,0,0,0,0,0,0,0,0,0)
    resi = []
    resi_list_atom = []
    resi_list_max = []
    for atom in pdb:
        if atom.seqid == atom_p.seqid:
            if atom.chainid == atom_p.chainid:
                resi += [atom]
        #Once it finishes going though all atoms of that residue
        else:
            if resi != []:
                resi_list_atom += [resi]
                resi = []
        #Reset so it compares with the previous
        atom_p = atom
    #add last residue:
    resi_list_atom += [resi]
    # Order the residue list by distance per residue and keep the max
    for residue in resi_list_atom:
        residue.sort(key=lambda x: x.xyz_change, reverse=True)
        resi_list_max += [residue[0]]
    return resi_list_max
    

import copy
file1 = "test.pdb"
pdb1 = get_atoms_from_pdb(file1)
pdb2 = copy.deepcopy(pdb1)
displace_xyz(pdb2)
compare_pdb_xyz(pdb1,pdb2)

resi_list = find_max_res(pdb1)
resi_list.sort(key=lambda x: x.xyz_change, reverse=True)
