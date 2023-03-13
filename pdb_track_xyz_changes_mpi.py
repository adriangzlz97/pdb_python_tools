#!/bin/env/python
import math
import sys
import numpy as np
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
    def __init__(self, atomid, element, altid, restyp, chainid, seqid, x, y, z, xyz_change):
        self.atomid = atomid
        self.element = element
        self.altid = altid
        self.restyp = restyp
        self.chainid = chainid
        self.seqid = seqid
        self.x = x
        self.y = y
        self.z = z
        #self.occ = occ
        #self.biso = biso
        self.xyz_change = xyz_change
    def print_info(self):
        """
        Prints the attributes of each atom tab separated. For testing purposes.
        
        """
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.atomid, self.element, self.altid, self.restyp, self.chainid, self.seqid, self.x, self.y, self.z, self.xyz_change))
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
        if line[:4] == "ATOM":
            line = line.split()
            if line[-1] != "H":
                if len(line[4]) > 2:
                    pdb += [Atom(line[1],line[-1], line[2], line[3], line[4][:1], line[4][1:], float(line[5]), float(line[6]), float(line[7]), 0 )]
                else:
                    pdb += [Atom(line[1],line[-1], line[2], line[3], line[4], line[5], float(line[6]), float(line[7]), float(line[8]), 0 )]
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
                print(atom1.atomid)
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
    atom_p = Atom(0,0,0,0,0,0,0,0,0,0)
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

def compare_pdb_mpi(pdb1,pdb2):
    """
    Compares two lists of Atoms (class). Implements mpi.

    Inputs
    ------
    pdb1, pdb2 : List of Atoms (class)

    Returns
    -------
    Modifies self.xyz_change from pdb1 atoms (class) based on the
    x, y, z change between pdb1 and pdb2.

    """
    # Scatter the lists
    sc_pdb1 = comm.scatter(pdb1, root=0)
    # Iterate through the part of the lists assigned to each rank
    for i in zip(sc_pdb1):
        for atom1 in i:
            for atom2 in pdb2:
                #Make sure it is the same atom being compared
                if atom1.chainid == atom2.chainid:
                    if atom1.seqid == atom2.seqid:
                        if atom1.altid == atom2.altid:
                            #Get coordinates from each atom
                            x1, y1, z1 = atom1.x, atom1.y, atom1.z
                            x2, y2, z2 = atom2.x, atom2.y, atom2.z
                            #Calculate vector distance
                            xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                            xyz = math.sqrt(xyz)
                            #Write distance to attribute
                            atom1.xyz_change = xyz
    # Gather results on rank 0
    pdb1 = comm.gather(sc_pdb1, root=0)
    return pdb1

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    pdb1 = get_atoms_from_pdb(sys.argv[1])
    df_pdb1 = np.array_split(pdb1,size)
    pdb2 = get_atoms_from_pdb(sys.argv[2])
else:
    df_pdb1 = None
    pdb2 = get_atoms_from_pdb(sys.argv[2])
pdb1 = compare_pdb_mpi(df_pdb1,pdb2)
if rank == 0:
    flat_pdb1 = [item for sublist in pdb1 for item in sublist]
    resi_list = find_max_res(flat_pdb1)
    resi_list.sort(key=lambda x: x.xyz_change, reverse=True)
    print("Chain\tResidue\tDistance")
    for i in resi_list:
        if i.xyz_change > 0:
            print("%s\t%s\t%s" % (i.chainid, i.seqid, i.xyz_change))
