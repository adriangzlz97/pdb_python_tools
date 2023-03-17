#!/bin/env/python
import math
import sys
import numpy as np
from mpi4py import MPI
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
def get_atoms_from_pdb(file, hetatm, hydrogens):
    """
    Parses through a pdb file and generates a list of atoms.

    Inputs
    ------
    String indicating the pdb file to parse.

    Returns
    -------
    List of atoms as Atom class.
    """
    # Read the file per line
    file = open(file, 'r')
    lines = file.readlines()
    pdb = []

    # Iterate through the lines
    for line in lines:
        # Check that the line has atom information
        if line[:4] == "ATOM":
            # Ignore hydrogens by default
            if line[-2] != "H":
                # Add atoms to list by getting the attributes through indexing the line
                pdb += [Atom(line[4:11].strip(),line[-2], line[11:17].strip(), line[17:21].strip(), line[21:22].strip(), line[22:31].strip(), float(line[31:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()), float(line[55:60].strip()), float(line[60:67].strip()), 0 )]
            # Gather hydrogens if -ignore-hydrogens-false flag is present
            elif line[-2] == "H" and hydrogens == "-ignore-hydrogens-false":
                pdb += [Atom(line[4:11].strip(),line[-2], line[11:17].strip(), line[17:21].strip(), line[21:22].strip(), line[22:31].strip(), float(line[31:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()), float(line[55:60].strip()), float(line[60:67].strip()), 0 )]
        # If the -HETATM flag is present, do the same for HETATMs
        if hetatm == "-HETATM":
            if line[:6] == "HETATM":
                if line[-2] != "H":
                    pdb += [Atom(line[4:11].strip(),line[-2], line[11:17].strip(), line[17:21].strip(), line[21:22].strip(), line[22:31].strip(), float(line[31:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()), float(line[55:60].strip()), float(line[60:67].strip()), 0 )]
                elif line[-2] == "H" and hydrogens == "-ignore-hydrogens-false":
                    pdb += [Atom(line[4:11].strip(),line[-2], line[11:17].strip(), line[17:21].strip(), line[21:22].strip(), line[22:31].strip(), float(line[31:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()), float(line[55:60].strip()), float(line[60:67].strip()), 0 )]
        # Ignore lines that do not include the atom information
        else:
            continue
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
    # Iterate through the first atom list
    for atom1 in pdb1:
        # Iterate through the second atom list
        for atom2 in pdb2:
            # Make sure it is the same atom being compared (same chain, seq number and atom name)
            if atom1.chainid == atom2.chainid:
                if atom1.seqid == atom2.seqid:
                    if atom1.altid == atom2.altid:
                        # Get coordinates from each atom
                        x1, y1, z1 = atom1.x, atom1.y, atom1.z
                        x2, y2, z2 = atom2.x, atom2.y, atom2.z
                        # Calculate vector distance
                        xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                        xyz = math.sqrt(xyz)
                        # Write distance to attribute on the first list
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
    # Prepare dummy variables
    atom_p = Atom(0,0,0,0,0,0,0,0,0,0,0,0)
    resi = []
    resi_list_atom = []
    resi_list_max = []
    # Iterate through the atoms in the pdb list
    for atom in pdb:
        # Check that the atom belongs to the same residue as the one before
        if atom.seqid == atom_p.seqid:
            if atom.chainid == atom_p.chainid:
                #Build the residue list with that atom
                resi += [atom]
        # Once it finishes going though all atoms of that residue (it considers that they are correctly ordered)
        else:
            # Double check that it is not empty
            if resi != []:
                # Build a list with the residues
                resi_list_atom += [resi]
                # Start a new residue
                resi = []
        # Reset so it compares with the previous
        atom_p = atom
    # add last residue:
    resi_list_atom += [resi]
    # Order the residue list by distance per residue and keep the max
    for residue in resi_list_atom:
        residue.sort(key=lambda x: x.xyz_change, reverse=True)
        if len(residue) > 0:
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
    # Set up MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # Scatter the lists
    sc_pdb1 = comm.scatter(pdb1, root=0)
    # Iterate through the part of the lists assigned to each rank
    for i in zip(sc_pdb1):
        for atom1 in i:
            # Slowest part - Possibility for improving performance
            for atom2 in pdb2:
                # Make sure it is the same atom being compared
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

def find_contacts(pdb, distance, chain):
    """
    Finds all atoms from different chains within a specific distance and returns a list of pairs.
    """
    atom_pairs = []
    # Iterate through atoms
    for atom1 in pdb:
        # Only consider atoms for a given chain
        if atom1.chainid == chain:
        # Compare with all other atoms. Slowest part - possibility for improvement
            for atom2 in pdb:
                # Consider only atoms from different chain for the comparison
                if atom1.chainid != atom2.chainid:
                    # Get coordinates from each atom
                    x1, y1, z1 = atom1.x, atom1.y, atom1.z
                    x2, y2, z2 = atom2.x, atom2.y, atom2.z
                    # Calculate vector distance
                    xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                    xyz = math.sqrt(xyz)
                    if xyz <= distance:
                        # Ignore duplicate (the opposite pair, reversed)
                        if [atom2,atom1,xyz] not in atom_pairs:
                            atom_pairs += [[atom1,atom2,xyz]]
    return(atom_pairs)

def find_contacts_mpi(pdb, df_pdb, distance, chain):
    """
    Finds all atoms from different chains within a specific distance and returns a list of pairs.
    """
    atom_pairs = []
    # Set up MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # Scatter the lists
    sc_pdb = comm.scatter(df_pdb, root=0)
    # Iterate through the part of the lists assigned to each rank
    for atom1 in pdb:
        # Only consider atoms for a given chain
        if atom1.chainid == chain:
            # Compare with all other atoms. Distributed over mpi. Slowest part - possibility for improvement
            for i in zip(sc_pdb):
                for atom2 in i:
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
    # Remove empty
    if rank == 0:
        while [] in atom_pairs:
                atom_pairs.remove([])
        # Clean list
        flat_atom_pairs = [item for sublist in atom_pairs for item in sublist]
        # Remove duplicates. The opposite pair, which is reversed.
        for i in flat_atom_pairs:
            for j in flat_atom_pairs:
                if [[i[1].atomid,i[0].atomid,i[2]]] == [[j[0].atomid,j[1].atomid,j[2]]]:
                            flat_atom_pairs.remove(i)
        return(flat_atom_pairs)
    
def get_atoms_from_cif(file, hetatm, hydrogens):
    """
    Parses through a cif file and generates a list of atoms.

    Inputs
    ------
    String indicating the pdb file to parse.

    Returns
    -------
    List of atoms as Atom class.
    """
    # Read file by lines
    file = open(file, 'r')
    lines = file.readlines()
    # Set up dummy variable
    cif = []
    # Set up count to get the column order
    count= -1
    # Iterate through the lines
    for line in lines:
        # Add counts each line
        count += 1
        # Reset counts at the start of a loop_, to get the column order
        if "loop_" in line:
            count = -1
        # Get the order for each attribute
        if "_atom_site.id" in line:
            atomid = count
        if "_atom_site.type_symbol" in line:
            element = count
        if "_atom_site.label_atom_id" in line:
            altid = count
        if "_atom_site.label_comp_id" in line:
            restyp = count
        if "_atom_site.auth_asym_id" in line:
            chainid = count
        if "_atom_site.auth_seq_id" in line:
            seqid = count
        if "_atom_site.Cartn_x" in line:
            x = count
        if "_atom_site.Cartn_y" in line:
            y = count
        if "_atom_site.Cartn_z" in line:
            z = count
        if "_atom_site.occupancy" in line:
            occ = count
        if "_atom_site.B_iso" in line:
            biso = count
        # Get atom attributes with the obtained order
        if line[:4] == "ATOM":
            line = line.split()
            if line[2] != "H":
                cif += [Atom(line[atomid], line[element], line[altid], line[restyp], line[chainid], line[seqid], float(line[x]), float(line[y]), float(line[z]), float(line[occ]), float(line[biso]), 0)]
            # Get hydrogens if flag is present
            elif line[2] == "H" and hydrogens == "-ignore-hydrogens-false":
                cif += [Atom(line[atomid], line[element], line[altid], line[restyp], line[chainid], line[seqid], float(line[x]), float(line[y]), float(line[z]), float(line[occ]), float(line[biso]), 0)]
        # Do the same for HETATM if flag is present
        if hetatm == "-HETATM":
            if line[:6] == "HETATM":
                line = line.split()
                if line[2] != "H":
                    cif += [Atom(line[atomid], line[element], line[altid], line[restyp], line[chainid], line[seqid], float(line[x]), float(line[y]), float(line[z]), float(line[occ]), float(line[biso]), 0)]
                elif line[2] == "H" and hydrogens == "-ignore-hydrogens-false":
                    cif += [Atom(line[atomid], line[element], line[altid], line[restyp], line[chainid], line[seqid], float(line[x]), float(line[y]), float(line[z]), float(line[occ]), float(line[biso]), 0)]
        else:
            continue
    return(cif)