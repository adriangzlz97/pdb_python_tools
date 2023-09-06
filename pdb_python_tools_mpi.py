#!/bin/env/python
import math
import numpy as np
from mpi4py import MPI
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
                        if atom1.altid == atom2.altid and isinstance(atom1.xyz_change, int):
                            #Get coordinates from each atom
                            x1, y1, z1 = atom1.x, atom1.y, atom1.z
                            x2, y2, z2 = atom2.x, atom2.y, atom2.z
                            #Calculate vector distance
                            xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                            xyz = math.sqrt(xyz)
                            #Write distance to attribute
                            atom1.xyz_change = xyz
                        if atom1.restyp == "TYR" or atom1.restyp == "PHE":
                            if "CE" in atom1.altid or "CD" in atom1.altid:
                                if "CE" in atom2.altid or "CD" in atom2.altid:
                                    x1, y1, z1 = atom1.x, atom1.y, atom1.z
                                    x2, y2, z2 = atom2.x, atom2.y, atom2.z
                                    xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                                    xyz = math.sqrt(xyz)
                                    if xyz < atom1.xyz_change or isinstance(atom1.xyz_change, int):
                                        atom1.xyz_change = xyz

    # Gather results on rank 0
    pdb1 = comm.gather(sc_pdb1, root=0)
    return pdb1

def compare_resi_pdb_mpi(pdb1,pdb2):
    """
    Compares two lists of Residues (class). Implements mpi.

    Inputs
    ------
    pdb1, pdb2 : List of Residues (class)

    Returns
    -------
    Modifies self.xyz_change from pdb1 Atoms within the list in the Residue (class) based on the
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
        for resi1 in i:
            for resi2 in pdb2:
                if resi1.chainid == resi2.chainid and resi1.seqid == resi2.seqid:
                    for atom1 in resi1.atom_list:
                # Slowest part - Possibility for improving performance
                        for atom2 in resi2.atom_list:
                        # Make sure it is the same atom being compared
                            if atom1.chainid == atom2.chainid:
                                if atom1.seqid == atom2.seqid:
                                    if atom1.altid == atom2.altid and isinstance(atom1.xyz_change, int):
                                        #Get coordinates from each atom
                                        x1, y1, z1 = atom1.x, atom1.y, atom1.z
                                        x2, y2, z2 = atom2.x, atom2.y, atom2.z
                                        #Calculate vector distance
                                        xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                                        xyz = math.sqrt(xyz)
                                        #Write distance to attribute
                                        atom1.xyz_change = xyz
                                    if atom1.altid == "CA" or atom1.altid == "C1'":
                                        resi1.CA.xyz_change = xyz
                                    if atom1.restyp == "TYR" or atom1.restyp == "PHE":
                                        if "CE" in atom1.altid or "CD" in atom1.altid:
                                            if "CE" in atom2.altid or "CD" in atom2.altid:
                                                x1, y1, z1 = atom1.x, atom1.y, atom1.z
                                                x2, y2, z2 = atom2.x, atom2.y, atom2.z
                                                xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                                                xyz = math.sqrt(xyz)
                                                if xyz < atom1.xyz_change or isinstance(atom1.xyz_change, int):
                                                    atom1.xyz_change = xyz
                                    elif atom1.restyp == "GLU":
                                        if "OE1" in atom1.altid or "OE2" in atom1.altid:
                                            x1, y1, z1 = atom1.x, atom1.y, atom1.z
                                            x2, y2, z2 = atom2.x, atom2.y, atom2.z
                                            xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                                            xyz = math.sqrt(xyz)
                                            if xyz < atom1.xyz_change or isinstance(atom1.xyz_change, int):
                                                atom1.xyz_change = xyz
                                    elif atom1.restyp == "ASP":
                                        if "OD1" in atom1.altid or "OD2" in atom1.altid:
                                            x1, y1, z1 = atom1.x, atom1.y, atom1.z
                                            x2, y2, z2 = atom2.x, atom2.y, atom2.z
                                            xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                                            xyz = math.sqrt(xyz)
                                            if xyz < atom1.xyz_change or isinstance(atom1.xyz_change, int):
                                                atom1.xyz_change = xyz
                                    elif atom1.restyp == "ARG":
                                        if "NH1" in atom1.altid or "NH2" in atom1.altid:
                                            x1, y1, z1 = atom1.x, atom1.y, atom1.z
                                            x2, y2, z2 = atom2.x, atom2.y, atom2.z
                                            xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                                            xyz = math.sqrt(xyz)
                                            if xyz < atom1.xyz_change or isinstance(atom1.xyz_change, int):
                                                atom1.xyz_change = xyz

    # Gather results on rank 0
    pdb1 = comm.gather(sc_pdb1, root=0)
    return pdb1

def find_contacts_mpi(pdb, df_pdb, distance, chain, polar):
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
            if polar == True:
                if atom1.element == "O" or atom1.element == "N" or atom1.element == "P" or atom1.element == "S":
                    for i in zip(sc_pdb):
                        for atom2 in i:
                            if atom1.chainid != atom2.chainid:
                                if atom2.element == "O" or atom2.element == "N" or atom2.element == "P" or atom2.element == "S":
                                    #Get coordinates from each atom
                                    x1, y1, z1 = atom1.x, atom1.y, atom1.z
                                    x2, y2, z2 = atom2.x, atom2.y, atom2.z
                                    #Calculate vector distance
                                    xyz = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                                    xyz = math.sqrt(xyz)
                                    if xyz <= distance:
                                        atom_pairs += [[atom1,atom2,xyz]]
            # Compare with all other atoms. Distributed over mpi. Slowest part - possibility for improvement
            else:
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
  
def find_contacts_resi_mpi(pdb, df_pdb, distance, chain, polar):
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
    for resi1 in pdb:
        # Only consider atoms for a given chain
        if resi1.chainid == chain:
            if polar == True:
                for i in zip(sc_pdb):
                    for resi2 in i:
                        # Consider only atoms from different chain for the comparison
                        if resi1.chainid != resi2.chainid:
                            #iterate over residue
                            for atom1 in resi1.atom_list:
                                if atom1.element == "O" or atom1.element == "N" or atom1.element == "P" or atom1.element == "S":
                                    for atom2 in resi2.atom_list:
                                        if atom2.element == "O" or atom2.element == "N" or atom2.element == "P" or atom2.element == "S":
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
            # Compare with all other atoms. Distributed over mpi. Slowest part - possibility for improvement
            else:
                for i in zip(sc_pdb):
                    for resi2 in i:
                        # Consider only atoms from different chain for the comparison
                        if resi1.chainid != resi2.chainid:
                            #iterate over residue
                            for atom1 in resi1.atom_list:
                                for atom2 in resi2.atom_list:
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