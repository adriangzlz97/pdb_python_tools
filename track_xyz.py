#!/bin/env/python
import sys
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
from pdb_python_tools import get_atoms_from_cif
from pdb_python_tools import compare_pdb_xyz
from pdb_python_tools import find_max_res
from pdb_python_tools import Residue
from pdb_python_tools import get_resi_from_cif
from pdb_python_tools import get_resi_from_pdb
from pdb_python_tools import compare_pdb_resi_xyz
hetatm = 0
hydrogens = 0
info = False
# Check for flags
for i in sys.argv:
    if i == "-HETATM":
        hetatm = i
    if i == "-ignore-hydrogens-false":
        hydrogens = i
    if i == "-more_info":
        info = True
if info == False:
    # Check format and parse with appropriate function
    if ".pdb" in sys.argv[1]:
        pdb1 = get_atoms_from_pdb(sys.argv[1], hetatm, hydrogens)
    elif ".cif" in sys.argv[1]:
        pdb1 = get_atoms_from_cif(sys.argv[1], hetatm, hydrogens)
    if ".pdb" in sys.argv[2]:
        pdb2 = get_atoms_from_pdb(sys.argv[2], hetatm, hydrogens)
    elif ".cif" in sys.argv[2]:
        pdb2 = get_atoms_from_cif(sys.argv[2], hetatm, hydrogens)

    # Compare both pdbs 
    compare_pdb_xyz(pdb1,pdb2)

    # Find maximum xyz change per residue and produce a list of those atoms
    resi_list = find_max_res(pdb1)

    # Sort by the xyz change
    resi_list.sort(key=lambda x: x.xyz_change, reverse=True)

    # Print table
    print("Chain\tResidue\tResidue name\tDistance\tAtom")
    for i in resi_list:
        if i.xyz_change > 0.01:
            print("%s\t%s\t%s\t%s\t%s" % (i.chainid, i.seqid, i.restyp, i.xyz_change, i.altid))

if info == True:
    # Check format and parse with appropriate function
    if ".pdb" in sys.argv[1]:
        pdb1 = get_atoms_from_pdb(sys.argv[1], hetatm, hydrogens)
    elif ".cif" in sys.argv[1]:
        pdb1 = get_atoms_from_cif(sys.argv[1], hetatm, hydrogens)
    if ".pdb" in sys.argv[2]:
        pdb2 = get_atoms_from_pdb(sys.argv[2], hetatm, hydrogens)
    elif ".cif" in sys.argv[2]:
        pdb2 = get_atoms_from_cif(sys.argv[2], hetatm, hydrogens)

    # Compare both pdbs 
    compare_pdb_resi_xyz(pdb1,pdb2)
    for resi in pdb1:
        resi.max_xyz = max(resi.atom_list, key=lambda x: x.xyz_change)
        resi.average_xyz = sum(atom.xyz_change for atom in resi.atom_list) / len(resi.atom_list)

    pdb1.sort(key=lambda x: x.max_xyz.xyz_change, reverse=True)
    # Print table
    print("Chain\tResidue\tResidue name\tMax_Distance\tMax_atom\tAverage_distance\tCA/C1'_distance")
    for resi in pdb1:
        if resi.max_xyz.xyz_change > 0.01:
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (resi.chainid, resi.seqid, resi.restyp, resi.max_xyz.xyz_change, resi.max_xyz.altid, resi.average_xyz, resi.CA_xyz))