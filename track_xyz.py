#!/bin/env/python
import sys
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
from pdb_python_tools import get_atoms_from_cif
from pdb_python_tools import compare_pdb_xyz
from pdb_python_tools import find_max_res
hetatm = 0
hydrogens = 0
# Check for flags
for i in sys.argv:
    if i == "-HETATM":
        hetatm = i
    if i == "-ignore-hydrogens-false":
        hydrogens = i

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
print("Chain\tResidue\tResidue name\tDistance")
for i in resi_list:
    if i.xyz_change > 0.01:
        print("%s\t%s\t%s\t%s" % (i.chainid, i.seqid, i.restyp, i.xyz_change))