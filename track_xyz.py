#!/bin/env/python
import sys
from pdb_python_tools import Atom
from pdb_python_tools import get_atoms_from_pdb
from pdb_python_tools import compare_pdb_xyz
from pdb_python_tools import find_max_res



pdb1 = get_atoms_from_pdb(sys.argv[1])
pdb2 = get_atoms_from_pdb(sys.argv[2])

compare_pdb_xyz(pdb1,pdb2)

resi_list = find_max_res(pdb1)
resi_list.sort(key=lambda x: x.xyz_change, reverse=True)
print("Chain\tResidue\tResidue name\tDistance")
for i in resi_list:
    if i.xyz_change > 0:
        print("%s\t%s\t%s\t%s" % (i.chainid, i.seqid, i.restyp, i.xyz_change))