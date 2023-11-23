

import os
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.outputs import Poscar

MY_API_KEY =os.getenv("MAPI_KEY")
a = MPRester("xJz5FSRC3KfibTth")

print ("To get the structure of an entry from materialsproject.org,")
print ("in POSCAR format, enter its material_id below.")
mpID = input('material_id: ')
entry = a.get_structure_by_material_id(mpID)

tmp_poscar = Poscar(entry)
tmp_poscar.write_file(mpID + '.POSCAR')
entry.to(filename=mpID + '.cif')
