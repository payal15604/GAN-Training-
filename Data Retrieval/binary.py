import os
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifWriter

api_key = 'LmTAf10d5yXFnT9XyGQ'

mp = MPRester(api_key)

from pymatgen.core.periodic_table import Element
elements = [e.symbol for e in Element]

element_pairs = [(element1, element2) for element1 in elements for element2 in elements]

if not os.path.exists('cif_file_binary'):
    os.makedirs('cif_file_binary')

for element_pair in element_pairs:
    criteria = {'elements': {'$all': element_pair}}
    properties = ['structure']

    results = mp.query(criteria=criteria, properties=properties)

    if results is not None and len(results) > 0:
        print(f"Number of structures found for {element_pair[0]}-{element_pair[1]}: {len(results)}")

        for entry in results:
            structure = entry['structure']
            cif_writer = CifWriter(structure)
            cif_file_path = os.path.join('cif_files', f"{element_pair[0]}_{element_pair[1]}.cif")
            cif_writer.write_file(cif_file_path)
    else:
        print(f"No materials found for {element_pair[0]}-{element_pair[1]}.")

print("CIF files downloaded for all compounds in the 'cif_files' directory.")
