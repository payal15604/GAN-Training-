import os
from pymatgen.ext.matproj import MPRester

# Replace 'YOUR_API_KEY' with your Materials Project API key
mpr = MPRester('WZ4TpJFvnVDhdszu')

# Get a list of all elements in the Materials Project
all_elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K', 'Ar',
            'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
            'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
            'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
            'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
            'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
            'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

# Generate combinations of two elements
element_combinations = [(elem1, elem2) for elem1 in all_elements for elem2 in all_elements]

output_folder = 'BINARY'
os.makedirs(output_folder, exist_ok=True)

# Download CIF files for compounds with two elements
for elem1, elem2 in element_combinations:
    criteria = {'elements': {'$all': [elem1, elem2]}}
    properties = ['material_id', 'pretty_formula']

    compounds = mpr.query(criteria, properties)

    for compound in compounds:
        material_id = compound['material_id']

        # Check if the CIF file already exists
        cif_filename = os.path.join(output_folder, f"{material_id}.cif")
        if os.path.exists(cif_filename):
            print(f"CIF file already exists for {material_id}: {cif_filename}")
            continue  # Skip to the next iteration if the file exists

        # Get the CIF data
        cif_data = mpr.query(criteria={'task_id': material_id}, properties=['cif'])[0]['cif']

        # Save the CIF data to a file in the 'BINARY' folder
        with open(cif_filename, 'w') as cif_file:
            cif_file.write(cif_data)

        print(f"CIF file saved for {material_id}: {cif_filename}")
