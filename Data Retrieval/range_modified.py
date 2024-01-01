import os
from mp_api.client import MPRester

mpr = MPRester('cTAfFUbu7l325RV1')

all_elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K', 'Ar',
                'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
                'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
                'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
                'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
                'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
                'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

# Generate combinations of two elements
element_combinations = [(elem1, elem2) for elem1 in all_elements for elem2 in all_elements]

output_folder = 'MODIFIED_BINARY'
os.makedirs(output_folder, exist_ok=True)

# Track max and min values for various properties
max_values = {'atoms': float('-inf'), 'alpha': float('-inf'), 'beta': float('-inf'), 'gamma': float('-inf'),
              'band_gap': float('-inf'), 'youngs_modulus': float('-inf'), 'a': float('-inf'), 'b': float('-inf'),
              'c': float('-inf'), 'formation_energy_per_atom': float('-inf')}
min_values = {'atoms': float('inf'), 'alpha': float('inf'), 'beta': float('inf'), 'gamma': float('inf'),
              'band_gap': float('inf'), 'youngs_modulus': float('inf'), 'a': float('inf'), 'b': float('inf'),
              'c': float('inf'), 'formation_energy_per_atom': float('inf')}

# Define the range and chunk size
start_index = 0
end_index = len(element_combinations)
chunk_size = 100

# Download CIF files in chunks
for i in range(start_index, end_index, chunk_size):
    chunk_combinations = element_combinations[i:i+chunk_size]

    for elem1, elem2 in chunk_combinations:
        criteria = {'elements': {'$all': [elem1, elem2]}, 'e_above_hull': {'$lt': 0.0}}
        properties = ['material_id', 'pretty_formula', 'nelements', 'alpha', 'beta', 'gamma', 'band_gap', 'elasticity',
                      'volume', 'formation_energy_per_atom']

        compounds = mpr.query(criteria, properties)

        for compound in compounds:
            material_id = compound['material_id']
            atoms = compound['nelements']

            # Update max and min values
            for prop in ['atoms', 'alpha', 'beta', 'gamma', 'band_gap', 'youngs_modulus', 'a', 'b', 'c',
                         'formation_energy_per_atom']:
                value = compound.get(prop, None)
                if value is not None:
                    max_values[prop] = max(max_values[prop], value)
                    min_values[prop] = min(min_values[prop], value)

            # Get the CIF data
            cif_data = mpr.query(criteria={'task_id': material_id}, properties=['cif'])[0]['cif']

            # Save the CIF data to a file in the 'BINARY' folder
            cif_filename = os.path.join(output_folder, f"{material_id}.cif")
            with open(cif_filename, 'w') as cif_file:
                cif_file.write(cif_data)

            print(f"CIF file saved for {material_id}: {cif_filename}")

# Print max and min values
print("Max Values:")
print(max_values)
print("Min Values:")
print(min_values)
