import os
from pymatgen.ext.matproj import MPRester

# Set your Materials Project API key here.
api_key = 'LmTAf10d5yXFnT9XyGQ'  # Replace with your actual API key

# Initialize the MPRester with your API key.
mpr = MPRester(api_key)

# List of elements from the periodic table.
elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K', 'Ar',
            'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
            'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
            'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
            'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
            'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
            'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

# Create a directory to save CIF files.
output_directory = 'CIF_Files'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Loop through all possible binary combinations of elements.
for i in range(len(elements)):
    for j in range(i + 1, len(elements)):
        element1 = elements[i]
        element2 = elements[j]

        # Generate the formula of the binary compound.
        formula = f'{element1}{element2}'

        try:
            # Get the CIF data for the binary compound.
            structure = mpr.get_structure_by_material_id(formula)

            # Save the CIF file to the output directory.
            cif_file_path = os.path.join(output_directory, f'{formula}.cif')
            structure.to(filename=cif_file_path, fmt="cif")
            print(f'Downloaded {formula}.cif')
        except Exception as e:
            print(f'Error downloading {formula} CIF: {str(e)}')

print('Downloaded all CIF files.')

