import os
from time import sleep
from mp_api.client import MPRester
import logging

def extract_properties(compound):
    material_id = compound['material_id']
    atoms = compound['nelements']
    alpha = compound['elasticity']['G_Reuss']
    beta = compound['elasticity']['G_VRH']
    gamma = compound['elasticity']['G_Voigt']
    band_gap = compound['band_gap']
    a = compound['volume'] / compound['nelements']
    b = a
    c = a
    formation_energy_per_atom = compound['formation_energy_per_atom']

    # update the maximum and minimum values
    global max_values, min_values
    max_values['atoms'] = max(max_values['atoms'], atoms)
    max_values['alpha'] = max(max_values['alpha'], alpha)
    max_values['beta'] = max(max_values['beta'], beta)
    max_values['gamma'] = max(max_values['gamma'], gamma)
    max_values['band_gap'] = max(max_values['band_gap'], band_gap)
    max_values['a'] = max(max_values['a'], a)
    max_values['b'] = max(max_values['b'], b)
    max_values['c'] = max(max_values['c'], c)
    max_values['formation_energy_per_atom'] = max(max_values['formation_energy_per_atom'], formation_energy_per_atom)

    min_values['atoms'] = min(min_values['atoms'], atoms)
    min_values['alpha'] = min(min_values['alpha'], alpha)
    min_values['beta'] = min(min_values['beta'], beta)
    min_values['gamma'] = min(min_values['gamma'], gamma)
    min_values['band_gap'] = min(min_values['band_gap'], band_gap)
    min_values['a'] = min(min_values['a'], a)
    min_values['b'] = min(min_values['b'], b)
    min_values['c'] = min(min_values['c'], c)
    min_values['formation_energy_per_atom'] = min(min_values['formation_energy_per_atom'], formation_energy_per_atom)

    # store the data in a text file
    with open(os.path.join(output_folder, '{}.txt'.format(material_id)), 'w') as f:
        f.write('material_id: {}\n'.format(material_id))
        f.write('pretty_formula: {}\n'.format(compound['pretty_formula']))
        f.write('atoms: {}\n'.format(atoms))
        f.write('alpha: {}\n'.format(alpha))
        f.write('beta: {}\n'.format(beta))
        f.write('gamma: {}\n'.format(gamma))
        f.write('band_gap: {}\n'.format(band_gap))
        f.write('a: {}\n'.format(a))
        f.write('b: {}\n'.format(b))
        f.write('c: {}\n'.format(c))
        f.write('formation_energy_per_atom: {}\n'.format(formation_energy_per_atom))

def make_api_request(mpr, criteria, properties):
    max_retries = 3
    retries = 0
    while retries < max_retries:
        try:
            compounds = mpr.materials.summary.search(criteria, properties)
            return compounds
        except Exception as e:
            logging.error(f"Error in API request: {e}")
            retries += 1
            sleep(2 ** retries)  # Exponential backoff before retrying

    logging.error(f"Failed after {max_retries} retries. Check API key, network, and retry later.")
    raise Exception("Failed after multiple retries. Check logs for details.")

if __name__ == '__main__':
    try:
        # Replace 'your_api_key' with your actual Materials Project API key
        mpr = MPRester("iUNJ21RiTNEJUMudW6WXt855qtLvABfx")

        output_folder = 'MODIFIED_BINARY'
        os.makedirs(output_folder, exist_ok=True)

        max_values = {'atoms': float('-inf'), 'alpha': float('-inf'), 'beta': float('-inf'), 'gamma': float('-inf'),
                     'band_gap': float('-inf'), 'a': float('-inf'), 'b': float('-inf'), 'c': float('-inf'),
                     'formation_energy_per_atom': float('-inf')}
        min_values = {'atoms': float('inf'), 'alpha': float('inf'), 'beta': float('inf'), 'gamma': float('inf'),
                     'band_gap': float('inf'), 'a': float('inf'), 'b': float('inf'), 'c': float('inf'),
                     'formation_energy_per_atom': float('inf')}

        all_elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K', 'Ar',
                        'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
                        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
                        'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
                        'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
                        'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
                        'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

        # Generate combinations of two elements
        element_combinations = [(elem1, elem2) for elem1 in all_elements for elem2 in all_elements]

        for elem1, elem2 in element_combinations:
            # Use a list for range criteria values
            criteria = {
                'elements': {'$all': [elem1, elem2]},
                'e_above_hull': {'$gte': 0.0, '$lte': 0.2}  # Update for the specified range
            }
            properties = ['material_id', 'pretty_formula', 'nelements', 'elasticity', 'band_gap', 'volume', 'formation_energy_per_atom']

            compounds = make_api_request(mpr, criteria, properties)

            for compound in compounds:
                extract_properties(compound)

        # Print max and min values
        print("Max Values:")
        print(max_values)
        print("Min Values:")
        print(min_values)

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
