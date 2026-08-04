from rdkit import Chem
from rdkit.Chem import Draw
import itertools

mol_smiles = 'NCC(N)COCC(=O)NC1CC(=O)N1Cc2c(C(=O)O)c(C(=O)O)c(C(=O)O)c(C(=O)O)c2O'
mol_obj = Chem.MolFromSmiles(mol_smiles)
palette = [(1, 0, 0), (1, 0.541, 0), (1, 0.984, 0), (0, 1, 0.071), (0, 0.384, 1), (0.604, 0, 1)] #roygb

funky_dict = {
    "beta-lactam" : "N1C(=O)CC1",
    "sulfonamide" : "Nc1ccc(cc1)S(=O)(=O)N",
    "carboxylic acid" : "[CX3](=O)[OX2H1]",
    "primary amine" : "[NX3H2]",
    "ester" : "[#6][CX3](=O)[OX2H0][#6]",
    "aldehyde" : "[CX3H1](=O)"
}

# Compile patterns correctly from SMARTS
compiled_patterns = {name: Chem.MolFromSmarts(s) for name, s in funky_dict.items()}

def find_groups(mol, dict):
    # Directly iterate through key-value pairs (no index tracking needed)
    atoms_highlighted = {}
    bonds_highlighted = {}
    group_colors = {}
    for (name, obj), color in zip(dict.items(), palette):
        group_colors[name] = color
    print(group_colors)
    for name, obj in dict.items():
        if obj is None:
            print(f"Skipping {name}: Invalid SMARTS pattern")
            continue
            
        # Check if the molecule matches the compiled pattern
        if mol.HasSubstructMatch(obj):
            matches = mol.GetSubstructMatches(obj)
    
            for sub_tuple in matches:
                for atom_idx in sub_tuple:
                    atoms_highlighted[atom_idx] = group_colors[name]
    
            for match in matches:
                for a_1, a_2 in itertools.combinations(match, 2):
                    bond = mol.GetBondBetweenAtoms(a_1, a_2)
                    if bond is not None:
                        bonds_highlighted[bond.GetIdx()] = group_colors[name]
    
            num_matches = len(matches)
            if num_matches == 1:
                print(f"Contains {num_matches} {name} group")
            else:
                print(f"Contains {num_matches} {name} groups")
            #print(name, "atoms:", [idx for sub_tuple in matches for idx in sub_tuple])
    img = Draw.MolToImage(
    mol,
    highlightAtoms=list(atoms_highlighted.keys()),
    highlightBonds=list(bonds_highlighted.keys()),
    highlightAtomColors=atoms_highlighted,
    highlightBondColors=bonds_highlighted
    )    
    img.save('image.png')

# Run the fixed function
find_groups(mol_obj, compiled_patterns)