from rdkit import Chem
from rdkit.Chem import Draw

mol_smiles = 'NC(C(=O)O)C1CC(=O)N1CC(CC(=O)O)C2CC(=O)N2CC(=O)O'
mol_obj = Chem.MolFromSmiles(mol_smiles)

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
    matches_highlighted = []
    for name, obj in dict.items():
        if obj is None:
            print(f"Skipping {name}: Invalid SMARTS pattern")
            continue
            
        # Check if the molecule matches the compiled pattern
        if mol.HasSubstructMatch(obj):
            matches = mol.GetSubstructMatches(obj)
            matches_highlighted.extend([num for sub_tuple in matches for num in sub_tuple])
            num_matches = len(matches)
            if num_matches == 1:
                print(f"Contains {num_matches} {name} group")
            else:
                print(f"Contains {num_matches} {name} groups")

    img = Draw.MolToImage(mol, highlightAtoms = matches_highlighted)
    img.save('image.png')

# Run the fixed function
find_groups(mol_obj, compiled_patterns)