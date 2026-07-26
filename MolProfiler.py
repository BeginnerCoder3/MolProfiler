from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
import pubchempy as pcp
import pandas as pd
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import PandasTools
IPythonConsole.ipython_usePNG = True

def mol_summary(mol_name):

    # get the compound from PubChem
    mol = pcp.get_compounds(mol_name, 'name')

    # get the name of the compound
    molecule = mol[0].synonyms[0]

    # determine molecular weight
    mol_MW = mol[0].molecular_weight

    # determine logP
    mol_logP = mol[0].xlogp

    # determine HBD
    mol_HBD = mol[0].h_bond_donor_count

    # determine HBA
    mol_HBA = mol[0].h_bond_acceptor_count

    # determine TPSA
    mol_TPSA = mol[0].tpsa

    # determine number of rotatable bonds
    mol_RB = mol[0].rotatable_bond_count

    # determine ring count
    smiles = mol[0].smiles
    rdkit_mol = Chem.MolFromSmiles(smiles)
    mol_RC = rdMolDescriptors.CalcNumRings(rdkit_mol)

    data = (
        f"""
        {molecule}:

        MW: {mol_MW}
        logP: {mol_logP}
        HBD: {mol_HBD}
        HBA: {mol_HBA}
        TPSA: {mol_TPSA}
        Rotatable Bonds: {mol_RB}
        Ring count: {mol_RC}"""
        
    )

    print(data)

    return {
        "Name": mol_name,
        "SMILES": smiles,
        "MW": mol[0].molecular_weight,
        "logP": mol[0].xlogp,
        "HBD": mol[0].h_bond_donor_count,
        "HBA": mol[0].h_bond_acceptor_count,
        "TPSA": mol[0].tpsa,
        "Rotatable_Bonds": mol[0].rotatable_bond_count,
        "Ring_Count": rdMolDescriptors.CalcNumRings(rdkit_mol),
        "_Mol": rdkit_mol  # Needed by WriteSDF to draw 2D structures in the file
    }

df_name = pd.read_csv('lead compounds.csv')
lead_count = len(df_name)
mol_data = []

for u in range (0, lead_count):
    summary = mol_summary(df_name.iloc[u, 0])
    mol_data.append(summary)
    u = u + 1

# create dataframe
df = pd.DataFrame(mol_data)

# Step A: View table
print("\n--- MOLECULAR TABLE ---")
print(df.drop(columns=["_Mol"]))

# export to CSV
df.drop(columns=["Name"]).to_csv("molecules.csv", index=False)

for i in range(0, lead_count):

    # export to SDF
    PandasTools.WriteSDF(df.loc[[i]], f"""{df.loc[i, "Name"]}.sdf""", molColName="_Mol", properties=list(df.columns))

    # create images
    img = Draw.MolToImage(df.loc[i, "_Mol"], size=(300, 300), kekulize=True, wedgeBonds=True)
    img.save(f"""{df.loc[i, "Name"]}.png""")
    i = i + 1