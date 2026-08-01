from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
import sys
import pandas as pd
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import AllChem  # <-- Added for 3D coordinate generation
from pathlib import Path
IPythonConsole.ipython_usePNG = True

df_lead = pd.read_csv('lead compounds.csv', sep=';', header = None)
print("DataFrame Shape:", df_lead.shape)
print("DataFrame Columns:", df_lead.columns.tolist())
lead_count = len(df_lead)
mol_data = []

def mol_summary(lead_name, smiles):

    # convert smiles to mol object
    lead = Chem.MolFromSmiles(smiles)

    if lead is None:
        print(f"Skipping {lead_name}: Invalid SMILES format '{smiles}'")
        return None

    # determine molecular weight
    mol_MW = Descriptors.MolWt(lead)

    # determine logP
    mol_logP = Descriptors.MolLogP(lead)

    # determine HBD
    mol_HBD = Descriptors.NumHDonors(lead)
    
    # determine HBA
    mol_HBA = Descriptors.NumHAcceptors(lead)

    # determine TPSA
    mol_TPSA = Descriptors.TPSA(lead)

    # determine number of rotatable bonds
    mol_RB = Descriptors.NumRotatableBonds(lead)

    # determine ring count
    mol_RC = rdMolDescriptors.CalcNumRings(lead)

    #check if passes rule of 5
    conditions = [mol_MW <= 500, mol_HBA <= 10, mol_HBD <= 5, mol_logP <= 5]
    pass_ro5 = conditions.count(True) >= 3
    ro5_status = "Pass" if pass_ro5 else "Fail"
    if pass_ro5:
        print(f"{lead_name} passes the Rule of 5.")
    else:
        print(f"{lead_name} fails the Rule of 5.")

    data = (
        f"""
        {lead_name}:

        MW: {mol_MW}
        logP: {mol_logP}
        HBD: {mol_HBD}
        HBA: {mol_HBA}
        TPSA: {mol_TPSA}
        Rotatable Bonds: {mol_RB}
        Ring count: {mol_RC}"""
    )

    print(data)

    # --- 3D COORDINATE GENERATION BLOCK ---
    # 1. Hydrogens are mandatory for proper 3D geometry calculations
    lead_3d = Chem.AddHs(lead)
    
    # 2. Embed 3D coordinates using the standard ETKDG toolkit
    embed_status = AllChem.EmbedMolecule(lead_3d, AllChem.ETKDGv3())
    
    # 3. Optimize the structure with an energy force field if embedding succeeds
    if embed_status == 0:
        AllChem.MMFFOptimizeMolecule(lead_3d)
    else:
        # Fallback to a simpler 2D coordinate layout if 3D fails for a complex molecule
        AllChem.Compute2DCoords(lead_3d)
        print(f"Warning: Could not generate 3D coordinates for {lead_name}, falling back to 2D.")
    # --------------------------------------

    return {
        "Name": lead_name,
        "SMILES": smiles,
        "MW": mol_MW,
        "logP": mol_logP,
        "HBD": mol_HBD,
        "HBA": mol_HBA,
        "TPSA": mol_TPSA,
        "Rotatable_Bonds": mol_RB,
        "Ring_Count": mol_RC,
        "Passes RO5's?": ro5_status,
        "_Mol": lead_3d  # Pass the 3D-embedded molecule object to the dataframe
    }

for u in range(0, lead_count):
    summary = mol_summary(df_lead.iloc[u, 0], df_lead.iloc[u, 1])
    mol_data.append(summary)
    
mol_data = [x for x in mol_data if x is not None]

# create dataframe
df = pd.DataFrame(mol_data)

# Step A: View table
print("\n--- ADMET SUMMARY ---")
print(df)

# export to CSV
df.drop(columns=["Name"], errors="ignore").to_csv("ADMET SUMMARY.csv", index=False)

for i in range(0, len(df)):

    # make a new folder
    compound_name = df.loc[i, "Name"]
    path = Path(compound_name)
    path.mkdir(parents=True, exist_ok=True)

    # export to SDF
    path_sdf = path / f"{compound_name}.sdf"
    PandasTools.WriteSDF(df.loc[[i]], str(path_sdf), molColName="_Mol", properties=list(df.columns))

    # create images
    path_png = path / f"{compound_name}.png"
    
    # Generate 2D clean drawing coordinates just for the 2D PNG file slice
    img_mol = Chem.Mol(df.loc[i, "_Mol"])
    AllChem.Compute2DCoords(img_mol)
    img = Draw.MolToImage(img_mol, size=(300, 300), kekulize=True, wedgeBonds=True)
    img.save(path_png)
