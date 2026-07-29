from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
import pubchempy as pcp
import sys
import pandas as pd
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors, Lipinski
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
    print(f"{lead_name} passes the Rule of 5's.")

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
        "_Mol": lead  # Needed by WriteSDF to draw 2D structures in the file
    }

for u in range (0, lead_count):
    summary = mol_summary(df_lead.iloc[u, 0], df_lead.iloc[u, 1])
    mol_data.append(summary)
    u = u + 1

# create dataframe
df = pd.DataFrame(mol_data)

# Step A: View table
print("\n--- ADMET SUMMARY ---")
print(df)

# export to CSV
df.drop(columns=["Name"], errors="ignore").to_csv("ADMET SUMMARY.csv", index=False)

for i in range(0, lead_count):

    # export to SDF
    PandasTools.WriteSDF(df.loc[[i]], f"""{df.loc[i, "Name"]}.sdf""", molColName="_Mol", properties=list(df.columns))

    # create images
    img = Draw.MolToImage(df.loc[i, "_Mol"], size=(300, 300), kekulize=True, wedgeBonds=True)
    img.save(f"""{df.loc[i, "Name"]}.png""")
    i = i + 1