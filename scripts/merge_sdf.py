import os, sys
import glob as glob
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


path = "/Users/takabak/work/rna_bgsu/minimized"
files = glob.glob(path + "/*.sdf")

mols = []
#mols = [ Chem.MolFromFile(f) for f in files if Chem.MolFromFile(f) != None ]
for f in files:
    basename = os.path.basename(f)
    if not basename.startswith("Triple"):
        m = Chem.MolFromMolFile(f)
        if m != None:
            n = basename.strip('.sdf')
            m.SetProp('_Name', n)
            mols.append(m)
        else:
            print("failed: {}".format(basename))
    
print(len(files))
print(len(mols))

with Chem.SDWriter('../mols.sdf') as w:
    for mol in mols:
        w.write(mol)