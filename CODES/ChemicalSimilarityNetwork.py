#Calculation Tanimoto coeff. for quantifying chemical structure similarity between molecules in a library
#Following analysis requires RDKit python package
#Input is SDF file with muliple chemical structures and output file name
#python3 ChemicalSimilarityNetwork.py <input sdf file> <output file name>
#output file with similarity. col1: molecule1; col2: molecule2; col3: Tanimoto ECFC4; col4: Tanimoto Maccs Key.
#R.P. Vivek-Ananth, IMSc, May 2022

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
import sys

print (rdkit.__version__)


figprnts_ecfc4={}
figprnts_maccs={}
molnames=[]

#reading molecules from sdf file
edcs = Chem.SDMolSupplier(sys.argv[1])# 1 --> input sdf file
print ('Number of molecules {}'.format(len(edcs)))

#calculation of fingerprints
for mol in edcs:
    t1=mol.GetProp('_Name')
    molnames.append(t1)
    if mol is None:
        print('Error: {}'.format(t1))
    else:
        figprnts_ecfc4[t1]=AllChem.GetMorganFingerprint(mol,2)
        figprnts_maccs[t1]=MACCSkeys.GenMACCSKeys(mol)

#calculation of tanimoto coefficient
temp=molnames[:]
for i in molnames:
    temp.remove(i)
    for j in temp:
        tani_ecfc4 = DataStructs.TanimotoSimilarity(figprnts_ecfc4[i],figprnts_ecfc4[j])
        tani_maccs = DataStructs.TanimotoSimilarity(figprnts_maccs[i],figprnts_maccs[j])
        tani = list(map(str,[i,j,tani_ecfc4,tani_maccs]))

        fo_tani.write('\t'.join(tani) + '\n')

fo_tani.close()