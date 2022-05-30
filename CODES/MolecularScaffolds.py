#Compute molecular scaffolds for chemical structures in a library
#Following analysis requires RDKit python package
#MurckoScaffold.py code of RDKit package was edited to compute Scaffold at Graph/Node level
#Replace the default MurckoScaffold.py code in RDKit with the provided edited code before executing this code
#Input is SDF file with muliple chemical structures and output file name
#python3 MolecularScaffolds.py <input sdf file> <output file name>
#output file with similarity. col1: molecule1; col2: molecule2; col3: Tanimoto ECFC4; col4: Tanimoto Maccs Key.
#R.P. Vivek-Ananth, IMSc, May 2022

import rdkit, gzip
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import sys

print ('RDkit version: ',rdkit.__version__)

def carboncounter(mol):
        return len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#6]"), maxMatches=mol.GetNumAtoms()))

suppl = Chem.SDMolSupplier(sys.argv[1]) # sdf files with many structures
print (len(suppl))

fout = open(sys.argv[2],'w',encoding='utf-8') # output file name

fout.write('\t'.join(['Chemical identifier','Scaffold-Graph','Scaffold-Graph_Node','Scaffold-Graph_Node_Bond'])+'\n')

for Counter,mol in enumerate(suppl):
    if mol is None:
        print ('Error',Counter)
        continue
    else:
        if Counter%50000 == 0:
            print (Counter)
        iden = mol.GetProp('_Name')
        numcarbons = str(carboncounter(mol))
        core = MurckoScaffold.GetScaffoldForMol(mol)
        sf1 = Chem.MolToSmiles(core,isomericSmiles=False,canonical=True)
        try:
            sf3 = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(core),isomericSmiles=False,canonical=True)
            sf2 = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGenericOnlyBonds(core),isomericSmiles=False,canonical=True)
            if '.' not in sf1 and ('.' not in sf3 and '.' not in sf2):
                fout.write(str(iden) + '\t' + sf3 + '\t' + sf2 + '\t' + sf1 + '\t' + numcarbons + '\n')
            else:
                fout.write(str(iden) + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + numcarbons + '\n')
        except:
            print ('Any Exception', iden)
            if '.' not in sf1:
                fout.write(str(iden) + '\t' + '' + '\t' + '' + '\t' + sf1 + '\t' + numcarbons + '\n')
            else:
                fout.write(str(iden) + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + numcarbons + '\n')

fout.close()
