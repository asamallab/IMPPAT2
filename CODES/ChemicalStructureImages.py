#Create SVG or PNG images for chemical structures in a library
#Following analysis requires RDKit python package
#python3 ChemicalStructureImages.py <input sdf file> <output file format: svg or png>
#R.P. Vivek-Ananth, IMSc, May 2022

from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Draw
import rdkit

import sys

print (rdkit.__version__)


print (sys.argv[1], sys.argv[2])

#reading molecules from sdf file arg 1
edcs = Chem.SDMolSupplier(sys.argv[1])
#print ('Number of molecules {}'.format(len(edcs)))
if len(edcs) != 1:
    print ('Error!!!')
elif len(edcs) == 1:
    mol = edcs[0]
    t1=mol.GetProp('_Name')
    outn = sys.argv[1].split('/')[-1].replace('.sdf','.'+sys.argv[2])
    print (outn)
    
    if mol is None:
        print('Error: {}'.format(t1))
    else:
        #create images of the 2D molecules        
        Draw.MolToFile(mol,outn, size=(400, 400), fitImage=True)
