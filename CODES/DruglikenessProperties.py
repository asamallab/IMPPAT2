#Evaluate drug-likeness of chemical structures in a library
#Following analysis requires RDKit python package
#The output from MolecularProperties.py is used as input
#python3 DruglikenessProperties.py <output file from MolecularProperties.py> <input sdf file> <output file name>
#R.P. Vivek-Ananth, IMSc, May 2022

import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import QED

print ('RDkit version: ',rdkit.__version__)

def lipinski_ro5(mwt,logp,hbd,hba):
    lst = [mwt < 500, logp < 5, hbd <= 5,hba <= 10]
    vio_c = lst.count(False)
    if vio_c == 0 or vio_c == 1:
        return int(vio_c), 'Passed'
    else:
        return int(vio_c), 'Failed'

def ghose(logp,mr,mwt,tna):
    lst = [logp >= -0.4 and logp <= 5.6, mr >= 40 and mr <= 130, mwt >= 160 and mwt <= 480, tna >= 20 and tna <= 70]
    vio_c = lst.count(False)
    if vio_c == 0:
        return int(vio_c), 'Passed'
    else:
        return int(vio_c), 'Failed'

def veber(tpsa,nrb,hbd,hba):
    lst1 = [(tpsa <= 140 or hbd+hba <= 12) and nrb <= 10]
    if all(lst1):
        return 'Good'
    else:
        return 'Bad'

def egan(tpsa,logp):
    lst = [logp >= -1 and logp <= 6, tpsa >= 0 and tpsa <= 132]
    if all(lst):
        return 'Good'
    else:
        return 'Bad'

def gsk4_400(mwt,logp):
    lst = [mwt < 400, logp < 4]
    if all(lst):
        return 'Good'
    else:
        return 'Bad'

def pfizer3_75(logp,tpsa):
    lst = [logp < 3, tpsa > 75]
    if all(lst):
        return 'Good'
    else:
        return 'Bad'

def qed(mol):
    qed_w = QED.weights_mean(mol)
    qed_uw = QED.weights_none(mol)
    return qed_w,qed_uw


suppl = Chem.SDMolSupplier(sys.argv[2]) # sdf files with many structures
print (len(suppl))

mol_dict = {}
for Counter,mol in enumerate(suppl):
    if mol is None:
        print ('Error',Counter)
        continue
    else:
        iden = mol.GetProp('_Name')
        mol_dict[iden] = mol
print (len(mol_dict))

fout = open(sys.argv[3],'w') # output file name

fout.write('\t'.join(['Identifier', 'Lipinski\'s RO5 violations','Lipinski filter','Ghose filter violations','Ghose filter','Veber filter','Egan filter','GSK 4/400 filter','Pfizer 3/75 filter','QEDw','QEDuw','Fraction of SP3 hybridized Carbon atoms','Shape complexity','Stereochemical complexity'])+'\n')

for lin in open(sys.argv[1]): # input computed physicochemical properties
    if 'Identifier' in lin:
        continue
    tmp = lin.strip().split('\t')
    iden = tmp[0]
    mwt = float(tmp[2])
    logp = float(tmp[4])
    mr = float(tmp[5])
    tpsa = float(tmp[6])
    hbd = int(tmp[7])
    hba = int(tmp[8])
    tna = int(mol_dict[iden].GetNumAtoms(onlyExplicit=False))
    nrb = int(tmp[22])
    Fsp3 = float(tmp[16])
    shape_compx = float(tmp[15])
    fcc = float(tmp[17])
    qed_w, qed_uw = qed(mol_dict[iden])
    l_vio, l_s = lipinski_ro5(mwt,logp,hbd,hba)
    g_vio, g_s = ghose(logp,mr,mwt,tna)
    v_s = veber(tpsa,nrb,hbd,hba)
    e_s = egan(tpsa,logp)
    gsk_s = gsk4_400(mwt,logp)
    p_s = pfizer3_75(logp,tpsa)
    fout.write('\t'.join([iden,str(l_vio),str(l_s),str(g_vio),str(g_s),str(v_s),str(e_s),str(gsk_s),str(p_s),str(qed_w),str(qed_uw),str(Fsp3),str(shape_compx),str(fcc)])+'\n')

fout.close()
