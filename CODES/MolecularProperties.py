#Compute physicochemical properties for chemical structures in a library
#Following analysis requires RDKit python package
#python3 MolecularProperties.py <sdf or smiles> <input structure file> <output file name>
#R.P. Vivek-Ananth, IMSc, May 2022

from rdkit import Chem
from rdkit.Chem import Descriptors
import sys

def counter(mol,pat):
        return len(mol.GetSubstructMatches(pat, maxMatches=mol.GetNumAtoms()))

inptype = sys.argv[1]# either 'smiles' or 'sdf'

if inptype == 'smiles':

    suppl = Chem.SmilesMolSupplier(sys.argv[2],delimiter='\t',titleLine=False)#input file name

elif inptype == 'sdf':

    suppl = Chem.SDMolSupplier(sys.argv[2])
    print (len(suppl))

else:

    print ('Error', 'Please provide file type: `smiles` or `sdf`')

f = open(sys.argv[3],'w')#output file name

f.write('\t'.join(['Identifier','Chemical formula','Mol.wt','Exact Mol.Mwt','LogP','Molar Refractivity','Topological polar surface area ','Number of Hydrogen bond donors','Number of Hydrogen bond acceptors','Number of Chiral Carbons','Number of Carbon atoms','Number of SP hybridized Carbon atoms','Number of SP2 hybridized Carbon atoms','Number of SP2 hybridized Carbon atoms test','Number of SP3 hybridized Carbon atoms','Shape complexity ','Fraction of SP3 hybridized Carbon atoms','Fraction of Chiral Carbon atoms','Number of Heavy atom','Number of Nitrogen atoms','Number of Sulphur atoms','Number of Heteroatoms','Number of Rotatable bonds','Number of Aliphatic Carbon cycles','Number of Aliphatic Heterocycles','Number of Aliphatic rings','Number of Aromatic Carbon cycles','Number of Aromatic Heterocycles','Number of Aromatic rings','Total number of rings','Number of Saturated Carbon cycles','Number of Saturated Heterocycles','Number of Saturated rings','Smallest set of smallest rings (SSSR)'])+'\n')

for Counter,mol in enumerate(suppl):

    if mol is None:
        print ('Error',Counter)
        continue
    
    else:
        iden = mol.GetProp('_Name')
        #print (Counter, iden)

        formula = Chem.rdMolDescriptors.CalcMolFormula(mol) # newly added 2 July 2019

        Molwt = Descriptors.MolWt(mol)

        ExactMolMwt = Descriptors.ExactMolWt(mol) # newly added

        LogP = Descriptors.MolLogP(mol)

        MR = Descriptors.MolMR(mol) # newly added

        TPSA = Descriptors.TPSA(mol)

        HDonors = Descriptors.NumHDonors(mol)

        HAccep = Descriptors.NumHAcceptors(mol)

        ChiralCenters = Chem.FindMolChiralCenters(Chem.MolFromSmiles(Chem.MolToSmiles(mol,isomericSmiles=False)),includeUnassigned=True) # new modification; can  include N or S or P chiral centers
        carbon_ChiralCenters = 0
        tmp_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol,isomericSmiles=False))
        for cc in ChiralCenters:
            if tmp_mol.GetAtomWithIdx(cc[0]).GetSymbol() == 'C':
                carbon_ChiralCenters = carbon_ChiralCenters + 1  # newly added; count only carbon chiral centers
            else:
                continue

        NumCarbon = counter(mol,Chem.MolFromSmarts("[#6]"))

        NumSp2C = counter(mol,Chem.MolFromSmarts("[$([cX3](:*):*),$([cX2+](:*):*),$([CX3]=*),$([CX2+]=*)]")) # newly added
        
        NumSp2C_test = counter(mol,Chem.MolFromSmarts("[#6^2]")) 

        NumSp3C = counter(mol,Chem.MolFromSmarts("[#6^3]"))

        if NumSp2C + NumSp3C == 0:
            shape_complexity = 0
        else:
            shape_complexity = round(float(NumSp3C) / (NumSp2C + NumSp3C),2)
        
        NumSpC = counter(mol,Chem.MolFromSmarts("[#6^1]"))

        if NumCarbon == 0:
            FCC = 0
        else:
            FCC = round(float(carbon_ChiralCenters)/float(NumCarbon),2)

        FSP3 = round(Descriptors.FractionCSP3(mol),2)

        HetAtoms = Descriptors.NumHeteroatoms(mol)

        NumNitrogen = counter(mol,Chem.MolFromSmarts("[#7]"))

        NumSulphur = counter(mol,Chem.MolFromSmarts("[#16]"))

        RotBonds_step1 = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol,Chem.rdMolDescriptors.NumRotatableBondsOptions.NonStrict) - Chem.rdMolDescriptors.CalcNumAmideBonds(mol) # All rotatble bonds minus amide bonds as required for Veber
        if RotBonds_step1 < 0:
            RotBonds = 0
        else:
            RotBonds = RotBonds_step1
        
        
        AliCarboCycl = Descriptors.NumAliphaticCarbocycles(mol)

        AliHetCycl = Descriptors.NumAliphaticHeterocycles(mol)

        AliRings = Descriptors.NumAliphaticRings(mol)

        AroCarboCycl = Descriptors.NumAromaticCarbocycles(mol)

        AroHetCycl = Descriptors.NumAromaticHeterocycles(mol)

        AroRings = Descriptors.NumAromaticRings(mol)

        TotalRings = Descriptors.RingCount(mol)

        SatCarboCycl = Descriptors.NumSaturatedCarbocycles(mol)

        SatHetCycl = Descriptors.NumSaturatedHeterocycles(mol)

        SatRings = Descriptors.NumSaturatedRings(mol)

        HeavyAtom = Descriptors.HeavyAtomCount(mol)

        SSSR = Chem.GetSSSR(mol)

        out = [str(i) for i in [iden,formula,Molwt,ExactMolMwt,LogP,MR,TPSA,HDonors,HAccep,carbon_ChiralCenters,NumCarbon,NumSpC,NumSp2C,NumSp2C_test,NumSp3C,shape_complexity,FSP3,FCC,HeavyAtom,NumNitrogen,NumSulphur,HetAtoms,RotBonds,AliCarboCycl,AliHetCycl,AliRings,AroCarboCycl,AroHetCycl,AroRings,TotalRings,SatCarboCycl,SatHetCycl,SatRings,SSSR]]

        f.write('\t'.join(out)+'\n')

f.close()


