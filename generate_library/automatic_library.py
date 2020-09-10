#!/usr/bin/python

"""
Open source consensus docking and ranking protocol: application with SARS-CoV-2 main protease
Authors: Rodrigo Ochoa, Natalia Adler, Karen Palacio-Rodriguez and Camila M. Clemente


Third-party tools required (put them in the auxiliar/software folder):

RDKit: https://github.com/rdkit/rdkit/releases - Ubuntu package: python-biopython
Crem: https://github.com/DrrDom/crem
Open Babel: https://sourceforge.net/projects/openbabel/ - Ubuntu package: openbabel
"""

import os
from rdkit import Chem
from rdkit.Chem import AllChem
from crem.crem import mutate_mol, grow_mol, link_mols

# Read the fragments
fragments=[x.strip() for x in open("non-covalent_fragments.txt")]

report=open("new_molecules.txt","w")
molecules=[]

for f1 in fragments:
    for f2 in fragments:
        if f1!=f2:
            m = Chem.MolFromSmiles(f1) 
            m2 = Chem.MolFromSmiles(f2)
            # The db file can be downloaded from: http://www.qsar4u.com/files/cremdb/replacements02_sc2.5.db.gz
            mols = list(link_mols(m, m2, db_name='replacements_sc2.5.db'))
            if mols: 
                for m in mols:
                    if m not in molecules:
                        molecules.append(m)
                        report.write("{}\t{}\t{}\n".format(f1,f2,m))

report.close()
print(len(molecules))

# Read the smiles
total_smiles = [x.strip() for x in open("full_molecules.txt")]
smiles=[]

for s in total_smiles:
    info=s.split()
    smiles.append(info[2])

# Iterate over the smules and generate the pdb file
ids=[]
for i,smi in enumerate(smiles):
    mol=Chem.MolFromSmiles(smi)
    molfile=open('molecule_{}.mol'.format(i),'w')
    ids.append("molecule_{}".format(i))
    mol2=Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol2)
    molfile.write(Chem.MolToMolBlock(mol2))
    molfile.close()
    os.system("babel -imol molecule_{}.mol -opdb molecule_{}.pdb".format(i,i))
    os.system("rm molecule_{}.mol".format(i))
