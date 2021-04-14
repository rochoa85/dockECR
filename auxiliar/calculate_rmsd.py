#!/usr/bin/python

import os
import sys
import math
from statistics import mean

# rdkit imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import AlignMol

def GetBestRMSD(probe,ref,refConfId=-1,probeConfId=-1,maps=None):
    def RMSD(probe,ref,amap):
        rmsd = 0.0
        # print(amap)
        atomNum = ref.GetNumAtoms() + 0.0
        for (pi,ri) in amap:
            posp = probe.pos[pi]
            posf = ref.pos[ri]
            rmsd += dist_2(posp,posf)
        rmsd = math.sqrt(rmsd/atomNum)
        return rmsd

    def dist_2(atoma_xyz, atomb_xyz):
        dis2 = 0.0
        for i, j  in zip(atoma_xyz,atomb_xyz):
            dis2 += (i -j)**2
        return dis2

    def orginXYZ(mol):
        mol_pos={}
        for i in range(0,mol.GetNumAtoms()):
            pos = mol.GetConformer().GetAtomPosition(i)
            mol_pos[i] = pos
        return mol_pos

    # When mapping the coordinate of probe will changed!!!
    ref.pos = orginXYZ(ref)
    probe.pos = orginXYZ(probe)

    if not maps:
        matches = ref.GetSubstructMatches(probe,uniquify=False)
        maps = [list(enumerate(match)) for match in matches]

        bestRMSD = 1000.0
        rmsd=100.0
        for amap in maps:
            rmsd = RMSD(probe,ref,amap)
        if rmsd<bestRMSD:
            bestRMSD = rmsd

    return bestRMSD

def calculate_mean_RMSD(folder,ligand,methods,sigma):
    pairs=[]
    rmsd_values=[]
    for m1 in methods:
        for m2 in methods:
            if m1!=m2:
                if (m1,m2) not in pairs and (m2,m1) not in pairs:
                    ref = Chem.MolFromPDBFile("{}/{}_{}.pdb".format(folder,ligand,m1))
                    probe = Chem.MolFromPDBFile("{}/{}_{}.pdb".format(folder,ligand,m2))
                    #itry:
                    rmsd = GetBestRMSD(probe,ref)
                    pairs.append((m1,m2))
                    rmsd_values.append(rmsd)
                    #except:
                    #    pairs.append((m1,m2))
                    #    pass
    # ECR rmsd
    prob_rmsd=0
    for val in rmsd_values:
        prob_rmsd+=math.exp(-1.0*val/sigma)/sigma

    return mean(rmsd_values),prob_rmsd
