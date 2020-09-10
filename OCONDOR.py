#!/usr/bin/python

"""
OCONDOR: Open CONsensus DOcking and Ranking protocol for virtual screening of small molecules
Authors: Rodrigo Ochoa, Karen Palacio-Rodriguez, Natalia Adler, Camila M. Clemente and Pilar Cossio

Please also cite: 
Palacio-Rodriguez, K., Lans, I., Cavasotto, C. N., & Cossio, P. (2019).
Exponential consensus ranking improves the outcome in docking and receptor ensemble docking.
Scientific reports, 9(1), 1-14.

Third-party tools required (put them in the auxiliar/software folder):

AutoDock Vina: http://vina.scripps.edu/download.html - stand-alone tools: vina
Smina: https://sourceforge.net/projects/smina/files/ - stand-alone tools: smina.static
LeDock: http://www.lephar.com/download.htm - stand-alone tools: ledock_linux_x86
rDock: https://sourceforge.net/projects/rdock/files/ - stand-alone tools: rbcavity, rbdock and sdsort
MGL Tools: http://mgltools.scripps.edu/downloads - stand-alone tools: pythonsh

Additional tools:
Open Babel: https://sourceforge.net/projects/openbabel/
"""

########################################################################################
# Authorship
########################################################################################

__credits__ = ["Rodrigo Ochoa","Karen Palacio-Rodriguez","Natalia Adler","Camila Clemente","Pilar Cossio"]
__license__ = "MIT"
__version__ = "1.0"

########################################################################################
# Modules to import
########################################################################################

import argparse
import subprocess
import operator
import os
import yaml
import multiprocessing
from auxiliar import calculate_rmsd

########################################################################################
# Preparation functions
########################################################################################

def prepare_pdbqt(target,ligand):
    os.system("./auxiliar/software/pythonsh auxiliar/prepare_receptor4.py -r target/{tar}.pdb -o {tar}.pdbqt -U waters".format(tar=target))
    os.system("./auxiliar/software/pythonsh auxiliar/prepare_ligand4.py -l ligands/{lig}.pdb -U '' -B -o {lig}.pdbqt".format(lig=ligand))

########################################################################################

def prepare_rdock_cavity(target,ligand,center_x,center_y,center_z,size_x,size_y,size_z):
    # Create config 
    prepare_pdbqt(target,ligand)
    config=open("config/config_vina_{}.txt".format(target),"w")
    config.write("center_x={}\n".format(center_x))
    config.write("center_y={}\n".format(center_y))
    config.write("center_z={}\n".format(center_z))
    config.write("size_x={}\n".format(size_x))
    config.write("size_y={}\n".format(size_y))
    config.write("size_z={}\n".format(size_z))
    config.write("cpu=1\n")
    config.write("exhaustiveness=1\n")
    config.write("num_modes = 1\n")
    config.close()
    
    os.system("./auxiliar/software/vina --receptor {}.pdbqt --ligand {}.pdbqt --log score.log --out out.pdbqt --config config/config_vina_{}.txt".format(target,ligand,target))
    # Split and get the first model
    os.system("csplit out.pdbqt /MODEL/ {*}; mv xx01 ligand.pdbqt")
    os.system("babel -ipdbqt ligand.pdbqt -osdf ligand.sd")
    os.system("babel -ipdb target/{}.pdb -omol2 receptor.mol2".format(target))
    os.system("./auxiliar/software/rbcavity -was -d -r config/rdock.prm")
    os.system("rm xx* score.log *.pdbqt")

########################################################################################
# Docking functions
########################################################################################

########################################################################################
            
def run_vina(target,ligand,center_x,center_y,center_z,size_x,size_y,size_z):
    
    # Run vina and store results
    print("Docking with vina between target {} and ligand {} ...".format(target,ligand))
    os.system("./auxiliar/software/vina --receptor {}.pdbqt --ligand {}.pdbqt --log score_{}.log --out out_{}.pdbqt --config config/config_vina_{}.txt".format(target,ligand,ligand,ligand,target))
    os.system("mv out_{}.pdbqt results/vina/{}_{}_out.pdbqt".format(ligand,target,ligand))

########################################################################################

def run_smina(target,ligand,center_x,center_y,center_z,size_x,size_y,size_z):
    
    # Run smina and store results
    print("Docking with smina between target {} and ligand {} ...".format(target,ligand))
    os.system("./auxiliar/software/smina.static --receptor {}.pdbqt --ligand {}.pdbqt --log score_{}.log --out out_{}.pdbqt --config config/config_smina_{}.txt".format(target,ligand,ligand,ligand,target))
    os.system("mv out_{}.pdbqt results/smina/{}_{}_out.pdbqt".format(ligand,target,ligand))

########################################################################################

def run_rdock(target,ligand):
    
    # Convert ligand to sdf file format
    os.system("babel -ipdb ligands/{}.pdb -osdf {}.sd".format(ligand,ligand))
    
    # Run rdock and store results
    print("Docking with rdock between target {} and ligand {} ...".format(target,ligand))
    os.system("./auxiliar/software/rbdock -i {}.sd -o {}_{} -r rdock.prm -p dock.prm -n 10".format(ligand,target,ligand,target))
    os.system("./auxiliar/software/sdsort -n -f'SCORE.INTER' {}_{}.sd > {}_{}.sorted".format(target,ligand,target,ligand))
    os.system("mv {}_{}.sorted results/rdock/{}_{}_out.sd".format(target,ligand,target,ligand))
    #os.system("rm *.sd *.mol2 rdock* dock.prm")

########################################################################################

def run_ledock(target,ligand,xmin,xmax,ymin,ymax,zmin,zmax):
    
    os.system("babel -ipdb ligands/{}.pdb -omol2 {}.mol2".format(ligand,ligand))
    os.system("cp target/{}.pdb {}_{}.pdb".format(target,target,ligand))
    ligands_file=open("mol-list_{}".format(ligand),"w")
    ligands_file.write("{}.mol2".format(ligand))
    ligands_file.close()
    
    # Write config file
    config=open("config_ledock_{}_{}.txt".format(target,ligand),"w")
    config.write("Receptor\n")
    config.write("{}_{}.pdb\n".format(target,ligand))                       
    config.write("RMSD\n")                           
    config.write("1.0\n")                            
    config.write("Binding pocket\n")
    config.write("{} {}\n".format(xmin,xmax))
    config.write("{} {}\n".format(ymin,ymax))
    config.write("{} {}\n".format(zmin,zmax))
    config.write("Number of binding poses\n")
    config.write("10\n")
    config.write("Ligands list\n")
    config.write("mol-list_{}\n".format(ligand))                    
    config.write("END")
    config.close()
    
    # Run ledock and store results
    print("Docking with ledock between target {} and ligand {} ...".format(target,ligand))
    os.system("auxiliar/software/ledock_linux_x86 config_ledock_{}_{}.txt".format(target,ligand))
    os.system("mv {}.dok results/ledock/{}_{}_out.dok".format(ligand,target,ligand))

########################################################################################
# Ranking functions
########################################################################################

def generate_complex(target,ligand,software):
    # Generate complexes with the receptor and the docked ligands
    os.system("grep HETATM ligand.pdb | sed 's/LIG  /LIG B/g' > ligand_fixed.pdb")
    bash="wc -l ligand_fixed.pdb | awk '{print $1}'"
    total_lines=int(subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8"))
    
    # Fix the atom numeration
    for i in range(1,total_lines+1):
        bash="head -n {} ligand_fixed.pdb | tail -n 1 | awk '{{print $3}}'".format(i)
        atom=subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
        if i<10:
            os.system("head -n {} ligand_fixed.pdb | tail -n 1 | sed 's/{}   /{}{}  /g' >> ligand_fixed_ref.pdb".format(i,atom,atom,i))
        else:
            os.system("head -n {} ligand_fixed.pdb | tail -n 1 | sed 's/{}   /{}{} /g' >> ligand_fixed_ref.pdb".format(i,atom,atom,i))
    
    # Store the files
    os.system("mv ligand_fixed_ref.pdb ligand_fixed.pdb")
    os.system("cp ligand_fixed.pdb temp_ranking/{}_{}.pdb".format(ligand,software))
    os.system("cp target/{}.pdb .".format(target))
    os.system("cat {}.pdb ligand_fixed.pdb | sed 's/END/TER/g' > complex_{}_{}_{}.pdb".format(target,target,ligand,software))
    os.system("echo 'TER' >> complex_{}_{}_{}.pdb".format(target,ligand,software))

########################################################################################

def score_vina(target,ligand):
    # Split and get the first model
    os.system("csplit -s results/vina/{}_{}_out.pdbqt /MODEL/ {{*}}; mv xx01 ligand.pdbqt".format(target,ligand))
    bash="grep VINA ligand.pdbqt | awk '{print $4}'"
    score_v = float(subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8"))
    
    # Generate complex
    os.system("babel -ipdbqt ligand.pdbqt -opdb ligand.pdb")
    software="vina"
    generate_complex(target,ligand,software)
    os.system("mv complex_{}_{}_vina.pdb complexes/vina".format(target,ligand))  
    os.system("rm xx* *.pdb *.pdbqt")
    
    return score_v

########################################################################################

def score_smina(target,ligand):
    # Split and get the first model
    os.system("csplit -s results/smina/{}_{}_out.pdbqt /MODEL/ {{*}}; mv xx01 ligand.pdbqt".format(target,ligand))
    bash="grep minimizedAffinity ligand.pdbqt | awk '{print $3}'"
    score_s = float(subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8"))
    
    # Generate complex
    os.system("babel -ipdbqt ligand.pdbqt -opdb ligand.pdb")
    software="smina"
    generate_complex(target,ligand,software)
    os.system("mv complex_{}_{}_smina.pdb complexes/smina".format(target,ligand))    
    os.system("rm xx* *.pdb *.pdbqt")
    
    return score_s

########################################################################################

def score_ledock(target,ligand):
    # Split and get the first model
    os.system("csplit -s results/ledock/{}_{}_out.dok /REMARK/ {{*}}; mv xx02 ligand.pdb".format(target,ligand))
    bash="grep Cluster ligand.pdb | awk '{print $8}'"
    score_l = float(subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8"))
    
    # Generate complex
    os.system("grep -v REMARK ligand.pdb | sed 's/ATOM  /HETATM/g' > ligand_mod.pdb")
    os.system("mv ligand_mod.pdb ligand.pdb")
    software="ledock"
    generate_complex(target,ligand,software)
    os.system("mv complex_{}_{}_ledock.pdb complexes/ledock".format(target,ligand))    
    os.system("rm xx* *.pdb")
    
    return score_l

########################################################################################

def score_rdock(target,ligand):
 
    # Split and get the first model
    os.system("babel -isdf results/rdock/{}_{}_out.sd -opdb {}_{}_out.pdb".format(target,ligand,target,ligand))
    os.system("csplit -s {}_{}_out.pdb /MODEL/ {{*}}; mv xx01 ligand.pdb".format(target,ligand))
    
    os.system("csplit -s results/rdock/{}_{}_out.sd /SCORE/ {{*}}; mv xx01 score.txt".format(target,ligand))
    bash="head -n 2 score.txt | tail -n 1"
    score_r = float(subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8"))
    
    # Generate complex
    software="rdock"
    generate_complex(target,ligand,software)
    os.system("mv complex_{}_{}_rdock.pdb complexes/rdock".format(target,ligand))    
    os.system("rm xx* *.pdb score.txt")
    
    return score_r

########################################################################################
# Main execution
########################################################################################

if __name__ == '__main__':
    
    # Script arguments
    parser = argparse.ArgumentParser(description='OCONDOR: Open CONsensus DOcking and Ranking protocol')
    parser.add_argument('-l', dest='list_ligands', action='store',required=True,
                        help='File with the list of ligands names available in the ligand folder')
    parser.add_argument('-s', dest='list_software', action='store',required=True,
                        help='File with the list of software that will be included in the consensus')
    parser.add_argument('-m', dest='mode', action='store',required=True,
                        help='Mode the script will be run. Options: (i) docking, (ii) ranking')
    parser.add_argument('-t', dest='target', action='store',required=True,
                        help='Name of the PDB structure used as target')
    parser.add_argument('-c', dest='config_file', type=argparse.FileType(mode='r'),
                        help='File containing the center coordinates and size of the box. Note: Only for docking mode') 
    
    #####################################################################################
    # Assignment of parameters
    #####################################################################################
    args = parser.parse_args()
    
    # Map the main parameters
    list_software=[x.strip() for x in open(args.list_software)] # List of software
    list_ligands=[x.strip() for x in open(args.list_ligands)] # List of ligands
    target=args.target # Target structure
    mode=args.mode # Mode: docking or ranking
    
    # Check the coordinates
    if args.config_file:
        data = yaml.load(args.config_file)
        delattr(args, 'config_file')
        arg_dict = args.__dict__
        for key, value in data.items():
            if isinstance(value, list):
                for v in value:
                    arg_dict[key].extend(v)
            else:
                arg_dict[key] = value
    else:
        if mode=="docking":
            print("A config file is necessary to run the protocol. Exiting ...")
            exit()
    
    # Check the arguments
    try:
        if args.center_x: center_x=float(args.center_x)
    except:
        if mode=="docking":
            print("The parameter 'center_x' is required for the analysis. Exiting ...")
            exit()
    try:
        if args.center_y: center_y=float(args.center_y)
    except:
        if mode=="docking":
            print("The parameter 'center_y' is required for the analysis. Exiting ...")
            exit()
    try:
        if args.center_z: center_z=float(args.center_z)
    except:
        if mode=="docking":
            print("The parameter 'center_z' is required for the analysis. Exiting ...")
            exit()
    try:
        if args.size_x: size_x=int(args.size_x)
    except:
        if mode=="docking":
            print("The parameter 'size_x' is required for the analysis. Exiting ...")
            exit()
    try:
        if args.size_y: size_y=int(args.size_y)
    except:
        if mode=="docking":
            print("The parameter 'size_y' is required for the analysis. Exiting ...")
            exit()
    try:
        if args.size_z: size_z=int(args.size_z)
    except:
        if mode=="docking":
            print("The parameter 'size_z' is required for the analysis. Exiting ...")
            exit()

    ####################################################################################
    # Run docking protocols
    ####################################################################################
    
    if mode=="docking":
        
        # Size parameter for ledock
        xmin=center_x-(size_x/2); xmax=center_x+(size_x/2)
        ymin=center_y-(size_y/2); ymax=center_y+(size_y/2)
        zmin=center_z-(size_z/2); zmax=center_z+(size_z/2)
        
        # Iterate over the software list
        for software in list_software:
            if software == "vina":
                
                # Create config file
                config=open("config/config_vina_{}.txt".format(target),"w")
                config.write("center_x={}\n".format(center_x))
                config.write("center_y={}\n".format(center_y))
                config.write("center_z={}\n".format(center_z))
                config.write("size_x={}\n".format(size_x))
                config.write("size_y={}\n".format(size_y))
                config.write("size_z={}\n".format(size_z))
                config.write("cpu=1\n")
                config.write("exhaustiveness=1\n")
                config.write("num_modes = 10\n")
                config.close()
                
                pool=multiprocessing.Pool()
                for ligand in list_ligands:
                    prepare_pdbqt(target,ligand)
                    # Create parallel jobs
                    pool.apply_async(run_vina, args=(target,ligand,center_x,center_y,center_z,size_x,size_y,size_z,))
                    
                pool.close()
                pool.join()
                # Delete all temporal files
                os.system("rm *.pdbqt score_*.log")
            
            if software == "smina":
                # Create config file
                config=open("config/config_smina_{}.txt".format(target),"w")
                config.write("center_x={}\n".format(center_x))
                config.write("center_y={}\n".format(center_y))
                config.write("center_z={}\n".format(center_z))
                config.write("size_x={}\n".format(size_x))
                config.write("size_y={}\n".format(size_y))
                config.write("size_z={}\n".format(size_z))
                config.write("cpu = 1\n")
                config.write("exhaustiveness=1\n")
                config.write("num_modes = 10\n")
                config.write("scoring = vinardo\n")
                config.close()
                
                pool=multiprocessing.Pool()
                for ligand in list_ligands:
                    prepare_pdbqt(target,ligand)
                    # Create parallel jobs
                    pool.apply_async(run_smina, args=(target,ligand,center_x,center_y,center_z,size_x,size_y,size_z,))
                
                pool.close()
                pool.join()
                # Delete all temporal files
                os.system("rm *.pdbqt score_*.log")
                    
            if software == "ledock":
                pool=multiprocessing.Pool()
                for ligand in list_ligands:
                    # Create parallel jobs
                    pool.apply_async(run_ledock, args=(target,ligand,xmin,xmax,ymin,ymax,zmin,zmax,))
                
                pool.close()
                pool.join()
                # Delete all temporal files
                os.system("rm *.pdb *.mol2 mol-list_* config_ledock_*")
                    
            if software == "rdock":
                pool=multiprocessing.Pool()
                for i,ligand in enumerate(list_ligands):
                    if i==0:
                        prepare_rdock_cavity(target,ligand,center_x,center_y,center_z,size_x,size_y,size_z)
                        os.system("babel -ipdb target/{}.pdb -omol2 receptor.mol2".format(target))
                        os.system("cp config/rdock* .")
                        os.system("cp auxiliar/dock.prm .")
                        run_rdock(target,ligand)
                        os.system("cp config/rdock* .")
                        os.system("cp auxiliar/dock.prm .")
                    else:
                        # Create parallel jobs
                        pool.apply_async(run_rdock, args=(target,ligand,))
                    
                pool.close()
                pool.join()
                # Delete all temporal files
                os.system("rm *.sd *.mol2 rdock* dock.prm")
    
    ####################################################################################
    # Run ranking protocols
    ####################################################################################
    
    if mode=="ranking":
        
        os.system("mkdir temp_ranking")
        
        for software in list_software:
            if software == "vina":
                ranked_ligands_vina={}
                for ligand in list_ligands:
                    score_v=score_vina(target,ligand)
                    ranked_ligands_vina[ligand]=score_v
                
                print(ranked_ligands_vina)
                sorted_x = sorted(ranked_ligands_vina.items(), key=operator.itemgetter(1))
                
                rank_file_vina=open("ranks/rank_vina_{}.txt".format(target),"w")
                for j,element in enumerate(sorted_x):
                    rank_file_vina.write("{}\t{}\t{}\n".format(j+1,element[0],element[1]))
                rank_file_vina.close()
            
            if software == "smina":
                ranked_ligands_smina={}
                for ligand in list_ligands:
                    score_s=score_smina(target,ligand)
                    ranked_ligands_smina[ligand]=score_s
                
                print(ranked_ligands_smina)
                sorted_x = sorted(ranked_ligands_smina.items(), key=operator.itemgetter(1))
                
                rank_file_smina=open("ranks/rank_smina_{}.txt".format(target),"w")
                for j,element in enumerate(sorted_x):
                    rank_file_smina.write("{}\t{}\t{}\n".format(j+1,element[0],element[1]))
                rank_file_smina.close()
            
            if software == "ledock":
                ranked_ligands_ledock={}
                for ligand in list_ligands:
                    score_l=score_ledock(target,ligand)
                    ranked_ligands_ledock[ligand]=score_l
                
                print(ranked_ligands_ledock)
                sorted_x = sorted(ranked_ligands_ledock.items(), key=operator.itemgetter(1))
                
                rank_file_ledock=open("ranks/rank_ledock_{}.txt".format(target),"w")
                for j,element in enumerate(sorted_x):
                    rank_file_ledock.write("{}\t{}\t{}\n".format(j+1,element[0],element[1]))
                rank_file_ledock.close()
            
            if software == "rdock":
                ranked_ligands_rdock={}
                for ligand in list_ligands:
                    score_r=score_rdock(target,ligand)
                    ranked_ligands_rdock[ligand]=score_r
                
                print(ranked_ligands_rdock)
                sorted_x = sorted(ranked_ligands_rdock.items(), key=operator.itemgetter(1))
                
                rank_file_rdock=open("ranks/rank_rdock_{}.txt".format(target),"w")
                for j,element in enumerate(sorted_x):
                    rank_file_rdock.write("{}\t{}\t{}\n".format(j+1,element[0],element[1]))
                rank_file_rdock.close()
        
        # Calculate full RMSD
        ranked_ligands_rmsd={}
        for ligand in list_ligands:
            rmsd_m=calculate_rmsd.calculate_mean_RMSD("temp_ranking",ligand,list_software)
            ranked_ligands_rmsd[ligand]=rmsd_m
        
        print(ranked_ligands_rmsd)
        sorted_x = sorted(ranked_ligands_rmsd.items(), key=operator.itemgetter(1))
        
        rank_file_rmsd=open("ranks/rank_rmsd_{}.txt".format(target),"w")
        for j,element in enumerate(sorted_x):
            rank_file_rmsd.write("{}\t{}\t{}\n".format(j+1,element[0],element[1]))
        rank_file_rmsd.close()
        
        os.system("rm -r temp_ranking")
            
        # Run ECR ranking
        os.system("cp list_software.txt list_metrics.txt; echo rmsd >> list_metrics.txt")
        os.system("./calculate_ecr.sh -t {} -s list_metrics.txt -l list_ligands.txt".format(target))
        os.system("rm list_metrics.txt")
        print("ECR ranking completed")
