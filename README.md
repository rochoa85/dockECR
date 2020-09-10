# OCONDOR

**OCONDOR: Open CONsensus DOcking and Ranking protocol for virtual screening of small molecules**

**Authors:** Rodrigo Ochoa, Karen Palacio-Rodriguez, Natalia Adler, Camila M. Clemente and Pilar Cossio

### Purpose

The goal of OCONDOR is to implement multiple docking programs to predict scores between a single target and a library of compounds. The molecules are ranked based on the ECR ranking metric using as input the scores obtained by each program. The protocol is open source and requires to install third-party docking software and the target and ligand structures in PDB format

### Third-party tools

The script is written in python3 with system calls to the Ubuntu OS. The docking stand-alone programs can be put after installation in the folder `auxiliar/software`. For instructions about how to download and install each program follow these websites:

- AutoDock Vina: http://vina.scripps.edu/download.html - stand-alone tools: vina
- Smina: https://sourceforge.net/projects/smina/files/ - stand-alone tools: smina.static
- LeDock: http://www.lephar.com/download.htm - stand-alone tools: ledock_linux_x86
- rDock: https://sourceforge.net/projects/rdock/files/ - stand-alone tools: rbcavity, rbdock and sdsort
- MGL Tools: http://mgltools.scripps.edu/downloads - stand-alone tools: pythonsh

In addition to these packages, it is required to install OpenBabel (https://sourceforge.net/projects/openbabel/), which can be obtained from the main OS package repository.

### How to run the consensus docking and ranking script

The workspace requires of the following folders:

```
auxiliar: Folder containing scripts and the stand-alone programs
target: Folder with the targets PDB files
ligands: Folder with the ligands PDB files
config: Folder that will store the configuration files created per program
results: Docking results per program used
ranks: Ranks per each target and docking program used
complexes: Final PDB files with the receptor and the docked best ligand in complex
```

An example of the script syntax is as follows::

`python OCONDOR.py [-h] -l LIST_LIGANDS -s LIST_SOFTWARE -m MODE -t
                            TARGET [-c CONFIG_FILE]`
                                       
where the arguments are:

```
arguments:
  -h, --help        show this help message and exit
  -l LIST_LIGANDS   File with the list of ligands names available in the
                    ligand folder
  -s LIST_SOFTWARE  File with the list of software that will be included in
                    the consensus
  -m MODE           Mode the script will be run. Options: (i) docking, (ii)
                    ranking
  -t TARGET         Name of the PDB structure used as target
  -c CONFIG_FILE    File containing the center coordinates and size of the box. Note: This is required only for docking mode
 ```

The script has two main modes. One is for running docking with any of the four programs mentioned in the `list_software.txt` file, and the second is for ranking with same list of docking programs. In the following sections we show how to use the main python script, and an additional bash script for ECR ranking. The examples provided in the code were run using the SARS-CoV-2 main protease with PDB id 6y2e.

### Docking example

After preparing the "list of software" file, and the "list of ligands" file with the PDB structures located in the `ligands` folder, the command can be run as:

```
python OCONDOR.py -m docking -l list_ligands.txt -s list_software.txt -t 6y2e -c config_docking.txt
```

For docking, we require to prepare a configuration file with the box center coordinates and the box size dimensions. The following is an example of a box located on the active site of the MPro structure 6y2e provided in the `target` folder:

```
center_x: -11.6
center_y: 14.5
center_z: 68.3
size_x: 30
size_y: 30
size_z: 30
```
After running the protocol, the results with 10 poses per program are stored in the folder `results/[program]` as files with extensions depending on the docking program output.

### Ranking example

With the docking results available in the workspace, the same script can be called using the *ranking* mode to generate rankings per program from the pool of ligands used in the virtual screening. An example using the same receptor is provided here:

```
python OCONDOR.py -m ranking -l list_ligands.txt -s list_software.txt -t 6y2e
```

The ranking results are stored in the `ranks` folder with a three-column format file, where the first column is the ranking position, the second column is the ligand id, and the third column the docking score:

```
1	mol3	-12.0649
2	mol2	-11.8435
3	mol1	-10.6486
```

With the rankings per software, it is possible to use the ECR method to generate a final consensus ranking as follows.

### ECR script

For the ECR calculation, a bash script was prepared, which can be called for this system as:

```
./calculate_ecr.sh -t 6y2e -s list_software.txt -l list_ligands.txt
```

The inputs are the target id, the file with the list of software and the file with the list of ligands used on the previous script. The results are finally sorted by the calculated ECR ranking in a file located in the `ranks` folder with the following format:

```
1 mol1       1.05880939
2 mol3       0.94259210
3 mol2       0.66164496
```

This final list can be used to prioritize candidates from a large number of molecules in a virtual screening campaign, using the outputs of multiple docking programs at the same time.

### Additional: automatic library generation

An additional script is provided to generate linked fragments from the list of MPro fragments reported in the COVID Moonshot initiative (https://covid.postera.ai/covid). The script makes use of the RDKit (https://www.rdkit.org/) and the CReM packages (https://github.com/DrrDom/crem) to link fragments with plaussible chemical modifications available in the ChEMBL database (https://www.ebi.ac.uk/chembl/). The generated SMILES are used to predit conformers in PDB format that can be included as input in the consensus docking and ranking open source protocol.

### Support

In case the protocol is useful for other research projects and require some advice, please contact us to the email: rodrigo.ochoa@udea.edu.co. If this protocol is implemented, please also cite: 

- Palacio-Rodriguez, K., Lans, I., Cavasotto, C. N., & Cossio, P. (2019). Exponential consensus ranking improves the outcome in docking and receptor ensemble docking. Scientific reports, 9(1), 1-14.
