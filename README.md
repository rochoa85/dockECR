# dockECR

**dockECR: open consensus docking and ranking protocol for virtual screening of small molecules**

**Authors:** Rodrigo Ochoa, Karen Palacio-Rodriguez, Camila M. Clemente and Natalia Adler

### Purpose

The goal of dockECR is to implement multiple docking programs to predict scores between a single target or multiple targets (merging and shrinking aproach) and a library of compounds. The molecules are ranked based on the ECR ranking metric using as input the scores obtained by each program. The protocol is open source and requires to install third-party docking software, as well as information of the targets and ligands structures in PDB format

### Third-party tools

The script is written in python3 with system calls to the Ubuntu OS. The docking stand-alone programs can be put after installation in the folder `auxiliar/software` **(alternatively, the code can be modified to call the programs without local paths, in case they have been installed in the system path.)**. For instructions about how to download and install each program follow these websites:

- AutoDock Vina: http://vina.scripps.edu/download.html - stand-alone tools: vina
- Smina: https://sourceforge.net/projects/smina/files/ - stand-alone tools: smina.static
- LeDock: http://www.lephar.com/download.htm - stand-alone tools: ledock_linux_x86
- rDock: https://sourceforge.net/projects/rdock/files/ - stand-alone tools: rbcavity, rbdock and sdsort. **To compile the program it is required gcc > 3 and < 6, plus additional requirements mentioned in the project website**
- MGL Tools: http://mgltools.scripps.edu/downloads - stand-alone tools: pythonsh and the scripts *prepare_ligand4.py* and *prepare_receptor4.py*
- RDKit: http://rdkit.org - version for python3. 

**NOTE:** To facilitate the compilation/installation of rDock and RDKit, it is recommended to create a conda (https://docs.conda.io/en/latest/) virtual environment. The packages can be installed by adding the corresponding channels after the environment is created:

```
conda config --add channels bioconda
conda config --add channels conda-forge
```
Then the packages are installed with:
```
conda install -c bioconda rdock
conda install -c conda-forge rdkit
```

In addition to these packages, it is required to install OpenBabel (https://sourceforge.net/projects/openbabel/), which can be obtained from the main OS package repository. **The version 2.3 is required (the program called babel instead of obabel).**

### How to run the consensus docking and ranking script

The workspace requires of the following folders:

```
auxiliar: Folder containing scripts and the stand-alone programs
target: Folder with the targets PDB files
ligands: Folder with the ligands PDB files
```

After running docking and ranking campaigns, the following additional folders will be created:

```
config: Folder that will store the configuration files created per program
results: Docking results per program used
ranks: Ranks per each target and docking program used, as well as the final ECR method.
```

An example of the script syntax is as follows::

`python dockECR.py [-h] -l LIST_LIGANDS -s LIST_SOFTWARE -m MODE -t LIST_TARGETS`
                                       
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
  -t LIST_TARGETS   File with the list of the PDB structures used as targets
 ```

The script has two main modes. One is for running docking with any of the four programs mentioned in the `list_software.txt` file, and the second is for ranking with same list of docking programs. In the following sections we show how to use the main python script. The examples provided in the code were run using the SARS-CoV-2 main protease with PDB id 6y2e and 5re4.

### Docking example

After preparing the `list_software.txt` file, and the `list_ligands.txt` file with the names of the PDB structures located in the `ligands` folder, the command can be run as:

```
python dockECR.py -m docking -l list_ligands.txt -s list_software.txt -t list_targets.txt
```
**NOTE: Examples of all the input files are included in the code.**

For docking, we require to prepare a configuration file with the box center coordinates and the box size dimensions. The following is an example of a box located in the active site of the MPro structure 6y2e provided in the `target` folder:

```
center_x: -11.6
center_y: 14.5
center_z: 68.3
size_x: 30
size_y: 30
size_z: 30
```
**NOTE: The config file should be named config_[TARGET].txt, where the [TARGET] is the name given in the list_targets.txt file.**

After running the protocol, the results with 10 poses per program are stored in the folder `results/[program]` as files with extensions depending on the docking program output.

### Ranking example

With the docking results available in the workspace, the same script can be called using the *ranking* mode to generate rankings per program from the pool of ligands used in the virtual screening. An example is provided here:

```
python dockECR.py -m ranking -l list_ligands.txt -s list_software.txt -t list_targets.txt
```

The `list_targets.txt` file can contain the names of multiple targets, in order to do the ranking individually, and using a merging and shrinking approach (see Ref. Palacio-Rodriguez et. al. 2019) to combine the results of multiple version of the same target. The ECR ranking is calculated in both cases.

The ranking results are stored in the `ranks` folder with a three-column format file, where the first column is the ranking position, the second column is the ligand id, and the third column the docking score:

```
1	mol3	-12.0649
2	mol2	-11.8435
3	mol1	-10.6486
```

With the rankings per software, it is possible to use the ECR method to generate a final consensus ranking. The results are sorted by the calculated ECR ranking in a file located in the `ranks` folder with the following format:

```
1 mol1       1.05880939
2 mol3       0.94259210
3 mol2       0.66164496
```

**NOTE: If multiple targets are included, an ECR ranking per target and for the combination of them is generated.**

This final list can be used to prioritize candidates from a large number of molecules in a virtual screening campaign, using the outputs of multiple docking programs at the same time.

An example of the docking and ranking results for the two included systems is available in the folder: `example_outputs`

### Support

In case the protocol is useful for other research projects and require some advice, please contact us to the email: rodrigo.ochoa@udea.edu.co. If this protocol is implemented, please also cite: 

- Palacio-Rodriguez, K., Lans, I., Cavasotto, C. N., & Cossio, P. (2019). Exponential consensus ranking improves the outcome in docking and receptor ensemble docking. Scientific reports, 9(1), 1-14.
