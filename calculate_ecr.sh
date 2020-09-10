#!/bin/bash

#################################################################################################################
# Script for consensus docking merging and shrinking  
# Input files: 3 column files (rankings) with rank - name of molecule - score                 
#################################################################################################################

########################
# - Arguments
########################

function help(){
echo "-------------------------------------------"
echo "HELP SECTION"
echo "-------------------------------------------"
echo "./calculate_ecr.sh -t target -s list_software.txt -l list_ligands
where:
 	-t : id of the target used in the screening
    -s : file with the docking software used to create the rankings
    -l : file having the ligands used in the screening
    -h : show the help section"
}


flag1=false;flag2=false;flag3=false # Start the flag in false

while getopts "h:t:s:l:" option; do
	case $option in
		h) help
		   exit;;                
		t) target=${OPTARG}; flag1=true;;
        s) list_software=${OPTARG}; flag2=true;;
        l) list_ligands=${OPTARG}; flag3=true;;
		*) help
       	   exit;;
        esac
done

# Checking if any of the mandatory fields have not been assigned
if [ $flag1 == false ]; then
    echo "-t must be included when the program is called"
    help
    exit
fi
if [ $flag2 == false ]; then
    echo "-s must be included when the program is called"
    help
    exit
fi
if [ $flag3 == false ]; then
    echo "-l must be included when the program is called"
    help
    exit
fi

# Number of ligands
sig=`wc -l $list_ligands | awk '{print $1}'`

while read software; do
    
    sort -k2 ranks/rank_${software}_${target}.txt | awk '{print $2"\t"$1}' > s_${software}_${target}

done < $list_software

############ CONSENSUS ECR ###########

paste s_* | awk -v sig=$sig '{prob=0; for(x=2;x<=NF;x=x+2){prob+=exp(-$x/sig)/(sig)}; printf"%s %16.8f\n", $1,prob}' | sort -nrk2 | awk '{printf"%d %s %16.8f\n", NR,$1,$2}' > ranks/rank_ecr_${target}.txt
rm s_*