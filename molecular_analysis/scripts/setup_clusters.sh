#!/bin/bash

# This script it used to move all of the PQR files to the correct location
# and prep them for running the setup.sh
# No longer converts from pdb to pqr because of saving frame approach from cpptraj.
# Frequently this would be ran in the /dATP_multiscale_modeling/molecular_analysis/bd_adp/metastable_0 (ADP or dADP) 
# and adjusted for the metstable state. The cluster directory would then be:
# dATP_multiscale_modeling/molecular_analysis/bd_adp/metastable_0/aligned_structures

# Note: This script will prepare the simulations for a single metastable state with 15 structures. Run 3x for ADP metastable
# states and another 3x for dADP metastable states. 

# Proper usage: setup_clusters.sh -t [Numer of threads] -n [Number of Trajectories] -c [path to cluster directory]

# IMPORTANT: You must change the MOLECULAR_ROOT_DIR to match your system to the "dATP_multiscale_modeling/molecular_analysis"
# location on your machine. 

while getopts t:n:c:a: flag
do
    case "${flag}" in
        t) threads=${OPTARG};;
        n) sims=${OPTARG};;
        c) cluster_dir=${OPTARG};;
        a) apbs=${OPTARG};;
    esac
done

echo $cluster_dir

if (($# == 0))
then
    echo "No positional arguments specified"
    exit
fi

# Cluster directory is now pulled from command line argument 

MOLECULAR_ROOT_DIR="/crucial/dATP_multiscale_modeling/molecular_analysis"
STARTING_DIR=$(pwd)

# Molecule 0
mol0="myosin"
# Molecule 1
mol1="actin"
# Reaction Distance 
reaction_distance="3.5"

echo Creating subdirectories for each of the clusters...
for item in {1..15}
do
    # Create directories for each cluster structure
    mkdir cluster_$item

    # Copy myosin
    echo Copy myosin
    echo cp ${cluster_dir}/myosin_$item.pqr cluster_$item/myosin.pqr
    cp $cluster_dir/myosin_$item.pqr cluster_$item/myosin.pqr
    # Copy actin
    cp ${MOLECULAR_ROOT_DIR}/prep_template_documents/actin.pqr cluster_$item/.

    # Copy other input files 
    # No longer needed because of the apbs files calculated consistently for all the structures
    # cp ${MOLECULAR_ROOT_DIR}/prep_template_documents/pre_input.xml cluster_$item/.

    # Copy input.xml file 
    cp ${MOLECULAR_ROOT_DIR}/prep_template_documents/input.xml cluster_$item/.


    # Copy APBS input files. 
    cp ${MOLECULAR_ROOT_DIR}/prep_template_documents/myosin*.in cluster_$item/.
    cp ${MOLECULAR_ROOT_DIR}/prep_template_documents/actin*.in cluster_$item/.

    # Replace text with threads and sims
    sed -i "s/NUM_TRAJECTORIES/$sims/g" cluster_$item/input.xml
    sed -i "s/NUM_THREADS/$threads/g" cluster_$item/input.xml


    # Copy protein contacts file
    cp ${MOLECULAR_ROOT_DIR}/prep_template_documents/protein_protein_contacts.xml.bak cluster_$item/protein_protein_contacts.xml

    pqr2xml < cluster_$item/${mol0}.pqr > cluster_$item/${mol0}.xml
    pqr2xml < cluster_$item/${mol1}.pqr > cluster_$item/${mol1}.xml

    echo Creating reaction pairs in cluster_$item state 
    # Define reaction pairs
    make_rxn_pairs -mol0 cluster_$item/${mol0}.xml -mol1 cluster_$item/${mol1}.xml -ctypes cluster_$item/protein_protein_contacts.xml -dist ${reaction_distance} > cluster_$item/reaction_pairs.xml


done
echo Done creating files 

echo Running python to create a superset of reaction pairs. 
python ${MOLECULAR_ROOT_DIR}/scripts/combine_pairs.py

echo Creating one reaction file to use for all 
for item in {1..15}
do 
    make_rxn_file -pairs cluster_$item/new_reaction_pairs.xml -state_from before -state_to after -rxn association -mol0 ${mol0} ${mol0} -mol1 ${mol1} ${mol1} -distance 0.0 -nneeded 3 > cluster_$item/reactions.xml

done



if [ "$apbs" = "TRUE" ]; then
    echo "Run APBS set to true, running APBS for all the input files. Assuming actin is identical between clusters."   

    # Generate actin grid first 
    cd cluster_1
    apbs actin0.in 
    apbs actin1.in 
    cd $STARTING_DIR
    

    for item in {1..15}
    do 
        cd cluster_$item 
        apbs myosin0.in 
        apbs myosin1.in 
        echo Copying actin dx 
        cp ../cluster_1/actin*dx .
        cd $STARTING_DIR
    done 
else
    echo "The APBS variable is not TRUE."
    echo "Not running APBS" 
fi

echo "Copying the run_all.sh" script that can be used to run BD on all the directories sequentially..."
cp ${MOLECULAR_ROOT_DIR}/scripts/run_all.sh



echo Complete. 

