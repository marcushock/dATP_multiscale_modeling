# Molecular Analysis  
This directory contains the scripts and notebooks that are used to carry out the molecular dynamics analysis, Markov State model construction and analysis, as well as the setup and analysis of the Brownian dynamics simulations. 

This analysis primarily uses the python package PyEmma for the MSM analysis, as well as other common python packages. The `MSM_environment.yml` file can be used to create a new conda environment to match that which was used in this approach. 

If you would like to follow the entire workflow, you will need to have the trajectories downloaded as well, which are not part of this github repository. Then the paths in some of the notebooks and scripts will need to be adjusted to point to the correct directory location. 

## Directory organization: 
The following directories are organized as follows: 
- `msm_datafiles` contains binary files that help to follow the previously used MSM analysis and load files more quickly. 
- `prep_template_documents` contains files that are used in the setup of the Brownian dynamics (BD) simulations after selecting structures from the MSM analysis. 
- The `scripts` directory contains python and bash scripts that are useful in preparing files for the BD simulations. 
- `bd_adp` and `bd_dadp` contains the template files that have been set up for Brownian dynamics simulations, as well out the output binding curves that were estimated by the simulations. This data is necessary for the `BD_results_analysis.ipynb` notebook. 

## General Workflow: 
Prior to running, make sure the [Browndye2](https://browndye.ucsd.edu/) is installed. The Browndye version 2.0: Version of 6 Oct 2023 was used for these simulations, but newer versions should be fine as well. Make sure that once Browndye is installed, the `browndye2/bin` is in the system path. 

As a summary, this approach: A) Creates a Markov Model for ADP.Pi bound myosin from and MD simulation, and similarly a MSM of dADP.Pi bound myosin. B) Samples structures from the metastbale states of the ADP and dADP MSM. C) Uses those sampled structures to carry out Brownian Dynamics simulations to estiamte the on-rate of prepowerstroke myosin to actin. 
1. The `Myosin_MSM_Analysis.ipynb` file is first used the carry out featurization, tICA dimmentionality reduction, k-means clustering, and PCCA+ clustering and finally provide a sampling of conformations from the identified metastable states. More details inside of notebook. 
2. Run cpptraj with the commands outputted from the `Myosin_MSM_Analysis.ipynb` notebook. Alternatively, the `bd_adp/adp_msm_frames` directory has the files already selected, as well as the script used to extract the frames. 
3. Rename the PDB files as necessary for the organizational structure. The `scripts` directory has some bash scripts that can be useful for changing the name. The final directory per metastable state should have a set of aligned structures that are named `myosin_N.pqr` where n is an integer from 1 to 15. 
4. Set up the individual cluster simulations using the `scripts/setup_clusters.sh` script. When in the desired directory (i.e. `metastable_1`) `bash setup_cluster.sh -t 24 -n 50000 -c aligned_structures -a TRUE`. (note, must add scripts to path, or specify the full path to the script location)
    - -t is the number of threads to parallelize the BD simulation 
    - -n is the number of trajectories to run (example 50000)
    - -a (optional) TRUE or FALSE to specify if you want the APBS electrostatic maps to be calculated as you set up the subdirectories. If empty or FALSE, in each subdirectory, `apbs actin0.in`, `apbs myosin0.in` etc will need to be ran. 
    - -c The directory from which to pull the structures. 
    - The existing github has the output of this command with -a set to False, therefore the APBS grids still need to be created. The outputs are in `bd_(d)adp/metastable_*` and have created directories `cluster_1...cluster_15`. 
5. Run the `run_all.sh` script within each of the `bd_(d)adp/metastable_*` directories which will run the `bd_top` setup program followed by `nam_simulation` program from Browndye. This analysis by default does not use a specified reaction distance (see `reactions.xml` file) and therefore will output a distance file. 
6. Run the `compute_all_rates.sh` to generate a binding curve which estimates the binding rate against a specified binding reaction criteria distance. This script will search for all instances of `output_distances.xml` in any subdirectories, and then run the BD program `rates_of_distances`. 
7. Carry out final analysis of the BD simulations using the `BD_results_analysis.ipynb` to estimate the overall binding rates based on the sampled structures from the MSMs. 