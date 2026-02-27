# dATP_multiscale_modeling
Contains code for the manuscript "Multiscale modeling shows how 2'-deoxy-ATP rescues ventricular function in heart failure"


The three main scripts for running the simulations from the paper are plot_script_myocyte_new.m, plot_script_ventricular_new.m, and plot_script_filament. Running these scripts will run all simulations (except for spatially explicit sarcomere simulations) and generate all of the main and supplemental figures from the paper.
## MOLECULAR 
The `molecular_analysis` directory contains scripts and notebooks for performing molecular dynamics analysis, constructing and analyzing Markov State Models (MSM), and setting up and analyzing Brownian dynamics (BD) simulations to estimate actomyosin interaction, and the changes due to dATP or ATP.  A much more detailed README is within the directory. The MD trajectories are available at the UC San Diego Library Digital Collections, dataset linked at https://doi.org/10.6075/J00G3KCC.

## FILAMENT (Spatially explicit sarcomere model)
Note: This code is set up to run on an NVIDIA GeForce RTX 3080 Ti with CUDA 11.7; you may need to reconfigure the code if you have a different GPU/version
*Refer to McCabe et al. 2020 for more details
- Makefile: If you change anything in any of the source files, you will need to run "make clean" and then "make all"

### Building with Docker (Recommended for reproducibility)
A Dockerfile is provided with CUDA 12.2 and Boost dependencies pre-configured.

**Build the Docker image:**
```bash
docker build -t my-cuda-boost:12.2 .
```

**Run compilation inside Docker:**
```bash
docker run --gpus all -v $(pwd):/workspace -it my-cuda-boost:12.2 bash
cd /workspace && make clean && make all
```

**Note:** Docker GPU support requires:
- Linux or Mac with Docker installed
- NVIDIA GPU with Docker GPU support
- Returns exit code 127 on Windows (GPU passthrough not supported)

### Local Build Requirements (without Docker)
- cuda compiler nvcc, version 11.7 or 12.2+
- boost library (libboost-system-dev, libboost-thread-dev)
- cuda tool kit (developer level)
- *Make* utility

### Running Simulations
- sumatra # Version 0.8.1 for carrying out different experimental runs. Example:
     `smt run expData/mohran_force_pCa_normal_normalized.csv parameters/parameter_input.csv`
    - Note that the `-r` flag can be used to name the reason for the "experiment" 
    - Also note that the `smt web` tool allows for a useful web browser gui to sort through the data. 
- Complete guide to set-up in progress. 



Directories:
- bin: contains binary called by sweep.py (MCMC_CUDA_10States)
- expData: .csv files with input values for force pCa curve (Force_pCa_Optmz_morePts.csv was used for this paper)
- obj: contains CUDA object files
- sweep: contains error files
- src: contains source code for model
    - main.cu: reads input params (csv_reader.cu/h), calls particles.cu
	- also calls gpuErrchk.h, problemDefines.h (set number of RUs, states, reps, timestep, etc), setGPU.h/cu
    - particles.cu: calculates max steady state force/twitch for different parameter combinations, minimizes cost function, writes out force     and states, calls force_pCa_curve.cu
    - force_pCa_curve.cu: calculates parameters, calls rates_trans_matrix.cu, update_RUs.cu, repeat_simul.cu
    - rates_trans_matrix.cu: calculate transition matrices (cooperative coefficients)
    - update_RUs.cu: updates state of each RU at each Markov step (based on probability)
    - repeat_simul.cu: calls rates_trans_matrix.cu, genrand.cu, update_RUs.cu, lin_interp_ca.cu (for twitch simulations; loops through     simulations
- SA_krecruit: simulation data for plotting supplemental figure S9
- States: simulation data for plotting supplemental figures S7 and S10
- MCMC_simulation_results: This is where all subsequent data runs are stored. Each new run with generate a datestamped folder. 
- resultsOrg: Helper functions for viewing, analyzing and working with simulations carried out from the MCMC model. 


## MYOCYTE (Spatially implicit sarcomere model)
- Read in digitized experimental data
	- Read in digitized experimental data from Regnier et al. 2004 figure 1a (force pCa data)
	- Interpolate data smoothly across timespan to match model
    - ATP_points.csv: ATP force pCa data (data points)
	  - ATP_Hill.csv: ATP force pCa data (Hill curve fit)
	  - dATP_points.csv: dATP force pCa data (data points)
	  - dATP_Hill.csv: dATP force pCa data (Hill curve fits)

- myocyte_model.m: Runs myocyte model (called by plot_script_myocyte_new.m)
*note: code requires parallel computing toolbox for use of parallel for loops (alternatively, this can be removed and normal for loops can be used)
  - Flags: 
	  - Ca_flag: specifies ATP or dATP calcium transient
	  - XB_protocol: specifies which type of simulation you want to run (force pCa or twitch)
	  - dATP_percent: specifies dATP percentage
    - Calls dXdT_myocyte_mechanics.m

- dXdT_myocyte_mechanics.m: differential equations for myofilament mechanics model (Lopez et al. 2020) - refer to Lopez et al. 2020 supplement for more details
    - Read in parameters specified in myocyte_model.m
    - Set additional parameters for sarcomere model
    - Adjust rates for metabolite concentrations
    - Set sarcomere length and calcium based on specifications in myocyte_model.m
    - For shortening simulations, uses calcium function   where parameters a, b, c, and Ca0 were optimized to match ATP and dATP experimental calcium transients
    - Compute active and passive forces, update state  variables

- pCa_calculate.m: calculates ec50, hill coefficient for force pCa simulations 
- fit_Ca_params.m: utilized to fit a, b, c, and Ca0

## VENTRICULAR
- CardiovascularMechanics.m: Runs ventricular model (called by plot_script_ventricular.m) 
    - Flags:
	    - Ca_flag: specifies ATP or dATP calcium transient
	    - HF_protocol: specifies healthy or failing simulation; for failing simulation, swaps in failing metabolite concentrations (as in     Lopez 2020)
	    - dATP_percent: specifies dATP percentage
    - Read in experimental rat data - for this study, I used the mean sham rat from Lopez 2020 
	        - Adjustable_parameters_table_rest.xlsx: contains data on geometry parameters, ATPase parameters etc (for all rats originally in   Lopez 2020 study)
	    - Adjustable_parameters_table_swap.xlsx: same as above but with healthy and failing rat values swapped 
	    - data1.xlsx: contains data on rat weight, echo data, metabolites (for all rats originally in Lopez 2020 study)
    - Run model
	  - Run cardiovascular mechanics model to steady state without coupled energetics model
	  - Run triseg model to get initial geometry (calls TrisegEquations.m)
	  - Run mechanics model (calls dXdT_cardiovascular_mechanics.m)
	  - Store output, initial calculations
	  - Run coupled cardiovascular and energetics model
	  - Run energetics model every 3 beats (calls EnergeticsModelScript.m)
	  - Run mechanics model (calls dXdT_cardiovascular_mechanics.m)
	  - Store output for every beat, concatenate at end

- EnergeticsModelScript.m: runs energetics model; taken directly from Lopez et al. 2020
    - calls dXdT_energetics.m
    - Outputs metabolite concentrations for input into mechanics model

- dXdT_energetics.m: differential equations for energetics model

- TrisegEquations.m: calculates ventricular geometry based on triseg equations (Lumens et al. 2009) *note: these equations are also found in dXdT_cardiovascular_mechanics.m; this script is only run at the beginning of the simulation to get the initial geometry

- dXdT_cardiovascular_mechanics.m: differential equations for mechanics model
    - Read in parameters specified in cardiovascular_model.m
    - Adjust rates for metabolite concentrations
    - Set calcium based on specifications in CardiovascularMechanics.m
	  - Uses calcium function where parameters a, b, c, and Ca0 were optimized to match ATP and dATP experimental calcium transients
    - Compute active and passive forces, changes in geometry in each segment (RV, LV, septum)

- Sensitivity_analysis.m: Utilized for supplemental figures S23, S25, S26
