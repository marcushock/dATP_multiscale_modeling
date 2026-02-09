### Created on 7/15/2025 by MTH
# Designed with the intent to make it easier to run the binary for use with Sumatra, and just in general for better reproducibility. 
# Note that the script assumes everything is being called from the root directory of the repository.
# Additional changes must be made if the function is being called from other directories.
# Also note, that MCMC code is always going to write to the 



#
# Expected usage: 
# python Run_MCMC.py experimental_force_pCa.csv input_params.csv
#



import os
import sys
import shutil
from datetime import datetime
import subprocess
import time


def run_mcmc(binary,
            exp_data_input,
            parameter_filename,
            general_results_dir='MCMC_simulation_results/General_results',
            output_results_dir='MCMC_simulation_results/',
            code_src=None):
    """Run the MCMC binary and collect results.

    Parameters:
    - binary: path to the binary executable
    - exp_data_input: path to experimental data CSV
    - parameter_filename: path to parameter CSV
    - general_results_dir: directory where the binary writes general results (this might need to be a full path) but is defined in the MCMC cuda code
    - output_results_dir: base directory to store timestamped results - This is the thing that might need to be modified. 
    - code_src: path to source code directory to copy into results (defult none, and won't copy the source) but does in the main() function 
    """
    # Check if the current general directory is empty: 
    # If not empty, exit with an error message. 
    if os.listdir(general_results_dir):
        print(f"General results directory {general_results_dir} not is empty. Make sure to empty the directory and move files before running this script. ")
        sys.exit(1)
    else:
        print(f"General results directory {general_results_dir} is empty. Proceeding with the run.")
        run_dir = general_results_dir

    # Get the current time to make unique storage places. 
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M")
    print(timestamp)
    destination_dir = os.path.join(output_results_dir, timestamp)
    print(destination_dir)

    # Create a results directory that will have a unique name. 
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)
        print(f"Created directory {destination_dir} for results.")
    else:
        print(f"Directory {destination_dir} already exists. Please remove it or choose a different name.")
        print("Going to add in the seconds...")
        timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
        destination_dir = os.path.join(output_results_dir, timestamp)
        if not os.path.exists(destination_dir):
            os.makedirs(destination_dir)
        else:
            print(f"Directory {destination_dir} still exists. Please remove it or choose a different name.")
            sys.exit(1)

    # Copy the src directory to the destination directory
    if code_src is not None:
        shutil.copytree(code_src, os.path.join(destination_dir, os.path.basename(code_src)))
        print(f"Copied source code to {os.path.join(destination_dir, os.path.basename(code_src))}")

    # Copy the input experimental file and input parameters to the destination direcotry 
    shutil.copy(exp_data_input, destination_dir)
    shutil.copy(parameter_filename, destination_dir)
    print(f"Copied input files to {destination_dir}")

    # Create a raw_data directory inside the destination directory
    raw_data_dir = os.path.join(destination_dir, "raw_data")
    os.makedirs(raw_data_dir, exist_ok=True)
    print(f"Created raw_data directory at {raw_data_dir}")

    print("Running the binary at " + binary)
    subprocess.run([binary, exp_data_input, parameter_filename], check=True)
    print("Binary run completed.")

    # Move all the files from the general results directory to the raw_data directory
    for item in os.listdir(general_results_dir):
        src = os.path.join(general_results_dir, item)
        dst = os.path.join(raw_data_dir, item)
        shutil.move(src, dst)
    return raw_data_dir


def main():
    ### CHANGE THESE PARAMETERS IF NECESSARY
    # Currently using relative paths, but you can change these to absolute paths if you prefer.
    general_results_dir = 'MCMC_simulation_results/General_results'
    # Change this if yo you want to move the results somewhere else. 
    output_results_dir = 'MCMC_simulation_results/'
    code_src = 'src'
    BINARY = 'bin/MCMC_CUDA_10States'

    if len(sys.argv) == 3:
        exp_data_input, parameter_filename = sys.argv[1], sys.argv[2]
        inputyaml = None
    elif len(sys.argv) == 4:
        exp_data_input, parameter_filename, inputyaml = sys.argv[1], sys.argv[2], sys.argv[3]
    else:
        print("ERROR. Incorrect number of argument. \nUsage: python Run_MCMC.py experimental_force_pCa.csv input_params.csv")
        sys.exit(1)

    # Call the refactored function using the module-level defaults
    run_mcmc(BINARY, exp_data_input, parameter_filename, general_results_dir=general_results_dir,
             output_results_dir=output_results_dir, code_src=code_src)


if __name__ == "__main__":
    main()