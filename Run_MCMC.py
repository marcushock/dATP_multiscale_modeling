### Created on 7/15/2025 by MTH
# Designed with the intent to make it easier to run the binary for use with Sumatra, and just in general for better reproducibility. 

### CHANGE THESE PARAMETERS IF NECESSARY
# Currently using relative paths, but you can change these to absolute paths if you prefer.
general_results_dir = 'MCMC_simulation_results/General_results'
# Change this if yo you want to move the results somewhere else. 
output_results_dir = 'MCMC_simulation_results/'
code_src = 'src'
BINARY = 'bin/MCMC_CUDA_10States'

import os
import sys
import shutil
from datetime import datetime
import subprocess
import time


def main():
    if len(sys.argv) == 3:
         exp_data_input, parameter_filename = sys.argv[1], sys.argv[2]
    elif len(sys.argv) == 4:
        exp_data_input, parameter_filename, inputyaml = sys.argv[1], sys.argv[2], sys.argv[3]
    else:
        print("ERROR. Incorrect number of argument. \nUsage: python Run_MCMC.py experimental_force_pCa.csv input_params.csv")
        sys.exit(1)

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

    print("Runnign the binary at " + BINARY)
    args = [BINARY, exp_data_input, parameter_filename]
    subprocess.run([BINARY, exp_data_input, parameter_filename], check=True)
    # with open("stdout.out", 'w') as out:
    #         with open("stderr.err", 'w') as err:
    #             print("Running " + str(args))
    #             time_start = time.time()
    #             subprocess.call(args, stdout=out, stderr=err, shell=False, 
    #                             # cwd=output_dir
    #                             )
    #             print("Time (seconds): ", time.time() - time_start)
    print("Binary run completed.")

    # Move all the files from the general results directory to the raw_data directory
    for item in os.listdir(general_results_dir):
        src = os.path.join(general_results_dir, item)
        dst = os.path.join(raw_data_dir, item)
        shutil.move(src, dst)

    # shutil.move("output.csv", run_dir)


if __name__ == "__main__":
    main()