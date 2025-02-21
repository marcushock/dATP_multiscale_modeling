#### Python file for optimazation of selected parameters
# Created 02/20/2025


# TODO
# [ ] - Create a funciton to write csv file that has all the simulations that we want to run 
# [X] - Readin in the some starting values as a base location for the PSO algorithm 
# [ ] - Create a function to run the PSO algorithm 
# [ ] - Calculate the error that can be used when returning from the PSO function and have a minimzation target 
# [ ] - Load in the ATPase data 
# [ ] - Maybe have a temp csv file that gets overwritten every PSO run and then also have a master CSV that shows all the particle runs 
# [ ] - Select whether to use the ATPase data or the force-pCa data 
# [ ] - Decide on how the drug amount will be estimated in the simulations (probably start with )

# Import libraries
import numpy as np 
import pyswarms as ps
import pandas as pd
import os
import subprocess
import time
import numpy as np 
import datetime




# Original dict - uncomment to use [starting value, lower bound upper bound]
optimization_dict = {
#  'protocol': [1.0, 0.5, 1.5], # Should not ever be optimized. Either can be 0 (twich) or 1 (force pCa)
#  'k_force_baseline': [0.2, 0.1, 0.3], 
#  'k_force_drug': [779.0, 0, 10000], 
#  'k_plus_SR_baseline': [16.0, 5.0, 20.0], 
#  'k_plus_SR_drug': [16.0, 10.0, 20.0], 
#  'k_minus_SR': [15.0, 10.0, 20.0], 
#  'k_xb': [5.0, 1.0, 10.0], 
#  'k1_plus_ref_baseline': [0.0025, 0.001, 0.005], 
#  'k1_plus_ref_drug': [0.00478, 0.002, 0.01], 
#  'k2_plus_ref_baseline': [0.0015, 0.001, 0.005], 
#  'k2_plus_ref_drug': [0.0015, 0.001, 0.005], 
#  'k3_plus_baseline': [0.05, 0.01, 0.1], 
#  'k3_plus_drug': [0.08, 0.01, 0.1], 
#  'k4_plus_ref_baseline': [0.135, 0.1, 0.2], 
#  'k4_plus_ref_drug': [0.23, 0.1, 0.3], 
#  'kB_plus_ref': [13.0, 10.0, 20.0], 
#  'kB_minus_ref': [0.1, 0.05, 0.2], 
#  'kCa_plus_ref': [0.09, 0.05, 0.2], 
#  'kCa_minus_ref': [0.57, 0.5, 1.0], 
#  'percent_drug': [0.0, 0.0, 1.0], 
#  'lambda': [0.0, 0.0, 1.0], 
#  'gamma_B': [45.0, 30.0, 60.0], 
#  'gamma_M': [21.0, 10.0, 30.0], 
#  'mu_B': [21.0, 10.0, 30.0], 
#  'mu_M': [2.0, 1.0, 5.0], 
#  'q': [1.0, 0.5, 1.5], 
#  'r': [1.0, 0.5, 1.5], 
#  'x_preR': [0.0, -0.1, 0.1], 
#  'x_xb': [0.075, 0.05, 0.1], 
#  'conc_ADP': [30.0, 10.0, 50.0], 
#  'conc_ATP': [3000.0, 1000.0, 5000.0], 
#  'conc_Pi': [3000.0, 1000.0, 5000.0], 
#  'delta_G_ATP': [-13.0, -20.0, -5.0], 
#  'alpha': [0.28, 0.1, 0.5], 
#  'beta': [0.35, 0.1, 0.5], 
#  'eta': [0.68,]
#  'g_Cb': [0.0, 0, 0], 
#  'g_Ca': [0.0, 0, 0]
}


'''
General program flow: 
The user can define which variables they want to optimize and the starting values (and potentially a range of values to optimize within)
Then once the user selects these parametes to optimize, the program will create an initial CSV file with values for the particles
that are going to be optimized, and the rest of the values will be the same as the default values.
Then the PSO algorithm will start running the simulations in a linear fashion 
because the entire GPU is used up by the simulations.

First priority is to just create a function that can take some specified parameters, and run the simulation, then return the error 
from the ATPase data. 
 - Probably need to create a new directory for each run 
 - Store the output state files 
 - Read in the ATP
'''

def run_MCMC(input_dict, exp_data = "ATPase"):
    '''
    This function takes in a dictionary which is used to specify the input parameters that are going to be optimized. 
    The rest of the input arguments will be left untouched. 

    '''
    write_csv_file(input_dict, filename = "MCMC_temp_input.csv")
    # Run the MCMC simulation 
    # Copied below from other file 
    '''
    args = [
		"bin/" + str(binName),
		"expData/Force_pCa_Optmz_morePts.csv", #Force_pCa_3pts_2 #Force_pCa_Optmz_morePts
		str(parameters["gamma_B"]),
		str(parameters["gamma_M"]),
		str(parameters["mu_M"]),
		str(parameters["k2_plus_ref"]),
		str(parameters["k3_plus"]),
		str(parameters["k4_plus_ref"]),
		str(parameters["kB_plus_ref"]),
		str(parameters["kB_minus_ref"]),
		str(parameters["lambda"]),
		str(parameters["kCa_plus_ref"]),
		str(parameters["kCa_minus_ref"]),
		str(parameters["percent_drug"]),
		str(parameters["k_force"]),
		str(parameters["k_plus_SR_ref"]),
		str(parameters["k_minus_SR_ref"]),
		str(parameters["protocol"]),
		 ]
	try:
		os.makedirs(folder)
	except OSError:
		pass
	fname = folder+"/"
	for key, value in parameters.items():
		fname += str(key) + " " + str(value) + " "
	fname = fname.strip()
	with open(fname+".out", 'w') as out:
		with open(fname+".err", 'w') as err:
			print("Running " + str(args))
			time_start = time.time()
			subprocess.call(args, stdout=out, stderr=err, shell=False)
			print("Time (seconds): ", time.time() - time_start)
    '''

    return 


def write_csv_file(input_dict, filename = "MCMC_temp_input.csv"):
    '''
    This function will write a CSV file that has the input parameters that are going to be optimized. 
    The rest of the values will be the same as the default values. 

    '''

    return

def read_default_params(filename = "/crucial/modified_MCMC/dATP_multiscale_modeling/parameters/parameter_input.csv"):
    '''
    This function will read in the default parameters that are going to be used in the simulations. 
    The user can specify the default parameters in a CSV file that is read in by this function. 
    If the user does not specify a default parameter, then the default value will be used. 

    '''
    parameters = pd.read_csv(filename, na_filter=True, comment='#')
    parameters.dropna(inplace=True)
    # Convert the dataframe (series) to a dictionary
    parameter_dict = parameters.iloc[0].to_dict()

    return parameter_dict




def make_exp_dir(custom_name = "", storing_directory = "/crucial/modified_MCMC/dATP_multiscale_modeling/MCMC_simulation_results/PSO_optimizations/"):
    '''
    This function will create a new directory for the experiment that is being run. 
    The directory will be created in the storing_directory that is specified by the user. 
    The name of the directory will be the current time that the function is called. 
    '''
    
    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    new_directory = storing_directory + current_time + custom_name
    os.mkdir(new_directory)
    return new_directory



def main():

    default_parameters = read_default_params()
    print("Read in default parameters are: \n", default_parameters)
    output_direcoty = make_exp_dir()
    print("Output directory is: ", output_direcoty)
    
    # This will eventually need to be inside of the PSO alogirthm 
    run_MCMC(default_parameters)

    return 


if __name__ == "__main__":
    main()
