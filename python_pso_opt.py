#### Python file for optimazation of selected parameters
# Created 02/20/2025


# TODO
# [ ] - Create a funciton to write csv file that has all the simulations that we want to run (not sure actually about this...)
# [X] - Readin in the some starting values as a base location for the PSO algorithm 
# [ ] - Create a function to run the PSO algorithm 
# [ ] - Calculate the error that can be used when returning from the PSO function and have a minimzation target 
# [ ] - Maybe have a temp csv file that gets overwritten every PSO run and then also have a master CSV that shows all the particle runs 
# [ ] - Select whether to use the ATPase data or the force-pCa data 
# [ ] - Decide on how the drug amount will be estimated in the simulations (probably start with )
# [ ] - Made a changable normalization. Or normalization using multiple data sets. 
# [ ] - Load in the ATPase data 


# Import libraries
import numpy as np 
import pyswarms as ps
import pandas as pd
import os
import subprocess
import time
import numpy as np 
import datetime
from resultsOrg import helper_functions as hf
import yaml 
import sys 
from pyswarms.single.local_best import LocalBestPSO









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


# Define function to read in the yaml from the input line argument 
def read_yaml_file():
    if len(sys.argv) < 2:
        raise ValueError("No command line YML file path provided")

    file_path = sys.argv[1]
    # file_path = "/crucial/modified_MCMC/dATP_multiscale_modeling/resultsOrg/optiziation_dir/example_yaml_in.yml"
    """
    Load PSO configuration from a YAML file with compact [baseline, min, max] format
    Returns:
        dict: Dictionary with processed PSO configuration
        
    Raises:
        FileNotFoundError: If the specified file does not exist
        yaml.YAMLError: If there is an error parsing the YAML file
    """
    try:
        with open(file_path, 'r') as file:
            config = yaml.safe_load(file)
        
        if config['search_space'] is None: 
            config['search_space'] = {}
        else:
            for param_name, param_values in config['search_space'].items():
                if len(param_values) != 3:
                    raise ValueError(f"Parameter {param_name} should have exactly 3 values [baseline, min, max]")
                        

        return config
        
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file not found: {file_path}")
    except yaml.YAMLError as e:
        raise yaml.YAMLError(f"Error parsing YAML file: {e}")


def rename_undated(results_directory, extra_append = ""):
	file_names = os.listdir(results_directory)
	time_str = time.strftime("%y%m%d-%H%M_")+extra_append
	for name in file_names:
		if "States_" == name[0:7] or "Force_" == name[0:6]:
			new_filename = time_str+fix_trailing_zeros(name)
			print("Renaming '{}' \nto \n'{}'\n".format(name, new_filename))
			original_path = results_directory+'/'+name
			new_path = results_directory + '/' + new_filename
			try:
				os.rename(original_path, new_path)
			except:
				print('Could not rename. ')
	return 

def fix_trailing_zeros(filename_string):
    new_name = ''
    for filename_piece in filename_string.strip('.csv').split(' '):
        try:
            new_value = float(filename_piece)
            new_name += str(new_value)
            new_name += ' '
        except:
            new_name += filename_piece
            new_name += ' '
    if new_name[-1] == ' ':
        new_name = new_name[0:-1]
    new_name += '.csv'
    return new_name


def write_csv_file(input_dict, filename = "MCMC_temp_input.csv", append = False):
    '''
    This function will write a CSV file that has the input parameters that are going to be optimized. 
    The rest of the values will be the same as the default values. 

    '''
    if append == False:
        f = open(filename, "w")
        count = 0
        new_line = ""
        for key in input_dict.keys():
            new_line += str(key)
            new_line += ","
        new_line = new_line[:-1] + "\n"
        f.write(new_line)
        new_line = ""
        for key in input_dict.keys():
            new_line += str(input_dict[key])
            new_line += ","
        new_line = new_line[:-1] + "\n"
        f.write(new_line)
        f.close()
        print("Wrote temporary input file to: ", filename)

    elif append == True:
        try:
            f = open(filename, "r")
            f.close()
            write_header = False
        except:
            print("File does not exist, creating new file. ")
            write_header = True
        f = open(filename, "a")
        if write_header == True:
            new_line = ""
            for key in input_dict.keys():
                new_line += str(key)
                new_line += ","
            new_line = new_line[:-1] + "\n"
            f.write(new_line)
        new_line = ""
        for key in input_dict.keys():
            new_line += str(input_dict[key])
            new_line += ","
        new_line = new_line[:-1] + "\n"
        f.write(new_line)
        f.close()
        print("Appended input file to: ", filename)
    return filename

def read_default_params(filename):
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

    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    new_directory = storing_directory + current_time + custom_name
    try: 
        os.mkdir(new_directory)
    except:
        print("Directory already exists... given timestamp. ")
        print("Data may be overwritten. ")
    return new_directory


def run_MCMC_bin(input_dict, output_dir, bin_namepath, exp_data_file):
    '''
    This function takes in a dictionary which is used to specify the input parameters that are going to be optimized. 
    The rest of the input arguments will be left untouched. 
    '''
    temp_running_file = write_csv_file(input_dict, filename = output_dir+"/PSO_temp_input.csv")
    log_all_runs = write_csv_file(input_dict, filename = output_dir+"/PSO_all_runs.csv", append = True)
    # Run the MCMC simulation 
    # Copied below from other file 
    
    # Exp Data path: 
    # Should be fixed to be more easily modified later rather than hard coded
    # exp_data_path = "/crucial/modified_MCMC/dATP_multiscale_modeling/expData/4.csv"
    ### FIX LATER ### 
    # This file path should be changed thoughtfully 
    args = [
        bin_namepath,
        exp_data_file,
        temp_running_file]

    os.chdir(output_dir)

    try:
        os.makedirs("PSO_results")
    except OSError:
        pass
    results_path = output_dir + "/PSO_results"


    with open("stdout.out", 'w') as out:
        with open("stderr.err", 'w') as err:
            print("Running " + str(args))
            time_start = time.time()
            subprocess.call(args, stdout=out, stderr=err, shell=False, cwd=output_dir)
            print("Time (seconds): ", time.time() - time_start)
    return results_path



def run_and_evalutate(new_full_parameter_set, config_dirs):
    output_directory = config_dirs['output_directory']
    binary_name = config_dirs['bin_name']
    exp_data_type = config_dirs['exp_data_type']
    exp_data_file = config_dirs['exp_data_file']


    # Run the binary
    results_dir = run_MCMC_bin(new_full_parameter_set, output_directory, binary_name, exp_data_file)
    # results_dir = '/crucial/modified_MCMC/dATP_multiscale_modeling/MCMC_simulation_results/PSO_optimizations/2025-02-25_06-22'+"/PSO_results"
    # Load in the exp data 
    exp_data = pd.read_csv(exp_data_file, names = ['pCa','Normalized_data'])
    
    # Load in the simulation data 
    if exp_data_type == "ATPase":
    #     # Load in the ATPase data 
        for f in os.listdir(results_dir):
            if "Rep_0ATP_out" in f:
                simulation_data = pd.read_csv(results_dir+'/'+f, names = ['Timestep']+list(exp_data.pCa))
                break
    
    if exp_data_type == "Force_pCa":
        for f in os.listdir(results_dir):
            if "Rep_0Force_out" in f:
                simulation_data = pd.read_csv(results_dir+'/'+f, names = ['Timestep']+list(exp_data.pCa))
                break 

    selected_times = simulation_data.Timestep>simulation_data.Timestep.max()*1/3
    sim_steady_states_values = simulation_data[selected_times][list(exp_data.pCa)].mean()
    # Define new steady state DF with each column a different state
    # And each row a different pCa (7,4,0.1 steap)
    # Also creating a place for the pCa first, then renaming 

        
    # Calculate the error between the two 
    sim_normalized = normalize_array(sim_steady_states_values.values)
    print("Sim normalized: ", sim_normalized)
    exp_normalized = normalize_array(exp_data.Normalized_data.values)
    print("Exp normalized: ", exp_normalized)
    errors = sim_normalized - exp_normalized
    rmse = np.sum((errors)**2)
    print("RMSE is: ", rmse)
    print("For the parameters: ", new_full_parameter_set)
    rename_files(results_dir)
    # Rename any file in the directory with the start of the file of "Rep0" to Particle_N_ where N starts 
    # at one and increases by one if that file exists already

    return rmse
    # Return the error 

def rename_files(directory):
    files = [f for f in os.listdir(directory) if "Rep_0" in f]
    files.sort()
    
    existing_files = set(os.listdir(directory))
    
    particles = [int(f.split("particle_")[1].split("_")[0]) for f in os.listdir(directory) if "particle_" in f]
    try:
        counter = max(particles)+1
    except:
        counter = 1
    
    for file in files:

        old_path = os.path.join(directory, file)
        new_name = fix_trailing_zeros(f"particle_{counter:04d}_"+file.strip("Rep_0"))
        new_path = os.path.join(directory, new_name)
        
        os.rename(old_path, new_path)
        existing_files.add(new_name)
        
        print(f"Renamed {file} -> {new_name}")

def normalize_array(input_array):
    '''
    This function will normalize the input array to have a max of 1, and min of 0 '''
    data_min = np.min(input_array)
    data_max = np.max(input_array)
    normalized_array = (input_array - data_min) / (data_max - data_min)
    return normalized_array

def pyswarm_run(X_scaled, yaml_config):
    default_parameters = read_default_params(yaml_config['files_and_directories']['default_param_file'])


    parameter_space_bounds = np.zeros((len(yaml_config['search_space']), 2))
    i = 0
    for key in yaml_config['search_space']:
        parameter_space_bounds[i] = yaml_config['search_space'][key][1:3]
        i+=1

    X = log_denormalize(X_scaled, parameter_space_bounds)
    
    results_array = np.zeros(X.shape[0]) # .reshape(-1,len(yaml_config['search_space']))

    # X has the shape (n_particles, n_features/parameters)
    for i in range(X.shape[0]):
        parameters = X[i]
        optimization_range_dict = yaml_config['search_space']
        if optimization_range_dict is None:
            optimization_range_dict = {}

        optimization_values = {}
        k = 0
        for key, value in optimization_range_dict.items():
            optimization_values[key] = parameters[k]
            k+=1
        # Optimization values is going to be a small dict
        new_full_parameter_set = default_parameters.copy()
        new_full_parameter_set.update(optimization_values)
        results_array[i] = run_and_evalutate(new_full_parameter_set, yaml_config['files_and_directories'])
    return results_array


def write_best_params_to_csv(position, cost, yaml_config):
    default_parameters = read_default_params(yaml_config['files_and_directories']['default_param_file'])
    new_full_parameter_set = default_parameters.copy()

    for param, key in zip(position, yaml_config['search_space'].keys()):
        new_full_parameter_set[key] = param 
        
    write_csv_file(new_full_parameter_set, filename = yaml_config['files_and_directories']['output_directory']+"/PSO_best_params.csv")


def log_normalize(x, orig_bounds):
    """ 
    Convert parameters to log-space and normalize to [0,1].
    
    Parameters:
        x (array): Original parameter values.
        orig_bounds (array): 2D array of min/max bounds for each parameter. Shape (n_params, 2)
    
    Returns:
        array: Normalized log-space values in [0,1].
    """
    log_bounds = np.log10(orig_bounds)  # Convert bounds to log-space
    log_x = np.log10(x)  # Convert values to log-space
    return (log_x - log_bounds[:, 0]) / (log_bounds[:, 1] - log_bounds[:, 0])  # Normalize to [0,1]

def log_denormalize(x_norm, orig_bounds):
    """ 
    Convert normalized parameters [0,1] back to original scale.
    
    Parameters:
        x_norm (array): Normalized parameter values in [0,1].
        orig_bounds (array): 2D array of min/max bounds for each parameter. Shape (n_params, 2)
    
    Returns:
        array: Parameters converted back to original scale.
    """
    log_bounds = np.log10(orig_bounds)  # Convert bounds to log-space
    log_x = x_norm * (log_bounds[:, 1] - log_bounds[:, 0]) + log_bounds[:, 0]  # Denormalize from [0,1] to log-space
    return 10**log_x  # Convert back to original scale

def set_particle_1(x_init_array, yaml_config):
    '''
    This function will set the first particle to be the starting values that are specified in the yaml file. 
    '''
    x_init = x_init_array.copy()
    i = 0
    for key in yaml_config['search_space']:
        x_init[0][i] = yaml_config['search_space'][key][0]
        i+=1
    return x_init


def main():
    # Read in yaml file 
    yaml_config = (read_yaml_file())
    
    # Read in the default parameters CSV file 
    default_parameters = read_default_params(yaml_config['files_and_directories']['default_param_file'])
    print("Read in default parameters are: \n", default_parameters)

    # Make new directory for the simulation to occur 
    output_directory = make_exp_dir()
    print("Output directory is: ", output_directory)
    yaml_config['files_and_directories']['output_directory'] = output_directory 


    #### This chunk of code can probably go into a different function later ### 
    # Take the small dict of values that have been chosen to be optimized, and just get the first value in 
    # the list which is the starting value 
    # Read in the yaml search space 
    # optimization_range_dict = yaml_config['search_space']
    # if optimization_range_dict is None:
    #     optimization_range_dict = {}

    n_dim = len(yaml_config['search_space'])
    n_particles = yaml_config['pso_algorithm']['n_particles']

    x_min = np.zeros(n_dim)
    x_max = np.ones(n_dim)
    x_init = np.zeros((n_particles, n_dim))
    x_init_0 = np.zeros(n_dim)
    i = 0
    for key in yaml_config['search_space']:
        print("Key: ", key)
        print("Value: ", yaml_config['search_space'][key])
        x_init_0[i] = yaml_config['search_space'][key][0]
        # x_init_0[i] = log_normalize(yaml_config['search_space'][key][0], yaml_config['search_space'][key][1:3])
        x_max[i] = yaml_config['search_space'][key][2]
        x_min[i] = yaml_config['search_space'][key][1]
        x_init[:,i] = np.random.uniform(x_min[i], x_max[i], size = yaml_config['pso_algorithm']['n_particles']) 
        i+=1

    x_init[0] = x_init_0
    x_init_normalized = log_normalize(x_init, np.array([x_min, x_max]).T)
    assert x_init_normalized.shape == x_init.shape

    bounds = (np.zeros(n_dim), np.ones(n_dim))
    bounds_working = (x_min, x_max)

    # Instantiate the optimizer for a 1D problem
    options = yaml_config['pso_algorithm']['swarm_options']


    # optimizer = ps.single.LocalBestPSO(n_particles=yaml_config['pso_algorithm']['n_particles'],
    optimizer = ps.single.GlobalBestPSO(n_particles=yaml_config['pso_algorithm']['n_particles'],
                                        dimensions=n_dim,
                                        options = options,
                                        bounds = bounds, 
                                        init_pos = x_init_normalized)

    # now run the optimization, pass a=1 and b=100 as a tuple assigned to args

    cost, pos = optimizer.optimize(pyswarm_run, 
                                iters = yaml_config['pso_algorithm']['n_iterations'],
                                yaml_config = yaml_config)
    
    print("Cost: ", cost)
    pos = log_denormalize(pos, np.array([x_min, x_max]).T)
    print("Position: ", pos)
    write_best_params_to_csv(pos, cost, yaml_config)
    np.save("position_history.npy", optimizer.pos_history)  
    return 


if __name__ == "__main__":
    main()
