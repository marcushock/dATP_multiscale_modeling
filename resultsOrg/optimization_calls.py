### This file is defined to have the necessary functions to call the within the python context, and clean and move around data 
### Such that it can be used within other optimization frameworks like MCMC sampling or Ax/BoTorch Optimization 


# Note that this approach is designed for carrying out a single call at a time, and ultiamtely trying to store results in a clearner way. 
# i.e. Each call of run_MCMC will lead to one simulation (OR MAYBE 2) if there is also a twitch. This is actually imoprtant to consider, but 
# should be added later I think. 
# Also note that this is hard coded to run a force pCa simulation always. A twitch will be run if the twitch file is provided. 

## Other misc todos: 
# [ ] - Make sure I don't overwrite data anywhere, either general results 
# [ ] - And also in whathever position I decide to save it it
# Lower priority todos:
# [ ] - LATER If twitch is implemented LATER, then decide if it's going to have it's own results and return value, or if it's just going to be part of the error metric.
# [ ] - Better management of the twitch simulation runs. 



import os
import sys
import shutil
from datetime import datetime
import subprocess
import time
import pandas as pd



# from /crucial/modified_MCMC/dATP_multiscale_modeling/Run_MCMC.py import run_mcmc
# Add the directory containing the module to sys.path
module_path = os.path.abspath("/crucial/modified_MCMC/dATP_multiscale_modeling/")
if module_path not in sys.path:
    sys.path.append(module_path)

# Import the module
from Run_MCMC import run_mcmc
from resultsOrg import helper_functions as hf

def combine_params(default_params, new_params):
    # Default needs to have stored or read in a full set of baseline parameters. 
    # The new params could be a dictionary or dataframe that has the new parameters we want to replace things with. 
    # Can return a df 
    new_param_df = default_params.copy()
    for key in new_params.keys():
        new_param_df[key] = new_params[key]
    return new_param_df 

def read_states_output(directory, exp_file, num_states=7):
    filename = os.path.join(directory, 'Rep_0States_out.csv')
    simulation = hf.states_structure(filename, num_states = num_states, skip_params = True, exp_file = exp_file)
    return simulation

def calcuate_error_metric(simulation, exp_file, type = 'SSE', normalization = None):
    exp_data = pd.read_csv(exp_file, names = ['pCa', 'Force'])
    # Need to scale the force_pCa first 
    if normalization is None: 
        force_pCa = simulation.force_pCa.values / simulation.force_pCa.max()
    else:
        print('Not implemented yet. Error!')
        exit(1)

    residues = exp_data['Force'].values - force_pCa
    if type == 'SSE':
        error_metric = sum(residues**2)
    elif type == 'residuals':
        error_metric = residues
    else:
        print('Error metric type not recognized!')
        exit(1)
    return error_metric
        
    

def evaluate_cuda_fit(trial_parameters, settings_dict):
    '''
    A function that can easily be called and ultimately returns an error measurement. 
    :param trial_parameters: A dictionary of the actual parameter that are changed and tested in the optimization 
    :param settings_dict: A detailed dictionary that contains the necessary information for the rest of the run. Likely defined outside the optimization loop
        Should contain: 
        - 'binary_path' = Path to the binary 
        - 'general_outdata' = Path to the General results dir 
        - 'new_savedata' = Path to the savedata location 
        - 'exp_force_pCa' = Experimental force pCa data path
        - 'exp_twitch' = Experimental twitch data (optional)
        - 'default_params' = Default parameter set (Should be a df)
        - 'temporary_parameter_file' = Name of temporary parameter file to write the combined parameters to (could be optional) [Not implemented yet]
        - 'code_src' = Path to the source code (optional) and probably won't be used much 
    
    '''
    # Unpack the settings dict so that they can be used. 
    binary_path = settings_dict['binary_path'] # Make sure it's the full path 
    general_outdata = settings_dict['general_outdata'] # /crucial/modified_MCMC/dATP_multiscale_modeling/MCMC_simulation_results/General_results
    new_savedata = settings_dict['new_savedata'] # Make sure it's the full path 
    exp_force_pCa_file = settings_dict['exp_force_pCa'] # Full path to experimental force pCa data
    exp_twitch_file = settings_dict['exp_twitch'] # Full path to experimental twitch data (optional, and likely will be None)
    default_params_df = settings_dict['default_params'] # Full path to default parameter set (CSV)
    code_src_dir = settings_dict['code_src'] # Full path to source code (optional, likely None)


    # Combine the default parameters with the trial parameters to create a full parameter set
    new_parameters_df = combine_params(default_params_df, trial_parameters)

    if not os.path.exists(new_savedata):
        os.makedirs(new_savedata)

    # Assert force pCa always 
    new_parameters_df['protocol'] = 1
    if exp_twitch_file is not None:
        new_parameter
    # Write the new parameters to a temporary CSV file
    new_parameter_file = os.path.join(new_savedata, 'temp_parameters.csv')
    new_parameters_df.to_csv(new_parameter_file, index=False)

    # Call run MCMC 
    # First ensure that the output directory exists


    raw_data_dir_output = run_mcmc(binary_path,
            exp_force_pCa_file,
            new_parameter_file,
            general_results_dir=general_outdata,
            output_results_dir=new_savedata, # This is going to be the imoprtant part of where we move things to. By default, I think it's just gonna make a new folder each time. 
            code_src=None)
    # Output will go into the "General_results" dir 
    # Then something will read the data from there and read:
    # - the parameters use
    # - the states data and create a states structure 
    sim = read_states_output(raw_data_dir_output, exp_file = exp_force_pCa_file, num_states =7)
    print(sim.force_pCa)
    # read_simulation()
    new_parameters_df['simulation'] = [sim]

    error_metric = calcuate_error_metric(sim, exp_force_pCa_file, type = 'SSE', normalization = None)

    return error_metric, new_parameters_df



