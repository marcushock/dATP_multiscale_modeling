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
import numpy as np 
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

def read_states_output(directory, exp_file, file_name = 'Rep_0States_out.csv', num_states=7):
    filename = os.path.join(directory, file_name)
    simulation = hf.states_structure(filename, num_states = num_states, skip_params = True, exp_file = exp_file)
    return simulation

def calcuate_error_metric(simulation, exp_file, type = 'SSE', scaling_factor = None):
    '''
    Docstring for calcuate_error_metric
    
    :param simulation: Simulation object that has been created and contains the force pCa
    :param exp_file: Filepath to the experimental csv for which this should be compared against
    :param type: Type of error to return. Options are 'SSE' for sum of squared errors, or 'residuals' for just the residuals
    :param scaling_factor: Value to normalize the simulation force to. If none, will autoscale to one. It should be a factor that when multiplied by the raw force, it will give it on a scale of 0 to 1. 
    Generally this should be the 1/max(raw force) of the simulation without any drug present. 
    Ex: 
    - No afi present, raw force = 0.3, then enter the normalziation = 1/0.3. = 3.33 
    '''

    exp_data = pd.read_csv(exp_file, names = ['pCa', 'Force'])
    # Need to scale the force_pCa first 
    if scaling_factor is None: 
        scaling_factor = 1 / simulation.force_pCa.max()
        force_pCa = simulation.force_pCa.values * scaling_factor
        assert abs(force_pCa.max() - 1.0) < 1e6, "Scaling factor did not properly scale the force to max of 1.0"
    else:
        force_pCa = simulation.force_pCa.values * scaling_factor
    

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
    A function that can easily be called and ultimately returns an error measurement in the form of a dictionary. 
    :param trial_parameters: A dictionary of the actual parameter that are changed and tested in the optimization 
    :param settings_dict: A detailed dictionary that contains the necessary information for the rest of the run. Likely defined outside the optimization loop
        Should contain: 
        - 'binary_path' = Path to the binary (Default value is set in function, but can be overwritten)
        - 'general_outdata' = Path to the General results dir (Default value is set in function, but can be overwritten)
        - 'new_savedata' = Path to the savedata location (required)
        - 'exp_force_pCa' = Experimental force pCa data path (required)
        - 'exp_SuperSlow_curve' = Experimental SuperSlow curve data path (optional, but required for additional fitting of SuperSlow ATPase curve) [Implementation in progress 2/5/26]
        - 'exp_twitch' = Experimental twitch data path (optional) 
        - 'default_params' = Default parameter set (Must be a df)
        - 'code_src' = Path to the source code (optional) and probably won't be used much 
        - 'fpCa_scaling_factor' = Normalization approach for force pCa (optional, but useful for drug optimizations) (should be 1/max(no drug force pCa))
        - 'L1_lambda_term' = Regularization strength for L1 regularization (optional, default = 0, meaning no regularization) 
        
    
    '''
    # Unpack the settings dict so that they can be used. 
    binary_path = settings_dict.get('binary_path', '/crucial/modified_MCMC/dATP_multiscale_modeling/bin/MCMC_CUDA_10States') # Path to the binary, this is required
    general_outdata = settings_dict.get('general_outdata', '/crucial/modified_MCMC/dATP_multiscale_modeling/MCMC_simulation_results/General_results') # Path to the general results directory, this is required
    new_savedata = settings_dict.get('new_savedata', None)  # Make sure it's the full path 
    exp_force_pCa_file = settings_dict.get('exp_force_pCa', None)  # Full path to experimental force pCa data
    exp_twitch_file = settings_dict.get('exp_twitch', None)  # Full path to experimental twitch data (optional, and likely will be None)
    default_params_df = settings_dict.get('default_params', None) 
    code_src_dir = settings_dict.get('code_src', None) # Full path to source code (optional, likely None)
    fpCa_scaling_factor = settings_dict.get('fpCa_scaling_factor', None) # Normalization approach for force pCa (optional)
    SuperSlow_curve_file = settings_dict.get('exp_SuperSlow', None) # Full PATH to experimental SuperSlow curve data (optional, likely None)
    L1_lambda = settings_dict.get('L1_term_lambda', 0) # Regularization strength for L1 regularization (optional, default = 0, meaning no regularization)


    if new_savedata is None:
        print("Error: new_savedata path must be provided in settings_dict")
        exit(1)
    if exp_force_pCa_file is None:
        print("Error: exp_force_pCa path must be provided in settings_dict")
        exit(1)
    if default_params_df is None or not isinstance(default_params_df, pd.DataFrame):
        print("Error: default_params as a must be provided in settings_dict as a pandas dataframe")
        exit(1)
    

    # Combine the default parameters with the trial parameters to create a full parameter set
    new_parameters_df = combine_params(default_params_df, trial_parameters)
    
    # If carrying out the super slow curve fitting, we need to read in the data, and lengthen the parameters df 
    if SuperSlow_curve_file is not None:
        super_slow_data, where_1_uM = read_SuperSlow_data(SuperSlow_curve_file)
        len_superslow = len(super_slow_data)
        new_parameters_df = pd.concat([new_parameters_df]*len_superslow, ignore_index=True)
        new_parameters_df['percent_drug'] = super_slow_data['drug_conc'].values # Apologies for different keys, but same meaning. 
    
    
    if not os.path.exists(new_savedata):
        os.makedirs(new_savedata)

    # Assert force pCa always 
    new_parameters_df['protocol'] = 1
    if exp_twitch_file is not None:
        new_parameter = 0 # Not implememnted yet 
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
    if SuperSlow_curve_file is None:
        sim = read_states_output(raw_data_dir_output, exp_file = exp_force_pCa_file, num_states =7)
        print(sim.force_pCa)
        # read_simulation()
        new_parameters_df['simulation'] = [sim]

        error_metric_fpCa = calcuate_error_metric(sim, exp_force_pCa_file, type = 'SSE', scaling_factor = fpCa_scaling_factor)
        print("Error metric for force pCa curve: ", error_metric_fpCa)
        
        L1_term = compute_L1_term(new_parameters_df, lambda_L1 = L1_lambda)
        # error_metric += L1_term
        return_dict = {'force_pCa_SSE': error_metric_fpCa, 'results_df': new_parameters_df, 'L1_term': L1_term}
        return return_dict
    elif SuperSlow_curve_file is not None:
        # Iterate through the number of simulations that have been run
        simulation_list = []
        for i in range(len_superslow):
            sim = read_states_output(raw_data_dir_output, exp_file = exp_force_pCa_file, num_states =7, file_name = f'Rep_{i}States_out.csv')
            simulation_list.append(sim)
            new_parameters_df.loc[i, 'simulation'] = sim
        sim_1_uM = simulation_list[where_1_uM]
        error_metric_fpCa = calcuate_error_metric(sim_1_uM, exp_force_pCa_file, type = 'SSE', scaling_factor = fpCa_scaling_factor)
        error_metric_superslow = calcualte_error_superslow_percentage(simulation_list, super_slow_data, metric_type = 'SSE')


        L1_term = compute_L1_term(new_parameters_df , lambda_L1 = L1_lambda)

        return_dict = {'force_pCa_SSE': error_metric_fpCa,
                       'superslow_SSE': error_metric_superslow,
                       'results_df': new_parameters_df,
                       'L1_term': L1_term}
        return return_dict

def calcualte_error_superslow_percentage(simulation_list, super_slow_data, metric_type = 'SSE'):
    '''
    Docstring for calcuate_error_superslow_percentage
    
    :param simulation_list: A list of simulation objects that have been created and contains the necessary data to compare against the SuperSlow curve. The order should be the same as the order of the super slow data points. 
    :param super_slow_data: A dataframe that has the drug concentrations and percent superslow for each point. 
    :param metric_type: Type of error to return. Options are 'SSE' for sum of squared errors, or 'residuals' for just the residuals
    '''
    residuals = np.zeros(len(simulation_list))
    for i, sim in enumerate(simulation_list):
        # For each simulation, we need to extract the percent superslow at maximal calcium, which is the last point in the force pCa curve. 
        min_index = sim.steady_states_all.index.min()
        sim_percent_superslow = sim.steady_states_all.loc[min_index,'SuperSlow'] # Note that the percentages are actually 0 to 1
        exp_percent_superslow = super_slow_data.iloc[i]['percent_superslow'] # Note that the percentages are actually 0 to 1
        residuals[i] = exp_percent_superslow - sim_percent_superslow
    if metric_type == 'SSE':
        error_metric = sum(residuals**2)
    elif metric_type == 'residuals':
        error_metric = residuals
    else:
        print('Error metric type not recognized!')
        exit(1)
    return error_metric
    # For each simulation, we need to extract the percent superslow at maximal calcium, which is the last point in the force pCa curve. 
    

### Super Slow via ATPase Curve Function Call and optimization 
# End goal for function: 
# ATPase_optimization(some_inputs) -> error_metric_single_force_pCa, error_metric_ATPase_curve, new_parameters_df (with the simulation object included)
# Overall, I think it'll make sense to use the same "evaluate cuda_fit" function, but then if there's ATPase data, then we run more steps... 
# In this optimization code, we assume that 1 uM aficamten is always going to be used for one of the points in the ATPase curve, and also be used for the force pCa
# General scheme: 
'''
1. Call the "evaluate cuda function" but have a setting for the Super Slow curve via ATPase data.
2. If: 
    - the super_slow curve is None, continue with normal optimization. 
   Elif SuperSlow curve is something: 
   [X]  First read in the SuperSlow data
   [X]  Then create a dataframe that has all the necessary afimcanten (or drug) concentrations 
   [X]  Write that CSV in the same way that the other CSV has been written. 
   [X]  Call the MCMC function to carry out those simulations, it will by default run all 4 (or however many points we have). 
   [X]  Read in each of the super slow states at maximal calcium, and compare it against the SuperSlow curve experimetnal data. 
   [X]  Compute an SSE metric for the SuperSlow curve 
   [X]  Then find the simulation that has the 1 uM aficamten, and compute the error metric for the force pCa curve.
   [X]  Return both error metrics, and the new parameters df with the simulation object for the force pCa curve.
'''

def read_SuperSlow_data(SuperSlow_curve_file):
    '''
    Docstring for read_SuperSlow_data
    
    :param SuperSlow_curve_file: Description
    Notes: 
    - The SuperSlow curve data should be in a CSV with two columns, one for drug concentration and one for percent SuperSlow state. The percent superslow should be between 0 and 1. (not 0 and 100)
    - If any data is great than one, it will automatically convert it to be between 0 and 1 by dividing by 100.
    '''
    # Note that this is actually storing the percent that's not superslow. 
    cols = ["drug_conc", "percent_superslow"]

    # Read CSV assuming there are no header names. 
    SuperSlow_data = pd.read_csv(SuperSlow_curve_file, names=cols)
    if SuperSlow_data['percent_superslow'].max() > 1:
        SuperSlow_data['percent_superslow'] = SuperSlow_data['percent_superslow'] / 100.0
    # SuperSlow_data['percent_superslow'] = 1 - SuperSlow_data['percent_superslow']
    where_1_uM = np.argwhere(SuperSlow_data['drug_conc'] == 1.0).flatten()[0]
    return SuperSlow_data, int(where_1_uM)

def compute_L1_term(parameters_df, lambda_L1 = 1.0):
    '''
    Docstring for compute_L1_term
    
    :param parameters_df: DataFrame of parameters, note, that including the simulation is optional. However, only the first row will be used 
    :param lambda_L1: Regularization strength for L1 regularization. 
    '''
    baseline_parmaeter_names = ['k_force_baseline','k_plus_SR_baseline','k1_plus_ref_baseline','k2_plus_baseline','k3_plus_baseline','k4_plus_ref_baseline']
    # This is the list of possible parameters that we might be changing
    # For these parameters, we care about the difference between the drug parameter and the baseline parameter, so we want to regularize the difference.
    drug_parameter_names = ['k_force_drug', 'k_plus_SR_drug','k1_plus_ref_drug','k2_plus_drug','k3_plus_drug','k4_plus_ref_drug']
    # These parameters, technically the no drug value is 0, so we just care about the value itself. 
    drug_only_names = ['k_plus_SS', 'k_minus_SS', 'k_plus_alt', 'k_minus_alt']
    # These parameters are necessary, but we might allow them to vary, therefore we don't want to regularize them. 
    no_regularization_parameters = ['K_D', 'coop_N'] # Probabyl wont' use this list 
    L1_term = 0
    for baseline_param, drug_param in zip(baseline_parmaeter_names, drug_parameter_names):
        baseline_value = parameters_df.loc[0, baseline_param]
        drug_value = parameters_df.loc[0, drug_param]
        L1_term += abs((drug_value - baseline_value))
    for drug_only_param in drug_only_names:
        drug_value = parameters_df.loc[0, drug_only_param]
        L1_term += abs(drug_value)
    return lambda_L1 * L1_term


