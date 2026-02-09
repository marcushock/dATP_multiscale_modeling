import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize as opt
import os 
from tqdm import tqdm

### TODO #### 
# 1. Fix the other auxillary functions to work with 5 or 4 states rather than just 1 
#    [ ] Get twitch 
#    [ ] Get steady statep
#    [ ] Get plot_twitch_from_sim
#    [ ] Get get_twitch_from_sim


class states_structure:
    def __init__(self, input_filename:'CSV', num_states = 6, skip_params = False, exp_file = '/crucial/modified_MCMC/dATP_multiscale_modeling/expData/Force_pCa_Optmz_morePts.csv') -> None: # type: ignore
        '''
        By default this code will assume we are using the updated 6 state model that 
        has 3 intermediates for the XB (M1, M2, M3). Use num_states = 5 for the 5
        state model or numstates = 4 for the 4 state model.
        '''
        self.exp_file = exp_file
        self.twitches = None
        self.exp_data = pd.read_csv(self.exp_file, names = ['pCa', 'Measurement'], header = None)
        # print('The filename is: ', input_filename)
        self.file_name = input_filename
        if num_states == 6:
            self.state_list = ['M3','M2','M1','C','B','OFF']
        elif num_states == 7:
            self.state_list = ['M3','M2','M1','C','B','OFF','SuperSlow']
        elif num_states == 5:
            self.state_list = ['M2','M1','C','B','OFF']
        else:
            self.state_list = ['M2','M1','C','B']
        self.read_states()
        if not skip_params:
            try:
                self.extract_parameters()
            except Exception as e:
                print('Failed to extract parameters from the file name')
                print('Try running with the skip_params = True option')
                print(f'Error: {e}')

        self.states_steadystate()

        return 
    
    def read_states(self):
        pd.read_csv(self.file_name, header = None)
        columns_list = ['Time'] # unit unknown 
        for pCa in self.exp_data.pCa:
            for state in self.state_list:
                new_item = ('{} {}'.format(state,pCa))
                columns_list.append(new_item)

        self.states_df = pd.read_csv(self.file_name, names = columns_list)
        return 

    def extract_parameters(self):
        # This assumes that there is "_States_out " preceeding the first set of parameters 
        self.parameters = {}
        param_half_string = self.file_name.split('States_out ')[1].strip(".csv")
        params_list = param_half_string.split(' ')
        i = 0 
        while i < len(params_list):
            try:
                self.parameters[params_list[i]] = float(params_list[i+1])
            except:
                try:
                    print('Flipping order for parameters')
                    self.parameters[params_list[i+1]] = float(params_list[i])
                except:
                    print('Failed to make with values {} and {}'.format(params_list[i], params_list[i+1]))
                    i +=1
                    next 
            i += 2
        # print('Read in parameters', self.parameters)
        return 

    def get_twitch(self):
        '''
        This function returns the temporal force states from the simulation and does not 
        average them across the replciaets. 
        '''
        force_cols = [col for col in self.states_df.columns if 'M3' in col or 'M2' in col] # 'M3 7.0 and M2 7.0'
        list_of_nums = np.unique([element[1] for element in self.states_df.columns[1:].str.split(' ')]) # self.states_df.columns.str(lambda x: x.split(' ')[1])
        new_twitch_df = pd.DataFrame(np.zeros((len(self.states_df.index), len(list_of_nums))), 
                                                    index = self.states_df.index,
                                                    columns = list_of_nums)

        for number in list_of_nums: 
            specific_cols = [col for col in force_cols if number in col]
            new_twitch_df[number] = self.states_df[specific_cols].sum(axis = 1)
        self.twitches = new_twitch_df
        return new_twitch_df

    def states_steadystate(self, end_time_amount= 500): 
        temp_series = self.states_df[self.states_df.Time>end_time_amount].mean(axis = 0)
        std_series = self.states_df[self.states_df.Time>end_time_amount].std(axis = 0)
        # Define new steady state DF with each column a different state
        # And each row a different pCa (7,4,0.1 steap)
        # Also creating a place for the pCa first, then renaming 
        num_rows = len(self.exp_data.pCa)
        num_cols = len(self.state_list) + 1
        df = pd.DataFrame(np.zeros((num_rows,num_cols)), columns = ['pCa'] + self.state_list)
        df.pCa = self.exp_data.pCa
        df = df.set_index('pCa')

        df_std = pd.DataFrame(np.zeros((num_rows,num_cols)), columns = ['pCa'] + self.state_list)
        df_std.pCa = self.exp_data.pCa
        df_std = df_std.set_index('pCa')

        # Annoying rounding necessary for the python decimal storage necessity
        for pCa in self.exp_data.pCa:
            for state in self.state_list:
                value = temp_series.get(f'{state} {pCa}')
                df.loc[pCa, state] = value 

                value = std_series.get(f'{state} {pCa}')
                df_std.at[round(pCa,2), state] = value 
        
        self.steady_states_all  = df
        self.steady_states_all_std  = df_std
        if 'M3' not in self.steady_states_all.columns:
            self.force_pCa =  self.steady_states_all['M2']
        else: 
            self.force_pCa = self.steady_states_all['M3'] + self.steady_states_all['M2']
        max_force = self.force_pCa.max()
        min_force = self.force_pCa.min()
        half_force = (max_force + min_force) / 2
        self.pCa_50 = np.interp(half_force, self.force_pCa.values, df.index)


        if 'M3' in self.steady_states_all.columns:
            # Note, this curve is likely going to be off because M2 + M3 issue now. 
            upper_curve = self.force_pCa + self.steady_states_all_std['M3']
            upper_half = (upper_curve.max() + upper_curve.min())/2
            self.upper_pCa_50  = np.interp(upper_half, upper_curve.values, upper_curve.index)
            

            lower_curve = self.force_pCa - self.steady_states_all_std['M3']
            lower_half = (lower_curve.max() + lower_curve.min())/2
            self.lower_pCa_50 = np.interp(lower_half, lower_curve.values, lower_curve.index)

    def twitch_tension_integral(self, reference_integral = 1):
        if self.twitches is None:
            self.get_twitch()
        integral = np.trapz(self.twitches.mean(axis  = 1), self.states_df.Time)
        return integral / reference_integral

    def time_to_peak(self):
        '''Note that this will also call some other function to set variables with regards to the twitch'''
        if self.twitches is None:
            self.get_twitch()

        smooth_twitch = self.twitches.mean(axis = 1).rolling(50, center = True).mean()    
        index_max = smooth_twitch.idxmax()
        time_to_peak = self.states_df.Time[index_max]
        self.TTP = time_to_peak

        self.relaxation_time_50()

        self.max_twitch_force = smooth_twitch.max()
        self.smoothed_twitch = smooth_twitch

        return time_to_peak

    def relaxation_time_50(self):
        if self.twitches is None:
            self.get_twitch()


        smooth_twitch = self.twitches.mean(axis = 1).rolling(50, center = True).mean()
        self.smoothed_twitch = smooth_twitch
        index_max = smooth_twitch.idxmax()
        max_peak = smooth_twitch.max()
        # Find the first time point after the peak where the force is less than 50% of the peak
        half_peak = max_peak / 2
        relaxation_index = smooth_twitch.loc[index_max:].lt(half_peak).idxmax()
        if relaxation_index == index_max: 
            self.RT50 = np.nan
            return np.nan
        else:
            time_to_peak = self.states_df.Time[index_max]
            relaxation_time = self.states_df.Time[relaxation_index]
            self.RT50 = relaxation_time - time_to_peak
            return relaxation_time - time_to_peak
    
    def calculate_nH(self):
        """
        Calculate Hill coefficient (nH) fitting the Hill equation.
        """
        # Extract parameters
        nparams, options = opt.curve_fit(normalized_force,
                                self.force_pCa.index,
                                self.force_pCa.values/np.max(self.force_pCa.values),
                                p0=[1, 6.5])
        n_H = nparams[0]
        # pCa_50 = nparams[1]
        self.n_H = n_H
        return
    





def extract_parameters(file_name, custom_sep = None, verbose = False):
    # This assumes that there is "_States_out " or "_Force_out " or "_Force_pCa_Optmz " or "_Force_pCa_Optmz_Normalized "
    # as the separator at the beginning. 
    if custom_sep != None:
        if custom_sep == '':
            param_half_string = file_name.strip(".csv")
        else:
            param_half_string = file_name.split(custom_sep)[1].strip(".csv")
    else:
        sep_options = ["_States_out ", "_Force_out ", "_Force_pCa_Optmz ", "_Force_pCa_Optmz_Normalized ", None]
        for sep in sep_options:
            if sep in file_name:
                break 
            
        # 
        param_half_string = file_name.split(sep)[1].strip(".csv")
    
    parameters = {}
    params_list = param_half_string.split(' ')
    i = 0 
    while i < len(params_list):
        try:
            parameters[params_list[i]] = float(params_list[i+1])
        except:
            try:
                print('Flipping order for parameters') if verbose else None
                parameters[params_list[i+1]] = float(params_list[i])
            except:
                print('Failed to make with values {} and {}'.format(params_list[i], params_list[i+1])) if verbose else None
                i +=1
                next 
        i += 2
    print('Read in parameters', parameters) if verbose else None
    return parameters

def plot_twitch_from_sim(input_states_structure, shading = False):
    force_cols = [col for col in input_states_structure.states_df.columns if 'M3' in col]

    plt.plot(input_states_structure.states_df.Time, input_states_structure.states_df[force_cols].mean(axis  = 1), 'C0', alpha = 1, label = 'Old')
    if shading: 
        plt.plot(input_states_structure.states_df.Time, input_states_structure.states_df[force_cols], 'C0', alpha = 0.05)

def get_twitch_from_sim(input_states_structure, mean = True, time = False):
    force_cols = [col for col in input_states_structure.states_df.columns if 'M3' in col]
    if time == False:    
        if mean:
            return input_states_structure.states_df[force_cols].mean(axis  = 1)
        else:
            return input_states_structure.states_df[force_cols]
    else:
        if mean:
            return input_states_structure.states_df.Time, input_states_structure.states_df[force_cols].mean(axis  = 1)
        else:
            return input_states_structure.states_df.Time, input_states_structure.states_df[force_cols]

# my_var = states_instance('/crucial/temp_MCMC/dATP_multiscale_modeling/MCMC_simulation_results/241004-1555_MR_640_States_out k2_plus_ref 0.002500 k3_plus 0.050000 k4_plus_ref 0.135000 kB_plus_ref 13.000000 kB_minus_ref 0.100000 kCa_plus_ref 0.090000 dATP 0.250000 k_force 0.000200 k_plus_SR_ref 16.000000 k_minus_SR_ref 15.000000.csv')

def normalized_force(pca: np.ndarray, n_H: float, pCa_50: float) -> np.ndarray:
    return 1 / (1 + 10**(n_H * (pca - pCa_50)))







def get_nH(state_structure: states_structure):
    """
    Calculate Hill coefficient (nH) fitting the Hill equation.
    """
    # Extract parameters
    nparams, options = opt.curve_fit(normalized_force,
                                state_structure.force_pCa.index,
                                state_structure.force_pCa.values/np.max(state_structure.force_pCa.values),
                                p0=[1, 6.5])
    n_H = nparams[0]
    pCa_50 = nparams[1]
    return n_H, pCa_50


#### Section for Manipulation of Lots of Simulations 

def create_results_df(parameter_filename, simulation_directory, slice_sims = False, 
                      exp_file = 'mohran_force_pCa_normal_normalized.csv', 
                      num_states = 7, suppress_output = False):
    """
    Create a DataFrame with simulation results from a given parameter file and directory.
    
    Args:
        parameter_filename (str): Path to the parameter file.
        simulation_directory (str): Directory containing simulation files (but not the raw_data directory).
        slice_sims (bool): Whether to rename the columns, and move the simulation column o the force_sim and twitch_sim, assumes all force sims first, then all twitch sims
        exp_file (str): A file that is used to read in the results and matching pCa values
        num_states (int): Number of states that are in the model (7 expected for the super slow state simulations)
    Returns:
        pd.DataFrame: DataFrame containing simulation results.
    """
    # Read the parameter file
    params = pd.read_csv(parameter_filename, comment='#')
    params['simulation'] = None
    # Initialize an empty list to store results
    results = []
    
    # Iterate through each row in the parameters DataFrame
    for ind in tqdm(params.index, disable=suppress_output):
        # Construct the filename based on the index
        filename = os.path.join(simulation_directory, f'Rep_{ind}States_out.csv')

        # Check if the file exists
        if os.path.exists(filename):
            # Load the states structure from the file
            simulation_instance = states_structure(filename, skip_params=True, exp_file=exp_file, num_states=num_states)
            # Append the result to the list
            params.loc[ind, 'simulation']  = simulation_instance
        else:
            print(f"File {filename} does not exist.")

    if slice_sims:
        n_sims = len(params)
        force_sims = params.simulation.iloc[0:n_sims//2].values
        twitch_sims = params.simulation.iloc[n_sims//2:n_sims].values
        params_first_half = params.iloc[0:n_sims//2, 0:-1].reset_index(drop=True).copy()
        params_first_half['force_sim'] = force_sims
        params_first_half['twitch_sim'] = twitch_sims
        return params_first_half
    else: 
        return params

def add_data(original_df, new_path, parameter_file = '/parameter_input_temp.csv'):
    '''
    Make sure that the new_path points to the new larger directory 
    The parmaeter file is just the file name. It's joined to the new path 
    It's assumed that the data is in raw_data
    Does not work with slicing sims. 
    '''
    temp_df = pd.read_csv(new_path + parameter_file, comment = '#')
    temp_df

    for index in temp_df.index: 
        filename = new_path+"/raw_data/"+f'Rep_{index}States_out.csv'
        simulation_instance = states_structure(filename, 
                            skip_params = True, 
                            exp_file = 'mohran_force_pCa_normal_normalized.csv')

        temp_df.loc[index, 'simulation'] = simulation_instance

    return pd.concat([original_df, temp_df], ignore_index=True)


def average_twitch(twitch: np.ndarray, time: np.ndarray, period: float, skip_periods: int = 0) -> tuple[pd.Series, pd.Series]:
    full_twitch = pd.DataFrame({'Time': time, 'Twitch': twitch})
    num_twitches = int(np.ceil(full_twitch.Time.max()/period))


    split_df = pd.DataFrame({'Time': time[time<period]})
    twitch_column_names = []
    for i in range(skip_periods, num_twitches):
        start_time = i * period
        end_time = (i + 1) * period
        mask = (full_twitch['Time'] >= start_time) & (full_twitch['Time'] < end_time)
        split_df[f'Twitch_{i}'] = np.nan  # Initialize with NaN
        twitch_column_names.append(f'Twitch_{i}')
        mask_len = sum(mask)
        split_df.loc[0:mask_len-1, f'Twitch_{i}'] = full_twitch.loc[mask, 'Twitch'].values
    return split_df['Time'], split_df[twitch_column_names].mean(axis = 1)


def apply_parameteter_set(input_df, corresponding_parameter = 'k1_plus_ref_baseline', inplace = True):
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    if not inplace:
        input_df = input_df.copy()
    input_df['parameter_set'] = None
    param_dict = {}
    counter = 0
    for i in input_df.index:
        if input_df.loc[i, corresponding_parameter] not in param_dict:
            param_dict[input_df.loc[i, corresponding_parameter]] = alphabet[counter]
            input_df.loc[i, 'parameter_set'] = alphabet[counter]
            counter += 1
        else:
            input_df.loc[i, 'parameter_set'] = param_dict[input_df.loc[i, corresponding_parameter]]
    return input_df

def end_state_analysis(param_and_sim_df: pd.DataFrame, 
                       separate_by: str, 
                       x_axis_var: str = 'k2_plus_baseline', 
                       ) -> list: 
    '''
    separate_by = "pCa" or "percent_drug"

    '''
    return 0
    list_of_pCa = param_and_sim_df.simulation[0].steady_states_all.index.tolist()
    
    # This code still needs a lot of work. 
    # for pCa_ind in param_and_sim_df.simulation[0].steady_states_all.index:
    #     percent_drug = 0
    #     end_states = []


    #     filtered_df = param_and_sim_df.query('protocol == 1').sort_values(by='k2_plus_baseline')


    #     for sim in filtered_df.simulation:
    #         # sim.steady_states_all.plot(color = colors, ls = lines[i], ax = ax)
    #         end_states.append(sim.steady_states_all.loc[pCa_ind].to_frame().T)

    #     end_state_df = pd.concat(end_states)

    #     end_state_df['k2_plus_baseline'] = filtered_df.k2_plus_baseline.values

    #     melted = pd.melt(end_state_df, id_vars = 'k2_plus_baseline', var_name = 'State', value_name = 'Fraction')

    #     melted


    #     i = 0
    #     for state in melted.State.unique():
    #         title_str = f'k2_plus_drug varied @ afic = {percent_drug} uM'
    #         # plt.figure(figsize=(8,6))
    #         col, row = divmod(i, 4)
    #         ax[row, col].set_title(f'State {state} fractions')
    #         ax[row, col].set_xlabel('k2_plus_drug')
    #         ax[row, col].set_ylabel('Fraction')
    #         # ax[row, col].plot(melted.query('State == @state').k2_plus_baseline, melted.query('State == @state').Fraction, marker='o')
    #         # ax[row, col].legend(None)
    #         i +=1 
    #         if state == 'C':
    #             sns.lineplot(data = melted.query('State == @state'), x = 'k2_plus_baseline', y = 'Fraction', marker='o', ax = ax[row, col], label = f'pCa {pCa_ind:.1f}')
    #             ax[row, col].legend(loc = (1.05,0))
    #         else:
    #             sns.lineplot(data = melted.query('State == @state'), x = 'k2_plus_baseline', y = 'Fraction', marker='o', ax = ax[row, col])

    #             # ax[row, col].legend()# [f'pCa {pca:.1f}' for pca in list_of_pCa], loc = (1.05,-0.2))
    #             # ax[row, col].set_yscale('log')
    # # plt.scatter(melted.query('State == "M3"').k2_plus_drug, melted.query('State == "M3"').Fraction, color = 'red')
    # # Hide the plot in the last subplot (bottom right)
    # ax[3,1].axis('off')
    # # plt.suptitle(f'Changes at all aficamten concentrations at various calcium')
    # plt.tight_layout()
    # plt.show()