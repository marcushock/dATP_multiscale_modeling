import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize as opt

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