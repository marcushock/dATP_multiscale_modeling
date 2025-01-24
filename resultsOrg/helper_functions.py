import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

class states_structure:
    def __init__(self, input_filename:'CSV', OFFState = True) -> None: # type: ignore
        print('The filename is: ', input_filename)
        self.file_name = input_filename
        if OFFState:
            self.state_list = ['M2','M1','C','B','OFF']
        else:
            self.state_list = ['M2','M1','C','B']
        self.read_states()
        self.extract_parameters()
        self.states_steadystate()

        return 
    
    def read_states(self):
        pd.read_csv(self.file_name, header = None)
        columns_list = ['Time'] # unit unknown 
        for i in np.arange(7,3.9,-0.1):
            for j in self.state_list:
                new_item = ('{} {}'.format(j,round(i,1)))
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
        print('Read in parameters', self.parameters)
        return 

    def get_twitch(self):
        '''
        This function returns the temporal force states from the simulation and does not 
        average them across the replciaets. 
        '''
        force_cols = [col for col in self.states_df.columns if 'M2' in col]
        self.twitch_states = self.states_df[force_cols]
        return self.states_df[force_cols]

    def states_steadystate(self, end_time_amount= 500): 
        temp_series = self.states_df[self.states_df.Time>end_time_amount].mean()
        std_series = self.states_df[self.states_df.Time>end_time_amount].std()
        # Define new steady state DF with each column a different state
        # And each row a different pCa (7,4,0.1 steap)
        # Also creating a place for the pCa first, then renaming 
        df = pd.DataFrame(np.zeros((31,6)), columns = ['pCa','M2','M1','C','B','OFF'])
        df.pCa = np.arange(7,3.9,-0.1)
        for i in range(len(df.pCa)):
            df.pCa[i] = round(df.pCa[i],1)
        df = df.set_index('pCa')

        df_std = pd.DataFrame(np.zeros((31,6)), columns = ['pCa','M2','M1','C','B','OFF'])
        df_std.pCa = np.arange(7,3.9,-0.1)
        for i in range(len(df_std.pCa)):
            df_std.pCa[i] = round(df_std.pCa[i],1)
        df_std = df_std.set_index('pCa')

        # Annoying rounding necessary for the python decimal storage necessity
        pca = 7
        while pca >= 4:
            for state in self.state_list:
                value = temp_series.get(f'{state} {pca:.1f}')
                df.at[round(pca,2), state] = value 

                value = std_series.get(f'{state} {pca:.1f}')
                df_std.at[round(pca,2), state] = value 
            pca -= 0.1
        
        max_force = df['M2'].max()
        min_force = df['M2'].min()
        half_force = (max_force + min_force) / 2
        self.pCa_50 = np.interp(half_force, df['M2'].values, df.index)

        self.steady_states_all  = df
        self.steady_states_all_std  = df_std
        self.force_pCa = self.steady_states_all['M2']
        df_std['M2'].index, df_std['M2'].values

        upper_curve = self.force_pCa + self.steady_states_all_std['M2']
        upper_half = (upper_curve.max() + upper_curve.min())/2
        self.upper_pCa_50  = np.interp(upper_half, upper_curve.values, upper_curve.index)
        

        lower_curve = self.force_pCa - self.steady_states_all_std['M2']
        lower_half = (lower_curve.max() + lower_curve.min())/2
        self.lower_pCa_50 = np.interp(lower_half, lower_curve.values, lower_curve.index)





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
    force_cols = [col for col in input_states_structure.states_df.columns if 'M2' in col]

    plt.plot(input_states_structure.states_df.Time, input_states_structure.states_df[force_cols].mean(axis  = 1), 'C0', alpha = 1, label = 'Old')
    if shading: 
        plt.plot(input_states_structure.states_df.Time, input_states_structure.states_df[force_cols], 'C0', alpha = 0.05)

def get_twitch_from_sim(input_states_structure, mean = True, time = False):
    force_cols = [col for col in input_states_structure.states_df.columns if 'M2' in col]
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

    