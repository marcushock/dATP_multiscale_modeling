import numpy as np 
import pandas as pd

class states_instance:
    def __init__(self, input_filename:'CSV') -> None:
        print('The filename is: ', input_filename)
        self.file_name = input_filename
        self.read_states()
        self.extract_parameters()
        return 
    
    def read_states(self):
        pd.read_csv(self.file_name, header = None)
        columns_list = ['Time'] # unit unknown 
        for i in np.arange(7,3.9,-0.1):
            for j in ['M2','M1','C','B','OFF']:
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
                except:
                    print('Failed to make dict')
                    return 
            i += 2
        print('Read in dict', self.parameters)
        return 

states_instance('/crucial/temp_MCMC/dATP_multiscale_modeling/MCMC_simulation_results/241003-1020_States_out k2_plus_ref 0.004780 k3_plus 0.050000 k4_plus_ref 0.135000 kB_plus_ref 13.000000 kB_minus_ref 0.100000 kCa_plus_ref 0.090000 dATP 0.050000 k_force 0.200000 k_plus_SR_ref 16.000000 k_minus_SR_ref 15.000000.csv')
