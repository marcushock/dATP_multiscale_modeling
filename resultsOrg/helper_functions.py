import numpy as np 
import pandas as pd

class states_instance:
    def __init__(self, input_filename:'CSV') -> None: # type: ignore
        print('The filename is: ', input_filename)
        self.file_name = input_filename
        self.read_states()
        self.extract_parameters()
        self.states_steadystate()
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
                    self.parameters[params_list[i+1]] = float(params_list[i])
                except:
                    print('Failed to make with values {} and {}'.format(params_list[i], params_list[i+1]))
                    i +=1
                    next 
            i += 2
        print('Read in parameters', self.parameters)
        return 


    def states_steadystate(self, end_time_amount= 500): 
        temp_series = self.states_df[self.states_df.Time>end_time_amount].mean()
        # print(temp_series)
        # Define new steady state DF with each column a different state
        # And each row a different pCa (7,4,0.1 steap)
        # Also creating a place for the pCa first, then renaming 
        df = pd.DataFrame(np.zeros((31,6)), columns = ['pCa','M2','M1','C','B','OFF'])
        df.pCa = np.arange(7,3.9,-0.1)
        for i in range(len(df.pCa)):
            df.pCa[i] = round(df.pCa[i],1)
        df = df.set_index('pCa')

        # Annoying rounding necessary for the python decimal storage necessity
        pca = 7
        while pca >= 4:
            for state in ['M2','M1','C','B','OFF']:
                value = temp_series.get(f'{state} {pca:.1f}')
                print(value)
                df.at[round(pca,2), state] = value 
            pca -= 0.1

        self.steady_states_all  = df
        self.force_pCa = self.steady_states_all['M2']




my_var = states_instance('/crucial/temp_MCMC/dATP_multiscale_modeling/MCMC_simulation_results/241004-1555_MR_640_States_out k2_plus_ref 0.002500 k3_plus 0.050000 k4_plus_ref 0.135000 kB_plus_ref 13.000000 kB_minus_ref 0.100000 kCa_plus_ref 0.090000 dATP 0.250000 k_force 0.000200 k_plus_SR_ref 16.000000 k_minus_SR_ref 15.000000.csv')

