import matplotlib.pyplot as plt
import helper_functions as hf 
import argparse
import numpy as np 


# What would actually be more useful, basically to call the create_results_df function 
# Then plot the different force pCa curves 

# Need to have a flag for slice, force pCa, or twitch 
# By defualt maybe it should just sweep though, all the rows in the simulation run. 
# No normaliztion is probably fine, or normalization to the very first row, which is the new baseline
#  

def read_command_line_arguments():
    parser = argparse.ArgumentParser(description='Read in the command line arguments. The -c flag for the csv file name is required, and other flags are recommended. ')
    parser.add_argument('-t', type=str, default = 'slice', 
                        help='Type of simulation df to read in and analyze. Options are "slice","force", or "twitch"')
    parser.add_argument('-n', type=int, default = 0,
                        help='Normalization index. Default is 0. Enter -1 for no normalization.')
    parser.add_argument('-c', type=str, help='CSV filename of the simulation results df to read in.',
                        required = True)
    parser.add_argument('-d', type=str, default = 'raw_data',
                        help='Data directory to read in the simulation results df from. Default is "raw_data".')
    parser.add_argument('-nt', type=int, default = 0, 
                        help='Number of twitches run in the simulation. ')
    args = parser.parse_args()
    return args



def plot_force_pCa(results_df, normalization_index):
    unique_dict = get_unique_parameters(results_df)
    if 'force_sim' in results_df.columns:
        sim_name = 'force_sim'
    elif 'simulation' in results_df.columns:
        sim_name = 'simulation'
    else:
        print('No force_pCa data found in the results df.')
        return 
    
    # Begin plotting process 
    plt.figure(figsize = (8,6))

    if normalization_index == -1:
        scaling_factor = 1
    else: 
        scaling_factor = results_df.loc[normalization_index, sim_name].force_pCa.max()
    for index, row in results_df.iterrows():
        additive_label = ''
        for key in unique_dict.keys():
            # Get the first value from all the keys 
            unique_value = row[key]
            additive_label += f'{key}: {unique_value:.2e}, '
        plt.plot(row[sim_name].force_pCa / scaling_factor, label = additive_label)

    
    plt.gca().invert_xaxis()
    plt.xlabel('pCa')
    plt.ylabel('Normalized Force')
    plt.title('Force-pCa Curves')
    plt.legend(loc = (1.04, 0))
    plt.savefig('force_pCa_curves.png', bbox_inches='tight')
    plt.show()
    return

def get_unique_parameters(results_df):
    # Find which parameters are changing through the course of the simulation and return a dict 
    parameter_columns = [col for col in results_df.columns if col not in ['simulation', 'force_sim', 'twitch_sim']]
    unique_parameters = {}
    for col in parameter_columns:
        unique_values = results_df[col].unique()
        if len(unique_values) > 1:
            unique_parameters[col] = unique_values
    return unique_parameters

def plot_twitches(results_df, nt, normalization_index):
    # Check the name for the column 
    if 'twitch_sim' in results_df.columns:
        sim_name = 'twitch_sim'
    elif 'simulation' in results_df.columns:
        sim_name = 'simulation'
    else:
        print('No twitch data found in the results df.')
        return
    # Then get the unique parameters dict
    unique_dict = get_unique_parameters(results_df)
    # Get the scaling twitch: 
    if normalization_index == -1:
        scaling = 1
    else:
        sim = results_df.loc[normalization_index, sim_name]
        scaling = (sim.get_twitch().mean(axis = 1).values[0:-1].reshape(nt, -1)).mean(axis = 0).max()

    # Create figure and start working through it 
    fig, ax = plt.subplots(figsize = (8,6))
    for index, row in results_df.iterrows():
        sim = row[sim_name]
        mean_value = (sim.get_twitch().mean(axis = 1).values[0:-1].reshape(nt, -1)).mean(axis = 0)
        std_value = (sim.get_twitch().mean(axis = 1).values[0:-1].reshape(nt, -1)).std(axis = 0)
        mean_value /= scaling
        std_value /= scaling

        n_steps = len(mean_value)

        additive_label = ''
        for key in unique_dict.keys():
            # Get the first value from all the keys 
            unique_value = row[key]
            additive_label += f'{key}: {unique_value:.2e}, '

        plt.plot(np.linspace(0,1,n_steps), mean_value, label = additive_label)
        plt.fill_between(np.linspace(0,1,n_steps), mean_value - std_value, mean_value + std_value, alpha = 0.2)
        
    plt.legend(loc = (1.04, 0))
    plt.xlabel('Time (s)')
    plt.ylabel('Normalized Force')
    plt.title('Twitch Force Curves')
    plt.savefig('twitch_force_curves.png', bbox_inches='tight')
    plt.show()
    return



def main():
    # Read the command line arguments
    print("Running program... ")
    args = read_command_line_arguments()
    type_of_df = args.t
    normalization_index = args.n
    csv_filename = args.c
    data_dir = args.d
    num_twitches = args.nt

    if (type_of_df == 'twitch' or type_of_df == 'slice') and num_twitches <= 0:
        print('Number of twitches must be specified correctly for twitch or slice simulations. Please enter a valid number of twitches with the -nt flag.')
        return -1

    if type_of_df not in ['slice', 'force', 'twitch']:
        print('Invalid type of df. Options are "slice","force", or "twitch".')
        return -1
    results_df = hf.create_results_df(csv_filename, data_dir, slice_sims = (type_of_df == 'slice'))
    # Read the data from the files
    # sim_df = hf.create_results_df(sim_type, data_dir, csv_filename)
    if type_of_df == 'force' or type_of_df == 'slice':
        plot_force_pCa(results_df, normalization_index)

    if type_of_df == 'twitch' or type_of_df == 'slice':
        plot_twitches(results_df, num_twitches, normalization_index)

    print('Unique parameters changing through the course of the simulation:', get_unique_parameters(results_df))
    
    return 0



if __name__ == '__main__':
    main()  # Run the main function
