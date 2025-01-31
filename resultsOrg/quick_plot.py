import matplotlib.pyplot as plt
import helper_functions as hf 
import argparse

#### Example usage: python quick_plot.py -T 'twitch_states.csv' -F 'force_pCa.csv'

# Create an command line argument reader that reads two optional arguments with the flags -T or -F 
# and returns the values of the arguments as a tuple.
def read_command_line_arguments():
    parser = argparse.ArgumentParser(description='Read two optional arguments with the flags -T or -F')
    parser.add_argument('-T', type=str, help='Twitch States File Name')
    parser.add_argument('-F', type=str, help='Force pCa States File Name')
    args = parser.parse_args()
    return args.T, args.F


def plot_twitch(input_data):
    force_cols = [col for col in input_data.states_df.columns if 'M2' in col]
    plt.plot(input_data.states_df.Time, input_data.states_df[force_cols], 'C0', alpha = 0.05)
    plt.plot(input_data.states_df.Time, input_data.states_df[force_cols].mean(axis  = 1), 'C0', alpha = 1, label = 'New')
    plt.show()
    return

def plot_force_pCa(input_data):
    plt.plot(input_data.force_pCa)
    # Flip the x-axis 
    plt.gca().invert_xaxis()
    plt.show()
    return 
    

def main():
    # Read the command line arguments
    twitch_file, force_file = read_command_line_arguments()
    # Read the data from the files
    if twitch_file is not None:
        data = hf.states_structure(twitch_file)
        plot_twitch(data)
    if force_file is not None:
        data = hf.states_structure(force_file)
        plot_force_pCa(data)
    
    return 0

if __name__ == '__main__':
    main()  # Run the main function