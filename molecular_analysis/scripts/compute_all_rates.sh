#!/bin/bash

### This is script will run the rates_of_distances BD function (must be in path)
### for all the files named "output_distances.xml" in any current or sub directories
### using the find command 

### Edited by MTH 07/12/2024


# Save the current working directory
start_dir=$(pwd)

compute_rates_command="rates_of_distances -plain -res results.xml -dist output_distances.xml  > output_rates.txt"

# Use 'find' to locate all files named 'output_distances.xml'
# and execute a command within their directory
echo "Searching for 'output_distances.xml' in the current directory and its subdirectories..."

find "$start_dir" -name 'output_distances.xml' -exec sh -c 'cd "$(dirname {})" && echo "Executing compute_rates in $(dirname {})" && rates_of_distances -plain -res results.xml -dist output_distances.xml  > output_rates.txt' \;

# Return to the original working directory
cd "$start_dir"
