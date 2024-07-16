#!/bin/bash

# Loop over all subdirectories
for item in {1..15}
do
    dir=cluster_$item
    # Check if directory exists (in case of no subdirectories)
    [ -d "$dir" ] || continue

    # Change to the subdirectory
    echo Changing to directory "$dir"
    cd "$dir" || continue

    # Run your simulation commands
    # Uncomment the run_apbs_inputs line if the electrostatic grids still need to be calculated
    # run_apbs_inputs input.xml
    bd_top input.xml
    nam_simulation myosin_actin_simulation.xml

    # Return to the original directory
    cd - || exit
done

