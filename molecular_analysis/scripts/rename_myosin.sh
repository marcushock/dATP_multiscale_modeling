#!/bin/bash

### Auxillary script used to quickly rename myosin pqr files. 
### Looks for files matching 'myosin_*_*.pqr' and removes the middle 
### integer which generally corresponds to the metastable state, and just returns
### the sample ID (at the end of the file name)

### Edited by MTH on 07/12/2024

# Iterate over all files matching the pattern 'myosin_*_*.pqr'
for file in metastable_*_*.pqr; do
    # Extract the last integer part (X) from the file name
    integer_part=$(echo "$file" | sed 's/metastable_[0-9]*_\([0-9]*\).pqr/\1/')

    # Construct the new file name
    new_file_name="myosin_${integer_part}.pqr"

    # Rename the file
    mv "$file" "$new_file_name"
done

