#!/bin/bash

### Renames any pdb files to a pqr file format. 
### Useful after extracting pqr files from cpptraj that have the PDB extension

# Iterate over all .pdb files in the current directory
for file in *.pdb; do
    # Check if the file exists before attempting to rename
    if [ -e "$file" ]; then
        # Use the 'mv' command to rename the file
        mv "$file" "${file%.pdb}.pqr"
    fi
done

