#!/bin/bash

# Activate the conda environment
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/marcus/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/marcus/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/marcus/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/marcus/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
conda activate cuda_11_7

# Clean the build
make clean

# Build the project
make all

# Run the Python script
python sweep.py
