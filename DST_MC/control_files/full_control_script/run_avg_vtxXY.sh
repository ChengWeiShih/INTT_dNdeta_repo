#!/bin/bash

#todo : the script should be copied to the corresponding control folder

# note : full directory
control_folder_name=`pwd`
file_folder_name=""

# Add a condition to execute only part of the script
if [[ $1 == "run_condor" ]]; then
    echo "the script is going to run_condor"
    echo "the control_folder_name: $control_folder_name"
    
    
elif [[ $1 == "run_merged_result" ]]; then
    echo echo "the script is going to run_merged_result"

else 
    echo "Usage: $0 {run_condor|run_merged_result}"
fi

echo "This is the end of the script"
