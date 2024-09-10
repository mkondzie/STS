#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <NUMBER-OF-THE-SETUP>"
    exit 1
fi

setup_number=$1
filename="tester_febs_setup${setup_number}_arr.py"

# Use sed to replace the line containing test_list_check = 
sed -i '/test_list_check = /s/.*/test_list_check = ["set_trim_calib", "check_trim", "#turn_hv_off", "#conn_check"]/' "$filename"

# Uncomment the line containing test_list.extend(test_list_calib)
sed -i "/#test_list.extend(test_list_calib)/s/.*/test_list.extend(test_list_calib)/" "$filename"

echo "Lines have been modified for calibration in $filename"
