#!/bin/bash

filename="feb_test.py"

# Use sed to replace the line containing #active_downlinks = [1
sed -i '/#active_downlinks = \[1/s/.*/    active_downlinks = \[1,2\] #For modules PB/' "$filename"

# Comment out the line containing active_downlinks = [0,3] 
sed -i '/active_downlinks = \[0/s/.*/    #active_downlinks = \[0,3\]  #For modules PA/' "$filename"

echo "Lines have been modified for a PB module in $filename"

