# STS
Raw data from calibration of the STS modules, analysis of the STS data and automation of analysis-related tasks and activities for the STS.
# How to use
## Running tests 
In order to launch the calibration (LV ON, HW ON), modify the python script accordinly by running a shell script:
```
./check_trim_calib.sh <number_of_the_setup>
```
to call it you need to provide the number of the setup as the argument. For instance, if you want to modify the script for setup 2, you need to call:
`` ./check_trim_calib.sh 2
``
In order to launch the connectivity check (LV ON, HV OFF), modify the python script accordinly by running a shell script:
```
 ./conn_check.sh <number_of_the_setup>
```

## Analysis
In order to analyze data for one module (example data included), run the following command:
```
$ ./analyze_pscan.sh M4DL0B1001611B2/
```
# Description
## Bash scripts
## Python scripts
## Module files
### p-scan files
### trim files
### connectivity check files
## Analysis macros

# Future development ideas
* rewrite from scratch
  + python serial communication scripts 
  + ROOT analysis macros
* divide the codebase using OOP principles
* synchronize data taking and analysis to plot results dynamically
