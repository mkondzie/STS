# STS - Silicon Tracking System
Raw data from calibration of the STS modules, analysis of the STS data and automation of analysis-related tasks and activities for the STS.
# How to use
## Running tests 
In order to launch the calibration (LV ON, HW ON), modify the `tester_febs_setup<number_of_the_setup>_arr.py` python script accordinly by running a shell script:
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
In principle, all provided analysis macros require [ROOT to be installed](https://root.cern/install/).

In order to analyze data for one module (sample data included), run the following command:
```
$ ./analyze_pscan.sh <module_ID>
```
For instance, if you wish to run analysis for module M5UL5B0010180A2, call:
``
$ ./analyze_pscan.sh M5UL5B0010180A2
``
# Description
## Bash scripts
## Python scripts
## calibration
The `tester_febs_setup<number_of_the_setup>_arr.py` script serves as hardware communication tool. It is used to connect EMUs, find FEBs and sync ASICs. The script runs test sequence (calibration of the module, setting the calibration trim and connectivity check) and writes data to files. 
## burn-in
tester\_setup<number\_of\_the\_setup>-THE-SETUP\_arr.py 
## Module files
For each tested module a separate directory named after module ID is created. The directory consists of 3 main folders, a data file and log files. 
### p-scan files
The pscan_files directory contains one txt file for each ASIC. For each out of chosen discriminators, signals are recorded for every channel of an ASIC. Signals originate from a pulse generator injecting the channels with fixed-amplitude pulses. Based on this procedure, an individual threshold is assigned to each discriminator.
### trim files
The trim_files directory contains threshold correction data for each discriminator of each channel of each ASIC of a given polarity.
### connectivity check files
The conn_check_files directory contains txt files, one for each ASIC. The files include ADC values for each channel and the corresponding amplitude. 
## Analysis macros
### trim adc
The trim_adc.hxx and trim_adc.cxx files introduce the trim_adc class, which is responsible for reading raw data (txt files) and converting it into ROOT files for further analysis.
### execution
The execution.C macro acts as a proxy for reading files and utilizes the trim_adc class.
### plot 1024
The plot_1024.C macro is used to plot the calibration data for all ASICs, sorted including the position of each ASIC on the FEB. Noise levels predicted based on sensor length, microcable length and intrinsic ASIC noise are adjusted for each module, based on the module ID. ADC gain as well as ADC threshold analysis is visualized. Broken channels can be detected with the distinction of no analog response (NAR), broken bond at the ASIC (ASIC) and broken bond at the sensor (SENS).
### analysis conn check
In order to check connectivity when high voltage is turned off, the analysis_conn_check.C macro is used as an alternative method for detecting broken channels.
# References 
[Nuclear Instruments and Methods in Physics Research A 1058 (2024) 168813](https://doi.org/10.1016/j.nima.2023.168813)
# Future development ideas
* rewrite from scratch
  + python serial communication scripts 
  + ROOT analysis macros
* organize the codebase using OOP principles
* synchronize data taking and analysis to plot results dynamically
