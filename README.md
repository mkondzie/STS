# STS - Silicon Tracking System
Raw data from calibration of the STS modules, analysis of the STS data and automation of analysis-related tasks and activities for the STS.
# How to use
## Running tests 
### Calibration
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
### Burn-in
Before launching a thermal stress test, the hardware communication script (feb_test.py or feb_test2.py) needs to be modified based on the module type (PA or PB). 
In order to modify the `feb_test.py` script for BURN-IN setup 1 in case of a PA-type module, run the following:
```
 ./module_PA.sh
```
Similarly, to modify the `feb_test.py` script for BURN-IN setup 1 in case of a PB-type module, run the following:
```
 ./module_PB.sh
```
When working with BURN-IN setup 2, to adjust the script to a PA-type module, please run:
```
 ./module_PA2.sh
```
Finally, in case of BURN-IN setup 2, to adjust the script to a PB-type module, run:
```
 ./module_PB2.sh
```
After initiating the appropriate script, configure the module-specific details by editing the `test_parameters.json` file for BURN-IN setup 1 or `test_parameters2.json` for BURN-IN setup 2, which contains config data for the burn-in process. Finally, run the `run_exp.py` or `run_exp2.py` script to begin the burn-in test and collect the necessary data.
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
### Calibration
The `check_trim_calib.sh` and `conn_check.sh` scripts are used (see the *Running tests* section) to modify the python script, for calibration and connectivity test, accordingly.
### Burn-in
As described in the *Running tests* section, the `module_PA.sh`, `module_PB.sh`, `module_PA2.sh` and `module_PB2.sh` scripts help during thermal stress test launch.
### ASIC sorting
The `find_ASICs.sh` script checks if all ASIC files are present for each polarity. In case some files are missing, it modifies the `plot_1024.C` macro accordingly. Thanks to this scirpt, even the partial results can be visualized without encountering errors.
The `count_select_sort_files.sh` script, as the name suggests, is responsible for counting, selecting the most recent files and sorting them in a specific order depending on the module type PA or PB. Similarly, `count_select_sort_root_files.sh` script handles the same tasks but for ROOT files. These scripts are employed during `analyze_pscan.sh` execution, to properly list the relevant files for further analysis.
### p-scan analysis
The `analyze_pscan.sh` script comprises of several repetitive steps that need to be performed in order to visualize the most recent calibration results. Several techniques have been used to ensure smooth macros execution even on imperfect data.
### Analysis cleanup 
In case you wish to start over with the analysis, there is a fast way to remove all the files created during the execution of the analysis macros. Run the `clean.sh` script to clean the module files from generated files, keeping all of the raw data files untouched.
## Python scripts
### Calibration
The `tester_febs_setup<number_of_the_setup>_arr.py` script serves as hardware communication tool. It is used to connect EMUs, find FEBs and sync ASICs. The script runs test sequence (calibration of the module, setting the calibration trim and connectivity check) and writes data to files. 
### Burn-in
The `feb_test.py` or the `feb_test2.py` script (depending on the setup) is used to manage the thermal stress testing process. This script controls a climatic chamber, a chiller, power supplies, as well as synchronizes the ASICs. 
## Module files
For each tested module a separate directory named after module ID is created. The directory consists of 3 main folders, a data file and log files. 
### p-scan files
The pscan_files directory contains one .txt file for each ASIC. For each out of chosen discriminators, signals are recorded for every channel of an ASIC. Signals originate from a pulse generator injecting the channels with fixed-amplitude pulses. Based on this procedure, an individual threshold is assigned to each discriminator.
### trim files
The trim_files directory contains threshold correction data for each discriminator of each channel of each ASIC of a given polarity.
### connectivity check files
The conn_check_files directory contains .txt files, one for each ASIC. The files include ADC values for each channel and the corresponding amplitude. 
## Analysis macros
### trim adc
The trim_adc.hxx and trim_adc.cxx files introduce the trim_adc class, which is responsible for reading raw data (.txt files) and converting it into ROOT files for further analysis.
### execution
The execution.C macro acts as a proxy for reading files and utilizes the trim_adc class.
### plot 1024
The plot_1024.C macro is used to plot the calibration data for all ASICs, sorted including the position of each ASIC on the FEB. Noise levels predicted based on sensor length, microcable length and intrinsic ASIC noise are adjusted for each module, based on the module ID. ADC gain as well as ADC threshold for each out of 1024 channels is visualized. Broken channels can be detected with the distinction of no analog response (NAR), broken bond at the ASIC (ASIC) and broken bond at the sensor (SENS).
### analysis conn check
In order to check connectivity when high voltage is turned off, the analysis_conn_check.C macro is used as an alternative method for detecting broken channels. The running median algorithm is implemented with optimized neighborhood size and a dynamic threshold parameter effective across a vast majority of modules.
# References 
[Nuclear Instruments and Methods in Physics Research A 1058 (2024) 168813](https://doi.org/10.1016/j.nima.2023.168813)
# Future development ideas
* rewrite from scratch
  + python serial communication scripts 
  + ROOT analysis macros
* organize the codebase using OOP principles
* synchronize data taking and analysis to plot results dynamically
