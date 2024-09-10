#!/bin/bash

find ../../module_files/ -type f \( -name '*.root' -o -name '*.pdf' -o -name '*.exe' -o -name '*.C' -o -name '*.cxx' -o -name '*.hxx' -o -name 'a.txt' -o -name 'plot.txt' -o -name 'conn_check_parameters.txt'  -o -name 'conn_check_summary.txt' \) -exec rm -v {} \;
find ../../ module_files/ -type f -name 'module_test_*' -exec rm -v {} \;
