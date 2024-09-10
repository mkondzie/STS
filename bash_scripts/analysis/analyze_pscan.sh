#!/bin/bash

# Check if module_id is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <module_id>"
  exit 1
fi

module_id=$(basename "$1")


# Navigate to the module's pscan_files directory
cd "$module_id/pscan_files/" || exit

# Copy the necessary files from the source directory
cp ../../../analysis_macros/execution.C .
cp ../../../analysis_macros/trim_adc.cxx .
cp ../../../analysis_macros/trim_adc.hxx .
cp ../../../analysis_macros/plot_1024.C .

# List the relevant pscan files and save to a.txt
cp ../ASIC_sorting/find_ASICs.sh .
cp ../ASIC_sorting/count_select_sort_files.sh .
./count_select_sort_files.sh

# Run ROOT commands
root -l <<EOF
.L trim_adc.cxx
.L execution.C
execution()
.q
EOF

# List the generated ROOT files and save to plot.txt
cp ../ASIC_sorting/count_select_sort_root_files.sh .
./count_select_sort_root_files.sh

# Run the plotting macro
root -l <<EOF
.x plot_1024.C
.q
EOF

# Find and open the generated PDF file 
pdf_file="module_test_${module_id}.pdf"

if [ -n "$pdf_file" ] && [ -s "$pdf_file" ]; then
  # File exists and is not empty, open it with the default PDF viewer
  xdg-open "$pdf_file"
else
  # No valid PDF file found
  echo "No PDF file named $pdf_file found or the file is empty."
fi
