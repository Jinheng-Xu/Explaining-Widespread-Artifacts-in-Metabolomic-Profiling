Copyright reserved by Jinheng Xu, Stanford University

1. Demo data:
NR.csv: background peaks you want to remove.
Metabolites.csv: a list of different metabolite names.
MS-data.csv: a mass spectrum, in the form of list of m/z vs. intensity.

When using your own data, please keep the same format in each csv file: for example, in MS-data.csv, the first column should be m/z and the second column should be intensity.

2. To run the Python code, 
a. Please install following packages:
"requests" "rdkit" "pandas" "molmass".
b. Please follow the instructions from Line 7 to Line 25 in "code.py" to set paths and parameters.

3. When running the Python code,
a. Please connect to the internet when running the program.
b. The progress and status will display in the terminal.
c. Generated species list (based on metabolite names input), Generated standard peaks list and Results will be saved in the designated folder as csv files.