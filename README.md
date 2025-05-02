# tooth-evo-devo

Code steps:
1. python3 lhs.py 19000     # This will generate 19000 input .txt files using Latin Hypercube Sampling of the variable genetic and cellular parameters in ToothMaker
2. cp runt.e ./19000_input_files     #This will copy runt.e into the subdirectory with the input .txt files
3. python3 make_off_multi.py ./19000_input_files ./19000_input_files     #this makes a "multirun" code to run to generate 19000 teeth in the cluster 



#check_height
This script processes all .off files in the current directory, checks the height (difference between the maximum and minimum z-coordinates) of each file, and classifies them as either "tall enough" or "too flat" based on a threshold of 0.5. The results are written to two separate output files: 18985_tall_enough.txt and 18985_too_flat.txt. If the file is not in the correct OFF or COFF format, an error is raised.

