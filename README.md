# tooth-evo-devo

#check_height
This script processes all .off files in the current directory, checks the height (difference between the maximum and minimum z-coordinates) of each file, and classifies them as either "tall enough" or "too flat" based on a threshold of 0.5. The results are written to two separate output files: 18985_tall_enough.txt and 18985_too_flat.txt. If the file is not in the correct OFF or COFF format, an error is raised.

