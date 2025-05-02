# tooth-evo-devo

Code steps to generate uniform random sampling of teeth using OPCR and cusp counts as the measurement for complexity (figures 3a & 4a):
1. python3 lhs.py 19000     # This will generate 19000 input .txt files using Latin Hypercube Sampling of the variable genetic and cellular parameters in ToothMaker
2. cp runt.e ./19000_input_files     #This will copy runt.e into the subdirectory with the input .txt files
3. python3 make_off_multi.py ./19000_input_files ./19000_input_files     #this makes a "multirun" code to run to generate 19000 teeth in the cluster 
4. cd ./19000_input_files #navigate to directory with input files
5. sed -i '/\.\/runt\.e \.\/multirun\.txt \. multirun 9000 1/d' multirun.txt     #this deletes the last line of the multirun file, if neccessary, which may say ./runt.e ./multirun.txt . multirun 5000 1 in multirun.txt
6. addqueue -n 50 /usr/local/shared/bin/multirun ./multirun.txt     #for running calculations on a cluster, this will generate all the tooth .off files inside ./19000_input_files
7. ls -1 ./*.off | wc -l    #this will tell you how many .off files were generated
8. cp convert.off.to.ply.py ./19000_input_files     #this is a file which will convert all the .off files into .ply files
9. python3 convert.off.to.ply.py .    #this will convert all the .off files into .ply files
10. cd ./19000_input_files/outputs_ply    #navigate to directory with .ply tooth files
11. ls -1 ./*.ply | wc -l    #should be 19000!
12. mkdir opc    #this will be the directory where we save the orientation patch count rotated values for each .ply file
13.  (copy the following files to ./outputs_ply : calc_opc.py, topomesh.py, plython.py, DNE.py, implicitfair.py, normcore.py, OPC.py, RFI.py, 1opc.py, run_opc.py)
14.  python3 run_opc.py    #this will generate a .txt file called commands_opc.txt, which is a command to calculate the opc for each .ply tooth file in an efficient manner
15.  addqueue -c "5 minutes" -n 100 /usr/local/shared/bin/multirun ./commands_opc.txt    #this will calculate the OPC values for all the .ply files in the directory and save them in /opc
16.  

#check_height
This script processes all .off files in the current directory, checks the height (difference between the maximum and minimum z-coordinates) of each file, and classifies them as either "tall enough" or "too flat" based on a threshold of 0.5. The results are written to two separate output files: 18985_tall_enough.txt and 18985_too_flat.txt. If the file is not in the correct OFF or COFF format, an error is raised.

