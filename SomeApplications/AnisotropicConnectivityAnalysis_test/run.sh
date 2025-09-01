#!/bin/bash

# Check if at least one argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <argument>"
    exit 1
fi

# Access the first argument
setno=$1

#source ~/.bashrc
#module load cuda/12.0
executablePath=/home/tingchang/cuDFNsys/SomeApplications/AnisotropicConnectivityAnalysis_test

setfilename="set"$(printf "%05d" "$setno")

echo "set name is "$setfilename

mkdir -p ${setfilename}
cd ${setfilename}

${executablePath}/main \
    200 1 1 1 2 \
    25 0  0  1 32 32 34 0 1.5             1           100             0 0.13 1.7e-05 \
    45 1           0  6.1232e-17 48 48 34 0 2.5            1           30            0 0.12 1.6e-05
        
# IfContinue=1
# count_loop=1
# while [ $IfContinue -eq 1 ]
# do
#     echo "loop "${count_loop}"---------------------------"
#     ${executablePath}/main \
#         250 1 1 1 2 \
#         22 1           0  6.1232e-17 71 71 29 0 2.4            1           25            0 0.24 1.2e-08 \
#         14 0  0  1 165 165 29 0 1.3            1           50            0 0.31 3.3e-07 
# 
#     echo "start checking if finished. Time: $(date)"
#     for i in {1..30}
#     do
#         DFNFileName=DFN_$(printf "%03d" "$i")
#         cd $DFNFileName
#         filename="PTFinished"
#         # Check if the file exists
#         if [ -e "$filename" ]; then
#             #echo "File $filename exists in the current directory."
#             echo $DFNFileName" finished"
#             IfContinue=0
#         else
#             echo $DFNFileName" did not finish"
#             IfContinue=1
#             cd ..
#             break
#         fi
#         cd ..
#     done
#     echo "--------------------------"
#     echo " "
#     ((count_loop++))
# done

echo $setfilename': task done'
