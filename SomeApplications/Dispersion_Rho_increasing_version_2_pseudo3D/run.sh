#!/bin/bash

# Check if at least one argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <argument>"
    exit 1
fi

# Access the first argument
setno=$1

source ~/.bashrc
#module load cuda/12.0
executablePath=~/cuDFNsys/SomeApplications/Dispersion_Rho_increasing_version_2_pseudo3D/

setfilename="set"$(printf "%05d" "$setno")

echo "set name is "$setfilename

mkdir -p ${setfilename}
cd ${setfilename}

IfContinue=1
count_loop=1
while [ $IfContinue -eq 1 ]
do
    echo "loop "${count_loop}"---------------------------"
    ${executablePath}/main ${executablePath} 75  1 1 1 2 0 420 40 29  3 7.5 0 0 0 0.1 5.1e-4 1 2 300000 4000 Flux-weighted 0 10000 1 0 0 1 400 1  300  0  0.2  0.05 0.99 5 300000

    echo "start checking if finished. Time: $(date)"
    for i in {1..30}
    do
        DFNFileName=DFN_$(printf "%03d" "$i")
        cd $DFNFileName
        filename="PTFinished"
        # Check if the file exists
        if [ -e "$filename" ]; then
            #echo "File $filename exists in the current directory."
            echo $DFNFileName" finished"
            IfContinue=0
        else
            echo $DFNFileName" did not finish"
            IfContinue=1
            cd ..
            break
        fi
        cd ..
    done
    echo "--------------------------"
    echo " "
    ((count_loop++))
done

echo $setfilename': task done'
