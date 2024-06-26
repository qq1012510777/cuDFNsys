#!/bin/bash

setNO=30
DFNNo=30

current_dir=$(pwd)

dataFile=$current_dir/puredata
mkdir -p $dataFile
for ((i = 1; i <= $setNO; i++))
do
    setName=$current_dir/set$(printf "%05d" "$i")

    if [ -d "$setName" ]; then
        donothing=1
    else
        echo  "set$(printf "%05d" "$i") does not exist."
        continue
    fi

    mkdir -p $dataFile/set$(printf "%05d" "$i")

    for ((j = 1; j <= $DFNNo; j++))
    do
        DFNname=$setName/DFN_$(printf "%03d" "$j")
        if [ -d "$DFNname" ]; then
            donothing=1
        else
            continue
        fi

        targetDir=$dataFile/set$(printf "%05d" "$i")/DFN_$(printf "%03d" "$j")
        mkdir -p $targetDir

        cp -r $DFNname/Dispersion_MeanSquareDisplacement.h5  $targetDir
        cp -r $DFNname/ParticlePositionResult $targetDir

    done

done