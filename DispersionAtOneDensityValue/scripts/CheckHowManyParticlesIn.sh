#!/bin/bash

let begin_=1          ########here change this variable!!!!!!!
let end_=40          ########here change this variable!!!!!!!

for i in $( seq $begin_ $end_)
do
    DFN_i=./DFN_$i

    if  [ -d "$DFN_i" ]
    then

        if [ -d "$DFN_i/ParticlePositionResult" ]
        then
            s=1
        else
            echo $DFN_i"ParticlePositionResult is not existing ************"
            echo " " 
            continue
        fi
        
        cd $DFN_i/ParticlePositionResult

        echo $DFN_i" ##################"
        h5dump -A  ParticlePositionLastStep.h5
        echo " "
        cd ../../
    else
        echo $DFN_i" is not existing ************"
        echo " "
    fi

done 