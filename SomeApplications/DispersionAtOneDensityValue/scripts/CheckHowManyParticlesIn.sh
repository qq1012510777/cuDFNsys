#!/bin/bash

let begin_=21          ########here change this variable!!!!!!!
let end_=30          ########here change this variable!!!!!!!

for i in $( seq $begin_ $end_)
do
    DFN_i=./DFN_$i

    if  [ -d "$DFN_i" ]
    then
        cd $DFN_i/ParticlePositionResult

        echo $DFN_i" ##################"
        str=$(h5dump -A  ParticlePositionLastStep.h5)
        str1=$(echo $str | cut -d '{' -f5 )
        str2=$(echo $str1 | cut -d ' ' -f3)
        echo $str2" particles still in"
        echo " "
        cd ../../
    else
        echo $DFN_i" is not existing ************"
        echo " "
    fi

done 