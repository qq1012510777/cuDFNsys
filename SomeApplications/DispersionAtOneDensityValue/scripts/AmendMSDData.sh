#!/bin/bash
#MyArray=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40)
MyArray=($(seq 1 1 40))
for i in "${MyArray[@]}"
do
    DFN_i=./DFN_$i
    if  [ -d "$DFN_i" ]
    then
        cd $DFN_i
        echo $DFN_i" ##########"
        str=$(python3 -V)
        if [ -z "$str" ];
        then
            echo "PythonVersion Problem"
            exit
        fi
        python3 ../DelDataSet.py 
        cd ..
    else
        echo $DFN_i" is not existing ************"
        echo " "
    fi
done 