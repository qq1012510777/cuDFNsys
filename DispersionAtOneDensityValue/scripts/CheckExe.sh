#!/bin/bash
# source ~/.bashrc

let begin_=21          ########here change this variable!!!!!!!
let end_=30          ########here change this variable!!!!!!!

for i in $( seq $begin_ $end_)
do
    DFN_i=./DFN_$i

    cd $DFN_i

    echo $DFN_i" #############################################"

    log_i=log_$i.a.txt

    log_i_S=log_$i.s.txt

    echo $log_i

    tail -10  $log_i

    echo " "

    echo $log_i_S

    tail -10  $log_i_S

    echo " "

    cd ..
done