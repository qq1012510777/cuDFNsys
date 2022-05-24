#!/bin/bash
let begin_=1
let end_=10000
filename='log.txt '
for i in $( seq $begin_ $end_)
do
    ./bin/main 1000 100 0.01 0.5 > ${filename} 
    result=0;
    result=$(grep "drop" ${filename})

    if [ -z "${result}" ]; then
        echo "I did not find warning sign!"
    else
        echo "Warning released!"
        break
    fi
done