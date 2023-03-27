#!/bin/bash

let begin_=1          ########here change this variable!!!!!!!
let end_=8          ########here change this variable!!!!!!!

commands="../../main 230 50 0 0.1 5.1e-4 0 1.5 1 15 0 1 0 1 3 100 20 2000 5000 30  0 20"

sign=1
CountTimes=0

Myarray=()
for i in $( seq $begin_ $end_); do
    Myarray+=(0)
done

while [ $sign -gt 0 ]; do
    CountTimes=`expr $CountTimes + 1`
    for i in $( seq $begin_ $end_)
    do
        sleep 0.25
        DFN_i=./DFN_$i

        if  [ -d "$DFN_i" ]
        then
            cd DFN_$i
        else
            echo "create DFN_$i"
            mkdir DFN_$i
            cd DFN_$i
        fi

        log_i=log_$i.txt

        # rm -rf ./SimulationFailed.txt ./NoPercolation.txt ./SimulationFinished.txt

        if [ -f "$log_i" ]
        then
            sign_d=1

            #if grep -wq "No percolation happens" $log_i ## the program has been finished
            if grep -wq "No percolation happens" $log_i
            then
                # run code
                echo "No percolation happens in $DFN_i"
                echo "Run a new DFN"
                rm -rf ./*.h5 ./ParticlePositionResult
                $commands > $log_i &
                Myarray[$i]=1
                sign_d=0
            fi

            if [ $sign_d -gt 0 ]
            then    
                if grep -wq "This simulation consumes" $log_i ## the program has been finished
                then
                    # run code
                    
                    echo "Finished simulation in $DFN_i"
                    echo "Run more steps"
                    $commands > $log_i &
                    Myarray[$i]=`expr ${Myarray[$i]} + 1`
                    
                fi
            fi

            if grep -wq "Failed simulation" $log_i  ## the program has failed
            then
                # run code
                echo "Failed simulation in $DFN_i"
                echo "Run a new DFN"
                rm -rf ./*.h5 ./ParticlePositionResult
                $commands > $log_i &
                Myarray[$i]=1
            fi

        else # no simulation has been done so far
            echo "Start simulation in $DFN_i"
            $commands > $log_i &
            Myarray[$i]=`expr ${Myarray[$i]} + 1`
        fi

        cd .. 

    done    

    for j in $( seq $begin_ $end_); do
        echo "DFN_$j has run "${Myarray[$j]}" times"
    done
    echo 'WhileTimes='$CountTimes
    sleep 1
done