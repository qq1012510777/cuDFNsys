#!/bin/bash
source ~/.bashrc

if [ ! -n "$1" ]; then
    echo "Please define the model NO as the 1st input"
    exit 0
fi 

if [ ! -n "$2" ]; then
    echo "Please define the model NO as the 1st input"
    exit 0
fi 

let begin_=$1           ########here change this variable!!!!!!!
let end_=$2          ########here change this variable!!!!!!!

NUM_fracs=740
ModelSize=75
Kappa=0
Beta=0.1
Gamma=5.1e-4
SizeFracMode=0
ParaSizeDistri=(1.5 1 15 0)
DomainDimensionRatioZ=1
IfRemoveDeadEnds=0
MinGridSize=1
MaxGridSize=8
NumStepsDispersion=200
NumParticlesRandomWalk=700000
FactorMeanTimeInGrid=0.3  # this important # change this # change this
Lengthscale=5.4772
Pe=100 # change this # change this
ControlPlaneSpacing=1125
IfoutputMsd=1
IfoutputParticleInfoAllsteps=0
# change above
MaxExeTimes=1000
SpecialChar="L75s_Pe100" # change this # change this  # change this # change this
LimitationOfArrivedParticles=2000 # if there are less than $LimitationOfArrivedParticles particles remain in the model, then discontinue simulations
###################################3

#commands="/storage/torresLab/yintingchang/cuDFNsys/DispersionAtOneDensityValue/main 230 50 0 0.1 5.1e-4 0 1.5 1 15 0 1 0 1 3 100 20 200 500000 25  0 10000"
commands="/storage/torresLab/yintingchang/cuDFNsys/DispersionAtOneDensityValue/main ${NUM_fracs} ${ModelSize} ${Kappa} ${Beta} ${Gamma} ${SizeFracMode} ${ParaSizeDistri[0]} ${ParaSizeDistri[1]} ${ParaSizeDistri[2]} ${ParaSizeDistri[3]} ${DomainDimensionRatioZ} ${IfRemoveDeadEnds} ${MinGridSize} ${MaxGridSize} ${NumStepsDispersion} ${NumParticlesRandomWalk} ${FactorMeanTimeInGrid} ${Lengthscale} ${Pe} ${ControlPlaneSpacing} ${IfoutputMsd} ${IfoutputParticleInfoAllsteps}"

sign=1
CountTimes=0

LLog_Errors="LLog_Errors.txt"
echo "" > $LLog_Errors

Myarray=()
ProgramErrorCountor=()
ShellErrorCounter=()
i_a=1
for i in $( seq ${i_a} `expr $end_ + 1`); do
    Myarray+=(0)
    ProgramErrorCountor+=(0)
    ShellErrorCounter+=(0)
done

while [ $sign -gt 0 ]; do
    CountTimes=`expr $CountTimes + 1`

    echo " "

    for i in $( seq $begin_ $end_)
    do

        sleep 0.25

        if [ ${Myarray[$i]} -gt ${MaxExeTimes} ];
        then
            continue
        fi

        DFN_i=./DFN_$i

        if  [ -d "$DFN_i" ]
        then
            cd DFN_$i
        else
            echo "create DFN_$i"
            mkdir DFN_$i
            cd DFN_$i
        fi

        #########h5dump

        if [ -f "./ParticlePositionResult/ParticlePositionLastStep.h5" ]
        then
            strsd=$(h5dump -A  ./ParticlePositionResult/ParticlePositionLastStep.h5)
            str1sd=$(echo $strsd | cut -d '{' -f5 )
            str2sd=$(echo $str1sd | cut -d ' ' -f3)
            if [[ "$str2sd" == *")"*  ]]
            then
                echo "1 particles still in "$DFN_i
                str2sd="1"
            else
                if [ -z "$str2sd" ]
                then
                    echo "0 particles still in "$DFN_i
                    str2sd="0"
                else
                    echo "$str2sd particles still in "$DFN_i
                fi
            fi
        
            NumParticlesInDFN=$(($str2sd))

            if [ ${NumParticlesInDFN} -gt $LimitationOfArrivedParticles ]
            then
                donothing=1 
            else
                Myarray[$i]=`expr ${MaxExeTimes} + 1`
                cd ..
                continue # do not need more steps
            fi
        else
            echo $DFN_i" has no ./ParticlePositionResult/ParticlePositionLastStep.h5"
        fi
        #################

        log_i=log_$i.a.txt

        log_i_S=log_$i.s.txt

        zeroPadding_i=`printf %05d $i`
        ExeName_i=MyMainDispersion_${zeroPadding_i}_$SpecialChar

        if [ -f "$log_i" ] ## the program has been run once
        then

            # first check if the program is running
            str_nvidia_smi=$(nvidia-smi --query-compute-apps=process_name --format=csv)
            #str_nvidia_smi=$(nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv)

            if [[ "$str_nvidia_smi" == *"$ExeName_i"* ]];
            then
                echo $ExeName_i" is running ..."
                cd ..
                continue
            fi

            # so now the program is not running
            # we check if the program is terminated due to shell's standard error like segmentation fault
            # this kind of error cannot be captured by the program 
            if [ -f "$log_i_S" ]
            then
                strds=$(cat ${log_i_S})
                
                if [ -z "$strds" ] # empty
                then
                    donothing=1
                else # not empty, meaning that a shell standard error happens
                    echo "Shell standard errors happen in $DFN_i"
                    echo $strds
                    echo $strds >> ../$LLog_Errors
                    echo " " >> ../$LLog_Errors
                    echo "Run a new DFN (shell errors)"
                    #exit
                    rm -rf ./*.h5 ./ParticlePositionResult
                    { bash -c "exec -a $ExeName_i $commands > $log_i"; } &> $log_i_S &
                    Myarray[$i]=1
                    ShellErrorCounter[$i]=`expr ${ShellErrorCounter[$i]} + 1`
                    cd ..
                    continue    
                fi
            fi    

            ## the program has been finished
            if grep -wq "No percolation happens" $log_i
            then
                # run code
                echo "No percolation happens in $DFN_i"
                echo "Run a new DFN (former no percolative)"
                rm -rf ./*.h5 ./ParticlePositionResult
                { bash -c "exec -a $ExeName_i $commands > $log_i"; } &> $log_i_S & 
                Myarray[$i]=1
                cd ..
                continue
            fi

            ## the program has been finished, and the percolation happened already
            if grep -wq "This simulation consumes" $log_i
            then
                # run code
                echo "Finished simulation in $DFN_i"
                echo "Run more steps"
                { bash -c "exec -a $ExeName_i $commands > $log_i"; } &> $log_i_S & 
                Myarray[$i]=`expr ${Myarray[$i]} + 1`
                cd ..
                continue
            fi

            ## the program has failed due to an error captured by the program    
            if grep -wq "Failed simulation" $log_i
            then

                # run code
                echo "Failed simulation in $DFN_i -------"
                tail -50 ${log_i}
                tail -50 ${log_i} >> ../$LLog_Errors
                echo " " >> ../$LLog_Errors
                echo "Run a new DFN (former failed)"
                rm -rf ./*.h5 ./ParticlePositionResult
                { bash -c "exec -a $ExeName_i $commands > $log_i"; } &> $log_i_S & 
                Myarray[$i]=1
                ProgramErrorCountor[$i]=`expr ${ProgramErrorCountor[$i]} + 1`
                cd ..
                continue
               
            fi

            # if log_i exists but the program is discontinued unexceptionally
            # e.g., it is interrupted by something
            # so, no substrings like "No percolation happens", "This simulation consumes" and "Failed simulation"
            # we can try to re-run the program
            echo "Unknown error in $DFN_i"
            echo "Run more steps"
            { bash -c "exec -a $ExeName_i $commands > $log_i"; } &> $log_i_S & 
            Myarray[$i]=`expr ${Myarray[$i]} + 1`

        else # no simulation has been done so far

            echo "Start simulation in $DFN_i"
            #echo ${commands}
            { bash -c "exec -a $ExeName_i $commands > $log_i"; } &> $log_i_S &  
            Myarray[$i]=`expr ${Myarray[$i]} + 1`    
        fi

        cd .. 

    done    

    for j in $( seq $begin_ $end_); do
        echo "DFN_$j has run "${Myarray[$j]}" times, ShellErrorCounter: "${ShellErrorCounter[$j]}", ProgramErrorCountor: "${ProgramErrorCountor[$j]}"."
    done

    echo 'WhileTimes='$CountTimes

    sleep 1

    #### if the number of particles retained in the model is smallar than a value, then exit
    #### if the number of particles retained in the model is smallar than a value, then exit
    #### if the number of particles retained in the model is smallar than a value, then exit
    #### if the number of particles retained in the model is smallar than a value, then exit
    Iuy=0
    for i in $( seq $begin_ $end_)
    do

        DFN_i=./DFN_$i

        #########h5dump

        if [ -f "./${DFN_i}/ParticlePositionResult/ParticlePositionLastStep.h5" ]
        then
            strsd=$(h5dump -A  ./${DFN_i}/ParticlePositionResult/ParticlePositionLastStep.h5)
            str1sd=$(echo $strsd | cut -d '{' -f5 )
            str2sd=$(echo $str1sd | cut -d ' ' -f3)
            #echo $str2sd" particles still in "$DFN_i

            if [[ "$str2sd" == *")"*  ]]
            then
                str2sd="1"
            fi
            if [ -z "$str2sd" ]
            then
                str2sd="0"                
            fi

            NumParticlesInDFN=$(($str2sd))

            if [ ${NumParticlesInDFN} -gt $LimitationOfArrivedParticles ]
            then 
                Iuy=1
            else
                donothing=1 ###   Iuy is still zero 
            fi
        else
            #echo $DFN_i" has no ./${DFN_i}/ParticlePositionResult/ParticlePositionLastStep.h5"
            #donothing=1 
            break
        fi

        if [ $i -eq $end_ ]
        then
            if [ $Iuy -eq 0 ]
            then
                echo " "
                echo "Enough particles arrived!"
                echo " "
                exit # 
            fi
        fi
    done

    #### random walk steps reached the maximum values
    for i in $( seq $begin_ $end_)
    do
        if [ ${Myarray[$i]} -gt ${MaxExeTimes} ];
        then
            if [ $i -eq $end_ ]
            then
                echo " "
                echo "Random walk steps reached the maximum values!"
                echo " "
                exit # 
            fi
        else
            break
        fi
    done

done
