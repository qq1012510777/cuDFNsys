#!/bin/bash
source ~/.bashrc

let begin_=1          ########here change this variable!!!!!!!
let end_=5          ########here change this variable!!!!!!!

NUM_fracs=3400
ModelSize=125
Kappa=0
Beta=0.1
Gamma=5.1e-4
SizeFracMode=0
ParaSizeDistri=(1.5 1 15 0)
DomainDimensionRatioZ=1
IfRemoveDeadEnds=0
MinGridSize=1
MaxGridSize=3
P_in=100
P_out=20
NumStepsDispersion=200
NumParticlesRandomWalk=500000
FactorMeanTimeInGrid=10  # this important
MolecuDiffu=0
ControlPlaneSpacing=25
IfoutputMsd=0
# change above
MaxExeTimes=500
SpecialChar="L125_I" # change this # change this  # change this # change this

###################################3

#commands="/storage/torresLab/yintingchang/cuDFNsys/DispersionAtOneDensityValue/main 230 50 0 0.1 5.1e-4 0 1.5 1 15 0 1 0 1 3 100 20 200 500000 25  0 10000"
commands="/storage/torresLab/yintingchang/cuDFNsys/DispersionAtOneDensityValue/main ${NUM_fracs} ${ModelSize} ${Kappa} ${Beta} ${Gamma} ${SizeFracMode} ${ParaSizeDistri[0]} ${ParaSizeDistri[1]} ${ParaSizeDistri[2]} ${ParaSizeDistri[3]} ${DomainDimensionRatioZ} ${IfRemoveDeadEnds} ${MinGridSize} ${MaxGridSize} ${P_in} ${P_out} ${NumStepsDispersion} ${NumParticlesRandomWalk} ${FactorMeanTimeInGrid} ${MolecuDiffu} ${ControlPlaneSpacing} ${IfoutputMsd}"

sign=1
CountTimes=0

LLog_Errors="LLog_Errors.txt"
echo "" > $LLog_Errors

Myarray=()
ProgramErrorCountor=()
ShellErrorCounter=()
for i in $( seq $1 `expr $end_ + 1`); do
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
            echo $str2sd" particles still in "$DFN_i
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

done