#!/bin/bash

### Run it in your (supercomputer) terminal (under slurm): bash CheckCpuGpuNodesResources.sh

RED='\033[0;31m'
NC='\033[0m' # No Color
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BOLD=$(tput bold)
NORMAL=$(tput sgr0)

echo " "
echo -e "${RED}${BOLD}Script's name: CheckCpuGpuNodesResources.sh${NC}"
echo -e "${RED}${BOLD}Purpose: output the resource information of CPU/GPU nodes (under slurm)${NC}"
echo -e "${RED}${BOLD}Author: Tingchang YIN, Westlake University, yintingchang@foxmail.com${NC}"
echo " "

### Check if slurm is installed or not
### Check if slurm is installed or not
### Check if slurm is installed or not
logh=log_CheckCpuGpuNodesResources.log
strE=$(sinfo -V)
if [ -z "${strE}" ]; then
        echo -e "${RED}${BOLD}Slurm is not installed${NC}"
        exit
fi
sinfo -V >${logh}
numLines=$(awk 'END {print NR}' ${logh})
numLines=$(($numLines))
if [ $numLines -gt 1 ]; then
        echo -e "${RED}${BOLD}Slurm is not installed${NC}"
        exit
fi
if grep -wq "not found" ${logh}; then
        echo -e "${RED}${BOLD}Slurm is not installed${NC}"
        exit
fi
if grep -wq "bash" ${logh}; then
        echo -e "${RED}${BOLD}Slurm is not installed${NC}"
        exit
fi

#node names
#node names
#node names
# NodeName=('gvno02')
NodeName=('xnode09'
        'xnode24'
        'grtq14'
        'grtq15')
for ((i = 2; i < 13; i++)); do
        TU=$(printf %02d $i)
        NodeName+=('gvnq'$TU)
done
NodeName+=('gvnq14')
NodeName+=('gvno01' 'gvno02' 'gvno03' 'grtq14' 'grtq15')
for ((i = 1; i < 20; i++)); do
        TU=$(printf %02d $i)
        NodeName+=('ga40q'$TU)
done
NodeName+=('gvna01' 'gvna02')
#-------------------------------------
# NodeName=('xnode09')
# sinfo -N --noheader >${logh}
# numLines=$(awk 'END {print NR}' ${logh})
# numLines=$(($numLines))
# NodeName=('grtq14')
# for ((i = 1; i <= $numLines; i++)); do
#         strUI="${i}q;d"
#         strTY=$(sed ${strUI} ${logh})
#         strTY=$(echo $strTY | cut -d ' ' -f 1)
#         NodeName+=(${strTY})
# done

# start checking
# start checking
# start checking
len=${#NodeName[@]}
for ((i = 0; i < $len; i++)); do
        Node_i=${NodeName[$i]}

        scontrol show node $Node_i >${logh}

        if grep -wq "not found" ${logh}; then
                echo -e ${RED}$BOLD"Node: "$Node_i" is not found"${NC}
                continue
        fi

        IfGPU=0
        if grep -wq "gpu" ${logh}; then
                IfGPU=1
        fi

        # cpu
        line1=$(sed '2q;d' ${logh})
        str1=$(echo $line1 | cut -d '=' -f 2)
        str2=$(echo $str1 | cut -d ' ' -f 1)

        str3=$(echo $line1 | cut -d '=' -f 4)
        str4=$(echo $str3 | cut -d ' ' -f 1)

        # gpu
        str5=""
        str6=""

        if [ $IfGPU -gt 0 ]; then
                line5=$(sed '5q;d' ${logh})
                str5=$(echo $line5 | cut -d ':' -f 2)
                line14=$(sed '14q;d' ${logh})
                str6=$(echo $line14 | cut -d '=' -f 5)
                if [ -z "${str6}" ];
                then
                        str6="0"       
                fi
        else
                str5="0"
                str6="0"
        fi

        # memory
        line8=$(sed '8q;d' ${logh})
        str7=$(echo $line8 | cut -d ' ' -f 1)
        str8=$(echo $str7 | cut -d '=' -f 2)
        RealMemory=$(($str8))
        RealMemory=$(printf "%f\n" $((RealMemory / 1024)))

        str9=$(echo $line8 | cut -d ' ' -f 2)
        str10=$(echo $str9 | cut -d '=' -f 2)
        OccupiedMemory=$(($str10))
        OccupiedMemory=$(printf "%f\n" $((OccupiedMemory / 1024)))

        # under which partition
        line10=$(sed '10q;d' ${logh})
        str11=$(echo $line10 | cut -d '=' -f 2)

        if [ $IfGPU -gt 0 ]; then
                echo -e ${GREEN}$BOLD$Node_i", partition: "$str11", CPU: "$str2"/"$str4", ${YELLOW}GPU: "$str6"/"$str5", ${GREEN}${BOLD}memory: "$OccupiedMemory"G/"$RealMemory"G"${NC}
        else
                echo -e $BOLD$Node_i", partition: "$str11", CPU: "$str2"/"$str4", GPU: "$str6"/"$str5", memory: "$OccupiedMemory"G/"$RealMemory"G"
        fi
done

echo " "$NORMAL
rm -rf ./$logh