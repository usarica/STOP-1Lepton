#!/bin/sh


INLIST=$1
OUTLIST=$2
DATE=$3
if [[ "$DATE" == "" ]];then
  DATE=$(date +%y%m%d)
fi
tplline=""

if [ -z $INLIST ] || [ -z $OUTLIST ];then
  echo "Either input or output file is not determined..."
  exit 1
fi
rm -f $OUTLIST

let isFirstLine=1
while IFS='' read -r line || [[ -n "$line" ]]; do
  if [[ "$line" == "#"* ]];then
    continue
  fi
  if [ $isFirstLine -eq 1 ] || [[ "$line" == "template"* ]];then
    tplline=$line
    if [ $isFirstLine -eq 1 ];then
      let isFirstLine=0
    fi
    if [[ "$tplline" != "template"* ]];then
      echo "Input file does not contain a template line to write. Exiting..."
      exit 1
    fi
    tplline=${tplline//"template "}
    tplline=${tplline//".oO[DATE]Oo."/$DATE}
    echo "NEW TEMPLATE FOUND"
    echo "=> ${tplline}"

    stropts=($(echo $tplline))
    for stropt in ${stropts[*]}; do
      if [[ "$stropt" == "condoroutdir="* ]];then
        condoroutdir=${stropt//"condoroutdir="}
        if [[ "$condoroutdir" != *".oO"* ]];then
          echo "Creating the directory ${condoroutdir}"
          mkdir -p $condoroutdir
        fi
      fi
    done
  else
    arrIN=(${line//\// })
    OUTFILE=""
    if [[ "$line" == *"MINIAODSIM"* ]];then
      OUTFILE=${arrIN[0]}".root"
    elif [[ "$line" == *"MINIAOD"* ]];then
      OUTFILE=${arrIN[1]}"_"${arrIN[0]}".root"
    else
      echo "Cannot determine if the sample is data or MC!"
      exit 1
    fi
    OUTFILECORE=${OUTFILE%.*}
    SAMPLE=$line

    outline=$tplline
    outline=${outline//".oO[SAMPLE]Oo."/$SAMPLE}
    outline=${outline//".oO[OUTFILE]Oo."/$OUTFILE}
    outline=${outline//".oO[OUTFILECORE]Oo."/$OUTFILECORE}
    #echo $outline |& tee -a $OUTLIST

    outline=${outline//\"}
    larr=($(echo $outline))
    let larrsize=${#larr[@]}
    lsub=${outline//" "${larr[larrsize-1]}}
    SplitFrameworkJob "$lsub" nfiles="${larr[larrsize-1]}" outfile="${OUTLIST}"
  fi

done < "$INLIST"
