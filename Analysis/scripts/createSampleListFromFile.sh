#!/bin/sh


INLIST=$1
OUTLIST=$2
DATE=$3
ISDATA=$4
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
    echo "NEW TEMPLATE FOUND"
    echo "\t=> ${tplline}"
  else
    arrIN=(${line//\// })
    OUTFILE=${arrIN[0]}".root"
    if [[ "$ISDATA" == "data" ]];then
      OUTFILE=${arrIN[1]}"_"${arrIN[0]}".root"
    fi
    SAMPLE=$line
    outline=$tplline
    outline=${outline//".oO[SAMPLE]Oo."/$SAMPLE}
    outline=${outline//".oO[OUTFILE]Oo."/$OUTFILE}
    outline=${outline//".oO[DATE]Oo."/$DATE}
    #echo $outline |& tee -a $OUTLIST

    outline=${outline//\"}
    larr=($(echo $outline))
    let larrsize=${#larr[@]}
    lsub=${outline//" "${larr[larrsize-1]}}
    SplitFrameworkJob "$lsub" nfiles="${larr[larrsize-1]}" outfile="${OUTLIST}"
  fi

done < "$INLIST"
