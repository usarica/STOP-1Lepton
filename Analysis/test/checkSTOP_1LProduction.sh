#!/bin/bash

chkdir=$1

let nOK=0
let nCOPYFAIL=0
let nFAIL=0

cd $chkdir

for d in $(ls ./); do
  logfilename=""
  for ldd in $(ls ./$d/Logs); do
    if [[ "$ldd" == "log_"* ]];then
      logfilename=$ldd
    fi
  done
  if [[ -z $logfilename ]];then
    continue
  fi
  res=$(grep -e "File generation was successful" $d/Logs/$logfilename)
  res2=$(grep -e "Copied successfully" $d/Logs/$logfilename)
  if [[ ! -z $res ]];then
    if [[ ! -z $res2 ]];then
      echo $d" ran successfully"
      let nOK=$nOK+1
    else
      echo $d" file ok but copy failed"
      let nCOPYFAIL=$nCOPYFAIL+1
    fi
  elif [ -s $d/Logs/$logfilename ];then
    echo $d" failed"
    let nFAIL=$nFAIL+1
  fi
done

cd -

echo "(OK:COPY_FAIL:FAIL) = (${nOK}:${nCOPYFAIL}:${nFAIL}"
