#!/bin/bash

chkdir=$1

let nOK=0
let nCOPYFAIL=0
let nFAIL=0
let nFILEDNE=0
let nUNKNOWN=0

cd $chkdir

for d in $(ls ./); do
  if [[ ! -d $d ]];then
    continue
  fi

  let dirok=1
  for logfilename in $(ls ./$d/Logs | grep -e "log_"); do
    res=$(grep -e "File generation was successful" $d/Logs/$logfilename)
    res2=$(grep -e "Copied successfully" $d/Logs/$logfilename)

    if [[ ! -z $res ]];then
      if [[ ! -z $res2 ]];then
        res3=$(grep -e "Running: env -i " $d/Logs/$logfilename)
        res3=${res3//*"gsiftp://gftp.t2.ucsd.edu"}
        if [[ -s $res3 ]];then
          echo $d" ran successfully"
          let nOK=$nOK+1
        else
          echo $d" ran successfully, but the file does not exist!"
          let nFILEDNE=$nFILEDNE+1
          let dirok=0
        fi
      else
        echo $d" file ok but copy failed"
        let nCOPYFAIL=$nCOPYFAIL+1
        let dirok=0
      fi
    elif [ -s $d/Logs/$logfilename ];then
      echo $d" failed"
      let nFAIL=$nFAIL+1
      let dirok=0
    else
      echo $d" status is not yet determined"
      let nUNKNOWN=$nUNKNOWN+1
      let dirok=0
    fi

  done

  if [[ $dirok -eq 1 ]];then
    TARFILE="${d}.tar"
    rm -f $TARFILE
    tar Jcf ${TARFILE} $d --exclude={*.tar}
    if [[ $? -eq 0 ]];then
      echo "- Compressed successfully, so removing the directory"
      rm -rf $d
    fi
  fi

done

cd -

echo "(OK:COPY_FAIL:FILE_DNE:FAIL:UNKNOWN) = (${nOK}:${nCOPYFAIL}:${nFILEDNE}:${nFAIL}:${nUNKNOWN})"
