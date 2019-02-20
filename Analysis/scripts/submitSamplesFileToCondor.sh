#!/bin/bash

INFILE=$1
OUTPUTDIR=$2
if [[ "$OUTPUTDIR" == "" ]];then
  OUTPUTDIR="./output"
fi
SCRIPTNAME="doFullProduction.cc"
FCN="doFullProduction"
QUEUE="default"
DATE=$(date +%y%m%d)

OUTDIR="${OUTPUTDIR}/${DATE}"

mkdir -p $OUTDIR

TARFILE="stop_1lepton.tar"
if [ ! -e ${OUTDIR}/${TARFILE} ];then
  createSTOP1LTarball.sh
  mv ${TARFILE} ${OUTDIR}/
fi


while IFS='' read -r line || [[ -n "$line" ]]; do
  fcnarglist=($(echo $line))
  fcnargname=""
  for farg in "${fcnarglist[@]}";do
    farg=${farg//"\""}
    if [[ "$farg" == "outfile="* ]];then
      fcnargname=$farg
      fcnargname=${fcnargname//"outfile="}
      fcnargname=${fcnargname//".root"}
      break
    fi
  done
  theOutdir="${OUTDIR}/${FCN}/${fcnargname}"
  mkdir -p $theOutdir
  ln -sf ${PWD}/${OUTDIR}/${TARFILE} ${PWD}/${theOutdir}/

  submitSTOP1LGenericJob.sh "$SCRIPTNAME" "$FCN" "$line" "$QUEUE" "$theOutdir"
done < "$INFILE"
