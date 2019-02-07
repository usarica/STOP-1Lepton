#!/bin/sh

SCRIPTNAME=$1
FCN=$2
FCNARGS=$3
QUEUE=$4
OUTDIR=$5

if [[ "$OUTDIR" == "" ]];then
  echo "You must set the output directory!"
  exit 1
fi

SCRIPTSO=${SCRIPTNAME%.*}"_"${SCRIPTNAME##*.}".so"

LOGSDIR=$OUTDIR"/Logs"

CMSENVDIR=$CMSSW_BASE
if [[ "$CMSENVDIR" == "" ]];then
  echo "Set up CMSSW first!"
  exit 1
fi


mkdir -p $LOGSDIR

extLog=$FCN
if [[ "$FCNARGS" != "" ]];then
  fcnargname=${FCNARGS//\"}
  fcnargname=${FCNARGS//\\}
  fcnargname=${fcnargname//"("}
  fcnargname=${fcnargname//")"}
  fcnargname=${fcnargname//","/"_"}
  extLog=$FCN"_"$fcnargname
fi


if [[ -f $SCRIPTNAME ]]; then
  echo "File "$SCRIPTNAME" exists."

  if [[ ! -f $SCRIPTSO ]]; then
    echo "Compiling "$SCRIPTNAME"..."
    root -l -b -q -e "gROOT->ProcessLine(\".x loadLib.C\"); gROOT->ProcessLine(\".L "$SCRIPTNAME"+\");"
  fi

  hname=$(hostname)
  echo $hname
  if [[ "$hname" == *"lxplus"* ]] || [[ "$hname" == *"ucsd"* ]];then
    echo "Host is on LXPLUS or UCSD, so need to use HTCONDOR"
    THEQUEUE="vanilla"
    if [[ "$QUEUE" != "default" ]];then
      THEQUEUE=$QUEUE
    fi
    checkGridProxy.sh
    createSTOP1LTarball.sh
    configureSTOP1LCondorJob.py --batchqueue="$THEQUEUE" --outdir="$OUTDIR" --outlog="Logs/log_$extLog.txt" --errlog="Logs/err_$extLog.err" --batchscript="submitSTOP1LGenericJob.condor.sh" --script="$SCRIPTNAME" --fcn="$FCN" --fcnargs="$FCNARGS"
  elif [[ "$hname" == *"login-node"* ]] || [[ "$hname" == *"bc-login"* ]]; then
    echo "Host is on MARCC, so need to use SLURM batch"
    THEQUEUE="lrgmem"
    if [[ "$QUEUE" != "default" ]];then
      THEQUEUE=$QUEUE
    fi
    sbatch --output=$LOGSDIR"/log_"$extLog".txt" --error=$LOGSDIR"/err_"$extLog".err" --partition=$THEQUEUE submitSTOP1LGenericJob.slurm.sh "$CMSENVDIR" "$SCRIPTNAME" "$FCN" "$FCNARGS"
  fi

fi
