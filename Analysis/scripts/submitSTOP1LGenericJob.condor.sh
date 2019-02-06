#!/bin/sh


if [ -z ${_CONDOR_SCRATCH_DIR+x} ]; then
  #running locally
  runninglocally=true
  _CONDOR_SCRATCH_DIR=$(mktemp -d)
  SUBMIT_DIR=$(pwd)
else
  runninglocally=false
  SUBMIT_DIR=$1
fi


cd $SUBMIT_DIR
eval `scram runtime -sh`
echo "CMSSW VERSION: "$CMSSW_VERSION


runfile=$2
fcn=$3
extcmd=$4

cp $runfile $_CONDOR_SCRATCH_DIR
cd $_CONDOR_SCRATCH_DIR

echo "Job running in "$(pwd)" at "$(date)


runGenericROOTCommand.py --loadlib="loadLib.C" --script="$runfile" --function="$fcn" --command="$extcmd"
