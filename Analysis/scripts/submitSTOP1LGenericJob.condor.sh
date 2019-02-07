#!/bin/sh


if [ -z ${_CONDOR_SCRATCH_DIR+x} ]; then
  #running locally
  runninglocally=true
  _CONDOR_SCRATCH_DIR=$(mktemp -d)
  SUBMIT_DIR=$(pwd)
else
  runninglocally=false
  SUBMIT_DIR=$3
fi

version=$4
runfile=$5
fcn=$6
extcmd=$7

echo "Host: "$(hostname)
echo "User: "$USER
echo "Cluster id: "$1
echo "Process id: "$2
echo "Condor scratch directory: "$_CONDOR_SCRATCH_DIR
echo "Submission directory: "$SUBMIT_DIR
echo "Script command to run: "$runfile"::""$fcn""(""$extcmd"")"
echo "Job running in "$(pwd)" at "$(date)
echo "Can I see hadoop?"
ls /hadoop/cms/store/user/usarica

mkdir $version
tar -xvf stop_1lepton.tar -C $version
cd $version
#cd $SUBMIT_DIR
eval `scram runtime -sh`
echo "CMSSW VERSION: "$CMSSW_VERSION
cd -

cp $runfile $_CONDOR_SCRATCH_DIR
cd $_CONDOR_SCRATCH_DIR

runGenericROOTCommand.py --loadlib="loadLib.C" --script="$runfile" --function="$fcn" --command="$extcmd"
