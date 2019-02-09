#!/bin/sh


CMSSWVERSION=$1
SCRAMARCH=$2
SUBMIT_DIR=$3 # Must be within $CMSSW_BASE/src/
TARFILE=$4
RUNFILE=$5
FCN=$6
FCNARGS=$7

# Make sure OUTPUTNAME doesn't have .root since we add it manually
#OUTPUTNAME=$(echo $OUTPUTNAME | sed 's/\.root//')

export SCRAM_ARCH=${SCRAMARCH}

echo -e "\n--- begin header output ---\n" #                     <----- section division
echo "CMSSWVERSION: $CMSSWVERSION"
echo "SCRAMARCH: $SCRAMARCH"
echo "SUBMIT_DIR: $SUBMIT_DIR"
echo "TARFILE: $TARFILE"
echo "RUNFILE: $RUNFILE"
echo "FCN: $FCN"
echo "FCNARGS: $FCNARGS"

echo "GLIDEIN_CMSSite: $GLIDEIN_CMSSite"
echo "hostname: $(hostname)"
echo "uname -a: $(uname -a)"
echo "time: $(date +%s)"
echo "args: $@"
echo "tag: $(getjobad tag)"
echo "taskname: $(getjobad taskname)"
echo -e "\n--- end header output ---\n" #                       <----- section division


if [ -r "$OSGVO_CMSSW_Path"/cmsset_default.sh ]; then
    echo "sourcing environment: source $OSGVO_CMSSW_Path/cmsset_default.sh"
    source "$OSGVO_CMSSW_Path"/cmsset_default.sh
elif [ -r "$OSG_APP"/cmssoft/cms/cmsset_default.sh ]; then
    echo "sourcing environment: source $OSG_APP/cmssoft/cms/cmsset_default.sh"
    source "$OSG_APP"/cmssoft/cms/cmsset_default.sh
elif [ -r /cvmfs/cms.cern.ch/cmsset_default.sh ]; then
    echo "sourcing environment: source /cvmfs/cms.cern.ch/cmsset_default.sh"
    source /cvmfs/cms.cern.ch/cmsset_default.sh
else
    echo "ERROR! Couldn't find $OSGVO_CMSSW_Path/cmsset_default.sh or /cvmfs/cms.cern.ch/cmsset_default.sh or $OSG_APP/cmssoft/cms/cmsset_default.sh"
    exit 1
fi


# holy crap this is a mess. :( why does PAT code have to do such insane
# things with paths?
# if the first file in the tarball filelist starts with CMSSW, then it is
# a tarball made outside of the full CMSSW directory, and must be handled
# differently
if [ ! -z $(tar -tf ${TARFILE} | head -n 1 | grep "^CMSSW") ]; then
    echo "this is a full cmssw tar file"
    tar xf ${TARFILE}
    cd $CMSSWVERSION
    echo $PWD
    echo "Running ProjectRename"
    scramv1 b ProjectRename
    echo "Running `scramv1 runtime -sh`"
    eval `scramv1 runtime -sh`
    mv ../${TARFILE} .
else
    echo "this is a selective cmssw tar file"
    eval `scramv1 project CMSSW $CMSSWVERSION`
    cd $CMSSWVERSION
    eval `scramv1 runtime -sh`
    if [ -e ../${TARFILE} ]; then
        mv ../${TARFILE} ${TARFILE}
        tar xf ${TARFILE}
        rm ${TARFILE}
    fi
    echo "================================================================"
    echo "Before cleaning everything, need to check lib/${SCRAM_ARCH}"
    echo "================================================================"
    ls lib/${SCRAM_ARCH}
    cd src/cmstas/CORE
    make clean
    make -j 12
    cd ../../../
    scramv1 b clean
    echo "================================================================"
    echo "After cleaning everything, need to check lib/${SCRAM_ARCH} again"
    echo "================================================================"
    ls lib/${SCRAM_ARCH}
    scramv1 b -j 12
    echo "CMSSW_BASE: ${CMSSW_BASE}"
    #[ -e ${TARFILE} ] && tar xf ${TARFILE}
    # Needed or else cmssw can't find libmcfm_705.so
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/${SCRAM_ARCH}
    # This is nicer than above. both work, and both have scary but benign warnings/printouts
    # cp ${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/${SCRAM_ARCH}/*.so ${CMSSW_BASE}/lib/${SCRAM_ARCH}/
    # "Needed" to get rid of benign warnings/printouts
    export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${CMSSW_BASE}/src/ZZMatrixElement/MELA/interface
fi

# Go into the submission directory within $CMSSW_BASE/src
cd src/$SUBMIT_DIR


echo "before running: ls -lrth"
ls -lrth

echo -e "\n--- begin running ---\n" #                           <----- section division


##############
# ACTUAL RUN #
##############
# Transfer needs to be done through the script run.
runGenericROOTCommand.py --loadlib="loadLib.C" --script="$RUNFILE" --function="$FCN" --command="$EXTCMD" --recompile # Must force recompilation
RUN_STATUS=$?

echo "after running: ls -lrth"
ls -lrth

if [[ $RUN_STATUS != 0 ]]; then
    echo "Run has crashed with exit code ${RUN_STATUS}"
    exit 1
fi

echo -e "\n--- end running ---\n" #                             <----- section division

echo "time at end: $(date +%s)"
