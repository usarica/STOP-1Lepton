#!/bin/sh


TARFILE="stop_1lepton.tar"
echo "SCRAM_ARCH: ${SCRAM_ARCH}"

HERE=$(pwd)

pushd $CMSSW_BASE

# Check metis/UserTarball for which directories to include from CMSSW_BASE
# Or just use it
tar Jcvf ${TARFILE} \
lib \
biglib \
cfipython \
config \
external \
bin \
src/ZZMatrixElement \
src/HiggsAnalysis \
src/TopTagger \
src/CMSDataTools \
src/cmstas \
src/cms-jet/JRDatabase \
src/STOP_1Lepton \
--exclude=lib/${SCRAM_ARCH}/* \
--exclude=src/ZZMatrixElement/MELA/COLLIER/*.so \
--exclude=src/ZZMatrixElement/MELA/data/Pdfdata/NNPDF30_lo_as_0130.LHgrid \
--exclude=src/ZZMatrixElement/MELA/test/reference \
--exclude=src/HiggsAnalysis/CombinedLimit/data \
--exclude=src/cmstas/CORE/*.o \
--exclude=src/cmstas/CORE/*.so \
--exclude=src/STOP_1Lepton/Analysis/test/output \
--exclude=src/STOP_1Lepton/Analysis/test/*.root \
--exclude=src/STOP_1Lepton/Analysis/test/*.pcm \
--exclude=src/STOP_1Lepton/Analysis/test/*.so \
--exclude=src/STOP_1Lepton/Analysis/test/*.a \
--exclude=src/STOP_1Lepton/Analysis/test/*.o \
--exclude=src/STOP_1Lepton/Analysis/test/*.d \
--exclude={.git,.gitignore,*.tar,libmcfm*}

mv stop_1lepton.tar $HERE/

popd
