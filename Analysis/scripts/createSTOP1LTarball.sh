#!/bin/sh


HERE=$(pwd)

pushd $CMSSW_BASE

tar -cvf stop_1lepton.tar \
lib \
bin \
src/ZZMatrixElement \
src/HiggsAnalysis \
src/TopTagger \
src/CMSDataTools \
src/cmstas \
src/cms-jet \
src/STOP_1Lepton \
--exclude=src/ZZMatrixElement/MELA/test/reference \
--exclude=src/HiggsAnalysis/CombinedLimit/data \
--exclude=src/STOP_1Lepton/Analysis/test/output \
--exclude={.git,.gitignore,*.tar}

mv stop_1lepton.tar $HERE/

popd
