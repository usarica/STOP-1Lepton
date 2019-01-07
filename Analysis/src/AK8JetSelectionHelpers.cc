#include <cassert>
#include "Samples.h"
#include "HelperFunctions.h"
#include "AK8JetSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


bool AK8JetSelectionHelpers::testSkimPtEta(AK8JetObject const& part, int icorr){
  CMSLorentzVector p4corr = part.getCorrectedMomentum(icorr);
  float pt = p4corr.Pt();
  float abs_eta = fabs(p4corr.Eta());
  return (pt>=ptThr_skim_preselection && abs_eta<etaThr_skim_preselection);
}
bool AK8JetSelectionHelpers::testPreselection(AK8JetObject const& part, int icorr){
  return (testSkimPtEta(part, icorr)/* && testTightSelection(part)*/);
}

void AK8JetSelectionHelpers::setSelectionBits(AK8JetObject& part){
  if (testSkimPtEta(part, 0)) part.setSelectionBit(kSkimPtEta);
  if (testPreselection(part, 0)) part.setSelectionBit(kPreselection);

  if (testSkimPtEta(part, 1)) part.setSelectionBit(kSkimPtEta_JECUp);
  if (testPreselection(part, 1)) part.setSelectionBit(kPreselection_JECUp);

  if (testSkimPtEta(part, -1)) part.setSelectionBit(kSkimPtEta_JECDn);
  if (testPreselection(part, -1)) part.setSelectionBit(kPreselection_JECDn);

  if (testSkimPtEta(part, 2)) part.setSelectionBit(kSkimPtEta_JERUp);
  if (testPreselection(part, 2)) part.setSelectionBit(kPreselection_JERUp);

  if (testSkimPtEta(part, -2)) part.setSelectionBit(kSkimPtEta_JERDn);
  if (testPreselection(part, -2)) part.setSelectionBit(kPreselection_JERDn);
}
