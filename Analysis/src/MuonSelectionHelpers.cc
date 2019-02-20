#include <cassert>
#include "Samples.h"
#include "HelperFunctions.h"
#include "MuonSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


namespace MuonSelectionHelpers{
  int muonEffAreaVersion = MuonSelectionHelpers::setMuonEffAreaVersion();
}

bool MuonSelectionHelpers::checkPOGSelectorBit(MuonObject const& part, POGSelectorBits ibit){ return test_bit(part.extras.POGSelectorBit, (unsigned int) ibit); }

int MuonSelectionHelpers::setMuonEffAreaVersion(){
  // From CORE/Config.cc
  if (SampleHelpers::theDataYear == 2016) return 1;
  else if (SampleHelpers::theDataYear == 2017) return 4;
  else if (SampleHelpers::theDataYear == 2018) return 4; // FIXME: Needs new version
  else return -1;
}
float MuonSelectionHelpers::muonEffArea(MuonObject const& part){
  if (muonEffAreaVersion<0) muonEffAreaVersion = MuonSelectionHelpers::setMuonEffAreaVersion();

  float eta = part.eta();
  float ea = 0.;
  switch (muonEffAreaVersion){
  case 0:
  {
    // PHYS14 muonEffAreaVersion
    if (fabs(eta)<=0.800) ea = 0.0913;
    else if (fabs(eta)<=1.300) ea = 0.0765;
    else if (fabs(eta)<=2.000) ea = 0.0546;
    else if (fabs(eta)<=2.200) ea = 0.0728;
    else if (fabs(eta)<=2.500) ea = 0.1177;
    break;
  }
  case 1:
  case 2:
  {
    // Spring15 muonEffAreaVersion
    // For version 2 it is same as Spring 15. Why? No idea.
    if (fabs(eta)<=0.800) ea = 0.0735;
    else if (fabs(eta)<=1.300) ea = 0.0619;
    else if (fabs(eta)<=2.000) ea = 0.0465;
    else if (fabs(eta)<=2.200) ea = 0.0433;
    else if (fabs(eta)<=2.500) ea = 0.0577;
    break;
  }
  case 3:
  case 4:
  {
    // Fall17 https://github.com/cms-data/PhysicsTools-NanoAOD/blob/master/effAreaMuons_cone03_pfNeuHadronsAndPhotons_94X.txt
    if (fabs(eta)<=0.800) ea = 0.0566;
    else if (fabs(eta)<=1.300) ea = 0.0562;
    else if (fabs(eta)<=2.000) ea = 0.0363;
    else if (fabs(eta)<=2.200) ea = 0.0119;
    else if (fabs(eta)<=2.500) ea = 0.0064;
    break;
  }
  default:
  {
    MELAerr << "MuonSelectionHelpers::muonEffArea: ERROR! muonEffAreaVersion=" << muonEffAreaVersion << " is unknown!" << endl;
    assert(0);
    break;
  }
  }
  return ea;
}

float MuonSelectionHelpers::miniAbsIso_DR0p3(MuonObject const& part){
  MuonVariables const& extras = part.extras;
  float pt = part.pt();
  float dr = 0.2;
  if (pt>200.f) dr = 0.05;
  else if (pt>50.f) dr = 10.f/pt;
  dr /= 0.3f;
  float correction = extras.rho * muonEffArea(part) * pow(dr, 2);
  float res = extras.miniIso_ch + std::max(0.f, extras.miniIso_nh + extras.miniIso_em - correction);
  return res;
}
float MuonSelectionHelpers::miniRelIso_DR0p3(MuonObject const& part){ float pt = part.pt(); return (pt>0. ? miniAbsIso_DR0p3(part)/pt : 0.f); }

bool MuonSelectionHelpers::isLooseMuonPOG(MuonObject const& part){
  MuonVariables const& extras = part.extras;

  if (!extras.isPFMuon) return false;
  bool isGlobal  = !(((extras.type) & (1<<1)) == 0);
  bool isTracker = !(((extras.type) & (1<<2)) == 0);

  return (isGlobal || isTracker);
}
bool MuonSelectionHelpers::isMediumMuonPOG(MuonObject const& part){
  if (!isLooseMuonPOG(part)) return false;

  MuonVariables const& extras = part.extras;

  bool isGlobal  = !(((extras.type) & (1<<1)) == 0);
  bool goodGlobal = (isGlobal && (extras.GlobalFit_Chisq/((float)extras.GlobalFit_Ndof))<3.f && extras.LocalPos_Chisq<12.f && extras.TrkKink<20.f);
  float validFraction = ((float) extras.validHits)/((float) (extras.validHits+extras.lostHits+extras.expectedMissingInnerHits+extras.expectedMissingOuterHits));

  return (validFraction > 0.8f && extras.SegComp >= (goodGlobal ? 0.303f : 0.451f));
}

bool MuonSelectionHelpers::testVetoCutBasedId(MuonObject const& part){ return testLooseCutBasedId(part); }
bool MuonSelectionHelpers::testVetoSelection(MuonObject const& part){ return testLooseSelection(part); }

bool MuonSelectionHelpers::testLooseCutBasedId(MuonObject const& part){
  MuonVariables const& extras = part.extras;

  if (SampleHelpers::theDataVersion>=SampleHelpers::kCMSSW_9_4_X){
    if (!checkPOGSelectorBit(part, CutBasedIdLoose)) return false;
  }
  else{
    if (!isLooseMuonPOG(part)) return false;
  }
  if (fabs(extras.dxyPV) > 0.1) return false;
  if (fabs(extras.dzPV) > 0.5) return false;

  return true;
}
bool MuonSelectionHelpers::testLooseSelection(MuonObject const& part){
  // Id cut
  if (!testLooseCutBasedId(part)) return false;
  // Iso cut
  float reliso = miniRelIso_DR0p3(part);
  if (reliso>0.2) return false;
  return true;
}


bool MuonSelectionHelpers::testMediumCutBasedId(MuonObject const& part){
  MuonVariables const& extras = part.extras;

  if (SampleHelpers::theDataVersion>=SampleHelpers::kCMSSW_9_4_X){
    if (!checkPOGSelectorBit(part, CutBasedIdMedium)) return false;
  }
  else{
    if (!isMediumMuonPOG(part)) return false;
  }
  if (fabs(extras.dxyPV) > 0.02) return false;
  if (fabs(extras.dzPV) > 0.1) return false;

  return true;
}
bool MuonSelectionHelpers::testMediumSelection(MuonObject const& part){
  // Id cut
  if (!testMediumCutBasedId(part)) return false;
  // Iso cut
  float reliso = miniRelIso_DR0p3(part);
  if (reliso>0.1) return false;
  return true;
}

bool MuonSelectionHelpers::testPtEtaGen(MuonObject const& part){
  // pT and eta gen cut
  return (part.pt()>=ptThr_gen && fabs(part.eta())<etaThr_gen);
}
bool MuonSelectionHelpers::testPtEtaSkim(MuonObject const& part){
  // pT and eta skim cut
  /*if (testTightSelection(part)) return (part.pt()>=ptThr_skim_tight && fabs(part.eta())<etaThr_skim_tight);
  else */if (testMediumSelection(part)) return (part.pt()>=ptThr_skim_medium && fabs(part.eta())<etaThr_skim_medium);
  else if (testLooseSelection(part)) return (part.pt()>=ptThr_skim_loose && fabs(part.eta())<etaThr_skim_loose);
  else if (testVetoSelection(part)) return (part.pt()>=ptThr_skim_veto && fabs(part.eta())<etaThr_skim_veto);
  else return false;
}
bool MuonSelectionHelpers::testPreselection(MuonObject const& part){
  return (
    (
    (bit_preselection_idisoreco == kVetoIDReco && testVetoSelection(part))
      ||
      (bit_preselection_idisoreco == kLooseIDReco && testLooseSelection(part))
      ||
      (bit_preselection_idisoreco == kMediumIDReco && testMediumSelection(part))
      //||
      //(bit_preselection_idisoreco == kTightIDReco && testTightSelection(part))
      ) && testPtEtaSkim(part)
    );
}

void MuonSelectionHelpers::setSelectionBits(MuonObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  if (testPtEtaGen(part)) part.setSelectionBit(kGenPtEta);
  if (testVetoCutBasedId(part)) part.setSelectionBit(kVetoID);
  if (testVetoSelection(part)) part.setSelectionBit(kVetoIDReco);
  if (testLooseCutBasedId(part)) part.setSelectionBit(kLooseID);
  if (testLooseSelection(part)) part.setSelectionBit(kLooseIDReco);
  if (testMediumCutBasedId(part)) part.setSelectionBit(kMediumID);
  if (testMediumSelection(part)) part.setSelectionBit(kMediumIDReco);
  //if (testTightCutBasedId(part)) part.setSelectionBit(kTightID);
  //if (testTightSelection(part)) part.setSelectionBit(kTightIDReco);
  if (testPtEtaSkim(part)) part.setSelectionBit(kSkimPtEta);
  if (testPreselection(part)) part.setSelectionBit(kPreselection);
}
