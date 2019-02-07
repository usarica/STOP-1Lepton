#include <cassert>
#include "Samples.h"
#include "HelperFunctions.h"
#include "AK4JetSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


bool AK4JetSelectionHelpers::isLooseAK4JetPOG(AK4JetObject const& part){
  if (SampleHelpers::theDataYear>2016) return true; // Loose id no longer needed in 2017 and 2018

  AK4JetVariables const& extras = part.extras;

  float uncorrE = (extras.undoJEC*part.energy());
  float nhf = extras.neutralHadronE / uncorrE;
  float nef = extras.neutralEmE / uncorrE;
  float chf = extras.chargedHadronE / uncorrE;
  float cef = extras.chargedEmE / uncorrE;
  float muf = extras.muonE / uncorrE;
  int const& cm = extras.chargedMultiplicity;
  int const& nm = extras.neutralMultiplicity;
  float eta = fabs(part.eta());

  if (cm + nm < 2) return false;
  if (nef >= 0.99) return false;
  if (nhf >= 0.99) return false;
  if (muf >= 0.8) return false;
  if (eta < 2.4){
    if (cm < 1) return false;
    if (!(chf > 0.f)) return false;
    if (cef >= 0.99) return false;
  }

  return true;
}
bool AK4JetSelectionHelpers::isTightAK4JetPOG(AK4JetObject const& part){
  if (!isLooseAK4JetPOG(part)) return false;

  AK4JetVariables const& extras = part.extras;

  float uncorrE = (extras.undoJEC*part.energy());
  float nhf = extras.neutralHadronE / uncorrE;
  float nef = extras.neutralEmE / uncorrE;
  float chf = extras.chargedHadronE / uncorrE;
  float cef = extras.chargedEmE / uncorrE;
  int const& ncands = extras.npfcands;
  int const& cm = extras.chargedMultiplicity;
  int const& nm = extras.neutralMultiplicity;
  float eta = fabs(part.eta());

  if (SampleHelpers::theDataYear==2016){
    if (nef >= 0.90) return false;
    if (nhf >= 0.90) return false;
    if (eta < 2.4 && cef >= 0.90) return false;
  }
  else if (SampleHelpers::theDataYear==2017){
    if (eta <= 2.4){
      if (chf == 0.00) return false;
      if (cm == 0) return false;
    }
    if (eta <= 2.7){
      if (nhf >= 0.90) return false;
      if (nef >= 0.90) return false;
      if (ncands <= 1) return false;
    }
    if (eta > 2.7 && eta <= 3.0){
      if (nef <= 0.02 || nef >= 0.99) return false;
      if (nm <= 2) return false;
    }
    if (eta > 3.0){
      if (nhf <= 0.02) return false;
      if (nef >= 0.90) return false;
      if (nm <= 10) return false;
    }
  }
  else if (SampleHelpers::theDataYear==2018){
    // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2018
    if (eta <= 2.6){
      if (nhf >= 0.90) return false;
      if (nef >= 0.90) return false;
      if (ncands <= 1) return false;
      if (chf <= 1e-6) return false;
      if (cm == 0) return false;
    }
    if (eta > 2.6 && eta <= 2.7){
      if (nhf >= 0.90) return false;
      if (nef >= 0.99) return false;
      if (cm == 0) return false;
    }
    if (eta > 2.7 && eta <= 3.0){
      if (nef <= 0.02 || nef >= 0.99) return false;
      if (nm <= 2) return false;
    }
    if (eta > 3.0){
      if (nhf <= 0.02) return false;
      if (nef >= 0.90) return false;
      if (nm <= 10) return false;
    }
  }

  return true;
}

bool AK4JetSelectionHelpers::testLooseId(AK4JetObject const& part){
  //AK4JetVariables const& extras = part.extras;
  if (!isLooseAK4JetPOG(part)) return false;
  return true;
}
bool AK4JetSelectionHelpers::testLooseSelection(AK4JetObject const& part){
  // Id cut
  if (!testLooseId(part)) return false;
  return true;
}

bool AK4JetSelectionHelpers::testTightId(AK4JetObject const& part){
  //AK4JetVariables const& extras = part.extras;
  if (!isTightAK4JetPOG(part)) return false;
  return true;
}
bool AK4JetSelectionHelpers::testTightSelection(AK4JetObject const& part){
  // Id cut
  if (!testTightId(part)) return false;
  return true;
}

bool AK4JetSelectionHelpers::isBadMuonJet(AK4JetObject& part, METObject const& obj){
  float dPhi = part.deltaPhi(obj.extras.phi);
  float uncorrE = (part.extras.undoJEC*part.energy());
  float muf = part.extras.muonE / uncorrE;
  bool res = (part.extras.undoJEC*part.pt()>200.f && dPhi>(TMath::Pi()-0.4f) && muf>0.5f);
  if (!res) part.setSelectionBit(kNotBadMuonJet);
  return res;
}

bool AK4JetSelectionHelpers::testSkimPtEta(AK4JetObject const& part, int icorr){
  CMSLorentzVector p4corr = part.getCorrectedMomentum(icorr);
  float pt = p4corr.Pt();
  float abs_eta = fabs(p4corr.Eta());
  return (pt>=ptThr_skim_preselection && abs_eta<etaThr_skim_preselection);
}
bool AK4JetSelectionHelpers::testPreselection(AK4JetObject const& part, int icorr){
  return (testSkimPtEta(part, icorr) && testTightSelection(part));
}

void AK4JetSelectionHelpers::setSelectionBits(AK4JetObject& part){
  if (testLooseId(part)) part.setSelectionBit(kLooseID);
  if (testTightId(part)) part.setSelectionBit(kTightID);

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
