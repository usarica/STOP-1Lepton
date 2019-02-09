#ifdef NOT_YET

#include <cassert>
#include "Samples.h"
#include "PhotonSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


float PhotonSelectionHelpers::photonCHEffArea_DR0p3(PhotonObject const& part){
  float aeta = fabs(part.eta());
  float EA = -999;

  if (SampleHelpers::theDataYear == 2016){
    // Spring16: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details
    if (aeta < 1.0) EA = 0.0360;
    else if (aeta < 1.479) EA = 0.0377;
    else if (aeta < 2.0) EA = 0.0306;
    else if (aeta < 2.2) EA = 0.0283;
    else if (aeta < 2.3) EA = 0.0254;
    else if (aeta < 2.4) EA = 0.0217;
    else EA = 0.0167;
  }
  else if (SampleHelpers::theDataYear == 2017 || SampleHelpers::theDataYear == 2018){
    // Fall17 V2: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
    if (aeta < 1.0) EA = 0.0112;
    else if (aeta < 1.479) EA = 0.0108;
    else if (aeta < 2.0) EA = 0.0106;
    else if (aeta < 2.2) EA = 0.01002;
    else if (aeta < 2.3) EA = 0.0098;
    else if (aeta < 2.4) EA = 0.0089;
    else EA = 0.0087;
  }

  return EA;
}
float PhotonSelectionHelpers::photonNHEffArea_DR0p3(PhotonObject const& part){
  float aeta = fabs(part.eta());
  float EA = -999;

  if (SampleHelpers::theDataYear == 2016){
    // Spring16: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details
    if (aeta < 1.0) EA = 0.0597;
    else if (aeta < 1.479) EA = 0.0807;
    else if (aeta < 2.0) EA = 0.0629;
    else if (aeta < 2.2) EA = 0.0197;
    else if (aeta < 2.3) EA = 0.0184;
    else if (aeta < 2.4) EA = 0.0284;
    else EA = 0.0591;
  }
  else if (SampleHelpers::theDataYear == 2017 || SampleHelpers::theDataYear == 2018){
    // Fall17 V2: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
    if (aeta < 1.0) EA = 0.0668;
    else if (aeta < 1.479) EA = 0.1054;
    else if (aeta < 2.0) EA = 0.0786;
    else if (aeta < 2.2) EA = 0.0233;
    else if (aeta < 2.3) EA = 0.0078;
    else if (aeta < 2.4) EA = 0.0028;
    else EA = 0.0137;
  }

  return EA;
}
float PhotonSelectionHelpers::photonEMEffArea_DR0p3(PhotonObject const& part){
  float aeta = fabs(part.eta());
  float EA = -999;

  if (SampleHelpers::theDataYear == 2016){
    // Spring16: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details
    if (aeta < 1.0) EA = 0.1210;
    else if (aeta < 1.479) EA = 0.1107;
    else if (aeta < 2.0) EA = 0.0699;
    else if (aeta < 2.2) EA = 0.1056;
    else if (aeta < 2.3) EA = 0.1457;
    else if (aeta < 2.4) EA = 0.1719;
    else EA = 0.1998;
  }
  else if (SampleHelpers::theDataYear == 2017 || SampleHelpers::theDataYear == 2018){
    // Fall17 V2: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
    if (aeta < 1.0) EA = 0.1113;
    else if (aeta < 1.479) EA = 0.0953;
    else if (aeta < 2.0) EA = 0.0619;
    else if (aeta < 2.2) EA = 0.0837;
    else if (aeta < 2.3) EA = 0.1070;
    else if (aeta < 2.4) EA = 0.1212;
    else EA = 0.1466;
  }

  return EA;
}

float PhotonSelectionHelpers::photonCHIsoCorr_DR0p3(PhotonObject const& part){
  float iso = part.extras.recoChargedHadronIso;
  float ea = photonCHEffArea_DR0p3(part);
  return std::max(float(0.0), iso - part.rho*ea);
}
float PhotonSelectionHelpers::photonNHIsoCorr_DR0p3(PhotonObject const& part){
  float iso = part.extras.recoNeutralHadronIso;
  float ea = photonNHEffArea_DR0p3(part);
  return std::max(float(0.0), iso - part.rho*ea);
}
float PhotonSelectionHelpers::photonEMIsoCorr_DR0p3(PhotonObject const& part){
  float iso = part.extras.recoPhotonIso;
  float ea = photonEMEffArea_DR0p3(part);
  return std::max(float(0.0), iso - part.rho*ea);
}

bool PhotonSelectionHelpers::testLooseCutBasedId(PhotonObject const& part){
  PhotonVariables const& extras = part.extras;

  float const pt  = part.pt();
  float const eta = part.eta();
  float const f_eta = fabs(eta);

  float const& hovere = extras.hOverE_full5x5;
  float const& sieie  = extras.sigmaIEtaIEta_full5x5;

  float const chiso = photonCHIsoCorr_DR0p3(part);
  float const nhiso = photonNHIsoCorr_DR0p3(part);
  float const emiso = photonEMIsoCorr_DR0p3(part);

  if (SampleHelpers::theDataYear == 2016){ // Same as CORE/PhotonSelections::isLoosePhotonPOG_Spring16
    if (f_eta < 1.479){
      if (hovere > 0.0597) return false;
      if (sieie > 0.01031) return false;
      if (chiso > 1.295) return false;
      if (nhiso > 10.910 + 0.0148*pt + 0.000017*pt*pt) return false;
      if (emiso > 3.630 + 0.0047*pt) return false;
    }
    else{
      if (hovere > 0.0481) return false;
      if (sieie > 0.03013) return false;
      if (chiso > 1.011) return false;
      if (nhiso > 5.931 + 0.0163*pt + 0.000014*pt*pt) return false;
      if (emiso > 6.641 + 0.0034*pt) return false;
    }
  }
  else if (SampleHelpers::theDataYear == 2017 || SampleHelpers::theDataYear == 2018){ // Same as CORE/PhotonSelections::isLoosePhotonPOG_Fall17V2
    if (f_eta < 1.479){
      if (hovere > 0.04596) return false;
      if (sieie > 0.0106) return false;
      if (chiso > 1.694) return false;
      if (nhiso > 24.032 + 0.01512*pt + 2.259e-5*pt*pt) return false;
      if (emiso > 2.876 + 0.004017*pt) return false;
    }
    else{
      if (hovere > 0.0590) return false;
      if (sieie > 0.0272) return false;
      if (chiso > 2.089) return false;
      if (nhiso > 19.722 + 0.0117*pt + 2.3e-5*pt*pt) return false;
      if (emiso > 4.162 + 0.0037*pt) return false;
    }
  }

  return true;
}
bool PhotonSelectionHelpers::testLooseSelection(PhotonObject const& part){
  // Id cut
  if (!testLooseCutBasedId(part)) return false;
  // Iso cut
  float reliso = miniRelIso_DR0p3(part);
  if (reliso>0.2) return false;
  return true;
}

bool PhotonSelectionHelpers::testMediumCutBasedId(PhotonObject const& part){
  PhotonVariables const& extras = part.extras;

  float const pt  = part.pt();
  float const eta = part.eta();
  float const f_eta = fabs(eta);

  float const& hovere = extras.hOverE_full5x5;
  float const& sieie  = extras.sigmaIEtaIEta_full5x5;

  float const chiso = photonCHIsoCorr_DR0p3(part);
  float const nhiso = photonNHIsoCorr_DR0p3(part);
  float const emiso = photonEMIsoCorr_DR0p3(part);

  if (SampleHelpers::theDataYear == 2016){ // Same as CORE/PhotonSelections::isMediumPhotonPOG_Spring16
    if (f_eta < 1.479){
      if (hovere > 0.0396) return false;
      if (sieie > 0.01022) return false;
      if (chiso > 0.441) return false;
      if (nhiso > 2.725 + 0.0148*pt + 0.000017*pt*pt) return false;
      if (emiso > 2.571 + 0.0047*pt) return false;
    }
    else{
      if (hovere > 0.0219) return false;
      if (sieie > 0.03001) return false;
      if (chiso > 0.442) return false;
      if (nhiso > 1.715 + 0.0163*pt + 0.000014*pt*pt) return false;
      if (emiso > 3.863 + 0.0034*pt) return false;
    }
  }
  else if (SampleHelpers::theDataYear == 2017 || SampleHelpers::theDataYear == 2018){ // Same as CORE/PhotonSelections::isMediumPhotonPOG_Fall17V2
    if (f_eta < 1.479){
      if (hovere > 0.02197) return false;
      if (sieie > 0.01015) return false;
      if (chiso > 1.141) return false;
      if (nhiso > 1.189 + 0.01512*pt + 2.259e-05*pt*pt) return false;
      if (emiso > 2.08 + 0.004017*pt) return false;
    }
    else{
      if (hovere > 0.0326) return false;
      if (sieie > 0.0272) return false;
      if (chiso > 1.051) return false;
      if (nhiso > 2.718 + 0.0117*pt + 2.3e-05*pt*pt) return false;
      if (emiso > 3.867 + 0.0037*pt) return false;
    }
  }

  return true;
}
bool PhotonSelectionHelpers::testMediumSelection(PhotonObject const& part){
  // Id cut
  if (!testMediumCutBasedId(part)) return false;
  // Iso cut
  float reliso = miniRelIso_DR0p3(part);
  if (reliso>0.1) return false;
  return true;
}

bool PhotonSelectionHelpers::testTightCutBasedId(PhotonObject const& part){
  PhotonVariables const& extras = part.extras;

  float const pt  = part.pt();
  float const eta = part.eta();
  float const f_eta = fabs(eta);

  float const& hovere = extras.hOverE_full5x5;
  float const& sieie  = extras.sigmaIEtaIEta_full5x5;

  float const chiso = photonCHIsoCorr_DR0p3(part);
  float const nhiso = photonNHIsoCorr_DR0p3(part);
  float const emiso = photonEMIsoCorr_DR0p3(part);

  if (SampleHelpers::theDataYear == 2016){ // Same as CORE/PhotonSelections::isTightPhotonPOG_Spring16
    if (f_eta < 1.479){
      if (hovere > 0.0269) return false;
      if (sieie > 0.00994) return false;
      if (chiso > 0.202) return false;
      if (nhiso > 0.264 + 0.0148*pt + 0.000017*pt*pt) return false;
      if (emiso > 2.362 + 0.0047*pt) return false;
    }
    else{
      if (hovere > 0.0213) return false;
      if (sieie > 0.03000) return false;
      if (chiso > 0.034) return false;
      if (nhiso > 0.586 + 0.0163*pt + 0.000014*pt*pt) return false;
      if (emiso > 2.617 + 0.0034*pt) return false;
    }
  }
  else if (SampleHelpers::theDataYear == 2017 || SampleHelpers::theDataYear == 2018){ // Same as CORE/PhotonSelections::isTightPhotonPOG_Fall17V2
    if (f_eta < 1.479){
      if (hovere > 0.02148) return false;
      if (sieie > 0.00996) return false;
      if (chiso > 0.65) return false;
      if (nhiso > 0.317 + 0.01512*pt + 2.259e-05*pt*pt) return false;
      if (emiso > 2.044 + 0.004017*pt) return false;
    }
    else{
      if (hovere > 0.0321) return false;
      if (sieie > 0.0271) return false;
      if (chiso > 0.517) return false;
      if (nhiso > 2.716 + 0.0117*pt + 2.3e-05*pt*pt) return false;
      if (emiso > 3.032 + 0.0037*pt) return false;
    }
  }

  return true;
}
bool PhotonSelectionHelpers::testTightSelection(PhotonObject const& part){
  // Id cut
  if (!testTightCutBasedId(part)) return false;
  // Iso cut
  float reliso = miniRelIso_DR0p3(part);
  if (reliso>0.1) return false;
  return true;
}

bool PhotonSelectionHelpers::testPtEtaGen(PhotonObject const& part){
  // pT and eta gen cut
  return (part.pt()>=ptThr_gen && fabs(part.eta())<etaThr_gen);
}
bool PhotonSelectionHelpers::testPtEtaSkim(PhotonObject const& part){
  // pT and eta skim cut
  if (
    (testTightSelection(part) && (part.pt()>=ptThr_skim_tight && fabs(part.eta())<etaThr_skim_tight))
    ||
    (testMediumSelection(part) && (part.pt()>=ptThr_skim_medium && fabs(part.eta())<etaThr_skim_medium))
    ||
    (testLooseSelection(part) && (part.pt()>=ptThr_skim_loose && fabs(part.eta())<etaThr_skim_loose))
    ||
    (testVetoSelection(part) && (part.pt()>=ptThr_skim_veto && fabs(part.eta())<etaThr_skim_veto))
    ) return true;
  return false;
}
bool PhotonSelectionHelpers::testPreselection(PhotonObject const& part){
  return (
    (
    (bit_preselection_idisoreco == kLooseIDReco && testLooseSelection(part))
      ||
      (bit_preselection_idisoreco == kMediumIDReco && testMediumSelection(part))
      ||
      (bit_preselection_idisoreco == kTightIDReco && testTightSelection(part))
      ) && testPtEtaSkim(part)
    );
}

void PhotonSelectionHelpers::setSelectionBits(PhotonObject& part){
  if (testPtEtaGen(part)) part.setSelectionBit(kGenPtEta);
  if (testLooseCutBasedId(part)) part.setSelectionBit(kLooseID);
  if (testLooseSelection(part)) part.setSelectionBit(kLooseIDReco);
  if (testMediumCutBasedId(part)) part.setSelectionBit(kMediumID);
  if (testMediumSelection(part)) part.setSelectionBit(kMediumIDReco);
  if (testTightCutBasedId(part)) part.setSelectionBit(kTightID);
  if (testTightSelection(part)) part.setSelectionBit(kTightIDReco);
  if (testPtEtaSkim(part)) part.setSelectionBit(kSkimPtEta);
  if (testPreselection(part)) part.setSelectionBit(kPreselection);
}

#endif // NOT_YET
