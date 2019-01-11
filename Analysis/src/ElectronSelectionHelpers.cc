#include <cassert>
#include "Samples.h"
#include "ElectronSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


namespace ElectronSelectionHelpers{
  int eleEffAreaVersion = ElectronSelectionHelpers::setEleEffAreaVersion();
}

int ElectronSelectionHelpers::setEleEffAreaVersion(){
  // From CORE/Config.cc
  if (SampleHelpers::theDataYear == 2016) return 1;
  else if (SampleHelpers::theDataYear == 2017) return 4;
  else if (SampleHelpers::theDataYear == 2018) return 4; // FIXME: Needs new version
  else return -1;
}
float ElectronSelectionHelpers::electronEffArea(ElectronObject const& part){
  if (eleEffAreaVersion<0) eleEffAreaVersion = ElectronSelectionHelpers::setEleEffAreaVersion();

  float eta = part.eta();
  float const& etaSC = part.extras.etaSC;
  float ea = 0.;
  switch (eleEffAreaVersion){
  case 0:
  {
    //PHYS14 eleEffAreaVersion
    if (fabs(eta)<=0.800) ea = 0.1013;
    else if (fabs(eta)<=1.300) ea = 0.0988;
    else if (fabs(eta)<=2.000) ea = 0.0572;
    else if (fabs(eta)<=2.200) ea = 0.0842;
    else if (fabs(eta)<=2.500) ea = 0.1530;
    break;
  }
  case 1:
  {
    //Spring15 eleEffAreaVersion
    if (fabs(etaSC)<=1.000) ea = 0.1752;
    else if (fabs(etaSC)<=1.479) ea = 0.1862;
    else if (fabs(etaSC)<=2.000) ea = 0.1411;
    else if (fabs(etaSC)<=2.200) ea = 0.1534;
    else if (fabs(etaSC)<=2.300) ea = 0.1903;
    else if (fabs(etaSC)<=2.400) ea = 0.2243;
    else if (fabs(etaSC)<=2.500) ea = 0.2687;
    break;
  }
  case 2:
  {
    //Spring16 eleEffAreaVersion
    if (fabs(etaSC)<=1.000) ea = 0.1703;
    else if (fabs(etaSC)<=1.479) ea = 0.1715;
    else if (fabs(etaSC)<=2.000) ea = 0.1213;
    else if (fabs(etaSC)<=2.200) ea = 0.1230;
    else if (fabs(etaSC)<=2.300) ea = 0.1635;
    else if (fabs(etaSC)<=2.400) ea = 0.1937;
    else if (fabs(etaSC)<=2.500) ea = 0.2393;
    break;
  }
  case 3:
  {
    //Fall17 92X eleEffAreaVersion https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt
    if (fabs(etaSC)<=1.000) ea = 0.1566;
    else if (fabs(etaSC)<=1.479) ea = 0.1626;
    else if (fabs(etaSC)<=2.000) ea = 0.1073;
    else if (fabs(etaSC)<=2.200) ea = 0.0854;
    else if (fabs(etaSC)<=2.300) ea = 0.1051;
    else if (fabs(etaSC)<=2.400) ea = 0.1204;
    else if (fabs(etaSC)<=2.500) ea = 0.1524;
    break;
  }
  case 4:
  {
    //Fall17 94X eleEffAreaVersion https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
    if (fabs(etaSC)<=1.000) ea = 0.1440;
    else if (fabs(etaSC)<=1.479) ea = 0.1562;
    else if (fabs(etaSC)<=2.000) ea = 0.1032;
    else if (fabs(etaSC)<=2.200) ea = 0.0859;
    else if (fabs(etaSC)<=2.300) ea = 0.1116;
    else if (fabs(etaSC)<=2.400) ea = 0.1321;
    else if (fabs(etaSC)<=2.500) ea = 0.1654;
    break;
  }
  default:
  {
    MELAerr << "ElectronSelectionHelpers::electronEffArea: ERROR! eleEffAreaVersion=" << eleEffAreaVersion << " is unknown!" << endl;
    assert(0);
    break;
  }
  }
  return ea;
}

float ElectronSelectionHelpers::miniAbsIso_DR0p3(ElectronObject const& part){
  ElectronVariables const& extras = part.extras;
  float pt = part.pt();
  float dr = 0.2;
  if (pt>200.f) dr = 0.05;
  else if (pt>50.f) dr = 10.f/pt;
  dr /= 0.3f;
  float correction = extras.rho * electronEffArea(part) * pow(dr, 2);
  float res = extras.miniIso_ch + std::max(0.f, extras.miniIso_nh + extras.miniIso_em - correction);
  return res;
}
float ElectronSelectionHelpers::miniRelIso_DR0p3(ElectronObject const& part){ float pt = part.pt(); return (pt>0. ? miniAbsIso_DR0p3(part)/pt : 0.f); }

bool ElectronSelectionHelpers::testVetoCutBasedId(ElectronObject const& part){
  ElectronVariables const& extras = part.extras;

  // Id cuts
  if (extras.conv_vtx_flag) return false;
  if (SampleHelpers::theDataYear == 2016){ // Same as CORE/ElectronSelections::isVetoElectronPOGspring15noIso_v1
    if (fabs(extras.etaSC)>2.5) return false;
    else if (fabs(extras.etaSC) > 1.479){ // Endcap cuts
      if (extras.sigmaIEtaIEta_full5x5 >= 0.0352) return false;
      if (fabs(extras.dEtaIn) >= 0.0113) return false;
      if (fabs(extras.dPhiIn) >= 0.237) return false;
      if (extras.hOverE >= 0.116) return false;
      if (fabs(part.EinvMinusPinv()) >= 0.174) return false;
      if (fabs(extras.dxyPV) >= 0.222) return false;
      if (fabs(extras.dzPV) >= 0.921) return false;
      if (extras.expectedMissingInnerHits > 3) return false;
    }
    else/* if (fabs(extras.etaSC) <= 1.479)*/{ // Barrel cuts
      if (extras.sigmaIEtaIEta_full5x5 >= 0.0114) return false;
      if (fabs(extras.dEtaIn) >= 0.0152) return false;
      if (fabs(extras.dPhiIn) >= 0.216) return false;
      if (extras.hOverE >= 0.181) return false;
      if (fabs(part.EinvMinusPinv()) >= 0.207) return false;
      if (fabs(extras.dxyPV) >= 0.0564) return false;
      if (fabs(extras.dzPV) >= 0.472) return false;
      if (extras.expectedMissingInnerHits > 2) return false;
    }
  }
  else if (SampleHelpers::theDataYear == 2017 || SampleHelpers::theDataYear == 2018){ // Same as CORE/ElectronSelections::isVetoElectronPOGfall17noIso_v2. FIXME: 2018 will be updated!
    if (extras.conv_vtx_flag) return false;
    if (fabs(extras.etaSC)>2.5) return false;
    else if (fabs(extras.etaSC) > 1.479){ // Endcap cuts
      if (extras.sigmaIEtaIEta_full5x5 >= 0.0457) return false;
      if (fabs(extras.dEtaIn - extras.etaSC + extras.etaSeedSC) >= 0.00814) return false;
      if (fabs(extras.dPhiIn) >= 0.19) return false;
      if (extras.hOverE >= 0.05 + (2.54 + 0.183*extras.rho) / extras.energySC) return false;
      if (fabs(part.EinvMinusPinv()) >= 0.132) return false;
      if (extras.expectedMissingInnerHits > 3) return false;
    }
    else/* if (fabs(extras.etaSC) <= 1.479)*/{ // Barrel cuts
      if (extras.sigmaIEtaIEta_full5x5 >= 0.0126) return false;
      if (fabs(extras.dEtaIn - extras.etaSC + extras.etaSeedSC) >= 0.00463) return false;
      if (fabs(extras.dPhiIn) >= 0.148) return false;
      if (extras.hOverE >= 0.05 + (1.16 + 0.0324*extras.rho) / extras.energySC) return false;
      if (fabs(part.EinvMinusPinv()) >= 0.209) return false;
      if (extras.expectedMissingInnerHits > 2) return false;
    }
  }

  return true;
}
bool ElectronSelectionHelpers::testVetoSelection(ElectronObject const& part){
  // Id cut
  if (!testVetoCutBasedId(part)) return false;
  // Iso cut
  float reliso = miniRelIso_DR0p3(part);
  if (reliso>0.2) return false;
  return true;
}

bool ElectronSelectionHelpers::testLooseCutBasedId(ElectronObject const& part){
  ElectronVariables const& extras = part.extras;

  if (extras.conv_vtx_flag) return false;
  if (SampleHelpers::theDataYear == 2016){ // Same as CORE/ElectronSelections::isLooseElectronPOGspring15noIso_v1
    if (fabs(extras.etaSC)>2.5) return false;
    else if (fabs(extras.etaSC) > 1.479){ // Endcap cuts
      if (extras.sigmaIEtaIEta_full5x5 >= 0.0301) return false;
      if (fabs(extras.dEtaIn) >= 0.00814) return false;
      if (fabs(extras.dPhiIn) >= 0.182) return false;
      if (extras.hOverE >= 0.0897) return false;
      if (fabs(part.EinvMinusPinv()) >= 0.126) return false;
      if (fabs(extras.dxyPV) >= 0.118) return false;
      if (fabs(extras.dzPV) >= 0.822) return false;
      if (extras.expectedMissingInnerHits > 1) return false;
    }
    else/* if (fabs(extras.etaSC) <= 1.479)*/{ // Barrel cuts
      if (extras.sigmaIEtaIEta_full5x5 >= 0.0103) return false;
      if (fabs(extras.dEtaIn) >= 0.0105) return false;
      if (fabs(extras.dPhiIn) >= 0.115) return false;
      if (extras.hOverE >= 0.104) return false;
      if (fabs(part.EinvMinusPinv()) >= 0.102) return false;
      if (fabs(extras.dxyPV) >= 0.0261) return false;
      if (fabs(extras.dzPV) >= 0.41) return false;
      if (extras.expectedMissingInnerHits > 2) return false;
    }
  }
  else if (SampleHelpers::theDataYear == 2017 || SampleHelpers::theDataYear == 2018){ // Same as CORE/ElectronSelections::isLooseElectronPOGfall17noIso_v2. FIXME: 2018 will be updated!
    if (extras.conv_vtx_flag) return false;
    if (fabs(extras.etaSC)>2.5) return false;
    else if (fabs(extras.etaSC) > 1.479){ // Endcap cuts
      if (extras.sigmaIEtaIEta_full5x5 >= 0.0425) return false;
      if (fabs(extras.dEtaIn - extras.etaSC + extras.etaSeedSC) >= 0.00674) return false;
      if (fabs(extras.dPhiIn) >= 0.169) return false;
      if (extras.hOverE >= 0.0441 + (2.54 + 0.183*extras.rho) / extras.energySC) return false;
      if (fabs(part.EinvMinusPinv()) >= 0.111) return false;
      if (extras.expectedMissingInnerHits > 1) return false;
    }
    else/* if (fabs(extras.etaSC) <= 1.479)*/{ // Barrel cuts
      if (extras.sigmaIEtaIEta_full5x5 >= 0.0112) return false;
      if (fabs(extras.dEtaIn - extras.etaSC + extras.etaSeedSC) >= 0.00377) return false;
      if (fabs(extras.dPhiIn) >= 0.0884) return false;
      if (extras.hOverE >= 0.05 + (1.16 + 0.0324*extras.rho) / extras.energySC) return false;
      if (fabs(part.EinvMinusPinv()) >= 0.193) return false;
      if (extras.expectedMissingInnerHits > 1) return false;
    }
  }

  return true;
}
bool ElectronSelectionHelpers::testLooseSelection(ElectronObject const& part){
  // Id cut
  if (!testLooseCutBasedId(part)) return false;
  // Iso cut
  float reliso = miniRelIso_DR0p3(part);
  if (reliso>0.2) return false;
  return true;
}


bool ElectronSelectionHelpers::testMediumCutBasedId(ElectronObject const& part){
  ElectronVariables const& extras = part.extras;

  if (extras.conv_vtx_flag) return false;
  if (SampleHelpers::theDataYear == 2016){ // Same as CORE/ElectronSelections::isMediumElectronPOGspring15noIso_v1
    if (fabs(extras.etaSC)>2.5) return false;
    else if (fabs(extras.etaSC) > 1.479){ // Endcap cuts
      if (extras.sigmaIEtaIEta_full5x5 >= 0.0283) return false;
      if (fabs(extras.dEtaIn) >= 0.00733) return false;
      if (fabs(extras.dPhiIn) >= 0.114) return false;
      if (extras.hOverE >= 0.0678) return false;
      if (fabs(part.EinvMinusPinv()) >= 0.0898) return false;
      if (fabs(extras.dxyPV) >= 0.0739) return false;
      if (fabs(extras.dzPV) >= 0.602) return false;
      if (extras.expectedMissingInnerHits > 1) return false;
    }
    else/* if (fabs(extras.etaSC) <= 1.479)*/{ // Barrel cuts
      if (extras.sigmaIEtaIEta_full5x5 >= 0.0101) return false;
      if (fabs(extras.dEtaIn) >= 0.0103) return false;
      if (fabs(extras.dPhiIn) >= 0.0336) return false;
      if (extras.hOverE >= 0.0876) return false;
      if (fabs(part.EinvMinusPinv()) >= 0.0174) return false;
      if (fabs(extras.dxyPV) >= 0.0118) return false;
      if (fabs(extras.dzPV) >= 0.373) return false;
      if (extras.expectedMissingInnerHits > 2) return false;
    }
  }
  else if (SampleHelpers::theDataYear == 2017 || SampleHelpers::theDataYear == 2018){ // Same as CORE/ElectronSelections::isMediumElectronPOGfall17noIso_v2. FIXME: 2018 will be updated!
    if (extras.conv_vtx_flag) return false;
    if (fabs(extras.etaSC)>2.5) return false;
    else if (fabs(extras.etaSC) > 1.479){ // Endcap cuts
      if (extras.sigmaIEtaIEta_full5x5 >= 0.0387) return false;
      if (fabs(extras.dEtaIn - extras.etaSC + extras.etaSeedSC) >= 0.00632) return false;
      if (fabs(extras.dPhiIn) >= 0.0394) return false;
      if (extras.hOverE >= 0.0275 + (2.52 + 0.183*extras.rho) / extras.energySC) return false;
      if (fabs(part.EinvMinusPinv()) >= 0.0721) return false;
      if (extras.expectedMissingInnerHits > 1) return false;
    }
    else/* if (fabs(extras.etaSC) <= 1.479)*/{ // Barrel cuts
      if (extras.sigmaIEtaIEta_full5x5 >= 0.0106) return false;
      if (fabs(extras.dEtaIn - extras.etaSC + extras.etaSeedSC) >= 0.0032) return false;
      if (fabs(extras.dPhiIn) >= 0.0547) return false;
      if (extras.hOverE >= 0.046 + (1.16 + 0.0324*extras.rho) / extras.energySC) return false;
      if (fabs(part.EinvMinusPinv()) >= 0.184) return false;
      if (extras.expectedMissingInnerHits > 1) return false;
    }
  }

  return true;
}
bool ElectronSelectionHelpers::testMediumSelection(ElectronObject const& part){
  // Id cut
  if (!testMediumCutBasedId(part)) return false;
  // Iso cut
  float reliso = miniRelIso_DR0p3(part);
  if (reliso>0.1) return false;
  return true;
}

bool ElectronSelectionHelpers::testPtEtaGen(ElectronObject const& part){
  // pT and eta gen cut
  return (part.pt()>=ptThr_gen && fabs(part.eta())<etaThr_gen);
}
bool ElectronSelectionHelpers::testPtEtaSkim(ElectronObject const& part){
  // pT and eta skim cut
  if (
    (testMediumSelection(part) && (part.pt()>=ptThr_skim_medium && fabs(part.eta())<etaThr_skim_medium))
    ||
    (testLooseSelection(part) && (part.pt()>=ptThr_skim_loose && fabs(part.eta())<etaThr_skim_loose))
    ||
    (testVetoSelection(part) && (part.pt()>=ptThr_skim_veto && fabs(part.eta())<etaThr_skim_veto))
    ) return true;
  return false;
}
bool ElectronSelectionHelpers::testPreselection(ElectronObject const& part){ return (testMediumSelection(part) && testPtEtaSkim(part)); }

void ElectronSelectionHelpers::setSelectionBits(ElectronObject& part){
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
