#include <cassert>
#include <cstdio>
#include "JECScaleFactorHandler.h"
#include "MELAStreamHelpers.hh"
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>


using namespace std;
using namespace SampleHelpers;
using namespace JECJERHelpers;
using namespace MELAStreamHelpers;


JECScaleFactorHandler::JECScaleFactorHandler(JECJERHelpers::JECJERType type_) :
  ScaleFactorHandlerBase(),
  type(type_),
  corrector_data(nullptr),
  corrector_MC_noFS(nullptr),
  corrector_MC_FS(nullptr),
  uncertaintyEstimator_MC_noFS(nullptr),
  uncertaintyEstimator_MC_FS(nullptr)
{
  setup();
}

JECScaleFactorHandler::~JECScaleFactorHandler(){ this->reset(); }

FactorizedJetCorrector* JECScaleFactorHandler::makeCorrector(std::vector<TString> const& fnames){
  if (fnames.empty()) return nullptr;

  for (auto const& fname:fnames){
    if (!HostHelpers::FileExists(fname)){
      MELAerr << "JECScaleFactorHandler::makeCorrector: File " << fname << " does not exists! Aborting..." << endl;
      assert(0);
    }
  }

  std::vector<JetCorrectorParameters> vParam; vParam.reserve(fnames.size());
  for (auto const& fname:fnames){
    const TString cmd = "echo ";
    FILE* f = popen((cmd + fname).Data(), "r");
    if (!f){
      MELAerr << "JECScaleFactorHandler::makeCorrector: Error opening pipe to execute " << cmd << fname << endl;
      return nullptr;
    }
    char corr_name[1024];
    int s = fscanf(f, " %1024s\n", corr_name);
    if (s != 1){
      MELAerr << "JECScaleFactorHandler::makeCorrector: Error reading file " << fname << endl;
      assert(0);
    }
    vParam.emplace_back(corr_name);
    pclose(f);
  }

  return new FactorizedJetCorrector(vParam);
}
JetCorrectionUncertainty* JECScaleFactorHandler::makeUncertaintyEstimator(TString const& fname){ return new JetCorrectionUncertainty(fname.Data()); }


bool JECScaleFactorHandler::setup(){
  bool res = true;
  this->reset();

  TDirectory* curdir = gDirectory;

  std::vector<TString> correctornames_data = getJECFileNames(type, false, false);
  corrector_data = makeCorrector(correctornames_data);

  std::vector<TString> correctornames_MC_noFS = getJECFileNames(type, true, false);
  corrector_MC_noFS = makeCorrector(correctornames_MC_noFS);

  std::vector<TString> correctornames_MC_FS = getJECFileNames(type, true, true);
  corrector_MC_FS = makeCorrector(correctornames_MC_FS);

  TString uncname_MC_noFS = getJECUncertaintyFileName(type, true, false);
  uncertaintyEstimator_MC_noFS = makeUncertaintyEstimator(uncname_MC_noFS);

  TString uncname_MC_FS = getJECUncertaintyFileName(type, true, true);
  uncertaintyEstimator_MC_FS = makeUncertaintyEstimator(uncname_MC_FS);

  curdir->cd();

  return res;
}
void JECScaleFactorHandler::reset(){
  delete corrector_data; corrector_data = nullptr;
  delete corrector_MC_noFS; corrector_MC_noFS = nullptr;
  delete corrector_MC_FS; corrector_MC_FS = nullptr;
  delete uncertaintyEstimator_MC_noFS; uncertaintyEstimator_MC_noFS = nullptr;
  delete uncertaintyEstimator_MC_FS; uncertaintyEstimator_MC_FS = nullptr;
}

void JECScaleFactorHandler::applyJEC(AK4JetObject* obj, bool isMC, bool isFastSim){
  if (type!=kAK4){
    MELAerr << "JECScaleFactorHandler::applyJEC: Mismatch between AK4 jets and type " << type << "! Aborting..." << endl;
    assert(0);
  }

  FactorizedJetCorrector* corrector = nullptr;
  if (!isMC) corrector = corrector_data;
  else if (!isFastSim) corrector = corrector_MC_noFS;
  else corrector = corrector_MC_FS;

  if (!corrector) return;

  JetCorrectionUncertainty* uncEst = nullptr;
  if (isMC){
    if (!isFastSim) uncEst = uncertaintyEstimator_MC_noFS;
    else uncEst = uncertaintyEstimator_MC_FS;
  }

  // Get nominal JEC
  float JEC = 1;
  double jeta = obj->eta();
  double jpt_unc = obj->extras.undoJEC*obj->momentum.Pt();
  if (corrector){
    corrector->setRho(obj->extras.rho);
    corrector->setJetA(obj->extras.area);
    corrector->setJetPt(jpt_unc);
    corrector->setJetEta(jeta);
    JEC = corrector->getCorrection();
  }
  if (JEC<0.f){
    MELAerr << "JECScaleFactorHandler::setJECs: Nominal JEC value " << JEC << "< 0!" << endl;
    JEC=0.;
  }

  // Get up/dn variations
  float JECunc = 0;
  if (uncEst){
    uncEst->setJetEta(jeta);
    uncEst->setJetPt(jpt_unc*JEC); // Must use corrected pT
    JECunc = uncEst->getUncertainty(true);
  }
  float JECup = JEC + JECunc;
  float JECdn = std::max(0.f, JEC - JECunc);

  obj->extras.JEC = JEC*obj->extras.undoJEC;
  obj->extras.JECup = JECup*obj->extras.undoJEC;
  obj->extras.JECdn = JECdn*obj->extras.undoJEC;
}
void JECScaleFactorHandler::applyJEC(AK8JetObject* obj, bool isMC, bool isFastSim){
  if (type!=kAK8){
    MELAerr << "JECScaleFactorHandler::applyJEC: Mismatch between AK8 jets and type " << type << "! Aborting..." << endl;
    assert(0);
  }

  FactorizedJetCorrector* corrector = nullptr;
  if (!isMC) corrector = corrector_data;
  else if (!isFastSim) corrector = corrector_MC_noFS;
  else corrector = corrector_MC_FS;

  if (!corrector) return;

  JetCorrectionUncertainty* uncEst = nullptr;
  if (isMC){
    if (!isFastSim) uncEst = uncertaintyEstimator_MC_noFS;
    else uncEst = uncertaintyEstimator_MC_FS;
  }

  // Get nominal JEC
  float JEC = 1;
  double jeta = obj->eta();
  double jpt_unc = obj->extras.undoJEC*obj->momentum.Pt();
  if (corrector){
    corrector->setRho(obj->extras.rho);
    corrector->setJetA(obj->extras.area);
    corrector->setJetPt(jpt_unc);
    corrector->setJetEta(jeta);
    JEC = corrector->getCorrection();
  }
  if (JEC<0.f){
    MELAerr << "JECScaleFactorHandler::setJECs: Nominal JEC value " << JEC << "< 0!" << endl;
    JEC=0.;
  }

  // Get up/dn variations
  float JECunc = 0;
  if (uncEst){
    uncEst->setJetEta(jeta);
    uncEst->setJetPt(jpt_unc*JEC); // Must use corrected pT
    JECunc = uncEst->getUncertainty(true);
  }
  float JECup = JEC + JECunc;
  float JECdn = std::max(0.f, JEC - JECunc);

  obj->extras.JEC = JEC*obj->extras.undoJEC;
  obj->extras.JECup = JECup*obj->extras.undoJEC;
  obj->extras.JECdn = JECdn*obj->extras.undoJEC;
}
