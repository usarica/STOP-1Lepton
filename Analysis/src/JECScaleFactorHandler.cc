#include <cassert>
#include <cstdio>
#include "JECScaleFactorHandler.h"
#include "MELAStreamHelpers.hh"
#include <cmstas/CORE/Tools/jetcorr/JetCorrectorParameters.h>


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

  for (TString const& fname:fnames){
    if (!HostHelpers::FileExists(fname)){
      MELAerr << "JECScaleFactorHandler::makeCorrector: File " << fname << " does not exists! Aborting..." << endl;
      assert(0);
    }
  }

  std::vector<JetCorrectorParameters> vParam; vParam.reserve(fnames.size());
  for (TString const& fname:fnames) vParam.emplace_back(fname.Data());

  return new FactorizedJetCorrector(vParam);
}
JetCorrectionUncertainty* JECScaleFactorHandler::makeUncertaintyEstimator(TString const& fname){ return new JetCorrectionUncertainty(fname.Data()); }


bool JECScaleFactorHandler::setup(){
  bool res = true;
  this->reset();

  TDirectory* curdir = gDirectory;

  std::vector<TString> correctornames_data = getJECFileNames(type, false, false);
  if (!correctornames_data.empty()) corrector_data = makeCorrector(correctornames_data);

  std::vector<TString> correctornames_MC_noFS = getJECFileNames(type, true, false);
  if (!correctornames_MC_noFS.empty()) corrector_MC_noFS = makeCorrector(correctornames_MC_noFS);

  std::vector<TString> correctornames_MC_FS = getJECFileNames(type, true, true);
  if (!correctornames_MC_FS.empty()) corrector_MC_FS = makeCorrector(correctornames_MC_FS);

  TString uncname_MC_noFS = getJECUncertaintyFileName(type, true, false);
  if (uncname_MC_noFS!="") uncertaintyEstimator_MC_noFS = makeUncertaintyEstimator(uncname_MC_noFS);

  TString uncname_MC_FS = getJECUncertaintyFileName(type, true, true);
  if (uncname_MC_FS!="") uncertaintyEstimator_MC_FS = makeUncertaintyEstimator(uncname_MC_FS);

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

  // Unfortunately, we also need to record the version of JEC without the muons, so we need to create new vectors
  CMSLorentzVector const jet_p4_uncor = obj->momentum*obj->extras.undoJEC;
  CMSLorentzVector const& jet_p4_nomus_uncor = obj->extras.momentum_nomus_uncor;

  // Get nominal JECs
  float JEC = 1;
  double jpt_unc = jet_p4_uncor.Pt();
  double jeta_unc = jet_p4_uncor.eta();

  float JEC_nomus = 1; // L3 in MC, L3 or L3+residuals in data
  float JEC_prev_nomus = 1;
  float JEC_L1_nomus = 1; // We need to store this when computing MET corrections according to the L1L2L3-L1 recipe
  size_t size_of_corrections_nomus = 0; // 3 or 4 hopefully
  double jpt_unc_nomus = jet_p4_nomus_uncor.Pt();
  double jeta_unc_nomus = jet_p4_nomus_uncor.eta();
  if (corrector){
    corrector->setRho(obj->extras.rho);
    corrector->setJetA(obj->extras.area);
    corrector->setJetPt(jpt_unc);
    corrector->setJetEta(jeta_unc);
    JEC = corrector->getCorrection();

    corrector->setRho(obj->extras.rho);
    corrector->setJetA(obj->extras.area);
    corrector->setJetPt(jpt_unc_nomus);
    corrector->setJetEta(jeta_unc_nomus);
    std::vector<float> const corr_vals_nomus = corrector->getSubCorrections(); // Subcorrections are stored with corr_vals(N) = corr(N)*corr(N-1)*...*corr(1)
    JEC_nomus = corr_vals_nomus.back();
    JEC_L1_nomus = corr_vals_nomus.front();
    size_of_corrections_nomus = corr_vals_nomus.size();
    if (!isMC && size_of_corrections_nomus==4) JEC_prev_nomus = corr_vals_nomus.at(corr_vals_nomus.size()-2); // We will use this L3 value in the uncertainty estimation below.
  }

  if (JEC<0.f){
    MELAerr << "JECScaleFactorHandler::setJECs: Nominal JEC value " << JEC << "< 0!" << endl;
    JEC=0.;
  }
  if (JEC_nomus<0.f){
    MELAerr << "JECScaleFactorHandler::setJECs: Nominal JEC_nomus value " << JEC_nomus << "< 0!" << endl;
    JEC_nomus=0.;
  }
  if (JEC_L1_nomus<0.f){
    MELAerr << "JECScaleFactorHandler::setJECs: Nominal JEC_L1_nomus value " << JEC_L1_nomus << "< 0!" << endl;
    JEC_L1_nomus=0.;
  }

  // Get up/dn variations
  float JECunc = 0;
  float JECunc_nomus = 0;
  if (uncEst){
    uncEst->setJetPt(jpt_unc*JEC); // Must use corrected pT
    uncEst->setJetEta(jeta_unc);
    JECunc = uncEst->getUncertainty(true);

    uncEst->setJetPt(jpt_unc_nomus*JEC_nomus);
    uncEst->setJetEta(jeta_unc_nomus);
    JECunc_nomus = uncEst->getUncertainty(true);
    if (!isMC && size_of_corrections_nomus==4) JECunc_nomus = sqrt(pow(JECunc_nomus, 2) + pow((JEC_nomus/JEC_prev_nomus-1.f), 2));
  }
  float JECup = JEC + JECunc;
  float JECdn = std::max(0.f, JEC - JECunc);

  obj->extras.JEC = JEC*obj->extras.undoJEC;
  obj->extras.JECup = JECup*obj->extras.undoJEC;
  obj->extras.JECdn = JECdn*obj->extras.undoJEC;

  obj->extras.JEC_raw_nomus = JEC_nomus;
  obj->extras.JEC_L1_raw_nomus = JEC_L1_nomus;
  obj->extras.JEC_raw_unc_nomus = JECunc_nomus;
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
