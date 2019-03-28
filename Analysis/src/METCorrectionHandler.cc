#include <cassert>
#include <string>
#include <fstream>
#include <sstream>
#include "METCorrectionHandler.h"
#include "MELAStreamHelpers.hh"
#include "TDirectory.h"
#include "TRandom3.h"


using namespace std;
using namespace SampleHelpers;
using namespace MELAStreamHelpers;


METCorrectionHandler::METCorrectionHandler() :
  ScaleFactorHandlerBase(),
  applyCorrection(false)
{
  setup();
}

METCorrectionHandler::~METCorrectionHandler(){ this->reset(); }


bool METCorrectionHandler::setup(){
  using namespace SystematicsHelpers;

  bool res = true;
  TDirectory* curdir = gDirectory;

  this->reset();

  TString strdatacore = Form("%s%s/%i/", STOP1LPKGDATAPATH.Data(), "METSFs", theDataYear);
  if (theDataPeriod == ""){
    MELAerr << "METCorrectionHandler::setup: Need to have a valid data period set." << endl;
    assert(0);
  }

  std::vector<TString> strValidPeriods = getValidDataPeriods();
  const size_t nValidPeriods = strValidPeriods.size();

  if (theDataYear == 2016){
    applyCorrection = true;
    // Use reserve since the push is done at the reading stage.
    lumilist.reserve(nValidPeriods);
    values_data_map.reserve(nValidPeriods);
    values_MC_map[sNominal] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[sNominal].reserve(nValidPeriods);
    values_MC_map[eJECDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJECDn].reserve(nValidPeriods);
    values_MC_map[eJECUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJECUp].reserve(nValidPeriods);
    values_MC_map[eJERDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJERDn].reserve(nValidPeriods);
    values_MC_map[eJERUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJERUp].reserve(nValidPeriods);
    values_MC_map[ePUDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[ePUDn].reserve(nValidPeriods);
    values_MC_map[ePUUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[ePUUp].reserve(nValidPeriods);
    for (TString const& period:strValidPeriods){
      if (theDataPeriod == "2016" || theDataPeriod == period){
        if (verbosity >= TVar::DEBUG) MELAout << "METCorrectionHandler::setup: Adding data period " << period << "..." << endl;
        if (period != "2016B") this->readFile(strdatacore+Form("fitparameters_Run%s-17Jul2018_MC.txt", period.Data()));
        else this->readFile(strdatacore+Form("fitparameters_Run%s-17Jul2018_ver2_MC.txt", period.Data()));
        lumilist.push_back(getIntegratedLuminosity(period));
      }
    }
    // Accumulate integrated luminosity and divide by the last one
    for (unsigned int il=1; il<lumilist.size(); il++) lumilist.at(il) += lumilist.at(il-1);
    for (unsigned int il=0; il<lumilist.size(); il++) lumilist.at(il) /= lumilist.back();
  }
  else if (theDataYear == 2017){
    applyCorrection = true;
    // Use reserve since the push is done at the reading stage.
    lumilist.reserve(nValidPeriods);
    values_data_map.reserve(nValidPeriods);
    values_MC_map[sNominal] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[sNominal].reserve(nValidPeriods);
    values_MC_map[eJECDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJECDn].reserve(nValidPeriods);
    values_MC_map[eJECUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJECUp].reserve(nValidPeriods);
    values_MC_map[eJERDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJERDn].reserve(nValidPeriods);
    values_MC_map[eJERUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJERUp].reserve(nValidPeriods);
    values_MC_map[ePUDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[ePUDn].reserve(nValidPeriods);
    values_MC_map[ePUUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[ePUUp].reserve(nValidPeriods);
    for (TString const& period:strValidPeriods){
      if (theDataPeriod == "2017" || theDataPeriod == period){
        if (verbosity >= TVar::DEBUG) MELAout << "METCorrectionHandler::setup: Adding data period " << period << "..." << endl;
        this->readFile(strdatacore+Form("fitparameters_Run%s-31Mar2018-v1_MC.txt", period.Data())); lumilist.push_back(getIntegratedLuminosity(period));
      }
    }
    if (theDataPeriod == "2017F-09May2018"){
      if (verbosity >= TVar::DEBUG) MELAout << "METCorrectionHandler::setup: Adding data period Run2017F-09May2018..." << endl;
      this->readFile(strdatacore+"fitparameters_Run2017F-09May2018-v1_MC.txt"); lumilist.push_back(getIntegratedLuminosity("2017F"));
    }
    // Accumulate integrated luminosity and divide by the last one
    for (unsigned int il=1; il<lumilist.size(); il++) lumilist.at(il) += lumilist.at(il-1);
    for (unsigned int il=0; il<lumilist.size(); il++) lumilist.at(il) /= lumilist.back();
  }
  else if (theDataYear == 2018){
    applyCorrection = true;
    // Use reserve since the push is done at the reading stage.
    lumilist.reserve(nValidPeriods);
    values_data_map.reserve(nValidPeriods);
    values_MC_map[sNominal] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[sNominal].reserve(nValidPeriods);
    values_MC_map[eJECDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJECDn].reserve(nValidPeriods);
    values_MC_map[eJECUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJECUp].reserve(nValidPeriods);
    values_MC_map[eJERDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJERDn].reserve(nValidPeriods);
    values_MC_map[eJERUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJERUp].reserve(nValidPeriods);
    values_MC_map[ePUDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[ePUDn].reserve(nValidPeriods);
    values_MC_map[ePUUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[ePUUp].reserve(nValidPeriods);
    for (TString const& period:strValidPeriods){
      if (theDataPeriod == "2018" || theDataPeriod == period){
        if (verbosity >= TVar::DEBUG) MELAout << "METCorrectionHandler::setup: Adding data period " << period << "..." << endl;
        if (period != "2018D") this->readFile(strdatacore+Form("fitparameters_Run%s-17Sep2018_MC.txt", period.Data()));
        else this->readFile(strdatacore+Form("fitparameters_Run%s-PromptReco_MC.txt", period.Data()));
        lumilist.push_back(getIntegratedLuminosity(period));
      }
    }
    // Accumulate integrated luminosity and divide by the last one
    for (unsigned int il=1; il<lumilist.size(); il++) lumilist.at(il) += lumilist.at(il-1);
    for (unsigned int il=0; il<lumilist.size(); il++) lumilist.at(il) /= lumilist.back();
  }

  curdir->cd();

  return res;
}
void METCorrectionHandler::reset(){
  applyCorrection = false;
  lumilist.clear();
  values_data_map.clear();
  values_MC_map.clear();
}

void METCorrectionHandler::readFile(TString const& strinput){
  using namespace SystematicsHelpers;

  // Notice that the systematics list does not contain eMETUp and eMETDn. This is because MET up and down variations involve varying the nominal values by the given errors.
  std::vector<SystematicVariationTypes> const validSysts={
    sNominal, eJECUp, eJECDn, eJERUp, eJERDn, ePUUp, ePUDn
  };
  const size_t nCopies = validSysts.size()+1;
  vector<pair<float, float>> value_error_S1; value_error_S1.reserve(nCopies);
  vector<pair<float, float>> value_error_S2; value_error_S2.reserve(nCopies);
  vector<pair<float, float>> value_error_S3; value_error_S3.reserve(nCopies);
  vector<pair<float, float>> value_error_F1; value_error_F1.reserve(nCopies);
  vector<pair<float, float>> value_error_F2; value_error_F2.reserve(nCopies);

  ifstream fin;
  fin.open(strinput.Data(), std::ios_base::in);
  if (fin.good()){
    string str_in;
    while (!fin.eof()){
      getline(fin, str_in);
      vector<pair<float, float>>* input_par_pair=nullptr;
      if (str_in.find("Sigma 1")!=string::npos) input_par_pair = &value_error_S1;
      else if (str_in.find("Sigma 2")!=string::npos) input_par_pair = &value_error_S2;
      else if (str_in.find("Sigma 3")!=string::npos) input_par_pair = &value_error_S3;
      else if (str_in.find("Frac 1")!=string::npos) input_par_pair = &value_error_F1;
      else if (str_in.find("Frac 2")!=string::npos) input_par_pair = &value_error_F2;
      if (input_par_pair){
        string dummy;
        float val, err;
        stringstream ss(str_in);
        //MELAout << str_in << endl;
        ss >> dummy >> dummy >> dummy >> val >> dummy >> err;
        input_par_pair->emplace_back(val, err);
      }
    }
  }
  fin.close();

  if (value_error_S1.size()!=nCopies || value_error_S2.size()!=nCopies || value_error_S3.size()!=nCopies || value_error_F1.size()!=nCopies || value_error_F2.size()!=nCopies){
    MELAerr << "METCorrectionHandler::readFile(" << strinput << "): Some values do not have size " << nCopies << endl;
    assert(0);
  }

  values_data_map.push_back(std::vector<std::pair<float, float>>());
  values_data_map.back().reserve(5);
  values_data_map.back().emplace_back(value_error_S1.front());
  values_data_map.back().emplace_back(value_error_S2.front());
  values_data_map.back().emplace_back(value_error_S3.front());
  values_data_map.back().emplace_back(value_error_F1.front());
  values_data_map.back().emplace_back(value_error_F2.front());

  size_t whichElement = 1;
  for (auto& syst:validSysts){
    auto* it = &(values_MC_map.find(syst)->second); it->push_back(std::vector<std::pair<float, float>>());
    it->back().reserve(5);
    it->back().emplace_back(value_error_S1.at(whichElement));
    it->back().emplace_back(value_error_S2.at(whichElement));
    it->back().emplace_back(value_error_S3.at(whichElement));
    it->back().emplace_back(value_error_F1.at(whichElement));
    it->back().emplace_back(value_error_F2.at(whichElement));
    whichElement++;
  }
}

void METCorrectionHandler::correctMET(float const& genMET, float const& genMETPhi, METObject* obj, bool /*useFastSim*/) const{
  using namespace SystematicsHelpers;

  if (!applyCorrection) return;
  if (!obj) return;

  std::vector<SystematicVariationTypes> const validSysts={
    sNominal, eMETDn, eMETUp, eJECDn, eJECUp, eJERDn, eJERUp, ePUDn, ePUUp
  };
  const TLorentzVector v_genMET(genMET*cos(genMETPhi), genMET*sin(genMETPhi), 0, genMET);

  TRandom3 rand;
  rand.SetSeed(std::abs(static_cast<int>(sin(genMETPhi)*100000)));
  float era_x = rand.Uniform();
  float frac_x = rand.Uniform();
  int i_era = -1;
  int i_frac = -1;
  for (auto const& lumi_era:lumilist){
    i_era++;
    if (era_x<=lumi_era) break;
  }
  if (i_era<0){
    MELAerr << "METCorrectionHandler::correctMET: i_era = " << i_era << " < 0." << endl;
    assert(0);
  }
  else if (i_era>=(int) values_data_map.size()){
    MELAerr << "METCorrectionHandler::correctMET: i_era = " << i_era << " >= size of eras = " << values_data_map.size() << endl;
    assert(0);
  }
  std::vector<std::pair<float, float>> const& data_values = values_data_map.at(i_era);
  float const& F1 = data_values.at(3).first;
  float const& F2 = data_values.at(4).first;
  if (frac_x<=F1) i_frac = 0;
  else if (frac_x<=(1.f - F1)*F2+F1) i_frac = 1;
  else i_frac = 2;

  if (verbosity>=TVar::DEBUG) MELAout << "METCorrectionHandler::correctMET: The i_era and i_frac random numbers: " << i_era << ", " << i_frac << endl;

  bool correctUpperError=false;
  float* met_value=nullptr;
  float* phi_value=nullptr;
  for (auto& syst:validSysts){
    auto systEff = syst;
    switch (syst){
    case sNominal:
      met_value = &(obj->extras.met);
      phi_value = &(obj->extras.phi);
      break;
    case eMETDn:
      systEff = sNominal;
      met_value = &(obj->extras.met_METdn);
      phi_value = &(obj->extras.phi_METdn);
      break;
    case eMETUp:
      systEff = sNominal;
      met_value = &(obj->extras.met_METup);
      phi_value = &(obj->extras.phi_METup);
      break;
    case eJECDn:
      met_value = &(obj->extras.met_JECdn);
      phi_value = &(obj->extras.phi_JECdn);
      break;
    case eJECUp:
      met_value = &(obj->extras.met_JECup);
      phi_value = &(obj->extras.phi_JECup);
      break;
    case eJERDn:
      met_value = &(obj->extras.met_JERdn);
      phi_value = &(obj->extras.phi_JERdn);
      break;
    case eJERUp:
      met_value = &(obj->extras.met_JERup);
      phi_value = &(obj->extras.phi_JERup);
      break;
    case ePUDn:
      met_value = &(obj->extras.met_PUdn);
      phi_value = &(obj->extras.phi_PUdn);
      break;
    case ePUUp:
      met_value = &(obj->extras.met_PUup);
      phi_value = &(obj->extras.phi_PUup);
      break;
    default:
      MELAout << "METCorrectionHandler::correctMET: Unsupported systematic " << syst << " in switch statement." << endl;
      assert(0);
      break;
    }

    TLorentzVector v_MET(
      *met_value * cos(*phi_value),
      *met_value * sin(*phi_value),
      0,
      *met_value
    );
    if (verbosity>=TVar::DEBUG) MELAout << "METCorrectionHandler::correctMET: The old MET vector: " << v_MET << endl;

    v_MET = v_MET - v_genMET;
    //v_MET = v_MET - TLorentzVector(genMET*cos(*phi_value), genMET*sin(*phi_value), 0, genMET);
    if (verbosity>=TVar::DEBUG) MELAout << "METCorrectionHandler::correctMET: The MET difference vector: " << v_MET << endl;

    auto const& mc_values = values_MC_map.find(systEff)->second;
    float data_sigma = data_values.at(i_frac).first;
    float mc_sigma = mc_values.at(i_era).at(i_frac).first;
    if (syst == eMETDn){
      data_sigma -= data_values.at(i_frac).second;
      mc_sigma += mc_values.at(i_era).at(i_frac).second;
    }
    else if (syst == eMETUp){
      data_sigma += data_values.at(i_frac).second;
      mc_sigma -= mc_values.at(i_era).at(i_frac).second;
    }

    float scale = data_sigma/mc_sigma;
    if (syst == sNominal) correctUpperError = (scale<1.f);
    if (correctUpperError && syst == eMETUp && scale<1.f) scale = 1.f;

    if (verbosity>=TVar::DEBUG) MELAout << "METCorrectionHandler::correctMET: The MET scale for systematic " << syst << ": " << data_sigma << " / " << mc_sigma << " = " << scale << endl;

    v_MET *= scale;
    v_MET += v_genMET;
    *met_value = v_MET.Pt();
    *phi_value = v_MET.Phi();

    if (verbosity>=TVar::DEBUG) MELAout << "METCorrectionHandler::correctMET: The new MET vector: " << v_MET << endl;
  }
}

void METCorrectionHandler::printParameters() const{
  using namespace SystematicsHelpers;

  if (!applyCorrection) return;

  MELAout << "Begin METCorrectionHandler::printParameters" << endl;
  MELAout << "Normalized lumi list: " << lumilist << endl;
  unsigned int i_era = 0;
  for (auto const& ed:values_data_map){
    MELAout << "Era " << i_era << " values in data: " << ed << endl;
    i_era++;
  }

  std::vector<SystematicVariationTypes> const validSysts={
    sNominal, eJECDn, eJECUp, eJERDn, eJERUp, ePUDn, ePUUp
  };
  for (auto const& syst:validSysts){
    i_era=0;
    for (auto const& ed:values_MC_map.find(syst)->second){
      MELAout << "Era " << i_era << " values in MC[" << syst << "]: " << ed << endl;
      i_era++;
    }
  }
  MELAout << "End METCorrectionHandler::printParameters" << endl;
}
