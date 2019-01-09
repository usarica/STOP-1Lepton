#include <cassert>
#include "BtagScaleFactorHandler.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace BtagHelpers;
using namespace MELAStreamHelpers;


BtagScaleFactorHandler::BtagScaleFactorHandler(BtagHelpers::BtagWPType type_, bool isFastSim_) :
  ScaleFactorHandlerBase(),
  type(type_),
  isFastSim(isFastSim_),
  m_calib(nullptr),
  m_fileEff(nullptr),
  WPval(-1)
{
  setup();
}

BtagScaleFactorHandler::~BtagScaleFactorHandler(){ this->reset(); }


bool BtagScaleFactorHandler::setup(){
  bool res = true;
  this->reset();

  TDirectory* curdir = gDirectory;

  WPval = BtagHelpers::getBtagWP(type);
  TString sf_fname = BtagHelpers::getBtagSFFileName(type, isFastSim);
  TString eff_fname = BtagHelpers::getBtagEffFileName(type, isFastSim);

  BTagEntry::OperatingPoint opPoint;
  TString calibname;
  switch (type){
  case kCSVv2_Loose:
  case kCSVv2_Medium:
  case kCSVv2_Tight:
    calibname = (isFastSim ? "deepcsv" : "DeepCSV");
    break;
  case kDeepCSV_Loose:
  case kDeepCSV_Medium:
  case kDeepCSV_Tight:
    calibname = "CSVv2";
    break;
  default:
    MELAerr << "BtagScaleFactorHandler::setup: No implementation for b tag WP " << type << ". Aborting..." << endl;
    assert(0);
    break;
  }
  switch (type){
  case kCSVv2_Loose:
  case kDeepCSV_Loose:
    opPoint = BTagEntry::OP_LOOSE;
    break;
  case kCSVv2_Medium:
  case kDeepCSV_Medium:
    opPoint = BTagEntry::OP_MEDIUM;
    break;
  case kCSVv2_Tight:
  case kDeepCSV_Tight:
    opPoint = BTagEntry::OP_TIGHT;
    break;
  default:
    MELAerr << "BtagScaleFactorHandler::setup: No implementation for b tag WP " << type << ". Aborting..." << endl;
    assert(0);
    break;
  }
  m_calib = new BTagCalibration(calibname.Data(), sf_fname.Data());
  BTagCalibrationReader* m_reader = new BTagCalibrationReader(opPoint, "central");
  BTagCalibrationReader* m_reader_up = new BTagCalibrationReader(opPoint, "up");
  BTagCalibrationReader* m_reader_down = new BTagCalibrationReader(opPoint, "down");
  m_readers = std::vector<BTagCalibrationReader*>{ m_reader, m_reader_up, m_reader_down };
  for (BTagCalibrationReader*& mr:m_readers){
    mr->load(*m_calib, BTagEntry::FLAV_B, "comb");
    mr->load(*m_calib, BTagEntry::FLAV_C, "comb");
    mr->load(*m_calib, BTagEntry::FLAV_UDSG, "incl");
  }

  std::vector<TString> hnames = BtagHelpers::getBtagEffHistogramNames(type, isFastSim);
  assert(hnames.size()==3);
  m_fileEff = TFile::Open(eff_fname, "read");
  for (unsigned int i=0; i<3; i++) m_hEff[i] = (TH1*) m_fileEff->Get(hnames.at(i));
  curdir->cd();

  return res;
}
void BtagScaleFactorHandler::reset(){
  ScaleFactorHandlerBase::closeFile(m_fileEff);
  for (auto*& m_reader:m_readers) delete m_reader;
  m_readers.clear();
  delete m_calib; m_calib = nullptr;
  WPval = -1;
}

float BtagScaleFactorHandler::getSF(int syst, int jetFlavor, float pt, float eta){
  float SF = 1.0;

  unsigned int isyst=0;
  if (syst==1) isyst=1;
  else if (syst==-1) isyst=2;

  BTagEntry::JetFlavor flav;
  if (abs(jetFlavor)==5) flav = BTagEntry::FLAV_B;
  else if (abs(jetFlavor)==4) flav = BTagEntry::FLAV_C;
  else flav = BTagEntry::FLAV_UDSG;

  float myPt = pt;
  float MaxJetEta = 2.5; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
  if (std::abs(eta) > MaxJetEta) return SF; // Do not apply SF for jets with eta higher than the threshold

  std::pair<float, float> pt_min_max = m_readers[isyst]->min_max_pt(flav, eta);
  if (pt_min_max.second<0.) return SF;
  bool DoubleUncertainty = false;
  if (myPt<pt_min_max.first){
    myPt = pt_min_max.first+1e-5;
    DoubleUncertainty = true;
  }
  else if (myPt>pt_min_max.second){
    myPt = pt_min_max.second-1e-5;
    DoubleUncertainty = true;
  }

  SF = m_readers[isyst]->eval(flav, eta, myPt);
  if (DoubleUncertainty && syst!=0){
    float SFcentral = m_readers[0]->eval(flav, eta, myPt);
    SF = 2.f*(SF - SFcentral) + SFcentral;
  }

  return SF;
}
float BtagScaleFactorHandler::getEff(int jetFlavor, float pt, float eta){
  int flav;
  if (abs(jetFlavor)==5) flav = 0;
  else if (abs(jetFlavor)==4) flav = 1;
  else flav = 2;

  TH1* h = m_hEff[flav];
  assert(h!=nullptr);

  float aEta = std::abs(eta);
  int binglobal = h->FindBin(pt, aEta);
  int binx, biny, binz;
  h->GetBinXYZ(binglobal, binx, biny, binz); // converts to x, y bins
  int nx = h->GetNbinsX();
  int ny = h->GetNbinsY();

  // under-overflows
  if (binx < 1) binx = 1;
  if (biny < 1) biny = 1;
  if (binx > nx) binx = nx;
  if (biny > ny) biny = ny;

  float eff = h->GetBinContent(binx, biny);

  // protection against wrongly measured efficiencies (low stat) --> reduce pT bin
  while (eff < 0.00000000001 && binx > 1){
    binx--;
    eff = h->GetBinContent(binx, biny);
  }

  return eff;
}
