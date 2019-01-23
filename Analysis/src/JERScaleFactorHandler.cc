#include <cassert>
#include "JERScaleFactorHandler.h"
#include "MELAStreamHelpers.hh"
#include "TRandom3.h"


using namespace std;
using namespace SampleHelpers;
using namespace JECJERHelpers;
using namespace MELAStreamHelpers;


JERScaleFactorHandler::JERScaleFactorHandler(JECJERHelpers::JECJERType type_) :
  ScaleFactorHandlerBase(),
  type(type_)
{
  setup();
}

JERScaleFactorHandler::~JERScaleFactorHandler(){ this->reset(); }


bool JERScaleFactorHandler::setup(){
  bool res = true;
  this->reset();

  TDirectory* curdir = gDirectory;

  resolution_pt_data = JME::JetResolution(getJERPtFileName(type, false, false).Data());
  resolution_phi_data = JME::JetResolution(getJERPhiFileName(type, false, false).Data());
  resolution_pt_mc = JME::JetResolution(getJERPtFileName(type, true, false).Data());
  resolution_phi_mc = JME::JetResolution(getJERPhiFileName(type, true, false).Data());
  resolution_sf = JME::JetResolutionScaleFactor(getJERSFFileName(type, false).Data());

  curdir->cd();

  return res;
}
void JERScaleFactorHandler::reset(){
  resolution_pt_data = JME::JetResolution();
  resolution_phi_data = JME::JetResolution();
  resolution_pt_mc = JME::JetResolution();
  resolution_phi_mc = JME::JetResolution();
  resolution_sf = JME::JetResolutionScaleFactor();
}

void JERScaleFactorHandler::smear(std::vector<AK4JetObject*>& jets, float rho, bool isMC){
  if (type!=JECJERHelpers::kAK4){
    MELAerr << "JERScaleFactorHandler::smear: The JECJERType is not AK4, but smearing for AK4 jets is called!" << endl;
    return;
  }

  /****************************************************************************/
  /* From CJLST/ZZAnalysis/blob/miniAOD_80X/AnalysisStep/plugins/JetFiller.cc */
  /****************************************************************************/
  JME::JetResolution& resolution_pt = (isMC ? resolution_pt_mc : resolution_pt_data);
  for (auto& jet:jets){

    CMSLorentzVector jetMomentum_onlyjec = jet->getFinalMomentum();
    float jpt = jetMomentum_onlyjec.Pt();
    float jeta = jetMomentum_onlyjec.Eta();
    float jphi = jetMomentum_onlyjec.Phi();

    JME::JetParameters res_parameters ={ { JME::Binning::JetPt, jpt }, { JME::Binning::JetEta, jeta }, { JME::Binning::Rho, rho } };
    float res_pt  = resolution_pt.getResolution(res_parameters); // Resolution/pT
    //float res_phi = resolution_phi.getResolution(res_parameters); //not used
    
    float sf=1, sf_up=1, sf_dn=1;
    if (isMC){
      float pt_jer=jpt, pt_jerup=jpt, pt_jerdn=jpt;
      auto* genjet = jet->associatedGenJet;
      JME::JetParameters sf_parameters ={ { JME::Binning::JetEta, jeta },{ JME::Binning::Rho, rho } };
      sf    = resolution_sf.getScaleFactor(sf_parameters);
      sf_up = resolution_sf.getScaleFactor(sf_parameters, Variation::UP);
      sf_dn = resolution_sf.getScaleFactor(sf_parameters, Variation::DOWN);

      bool hasMatched = (genjet!=nullptr);
      bool isMatched = hasMatched;
      float gen_pt=0;
      if (hasMatched){
        float deltaR = reco::deltaR(genjet->momentum, jetMomentum_onlyjec);
        gen_pt = genjet->pt();
        float diff_pt = fabs(jpt - gen_pt);
        isMatched = (deltaR < AK4JetObject::ConeRadiusConstant/2. && diff_pt < 3.f*res_pt*jpt);
      }
      if (isMatched){
        //- apply scaling
        pt_jer   = max(0.f, gen_pt + sf   *(jpt-gen_pt));
        pt_jerup = max(0.f, gen_pt + sf_up*(jpt-gen_pt));
        pt_jerdn = max(0.f, gen_pt + sf_dn*(jpt-gen_pt));
      }
      else{
        //- apply smearing
        TRandom3 rand;
        rand.SetSeed(abs(static_cast<int>(sin(jphi)*100000)));
        float smear = rand.Gaus(0, 1.);
        float sigma   = sqrt(sf   *sf   -1.) * res_pt*jpt;
        float sigmaup = sqrt(sf_up*sf_up-1.) * res_pt*jpt;
        float sigmadn = sqrt(sf_dn*sf_dn-1.) * res_pt*jpt;
        pt_jer   = max(0.f, smear*sigma   + jpt);
        pt_jerup = max(0.f, smear*sigmaup + jpt);
        pt_jerdn = max(0.f, smear*sigmadn + jpt);
      }
      jet->extras.JER = (pt_jer/jpt);
      jet->extras.JERup = (pt_jerup/jpt);
      jet->extras.JERdn = (pt_jerdn/jpt);
    }
    jet->extras.estimatedPtResolution = res_pt*jpt*sf;
  }
}
void JERScaleFactorHandler::smear(std::vector<AK8JetObject*>& jets, float rho, bool isMC){
  if (type!=JECJERHelpers::kAK8){
    MELAerr << "JERScaleFactorHandler::smear: The JECJERType is not AK8, but smearing for AK8 jets is called!" << endl;
    return;
  }

  /****************************************************************************/
  /* From CJLST/ZZAnalysis/blob/miniAOD_80X/AnalysisStep/plugins/JetFiller.cc */
  /****************************************************************************/
  JME::JetResolution& resolution_pt = (isMC ? resolution_pt_mc : resolution_pt_data);
  for (auto& jet:jets){
    CMSLorentzVector jetMomentum_onlyjec = jet->getFinalMomentum();
    float jpt = jetMomentum_onlyjec.Pt();
    float jeta = jetMomentum_onlyjec.Eta();
    float jphi = jetMomentum_onlyjec.Phi();

    JME::JetParameters res_parameters ={ { JME::Binning::JetPt, jpt },{ JME::Binning::JetEta, jeta },{ JME::Binning::Rho, rho } };
    float res_pt  = resolution_pt.getResolution(res_parameters);
    //float res_phi = resolution_phi.getResolution(res_parameters); //not used

    float sf=1, sf_up=1, sf_dn=1;
    if (isMC){
      float pt_jer=jpt, pt_jerup=jpt, pt_jerdn=jpt;
      auto* genjet = jet->associatedGenJet;
      JME::JetParameters sf_parameters ={ { JME::Binning::JetEta, jeta },{ JME::Binning::Rho, rho } };
      sf    = resolution_sf.getScaleFactor(sf_parameters);
      sf_up = resolution_sf.getScaleFactor(sf_parameters, Variation::UP);
      sf_dn = resolution_sf.getScaleFactor(sf_parameters, Variation::DOWN);

      bool hasMatched = (genjet!=nullptr);
      bool isMatched = hasMatched;
      float gen_pt=0;
      if (hasMatched){
        float deltaR = reco::deltaR(genjet->momentum, jetMomentum_onlyjec);
        gen_pt = genjet->pt();
        float diff_pt = fabs(jpt - gen_pt);
        isMatched = (deltaR < AK8JetObject::ConeRadiusConstant/2. && diff_pt < 3.f*res_pt*jpt);
      }
      if (isMatched){
        //- apply scaling
        pt_jer   = max(0.f, gen_pt + sf   *(jpt-gen_pt));
        pt_jerup = max(0.f, gen_pt + sf_up*(jpt-gen_pt));
        pt_jerdn = max(0.f, gen_pt + sf_dn*(jpt-gen_pt));
      }
      else{
        //- apply smearing
        TRandom3 rand;
        rand.SetSeed(abs(static_cast<int>(sin(jphi)*100000)));
        float smear = rand.Gaus(0, 1.);
        float sigma   = sqrt(sf   *sf   -1.) * res_pt*jpt;
        float sigmaup = sqrt(sf_up*sf_up-1.) * res_pt*jpt;
        float sigmadn = sqrt(sf_dn*sf_dn-1.) * res_pt*jpt;
        pt_jer   = max(0.f, smear*sigma   + jpt);
        pt_jerup = max(0.f, smear*sigmaup + jpt);
        pt_jerdn = max(0.f, smear*sigmadn + jpt);
      }
      jet->extras.JER = (pt_jer/jpt);
      jet->extras.JERup = (pt_jerup/jpt);
      jet->extras.JERdn = (pt_jerdn/jpt);
    }
    jet->extras.estimatedPtResolution = res_pt*jpt*sf;
  }
}
