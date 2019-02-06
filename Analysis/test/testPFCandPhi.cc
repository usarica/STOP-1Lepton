#include <cassert>
#include "common_includes.h"
#include "PDGHelpers.h"
#include "MELAParticle.h"
#include "RooMsgService.h"
#include "RooRelBW2ProngPdf.h"
#include "RooRelBW3ProngPdf.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "TNumericUtil.hh"
#include <HiggsAnalysis/CombinedLimit/interface/AsymPow.h>
#include <limits>
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"


#define DATA_CLASS(t, name, ...) \
  t* name=nullptr; \
  tree->SetBranchStatus(#name, 1); \
  tree->SetBranchAddress(#name, &name);
#define DATA_SIMPLE(t, name, ...) \
  t name=0; \
  tree->SetBranchStatus(#name, 1); \
  tree->SetBranchAddress(#name, &name);
#define BOOK_BRANCH(t, name, outtree, ...) \
  t name; (outtree).Branch(#name, &name);


bool getBestZCandidate(std::vector<MELAParticle>& leptons, MELAParticle& bestZ){
  std::vector<MELAParticle*> leps, aleps;
  for (auto& l:leptons){
    if (l.id>0) leps.push_back(&l);
    else if (l.id<0) aleps.push_back(&l);
  }
  std::pair<MELAParticle*, MELAParticle*> Zpair=std::pair<MELAParticle*, MELAParticle*>(nullptr, nullptr);
  for (auto& lep:leps){
    for (auto& alep:aleps){
      if (PDGHelpers::getCoupledVertex(lep->id, alep->id)!=23) continue;
      if (
        (!Zpair.first || !Zpair.second)
        ||
        fabs((lep->p4 + alep->p4).M()-PDGHelpers::Zmass)<fabs((Zpair.first->p4 + Zpair.second->p4).M()-PDGHelpers::Zmass)
        ){ Zpair.first = lep; Zpair.second = alep; }
    }
  }
  if (Zpair.first && Zpair.second){
    MELAParticle tmpZ(23, Zpair.first->p4 + Zpair.second->p4);
    tmpZ.addDaughter(Zpair.first);
    tmpZ.addDaughter(Zpair.second);
    bestZ.swap(tmpZ);
    return true;
  }
  else return false;
}

float properDPhi(float phi1, float phi2){
  float res = phi1-phi2;
  if (std::abs(res)>=TMath::Pi()) res += (res>0. ? -2. : +2.)*TMath::Pi();
  return res;
}

void matchRecoToGen(std::vector<MELAParticle> const& recolist, std::vector<MELAParticle*> const& genlist, std::vector<MELAParticle*>& reco_matchedGenParts){
  if (recolist.empty()) return;
  reco_matchedGenParts.clear();
  reco_matchedGenParts.assign(recolist.size(), nullptr);
  if (genlist.empty()) return;

  std::vector<unsigned int> remaining_recoparts; if (recolist.size()>0) remaining_recoparts.reserve(recolist.size());
  for (unsigned int ipart=0; ipart<recolist.size(); ipart++) remaining_recoparts.push_back(ipart);
  std::vector<unsigned int> remaining_genparts; remaining_genparts.reserve(genlist.size());
  for (unsigned int ipart=0; ipart<genlist.size(); ipart++) remaining_genparts.push_back(ipart);
  while (!remaining_recoparts.empty() && !remaining_genparts.empty()){
    int chosenRecoPart=-1;
    int chosenGenPart=-1;
    float minDeltaR=-1;
    for (unsigned int const& rpart:remaining_recoparts){
      TLorentzVector const& pReco = recolist.at(rpart).p4;
      for (unsigned int const& gpart:remaining_genparts){
        TLorentzVector const& pGen = genlist.at(gpart)->p4;
        float deltaR = pGen.DeltaR(pReco);
        if (minDeltaR==-1. || deltaR<minDeltaR){
          minDeltaR=deltaR;
          chosenRecoPart=rpart;
          chosenGenPart=gpart;
        }
      }
    }

    if (chosenRecoPart>=0 && chosenGenPart>=0) reco_matchedGenParts.at(chosenRecoPart) = genlist.at(chosenGenPart);
    for (auto it=remaining_recoparts.begin(); it!=remaining_recoparts.end(); it++){ if ((int) *it == chosenRecoPart){ remaining_recoparts.erase(it); break; } }
    for (auto it=remaining_genparts.begin(); it!=remaining_genparts.end(); it++){ if ((int) *it == chosenGenPart){ remaining_genparts.erase(it); break; } }
  }
}

void transformElectricAndMagneticFields(TVector3 const& boostVector, TVector3& vecE, TVector3& vecB){
  double betaSq = boostVector.Mag2();
  if (betaSq>=1.){
    MELAerr << "transformElectricAndMagneticFields: Boost vector has magnitude >= 1" << endl;
    return;
  }

  double gamma = 1./sqrt(1.-betaSq);
  TVector3 newE = (vecE - boostVector.Cross(vecB))*gamma - boostVector.Unit()*(gamma-1.)*vecE.Dot(boostVector.Unit());
  TVector3 newB = (vecB + boostVector.Cross(vecE))*gamma - boostVector.Unit()*(gamma-1.)*vecB.Dot(boostVector.Unit());

  vecE = newE;
  vecB = newB;
  return;
}

void getTransformedAxisCoordinates(TVector3 const& vecE, TVector3 const& vecB, TLorentzVector const& part, TVector3& vecX, TVector3& vecY, TVector3& vecZ){
  TVector3 betaVec = part.BoostVector();
  TVector3 accVec = (vecE + betaVec.Cross(vecB));
  TVector3 dimX = betaVec;
  TVector3 dimZ = dimX.Cross(accVec);
  TVector3 dimY = dimZ.Cross(dimX); // In the direction of accVec perp. to betaVec
  /*
  TVector3 dimZ = accVec.Cross(betaVec);
  TVector3 dimY = dimZ.Cross(betaVec);
  TVector3 dimX = dimY.Cross(dimZ);
  */
  vecX = dimX.Unit();
  vecY = dimY.Unit();
  vecZ = dimZ.Unit();
}

void getTransformedEtaPhi(TLorentzVector const& p4, TVector3 const& vecX, TVector3 const& vecY, TVector3 const& vecZ, float& eta, float& phi){
  TVector3 p3 = p4.Vect();
  TVector3 curlVec = vecZ.Cross(p3);
  float px = curlVec.Dot(vecY);
  float py = -curlVec.Dot(vecX);
  float pt = sqrt(pow(px, 2)+pow(py, 2));
  float pz = p3.Dot(vecZ);
  phi = atan2(py, px);
  eta = asinh(pz/pt);
}

void testPFCandPhi(TString fname, bool useTransformedFrame=true){
  gStyle->SetOptStat(0);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  SampleHelpers::theDataYear=2018;
  SampleHelpers::theDataPeriod="2018";
  SampleHelpers::theDataVersion=SampleHelpers::kCMSSW_10_X;
  float const btagWP = BtagHelpers::getBtagWP(BtagHelpers::kDeepCSV_Medium);
  float const totalLumi=150.f*1000.f;

  bool isWZ = (fname.Contains("WZTo"));
  bool isWZZ = (fname.Contains("WZZ"));
  bool isWWW = (fname.Contains("WWW"));
  bool isDYLLJets = (fname.Contains("DYJetsToLL"));
  bool isTTWLNu = (fname.Contains("TTWJetsToLNu"));
  bool isTTZLLNuNu = (fname.Contains("TTZToLLNuNu"));
  bool isTT2L2Nu = (fname.Contains("TTTo2L2Nu"));
  bool isTT2L2Nu_0j_special = (fname.Contains("TT_DiLept"));
  bool isTT2L2Nu_1j_special = (fname.Contains("TTPlus1Jet_DiLept"));
  bool isZZ4L = (fname.Contains("ZZTo4L"));
  if (!fname.Contains(".root")) fname += ".root";
  const TString strindir="output/test_GenWeights/2018/";
  TString stroutputcore=strindir+"/WZtagging/";
  TString stroutput=stroutputcore;
  TString sampleLabel;
  if (isWZ) sampleLabel="WZTo3LNu_comparison";
  if (isWZZ) sampleLabel="WZZ_comparison";
  if (isWWW) sampleLabel="WWW_comparison";
  if (isDYLLJets) sampleLabel="DYJetsToLL_comparison";
  if (isTTWLNu) sampleLabel="TTWLNu_comparison";
  if (isTTZLLNuNu) sampleLabel="TTZLLNuNu_comparison";
  if (isTT2L2Nu) sampleLabel="TTTo2L2Nu_comparison";
  if (isTT2L2Nu_0j_special) sampleLabel="TTTo2L2Nu_0j_special_comparison";
  if (isTT2L2Nu_1j_special) sampleLabel="TTTo2L2Nu_1j_special_comparison";
  if (isZZ4L) sampleLabel="ZZTo4L_comparison";
  stroutput += sampleLabel + ".root";

  TFile* finput = TFile::Open(strindir+fname, "read");
  TTree* tree = (TTree*) finput->Get("test");
  tree->SetBranchStatus("*", 0);

  // PFCands
  DATA_CLASS(std::vector<float>, pfcands_pt)
  DATA_CLASS(std::vector<float>, pfcands_eta)
  DATA_CLASS(std::vector<float>, pfcands_phi)
  DATA_CLASS(std::vector<float>, pfcands_mass)
  DATA_CLASS(std::vector<bool>, pfcands_trackHighPurity)
  DATA_CLASS(std::vector<int>, pfcands_id)
  DATA_CLASS(std::vector<int>, pfcands_charge)
  DATA_CLASS(std::vector<float>, pfcands_dxy)
  DATA_CLASS(std::vector<float>, pfcands_dz)
  DATA_CLASS(std::vector<float>, pfcands_dxyError)
  DATA_CLASS(std::vector<float>, pfcands_dzError)

  // AK4Jets
  DATA_CLASS(std::vector<float>, ak4jets_pt)
  DATA_CLASS(std::vector<float>, ak4jets_eta)
  DATA_CLASS(std::vector<float>, ak4jets_phi)
  DATA_CLASS(std::vector<float>, ak4jets_mass)
  DATA_CLASS(std::vector<unsigned int>, ak4jets_npfcands)
  DATA_CLASS(std::vector<long long>, ak4jets_selectionBits)
  DATA_CLASS(std::vector<float>, ak4jets_estimatedPtResolution)
  DATA_CLASS(std::vector<float>, ak4jets_deepCSVl)
  DATA_CLASS(std::vector<float>, ak4jets_deepCSVc)
  DATA_CLASS(std::vector<float>, ak4jets_deepCSVb)
  DATA_CLASS(std::vector<float>, ak4jets_deepCSVbb)

  // Gen. particles
  //DATA_CLASS(std::vector<bool>, genparticles_isPromptFinalState)
  //DATA_CLASS(std::vector<bool>, genparticles_isPromptDecayed)
  //DATA_CLASS(std::vector<bool>, genparticles_isDirectPromptTauDecayProductFinalState)
  DATA_CLASS(std::vector<bool>, genparticles_isHardProcess)
  //DATA_CLASS(std::vector<bool>, genparticles_fromHardProcessFinalState)
  //DATA_CLASS(std::vector<bool>, genparticles_fromHardProcessDecayed)
  //DATA_CLASS(std::vector<bool>, genparticles_isDirectHardProcessTauDecayProductFinalState)
  //DATA_CLASS(std::vector<bool>, genparticles_fromHardProcessBeforeFSR)
  //DATA_CLASS(std::vector<bool>, genparticles_isLastCopy)
  //DATA_CLASS(std::vector<bool>, genparticles_isLastCopyBeforeFSR)
  DATA_CLASS(std::vector<int>, genparticles_id)
  DATA_CLASS(std::vector<int>, genparticles_status)
  DATA_CLASS(std::vector<float>, genparticles_px)
  DATA_CLASS(std::vector<float>, genparticles_py)
  DATA_CLASS(std::vector<float>, genparticles_pz)
  DATA_CLASS(std::vector<float>, genparticles_mass)


  TFile* foutput = TFile::Open(stroutput, "recreate");
  TTree tout("Wtagging", "");
  BOOK_BRANCH(int, progenitorVid, tout);
  BOOK_BRANCH(float, jet1_pt, tout);
  BOOK_BRANCH(float, jet1_energy, tout);
  BOOK_BRANCH(float, jet1_eta, tout);
  BOOK_BRANCH(float, jet1_phi, tout);
  BOOK_BRANCH(float, jet1_mass, tout);
  BOOK_BRANCH(float, jet1_prob_udsg, tout);
  BOOK_BRANCH(float, jet1_prob_c, tout);
  BOOK_BRANCH(float, jet1_prob_b, tout);
  BOOK_BRANCH(float, jet1_prob_bb, tout);
  BOOK_BRANCH(float, jet1_genDeltaR, tout);
  BOOK_BRANCH(int, jet1_genId, tout);
  BOOK_BRANCH(unsigned int, jet1_npfcands_allpt, tout);
  BOOK_BRANCH(unsigned int, jet1_npfcands, tout);
  BOOK_BRANCH(unsigned int, jet1_npfcands_pos_dphiplus, tout);
  BOOK_BRANCH(unsigned int, jet1_npfcands_pos_dphiminus, tout);
  BOOK_BRANCH(unsigned int, jet1_npfcands_neg_dphiplus, tout);
  BOOK_BRANCH(unsigned int, jet1_npfcands_neg_dphiminus, tout);
  BOOK_BRANCH(unsigned int, jet1_npfcands_neut_dphiplus, tout);
  BOOK_BRANCH(unsigned int, jet1_npfcands_neut_dphiminus, tout);

  BOOK_BRANCH(float, jet2_pt, tout);
  BOOK_BRANCH(float, jet2_energy, tout);
  BOOK_BRANCH(float, jet2_eta, tout);
  BOOK_BRANCH(float, jet2_phi, tout);
  BOOK_BRANCH(float, jet2_mass, tout);
  BOOK_BRANCH(float, jet2_prob_udsg, tout);
  BOOK_BRANCH(float, jet2_prob_c, tout);
  BOOK_BRANCH(float, jet2_prob_b, tout);
  BOOK_BRANCH(float, jet2_prob_bb, tout);
  BOOK_BRANCH(float, jet2_genDeltaR, tout);
  BOOK_BRANCH(int, jet2_genId, tout);
  BOOK_BRANCH(unsigned int, jet2_npfcands_allpt, tout);
  BOOK_BRANCH(unsigned int, jet2_npfcands, tout);
  BOOK_BRANCH(unsigned int, jet2_npfcands_pos_dphiplus, tout);
  BOOK_BRANCH(unsigned int, jet2_npfcands_pos_dphiminus, tout);
  BOOK_BRANCH(unsigned int, jet2_npfcands_neg_dphiplus, tout);
  BOOK_BRANCH(unsigned int, jet2_npfcands_neg_dphiminus, tout);
  BOOK_BRANCH(unsigned int, jet2_npfcands_neut_dphiplus, tout);
  BOOK_BRANCH(unsigned int, jet2_npfcands_neut_dphiminus, tout);

  BOOK_BRANCH(std::vector<float>, jet1_pfcands_pt, tout);
  BOOK_BRANCH(std::vector<float>, jet1_pfcands_energy, tout);
  BOOK_BRANCH(std::vector<float>, jet1_pfcands_eta, tout);
  BOOK_BRANCH(std::vector<float>, jet1_pfcands_phi, tout);
  BOOK_BRANCH(std::vector<float>, jet1_pfcands_dphi, tout);
  BOOK_BRANCH(std::vector<float>, jet1_pfcands_mass, tout);
  BOOK_BRANCH(std::vector<int>, jet1_pfcands_id, tout);
  BOOK_BRANCH(std::vector<int>, jet1_pfcands_charge, tout);
  BOOK_BRANCH(std::vector<float>, jet2_pfcands_pt, tout);
  BOOK_BRANCH(std::vector<float>, jet2_pfcands_energy, tout);
  BOOK_BRANCH(std::vector<float>, jet2_pfcands_eta, tout);
  BOOK_BRANCH(std::vector<float>, jet2_pfcands_phi, tout);
  BOOK_BRANCH(std::vector<float>, jet2_pfcands_dphi, tout);
  BOOK_BRANCH(std::vector<float>, jet2_pfcands_mass, tout);
  BOOK_BRANCH(std::vector<int>, jet2_pfcands_id, tout);
  BOOK_BRANCH(std::vector<int>, jet2_pfcands_charge, tout);


  float sumCtrWgt[WeightVariables::nWeightTypes]={ 0 };
  int const nEntries = tree->GetEntries();
  for (int ev=0; ev<nEntries; ev++){
    tree->GetEntry(ev);
    HelperFunctions::progressbar(ev, nEntries);

    // Gen particles
    //MELAout << "Constructing gen. parts." << endl;
    size_t const ngenparticles = genparticles_px->size();
    std::vector<MELAParticle> genparticles; genparticles.reserve(ngenparticles);
    std::vector<MELAParticle*> genZlhe, genWlhe, genToplhe, genHardParticles;
    for (size_t ip=0; ip<ngenparticles; ip++){
      const int& id = genparticles_id->at(ip);
      double E = pow(genparticles_px->at(ip), 2)+pow(genparticles_py->at(ip), 2)+pow(genparticles_pz->at(ip), 2);
      E += (genparticles_mass->at(ip)>=0.f ? 1.f : -1.f) * pow(genparticles_mass->at(ip), 2);
      E = sqrt(E);
      TLorentzVector v(genparticles_px->at(ip), genparticles_py->at(ip), genparticles_pz->at(ip), E);
      genparticles.emplace_back(id, v);
      genparticles.back().setGenStatus(genparticles_status->at(ip));
      if (genparticles_isHardProcess->at(ip)){
        genHardParticles.push_back(&(genparticles.back()));
        if (std::abs(id)==6) genToplhe.push_back(&(genparticles.back()));
        else if (std::abs(id)==23) genZlhe.push_back(&(genparticles.back()));
        else if (std::abs(id)==24) genWlhe.push_back(&(genparticles.back()));
      }
    }

    // Find the mother/daughter relationship for the Z, W and tops
    //MELAout << "Constructing gen. Ws" << endl;
    for (auto Wpart:genWlhe){
      int const& idW = Wpart->id;
      TLorentzVector const& pW = Wpart->p4;
      const float WmassSq = pW.M2();
      MELAParticle* pDau1=nullptr;
      MELAParticle* pDau2=nullptr;
      float bestMassDiff=-1;
      for (auto ipart=genHardParticles.begin(); ipart!=genHardParticles.end(); ipart++){
        if (!((*ipart)->genStatus==23 || (*ipart)->genStatus==1 || (*ipart)->genStatus==2)) continue;
        for (auto jpart=ipart+1; jpart!=genHardParticles.end(); jpart++){
          if (!((*jpart)->genStatus==23 || (*jpart)->genStatus==1 || (*jpart)->genStatus==2)) continue;

          if (idW!=PDGHelpers::getCoupledVertex((*ipart)->id, (*jpart)->id)) continue;
          TLorentzVector pTotal = (*ipart)->p4 + (*jpart)->p4;
          float const tmpMassDiff = fabs(pTotal.M2() - WmassSq);
          if (bestMassDiff<0.f || bestMassDiff>tmpMassDiff){
            bestMassDiff = tmpMassDiff;
            pDau1 = *ipart;
            pDau2 = *jpart;
          }
        }
      }
      if (pDau1 && pDau2){
        if (pDau1->id<0 && pDau2->id>0) std::swap(pDau1, pDau2);
        Wpart->addDaughter(pDau1);
        Wpart->addDaughter(pDau2);
        pDau1->addMother(Wpart);
        pDau2->addMother(Wpart);
        /*
        MELAout << "\t- W particle:" << endl;
        MELAout << *Wpart << endl;
        MELAout << "\t\t- Daughters:" << endl;
        MELAout << *pDau1 << endl;
        MELAout << *pDau2 << endl;
        */
      }
    }
    //MELAout << "Constructing gen. Zs" << endl;
    for (auto Zpart:genZlhe){
      int const& idZ = Zpart->id;
      TLorentzVector const& pZ = Zpart->p4;
      const float ZmassSq = pZ.M2();
      MELAParticle* pDau1=nullptr;
      MELAParticle* pDau2=nullptr;
      float bestMassDiff=-1;
      for (auto ipart=genHardParticles.begin(); ipart!=genHardParticles.end(); ipart++){
        if (!((*ipart)->genStatus==23 || (*ipart)->genStatus==1 || (*ipart)->genStatus==2)) continue;
        for (auto jpart=ipart+1; jpart!=genHardParticles.end(); jpart++){
          if (!((*jpart)->genStatus==23 || (*jpart)->genStatus==1 || (*jpart)->genStatus==2)) continue;

          if (idZ!=PDGHelpers::getCoupledVertex((*ipart)->id, (*jpart)->id)) continue;
          TLorentzVector pTotal = (*ipart)->p4 + (*jpart)->p4;
          float const tmpMassDiff = fabs(pTotal.M2() - ZmassSq);
          if (bestMassDiff<0.f || bestMassDiff>tmpMassDiff){
            bestMassDiff = tmpMassDiff;
            pDau1 = *ipart;
            pDau2 = *jpart;
          }
        }
      }
      if (pDau1 && pDau2){
        if (pDau1->id<0 && pDau2->id>0) std::swap(pDau1, pDau2);
        Zpart->addDaughter(pDau1);
        Zpart->addDaughter(pDau2);
        pDau1->addMother(Zpart);
        pDau2->addMother(Zpart);
        /*
        MELAout << "\t- Z particle:" << endl;
        MELAout << *Zpart << endl;
        MELAout << "\t\t- Daughters:" << endl;
        MELAout << *pDau1 << endl;
        MELAout << *pDau2 << endl;
        */
      }
    }
    //MELAout << "Constructing gen. tops" << endl;
    for (auto Toppart:genToplhe){
      int const& idTop = Toppart->id;
      int const bsign = (idTop>0 ? 1 : -1);
      int const Widmatch = 24*bsign;

      TLorentzVector const& pTop = Toppart->p4;
      const float TopmassSq = pTop.M2();
      MELAParticle* pDau1=nullptr;
      MELAParticle* pDau2=nullptr;
      MELAParticle* pDau3=nullptr;
      float bestMassDiff=-1;
      for (auto ipart=genHardParticles.begin(); ipart!=genHardParticles.end(); ipart++){
        if (!((*ipart)->genStatus==23 || (*ipart)->genStatus==1 || (*ipart)->genStatus==2)) continue;
        for (auto jpart=ipart+1; jpart!=genHardParticles.end(); jpart++){
          if (!((*jpart)->genStatus==23 || (*jpart)->genStatus==1 || (*jpart)->genStatus==2)) continue;
          for (auto kpart=jpart+1; kpart!=genHardParticles.end(); kpart++){
            if (!((*kpart)->genStatus==23 || (*kpart)->genStatus==1 || (*kpart)->genStatus==2)) continue;

            int const id12 = PDGHelpers::getCoupledVertex((*ipart)->id, (*jpart)->id);
            int const id13 = PDGHelpers::getCoupledVertex((*ipart)->id, (*kpart)->id);
            int const id23 = PDGHelpers::getCoupledVertex((*jpart)->id, (*kpart)->id);
            if (!((id12==Widmatch && PDGHelpers::isAQuark((*kpart)->id) && (*kpart)->id*bsign>0) || (id13==Widmatch && PDGHelpers::isAQuark((*jpart)->id) && (*jpart)->id*bsign>0) || (id23==Widmatch && PDGHelpers::isAQuark((*ipart)->id) && (*ipart)->id*bsign>0))) continue;
            TLorentzVector pTotal = (*ipart)->p4 + (*jpart)->p4 + (*kpart)->p4;
            float const tmpMassDiff = fabs(pTotal.M2() - TopmassSq);
            if (bestMassDiff<0.f || bestMassDiff>tmpMassDiff){
              bestMassDiff = tmpMassDiff;
              pDau1 = *ipart;
              pDau2 = *jpart;
              pDau3 = *kpart;
            }
          }
        }
      }
      if (pDau1 && pDau2 && pDau3){
        /*
        MELAout << "\t- Top particle:" << endl;
        MELAout << *Toppart << endl;
        MELAout << "\t\t- Potential daughters:" << endl;
        MELAout << *pDau1 << endl;
        MELAout << *pDau2 << endl;
        MELAout << *pDau3 << endl;
        */
        int const id12 = PDGHelpers::getCoupledVertex(pDau1->id, pDau2->id);
        int const id13 = PDGHelpers::getCoupledVertex(pDau1->id, pDau3->id);
        int const id23 = PDGHelpers::getCoupledVertex(pDau2->id, pDau3->id);
        std::pair<MELAParticle*, MELAParticle*> bestWpair(nullptr, nullptr);
        MELAParticle* bestB=nullptr;
        float bestWpairMassDiff=-1;
        if (id12==Widmatch && PDGHelpers::isAQuark(pDau3->id) && pDau3->id*bsign>0 && (bestWpairMassDiff<0.f || bestWpairMassDiff>fabs((pDau1->p4+pDau2->p4).M()-PDGHelpers::Wmass))){
          bestWpairMassDiff = fabs((pDau1->p4+pDau2->p4).M()-PDGHelpers::Wmass);
          bestB=pDau3;
          bestWpair.first=pDau1;
          bestWpair.second=pDau2;
        }
        if (id13==Widmatch && PDGHelpers::isAQuark(pDau2->id) && pDau2->id*bsign>0 && (bestWpairMassDiff<0.f || bestWpairMassDiff>fabs((pDau1->p4+pDau3->p4).M()-PDGHelpers::Wmass))){
          bestWpairMassDiff = fabs((pDau1->p4+pDau3->p4).M()-PDGHelpers::Wmass);
          bestB=pDau2;
          bestWpair.first=pDau1;
          bestWpair.second=pDau3;
        }
        if (id23==Widmatch && PDGHelpers::isAQuark(pDau1->id) && pDau1->id*bsign>0 && (bestWpairMassDiff<0.f || bestWpairMassDiff>fabs((pDau2->p4+pDau3->p4).M()-PDGHelpers::Wmass))){
          bestWpairMassDiff = fabs((pDau2->p4+pDau3->p4).M()-PDGHelpers::Wmass);
          bestB=pDau1;
          bestWpair.first=pDau2;
          bestWpair.second=pDau3;
        }
        if (bestWpair.first->id<0 && bestWpair.second->id>0) std::swap(bestWpair.first, bestWpair.second);
        Toppart->addDaughter(bestB);
        Toppart->addDaughter(bestWpair.first);
        Toppart->addDaughter(bestWpair.second);
        if (!genWlhe.empty()){
          for (auto& genW:genWlhe){
            if (genW->hasDaughter(bestWpair.first) && genW->hasDaughter(bestWpair.second)){
              genW->addMother(Toppart);
              break;
            }
          }
        }
      }
      if (Toppart->getNDaughters()!=3){
        MELAerr << "LHE top particle does not have daughters! Summary of the hard process particles:" << endl;
        for (auto part:genHardParticles) MELAerr << *part << " | status: " << part->genStatus << endl;
        exit(1);
      }
    }

    //MELAout << "Constructing pfcands" << endl;
    size_t const npfcands = pfcands_pt->size();
    std::vector<MELAParticle> pfcands; pfcands.reserve(npfcands);
    std::vector<int> pfcandcharges; pfcandcharges.reserve(npfcands);
    for (size_t ip=0; ip<npfcands; ip++){
      if (std::abs(pfcands_dz->at(ip))>=0.1) continue;
      if (std::abs(pfcands_dxy->at(ip))>=0.2) continue;
      if (pfcands_charge->at(ip)!=0 && !pfcands_trackHighPurity->at(ip)) continue;
      TLorentzVector v; v.SetPtEtaPhiM(pfcands_pt->at(ip), pfcands_eta->at(ip), pfcands_phi->at(ip), pfcands_mass->at(ip));
      pfcands.emplace_back(pfcands_id->at(ip), v);
      pfcandcharges.push_back(pfcands_charge->at(ip));
    }
    size_t const npfcands_preselected = pfcands.size();


    unsigned int nak4jets_preselected;
    std::vector<unsigned int> ak4jets_preselected_npfcands;
    std::vector<int> ak4jets_preselected_matchedParton_Id;
    std::vector<float> ak4jets_preselected_matchedParton_DeltaR;
    std::vector<float> ak4jets_preselected_pt;
    std::vector<float> ak4jets_preselected_eta;
    std::vector<float> ak4jets_preselected_phi;
    std::vector<float> ak4jets_preselected_mass;
    std::vector<float> ak4jets_preselected_relPtResolution;
    std::vector<float> ak4jets_preselected_deepCSVl;
    std::vector<float> ak4jets_preselected_deepCSVc;
    std::vector<float> ak4jets_preselected_deepCSVb;
    std::vector<float> ak4jets_preselected_deepCSVbb;

    //MELAout << "Constructing ak4jets" << endl;
    size_t const nak4jets = ak4jets_pt->size();
    std::vector<MELAParticle> ak4jets; ak4jets.reserve(nak4jets);
    std::vector<std::pair<MELAParticle, float>> ak4jet_ptres_pairs; ak4jet_ptres_pairs.reserve(nak4jets);
    for (size_t ip=0; ip<nak4jets; ip++){
      if (!HelperFunctions::test_bit(ak4jets_selectionBits->at(ip), AK4JetSelectionHelpers::kTightID) || ak4jets_pt->at(ip)<AK4JetSelectionHelpers::ptThr_skim_preselection || fabs(ak4jets_eta->at(ip))>=4.7) continue;
      TLorentzVector v; v.SetPtEtaPhiM(ak4jets_pt->at(ip), ak4jets_eta->at(ip), ak4jets_phi->at(ip), ak4jets_mass->at(ip));
      ak4jets.emplace_back(0, v);
      ak4jet_ptres_pairs.emplace_back(ak4jets.back(), ak4jets_estimatedPtResolution->at(ip));

      ak4jets_preselected_pt.push_back(v.Pt());
      ak4jets_preselected_eta.push_back(v.Eta());
      ak4jets_preselected_phi.push_back(v.Phi());
      ak4jets_preselected_mass.push_back(v.M());
      ak4jets_preselected_relPtResolution.push_back(ak4jets_estimatedPtResolution->at(ip)/v.Pt());
      ak4jets_preselected_deepCSVl.push_back(ak4jets_deepCSVl->at(ip));
      ak4jets_preselected_deepCSVc.push_back(ak4jets_deepCSVc->at(ip));
      ak4jets_preselected_deepCSVb.push_back(ak4jets_deepCSVb->at(ip));
      ak4jets_preselected_deepCSVbb.push_back(ak4jets_deepCSVbb->at(ip));
      ak4jets_preselected_npfcands.push_back(ak4jets_npfcands->at(ip));
    }
    nak4jets_preselected = ak4jets.size();


    std::vector<MELAParticle*> ak4jets_preselected_progenitorV;
    std::vector<MELAParticle*> ak4jets_preselected_bestGenParticles;
    matchRecoToGen(ak4jets, genHardParticles, ak4jets_preselected_bestGenParticles);
    for (size_t ip=0; ip<nak4jets_preselected; ip++){
      auto const& part=ak4jets_preselected_bestGenParticles.at(ip);
      if (!part){
        ak4jets_preselected_progenitorV.push_back(nullptr);
        ak4jets_preselected_matchedParton_Id.push_back(-9000);
        ak4jets_preselected_matchedParton_DeltaR.push_back(-1);
      }
      else{
        ak4jets_preselected_matchedParton_Id.push_back(part->id);
        ak4jets.at(ip).id = part->id;
        ak4jets_preselected_matchedParton_DeltaR.push_back(part->deltaR(ak4jets.at(ip)));
        if (part->getNMothers()>0) ak4jets_preselected_progenitorV.push_back(part->getMother(0));
        else ak4jets_preselected_progenitorV.push_back(nullptr);
      }
    }

    //MELAout << "Starting to match gen Vs" << endl;
    std::vector<MELAParticle*>* genVlhe[2]={ &genWlhe, &genZlhe };
    for (unsigned int iv=0; iv<2; iv++){
      if (genVlhe[iv]->empty()) continue;
      //MELAout << "\t- gen V case " << iv << endl;
      std::vector<MELAParticle> recoWmatchedcands; recoWmatchedcands.reserve(genVlhe[iv]->size());
      // Info about jets
      std::vector<std::pair<MELAParticle, MELAParticle>> recoWmatchedcanddaus; recoWmatchedcanddaus.reserve(genVlhe[iv]->size());
      std::vector<std::pair<float, float>> recoWmatchedcanddaus_dR; recoWmatchedcanddaus_dR.reserve(genVlhe[iv]->size());
      std::vector<std::pair<unsigned int, unsigned int>> recoWmatchedcanddaus_npfcands; recoWmatchedcanddaus_npfcands.reserve(genVlhe[iv]->size());
      std::vector<std::pair<float, float>> recoWmatchedcanddaus_prob_udsg; recoWmatchedcanddaus_prob_udsg.reserve(genVlhe[iv]->size());
      std::vector<std::pair<float, float>> recoWmatchedcanddaus_prob_c; recoWmatchedcanddaus_prob_c.reserve(genVlhe[iv]->size());
      std::vector<std::pair<float, float>> recoWmatchedcanddaus_prob_b; recoWmatchedcanddaus_prob_b.reserve(genVlhe[iv]->size());
      std::vector<std::pair<float, float>> recoWmatchedcanddaus_prob_bb; recoWmatchedcanddaus_prob_bb.reserve(genVlhe[iv]->size());

      for (size_t iw=0; iw<genVlhe[iv]->size(); iw++){
        MELAParticle* p1=nullptr; MELAParticle* p2=nullptr;
        float* p1dr=nullptr; float* p2dr=nullptr;
        unsigned int* p1npfcands = nullptr; unsigned int* p2npfcands = nullptr;
        float* p1probudsg=nullptr; float* p2probudsg=nullptr;
        float* p1probc=nullptr; float* p2probc=nullptr;
        float* p1probb=nullptr; float* p2probb=nullptr;
        float* p1probbb=nullptr; float* p2probbb=nullptr;

        for (size_t ip=0; ip<nak4jets_preselected; ip++){
          if (ak4jets_preselected_progenitorV.at(ip)==genVlhe[iv]->at(iw)){
            if (!p1){
              p1=&(ak4jets.at(ip));
              p1dr=&(ak4jets_preselected_matchedParton_DeltaR.at(ip));
              p1npfcands=&(ak4jets_preselected_npfcands.at(ip));
              p1probudsg=&(ak4jets_preselected_deepCSVl.at(ip));
              p1probc=&(ak4jets_preselected_deepCSVc.at(ip));
              p1probb=&(ak4jets_preselected_deepCSVb.at(ip));
              p1probbb=&(ak4jets_preselected_deepCSVbb.at(ip));
            }
            else if (!p2){
              p2=&(ak4jets.at(ip));
              p2dr=&(ak4jets_preselected_matchedParton_DeltaR.at(ip));
              p2npfcands=&(ak4jets_preselected_npfcands.at(ip));
              p2probudsg=&(ak4jets_preselected_deepCSVl.at(ip));
              p2probc=&(ak4jets_preselected_deepCSVc.at(ip));
              p2probb=&(ak4jets_preselected_deepCSVb.at(ip));
              p2probbb=&(ak4jets_preselected_deepCSVbb.at(ip));
            }
            else{
              MELAerr << "More than 3 jets matched a W boson!" << endl;
              exit(1);
            }
          }
        }
        if (p1 && p2){
          recoWmatchedcands.emplace_back(genVlhe[iv]->at(iw)->id, (p1->p4+p2->p4));
          recoWmatchedcanddaus.emplace_back(*p1, *p2); // Make a separate copy
          recoWmatchedcanddaus_dR.emplace_back(*p1dr, *p2dr);
          recoWmatchedcanddaus_npfcands.emplace_back(*p1npfcands, *p2npfcands);
          recoWmatchedcanddaus_prob_udsg.emplace_back(*p1probudsg, *p2probudsg);
          recoWmatchedcanddaus_prob_c.emplace_back(*p1probc, *p2probc);
          recoWmatchedcanddaus_prob_b.emplace_back(*p1probb, *p2probb);
          recoWmatchedcanddaus_prob_bb.emplace_back(*p1probbb, *p2probbb);
          recoWmatchedcands.back().addDaughter(&(recoWmatchedcanddaus.back().first));
          recoWmatchedcands.back().addDaughter(&(recoWmatchedcanddaus.back().second));
        }
      }

      for (unsigned int iw=0; iw<recoWmatchedcands.size(); iw++){
        MELAParticle& recoW = recoWmatchedcands.at(iw);
        MELAParticle* p1 = recoWmatchedcands.at(iw).getDaughter(0);
        MELAParticle* p2 = recoWmatchedcands.at(iw).getDaughter(1);

        jet1_genDeltaR = recoWmatchedcanddaus_dR.at(iw).first;
        jet2_genDeltaR = recoWmatchedcanddaus_dR.at(iw).second;
        jet1_prob_udsg = recoWmatchedcanddaus_prob_udsg.at(iw).first;
        jet2_prob_udsg = recoWmatchedcanddaus_prob_udsg.at(iw).second;
        jet1_prob_c = recoWmatchedcanddaus_prob_c.at(iw).first;
        jet2_prob_c = recoWmatchedcanddaus_prob_c.at(iw).second;
        jet1_prob_b = recoWmatchedcanddaus_prob_b.at(iw).first;
        jet2_prob_b = recoWmatchedcanddaus_prob_b.at(iw).second;
        jet1_prob_bb = recoWmatchedcanddaus_prob_bb.at(iw).first;
        jet2_prob_bb = recoWmatchedcanddaus_prob_bb.at(iw).second;

        jet1_pfcands_pt.clear();
        jet1_pfcands_energy.clear();
        jet1_pfcands_eta.clear();
        jet1_pfcands_phi.clear();
        jet1_pfcands_dphi.clear();
        jet1_pfcands_mass.clear();
        jet1_pfcands_id.clear();
        jet1_pfcands_charge.clear();

        jet2_pfcands_pt.clear();
        jet2_pfcands_energy.clear();
        jet2_pfcands_eta.clear();
        jet2_pfcands_phi.clear();
        jet2_pfcands_dphi.clear();
        jet2_pfcands_mass.clear();
        jet2_pfcands_id.clear();
        jet2_pfcands_charge.clear();

        jet1_npfcands_allpt = recoWmatchedcanddaus_npfcands.at(iw).first;
        jet1_npfcands = 0;
        jet1_npfcands_pos_dphiplus = 0;
        jet1_npfcands_pos_dphiminus = 0;
        jet1_npfcands_neg_dphiplus = 0;
        jet1_npfcands_neg_dphiminus = 0;
        jet1_npfcands_neut_dphiplus = 0;
        jet1_npfcands_neut_dphiminus = 0;

        jet2_npfcands_allpt = recoWmatchedcanddaus_npfcands.at(iw).second;
        jet2_npfcands = 0;
        jet2_npfcands_pos_dphiplus = 0;
        jet2_npfcands_pos_dphiminus = 0;
        jet2_npfcands_neg_dphiplus = 0;
        jet2_npfcands_neg_dphiminus = 0;
        jet2_npfcands_neut_dphiplus = 0;
        jet2_npfcands_neut_dphiminus = 0;

        std::vector<MELAParticle> pfcands_j1, pfcands_j2;
        std::vector<int> pfcandcharges_j1, pfcandcharges_j2;
        pfcands_j1.reserve(npfcands_preselected);
        pfcands_j2.reserve(npfcands_preselected);
        for (size_t ipf=0; ipf<npfcands_preselected; ipf++){
          MELAParticle const& pfcand = pfcands.at(ipf);
          if (pfcand.deltaR(*p1)<0.4 && (pfcand.deltaR(*p2)>=0.4 || pfcand.deltaR(*p2)>pfcand.deltaR(*p1))){
            if (pfcand.t()>p1->t()) continue;
            // Make a copy of the PFCand
            pfcands_j1.emplace_back(pfcand); p1->addDaughter(&(pfcands_j1.back()));
            pfcandcharges_j1.push_back(pfcandcharges.at(ipf));
          }
          else if (pfcand.deltaR(*p2)<0.4){
            if (pfcand.t()>p2->t()) continue;
            // Make a copy of the PFCand
            pfcands_j2.emplace_back(pfcand); p2->addDaughter(&(pfcands_j2.back()));
            pfcandcharges_j2.push_back(pfcandcharges.at(ipf));
          }
        }

        TVector3 bv = -recoW.p4.BoostVector();
        recoW.boost(bv, true);

        TVector3 vecX(1, 0, 0);
        TVector3 vecY(0, 1, 0);
        TVector3 vecZ(0, 0, 1);
        if (useTransformedFrame){
          TVector3 vecE(0, 0, 0);
          TVector3 vecB(0, 0, 1);
          transformElectricAndMagneticFields(bv, vecE, vecB);
          getTransformedAxisCoordinates(vecE, vecB, recoW.getDaughter(0)->p4, vecX, vecY, vecZ);
        }

        progenitorVid = recoW.id;

        getTransformedEtaPhi(p1->p4, vecX, vecY, vecZ, jet1_eta, jet1_phi);
        jet1_pt = p1->pt();
        jet1_energy = p1->t();
        jet1_mass = p1->m();
        jet1_genId = p1->id;

        getTransformedEtaPhi(p2->p4, vecX, vecY, vecZ, jet2_eta, jet2_phi);
        jet2_pt = p2->pt();
        jet2_energy = p2->t();
        jet2_mass = p2->m();
        jet2_genId = p2->id;

        jet1_npfcands = pfcands_j1.size();
        for (unsigned int ipf=0; ipf<jet1_npfcands; ipf++){
          MELAParticle const& pfc = pfcands_j1.at(ipf);
          float pfc_eta, pfc_phi;
          getTransformedEtaPhi(pfc.p4, vecX, vecY, vecZ, pfc_eta, pfc_phi);
          jet1_pfcands_pt.push_back(pfc.pt());
          jet1_pfcands_energy.push_back(pfc.t());
          jet1_pfcands_eta.push_back(pfc_eta);
          jet1_pfcands_phi.push_back(pfc_phi);
          jet1_pfcands_dphi.push_back(properDPhi(jet1_pfcands_phi.back(), jet1_phi));
          jet1_pfcands_mass.push_back(pfc.m());
          jet1_pfcands_id.push_back(pfc.id);
          jet1_pfcands_charge.push_back(pfcandcharges_j1.at(ipf));

          if (jet1_pfcands_dphi.back()>=0.f && jet1_pfcands_charge.back()>0) jet1_npfcands_pos_dphiplus++;
          else if (jet1_pfcands_dphi.back()>=0.f && jet1_pfcands_charge.back()==0) jet1_npfcands_neut_dphiplus++;
          else if (jet1_pfcands_dphi.back()>=0.f && jet1_pfcands_charge.back()<0) jet1_npfcands_neg_dphiplus++;
          else if (jet1_pfcands_dphi.back()<0.f && jet1_pfcands_charge.back()>0) jet1_npfcands_pos_dphiminus++;
          else if (jet1_pfcands_dphi.back()<0.f && jet1_pfcands_charge.back()==0) jet1_npfcands_neut_dphiminus++;
          else if (jet1_pfcands_dphi.back()<0.f && jet1_pfcands_charge.back()<0) jet1_npfcands_neg_dphiminus++;
        }
        jet2_npfcands = pfcands_j2.size();
        for (unsigned int ipf=0; ipf<jet2_npfcands; ipf++){
          MELAParticle const& pfc = pfcands_j2.at(ipf);
          float pfc_eta, pfc_phi;
          getTransformedEtaPhi(pfc.p4, vecX, vecY, vecZ, pfc_eta, pfc_phi);
          jet2_pfcands_pt.push_back(pfc.pt());
          jet2_pfcands_energy.push_back(pfc.t());
          jet2_pfcands_eta.push_back(pfc_eta);
          jet2_pfcands_phi.push_back(pfc_phi);
          jet2_pfcands_dphi.push_back(properDPhi(jet2_pfcands_phi.back(), jet2_phi));
          jet2_pfcands_mass.push_back(pfc.m());
          jet2_pfcands_id.push_back(pfc.id);
          jet2_pfcands_charge.push_back(pfcandcharges_j2.at(ipf));

          if (jet2_pfcands_dphi.back()>=0.f && jet2_pfcands_charge.back()>0) jet2_npfcands_pos_dphiplus++;
          else if (jet2_pfcands_dphi.back()>=0.f && jet2_pfcands_charge.back()==0) jet2_npfcands_neut_dphiplus++;
          else if (jet2_pfcands_dphi.back()>=0.f && jet2_pfcands_charge.back()<0) jet2_npfcands_neg_dphiplus++;
          else if (jet2_pfcands_dphi.back()<0.f && jet2_pfcands_charge.back()>0) jet2_npfcands_pos_dphiminus++;
          else if (jet2_pfcands_dphi.back()<0.f && jet2_pfcands_charge.back()==0) jet2_npfcands_neut_dphiminus++;
          else if (jet2_pfcands_dphi.back()<0.f && jet2_pfcands_charge.back()<0) jet2_npfcands_neg_dphiminus++;
        }

        tout.Fill();
      }
      //MELAout << "\t- End case" << endl;
    }
    //MELAout << "End event" << endl;
  }
  foutput->WriteTObject(&tout); // Write the output tree
  foutput->Close();
  finput->Close();
}
