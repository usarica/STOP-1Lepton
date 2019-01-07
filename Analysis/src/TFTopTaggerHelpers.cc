#include "TFTopTaggerHelpers.h"
#include <TopTagger/TopTagger/interface/TopTaggerResults.h>
#include <TopTagger/TopTagger/interface/TopTaggerUtilities.h>


namespace TFTopTaggerHelpers{
  TopTagger tagger;
  bool TopTaggerCfg_firsttime = setTFTopTaggerConfigFile();
}

bool TFTopTaggerHelpers::setTFTopTaggerConfigFile(){
  // For now, there is a single config file to set
  tagger.setCfgFile(TString(TFTOPTAGGERPKGDATAPATH+"TopTagger.cfg").Data());
  return true;
}

std::vector<TFTopObject*> TFTopTaggerHelpers::getTopsFromResolvedJets(std::vector<AK4JetObject*> const& jets){
  std::vector<TFTopObject*> res;
  if (!TopTaggerCfg_firsttime) return res;
  if (jets.empty()) return res;
  const auto njets = jets.size();

  // Originally from StopAnalysis/StopBabyMaker/JetTree.cc
  std::vector<TLorentzVector> ak4jets_TLV; ak4jets_TLV.reserve(njets);
  std::vector<float> ak4pfjets_CSV; ak4pfjets_CSV.reserve(njets);
  std::vector<float> ak4pfjets_ptD; ak4pfjets_ptD.reserve(njets);
  std::vector<float> ak4pfjets_axis1; ak4pfjets_axis1.reserve(njets);
  std::vector<float> ak4pfjets_axis2; ak4pfjets_axis2.reserve(njets);

  std::vector<int> ak4pfjets_mult; ak4pfjets_mult.reserve(njets);

  std::vector<float> ak4pfjets_chf; ak4pfjets_chf.reserve(njets);
  std::vector<float> ak4pfjets_cef; ak4pfjets_cef.reserve(njets);
  std::vector<float> ak4pfjets_nef; ak4pfjets_nef.reserve(njets);
  std::vector<float> ak4pfjets_muf; ak4pfjets_muf.reserve(njets);
  std::vector<float> ak4pfjets_hhf; ak4pfjets_hhf.reserve(njets);
  std::vector<float> ak4pfjets_hef; ak4pfjets_hef.reserve(njets);
  std::vector<float> ak4pfjets_nhf; ak4pfjets_nhf.reserve(njets);
  std::vector<float> ak4pfjets_phf; ak4pfjets_phf.reserve(njets);
  std::vector<float> ak4pfjets_elf; ak4pfjets_elf.reserve(njets);
  std::vector<float> ak4pfjets_chm; ak4pfjets_chm.reserve(njets);
  std::vector<float> ak4pfjets_nhm; ak4pfjets_nhm.reserve(njets);
  std::vector<float> ak4pfjets_pm; ak4pfjets_pm.reserve(njets);
  std::vector<float> ak4pfjets_em; ak4pfjets_em.reserve(njets);
  std::vector<float> ak4pfjets_mm; ak4pfjets_mm.reserve(njets);
  std::vector<float> ak4pfjets_deepCSVb; ak4pfjets_deepCSVb.reserve(njets);
  std::vector<float> ak4pfjets_deepCSVc; ak4pfjets_deepCSVc.reserve(njets);
  std::vector<float> ak4pfjets_deepCSVl; ak4pfjets_deepCSVl.reserve(njets);
  std::vector<float> ak4pfjets_deepCSVbb; ak4pfjets_deepCSVbb.reserve(njets);
  //std::vector<float> ak4pfjets_deepCSVcc; ak4pfjets_deepCSVcc.reserve(njets);

  for (auto const* jet:jets){
    AK4JetVariables const& extras = jet->extras;

    ak4jets_TLV.emplace_back(jet->x(), jet->y(), jet->z(), jet->t());
    ak4pfjets_mult.emplace_back(extras.totalMultiplicity);

    ak4pfjets_CSV.emplace_back(extras.pfCombinedInclusiveSecondaryVertexV2BJetTag);
    ak4pfjets_ptD.emplace_back(extras.ptDistribution);
    ak4pfjets_axis1.emplace_back(extras.axis1);
    ak4pfjets_axis2.emplace_back(extras.axis2);

    float jetenergy_inv = 1./(extras.undoJEC*jet->energy());
    ak4pfjets_chf.push_back(extras.chargedHadronE * jetenergy_inv);
    ak4pfjets_cef.push_back(extras.chargedEmE * jetenergy_inv);
    ak4pfjets_nef.push_back(extras.neutralEmE * jetenergy_inv);
    ak4pfjets_muf.push_back(extras.muonE * jetenergy_inv);
    ak4pfjets_hhf.push_back(extras.hfHadronE * jetenergy_inv);
    ak4pfjets_hef.push_back(extras.hfEmE * jetenergy_inv);
    ak4pfjets_nhf.push_back(extras.neutralHadronE * jetenergy_inv);
    ak4pfjets_phf.push_back(extras.photonE * jetenergy_inv);
    ak4pfjets_elf.push_back(extras.electronE * jetenergy_inv);

    ak4pfjets_chm.push_back(extras.chargedHadronMultiplicity);
    ak4pfjets_nhm.push_back(extras.neutralHadronMultiplicity);
    ak4pfjets_pm.push_back(extras.photonMultiplicity);
    ak4pfjets_em.push_back(extras.electronMultiplicity);
    ak4pfjets_mm.push_back(extras.muonMultiplicity);

    ak4pfjets_deepCSVb.emplace_back(extras.deepCSVb);
    ak4pfjets_deepCSVc.emplace_back(extras.deepCSVc);
    ak4pfjets_deepCSVl.emplace_back(extras.deepCSVl);
    ak4pfjets_deepCSVbb.emplace_back(extras.deepCSVbb);
    //ak4pfjets_deepCSVcc.emplace_back(extras.deepCSVcc);
  }

  ttUtility::ConstAK4Inputs<float> AK4Inputs(ak4jets_TLV, ak4pfjets_CSV);
  AK4Inputs.addSupplamentalVector("qgPtD", ak4pfjets_ptD);
  AK4Inputs.addSupplamentalVector("qgAxis1", ak4pfjets_axis1);
  AK4Inputs.addSupplamentalVector("qgAxis2", ak4pfjets_axis2);
  auto ak4jets_mult = std::vector<float>(ak4pfjets_mult.begin(), ak4pfjets_mult.end()); // because everything has to be vector<float>
  AK4Inputs.addSupplamentalVector("qgMult", ak4jets_mult);
  AK4Inputs.addSupplamentalVector("recoJetschargedHadronEnergyFraction", ak4pfjets_chf);
  AK4Inputs.addSupplamentalVector("recoJetschargedEmEnergyFraction", ak4pfjets_cef);
  AK4Inputs.addSupplamentalVector("recoJetsneutralEmEnergyFraction", ak4pfjets_nef);
  AK4Inputs.addSupplamentalVector("recoJetsmuonEnergyFraction", ak4pfjets_muf);
  AK4Inputs.addSupplamentalVector("recoJetsHFHadronEnergyFraction", ak4pfjets_hhf);
  AK4Inputs.addSupplamentalVector("recoJetsHFEMEnergyFraction", ak4pfjets_hef);
  AK4Inputs.addSupplamentalVector("recoJetsneutralEnergyFraction", ak4pfjets_nhf);
  AK4Inputs.addSupplamentalVector("PhotonEnergyFraction", ak4pfjets_phf);
  AK4Inputs.addSupplamentalVector("ElectronEnergyFraction", ak4pfjets_elf);
  AK4Inputs.addSupplamentalVector("ChargedHadronMultiplicity", ak4pfjets_chm);
  AK4Inputs.addSupplamentalVector("NeutralHadronMultiplicity", ak4pfjets_nhm);
  AK4Inputs.addSupplamentalVector("PhotonMultiplicity", ak4pfjets_pm);
  AK4Inputs.addSupplamentalVector("ElectronMultiplicity", ak4pfjets_em);
  AK4Inputs.addSupplamentalVector("MuonMultiplicity", ak4pfjets_mm);
  AK4Inputs.addSupplamentalVector("DeepCSVb", ak4pfjets_deepCSVb);
  AK4Inputs.addSupplamentalVector("DeepCSVc", ak4pfjets_deepCSVc);
  AK4Inputs.addSupplamentalVector("DeepCSVl", ak4pfjets_deepCSVl);
  AK4Inputs.addSupplamentalVector("DeepCSVbb", ak4pfjets_deepCSVbb);
  auto dummy_deepCSVcc = std::vector<float>(ak4jets_TLV.size(), 0); // dealing with the fact that deepCSVcc not present in 94X
  AK4Inputs.addSupplamentalVector("DeepCSVcc", dummy_deepCSVcc);
  //AK4Inputs.addSupplamentalVector("DeepCSVcc", ak4pfjets_deepCSVcc);

  // Construct the constitutents and run the tagger
  std::vector<Constituent> constituents = ttUtility::packageConstituents(AK4Inputs);
  tagger.runTagger(constituents);
  // Retrieve the top tagger results object
  const TopTaggerResults& ttr = tagger.getResults();
  // Get the reconstructed tops
  const std::vector<TopObject*>& tftops = ttr.getTops();
  for (TopObject const* top:tftops){
    TFTopObject* obj = new TFTopObject(0, CMSLorentzVector(top->p().Px(), top->p().Py(), top->p().Pz(), top->p().E()));
    TFTopVariables& extras = obj->extras;

    extras.disc = top->getDiscriminator();
    extras.nSubjets = top->getNConstituents();
    for (Constituent const* subjet : top->getConstituents()) extras.subjets_momentum.push_back(subjet->p());

    res.push_back(obj);
  }

  return res;
}

