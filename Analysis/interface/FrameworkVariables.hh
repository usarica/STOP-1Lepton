#ifndef FRAMEWORKVARIABLES_HH
#define FRAMEWORKVARIABLES_HH

// This header defines the branch aliases

// Event variables
// Framework-dependent types
typedef unsigned int RunNumber_t;
typedef unsigned int Lumisection_t;
typedef unsigned long long EventNumber_t;
#define _event_RunNumber_ "evt_run"
#define _event_Lumisection_ "evt_lumiBlock"
#define _event_EventNumber_ "evt_event"

// Event weight variables
// Default weight list
// float
#define _genHEPMCweight_ "genHEPMCweight"
#define _genHEPMCweight_2016_ "genHEPMCweight_2016"
#define _LHEweight_QCDscale_muR1_muF1_ "gen_LHEweight_QCDscale_muR1_muF1"
#define _LHEweight_QCDscale_muR1_muF2_ "gen_LHEweight_QCDscale_muR1_muF2"
#define _LHEweight_QCDscale_muR1_muF0p5_ "gen_LHEweight_QCDscale_muR1_muF0p5"
#define _LHEweight_QCDscale_muR2_muF1_ "gen_LHEweight_QCDscale_muR2_muF1"
#define _LHEweight_QCDscale_muR2_muF2_ "gen_LHEweight_QCDscale_muR2_muF2"
#define _LHEweight_QCDscale_muR2_muF0p5_ "gen_LHEweight_QCDscale_muR2_muF0p5"
#define _LHEweight_QCDscale_muR0p5_muF1_ "gen_LHEweight_QCDscale_muR0p5_muF1"
#define _LHEweight_QCDscale_muR0p5_muF2_ "gen_LHEweight_QCDscale_muR0p5_muF2"
#define _LHEweight_QCDscale_muR0p5_muF0p5_ "gen_LHEweight_QCDscale_muR0p5_muF0p5"
#define _LHEweight_PDFVariation_Up_ "gen_LHEweight_PDFVariation_Up"
#define _LHEweight_PDFVariation_Dn_ "gen_LHEweight_PDFVariation_Dn"
#define _LHEweight_AsMZ_Up_ "gen_LHEweight_AsMZ_Up"
#define _LHEweight_AsMZ_Dn_ "gen_LHEweight_AsMZ_Dn"
#define _LHEweight_PDFVariation_Up_2016_ "gen_LHEweight_PDFVariation_Up_2016"
#define _LHEweight_PDFVariation_Dn_2016_ "gen_LHEweight_PDFVariation_Dn_2016"
#define _LHEweight_AsMZ_Up_2016_ "gen_LHEweight_AsMZ_Up_2016"
#define _LHEweight_AsMZ_Dn_2016_ "gen_LHEweight_AsMZ_Dn_2016"
// Alternative weight list
#define _genHEPMCweight_old_ "genps_weight"
// vfloat
#define _genweights_ "genweights"
// vstring
#define _genweightIDs_ "genweightsID"

// Electron variables
// float
#define _electrons_rho_ "evt_fixgridfastjet_all_rho"
// vbool
#define _electrons_conv_vtx_flag_ "els_conv_vtx_flag"
// vint
#define _electrons_charge_ "els_charge"
#define _electrons_expectedMissingInnerHits_ "els_exp_innerlayers"
// vfloat
#define _electrons_energySC_ "els_eSC"
#define _electrons_etaSC_ "els_etaSC"
#define _electrons_etaSeedSC_ "els_scSeedEta"
#define _electrons_sigmaIEtaIEta_full5x5_ "els_sigmaIEtaIEta_full5x5"
#define _electrons_dEtaIn_ "els_dEtaIn"
#define _electrons_dPhiIn_ "els_dPhiIn"
#define _electrons_hOverE_ "els_hOverE"
#define _electrons_ecalEnergy_ "els_ecalEnergy"
#define _electrons_eOverPIn_ "els_eOverPIn"
#define _electrons_dxyPV_ "els_dxyPV"
#define _electrons_dzPV_ "els_dzPV"
#define _electrons_miniIso_ch_ "els_miniIso_ch"
#define _electrons_miniIso_nh_ "els_miniIso_nh"
#define _electrons_miniIso_em_ "els_miniIso_em"
// vCMSLorentzVector
#define _electrons_momentum_ "els_p4"
// End electron variables


// Muon variables
// float
#define _muons_rho_ "evt_fixgridfastjet_all_rho"
// vuint
#define _muons_POGSelectorBit_ "mus_selectors"
// vint
#define _muons_charge_ "mus_charge"
#define _muons_isPFMuon_ "mus_pid_PFMuon" // Why is this an int???
#define _muons_type_ "mus_type"
#define _muons_validHits_ "mus_validHits"
#define _muons_lostHits_ "mus_lostHits"
#define _muons_expectedMissingInnerHits_ "mus_exp_innerlayers"
#define _muons_expectedMissingOuterHits_ "mus_exp_outerlayers"
#define _muons_GlobalFit_Ndof_ "mus_gfit_ndof"
// vfloat
#define _muons_GlobalFit_Chisq_ "mus_gfit_chi2"
#define _muons_LocalPos_Chisq_ "mus_chi2LocalPosition"
#define _muons_TrkKink_ "mus_trkKink"
#define _muons_SegComp_ "mus_segmCompatibility"
#define _muons_dxyPV_ "mus_dxyPV"
#define _muons_dzPV_ "mus_dzPV"
// vCMSLorentzVector
#define _muons_momentum_ "mus_p4"
// End muon variables


#endif