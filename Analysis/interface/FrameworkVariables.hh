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
// End event variables


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
// End event weight variables


// Vertex variables
// vint
// Why are these int anyway?!
#define _vtxs_isValid_ "vtxs_isValid"
#define _vtxs_isFake_ "vtxs_isFake"
// vfloat
#define _vtxs_ndof_ "vtxs_ndof"
// vCMSLorentzVector
#define _vtxs_position_ "vtxs_position"
// End vertex variables


// Pile-up variables
// vint
#define _puinfos_bunchCrossing_ "puInfo_bunchCrossing"
#define _puinfos_nPUVtxs_ "puInfo_nPUvertices"
// vfloat
// Why this is a vfloat is a big mystery (!)
#define _puinfos_nTrueVtxs_ "puInfo_trueNumInteractions"
// End pile-up variables


// PF candidate variables
// vbool
#define _pfcands_trackHighPurity_ "pfcands_trackHighPurity"
// vint
#define _pfcands_charge_ "pfcands_charge"
#define _pfcands_id_ "pfcands_particleId"
// vfloat
#define _pfcands_dxy_ "pfcands_dxy"
#define _pfcands_dz_ "pfcands_dz"
#define _pfcands_dxyError_ "pfcands_dxyError"
#define _pfcands_dzError_ "pfcands_dzError"
// vCMSLorentzVector
#define _pfcands_momentum_ "pfcands_p4"
// End PF candidate variables


// Iso. tracks variables
// vbool
#define _isotracks_isPFCand_ "isotracks_isPFCand"
#define _isotracks_hasLepOverlap_ "isotracks_lepOverlap"
// vint
#define _isotracks_charge_ "isotracks_charge"
#define _isotracks_id_ "isotracks_particleId"
// vfloat
#define _isotracks_pfIso_ch_ "isotracks_pfIso_ch"
#define _isotracks_dz_ "isotracks_dz"
// vCMSLorentzVector
#define _isotracks_momentum_ "isotracks_p4"
// End iso. tracks variables


// Tau variables
// vint
#define _taus_charge_ "taus_pf_charge"
// vCMSLorentzVector
#define _taus_momentum_ "taus_pf_p4"
// vTString
#define _taus_pfidnames_ "taus_pf_IDnames"
// vvfloat -> Why are these floats?!
#define _taus_pfids_ "taus_pf_IDs"
// End tau variables


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


// Photon variables
// float
#define _photons_rho_ "evt_fixgridfastjet_all_rho"
// vfloat
//#define _photons_etaSC_ "photons_etaSC" // FIXME: No such variable is defined in cmstas/CORE yet!
#define _photons_recoChargedHadronIso_ "photons_recoChargedHadronIso"
#define _photons_recoNeutralHadronIso_ "photons_recoNeutralHadronIso"
#define _photons_recoPhotonIso_ "photons_recoPhotonIso"
#define _photons_sigmaIEtaIEta_full5x5_ "photons_full5x5_sigmaIEtaIEta"
#define _photons_hOverE_full5x5_ "photons_full5x5_hOverE"
// vCMSLorentzVector
#define _photons_momentum_ "photons_p4"
// End photon variables


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
#define _muons_IP3D_ "mus_ip3d"
#define _muons_IP3Derr_ "mus_ip3derr"
#define _muons_miniIso_ch_ "mus_miniIso_ch"
#define _muons_miniIso_nh_ "mus_miniIso_nh"
#define _muons_miniIso_em_ "mus_miniIso_em"
// vCMSLorentzVector
#define _muons_momentum_ "mus_p4"
// End muon variables


// ak4 jet variables
// float
#define _ak4jets_rho_ "evt_fixgridfastjet_all_rho"
// vint
#define _ak4jets_npfcands_ "pfjets_npfcands"
#define _ak4jets_parton_flavor_ "pfjets_partonFlavour"
#define _ak4jets_hadron_flavor_ "pfjets_hadronFlavour"
#define _ak4jets_chargedHadronMultiplicity_ "pfjets_chargedHadronMultiplicity"
#define _ak4jets_neutralHadronMultiplicity_ "pfjets_neutralHadronMultiplicity"
#define _ak4jets_photonMultiplicity_ "pfjets_photonMultiplicity"
#define _ak4jets_electronMultiplicity_ "pfjets_electronMultiplicity"
#define _ak4jets_muonMultiplicity_ "pfjets_muonMultiplicity"
#define _ak4jets_chargedMultiplicity_ "pfjets_chargedMultiplicity"
#define _ak4jets_neutralMultiplicity_ "pfjets_neutralMultiplicity"
#define _ak4jets_totalMultiplicity_ "pfjets_totalMultiplicity"
// vfloat
#define _ak4jets_area_ "pfjets_area"
#define _ak4jets_undoJEC_ "pfjets_undoJEC"
#define _ak4jets_chargedHadronE_ "pfjets_chargedHadronE"
#define _ak4jets_chargedEmE_ "pfjets_chargedEmE"
#define _ak4jets_neutralHadronE_ "pfjets_neutralHadronE"
#define _ak4jets_neutralEmE_ "pfjets_neutralEmE"
#define _ak4jets_hfHadronE_ "pfjets_hfHadronE"
#define _ak4jets_hfEmE_ "pfjets_hfEmE"
#define _ak4jets_photonE_ "pfjets_photonE"
#define _ak4jets_electronE_ "pfjets_electronE"
#define _ak4jets_muonE_ "pfjets_muonE"
//
#define _ak4jets_pfCombinedInclusiveSecondaryVertexV2BJetTag_ "pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag"
#define _ak4jets_ptDistribution_ "pfjets_ptDistribution"
#define _ak4jets_axis1_ "pfjets_axis1"
#define _ak4jets_axis2_ "pfjets_axis2"
// vtstring
#define _ak4jets_bDiscriminatorNames_ "pfjets_bDiscriminatorNames"
// vvfloat
#define _ak4jets_bDiscriminators_ "pfjets_bDiscriminators"
// vCMSLorentzVector
#define _ak4jets_momentum_ "pfjets_p4"
// vvCMSLorentzVector
#define _ak4jets_mucands_momentum_ "pfjets_pfcandmup4"
// End ak4 jet variables


// ak8 jet variables
// float
#define _ak8jets_rho_ "evt_fixgridfastjet_all_rho"
// vint
#define _ak8jets_parton_flavor_ "ak8jets_partonFlavour"
// vfloat
#define _ak8jets_area_ "ak8jets_area"
#define _ak8jets_undoJEC_ "ak8jets_undoJEC"
#define _ak8jets_tau1_ "ak8jets_nJettinessTau1"
#define _ak8jets_tau2_ "ak8jets_nJettinessTau2"
#define _ak8jets_tau3_ "ak8jets_nJettinessTau3"
#define _ak8jets_deepdisc_qcd_ "ak8jets_deep_rawdisc_qcd"
#define _ak8jets_deepdisc_top_ "ak8jets_deep_rawdisc_top"
#define _ak8jets_deepdisc_w_ "ak8jets_deep_rawdisc_w"
#define _ak8jets_deepdisc_z_ "ak8jets_deep_rawdisc_z"
#define _ak8jets_deepdisc_zbb_ "ak8jets_deep_rawdisc_zbb"
#define _ak8jets_deepdisc_hbb_ "ak8jets_deep_rawdisc_hbb"
#define _ak8jets_deepdisc_h4q_ "ak8jets_deep_rawdisc_h4q"
// vCMSLorentzVector
#define _ak8jets_momentum_ "ak8jets_p4"
// End ak8 jet variables


// Secondary vertex objects (soft b)
// vint
#define _svs_nTracks_ "svs_nTrks"
// vfloat
#define _svs_IP2D_ "svs_distXYval"
#define _svs_SIP2D_ "svs_distXYsig"
#define _svs_IP3D_ "svs_dist3Dval"
#define _svs_SIP3D_ "svs_dist3Dsig"
#define _svs_anglePV_ "svs_anglePV"
// vCMSLorentzVector
#define _svs_position_ "svs_position"
#define _svs_momentum_ "svs_p4"
// End soft b objects

// Gen. jet variables
// vCMSLorentzVector
#define _genjets_momentum_ "genjets_p4NoMuNoNu"
// End gen. jet variables


// Gen. particle variables
// vbool
#define _genparticles_isPromptFinalState_ "genps_isPromptFinalState"
#define _genparticles_isPromptDecayed_ "genps_isPromptDecayed"
#define _genparticles_isDirectPromptTauDecayProductFinalState_ "genps_isDirectPromptTauDecayProductFinalState"
#define _genparticles_isHardProcess_ "genps_isHardProcess"
#define _genparticles_fromHardProcessFinalState_ "genps_fromHardProcessFinalState"
#define _genparticles_fromHardProcessDecayed_ "genps_fromHardProcessDecayed"
#define _genparticles_isDirectHardProcessTauDecayProductFinalState_ "genps_isDirectHardProcessTauDecayProductFinalState"
#define _genparticles_fromHardProcessBeforeFSR_ "genps_fromHardProcessBeforeFSR"
#define _genparticles_isLastCopy_ "genps_isLastCopy"
#define _genparticles_isLastCopyBeforeFSR_ "genps_isLastCopyBeforeFSR"
// vint
#define _genparticles_status_ "genps_status"
#define _genparticles_id_ "genps_id"
// vCMSLorentzVector
#define _genparticles_p4_ "genps_p4"
// End gen. particle variables


// Gen. info variables
// uint
#define _geninfo_processID_ "genps_signalProcessID"
// float
#define _geninfo_qscale_ "genps_qScale"
#define _geninfo_alphaS_ "genps_alphaQCD"
#define _gen_met_ "gen_met"
#define _gen_metPhi_ "gen_metPhi"
// End gen. info variables


// MET variables
// float
#define _pfmetraw_ "evt_pfmet_raw"
#define _pfmetrawPhi_ "evt_pfmetPhi_raw"
#define _pfmetraw_old_ "evt_old_pfmet_raw"
#define _pfmetrawPhi_old_ "evt_old_pfmetPhi_raw"
#define _pfmetraw_muegclean_ "evt_muegclean_pfmet_raw"
#define _pfmetrawPhi_muegclean_ "evt_muegclean_pfmetPhi_raw"
#define _pfmet_ "evt_pfmet"
#define _pfmetPhi_ "evt_pfmetPhi"


// Vertex variables
// vint
#define _vertex_isFake_ "vtxs_isFake"
// vfloat
#define _vertex_Ndof_ "vtxs_ndof"
// vCMSLorentzVector
#define _vertex_positions_ "vtxs_position"
// End vertex variables


// Event filter variables
// bool
#define _filt_hcalLaser_ "filt_hcalLaser"
#define _filt_ecalBoundaryEnergy_ "filt_ecalBoundaryEnergy"
#define _filt_hcalStrip_ "filt_hcalStrip"
#define _filt_cscBeamHaloTrkMuUnveto_ "filt_cscBeamHaloTrkMuUnveto"
#define _filt_met_filter_ "filt_metfilter"
#define _filt_badMuons_ "filt_badMuons"
#define _filt_BadPFMuon_filter_ "filt_BadPFMuonFilter"
#define _filt_BadChargedCandidate_filter_ "filt_BadChargedCandidateFilter"
#define _filt_duplicateMuons_ "filt_duplicateMuons"
#define _filt_noBadMuons_ "filt_noBadMuons"
#define _filt_trkPOG_filters_ "filt_trkPOGFilters"
#define _filt_chargedHadronTrackResolution_ "filt_chargedHadronTrackResolution"
#define _filt_ecalLaser_ "filt_ecalLaser"
#define _filt_goodVertices_ "filt_goodVertices"
#define _filt_hbheNoiseIso_ "filt_hbheNoiseIso"
#define _filt_hbheNoise_ "filt_hbheNoise"
#define _filt_cscBeamHalo2015_ "filt_cscBeamHalo2015"
#define _filt_ecalTP_ "filt_ecalTP"
#define _filt_cscBeamHalo_ "filt_cscBeamHalo"
#define _filt_globalSuperTightHalo2016_ "filt_globalSuperTightHalo2016"
#define _filt_trackingFailure_ "filt_trackingFailure"
#define _filt_trkPOG_manystripclus53X_ "filt_trkPOG_manystripclus53X"
#define _filt_trkPOG_toomanystripclus53X_ "filt_trkPOG_toomanystripclus53X"
#define _filt_eeBadSc_ "filt_eeBadSc"
#define _filt_ecalBadCalib_filter_ "filt_ecalBadCalibFilter"
#define _filt_ecalBadCalib_filterUpdate_ "filt_ecalBadCalibFilterUpdate"
#define _filt_trkPOG_logErrorTooManyClusters_ "filt_trkPOG_logErrorTooManyClusters"
#define _filt_muonBadTrack_ "filt_muonBadTrack"
#define _filt_globalTightHalo2016_ "filt_globalTightHalo2016"
// End event filter variables


// HLT trigger variables
// TBit
#define _hlt_bits_ "hlt_bits"
// vuint
#define _hlt_prescales_ "hlt_prescales"
#define _hlt_l1prescales_ "hlt_l1prescales"
// vTString
#define _hlt_trigNames_ "hlt_trigNames"
// End HLT trigger variables


#endif
