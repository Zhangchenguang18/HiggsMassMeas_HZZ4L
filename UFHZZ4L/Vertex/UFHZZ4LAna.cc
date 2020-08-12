// -*- C++ -*-
//
// Package:    UFHZZ4LAna
// Class:      UFHZZ4LAna
// 
//

// system include files
#include <memory>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iomanip>
#include <set>

#define PI 3.14159

// user include files 
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TSpline.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TClonesArray.h"
#include "TCanvas.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

//HTXS
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
//#include "SimDataFormats/HZZFiducial/interface/HZZFiducialVolume.h"

// PAT
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// Reco
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

// KD's
#include "ZZMatrixElement/MELA/interface/Mela.h"

//Helper
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LHelper.h"
//Muons
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LMuonAna.h"
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LMuonTree.h"
//Electrons
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LElectronTree.h"
//Photons
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LPhotonTree.h"
//Jets
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LJetTree.h"
//Final Leps
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LFinalLepTree.h"
//Sip
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LSipAna.h"
//PU
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LPileUp.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//GEN
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LGENAna.h"
//VBF Jets
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LJets.h"

// Jet energy correction
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include <vector>

// Kinematic Fit
#include "KinZfitter/KinZfitter/interface/KinZfitter.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"

// EWK corrections
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/EwkCorrections.h"

// JEC related
#include "PhysicsTools/PatAlgos/plugins/PATJetUpdater.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
 
//JER related
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

//BTag Calibration

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

//Muon MVA
//#include "MuonMVAReader/Reader/interface/MuonGBRForestReader.hpp"

//Track Re-Fit
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//Correct pt
//#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"
#include "UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/src/RoccoR.cc"

//
// class declaration
//
using namespace EwkCorrections;

class UFHZZ4LAna : public edm::EDAnalyzer {
public:
    explicit UFHZZ4LAna(const edm::ParameterSet&);
    ~UFHZZ4LAna();
  
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    static bool sortByPt( const reco::GenParticle &p1, const reco::GenParticle &p2 ){ return (p1.pt() > p2.pt()); };
  
private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
  
    virtual void beginRun(edm::Run const&, const edm::EventSetup& iSetup);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const& lumiSeg,edm::EventSetup const& eSetup);
  
    void findHiggsCandidate(std::vector< pat::Muon > &selectedMuons, std::vector< pat::Electron > &selectedElectrons, const edm::Event& iEvent);
    void findZ1LCandidate(const edm::Event& iEvent);

    //Helper Class
    HZZ4LHelper helper;
    //GEN
    HZZ4LGENAna genAna;
    //VBF
    HZZ4LJets jetHelper;
    //PU Reweighting
    edm::LumiReWeighting *lumiWeight;
    HZZ4LPileUp pileUp;
    //JES Uncertainties
    std::unique_ptr<JetCorrectionUncertainty> jecunc;
    // kfactors
    TSpline3 *kFactor_ggzz;
    std::vector<std::vector<float> > tableEwk;
    // data/MC scale factors
    TH2F *hElecScaleFac;
    TH2F *hElecScaleFac_Cracks;
    TH2F *hElecScaleFacGsf;
    TH2F *hElecScaleFacGsfLowET;
    TH2F *hMuScaleFac;
    TH2F *hMuScaleFacUnc;
    TH1D *h_pileup;
    TH1D *h_pileupUp;
    TH1D *h_pileupDn;
    std::vector<TH1F*> h_medians;
    TH2F *hbTagEffi;
    TH2F *hcTagEffi;
    TH2F *hudsgTagEffi;

    BTagCalibrationReader* reader;

    //Saved Events Trees
    TTree *passedEventsTree_All;

    void bookPassedEventTree(TString treeName, TTree *tree);
    void setTreeVariables( const edm::Event&, const edm::EventSetup&, 
                           std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons, 
                           std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons, 
                           std::vector<pat::Jet> goodJets, std::vector<float> goodJetQGTagger, 
                           std::vector<float> goodJetaxis2, std::vector<float> goodJetptD, std::vector<int> goodJetmult, 
                           std::vector<pat::Jet> selectedMergedJets);
    void setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles,
                         edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles,
                         edm::Handle<edm::View<reco::GenJet> > genJets);
    bool mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4, bool makeCuts);

    // -------------------------
    // RECO level information
    // -------------------------

    // Event Variables
    ULong64_t Run, Event, LumiSect;
    int nVtx, nInt;
    int finalState;
    std::string triggersPassed;
    bool passedTrig, passedFullSelection, passedZ4lSelection, passedQCDcut;
    bool passedZ1LSelection, passedZ4lZ1LSelection, passedZ4lZXCRSelection, passedZXCRSelection;
    int nZXCRFailedLeptons;   


    // Event Weights
    float genWeight, pileupWeight, pileupWeightUp, pileupWeightDn, dataMCWeight, eventWeight, prefiringWeight;
    float k_ggZZ, k_qqZZ_qcd_dPhi, k_qqZZ_qcd_M, k_qqZZ_qcd_Pt, k_qqZZ_ewk;
    // pdf weights                                                                   
    vector<float> qcdWeights;
    vector<float> nnloWeights;
    vector<float> pdfWeights;
    int posNNPDF;
    float pdfRMSup, pdfRMSdown, pdfENVup, pdfENVdown;
    // lepton variables
    vector<double> lep_pt; vector<double> lep_pterr; vector<double> lep_pterrold; 
    vector<double> lep_p; vector<double> lep_ecalEnergy; vector<int> lep_isEB; vector<int> lep_isEE;
    vector<double> lep_eta; vector<double> lep_phi; vector<double> lep_mass;
    vector<double> lepFSR_pt; vector<double> lepFSR_eta; vector<double> lepFSR_phi; vector<double> lepFSR_mass;
    int lep_Hindex[4];//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z2 sub
    float pTL1, pTL2, pTL3, pTL4;
    float etaL1, etaL2, etaL3, etaL4;
    int idL1, idL2, idL3, idL4;
    float mL1, mL2, mL3, mL4;
    float pTErrL1, pTErrL2, pTErrL3, pTErrL4;
    float phiL1, phiL2, phiL3, phiL4;
    float pTL1FSR, pTL2FSR, pTL3FSR, pTL4FSR;
    vector<float> lep_dataMC; vector<float> lep_dataMCErr;
    vector<int> lep_genindex; //position of lepton in GENlep_p4 (if gen matched, -1 if not gen matched)
    vector<int> lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId; // gen matching even if not in GENlep_p4
    vector<int> lep_id;
    vector<float> lep_mva; vector<int> lep_ecalDriven; 
    vector<int> lep_tightId; vector<int> lep_tightIdSUS; vector<int> lep_tightIdHiPt; //vector<int> lep_tightId_old;
    vector<float> lep_Sip; vector<float> lep_IP; vector<float> lep_isoNH; vector<float> lep_isoCH; vector<float> lep_isoPhot;
    vector<float> lep_isoPU; vector<float> lep_isoPUcorr; 
    vector<float> lep_RelIso; vector<float> lep_RelIsoNoFSR; vector<float> lep_MiniIso; 
    vector<float> lep_ptRatio; vector<float> lep_ptRel;
    vector<int> lep_missingHits;
    vector<string> lep_filtersMatched; // for each lepton, all filters it is matched to
    int nisoleptons;
    double muRho, elRho, rhoSUS;
	
	// electron scale unc study
	vector<double> eleUncorr_pt; vector<double> eleUncorr_mass;
	
	// no roccor muon
	vector<double> muonUncorr_pt; vector<double> muonUncorr_pterr; vector<double> muon_scale_sf;
	
	// vertex constraint study
	vector<double> track_vtx_pt; vector<double> track_vtx_pterr; vector<double> track_vtx_eta; vector<double> track_vtx_phi;
	vector<double> track_Hvtx_pt; vector<double> track_Hvtx_pterr; vector<double> track_Hvtx_eta; vector<double> track_Hvtx_phi;
//	vector<double> track_vtx_ptcorr; vector<double> track_vtx_pterrcorr;
//	vector<double> track_vtx_FSRpt; vector<double> track_vtx_FSRpterr; vector<double> track_vtx_FSReta; vector<double> track_vtx_FSRphi; vector<double> track_vtx_FSRmass;
//	vector<double> track_vtx_FSRptcorr; vector<double> track_vtx_FSRpterrcorr;
	
	vector<double> track_vtx_dxy;
	vector<double> track_vtx_dxycorr;	
	vector<double> track_Hvtx_dxycorr;
	
	double track_vtx_xcorr, track_vtx_ycorr, track_vtx_zcorr;
	double track_Hvtx_xcorr, track_Hvtx_ycorr, track_Hvtx_zcorr;
	double pv_x, pv_y, pv_z;
	
	
	vector<double> track_bs_pt; vector<double> track_bs_pterr; vector<double> track_bs_eta; vector<double> track_bs_phi;
	vector<double> track_Hbs_pt; vector<double> track_Hbs_pterr; vector<double> track_Hbs_eta; vector<double> track_Hbs_phi;
//	vector<double> track_bs_ptcorr; vector<double> track_bs_pterrcorr;
//	vector<double> track_bs_FSRpt; vector<double> track_bs_FSRpterr; vector<double> track_bs_FSReta; vector<double> track_bs_FSRphi; vector<double> track_bs_FSRmass;
//	vector<double> track_bs_FSRptcorr; vector<double> track_bs_FSRpterrcorr;
	
	vector<double> track_bs_dxy;
	vector<double> track_bs_dxycorr;
	vector<double> track_Hbs_dxycorr;
	
	double track_bs_xcorr, track_bs_ycorr, track_bs_zcorr;
	double track_Hbs_xcorr, track_Hbs_ycorr, track_Hbs_zcorr;
	double bs_x, bs_y, bs_z;

	vector<double> mu_pt_roc; vector<double> mu_eta_roc; vector<double> mu_phi_roc; vector<int> mu_q_roc; vector<int> mu_tracktype_roc; vector<int> mu_nl_roc; vector<int> mu_genparticle_roc; vector<double> mu_gen_pt_roc;
	vector<reco::TransientTrack> tracks_tmp;
	vector<reco::TransientTrack> H4Lrefittedtracks;
	
	vector<unsigned int> passLepindex;//for Z recovery
    
	// tau variables
    vector<int> tau_id;
    vector<double> tau_pt, tau_eta, tau_phi, tau_mass;

    // photon variables
    vector<double> pho_pt, pho_eta, pho_phi, photonCutBasedIDLoose;

    // Higgs candidate variables
    vector<double> H_pt; vector<double> H_eta; vector<double> H_phi; vector<double> H_mass;
    vector<double> H_noFSR_pt; vector<double> H_noFSR_eta; vector<double> H_noFSR_phi; vector<double> H_noFSR_mass;
    float mass4l, mass4l_noFSR, mass4e, mass4mu, mass2e2mu, pT4l, eta4l, phi4l, rapidity4l;
    float cosTheta1, cosTheta2, cosThetaStar, Phi, Phi1;
    float mass3l;

    // kin fit
    float mass4lREFIT, massZ1REFIT, massZ2REFIT, mass4lErr, mass4lErrREFIT;

    // Z candidate variables
    vector<double> Z_pt; vector<double> Z_eta; vector<double> Z_phi; vector<double> Z_mass;
    vector<double> Z_noFSR_pt; vector<double> Z_noFSR_eta; vector<double> Z_noFSR_phi; vector<double> Z_noFSR_mass;
    int Z_Hindex[2]; // position of Z1 and Z2 in Z_p4
    float massZ1, massZ1_Z1L, massZ2, pTZ1, pTZ2;

    // MET
    float met; float met_phi;
    float met_jesup, met_phi_jesup, met_jesdn, met_phi_jesdn;
    float met_uncenup, met_phi_uncenup, met_uncendn, met_phi_uncendn;

    // Jets
    vector<int>    jet_iscleanH4l;
    int jet1index, jet2index;
    vector<double> jet_pt; vector<double> jet_eta; vector<double> jet_phi; vector<double> jet_mass; vector<double> jet_pt_raw;
    vector<float>  jet_pumva, jet_csvv2, jet_csvv2_; vector<int> jet_isbtag;
    vector<int>    jet_hadronFlavour, jet_partonFlavour;
    vector<float>  jet_QGTagger, jet_QGTagger_jesup, jet_QGTagger_jesdn; 
    vector<float>  jet_axis2, jet_ptD; vector<int> jet_mult;
    vector<float>  jet_relpterr; vector<float>  jet_phierr;
    vector<float>  jet_bTagEffi;
    vector<float>  jet_cTagEffi;
    vector<float>  jet_udsgTagEffi;
    vector<int>    jet_jesup_iscleanH4l;
    vector<double> jet_jesup_pt; vector<double> jet_jesup_eta; 
    vector<double> jet_jesup_phi; vector<double> jet_jesup_mass;
    vector<int>    jet_jesdn_iscleanH4l;
    vector<double> jet_jesdn_pt; vector<double> jet_jesdn_eta; 
    vector<double> jet_jesdn_phi; vector<double> jet_jesdn_mass;
    vector<int>    jet_jerup_iscleanH4l;
    vector<double> jet_jerup_pt; vector<double> jet_jerup_eta; 
    vector<double> jet_jerup_phi; vector<double> jet_jerup_mass;
    vector<int>    jet_jerdn_iscleanH4l;
    vector<double> jet_jerdn_pt; vector<double> jet_jerdn_eta; 
    vector<double> jet_jerdn_phi; vector<double> jet_jerdn_mass;    
    int njets_pt30_eta4p7; int njets_pt30_eta4p7_jesup; int njets_pt30_eta4p7_jesdn; 
    int njets_pt30_eta4p7_jerup; int njets_pt30_eta4p7_jerdn;
    int njets_pt30_eta2p5; int njets_pt30_eta2p5_jesup; int njets_pt30_eta2p5_jesdn; 
    int njets_pt30_eta2p5_jerup; int njets_pt30_eta2p5_jerdn;
    int nbjets_pt30_eta4p7; int nvjets_pt40_eta2p4;
    float pt_leadingjet_pt30_eta4p7;
    float pt_leadingjet_pt30_eta4p7_jesup; float pt_leadingjet_pt30_eta4p7_jesdn;
    float pt_leadingjet_pt30_eta4p7_jerup; float pt_leadingjet_pt30_eta4p7_jerdn;
    float pt_leadingjet_pt30_eta2p5;
    float pt_leadingjet_pt30_eta2p5_jesup; float pt_leadingjet_pt30_eta2p5_jesdn;
    float pt_leadingjet_pt30_eta2p5_jerup; float pt_leadingjet_pt30_eta2p5_jerdn;
    float absrapidity_leadingjet_pt30_eta4p7;
    float absrapidity_leadingjet_pt30_eta4p7_jesup; float absrapidity_leadingjet_pt30_eta4p7_jesdn;
    float absrapidity_leadingjet_pt30_eta4p7_jerup; float absrapidity_leadingjet_pt30_eta4p7_jerdn;
    float absdeltarapidity_hleadingjet_pt30_eta4p7;
    float absdeltarapidity_hleadingjet_pt30_eta4p7_jesup; float absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn;
    float absdeltarapidity_hleadingjet_pt30_eta4p7_jerup; float absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn;
    float DijetMass, DijetDEta, DijetFisher;

    // merged jets
    vector<int>   mergedjet_iscleanH4l;
    vector<float> mergedjet_pt; vector<float> mergedjet_eta; vector<float> mergedjet_phi; vector<float> mergedjet_mass;
    
    vector<float> mergedjet_tau1; vector<float> mergedjet_tau2;
    vector<float> mergedjet_btag;

    vector<float> mergedjet_L1;
    vector<float> mergedjet_prunedmass; vector<float> mergedjet_softdropmass;
    
    vector<int> mergedjet_nsubjet;
    vector<vector<float> > mergedjet_subjet_pt; vector<vector<float> > mergedjet_subjet_eta;
    vector<vector<float> > mergedjet_subjet_phi; vector<vector<float> > mergedjet_subjet_mass;
    vector<vector<float> > mergedjet_subjet_btag;
    vector<vector<int> > mergedjet_subjet_partonFlavour, mergedjet_subjet_hadronFlavour;

    // FSR Photons
    int nFSRPhotons;
    vector<int> fsrPhotons_lepindex; 
    vector<double> fsrPhotons_pt; vector<double> fsrPhotons_pterr;
    vector<double> fsrPhotons_eta; vector<double> fsrPhotons_phi;
    vector<double> fsrPhotons_mass;
    vector<float> fsrPhotons_dR; vector<float> fsrPhotons_iso;
    vector<float> allfsrPhotons_dR; vector<float> allfsrPhotons_pt; vector<float> allfsrPhotons_iso;

    // Z4l? FIXME
    float theta12, theta13, theta14;  
    float minM3l, Z4lmaxP, minDeltR, m3l_soft;
    float minMass2Lep, maxMass2Lep;
    float thetaPhoton, thetaPhotonZ;

    // Event Category
    int EventCat;

    // -------------------------
    // GEN level information
    // -------------------------

    //Event variables
    int GENfinalState;
    bool passedFiducialSelection;

    // lepton variables
    vector<double> GENlep_pt; vector<double> GENlep_eta; vector<double> GENlep_phi; vector<double> GENlep_mass; 
    vector<int> GENlep_id; vector<int> GENlep_status; 
    vector<int> GENlep_MomId; vector<int> GENlep_MomMomId;
    int GENlep_Hindex[4];//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z3 sub
    vector<float> GENlep_isoCH; vector<float> GENlep_isoNH; vector<float> GENlep_isoPhot; vector<float> GENlep_RelIso; 

    // Higgs candidate variables (calculated using selected gen leptons)
    vector<double> GENH_pt; vector<double> GENH_eta; vector<double> GENH_phi; vector<double> GENH_mass; 
    float GENmass4l, GENmass4e, GENmass4mu, GENmass2e2mu, GENpT4l, GENeta4l, GENrapidity4l;
    float GENMH; //mass directly from gen particle with id==25
    float GENcosTheta1, GENcosTheta2, GENcosThetaStar, GENPhi, GENPhi1;

    // Z candidate variables
    vector<double> GENZ_pt; vector<double> GENZ_eta; vector<double> GENZ_phi; vector<double> GENZ_mass; 
    vector<int> GENZ_DaughtersId; vector<int> GENZ_MomId;
    float  GENmassZ1, GENmassZ2, GENpTZ1, GENpTZ2, GENdPhiZZ, GENmassZZ, GENpTZZ;

    // Higgs variables directly from GEN particle
    float GENHmass;

    // Jets
    vector<double> GENjet_pt; vector<double> GENjet_eta; vector<double> GENjet_phi; vector<double> GENjet_mass; 
    int GENnjets_pt30_eta4p7; float GENpt_leadingjet_pt30_eta4p7; 
    int GENnjets_pt30_eta2p5; float GENpt_leadingjet_pt30_eta2p5; 
    float GENabsrapidity_leadingjet_pt30_eta4p7; float GENabsdeltarapidity_hleadingjet_pt30_eta4p7;
    int lheNb, lheNj, nGenStatus2bHad;

    // STXS info
    int stage0cat;
    int stage1cat;
    int stage1p1cat;
    int stage1p2cat;
    // Fiducial Rivet
    int passedFiducialRivet;
    float GENpT4lRivet;
    int GENnjets_pt30_eta4p7Rivet;
    float GENpt_leadingjet_pt30_eta4p7Rivet;

    //KinZfitter
    KinZfitter *kinZfitter;

    // MEM
    Mela* mela;

    float me_0plus_JHU, me_qqZZ_MCFM, p0plus_m4l, bkg_m4l;
    float D_bkg_kin, D_bkg, D_g4, D_g1g4;

    float p0minus_VAJHU, Dgg10_VAMCFM, pg1g4_VAJHU;

    // old but working
    float phjj_VAJHU, pvbf_VAJHU;
    float pwh_hadronic_VAJHU, pzh_hadronic_VAJHU;
    float pAux_vbf_VAJHU, phj_VAJHU;

    float p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal;
    float p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal;
    float p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal;
    float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
    float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
    float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
    float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
    float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
    float p_JQCD_SIG_ghv1_1_JHUGen_JECNominal;
    float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;

    float p_JJVBF_BKG_MCFM_JECNominal;
    float p_HadZH_BKG_MCFM_JECNominal;
    float p_HadWH_BKG_MCFM_JECNominal;
    float p_JJQCD_BKG_MCFM_JECNominal;

    float p_HadZH_mavjj_JECNominal;
    float p_HadZH_mavjj_true_JECNominal;
    float p_HadWH_mavjj_JECNominal;
    float p_HadWH_mavjj_true_JECNominal;

    float pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal;
    float pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal;
    float pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal;
    float pConst_JJVBF_BKG_MCFM_JECNominal;
    float pConst_HadZH_BKG_MCFM_JECNominal;
    float pConst_HadWH_BKG_MCFM_JECNominal;
    float pConst_JJQCD_BKG_MCFM_JECNominal;
    float pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
    float pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
    float pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
    float pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;

    float D_VBF, D_VBF1j, D_HadWH, D_HadZH;
    float D_VBF_QG, D_VBF1j_QG, D_HadWH_QG, D_HadZH_QG;
    float D_bkg_VBFdec, D_bkg_VHdec;

    // a vector<float> for each vector<double>
    vector<float> lep_pt_float, lep_pterr_float, lep_pterrold_float;
    vector<float> lep_p_float, lep_ecalEnergy_float;
    vector<float> lep_eta_float, lep_phi_float, lep_mass_float;
    vector<float> lepFSR_pt_float, lepFSR_eta_float;
    vector<float> lepFSR_phi_float, lepFSR_mass_float;
    vector<float> tau_pt_float, tau_eta_float, tau_phi_float, tau_mass_float;
    vector<float> pho_pt_float, pho_eta_float, pho_phi_float, photonCutBasedIDLoose_float;
    vector<float> H_pt_float, H_eta_float, H_phi_float, H_mass_float;
    vector<float> H_noFSR_pt_float, H_noFSR_eta_float; 
    vector<float> H_noFSR_phi_float, H_noFSR_mass_float;
    vector<float> Z_pt_float, Z_eta_float, Z_phi_float, Z_mass_float;
    vector<float> Z_noFSR_pt_float, Z_noFSR_eta_float;
    vector<float> Z_noFSR_phi_float, Z_noFSR_mass_float;
    vector<float> jet_pt_float, jet_eta_float, jet_phi_float, jet_mass_float, jet_pt_raw_float;
    vector<float> jet_jesup_pt_float, jet_jesup_eta_float; 
    vector<float> jet_jesup_phi_float, jet_jesup_mass_float;
    vector<float> jet_jesdn_pt_float, jet_jesdn_eta_float;
    vector<float> jet_jesdn_phi_float, jet_jesdn_mass_float;
    vector<float> jet_jerup_pt_float, jet_jerup_eta_float;
    vector<float> jet_jerup_phi_float, jet_jerup_mass_float;
    vector<float> jet_jerdn_pt_float, jet_jerdn_eta_float;
    vector<float> jet_jerdn_phi_float, jet_jerdn_mass_float;
    vector<float> fsrPhotons_pt_float, fsrPhotons_pterr_float;
    vector<float> fsrPhotons_eta_float, fsrPhotons_phi_float, fsrPhotons_mass_float;
    vector<float> GENlep_pt_float, GENlep_eta_float;
    vector<float> GENlep_phi_float, GENlep_mass_float;
    vector<float> GENH_pt_float, GENH_eta_float;
    vector<float> GENH_phi_float, GENH_mass_float;
    vector<float> GENZ_pt_float, GENZ_eta_float;
    vector<float> GENZ_phi_float, GENZ_mass_float;
    vector<float> GENjet_pt_float, GENjet_eta_float;
    vector<float> GENjet_phi_float, GENjet_mass_float;
	
	vector<float> eleUncorr_pt_float, eleUncorr_mass_float;
	
	vector<float> track_vtx_pt_float, track_vtx_eta_float, track_vtx_phi_float, track_vtx_pterr_float;// track_vtx_ptcorr_float, track_vtx_pterrcorr_float;
//	vector<float> track_vtx_FSRpt_float, track_vtx_FSRpterr_float, track_vtx_FSRptcorr_float, track_vtx_FSRpterrcorr_float, track_vtx_FSReta_float, track_vtx_FSRphi_float, track_vtx_FSRmass_float;
	vector<float> track_vtx_dxy_float, track_vtx_dxycorr_float;
//	vector<float> track_vtx_xcorr_float, track_vtx_ycorr_float, track_vtx_zcorr_float;
//	vector<float> pv_x_float, pv_y_float, pv_z_float;//23
	
	vector<float> track_bs_pt_float, track_bs_eta_float, track_bs_phi_float, track_bs_pterr_float;//track_bs_ptcorr_float, track_bs_pterrcorr_float;
//	vector<float> track_bs_FSRpt_float, track_bs_FSRpterr_float, track_bs_FSRptcorr_float, track_bs_FSRpterrcorr_float, track_bs_FSReta_float, track_bs_FSRphi_float, track_bs_FSRmass_float;
	vector<float> track_bs_dxy_float, track_bs_dxycorr_float;
//	vector<float> track_bs_xcorr_float, track_bs_ycorr_float, track_bs_zcorr_float;
//	vector<float> bs_x_float, bs_y_float, bs_z_float;
	
	vector<float> track_Hvtx_pt_float, track_Hvtx_eta_float, track_Hvtx_phi_float, track_Hvtx_pterr_float;
	vector<float> track_Hvtx_dxycorr_float; 
	vector<float> track_Hbs_pt_float, track_Hbs_eta_float, track_Hbs_phi_float, track_Hbs_pterr_float;
	vector<float> track_Hbs_dxycorr_float;
	
	
    // Global Variables but not stored in the tree
    vector<double> lep_ptreco;
    vector<int> lep_ptid; vector<int> lep_ptindex;
    vector<pat::Muon> recoMuons; vector<pat::Electron> recoElectrons; vector<pat::Electron> recoElectronsUnS; vector<pat::Muon> recoMuonsUnS;
    vector<pat::Tau> recoTaus; vector<pat::Photon> recoPhotons;
    vector<pat::PFParticle> fsrPhotons; 
    TLorentzVector HVec, HVecNoFSR, Z1Vec, Z2Vec;
    TLorentzVector GENZ1Vec, GENZ2Vec;
    bool RecoFourMuEvent, RecoFourEEvent, RecoTwoETwoMuEvent, RecoTwoMuTwoEEvent;
    bool foundHiggsCandidate; bool foundZ1LCandidate; bool firstEntry;
    float jet1pt, jet2pt;

    // hist container
    std::map<std::string,TH1F*> histContainer_;

    //Input edm
    edm::EDGetTokenT<edm::View<pat::Electron> > elecSrc_;
    edm::EDGetTokenT<edm::View<pat::Electron> > elecUnSSrc_;
    edm::EDGetTokenT<edm::View<pat::Muon> > muonSrc_;
    edm::EDGetTokenT<edm::View<pat::Muon> > muonUnSSrc_;
    edm::EDGetTokenT<edm::View<pat::Tau> > tauSrc_;
    edm::EDGetTokenT<edm::View<pat::Photon> > photonSrc_;
    edm::EDGetTokenT<edm::View<pat::Jet> > jetSrc_;
    edm::EDGetTokenT<edm::ValueMap<float> > qgTagSrc_;
    edm::EDGetTokenT<edm::ValueMap<float> > axis2Src_;
    edm::EDGetTokenT<edm::ValueMap<int> > multSrc_;
    edm::EDGetTokenT<edm::ValueMap<float> > ptDSrc_;
    edm::EDGetTokenT<edm::View<pat::Jet> > mergedjetSrc_;
    edm::EDGetTokenT<edm::View<pat::MET> > metSrc_;
    //edm::InputTag triggerSrc_;
    edm::EDGetTokenT<edm::TriggerResults> triggerSrc_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
    edm::EDGetTokenT<std::vector<reco::Conversion> > conversionSrc_;
    edm::EDGetTokenT<double> muRhoSrc_;
    edm::EDGetTokenT<double> elRhoSrc_;
    edm::EDGetTokenT<double> rhoSrcSUS_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSrc_;
    edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandsSrc_;
    edm::EDGetTokenT<edm::View<pat::PFParticle> > fsrPhotonsSrc_;
    edm::EDGetTokenT<reco::GenParticleCollection> prunedgenParticlesSrc_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedgenParticlesSrc_;
    edm::EDGetTokenT<edm::View<reco::GenJet> > genJetsSrc_;
    edm::EDGetTokenT<GenEventInfoProduct> generatorSrc_;
    edm::EDGetTokenT<LHEEventProduct> lheInfoSrc_;
    edm::EDGetTokenT<LHERunInfoProduct> lheRunInfoToken_;
    edm::EDGetTokenT<HTXS::HiggsClassification> htxsSrc_;
    //edm::EDGetTokenT<HZZFid::FiducialSummary> fidRivetSrc_;
    edm::EDGetTokenT< double > prefweight_token_;

    // Configuration
    const float Zmass;
    float mZ1Low, mZ2Low, mZ1High, mZ2High, m4lLowCut;
    float jetpt_cut, jeteta_cut;
    std::string elecID;
    bool isMC, isSignal;
    float mH;
    float crossSection;
    bool weightEvents;
    float isoCutEl, isoCutMu; 
    double isoConeSizeEl, isoConeSizeMu;
    float sip3dCut, leadingPtCut, subleadingPtCut;
    float genIsoCutEl, genIsoCutMu;
    double genIsoConeSizeEl, genIsoConeSizeMu;
    float _elecPtCut, _muPtCut, _tauPtCut, _phoPtCut;
    float BTagCut;
    bool reweightForPU;
    std::string PUVersion;
    bool doFsrRecovery,bestCandMela, doMela, GENbestM4l;
    bool doPUJetID;
    int jetIDLevel;
    bool doJER;
    bool doJEC;
    bool doRefit;
    bool doTriggerMatching;
    bool checkOnlySingle;
    std::vector<std::string> triggerList;
    int skimLooseLeptons, skimTightLeptons;
    bool verbose;

    int year;///use to choose Muon BDT
	bool Zvtx;
	bool Hvtx;
	bool Checkvtx;
	
    // register to the TFileService
    edm::Service<TFileService> fs;

    // Counters
    float nEventsTotal;
    float sumWeightsTotal;
    float sumWeightsTotalPU;

    // JER
    JME::JetResolution resolution_pt, resolution_phi;
    JME::JetResolutionScaleFactor resolution_sf;

    string EleBDT_name_161718;
    string heepID_name_161718;
	
	//muon correct
//	KalmanMuonCalibrator *kalmanMuonCalibrator;
    RoccoR  *calibrator;

};


UFHZZ4LAna::UFHZZ4LAna(const edm::ParameterSet& iConfig) :
    histContainer_(),
    elecSrc_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc"))),
    elecUnSSrc_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("electronUnSSrc"))),
    muonSrc_(consumes<edm::View<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc"))),
	muonUnSSrc_(consumes<edm::View<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonUnSSrc"))),
    tauSrc_(consumes<edm::View<pat::Tau> >(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc"))),
    photonSrc_(consumes<edm::View<pat::Photon> >(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc"))),
    jetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc"))),
    qgTagSrc_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"))),
    axis2Src_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "axis2"))),
    multSrc_(consumes<edm::ValueMap<int>>(edm::InputTag("QGTagger", "mult"))),
    ptDSrc_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "ptD"))),
     mergedjetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("mergedjetSrc"))),
    metSrc_(consumes<edm::View<pat::MET> >(iConfig.getUntrackedParameter<edm::InputTag>("metSrc"))),
    triggerSrc_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerSrc"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
    vertexSrc_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc"))),
    beamSpotSrc_(consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc"))),
    conversionSrc_(consumes<std::vector<reco::Conversion> >(iConfig.getUntrackedParameter<edm::InputTag>("conversionSrc"))),
    muRhoSrc_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("muRhoSrc"))),
    elRhoSrc_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("elRhoSrc"))),
    rhoSrcSUS_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("rhoSrcSUS"))),
    pileupSrc_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>("pileupSrc"))),
    pfCandsSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfCandsSrc"))),
    fsrPhotonsSrc_(consumes<edm::View<pat::PFParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("fsrPhotonsSrc"))),
    prunedgenParticlesSrc_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("prunedgenParticlesSrc"))),
    packedgenParticlesSrc_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("packedgenParticlesSrc"))),
    genJetsSrc_(consumes<edm::View<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("genJetsSrc"))),
    generatorSrc_(consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("generatorSrc"))),
    lheInfoSrc_(consumes<LHEEventProduct>(iConfig.getUntrackedParameter<edm::InputTag>("lheInfoSrc"))),
    lheRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer",""))),
    htxsSrc_(consumes<HTXS::HiggsClassification>(edm::InputTag("rivetProducerHTXS","HiggsClassification"))),
    prefweight_token_(consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"))),//2017 2017 need this. not for 2018
    //fidRivetSrc_(consumes<HZZFid::FiducialSummary>(edm::InputTag("rivetProducerHZZFid","FiducialSummary"))),
    Zmass(91.1876),
    mZ1Low(iConfig.getUntrackedParameter<double>("mZ1Low",40.0)),
    mZ2Low(iConfig.getUntrackedParameter<double>("mZ2Low",12.0)), // was 12
    mZ1High(iConfig.getUntrackedParameter<double>("mZ1High",120.0)),
    mZ2High(iConfig.getUntrackedParameter<double>("mZ2High",120.0)),
    m4lLowCut(iConfig.getUntrackedParameter<double>("m4lLowCut",70.0)),
    jetpt_cut(iConfig.getUntrackedParameter<double>("jetpt_cut",10.0)),
    jeteta_cut(iConfig.getUntrackedParameter<double>("eta_cut",4.7)),
    elecID(iConfig.getUntrackedParameter<std::string>("elecID","NonTrig")),
    isMC(iConfig.getUntrackedParameter<bool>("isMC",true)),
    isSignal(iConfig.getUntrackedParameter<bool>("isSignal",false)),
    mH(iConfig.getUntrackedParameter<double>("mH",0.0)),
    crossSection(iConfig.getUntrackedParameter<double>("CrossSection",1.0)),
    weightEvents(iConfig.getUntrackedParameter<bool>("weightEvents",false)),
    isoCutEl(iConfig.getUntrackedParameter<double>("isoCutEl",9999.0)),
    isoCutMu(iConfig.getUntrackedParameter<double>("isoCutMu",0.35)),/////ios is applied to new Muon BDT //previous 0.35///Qianying
    //isoCutMu(iConfig.getUntrackedParameter<double>("isoCutMu",9999.0)),/////ios is applied to new Muon BDT //previous 0.35///Qianying
    isoConeSizeEl(iConfig.getUntrackedParameter<double>("isoConeSizeEl",0.3)),
    isoConeSizeMu(iConfig.getUntrackedParameter<double>("isoConeSizeMu",0.3)),
    sip3dCut(iConfig.getUntrackedParameter<double>("sip3dCut",4)),
    leadingPtCut(iConfig.getUntrackedParameter<double>("leadingPtCut",20.0)),
    subleadingPtCut(iConfig.getUntrackedParameter<double>("subleadingPtCut",10.0)),
    genIsoCutEl(iConfig.getUntrackedParameter<double>("genIsoCutEl",0.35)), 
    genIsoCutMu(iConfig.getUntrackedParameter<double>("genIsoCutMu",0.35)), 
    genIsoConeSizeEl(iConfig.getUntrackedParameter<double>("genIsoConeSizeEl",0.3)), 
    genIsoConeSizeMu(iConfig.getUntrackedParameter<double>("genIsoConeSizeMu",0.3)), 
    _elecPtCut(iConfig.getUntrackedParameter<double>("_elecPtCut",7.0)),
    _muPtCut(iConfig.getUntrackedParameter<double>("_muPtCut",5.0)),
    _tauPtCut(iConfig.getUntrackedParameter<double>("_tauPtCut",20.0)),
    _phoPtCut(iConfig.getUntrackedParameter<double>("_phoPtCut",10.0)),
    //BTagCut(iConfig.getUntrackedParameter<double>("BTagCut",0.4184)),/////2016: 0.6321; 2017: 0.4941; 2018: 0.4184
    reweightForPU(iConfig.getUntrackedParameter<bool>("reweightForPU",true)),
    PUVersion(iConfig.getUntrackedParameter<std::string>("PUVersion","Summer16_80X")),
    doFsrRecovery(iConfig.getUntrackedParameter<bool>("doFsrRecovery",true)),
    bestCandMela(iConfig.getUntrackedParameter<bool>("bestCandMela",true)),
    doMela(iConfig.getUntrackedParameter<bool>("doMela",true)),
    GENbestM4l(iConfig.getUntrackedParameter<bool>("GENbestM4l",false)),
    doPUJetID(iConfig.getUntrackedParameter<bool>("doPUJetID",true)),
    jetIDLevel(iConfig.getUntrackedParameter<int>("jetIDLevel",2)),
    doJER(iConfig.getUntrackedParameter<bool>("doJER",true)),
    doJEC(iConfig.getUntrackedParameter<bool>("doJEC",true)),
    doRefit(iConfig.getUntrackedParameter<bool>("doRefit",true)),
    doTriggerMatching(iConfig.getUntrackedParameter<bool>("doTriggerMatching",!isMC)),
    checkOnlySingle(iConfig.getUntrackedParameter<bool>("checkOnlySingle",false)),
    triggerList(iConfig.getUntrackedParameter<std::vector<std::string>>("triggerList")),
    skimLooseLeptons(iConfig.getUntrackedParameter<int>("skimLooseLeptons",2)),    
    skimTightLeptons(iConfig.getUntrackedParameter<int>("skimTightLeptons",2)),    
    verbose(iConfig.getUntrackedParameter<bool>("verbose",false)),
    year(iConfig.getUntrackedParameter<int>("year",2018)),////for year put 2016,2017, or 2018 to select correct training
	Zvtx(iConfig.getUntrackedParameter<bool>("Zvtx",false)),
	Hvtx(iConfig.getUntrackedParameter<bool>("Hvtx",false)),
	Checkvtx(iConfig.getUntrackedParameter<bool>("Checkvtx",false))

{
  
    if(!isMC){reweightForPU = false;}

    nEventsTotal=0.0;
    sumWeightsTotal=0.0;
    sumWeightsTotalPU=0.0;
    histContainer_["NEVENTS"]=fs->make<TH1F>("nEvents","nEvents in Sample",2,0,2);
    histContainer_["SUMWEIGHTS"]=fs->make<TH1F>("sumWeights","sum Weights of Sample",2,0,2);
    histContainer_["SUMWEIGHTSPU"]=fs->make<TH1F>("sumWeightsPU","sum Weights and PU of Sample",2,0,2);
    histContainer_["NVTX"]=fs->make<TH1F>("nVtx","Number of Vertices",36,-0.5,35.5);
    histContainer_["NVTX_RW"]=fs->make<TH1F>("nVtx_ReWeighted","Number of Vertices",36,-0.5,35.5);
    histContainer_["NINTERACT"]=fs->make<TH1F>("nInteractions","Number of True Interactions",61,-0.5,60.5);
    histContainer_["NINTERACT_RW"]=fs->make<TH1F>("nInteraction_ReWeighted","Number of True Interactions",61,-0.5,60.5);

    passedEventsTree_All = new TTree("passedEvents","passedEvents");

    edm::FileInPath kfacfileInPath("UFHZZAnalysisRun2/UFHZZ4LAna/data/Kfactor_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
    TFile *fKFactor = TFile::Open(kfacfileInPath.fullPath().c_str());
    kFactor_ggzz = (TSpline3*) fKFactor->Get("sp_Kfactor");
    fKFactor->Close();
    delete fKFactor;

    tableEwk = readFile_and_loadEwkTable("ZZBG");   
   
    kinZfitter = new KinZfitter(!isMC);
    
   if(doMela){
    mela = new Mela(13.0, 125.0, TVar::SILENT);
    mela->setCandidateDecayMode(TVar::CandidateDecay_ZZ); 
	}

    //string elec_scalefac_Cracks_name_161718[3] = {"egammaEffi.txt_EGM2D_cracks.root", "egammaEffi.txt_EGM2D_Moriond2018v1_gap.root", "egammaEffi.txt_EGM2D_Moriond2019_v1_gap.root"};
    string elec_scalefac_Cracks_name_161718[3] = {"ElectronSF_Legacy_2016_Gap.root", "ElectronSF_Legacy_2017_Gap.root", "ElectronSF_Legacy_2018_Gap.root"};
    edm::FileInPath elec_scalefacFileInPathCracks(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+elec_scalefac_Cracks_name_161718[year-2016]).c_str());
    TFile *fElecScalFacCracks = TFile::Open(elec_scalefacFileInPathCracks.fullPath().c_str());
    hElecScaleFac_Cracks = (TH2F*)fElecScalFacCracks->Get("EGamma_SF2D");    

    //string elec_scalefac_name_161718[3] = {"egammaEffi.txt_EGM2D.root", "egammaEffi.txt_EGM2D_Moriond2018v1.root", "egammaEffi.txt_EGM2D_Moriond2019_v1.root"};
    string elec_scalefac_name_161718[3] = {"ElectronSF_Legacy_2016_NoGap.root", "ElectronSF_Legacy_2017_NoGap.root", "ElectronSF_Legacy_2018_NoGap.root"};
    edm::FileInPath elec_scalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+elec_scalefac_name_161718[year-2016]).c_str());
    TFile *fElecScalFac = TFile::Open(elec_scalefacFileInPath.fullPath().c_str());
    hElecScaleFac = (TH2F*)fElecScalFac->Get("EGamma_SF2D");    

    //string elec_Gsfscalefac_name_161718[3] = {"egammaEffi.txt_EGM2D_GSF.root", "egammaEffi.txt_EGM2D_Moriond2018v1_runBCDEF_passingRECO.root", "Ele_Reco_2018.root"};//was previous;
    string elec_Gsfscalefac_name_161718[3] = {"Ele_Reco_2016.root", "Ele_Reco_2017.root", "Ele_Reco_2018.root"};
    edm::FileInPath elec_GsfscalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+elec_Gsfscalefac_name_161718[year-2016]).c_str());
    TFile *fElecScalFacGsf = TFile::Open(elec_GsfscalefacFileInPath.fullPath().c_str());
    hElecScaleFacGsf = (TH2F*)fElecScalFacGsf->Get("EGamma_SF2D");

    //string elec_GsfLowETscalefac_name_161718[3]= {"", "egammaEffi.txt_EGM2D_Moriond2018v1_runBCDEF_passingRECO_lowEt.root", "Ele_Reco_LowEt_2018.root"};//was previous
    string elec_GsfLowETscalefac_name_161718[3]= {"Ele_Reco_LowEt_2016.root", "Ele_Reco_LowEt_2017.root", "Ele_Reco_LowEt_2018.root"};
    edm::FileInPath elec_GsfLowETscalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+elec_GsfLowETscalefac_name_161718[year-2016]).c_str());
    TFile *fElecScalFacGsfLowET = TFile::Open(elec_GsfLowETscalefacFileInPath.fullPath().c_str());
    hElecScaleFacGsfLowET = (TH2F*)fElecScalFacGsfLowET->Get("EGamma_SF2D");

    //string mu_scalefac_name_161718[3] = {"final_HZZ_Moriond17Preliminary_v4.root", "ScaleFactors_mu_Moriond2018_final.root", "final_HZZ_muon_SF_2018RunA2D_ER_2702.root"};//was previous; 
    string mu_scalefac_name_161718[3] = {"final_HZZ_SF_2016_legacy_mupogsysts.root", "final_HZZ_SF_2017_rereco_mupogsysts_3010.root", "final_HZZ_SF_2018_rereco_mupogsysts_3010.root"};
    edm::FileInPath mu_scalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+mu_scalefac_name_161718[year-2016]).c_str());
    TFile *fMuScalFac = TFile::Open(mu_scalefacFileInPath.fullPath().c_str());
    hMuScaleFac = (TH2F*)fMuScalFac->Get("FINAL");
    hMuScaleFacUnc = (TH2F*)fMuScalFac->Get("ERROR");

    //string pileup_name_161718[3] = {"puWeightsMoriond17_v2.root", "puWeightsMoriond18.root", "pu_weights_2018.root"};///was previous
    string pileup_name_161718[3] = {"pu_weights_2016.root", "pu_weights_2017.root", "pu_weights_2018.root"};
    edm::FileInPath pileup_FileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+pileup_name_161718[year-2016]).c_str());
    TFile *f_pileup = TFile::Open(pileup_FileInPath.fullPath().c_str());
    h_pileup = (TH1D*)f_pileup->Get("weights");
    h_pileupUp = (TH1D*)f_pileup->Get("weights_varUp");
    h_pileupDn = (TH1D*)f_pileup->Get("weights_varDn");


    string bTagEffi_name_161718[3] = {"bTagEfficiencies_2016.root", "bTagEfficiencies_2017.root", "bTagEfficiencies_2018.root"};
    edm::FileInPath BTagEffiInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+bTagEffi_name_161718[year-2016]).c_str());
    TFile *fbTagEffi = TFile::Open(BTagEffiInPath.fullPath().c_str());
    hbTagEffi = (TH2F*)fbTagEffi->Get("eff_b_M_ALL");
    hcTagEffi = (TH2F*)fbTagEffi->Get("eff_c_M_ALL");
    hudsgTagEffi = (TH2F*)fbTagEffi->Get("eff_udsg_M_ALL");


    //BTag calibration
    string csv_name_161718[3] = {"DeepCSV_2016LegacySF_V1.csv", "DeepCSV_94XSF_V4_B_F.csv", "DeepCSV_102XSF_V1.csv"};
    edm::FileInPath btagfileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+csv_name_161718[year-2016]).c_str());

    BTagCalibration calib("DeepCSV", btagfileInPath.fullPath().c_str());
    reader = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,  // operating point
                                       "central",             // central sys type
                                       {"up", "down"});      // other sys types
   

    reader->load(calib,                // calibration instance
                BTagEntry::FLAV_B,    // btag flavour
                "comb");               // measurement type

    if(year==2018)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Autumn18IdIsoValues"; BTagCut=0.4184; heepID_name_161718 = "heepElectronID-HEEPV70";}
    if(year==2017)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Fall17IsoV2Values"; BTagCut=0.4941; heepID_name_161718 = "heepElectronID-HEEPV70";}
    if(year==2016)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Summer16IdIsoValues"; BTagCut=0.6321; heepID_name_161718 = "heepElectronID-HEEPV70";}
	

//	std::string name;
//	if(isMC) name = "MC_80X_13TeV";
//	if(!isMC) name = "DATA_80X_13TeV";
//	kalmanMuonCalibrator = new KalmanMuonCalibrator(name);
	std::string DATAPATH = std::getenv( "CMSSW_BASE" );
	if(year == 2018) DATAPATH+="/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v3/RoccoR2018.txt";
	if(year == 2017) DATAPATH+="/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v3/RoccoR2017.txt";
	if(year == 2016) DATAPATH+="/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v3/RoccoR2016.txt";
	calibrator = new RoccoR(DATAPATH);
	
}



UFHZZ4LAna::~UFHZZ4LAna()
{
    //destructor --- don't do anything here  
}



// ------------ method called for each event  ------------
void
UFHZZ4LAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    using namespace edm;
    using namespace std;
    using namespace pat;
    using namespace trigger;
    using namespace EwkCorrections;

    nEventsTotal += 1.0;

    Run = iEvent.id().run();
    Event = iEvent.id().event();
    LumiSect = iEvent.id().luminosityBlock();

    if (verbose) {
       cout<<"Run: " << Run << ",Event: " << Event << ",LumiSect: "<<LumiSect<<endl;
    }

    // ======= Get Collections ======= //
    if (verbose) {cout<<"getting collections"<<endl;}

    // trigger collection
    edm::Handle<edm::TriggerResults> trigger;
    iEvent.getByToken(triggerSrc_,trigger);
    const edm::TriggerNames trigNames = iEvent.triggerNames(*trigger);

    // trigger Objects
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects);

    // vertex collection
    edm::Handle<reco::VertexCollection> vertex;
    iEvent.getByToken(vertexSrc_,vertex);

    // electron collection
    edm::Handle<edm::View<pat::Electron> > electrons;
    iEvent.getByToken(elecSrc_,electrons);
    if (verbose) cout<<electrons->size()<<" total electrons in the collection"<<endl;

    // electron before scale/smearing corrections
    edm::Handle<edm::View<pat::Electron> > electronsUnS;
    iEvent.getByToken(elecUnSSrc_,electronsUnS);

    // muon collection
    edm::Handle<edm::View<pat::Muon> > muons;
    iEvent.getByToken(muonSrc_,muons);
    if (verbose) cout<<muons->size()<<" total muons in the collection"<<endl;
	
	// muon before roccor
	edm::Handle<edm::View<pat::Muon> > muonsUnS;
	iEvent.getByToken(muonUnSSrc_,muonsUnS);

    // tau collection
    edm::Handle<edm::View<pat::Tau> > taus;
    iEvent.getByToken(tauSrc_,taus);
    if (verbose) cout<<taus->size()<<" total taus in the collection"<<endl;

    // photon collection 
    edm::Handle<edm::View<pat::Photon> > photons;
    iEvent.getByToken(photonSrc_,photons);
    if (verbose) cout<<photons->size()<<" total photons in the collection"<<endl;
  
    // met collection 
    edm::Handle<edm::View<pat::MET> > mets;
    iEvent.getByToken(metSrc_,mets);
    
    // Rho Correction
    edm::Handle<double> eventRhoMu;
    iEvent.getByToken(muRhoSrc_,eventRhoMu);
    muRho = *eventRhoMu;

    edm::Handle<double> eventRhoE;
    iEvent.getByToken(elRhoSrc_,eventRhoE);
    elRho = *eventRhoE;

    edm::Handle<double> eventRhoSUS;
    iEvent.getByToken(rhoSrcSUS_,eventRhoSUS);
    rhoSUS = *eventRhoSUS;

    // Conversions
    edm::Handle< std::vector<reco::Conversion> > theConversions;
    iEvent.getByToken(conversionSrc_, theConversions);
 
    // Beam Spot
    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(beamSpotSrc_,beamSpot);
    const reco::BeamSpot BS = *beamSpot;

    // Particle Flow Cands
    edm::Handle<pat::PackedCandidateCollection> pfCands;
    iEvent.getByToken(pfCandsSrc_,pfCands);

    // FSR Photons
    edm::Handle<edm::View<pat::PFParticle> > photonsForFsr;
    iEvent.getByToken(fsrPhotonsSrc_,photonsForFsr);
  
    // Jets
    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByToken(jetSrc_,jets);

    if (!jecunc) {
        edm::ESHandle<JetCorrectorParametersCollection> jetCorrParameterSet;
        iSetup.get<JetCorrectionsRecord>().get("AK4PFchs", jetCorrParameterSet);
        const JetCorrectorParameters& jetCorrParameters = (*jetCorrParameterSet)["Uncertainty"];
        jecunc.reset(new JetCorrectionUncertainty(jetCorrParameters));
    }

    resolution_pt = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
    resolution_phi = JME::JetResolution::get(iSetup, "AK4PFchs_phi");
    resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");

    edm::Handle<edm::ValueMap<float> > qgHandle;
    iEvent.getByToken(qgTagSrc_, qgHandle);

    edm::Handle<edm::ValueMap<float> > axis2Handle;
    iEvent.getByToken(axis2Src_, axis2Handle);

    edm::Handle<edm::ValueMap<int> > multHandle;
    iEvent.getByToken(multSrc_, multHandle);

    edm::Handle<edm::ValueMap<float> > ptDHandle;
    iEvent.getByToken(ptDSrc_, ptDHandle);
 
    edm::Handle<edm::View<pat::Jet> > mergedjets;
    iEvent.getByToken(mergedjetSrc_,mergedjets);

    // GEN collections
    edm::Handle<reco::GenParticleCollection> prunedgenParticles;
    iEvent.getByToken(prunedgenParticlesSrc_, prunedgenParticles);

    edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles;
    iEvent.getByToken(packedgenParticlesSrc_, packedgenParticles);
    
    edm::Handle<edm::View<reco::GenJet> > genJets;
    iEvent.getByToken(genJetsSrc_, genJets);
    
    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken(generatorSrc_,genEventInfo);
    
    //vector<edm::Handle<LHEEventProduct> > lheInfos;
    //iEvent.getManyByType(lheInfos); // using this method because the label is not always the same (e.g. "source" in the ttH sample)

    edm::Handle<LHEEventProduct> lheInfo;
    iEvent.getByToken(lheInfoSrc_, lheInfo);

    // STXS info    
    if (isMC) {
        edm::Handle<HTXS::HiggsClassification> htxs;
        iEvent.getByToken(htxsSrc_,htxs);
        stage0cat = htxs->stage0_cat;
        stage1cat = htxs->stage1_cat_pTjet30GeV;
        //stage1p1cat = htxs->stage1p1_cat;
        stage1p1cat = htxs->stage1_1_cat_pTjet30GeV;
        stage1p2cat = htxs->stage1_2_cat_pTjet30GeV;
        if (verbose) cout<<"stage1cat "<<stage1cat<<endl;
    }

    // Fiducial Rivet
    /*
     passedFiducialRivet = false;
     GENpT4lRivet = -1.0;
     GENnjets_pt30_eta4p7Rivet = -1;
     GENpt_leadingjet_pt30_eta4p7Rivet = -1.0;
     edm::Handle<HZZFid::FiducialSummary> rivetfid;
     iEvent.getByToken(fidRivetSrc_,rivetfid);
     passedFiducialRivet = rivetfid->passedFiducial;
     GENpT4lRivet = rivetfid->higgs.Pt();
     GENnjets_pt30_eta4p7Rivet = rivetfid->jets.size();
     if (rivetfid->jets.size()>0) GENpt_leadingjet_pt30_eta4p7Rivet = rivetfid->jets[0].Pt();
    */

    if (isMC) {
        edm::Handle< double > theprefweight;
        iEvent.getByToken(prefweight_token_, theprefweight ) ;
        if (year == 2016 || year == 2017)
            prefiringWeight =(*theprefweight);
        else if (year == 2018)
            prefiringWeight =1.0;
    }
    else 
        prefiringWeight =1.0;

    // ============ Initialize Variables ============= //

    // Event Variables
    if (verbose) {cout<<"clear variables"<<endl;}
    nVtx = -1.0; nInt = -1.0;
    finalState = -1;
    triggersPassed="";
    passedTrig=false; passedFullSelection=false; passedZ4lSelection=false; passedQCDcut=false; 
    passedZ1LSelection=false; passedZ4lZ1LSelection=false; passedZ4lZXCRSelection=false; passedZXCRSelection=false;
    nZXCRFailedLeptons=0;


    // Event Weights
    genWeight=1.0; pileupWeight=1.0; pileupWeightUp=1.0; pileupWeightDn=1.0; dataMCWeight=1.0; eventWeight=1.0;
    k_ggZZ=1.0; k_qqZZ_qcd_dPhi = 1.0; k_qqZZ_qcd_M = 1.0; k_qqZZ_qcd_Pt = 1.0; k_qqZZ_ewk = 1.0;

    qcdWeights.clear(); nnloWeights.clear(); pdfWeights.clear();
    pdfRMSup=1.0; pdfRMSdown=1.0; pdfENVup=1.0; pdfENVdown=1.0;

    //lepton variables
    lep_pt.clear(); lep_pterr.clear(); lep_pterrold.clear(); 
    lep_p.clear(); lep_ecalEnergy.clear(); lep_isEB.clear(); lep_isEE.clear();
    lep_eta.clear(); lep_phi.clear(); lep_mass.clear(); 
    lepFSR_pt.clear(); lepFSR_eta.clear(); lepFSR_phi.clear(); lepFSR_mass.clear(); 
    for (int i=0; i<4; ++i) {lep_Hindex[i]=-1;}
    pTL1=-1.0; pTL2=-1.0; pTL3=-1.0; pTL4=-1.0;
    etaL1=9999.0; etaL2=9999.0; etaL3=9999.0; etaL4=9999.0;
    idL1=9999; idL2=9999; idL3=9999; idL4=9999;
    mL1=-1.0; mL2=-1.0; mL3=-1.0; mL4=-1.0;
    pTErrL1=-1.0; pTErrL2=-1.0; pTErrL3=-1.0; pTErrL4=-1.0;
    phiL1=9999.0; phiL2=9999.0; phiL3=9999.0; phiL4=9999.0;
    pTL1FSR=-1.0; pTL2FSR=-1.0; pTL3FSR=-1.0; pTL4FSR=-1.0;
    lep_genindex.clear(); lep_id.clear(); lep_dataMC.clear(); lep_dataMCErr.clear();
    lep_matchedR03_PdgId.clear(); lep_matchedR03_MomId.clear(); lep_matchedR03_MomMomId.clear();
    lep_mva.clear(); lep_ecalDriven.clear(); 
    lep_tightId.clear(); lep_tightIdSUS.clear(); lep_tightIdHiPt.clear(); //lep_tightId_old.clear();
    lep_Sip.clear(); lep_IP.clear(); 
    lep_isoNH.clear(); lep_isoCH.clear(); lep_isoPhot.clear(); lep_isoPU.clear(); lep_isoPUcorr.clear(); 
    lep_RelIso.clear(); lep_RelIsoNoFSR.clear(); lep_MiniIso.clear();
    lep_ptRatio.clear(); lep_ptRel.clear();
    lep_missingHits.clear();
    lep_filtersMatched.clear();    
    nisoleptons=0;
	
	eleUncorr_pt.clear(); eleUncorr_mass.clear();
	
	track_vtx_pt.clear(); track_vtx_eta.clear(); track_vtx_phi.clear(); track_vtx_pterr.clear();// track_vtx_ptcorr.clear(); track_vtx_pterrcorr.clear();
//	track_vtx_FSRpt.clear(); track_vtx_FSRpterr.clear(); track_vtx_FSRptcorr.clear(); track_vtx_FSRpterrcorr.clear(); track_vtx_FSReta.clear(); track_vtx_FSRphi.clear(); track_vtx_FSRmass.clear();
	track_vtx_dxy.clear(); track_vtx_dxycorr.clear();
	track_vtx_xcorr = -999.9; track_vtx_ycorr = -999.9; track_vtx_zcorr = -999.9;
	pv_x = -999.9; pv_y = -999.9; pv_z = -999.9;
	
	track_bs_pt.clear(); track_bs_eta.clear(); track_bs_phi.clear(); track_bs_pterr.clear(); //track_bs_ptcorr.clear(); track_bs_pterrcorr.clear();//6
//	track_bs_FSRpt.clear(); track_bs_FSRpterr.clear(); track_bs_FSRptcorr.clear(); track_bs_FSRpterrcorr.clear(); track_vtx_FSReta.clear(); track_vtx_FSRphi.clear(); track_vtx_FSRmass.clear();//7
	track_bs_dxycorr.clear(); track_bs_dxy.clear();//4
	track_bs_xcorr = -999.9; track_bs_ycorr = -999.9; track_bs_zcorr = -999.9;//3
	bs_x = -999.9; bs_y = -999.9; bs_z = -999.9;//3
	
	track_Hvtx_pt.clear(); track_Hvtx_eta.clear(); track_Hvtx_phi.clear(); track_Hvtx_pterr.clear();
	track_Hvtx_dxycorr.clear();
	track_Hbs_pt.clear(); track_Hbs_eta.clear(); track_Hbs_phi.clear(); track_Hbs_pterr.clear();
	track_Hbs_dxycorr.clear();
	
	track_Hvtx_xcorr = -999.0; track_Hvtx_ycorr = -999.0; track_Hvtx_zcorr = -999.0;
	track_Hbs_xcorr = -999.0; track_Hbs_ycorr = -999.0; track_Hbs_xcorr = -999.0;
	
	mu_pt_roc.clear(); mu_eta_roc.clear(); mu_phi_roc.clear(); mu_q_roc.clear(); mu_tracktype_roc.clear(); mu_nl_roc.clear(); mu_genparticle_roc.clear(); mu_gen_pt_roc.clear();

	passLepindex.clear();
	H4Lrefittedtracks.clear();
	tracks_tmp.clear();
	
    //tau variables
    tau_id.clear(); tau_pt.clear(); tau_eta.clear(); tau_phi.clear(); tau_mass.clear(); 

    // photon variables
    pho_pt.clear(); pho_eta.clear(); pho_phi.clear(); photonCutBasedIDLoose.clear(); 

    // Higgs candidate variables
    H_pt.clear(); H_eta.clear(); H_phi.clear(); H_mass.clear(); 
    H_noFSR_pt.clear(); H_noFSR_eta.clear(); H_noFSR_phi.clear(); H_noFSR_mass.clear(); 
    mass4l=-1.0; mass4l_noFSR=-1.0; mass4e=-1.0; mass4mu=-1.0; mass2e2mu=-1.0; pT4l=-1.0; eta4l=9999.0; phi4l=9999.0; rapidity4l=9999.0;
    cosTheta1=9999.0; cosTheta2=9999.0; cosThetaStar=9999.0; Phi=9999.0; Phi1=9999.0;
    mass3l=-1.0;

    // kin fitter
    mass4lREFIT = -999.0; massZ1REFIT = -999.0; massZ2REFIT = -999.0; mass4lErr = -999.0; mass4lErrREFIT = -999.0;

    // Z candidate variables
    Z_pt.clear(); Z_eta.clear(); Z_phi.clear(); Z_mass.clear(); 
    Z_noFSR_pt.clear(); Z_noFSR_eta.clear(); Z_noFSR_phi.clear(); Z_noFSR_mass.clear(); 
    for (int i=0; i<2; ++i) {Z_Hindex[i]=-1;}
    massZ1=-1.0; massZ1_Z1L=-1.0; massZ2=-1.0; pTZ1=-1.0; pTZ2=-1.0;

    // MET
    met=-1.0; met_phi=9999.0;
    met_jesup=-1.0; met_phi_jesup=9999.0; met_jesdn=-1.0; met_phi_jesdn=9999.0; 
    met_uncenup=-1.0; met_phi_uncenup=9999.0; met_uncendn=-1.0; met_phi_uncendn=9999.0; 

    // Jets
    jet_pt.clear(); jet_eta.clear(); jet_phi.clear(); jet_mass.clear(); jet_pt_raw.clear(); 
    jet_jesup_pt.clear(); jet_jesup_eta.clear(); jet_jesup_phi.clear(); jet_jesup_mass.clear(); 
    jet_jesdn_pt.clear(); jet_jesdn_eta.clear(); jet_jesdn_phi.clear(); jet_jesdn_mass.clear(); 
    jet_jerup_pt.clear(); jet_jerup_eta.clear(); jet_jerup_phi.clear(); jet_jerup_mass.clear(); 
    jet_jerdn_pt.clear(); jet_jerdn_eta.clear(); jet_jerdn_phi.clear(); jet_jerdn_mass.clear(); 
    jet_pumva.clear(); jet_csvv2.clear(); jet_isbtag.clear();
    jet_csvv2_.clear();
    jet_hadronFlavour.clear(); jet_partonFlavour.clear();
    jet_QGTagger.clear(); jet_QGTagger_jesup.clear(); jet_QGTagger_jesdn.clear(); 
    jet_relpterr.clear(); jet_phierr.clear();
    jet_bTagEffi.clear();
    jet_cTagEffi.clear();
    jet_udsgTagEffi.clear();
    jet_axis2.clear(); jet_ptD.clear(); jet_mult.clear();

    jet_iscleanH4l.clear();
    jet1index=-1; jet2index=-1;
    jet_jesup_iscleanH4l.clear(); jet_jesdn_iscleanH4l.clear(); 
    jet_jerup_iscleanH4l.clear(); jet_jerdn_iscleanH4l.clear();

    njets_pt30_eta4p7=0;
    njets_pt30_eta4p7_jesup=0; njets_pt30_eta4p7_jesdn=0;
    njets_pt30_eta4p7_jerup=0; njets_pt30_eta4p7_jerdn=0;

    njets_pt30_eta2p5=0;
    njets_pt30_eta2p5_jesup=0; njets_pt30_eta2p5_jesdn=0;
    njets_pt30_eta2p5_jerup=0; njets_pt30_eta2p5_jerdn=0;

    nbjets_pt30_eta4p7=0; nvjets_pt40_eta2p4=0;

    pt_leadingjet_pt30_eta4p7=-1.0;
    pt_leadingjet_pt30_eta4p7_jesup=-1.0; pt_leadingjet_pt30_eta4p7_jesdn=-1.0;
    pt_leadingjet_pt30_eta4p7_jerup=-1.0; pt_leadingjet_pt30_eta4p7_jerdn=-1.0;

    pt_leadingjet_pt30_eta2p5=-1.0;
    pt_leadingjet_pt30_eta2p5_jesup=-1.0; pt_leadingjet_pt30_eta2p5_jesdn=-1.0;
    pt_leadingjet_pt30_eta2p5_jerup=-1.0; pt_leadingjet_pt30_eta2p5_jerdn=-1.0;

    absrapidity_leadingjet_pt30_eta4p7=-1.0;
    absrapidity_leadingjet_pt30_eta4p7_jesup=-1.0; absrapidity_leadingjet_pt30_eta4p7_jesdn=-1.0;
    absrapidity_leadingjet_pt30_eta4p7_jerup=-1.0; absrapidity_leadingjet_pt30_eta4p7_jerdn=-1.0;

    absdeltarapidity_hleadingjet_pt30_eta4p7=-1.0;
    absdeltarapidity_hleadingjet_pt30_eta4p7_jesup=-1.0; absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn=-1.0;
    absdeltarapidity_hleadingjet_pt30_eta4p7_jerup=-1.0; absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn=-1.0;

    DijetMass=-1.0; DijetDEta=9999.0; DijetFisher=9999.0;
    
    mergedjet_iscleanH4l.clear();
    mergedjet_pt.clear(); mergedjet_eta.clear(); mergedjet_phi.clear(); mergedjet_mass.clear();
    mergedjet_L1.clear();
    mergedjet_softdropmass.clear(); mergedjet_prunedmass.clear();
    mergedjet_tau1.clear(); mergedjet_tau2.clear();
    mergedjet_btag.clear();

    mergedjet_nsubjet.clear();
    mergedjet_subjet_pt.clear(); mergedjet_subjet_eta.clear(); 
    mergedjet_subjet_phi.clear(); mergedjet_subjet_mass.clear();
    mergedjet_subjet_btag.clear();
    mergedjet_subjet_partonFlavour.clear(); mergedjet_subjet_hadronFlavour.clear();
    
    // FSR Photons
    nFSRPhotons=0;
    fsrPhotons_lepindex.clear(); fsrPhotons_pt.clear(); fsrPhotons_pterr.clear(); 
    fsrPhotons_eta.clear(); fsrPhotons_phi.clear();
    fsrPhotons_dR.clear(); fsrPhotons_iso.clear();
    allfsrPhotons_dR.clear(); allfsrPhotons_pt.clear(); allfsrPhotons_iso.clear();

    // Z4l? FIXME
    theta12=9999.0; theta13=9999.0; theta14=9999.0;
    minM3l=-1.0; Z4lmaxP=-1.0; minDeltR=9999.0; m3l_soft=-1.0;
    minMass2Lep=-1.0; maxMass2Lep=-1.0;
    thetaPhoton=9999.0; thetaPhotonZ=9999.0;

    // -------------------------
    // GEN level information
    // ------------------------- 

    //Event variables
    GENfinalState=-1;
    passedFiducialSelection=false;

    // lepton variables
    GENlep_pt.clear(); GENlep_eta.clear(); GENlep_phi.clear(); GENlep_mass.clear();
    GENlep_id.clear(); GENlep_status.clear(); GENlep_MomId.clear(); GENlep_MomMomId.clear();
    for (int i=0; i<4; ++i) {GENlep_Hindex[i]=-1;};//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z3 sub
    GENlep_isoCH.clear(); GENlep_isoNH.clear(); GENlep_isoPhot.clear(); GENlep_RelIso.clear();

    // Higgs candidate variables (calculated using selected gen leptons)
    GENH_pt.clear(); GENH_eta.clear(); GENH_phi.clear(); GENH_mass.clear();
    GENmass4l=-1.0; GENmassZ1=-1.0; GENmassZ2=-1.0; GENpT4l=-1.0; GENeta4l=9999.0; GENrapidity4l=9999.0; GENMH=-1.0;
    GENcosTheta1=9999.0; GENcosTheta2=9999.0; GENcosThetaStar=9999.0; GENPhi=9999.0; GENPhi1=9999.0;

    // Z candidate variables
    GENZ_DaughtersId.clear(); GENZ_MomId.clear();
    GENZ_pt.clear(); GENZ_eta.clear(); GENZ_phi.clear(); GENZ_mass.clear();
    GENmassZ1=-1.0; GENmassZ2=-1.0; GENpTZ1=-1.0; GENpTZ2=-1.0, GENdPhiZZ=9999.0, GENmassZZ=-1.0, GENpTZZ=-1.0;

    // Higgs variables directly from GEN particle
    GENHmass=-1.0;

    // Jets
    GENjet_pt.clear(); GENjet_eta.clear(); GENjet_phi.clear(); GENjet_mass.clear(); 
    GENnjets_pt30_eta4p7=0;
    GENnjets_pt30_eta2p5=0;
    GENpt_leadingjet_pt30_eta4p7=-1.0; GENabsrapidity_leadingjet_pt30_eta4p7=-1.0; GENabsdeltarapidity_hleadingjet_pt30_eta4p7=-1.0;
    GENpt_leadingjet_pt30_eta2p5=-1.0; 
    lheNb=0; lheNj=0; nGenStatus2bHad=0;

    // ME
    me_0plus_JHU=999.0; me_qqZZ_MCFM=999.0; p0plus_m4l=999.0; bkg_m4l=999.0; D_bkg_kin=999.0; D_bkg=999.0;   

    p0minus_VAJHU=999.0; pg1g4_VAJHU=999.0; Dgg10_VAMCFM=999.0, D_g4=999.0; D_g1g4=999.0;
    // OLD but working
    phjj_VAJHU=999.0; pvbf_VAJHU=999.0; pAux_vbf_VAJHU=999.0;
    pwh_hadronic_VAJHU=999.0; pwh_hadronic_VAJHU=999.0; 
    pzh_hadronic_VAJHU=999.0; pzh_hadronic_VAJHU=999.0; 

    p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal=999.0;
    p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal=999.0;
    p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal=999.0;
    p_JJVBF_BKG_MCFM_JECNominal=999.0;
    p_HadZH_BKG_MCFM_JECNominal=999.0;
    p_HadWH_BKG_MCFM_JECNominal=999.0;
    p_JJQCD_BKG_MCFM_JECNominal=999.0;
    p_HadZH_mavjj_JECNominal=999.0;
    p_HadZH_mavjj_true_JECNominal=999.0;
    p_HadWH_mavjj_JECNominal=999.0;
    p_HadWH_mavjj_true_JECNominal=999.0;
    pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal=999.0;
    pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal=999.0;
    pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal=999.0;
    pConst_JJVBF_BKG_MCFM_JECNominal=999.0;
    pConst_HadZH_BKG_MCFM_JECNominal=999.0;
    pConst_HadWH_BKG_MCFM_JECNominal=999.0;
    pConst_JJQCD_BKG_MCFM_JECNominal=999.0;

    p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal=999.0;
    pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal=999.0;
    p_HadWH_SIG_ghw1_1_JHUGen_JECNominal=999.0;
    pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal=999.0;
    p_HadZH_SIG_ghz1_1_JHUGen_JECNominal=999.0;
    pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal=999.0;

    D_HadWH=999.0; D_HadZH=999.0; 
    D_VBF=999.0; D_VBF1j=999.0; D_HadWH=999.0; D_HadZH=999.0;
    D_VBF_QG=999.0; D_VBF1j_QG=999.0; D_HadWH_QG=999.0; D_HadZH_QG=999.0;
    D_bkg_VBFdec=999.0; D_bkg_VHdec=999.0;

    if (verbose) {cout<<"clear other variables"<<endl; }
    // Resolution
    //massErrorUCSD=-1.0; massErrorUCSDCorr=-1.0; massErrorUF=-1.0; massErrorUFCorr=-1.0; massErrorUFADCorr=-1.0;

    // Event Category
    EventCat=-1;

    // Global variables not stored in tree
    lep_ptreco.clear(); lep_ptid.clear(); lep_ptindex.clear();
    recoMuons.clear(); recoElectrons.clear(); fsrPhotons.clear(); recoElectronsUnS.clear(); recoMuonsUnS.clear();
    HVec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    HVecNoFSR.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    Z1Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    Z2Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    GENZ1Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    GENZ2Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    RecoFourMuEvent = false; RecoFourEEvent = false;
    RecoTwoETwoMuEvent = false; RecoTwoMuTwoEEvent = false;
    foundHiggsCandidate = false; foundZ1LCandidate = false;
    jet1pt=-1.0; jet2pt=-1.0;

    // Float vectors
    lep_pt_float.clear(); lep_pterr_float.clear(); lep_pterrold_float.clear(); 
    lep_p_float.clear(); lep_ecalEnergy_float.clear();  
    lep_eta_float.clear(); lep_phi_float.clear(); lep_mass_float.clear();
    lepFSR_pt_float.clear(); lepFSR_eta_float.clear(); lepFSR_phi_float.clear(); lepFSR_mass_float.clear();
    tau_pt_float.clear(); tau_eta_float.clear(); tau_phi_float.clear(); tau_mass_float.clear();    
    pho_pt_float.clear(); pho_eta_float.clear(); pho_phi_float.clear(); photonCutBasedIDLoose_float.clear();
    H_pt_float.clear(); H_eta_float.clear(); H_phi_float.clear(); H_mass_float.clear();
    H_noFSR_pt_float.clear(); H_noFSR_eta_float.clear(); H_noFSR_phi_float.clear(); H_noFSR_mass_float.clear();
    Z_pt_float.clear(); Z_eta_float.clear(); Z_phi_float.clear(); Z_mass_float.clear();
    Z_noFSR_pt_float.clear(); Z_noFSR_eta_float.clear(); Z_noFSR_phi_float.clear(); Z_noFSR_mass_float.clear();
    jet_pt_float.clear(); jet_eta_float.clear(); jet_phi_float.clear(); jet_mass_float.clear(); jet_pt_raw_float.clear(); 
    jet_jesup_pt_float.clear(); jet_jesup_eta_float.clear(); jet_jesup_phi_float.clear(); jet_jesup_mass_float.clear();
    jet_jesdn_pt_float.clear(); jet_jesdn_eta_float.clear(); jet_jesdn_phi_float.clear(); jet_jesdn_mass_float.clear();
    jet_jerup_pt_float.clear(); jet_jerup_eta_float.clear(); jet_jerup_phi_float.clear(); jet_jerup_mass_float.clear();
    jet_jerdn_pt_float.clear(); jet_jerdn_eta_float.clear(); jet_jerdn_phi_float.clear();  jet_jerdn_mass_float.clear();
    fsrPhotons_pt_float.clear(); fsrPhotons_pterr_float.clear(); fsrPhotons_eta_float.clear(); fsrPhotons_phi_float.clear(); fsrPhotons_mass_float.clear();
	
	eleUncorr_pt_float.clear(); eleUncorr_mass_float.clear();
	muonUncorr_pt.clear(); muonUncorr_pterr.clear(); muon_scale_sf.clear();
	
	track_vtx_pt_float.clear(); track_vtx_eta_float.clear(); track_vtx_phi_float.clear(); track_vtx_pterr_float.clear();//track_vtx_ptcorr_float.clear(); track_vtx_pterrcorr_float.clear();//6
//	track_vtx_FSRpt_float.clear(); track_vtx_FSRpterr_float.clear(); track_vtx_FSRptcorr_float.clear(); track_vtx_FSRpterrcorr_float.clear(); track_vtx_FSReta_float.clear(); track_vtx_FSRphi_float.clear(); track_vtx_FSRmass_float.clear();//7
	track_vtx_dxy_float.clear(); track_vtx_dxycorr_float.clear();//4
		
	track_bs_pt_float.clear(); track_bs_eta_float.clear(); track_bs_phi_float.clear(); track_bs_pterr_float.clear();// track_bs_ptcorr_float.clear(); track_bs_pterrcorr_float.clear();//6
//	track_bs_FSRpt_float.clear(); track_bs_FSRpterr_float.clear(); track_bs_FSRptcorr_float.clear(); track_bs_FSRpterrcorr_float.clear(); track_bs_FSReta_float.clear(); track_bs_FSRphi_float.clear(); track_bs_FSRmass_float.clear();//7
	track_bs_dxy_float.clear(); track_bs_dxycorr_float.clear();//4
//	track_bs_xcorr_float.clear(); track_bs_ycorr_float.clear(); track_bs_zcorr_float.clear();//3
//	bs_x_float.clear(); bs_y_float.clear(); bs_z_float.clear();//3
	
	track_Hvtx_pt_float.clear(); track_Hvtx_eta_float.clear(); track_Hvtx_phi_float.clear(); track_Hvtx_pterr_float.clear();
	track_Hvtx_dxycorr_float.clear();
	track_Hbs_pt_float.clear(); track_Hbs_eta_float.clear(); track_Hbs_phi_float.clear(); track_Hbs_pterr_float.clear();
	track_Hbs_dxycorr_float.clear();
	
	
	
	edm::ESHandle<TransientTrackBuilder> trackbuilder;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackbuilder);
	
    // ====================== Do Analysis ======================== //

    std::map<int, TLorentzVector> fsrmap;
    vector<reco::Candidate*> selectedLeptons;
    std::map<unsigned int, TLorentzVector> selectedFsrMap;

    fsrmap.clear(); selectedFsrMap.clear(); selectedLeptons.clear();

    if (verbose) cout<<"start pileup reweighting"<<endl;
    // PU information
    if(isMC && reweightForPU) {        
        edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
        iEvent.getByToken(pileupSrc_, PupInfo);      

        if (verbose) cout<<"got pileup info"<<endl;

        std::vector<PileupSummaryInfo>::const_iterator PVI;      
        int npv = -1;
        for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
            int BX = PVI->getBunchCrossing();
            if(BX == 0) { npv = PVI->getTrueNumInteractions(); continue;}
        }        
        if (verbose) cout<<"N true interations = "<<npv<<endl;
        nInt = npv;
        //pileupWeight = pileUp.getPUWeight(npv,PUVersion);
        pileupWeight = pileUp.getPUWeight(h_pileup,npv);
        pileupWeightUp = pileUp.getPUWeight(h_pileupUp,npv);
        pileupWeightDn = pileUp.getPUWeight(h_pileupDn,npv);
        if (verbose) cout<<"pileup weight = "<<pileupWeight<<", filling histograms"<<endl;
        histContainer_["NINTERACT"]->Fill(npv);
        histContainer_["NINTERACT_RW"]->Fill(npv,pileupWeight);
    } else { pileupWeight = 1.0;}   

    if (verbose) {cout<<"finished pileup reweighting"<<endl; }
    
    if(isMC) {
        float tmpWeight = genEventInfo->weight();
        genWeight = (tmpWeight > 0 ? 1.0 : -1.0);
        if (verbose) {cout<<"tmpWeight: "<<tmpWeight<<"; genWeight: "<<genWeight<<endl;}
        double rms = 0.0;

        //std::cout<<"tmpWeight: "<<tmpWeight<<std::endl;

        if(lheInfo.isValid()){
            
            for(unsigned int i = 0; i < lheInfo->weights().size(); i++) {

                tmpWeight = genEventInfo->weight();
                tmpWeight *= lheInfo->weights()[i].wgt/lheInfo->originalXWGTUP();
                pdfWeights.push_back(tmpWeight);

                if (i<=8 or int(i)>=posNNPDF) {
                    tmpWeight = genEventInfo->weight();
                    tmpWeight *= lheInfo->weights()[i].wgt/lheInfo->originalXWGTUP();
                    if (int(i)<posNNPDF) {qcdWeights.push_back(tmpWeight);}
                }
                else {
                    tmpWeight = lheInfo->weights()[i].wgt;
                    tmpWeight /= lheInfo->originalXWGTUP();
                    //if (i==9) genWeight = tmpWeight;
                    if (int(i)<posNNPDF) {nnloWeights.push_back(tmpWeight);}
                }
                // NNPDF30 variations
                if (int(i)>=posNNPDF && int(i)<=(posNNPDF+100)) {
                    rms += tmpWeight*tmpWeight;
                    if (tmpWeight>pdfENVup) pdfENVup=tmpWeight;
                    if (tmpWeight<pdfENVdown) pdfENVdown=tmpWeight;
                }
            }
            pdfRMSup=sqrt(rms/100.0); pdfRMSdown=1.0/pdfRMSup;
            if (verbose) cout<<"pdfRMSup "<<pdfRMSup<<" pdfRMSdown "<<pdfRMSdown<<endl;
        
            const lhef::HEPEUP& lheEvent = lheInfo->hepeup();
            std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
            for ( size_t idxParticle = 0; idxParticle < lheParticles.size(); ++idxParticle ) {
                int id = std::abs(lheEvent.IDUP[idxParticle]);
                int status = lheEvent.ISTUP[idxParticle];
                if ( status == 1 && id==5 ) { 
                    lheNb += 1;
                }
                if ( status == 1 && ((id >= 1 && id <= 6) || id == 21) ) { 
                    lheNj += 1;
                }
            }
        
        }
        
        if (verbose) cout<<"setting gen variables"<<endl;       
        setGENVariables(prunedgenParticles,packedgenParticles,genJets); 
        if (verbose) { cout<<"finshed setting gen variables"<<endl;  }

        if (int(GENZ_pt.size()) == 2) {        
            GENZ1Vec.SetPtEtaPhiM(GENZ_pt[0], GENZ_eta[0], GENZ_phi[0], GENZ_mass[0]);
            GENZ2Vec.SetPtEtaPhiM(GENZ_pt[1], GENZ_eta[1], GENZ_phi[1], GENZ_mass[1]);
            k_qqZZ_ewk = getEwkCorrections(prunedgenParticles, tableEwk, genEventInfo, GENZ1Vec, GENZ2Vec);
            k_ggZZ = kFactor_ggzz->Eval((GENZ1Vec+GENZ2Vec).M());
            if (verbose) cout<<"ZZmass: "<< (GENZ1Vec+GENZ2Vec).M() <<"k factor qqZZ ewk: "<<k_qqZZ_ewk<<" ggZZ qcd"<<k_ggZZ << endl;

        }

    }
    sumWeightsTotal += genWeight;
    sumWeightsTotalPU += pileupWeight*genWeight;

    eventWeight = pileupWeight*genWeight;

    unsigned int _tSize = trigger->size();
    // create a string with all passing trigger names
    for (unsigned int i=0; i<_tSize; ++i) {
        std::string triggerName = trigNames.triggerName(i);
        if (strstr(triggerName.c_str(),"_step")) continue;
        if (strstr(triggerName.c_str(),"MC_")) continue;
        if (strstr(triggerName.c_str(),"AlCa_")) continue;
        if (strstr(triggerName.c_str(),"DST_")) continue;
        if (strstr(triggerName.c_str(),"HLT_HI")) continue;
        if (strstr(triggerName.c_str(),"HLT_Physics")) continue;
        if (strstr(triggerName.c_str(),"HLT_Random")) continue;
        if (strstr(triggerName.c_str(),"HLT_ZeroBias")) continue;
        if (strstr(triggerName.c_str(),"HLT_IsoTrack")) continue;
        if (strstr(triggerName.c_str(),"Hcal")) continue;
        if (strstr(triggerName.c_str(),"Ecal")) continue;
        if (trigger->accept(i)) triggersPassed += triggerName; 
    }
    if (firstEntry) cout<<"triggersPassed: "<<triggersPassed<<endl;
    firstEntry = false;
    // check if any of the triggers in the user list have passed
    bool passedSingleEl=false;
    bool passedSingleMu=false;
    bool passedAnyOther=false;
    for (unsigned int i=0; i<triggerList.size(); ++i) {
        if (strstr(triggersPassed.c_str(),triggerList.at(i).c_str())) {
            passedTrig=true;
            if (!isMC) {
                if (strstr(triggerList.at(i).c_str(),"_WP")) passedSingleEl=true;
                if (strstr(triggerList.at(i).c_str(),"HLT_Iso")) passedSingleMu=true;
                if (strstr(triggerList.at(i).c_str(),"CaloIdL")) passedAnyOther=true;
                if (strstr(triggerList.at(i).c_str(),"TrkIsoVVL")) passedAnyOther=true;
                if (strstr(triggerList.at(i).c_str(),"Triple")) passedAnyOther=true;
            }
        }
    }
    
    bool passedOnlySingle=((passedSingleEl && !passedAnyOther) || (passedSingleMu && !passedSingleEl && !passedAnyOther));
    //bool trigConditionData = ( passedTrig && (!checkOnlySingle || (checkOnlySingle && passedOnlySingle)) );
    bool trigConditionData = true;
    
    if (verbose) cout<<"checking PV"<<endl;       
    const reco::Vertex *PV = 0;
    int theVertex = -1;
    for (unsigned int i=0; i<vertex->size(); i++) {
        PV = &(vertex->at(i));        
        if (verbose) std::cout<<"isFake: "<<PV->isFake()<<" chi2 "<<PV->chi2()<<" ndof "<<PV->ndof()<<" rho "<<PV->position().Rho()<<" Z "<<PV->position().Z()<<endl; 
        //if (PV->chi2()==0 && PV->ndof()==0) continue;
        if (PV->isFake()) continue;
        if (PV->ndof()<=4 || PV->position().Rho()>2.0 || fabs(PV->position().Z())>24.0) continue;
        	pv_x = PV->x();
		pv_y = PV->y();
		pv_z = PV->z();
		theVertex=(int)i; break;//why do we only take into account one vertex?
    }        

    if (verbose) std::cout<<"vtx: "<<theVertex<<" trigConditionData "<<trigConditionData<<" passedTrig "<<passedTrig<<std::endl;
 
    if(theVertex >= 0 && (isMC || (!isMC && trigConditionData)) )  {

        if (verbose) cout<<"good PV "<<theVertex<<endl;
		bs_x = BS.x0();
		bs_y = BS.y0();
		bs_z = BS.z0();
		
        //N Vertex 
        if (verbose) {cout<<"fill nvtx histogram"<<endl;}
        nVtx = vertex->size();
        histContainer_["NVTX"]->Fill(nVtx);
        histContainer_["NVTX_RW"]->Fill(nVtx,pileupWeight);

        //MET
        if (verbose) {cout<<"get met value"<<endl;}
        if (!mets->empty()) {
            met = (*mets)[0].et();
            met_phi = (*mets)[0].phi();
            met_jesup = (*mets)[0].shiftedPt(pat::MET::JetEnUp);
            met_phi_jesup = (*mets)[0].shiftedPhi(pat::MET::JetEnUp);
            met_jesdn = (*mets)[0].shiftedPt(pat::MET::JetEnDown);
            met_phi_jesdn = (*mets)[0].shiftedPhi(pat::MET::JetEnDown);
            met_uncenup = (*mets)[0].shiftedPt(pat::MET::UnclusteredEnUp);
            met_phi_uncenup = (*mets)[0].shiftedPhi(pat::MET::UnclusteredEnUp);
            met_uncendn = (*mets)[0].shiftedPt(pat::MET::UnclusteredEnDown);
            met_phi_uncendn = (*mets)[0].shiftedPhi(pat::MET::UnclusteredEnDown);        
        }

        if (verbose) cout<<"start lepton analysis"<<endl;           
        vector<pat::Electron> AllElectrons; vector<pat::Muon> AllMuons; 
        vector<pat::Electron> AllElectronsUnS; vector<pat::Muon> AllMuonsUnS; ////uncorrected electron and muons 
        vector<pat::Tau> AllTaus; vector<pat::Photon> AllPhotons;
        AllElectrons = helper.goodLooseElectrons2012(electrons,_elecPtCut);
        AllElectronsUnS = helper.goodLooseElectrons2012(electrons,electronsUnS,_elecPtCut);
        AllMuons = helper.goodLooseMuons2012(muons,_muPtCut);
		AllMuonsUnS = helper.goodLooseMuons2012(muons,muonsUnS,_muPtCut);
        AllTaus = helper.goodLooseTaus2015(taus,_tauPtCut);
        AllPhotons = helper.goodLoosePhotons2015(photons,_phoPtCut);

        helper.cleanOverlappingLeptons(AllMuons,AllElectrons,PV);
        helper.cleanOverlappingLeptons(AllMuons,AllElectronsUnS,PV);
        recoMuons = helper.goodMuons2015_noIso_noPf(AllMuons,_muPtCut,PV,sip3dCut);
		recoMuonsUnS = helper.goodMuons2015_noIso_noPf(AllMuonsUnS,_muPtCut,PV,sip3dCut);
        recoElectrons = helper.goodElectrons2015_noIso_noBdt(AllElectrons,_elecPtCut,elecID,PV,iEvent,sip3dCut);
        recoElectronsUnS = helper.goodElectrons2015_noIso_noBdt(AllElectronsUnS,_elecPtCut,elecID,PV,iEvent,sip3dCut);
        helper.cleanOverlappingTaus(recoMuons,recoElectrons,AllTaus,isoCutMu,isoCutEl,muRho,elRho);
        recoTaus = helper.goodTaus2015(AllTaus,_tauPtCut);
        recoPhotons = helper.goodPhotons2015(AllPhotons,_phoPtCut,year);
        if (verbose) cout<<AllMuons.size()<<" loose muons "<<AllElectrons.size()<<" loose electrons"<<endl;

        //sort electrons and muons by pt
        if (verbose) cout<<recoMuons.size()<<" good muons and "<<recoElectrons.size()<<" good electrons to be sorted"<<endl;
        if (verbose) cout<<"start pt-sorting leptons"<<endl;
        if (verbose) cout<<"adding muons to sorted list"<<endl;

        if( (recoMuons.size() + recoElectrons.size()) >= (uint)skimLooseLeptons && (recoMuons.size()==recoMuonsUnS.size()) && (recoElectrons.size()==recoElectronsUnS.size()) ) {

            if (verbose) cout<<"found two leptons"<<endl;
//		cout<<Run<<" "<<Event<<" "<<LumiSect<<endl;
//	    cout<<"ele.size() "<<recoElectrons.size()<<" mu.size() "<<recoMuons.size()<<endl;
//	    cout<<"eleUnc.size() "<<recoElectronsUnS.size()<<" muUnc.size() "<<recoMuonsUnS.size()<<endl;
            
            for(unsigned int i = 0; i < recoMuons.size(); i++) {
                if (lep_ptreco.size()==0 || recoMuons[i].pt()<lep_ptreco[lep_ptreco.size()-1]) {
                    lep_ptreco.push_back(recoMuons[i].pt());
                    lep_ptid.push_back(recoMuons[i].pdgId());
                    lep_ptindex.push_back(i);
                    continue;
                }
                for (unsigned int j=0; j<lep_ptreco.size(); j++) {
                    if (recoMuons[i].pt()>lep_ptreco[j]) {
                        lep_ptreco.insert(lep_ptreco.begin()+j,recoMuons[i].pt());
                        lep_ptid.insert(lep_ptid.begin()+j,recoMuons[i].pdgId());
                        lep_ptindex.insert(lep_ptindex.begin()+j,i);
                        break;
                    }
                }
            }
            if (verbose) cout<<"adding electrons to sorted list"<<endl;           
            for(unsigned int i = 0; i < recoElectrons.size(); i++) {
                if (lep_ptreco.size()==0 || recoElectrons[i].pt()<lep_ptreco[lep_ptreco.size()-1]) {
                    lep_ptreco.push_back(recoElectrons[i].pt());
                    lep_ptid.push_back(recoElectrons[i].pdgId());
                    lep_ptindex.push_back(i);
                    continue;
                }
                for (unsigned int j=0; j<lep_ptreco.size(); j++) {
                    if (recoElectrons[i].pt()>lep_ptreco[j]) {
                        lep_ptreco.insert(lep_ptreco.begin()+j,recoElectrons[i].pt());
                        lep_ptid.insert(lep_ptid.begin()+j,recoElectrons[i].pdgId());
                        lep_ptindex.insert(lep_ptindex.begin()+j,i);
                        break;
                    }
                }
            }
            
            for(unsigned int i = 0; i < lep_ptreco.size(); i++) {
                if (verbose) cout<<"sorted lepton "<<i<<" pt "<<lep_ptreco[i]<<" id "<<lep_ptid[i]<<" index "<<lep_ptindex[i]<<endl;
                
                if (abs(lep_ptid[i])==11) {

					eleUncorr_pt.push_back(recoElectronsUnS[lep_ptindex[i]].pt());
					eleUncorr_mass.push_back(recoElectronsUnS[lep_ptindex[i]].mass());
					track_vtx_dxy.push_back(recoElectrons[lep_ptindex[i]].gsfTrack()->dxy(PV->position()));
					track_bs_dxy.push_back(recoElectrons[lep_ptindex[i]].gsfTrack()->dxy(BS.position()));
					
                    lep_id.push_back(recoElectrons[lep_ptindex[i]].pdgId());
                    lep_pt.push_back(recoElectrons[lep_ptindex[i]].pt());
                    lep_pterrold.push_back(recoElectrons[lep_ptindex[i]].p4Error(reco::GsfElectron::P4_COMBINATION));
                    lep_isEB.push_back(recoElectrons[lep_ptindex[i]].isEB());
                    lep_isEE.push_back(recoElectrons[lep_ptindex[i]].isEE());
                    lep_p.push_back(recoElectrons[lep_ptindex[i]].p());
                    lep_ecalEnergy.push_back(recoElectrons[lep_ptindex[i]].correctedEcalEnergy());
                    double perr = 0.0;
                    if (recoElectrons[lep_ptindex[i]].ecalDriven()) {
                        perr = recoElectrons[lep_ptindex[i]].p4Error(reco::GsfElectron::P4_COMBINATION);
                    }
                    else {
                        double ecalEnergy = recoElectrons[lep_ptindex[i]].correctedEcalEnergy();
                        double err2 = 0.0;
                        if (recoElectrons[lep_ptindex[i]].isEB()) {
                            err2 += (5.24e-02*5.24e-02)/ecalEnergy;
                            err2 += (2.01e-01*2.01e-01)/(ecalEnergy*ecalEnergy);
                            err2 += 1.00e-02*1.00e-02;
                        } else if (recoElectrons[lep_ptindex[i]].isEE()) {
                            err2 += (1.46e-01*1.46e-01)/ecalEnergy;
                            err2 += (9.21e-01*9.21e-01)/(ecalEnergy*ecalEnergy);
                            err2 += 1.94e-03*1.94e-03;
                        }
                        perr = ecalEnergy * sqrt(err2);
                    }
                    double pterr = perr*recoElectrons[lep_ptindex[i]].pt()/recoElectrons[lep_ptindex[i]].p();
                    lep_pterr.push_back(pterr);
                    lep_eta.push_back(recoElectrons[lep_ptindex[i]].eta());
                    lep_phi.push_back(recoElectrons[lep_ptindex[i]].phi());
                    lep_mass.push_back(recoElectrons[lep_ptindex[i]].mass());
                    lepFSR_pt.push_back(recoElectrons[lep_ptindex[i]].pt());
                    lepFSR_eta.push_back(recoElectrons[lep_ptindex[i]].eta());
                    lepFSR_phi.push_back(recoElectrons[lep_ptindex[i]].phi());
                    lepFSR_mass.push_back(recoElectrons[lep_ptindex[i]].mass());
                    if (isoConeSizeEl==0.4) {
                        lep_RelIso.push_back(helper.pfIso(recoElectrons[lep_ptindex[i]],elRho));
                        lep_RelIsoNoFSR.push_back(helper.pfIso(recoElectrons[lep_ptindex[i]],elRho));
                        lep_isoCH.push_back(recoElectrons[lep_ptindex[i]].chargedHadronIso());
                        lep_isoNH.push_back(recoElectrons[lep_ptindex[i]].neutralHadronIso());
                        lep_isoPhot.push_back(recoElectrons[lep_ptindex[i]].photonIso());
                        lep_isoPU.push_back(recoElectrons[lep_ptindex[i]].puChargedHadronIso());
                        lep_isoPUcorr.push_back(helper.getPUIso(recoElectrons[lep_ptindex[i]],elRho));
                    } else if (isoConeSizeEl==0.3) {
                        lep_RelIso.push_back(helper.pfIso03(recoElectrons[lep_ptindex[i]],elRho));
                        lep_RelIsoNoFSR.push_back(helper.pfIso03(recoElectrons[lep_ptindex[i]],elRho));
                        lep_isoCH.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumChargedHadronPt);
                        lep_isoNH.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumNeutralHadronEt);
                        lep_isoPhot.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumPhotonEt);
                        lep_isoPU.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumPUPt);
                        lep_isoPUcorr.push_back(helper.getPUIso03(recoElectrons[lep_ptindex[i]],elRho));
                    }
                    lep_MiniIso.push_back(helper.miniIsolation(pfCands, dynamic_cast<const reco::Candidate *>(&recoElectrons[lep_ptindex[i]]), 0.05, 0.2, 10., rhoSUS, false));
                    lep_Sip.push_back(helper.getSIP3D(recoElectrons[lep_ptindex[i]]));           
                    //lep_mva.push_back(recoElectrons[lep_ptindex[i]].userFloat("ElectronMVAEstimatorRun2Autumn18IdIsoValues")); 
                    //cout<<EleBDT_name_161718<<endl;
                    lep_mva.push_back(recoElectrons[lep_ptindex[i]].userFloat(EleBDT_name_161718.c_str())); 
                    lep_ecalDriven.push_back(recoElectrons[lep_ptindex[i]].ecalDriven()); 
                    //lep_tightId.push_back(helper.passTight_BDT_Id(recoElectrons[lep_ptindex[i]],recoElectrons[lep_ptindex[i]].userFloat("ElectronMVAEstimatorRun2Autumn18IdIsoValues"), year));           
                    lep_tightId.push_back(helper.passTight_BDT_Id(recoElectronsUnS[lep_ptindex[i]],year));
                    //lep_tightId_old.push_back(helper.passTight_BDT_Id(recoElectronsUnS[lep_ptindex[i]],year));
                    //lep_tightId_old.push_back(helper.passTight_BDT_Id(0);
                    //cout<<"old "<<recoElectrons[lep_ptindex[i]].userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") <<" new" <<recoElectrons[lep_ptindex[i]].userFloat("ElectronMVAEstimatorRun2Spring16HZZV1Values")<<endl;
                    lep_tightIdSUS.push_back(helper.passTight_Id_SUS(recoElectrons[lep_ptindex[i]],elecID,PV,BS,theConversions, year));           
                    //lep_tightIdHiPt.push_back(recoElectrons[lep_ptindex[i]].electronID("heepElectronID-HEEPV70"));
                    lep_tightIdHiPt.push_back(recoElectrons[lep_ptindex[i]].electronID(heepID_name_161718.c_str()));
                    lep_ptRatio.push_back(helper.ptRatio(recoElectrons[lep_ptindex[i]],jets,true)); // no L2L3 yet           
                    lep_ptRel.push_back(helper.ptRel(recoElectrons[lep_ptindex[i]],jets,true)); // no L2L3 yet           
                    lep_dataMC.push_back(helper.dataMC(recoElectrons[lep_ptindex[i]],hElecScaleFac,hElecScaleFac_Cracks,hElecScaleFacGsf,hElecScaleFacGsfLowET));
                    lep_dataMCErr.push_back(helper.dataMCErr(recoElectrons[lep_ptindex[i]],hElecScaleFac,hElecScaleFac_Cracks));
                    lep_genindex.push_back(-1.0);
                }
                if (abs(lep_ptid[i])==13) { 
		//		cout<<"uncorr eta: "<<recoMuonsUnS[lep_ptindex[i]].eta()<<" corr eta: "<<recoMuons[lep_ptindex[i]].eta()<<endl;	
					muonUncorr_pt.push_back(recoMuonsUnS[lep_ptindex[i]].pt());
					if(recoMuons[lep_ptindex[i]].hasUserFloat("scale_sf")){
						muon_scale_sf.push_back(recoMuons[lep_ptindex[i]].userFloat("scale_sf"));
					} else{
						muon_scale_sf.push_back(-999);
					}
					if(recoMuons[lep_ptindex[i]].hasUserFloat("mu_pt_roc")){
						mu_pt_roc.push_back(recoMuons[lep_ptindex[i]].userFloat("mu_pt_roc"));//uncorr pt 
					} else{
						mu_pt_roc.push_back(-999);
					}
					if(recoMuons[lep_ptindex[i]].hasUserFloat("mu_eta_roc")){
						mu_eta_roc.push_back(recoMuons[lep_ptindex[i]].userFloat("mu_eta_roc"));
					} else{
						mu_eta_roc.push_back(-999);
					}
					if(recoMuons[lep_ptindex[i]].hasUserFloat("mu_phi_roc")){
						mu_phi_roc.push_back(recoMuons[lep_ptindex[i]].userFloat("mu_phi_roc"));
					} else{
						mu_phi_roc.push_back(-999);
					}
					if(recoMuons[lep_ptindex[i]].hasUserFloat("mu_q_roc")){
						mu_q_roc.push_back(recoMuons[lep_ptindex[i]].userFloat("mu_q_roc"));
					} else{
						mu_q_roc.push_back(-999);
					}
					if(recoMuons[lep_ptindex[i]].hasUserFloat("mu_tracktype_roc")){
						mu_tracktype_roc.push_back(recoMuons[lep_ptindex[i]].userFloat("mu_tracktype_roc"));
					} else{
						mu_tracktype_roc.push_back(-999);
					}
					if(recoMuons[lep_ptindex[i]].hasUserFloat("mu_nl_roc")){
						mu_nl_roc.push_back(recoMuons[lep_ptindex[i]].userFloat("mu_nl_roc"));
					} else{
						mu_nl_roc.push_back(-999);
					}
					if(recoMuons[lep_ptindex[i]].hasUserFloat("mu_genparticle_roc")){
						mu_genparticle_roc.push_back(recoMuons[lep_ptindex[i]].userFloat("mu_genparticle_roc"));
					} else{
						mu_genparticle_roc.push_back(-999);
					}
					if(recoMuons[lep_ptindex[i]].hasUserFloat("mu_gen_pt_roc")){
						mu_gen_pt_roc.push_back(recoMuons[lep_ptindex[i]].userFloat("mu_gen_pt_roc"));
					} else{
						mu_gen_pt_roc.push_back(-999);
					}
                    lep_id.push_back(recoMuons[lep_ptindex[i]].pdgId());
                    lep_pt.push_back(recoMuons[lep_ptindex[i]].pt());
                    lep_pterrold.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->ptError());
                    lep_isEB.push_back(0);
                    lep_isEE.push_back(0);
                    lep_p.push_back(recoMuons[lep_ptindex[i]].p());
                    lep_ecalEnergy.push_back(0);

                    if (recoMuons[lep_ptindex[i]].hasUserFloat("correctedPtError")) {
                        lep_pterr.push_back(recoMuons[lep_ptindex[i]].userFloat("correctedPtError"));
						muon_scale_sf.push_back(recoMuons[lep_ptindex[i]].userFloat("scale_sf"));
                    } else {
                        lep_pterr.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->ptError());
                    }
                    lep_eta.push_back(recoMuons[lep_ptindex[i]].eta());
                    lep_phi.push_back(recoMuons[lep_ptindex[i]].phi());
                    if (recoMuons[lep_ptindex[i]].mass()<0.105) cout<<"muon mass: "<<recoMuons[lep_ptindex[i]].mass()<<endl;
                    lep_mass.push_back(recoMuons[lep_ptindex[i]].mass());
                    lepFSR_pt.push_back(recoMuons[lep_ptindex[i]].pt());
                    lepFSR_eta.push_back(recoMuons[lep_ptindex[i]].eta());
                    lepFSR_phi.push_back(recoMuons[lep_ptindex[i]].phi());
                    lepFSR_mass.push_back(recoMuons[lep_ptindex[i]].mass());
                    if (isoConeSizeMu==0.4) {
                        lep_RelIso.push_back(helper.pfIso(recoMuons[lep_ptindex[i]],muRho));
                        lep_RelIsoNoFSR.push_back(helper.pfIso(recoMuons[lep_ptindex[i]],muRho));
                        lep_isoCH.push_back(recoMuons[lep_ptindex[i]].chargedHadronIso());
                        lep_isoNH.push_back(recoMuons[lep_ptindex[i]].neutralHadronIso());
                        lep_isoPhot.push_back(recoMuons[lep_ptindex[i]].photonIso());
                        lep_isoPU.push_back(recoMuons[lep_ptindex[i]].puChargedHadronIso());
                        lep_isoPUcorr.push_back(helper.getPUIso(recoMuons[lep_ptindex[i]],muRho));
                    } else if (isoConeSizeMu==0.3) {
                        lep_RelIso.push_back(helper.pfIso03(recoMuons[lep_ptindex[i]],muRho));
                        lep_RelIsoNoFSR.push_back(helper.pfIso03(recoMuons[lep_ptindex[i]],muRho));
                        lep_isoCH.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumChargedHadronPt);
                        lep_isoNH.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumNeutralHadronEt);
                        lep_isoPhot.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumPhotonEt);
                        lep_isoPU.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumPUPt);
                        lep_isoPUcorr.push_back(helper.getPUIso03(recoMuons[lep_ptindex[i]],muRho));
                    }
                    lep_MiniIso.push_back(helper.miniIsolation(pfCands, dynamic_cast<const reco::Candidate *>(&recoMuons[lep_ptindex[i]]), 0.05, 0.2, 10., rhoSUS, false));
                    lep_Sip.push_back(helper.getSIP3D(recoMuons[lep_ptindex[i]]));            
                    lep_mva.push_back(recoMuons[lep_ptindex[i]].isPFMuon());  
                    //lep_mva.push_back(helper.get_Muon_MVA_Value(recoMuons[lep_ptindex[i]],vertex,muRho,year,PV));
                    lep_ecalDriven.push_back(0);  
                    lep_tightId.push_back(helper.passTight_Id(recoMuons[lep_ptindex[i]],PV));         
                    //lep_tightId.push_back(helper.passTight_BDT_Id(recoMuons[lep_ptindex[i]],vertex,muRho,year,PV));
                    //lep_tightId_old.push_back(helper.passTight_Id(recoMuons[lep_ptindex[i]],PV));
                    //lep_tightId_old.push_back(helper.passTight_Id(0);
                    lep_tightIdSUS.push_back(helper.passTight_Id_SUS(recoMuons[lep_ptindex[i]],PV));
                    lep_tightIdHiPt.push_back(recoMuons[lep_ptindex[i]].isHighPtMuon(*PV));
                    lep_ptRatio.push_back(helper.ptRatio(recoMuons[lep_ptindex[i]],jets,true)); // no L2L3 yet           
                    lep_ptRel.push_back(helper.ptRel(recoMuons[lep_ptindex[i]],jets,true)); // no L2L3 yet       
                    lep_dataMC.push_back(helper.dataMC(recoMuons[lep_ptindex[i]],hMuScaleFac));
                    lep_dataMCErr.push_back(helper.dataMCErr(recoMuons[lep_ptindex[i]],hMuScaleFacUnc));
                    lep_genindex.push_back(-1.0);
             	    //muonrefittedtracks.push_back(trackbuilder->build(recoMuons[lep_ptindex[i]].muonBestTrack()));
					//muontracks.push_back(recoMuons[lep_ptindex[i]].muonBestTrack());
					track_vtx_dxy.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->dxy(PV->position()));
					track_bs_dxy.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->dxy(BS.position()));
						
                }
                if (verbose) {cout<<" eta: "<<lep_eta[i]<<" phi: "<<lep_phi[i];
                              if(abs(lep_ptid[i])==11)  cout<<" eSuperClusterOverP: "<<recoElectrons[lep_ptindex[i]].eSuperClusterOverP()<<" ecalEnergy: "<<recoElectrons[lep_ptindex[i]].ecalEnergy()<<" p: "<<recoElectrons[lep_ptindex[i]].p();
                              cout<<" RelIso: "<<lep_RelIso[i]<<" isoCH: "<<lep_isoCH[i]<<" isoNH: "<<lep_isoNH[i]
                                  <<" isoPhot: "<<lep_isoPhot[i]<<" lep_isoPU: "<<lep_isoPU[i]<<" isoPUcorr: "<<lep_isoPUcorr[i]<<" Sip: "<<lep_Sip[i]
                                  <<" MiniIso: "<<lep_MiniIso[i]<<" ptRatio: "<<lep_ptRatio[i]<<" ptRel: "<<lep_ptRel[i]<<" lep_mva: "<<lep_mva[i];
                              if(abs(lep_ptid[i])==11)    cout<<" SCeta: "<<recoElectrons[lep_ptindex[i]].superCluster()->eta()<<" dxy: "<<recoElectrons[lep_ptindex[i]].gsfTrack()->dxy(PV->position())<<" dz: "<<recoElectrons[lep_ptindex[i]].gsfTrack()->dz(PV->position());
                              if(abs(lep_ptid[i])==11)    cout<<" Rho: "<<elRho<<" EleBDT_name: "<<EleBDT_name_161718<<" Uncorrected electron pt: "<<recoElectronsUnS[lep_ptindex[i]].pt();
                              if(abs(lep_ptid[i])==13)    cout<<" Rho: "<<muRho;
                              cout<<" dataMC: "<<lep_dataMC[i]<<" dataMCErr: "<<lep_dataMCErr[i];
                              cout<<" lep_pterr: "<<lep_pterr[i]<<" lep_pterrold: "<<lep_pterrold[i]<<" lep_tightIdHiPt: "<<lep_tightIdHiPt[i]<<endl;
                              if((abs(lep_ptid[i])==13)&&lep_pt[i]>200)    cout<<"Muon pt over 200 isTrackerHighPtID? "<<helper.isTrackerHighPt(recoMuons[lep_ptindex[i]],PV)<<endl;}
            
            }
			if(Zvtx){//vertex constraint for Z events
				for(unsigned int i = 0; i < lep_ptreco.size(); i++) {
					if(lep_ptreco.size()<2)break;
					if(abs(lep_ptid[i])==11){
						if(!helper.passTight_BDT_Id(recoElectronsUnS[lep_ptindex[i]],year))continue;
						if(isoConeSizeEl==0.4&&helper.pfIso(recoElectrons[lep_ptindex[i]],elRho)>0.35)continue;
						if(isoConeSizeEl==0.3&&helper.pfIso03(recoElectrons[lep_ptindex[i]],elRho)>0.35)continue;
						tracks_tmp.push_back(trackbuilder->build(recoElectrons[lep_ptindex[i]].gsfTrack()));
						passLepindex.push_back(i);
					}
					if(abs(lep_ptid[i])==13){
						if(!helper.passTight_Id(recoMuons[lep_ptindex[i]],PV))continue;
						if(isoConeSizeMu==0.4&&helper.pfIso(recoMuons[lep_ptindex[i]],muRho)>0.35)continue;
						if(isoConeSizeMu==0.3&&helper.pfIso03(recoMuons[lep_ptindex[i]],muRho)>0.35)continue;
						tracks_tmp.push_back(trackbuilder->build(recoMuons[lep_ptindex[i]].muonBestTrack()));
						passLepindex.push_back(i);
					}
				}
				if(passLepindex.size()==2&&tracks_tmp.size()==2&&(lep_id[passLepindex[0]]+lep_id[passLepindex[1]])==0){
					KalmanVertexFitter fitter(true);
					TransientVertex Zmumuvertex = fitter.vertex(tracks_tmp);
					if(Zmumuvertex.isValid()&&Zmumuvertex.hasRefittedTracks()&&Zmumuvertex.refittedTracks().size()==2){
						track_vtx_xcorr = Zmumuvertex.position().x();
						track_vtx_ycorr = Zmumuvertex.position().y();
						track_vtx_zcorr = Zmumuvertex.position().z();
						GlobalPoint vert(pv_x, pv_y, pv_z);
						for(unsigned int k = 0; k < Zmumuvertex.refittedTracks().size(); k++){
							TrajectoryStateClosestToPoint traj = Zmumuvertex.refittedTracks().at(k).trajectoryStateClosestToPoint(vert);
							track_vtx_dxycorr.push_back(traj.perigeeParameters().transverseImpactParameter());
							track_vtx_pt.push_back(Zmumuvertex.refittedTracks().at(k).track().pt());
							track_vtx_eta.push_back(Zmumuvertex.refittedTracks().at(k).track().eta());
							track_vtx_phi.push_back(Zmumuvertex.refittedTracks().at(k).track().phi());
							track_vtx_pterr.push_back(Zmumuvertex.refittedTracks().at(k).track().ptError());
						}
					}
				
					TransientVertex ZmumubeamSpot = fitter.vertex(tracks_tmp, BS);
					if(ZmumubeamSpot.isValid()&&ZmumubeamSpot.hasRefittedTracks()&&ZmumubeamSpot.refittedTracks().size()==2){
						track_bs_xcorr = ZmumubeamSpot.position().x();
						track_bs_ycorr = ZmumubeamSpot.position().y();
						track_bs_zcorr = ZmumubeamSpot.position().z();
						GlobalPoint vert(pv_x, pv_y, pv_z);
						for(unsigned int k = 0; k < ZmumubeamSpot.refittedTracks().size(); k++){
							TrajectoryStateClosestToPoint traj = ZmumubeamSpot.refittedTracks().at(k).trajectoryStateClosestToPoint(vert);
							track_bs_dxycorr.push_back(traj.perigeeParameters().transverseImpactParameter());
							track_bs_pt.push_back(ZmumubeamSpot.refittedTracks().at(k).track().pt());
							track_bs_eta.push_back(ZmumubeamSpot.refittedTracks().at(k).track().eta());
							track_bs_phi.push_back(ZmumubeamSpot.refittedTracks().at(k).track().phi());
							track_bs_pterr.push_back(ZmumubeamSpot.refittedTracks().at(k).track().ptError());
						}
					}

					if(Checkvtx&&Zmumuvertex.isValid()&&Zmumuvertex.hasRefittedTracks()&&Zmumuvertex.refittedTracks().size()==2){
						std::cout<<"!!!!!!!!!!!!!!!!!!"<<std::endl;
						std::cout<<"Zevent oldpt:"<<lep_pt[passLepindex[1]]<<std::endl;
						std::cout<<"Zevent newpt:"<<track_vtx_pt[1]<<std::endl;
					}
				}
			}
/*			//set rochester correction 
			for(unsigned int c = 0; c < 2; c++){
				if(c==0&&track_vtx_pt.size()!=2||passLepindex.size()!=2)break;
				if(c==1&&track_bs_pt.size()!=2||passLepindex.size()!=2)break;
				std::vector<double> lep_tmp_pt; std::vector<double> lep_tmp_eta; std::vector<double> lep_tmp_phi; std::vector<double> lep_tmp_pterr; 
				lep_tmp_pt.clear(); lep_tmp_eta.clear(); lep_tmp_phi.clear(); lep_tmp_pterr.clear();
				
				if(c==0){//for common vtx constraint 
					lep_tmp_pt.push_back(track_vtx_pt[0]); lep_tmp_pt.push_back(track_vtx_pt[1]);
					lep_tmp_eta.push_back(track_vtx_eta[0]); lep_tmp_eta.push_back(track_vtx_eta[1]);
					lep_tmp_phi.push_back(track_vtx_phi[0]); lep_tmp_phi.push_back(track_vtx_phi[1]);
					lep_tmp_pterr.push_back(track_vtx_pterr[0]); lep_tmp_pterr.push_back(track_vtx_pterr[1]);
				}
				if(c==1){//for beamspot constraint
					lep_tmp_pt.push_back(track_bs_pt[0]); lep_tmp_pt.push_back(track_bs_pt[1]);
					lep_tmp_eta.push_back(track_bs_eta[0]); lep_tmp_eta.push_back(track_bs_eta[1]);
					lep_tmp_phi.push_back(track_bs_phi[0]); lep_tmp_phi.push_back(track_bs_phi[1]);
					lep_tmp_pterr.push_back(track_bs_pterr[0]); lep_tmp_pterr.push_back(track_bs_pterr[1]);
				}
				for(unsigned int i = 0; i < 2; i++){
					TRandom3 rand;                                                                                                                                                           			            rand.SetSeed(abs(static_cast<int>(sin(lep_tmp_phi[i])*100000)));
					double u1 = rand.Uniform(1.0);
					double oldpt = lep_tmp_pt[i];
					double oldpterr = lep_tmp_pterr[i];
					double newpt = oldpt;
					double newpterr = oldpterr;
					double scale_factor = 1.0;
//					double smear_error = 0.0;
//					double scale_error = 0.0;
					if(year == 2018 || year == 2017 || year == 2016){//use RoccoR
						if(recoMuons[lep_ptindex[passLepindex[i]]].muonBestTrackType() == 1 && oldpt < 200){
							int nl = recoMuons[lep_ptindex[passLepindex[i]]].muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
							if(isMC && nl > 5){
								auto gen_particle = recoMuons[lep_ptindex[passLepindex[i]]].genParticle();
								if(gen_particle != 0 ){
									scale_factor = calibrator->kSpreadMC(recoMuons[lep_ptindex[passLepindex[i]]].charge(), oldpt, lep_tmp_eta[i], lep_tmp_phi[i], gen_particle->pt());
//									smear_error = calibrator->kSpreadMCerror(recoMuons[lep_ptindex[passLepindex[i]]].charge(), oldpt, lep_tmp_eta[i], lep_tmp_phi[i], gen_particle->pt());
								}
								else{
									scale_factor = calibrator->kSmearMC(recoMuons[lep_ptindex[passLepindex[i]]].charge(), oldpt, lep_tmp_eta[i], lep_tmp_phi[i], nl, u1);
//									smear_error = calibrator->kSmearMCerror(recoMuons[lep_ptindex[passLepindex[i]]].charge(), oldpt, lep_tmp_eta[i], lep_tmp_phi[i], nl, u1);
								}
//								scale_error = calibrator->kScaleDTerror(recoMuons[lep_ptindex[passLepindex[i]]].charge(), oldpt, lep_tmp_eta[i], lep_tmp_phi[i]);
								newpt = scale_factor*oldpt;
								newpterr = scale_factor*oldpterr;
							}
							else if(!isMC && nl > 5){
								if(oldpt>2.0 && fabs(lep_tmp_eta[i])<2.4){
									scale_factor = calibrator->kScaleDT(recoMuons[lep_ptindex[passLepindex[i]]].charge(), oldpt, lep_tmp_eta[i], lep_tmp_phi[i]);
//									scale_error = calibrator->kScaleDTerror(recoMuons[lep_ptindex[passLepindex[i]]].charge(), oldpt, lep_tmp_eta[i], lep_tmp_phi[i]);
//									smear_error = calibrator->kSmearMCerror(recoMuons[lep_ptindex[passLepindex[i]]].charge(), oldpt, lep_tmp_eta[i], lep_tmp_phi[i], nl, u1);
								}
								else{
									scale_factor = 1.0;
//									scale_error = 0.0;
//									smear_error = 0.0;
								}
								newpt = scale_factor*oldpt;
								newpterr = scale_factor*oldpt;
							}
						}
					}//use RoccoR
//					else if(year == 2016){//use kalmanCalibrator
//						if(recoMuons[lep_ptindex[passLepindex[i]]].muonBestTrackType() == 1 && oldpt < 200){
//							if(!isMC){
//								if(oldpt>2.0 && fabs(lep_tmp_eta[i])<2.4){
//									newpt = kalmanMuonCalibrator->getCorrectedPt(oldpt, lep_tmp_eta[i], lep_tmp_phi[i], recoMuons[lep_ptindex[passLepindex[i]]].charge());
//									newpterr = newpt*kalmanMuonCalibrator->getCorrectedError(newpt, lep_tmp_eta[i], oldpterr/newpt);
//								}
//							}
//							else{
//								double unsmearednewpt = kalmanMuonCalibrator->getCorrectedPt(oldpt, lep_tmp_eta[i], lep_tmp_phi[i], recoMuons[lep_ptindex[passLepindex[i]]].charge());
//								newpt = kalmanMuonCalibrator->smearForSync(unsmearednewpt, lep_tmp_eta[i]);
//								newpterr = newpt*kalmanMuonCalibrator->getCorrectedError(newpt, lep_tmp_eta[i], oldpterr/newpt );
//							}//isMC
//						}
//					}//use kalmanCalibrator
					
					if(c==0){
						track_vtx_ptcorr.push_back(newpt);
						track_vtx_pterrcorr.push_back(newpterr);
					}else{
						track_bs_ptcorr.push_back(newpt);
						track_bs_pterrcorr.push_back(newpterr);
					}			
				}//2 leptons
				
			}// vtx and bs
						
*/				
/*			
			//set FSR for Zmumu vtx constrainted and Zmumu beamspot constraint events
			for(unsigned int c = 0; c < 2; c++){
				if(track_vtx_pt.size()!=2||passLepindex.size()!=2)break;
				std::vector<double> lep_tmp_pt; std::vector<double> lep_tmp_eta; std::vector<double> lep_tmp_phi;
				lep_tmp_pt.clear(); lep_tmp_eta.clear(); lep_tmp_phi.clear();
				
				if(c==0){//for common vtx constraint 
					lep_tmp_pt.push_back(track_vtx_ptcorr[0]); lep_tmp_pt.push_back(track_vtx_ptcorr[1]);
					lep_tmp_eta.push_back(track_vtx_eta[0]); lep_tmp_eta.push_back(track_vtx_eta[1]);
					lep_tmp_phi.push_back(track_vtx_phi[0]); lep_tmp_phi.push_back(track_vtx_phi[1]);
				}
				if(c==1){//for beamspot constraint
					lep_tmp_pt.push_back(track_bs_ptcorr[0]); lep_tmp_pt.push_back(track_bs_ptcorr[1]);
					lep_tmp_eta.push_back(track_bs_eta[0]); lep_tmp_eta.push_back(track_bs_eta[1]);
					lep_tmp_phi.push_back(track_bs_phi[0]); lep_tmp_phi.push_back(track_bs_phi[1]);
				}
				
				for(unsigned int i = 0; i < 2; i++){
					TLorentzVector thisLep;
					thisLep.SetPtEtaPhiM(lep_tmp_pt[i], lep_tmp_eta[i], lep_tmp_phi[i], 0.1057);
					
					double minDrOEt2 = 999.0; double selectedPhotonIso=9999.0; double selectedPhotonDr=9999.0;
					bool selected=false; pat::PFParticle selectedPhoton; 
				
					for(edm::View<pat::PFParticle>::const_iterator phot=photonsForFsr->begin(); phot!=photonsForFsr->end(); ++phot){
						if(fabs(phot->eta()) > 2.4)continue;
						if(phot->pt() < 2.0)continue;
						double fsrDr = deltaR(thisLep.Eta(), thisLep.Phi(), phot->eta(), phot->phi());
						if (fsrDr>0.5) continue;
						bool closest=true;
						for (unsigned int j=0; j<2; j++) {
							TLorentzVector otherLep;
							otherLep.SetPtEtaPhiM(lep_tmp_pt[j],lep_tmp_eta[j],lep_tmp_phi[j],0.1057);
							double fsrDrOther = deltaR(otherLep.Eta(), otherLep.Phi(), phot->eta(), phot->phi());
							if (j!=i && fsrDrOther<fsrDr) {closest=false;}
						}
						if(!closest)continue;
						double photoniso = helper.photonPfIso03(*phot,pfCands)/phot->pt();
						double fsrDrOEt2 = fsrDr/(phot->pt()*phot->pt());
						if (photoniso>1.8) continue;
						if (fsrDrOEt2>0.012) continue;
						if( fsrDrOEt2 < minDrOEt2 ){
							selected = true;
							selectedPhoton=(*phot);
						}
					}
					if(selected){
						TLorentzVector vtx_fsr;
						vtx_fsr.SetPtEtaPhiM(selectedPhoton.pt(),selectedPhoton.eta(),selectedPhoton.phi(),0);
						TLorentzVector vtx_lepfsr;
						vtx_lepfsr = thisLep + vtx_lepfsr;
						if(c==0){
							track_vtx_FSRpt.push_back(vtx_lepfsr.Pt());
							track_vtx_FSReta.push_back(vtx_lepfsr.Eta());
							track_vtx_FSRphi.push_back(vtx_lepfsr.Phi());
							track_vtx_FSRmass.push_back(vtx_lepfsr.M());
						}
						else{
							track_bs_FSRpt.push_back(vtx_lepfsr.Pt());
							track_bs_FSReta.push_back(vtx_lepfsr.Eta());
							track_bs_FSRphi.push_back(vtx_lepfsr.Phi());
							track_bs_FSRmass.push_back(vtx_lepfsr.M());
						}
					}
					else{
						if(c==0){
							track_vtx_FSRpt.push_back(track_vtx_ptcorr[i]);
							track_vtx_FSReta.push_back(track_vtx_eta[i]);
							track_vtx_FSRphi.push_back(track_vtx_phi[i]);
							track_vtx_FSRmass.push_back(0.1057);
						}
						else{
							track_bs_FSRpt.push_back(track_vtx_ptcorr[i]);
							track_bs_FSReta.push_back(track_vtx_eta[i]);
							track_bs_FSRphi.push_back(track_vtx_phi[i]);
							track_bs_FSRmass.push_back(0.1057);
						}	
					}//select
				}//all lepton
			}//vtx and bs
*/           
			
			
            if (verbose) cout<<"adding taus to sorted list"<<endl;           
            for(int i = 0; i < (int)recoTaus.size(); i++) {
                tau_id.push_back(recoTaus[i].pdgId());
                tau_pt.push_back(recoTaus[i].pt());
                tau_eta.push_back(recoTaus[i].eta());
                tau_phi.push_back(recoTaus[i].phi());            
                tau_mass.push_back(recoTaus[i].mass());
            }
            
            if (verbose) cout<<"adding photons to sorted list"<<endl;           
            for(int i = 0; i < (int)recoPhotons.size(); i++) {
                pho_pt.push_back(recoPhotons[i].pt());
                pho_eta.push_back(recoPhotons[i].eta());
                pho_phi.push_back(recoPhotons[i].phi());            
                photonCutBasedIDLoose.push_back(recoPhotons[i].photonID("cutBasedPhotonID-Fall17-94X-V2-loose"));
            }
            
            if (doTriggerMatching) {
                if (verbose) cout<<"start trigger matching"<<endl;        
                // trigger Matching
                for(unsigned int i = 0; i < lep_pt.size(); i++) {
                    
                    TLorentzVector reco;
                    reco.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
                    
                    double reco_eta = reco.Eta();
                    double reco_phi = reco.Phi();
                    
                    std::string filtersMatched = "";
                    
                    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                        double hlt_eta = obj.eta();
                        double hlt_phi = obj.phi();
                        double dR =  deltaR(reco_eta,reco_phi,hlt_eta,hlt_phi); 
                        if (dR<0.5) {
                            obj.unpackFilterLabels(iEvent, *trigger); 
                            for (unsigned h = 0; h < obj.filterLabels().size(); ++h) filtersMatched += obj.filterLabels()[h];
                        }
                    }
                    
                    if (verbose) cout<<"Trigger matching lep id: "<<lep_id[i]<<" pt: "<<reco.Pt()<<" filters: "<<filtersMatched<<endl;
                    lep_filtersMatched.push_back(filtersMatched);
                    
                }
            }
            
            // GEN matching
            if(isMC) {
                if (verbose) cout<<"begin gen matching"<<endl;
                // for each reco lepton find the nearest gen lepton with same ID
                for(unsigned int i = 0; i < lep_pt.size(); i++) {
                    
                    double minDr=9999.0;
                    
                    TLorentzVector reco, gen;
                    reco.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
                    
                    for (unsigned int j = 0; j < GENlep_id.size(); j++) {
                        
                        if (GENlep_id[j]!=lep_id[i]) continue;
                        
                        gen.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);
                        double thisDr = deltaR(reco.Eta(),reco.Phi(),gen.Eta(),gen.Phi());
                        
                        if (thisDr<minDr && thisDr<0.5) {
                            lep_genindex[i]=j;
                            minDr=thisDr;
                        }
                        
                    } // all gen leptons
                    
                } // all reco leptons
                
                // for each reco lepton find the nearest gen particle and save its ID and mothers 
                for(unsigned int i = 0; i < lep_pt.size(); i++) {
                    
                    double minDr=9999.0;
                    
                    TLorentzVector reco;
                    reco.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
                    
                    reco::GenParticleCollection::const_iterator genPart;
                    int j = -1;
                    int tmpPdgId = 0;
                    int tmpMomId = 0;
                    int tmpMomMomId = 0;
                    for(genPart = prunedgenParticles->begin(); genPart != prunedgenParticles->end(); genPart++) {
                        j++;
                        double thisDr = deltaR(reco.Eta(),reco.Phi(),genPart->eta(),genPart->phi());
                        
                        if (thisDr<minDr && thisDr<0.3) {
                            tmpPdgId=genPart->pdgId();
                            tmpMomId=genAna.MotherID(&prunedgenParticles->at(j));
                            tmpMomMomId=genAna.MotherMotherID(&prunedgenParticles->at(j));
                            
                            minDr=thisDr;
                        }
                        
                    } // all gen particles
                    // storing the matches                                                                                                            
                    lep_matchedR03_PdgId.push_back(tmpPdgId);
                    lep_matchedR03_MomId.push_back(tmpMomId);
                    lep_matchedR03_MomMomId.push_back(tmpMomMomId);
                    
                } // all reco leptons
                
                if (verbose) {cout<<"finished gen matching"<<endl;}
            } //isMC
            
            unsigned int Nleptons = lep_pt.size();
            
            // FSR Photons
            if(doFsrRecovery) {
                
                if (verbose) cout<<"checking "<<photonsForFsr->size()<<" fsr photon candidates"<<endl;
                
                // try to find an fsr photon for each lepton
                for (unsigned int i=0; i<Nleptons; i++) {
                        
                    TLorentzVector thisLep;
                    thisLep.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

                    double minDrOEt2 = 999.0; double selectedPhotonIso=9999.0; double selectedPhotonDr=9999.0;
                    bool selected=false; pat::PFParticle selectedPhoton; 
                    for(edm::View<pat::PFParticle>::const_iterator phot=photonsForFsr->begin(); phot!=photonsForFsr->end(); ++phot) {

                        // preselection
                        if (fabs(phot->eta()) > 2.4) continue;
                        if (phot->pt()<2.0) continue;
                        double fsrDr = deltaR(thisLep.Eta(), thisLep.Phi(), phot->eta(), phot->phi());
                        if (fsrDr>0.5) continue;

                        // check super cluster veto against all electrons for each photon
                        // at the same time check that this is the closest lepton for this photon
                        bool matched=false;
                        bool closest=true;                    

                        for (unsigned int j=0; j<Nleptons; j++) {                                            

                            TLorentzVector otherLep;
                            otherLep.SetPtEtaPhiM(lep_pt[j],lep_eta[j],lep_phi[j],lep_mass[j]);

                            double fsrDrOther = deltaR(otherLep.Eta(), otherLep.Phi(), phot->eta(), phot->phi());
                            if (j!=i && fsrDrOther<fsrDr) {closest=false;}

                            if ( abs(lep_id[(int)j])==11) {                                                               
 
                                for(size_t ecand = 0; ecand < recoElectrons[lep_ptindex[j]].associatedPackedPFCandidates().size(); ecand++){


                                    double ecandpt = recoElectrons[lep_ptindex[j]].associatedPackedPFCandidates()[ecand]->pt();
                                    double ecandeta = recoElectrons[lep_ptindex[j]].associatedPackedPFCandidates()[ecand]->eta();
                                    double ecandphi = recoElectrons[lep_ptindex[j]].associatedPackedPFCandidates()[ecand]->phi();
                                    if (abs(ecandpt-phot->pt())<1e-10 && abs(ecandeta-phot->eta())<1e-10 && abs(ecandphi-phot->phi())<1e-10) matched=true;

                                }                                
                            }
                        }
                        if (matched) continue;
                        if (!closest) continue;

                        // comput iso, dR/ET^2
                        double photoniso = helper.photonPfIso03(*phot,pfCands)/phot->pt();
                        double fsrDrOEt2 = fsrDr/(phot->pt()*phot->pt());

                        // fill all fsr photons before iso, dR/Et^2 cuts
                        allfsrPhotons_dR.push_back(fsrDr);
                        allfsrPhotons_iso.push_back(photoniso);
                        allfsrPhotons_pt.push_back(phot->pt());

                        // require photon iso, dR/Et^2
                        if (photoniso>1.8) continue;
                        if (fsrDrOEt2>0.012) continue;

                        // this photon is now a good one, check if it is the best one
                        
                        if ( verbose) cout<<"fsr photon cand, pt: "<<phot->pt()<<" eta: "<<phot->eta()<<" phi: "<<phot->phi()
                                          //<<" isoCHPUNoPU: "<<phot->userFloat("fsrPhotonPFIsoChHadPUNoPU03pt02")
                                          //<<" isoNHPhoton: "<<phot->userFloat("fsrPhotonPFIsoNHadPhoton03")
                                          <<" photoniso: "<<photoniso<<" DrOEt2: "<< fsrDrOEt2 <<endl;                                       

                        if( fsrDrOEt2 < minDrOEt2 ) {
                            selected = true;
                            selectedPhoton=(*phot);
                            selectedPhotonIso=photoniso;
                            selectedPhotonDr=fsrDr;
                            minDrOEt2 = fsrDrOEt2;
                            if (verbose) cout<<"****selected fsr: "<<i<<endl;
                            if (verbose) cout<<"photoniso: "<<photoniso<<" fsr dR over Et^2: "<<fsrDrOEt2<<" fsr dR: "<<fsrDr<<endl;
                        }
                        
                    } // all photons

                    if (selected) {
                        nFSRPhotons++;
                        fsrPhotons.push_back(selectedPhoton); 
                        fsrPhotons_pt.push_back(selectedPhoton.pt());
                        double perr = PFEnergyResolution().getEnergyResolutionEm(selectedPhoton.energy(), selectedPhoton.eta());
                        double pterr = perr*selectedPhoton.pt()/selectedPhoton.p();
                        fsrPhotons_pterr.push_back(pterr);
                        fsrPhotons_eta.push_back(selectedPhoton.eta());
                        fsrPhotons_phi.push_back(selectedPhoton.phi());
                        fsrPhotons_lepindex.push_back((int)i);
                        fsrPhotons_dR.push_back(selectedPhotonDr);
                        fsrPhotons_iso.push_back(selectedPhotonIso);
                        TLorentzVector phofsr;
                        phofsr.SetPtEtaPhiM(selectedPhoton.pt(),selectedPhoton.eta(),selectedPhoton.phi(),0.0);
                        TLorentzVector lepfsr;
                        lepfsr = thisLep+phofsr;
                        lepFSR_pt[i] = lepfsr.Pt();
                        lepFSR_eta[i] = lepfsr.Eta();
                        lepFSR_phi[i] = lepfsr.Phi();
                        lepFSR_mass[i] = lepfsr.M();

                        fsrmap[i] = phofsr;
                        if (verbose) cout<<"****selected fsr: "<<i<<endl;
                        if (verbose) cout<<"phofsr pt: "<<phofsr.Pt()<<" eta: "<<phofsr.Eta()<<" phi: "<<phofsr.Phi()<<" mass: "<<phofsr.M()<<endl;
                        if (verbose) cout<<"lep pt: "<<thisLep.Pt()<<" eta: "<<thisLep.Eta()<<" phi: "<<thisLep.Phi()<<" mass: "<<thisLep.M()<<endl;
                        if (verbose) cout<<"lep+fsr pt: "<<lepFSR_pt[i]<<" eta: "<<lepFSR_eta[i]<<" phi: "<<lepFSR_phi[i]<<" mass: "<<lepfsr.M()<<endl;
                        
 
                    }

                } // all leptons

                // subtract selected photons from all leptons isolations
                for (unsigned int i=0; i<Nleptons; i++) {
                        
                    TLorentzVector lep_nofsr;
                    lep_nofsr.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

                    double isoFSR=0.0;
                    for (unsigned int j=0; j<fsrPhotons.size(); j++) {
                        double fsrDr = deltaR(lep_nofsr.Eta(),lep_nofsr.Phi(),fsrPhotons[j].eta(),fsrPhotons[j].phi());
                        bool isoVeto=true;
                        if (abs(lep_id[i])==13 && fsrDr>0.01) isoVeto=false;
                        if (abs(lep_id[i])==11 && (abs(recoElectrons[lep_ptindex[i]].superCluster()->eta())<1.479 || fsrDr>0.08)) isoVeto=false;
                        if (fsrDr<((abs(lep_id[i])==11)?isoConeSizeEl:isoConeSizeMu) && !isoVeto) isoFSR += fsrPhotons[j].pt();
                    }

                    double RelIsoNoFSR = (lep_isoCH[i]+std::max(lep_isoNH[i]+lep_isoPhot[i]-lep_isoPUcorr[i]-isoFSR,0.0))/lep_nofsr.Pt();
                    lep_RelIsoNoFSR[i] = RelIsoNoFSR;               
                    if (verbose) cout<<"lep pt: "<<lep_nofsr.Pt()<<" eta: "<<lep_nofsr.Eta()<<" phi: "<<lep_nofsr.Phi()<<" RelIsoNoFSR: "<<RelIsoNoFSR<<" lep mva: "<<lep_mva[i]<<" tightId? "<<lep_tightId[i]<<endl;
                }
                if (verbose) {cout<<"finished filling fsr photon candidates"<<endl;}
            } // doFsrRecovery

            // count tight ID iso leptons
            uint ntight=0;
            for (unsigned int i=0; i<Nleptons; i++) {
                if (abs(lep_id[i])==11 && lep_RelIsoNoFSR[i]<isoCutEl && lep_tightId[i]==1) ntight+=1;
                if (abs(lep_id[i])==13 && lep_RelIsoNoFSR[i]<isoCutMu && lep_tightId[i]==1) ntight+=1;
            }

            if ( ntight >= (uint)skimTightLeptons ) {

                // Fake Rate Study (Z+1L Control Region)
                if (verbose) cout<<"begin Z+1L fake rate study"<<endl;
                // Z+1L selection
                findZ1LCandidate(iEvent);
                if (foundZ1LCandidate) { 
                    passedZ4lZ1LSelection = true;
                    if (passedTrig) passedZ1LSelection = true;
                }
                if (verbose) {cout<<"finished Z+1L fake rate study"<<endl;}

                // creat vectors for selected objects
                vector<pat::Muon> selectedMuons;
                vector<pat::Electron> selectedElectrons;
                
                if (verbose) cout<<"begin looking for higgs candidate"<<endl;                    
                findHiggsCandidate(selectedMuons,selectedElectrons,iEvent);
                if (verbose) {cout<<"found higgs candidate? "<<foundHiggsCandidate<<endl; }
                
                // Jets
                if (verbose) cout<<"begin filling jet candidates"<<endl;
                
                vector<pat::Jet> goodJets;
                vector<float> patJetQGTagger, patJetaxis2, patJetptD;
                vector<float> goodJetQGTagger, goodJetaxis2, goodJetptD; 
                vector<int> patJetmult, goodJetmult;
                
                for(auto jet = jets->begin();  jet != jets->end(); ++jet){
                    edm::RefToBase<pat::Jet> jetRef(edm::Ref<edm::View<pat::Jet> >(jets, jet - jets->begin()));
                    float qgLikelihood = (*qgHandle)[jetRef];
                    float axis2 = (*axis2Handle)[jetRef];
                    float ptD = (*ptDHandle)[jetRef];
                    int mult = (*multHandle)[jetRef];
                    patJetQGTagger.push_back(qgLikelihood);  
                    patJetaxis2.push_back(axis2);  
                    patJetmult.push_back(mult);  
                    patJetptD.push_back(ptD);  
                }
           
                
                for(unsigned int i = 0; i < jets->size(); ++i) {
                    
                    const pat::Jet & jet = jets->at(i);
                    
                    //JetID ID
                    if (verbose) cout<<"checking jetid..."<<endl;
                    float jpumva=0.;
                    bool passPU;
                    if (doJEC && (year==2017 || year==2018)) {
                        passPU = bool(jet.userInt("pileupJetIdUpdated:fullId") & (1 << 0));
                        jpumva=jet.userFloat("pileupJetIdUpdated:fullDiscriminant");
                    } else {
                        passPU = bool(jet.userInt("pileupJetId:fullId") & (1 << 0));
                        jpumva=jet.userFloat("pileupJetId:fullDiscriminant");
                    }
                    if (verbose) cout<< " jet pu mva  "<<jpumva <<endl;
                    /*
                    if(jet.pt()>20){
                        if(abs(jet.eta())>3.){
                            if(jpumva<=-0.21)passPU=false;
                        }else if(abs(jet.eta())>2.75){
                            if(jpumva<=-0.26)passPU=false;
                        }else if(abs(jet.eta())>2.5){
                            if(jpumva<=-0.35)passPU=false;
                        }else if(jpumva<=0.69)passPU=false;
                    }else{
                        if(abs(jet.eta())>3.){
                            if(jpumva<=-0.01)passPU=false;
                        }else if(abs(jet.eta())>2.75){
                            if(jpumva<=-0.05)passPU=false;
                        }else if(abs(jet.eta())>2.5){
                            if(jpumva<=-0.10)passPU=false;
                        }else if(jpumva<=0.86)passPU=false;
                    }
                    */                   
                        
                    if (verbose) cout<<"pt: "<<jet.pt()<<" eta: "<<jet.eta()<<" passPU: "<<passPU
                                     <<" jetid: "<<jetHelper.patjetID(jet,year)<<endl;
                    
                    ///cut is seperated into two part since we will apply corrected pt into passPU|| jet.pt()>50; another part in L3916
                    //cut = cms.string("pt>20 && abs(eta)<4.7 && userFloat('JetID') && (userFloat('PUjetID') || pt>50)"),
                    //if( jetHelper.patjetID(jet,year)>=jetIDLevel && (passPU || !doPUJetID || jet.pt()>50) ) {
                    if( jetHelper.patjetID(jet,year)>=jetIDLevel ) {  
                        
                        if (verbose) cout<<"passed pf jet id and pu jet id"<<endl;
                        
                        // apply loose pt cut here (10 GeV cut is already applied in MINIAOD) since we are before JES/JER corrections
                        //if(jet.pt() > 10.0 && fabs(jet.eta()) < jeteta_cut) {
                        if(fabs(jet.eta()) < jeteta_cut) { ///move all pt cut after JES/JER corrections
                            
                            // apply scale factor for PU Jets by demoting 1-data/MC % of jets jets in certain pt/eta range 
                            // Configured now that SF is 1.0
                            if (verbose) cout<<"adding pu jet scale factors..."<<endl;       
                            bool dropit=false;
                            if (abs(jet.eta())>3.0 && isMC) {
                                TRandom3 rand;
                                rand.SetSeed(abs(static_cast<int>(sin(jet.phi())*100000)));
                                float coin = rand.Uniform(1.); 
                                if (jet.pt()>=20.0 && jet.pt()<36.0 && coin>1.0) dropit=true;
                                if (jet.pt()>=36.0 && jet.pt()<50.0 && coin>1.0) dropit=true;
                                if (jet.pt()>=50.0 && coin>1.0) dropit=true;
                            }                                        
                            
                            if (!dropit) {
                                if (verbose) cout<<"adding jet candidate, pt: "<<jet.pt()<<" eta: "<<jet.eta()<<endl;
                                goodJets.push_back(jet);
                                goodJetQGTagger.push_back(patJetQGTagger[i]);
                                goodJetaxis2.push_back(patJetaxis2[i]);
                                goodJetptD.push_back(patJetptD[i]);
                                goodJetmult.push_back(patJetmult[i]);
                            } // pu jet scale factor
                            
                        } // pass loose pt cut 
                    } // pass loose pf jet id and pu jet id
                } // all jets
                
                vector<pat::Jet> selectedMergedJets;
                /*
                for(unsigned int i = 0; i < mergedjets->size(); ++i) { // all merged jets
                    
                    const pat::Jet & mergedjet = mergedjets->at(i);
                    double pt = double(mergedjet.pt());
                    double eta = double(mergedjet.eta());
                    
                    if(pt>60 && abs(eta)<2.5) selectedMergedJets.push_back(mergedjet);
                    
                } // all merged jets

                */

                
                if( foundHiggsCandidate ){
                    
                    for(unsigned int i = 0; i<4;i++){
                        
                        int index = lep_Hindex[i];
                        if(fsrmap[index].Pt()!=0){
                            if (verbose) cout<<"find a fsr photon for "<<i<<" th Higgs lepton"<<endl;
                            selectedFsrMap[i] = fsrmap[index];
                        }
                    }
                    
                    
                    if (verbose) cout<<"storing H_p4_noFSR"<<endl; 
                    math::XYZTLorentzVector tmpHVec;
                    if (verbose) cout<<"selectedMuons "<<selectedMuons.size()<<" selectedElectrons "<<selectedElectrons.size()<<endl;
                    if(RecoFourMuEvent) {
                        tmpHVec = selectedMuons[0].p4() + selectedMuons[1].p4() + selectedMuons[2].p4() + selectedMuons[3].p4();
                        HVecNoFSR.SetPtEtaPhiM(tmpHVec.Pt(),tmpHVec.Eta(),tmpHVec.Phi(),tmpHVec.M());
                        
                        reco::Candidate *c1 = dynamic_cast<reco::Candidate* >(&selectedMuons[0]);
                        selectedLeptons.push_back(c1);
                        reco::Candidate *c2 = dynamic_cast<reco::Candidate* >(&selectedMuons[1]);
                        selectedLeptons.push_back(c2);
                        reco::Candidate *c3 = dynamic_cast<reco::Candidate* >(&selectedMuons[2]);
                        selectedLeptons.push_back(c3);
                        reco::Candidate *c4 = dynamic_cast<reco::Candidate* >(&selectedMuons[3]);
                        selectedLeptons.push_back(c4);
                        
                    }
                    else if(RecoFourEEvent) {
                        tmpHVec = selectedElectrons[0].p4() + selectedElectrons[1].p4() + selectedElectrons[2].p4() + selectedElectrons[3].p4();
                        HVecNoFSR.SetPtEtaPhiM(tmpHVec.Pt(),tmpHVec.Eta(),tmpHVec.Phi(),tmpHVec.M());
                        
                        reco::Candidate *c1 = dynamic_cast<reco::Candidate* >(&selectedElectrons[0]);
                        selectedLeptons.push_back(c1);
                        reco::Candidate *c2 = dynamic_cast<reco::Candidate* >(&selectedElectrons[1]);
                        selectedLeptons.push_back(c2);
                        reco::Candidate *c3 = dynamic_cast<reco::Candidate* >(&selectedElectrons[2]);
                        selectedLeptons.push_back(c3);
                        reco::Candidate *c4 = dynamic_cast<reco::Candidate* >(&selectedElectrons[3]);
                        selectedLeptons.push_back(c4);
                        
                    }
                    else if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){
                        tmpHVec = selectedMuons[0].p4() + selectedMuons[1].p4() + selectedElectrons[0].p4() + selectedElectrons[1].p4();
                        HVecNoFSR.SetPtEtaPhiM(tmpHVec.Pt(),tmpHVec.Eta(),tmpHVec.Phi(),tmpHVec.M());
                        
                        if(RecoTwoETwoMuEvent){
                            
                            reco::Candidate *c1 = dynamic_cast<reco::Candidate* >(&selectedElectrons[0]);
                            selectedLeptons.push_back(c1);
                            reco::Candidate *c2 = dynamic_cast<reco::Candidate* >(&selectedElectrons[1]);
                            selectedLeptons.push_back(c2);
                            reco::Candidate *c3 = dynamic_cast<reco::Candidate* >(&selectedMuons[0]);
                            selectedLeptons.push_back(c3);
                            reco::Candidate *c4 = dynamic_cast<reco::Candidate* >(&selectedMuons[1]);
                            selectedLeptons.push_back(c4);
                        }
                        else{
                            
                            reco::Candidate *c1 = dynamic_cast<reco::Candidate* >(&selectedMuons[0]);
                            selectedLeptons.push_back(c1);
                            reco::Candidate *c2 = dynamic_cast<reco::Candidate* >(&selectedMuons[1]);
                            selectedLeptons.push_back(c2);
                            reco::Candidate *c3 = dynamic_cast<reco::Candidate* >(&selectedElectrons[0]);
                            selectedLeptons.push_back(c3);
                            reco::Candidate *c4 = dynamic_cast<reco::Candidate* >(&selectedElectrons[1]);
                            selectedLeptons.push_back(c4);
                        }
                        
                    }
                    
					//set vertex constraint for 4mu events
					if(Hvtx){
						if(RecoFourMuEvent){
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedMuons[0].muonBestTrack()));
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedMuons[1].muonBestTrack()));
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedMuons[2].muonBestTrack()));
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedMuons[3].muonBestTrack()));
						}
						if(RecoTwoMuTwoEEvent){
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedMuons[0].muonBestTrack()));
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedMuons[1].muonBestTrack()));
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedElectrons[0].gsfTrack()));
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedElectrons[1].gsfTrack()));
						}
						if(RecoTwoETwoMuEvent){
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedElectrons[0].gsfTrack()));
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedElectrons[1].gsfTrack()));
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedMuons[0].muonBestTrack()));
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedMuons[1].muonBestTrack()));
						}
						if(RecoFourEEvent){
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedElectrons[0].gsfTrack()));
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedElectrons[1].gsfTrack()));
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedElectrons[2].gsfTrack()));
							H4Lrefittedtracks.push_back(trackbuilder->build(selectedElectrons[3].gsfTrack()));//
						}
						if(H4Lrefittedtracks.size()==4){
							KalmanVertexFitter H4Lfitter(true);
							TransientVertex H4Lvertex = H4Lfitter.vertex(H4Lrefittedtracks);
							if(H4Lvertex.isValid()&&H4Lvertex.hasRefittedTracks()&&H4Lvertex.refittedTracks().size()==4){
								track_Hvtx_xcorr = H4Lvertex.position().x();
								track_Hvtx_ycorr = H4Lvertex.position().y();
								track_Hvtx_zcorr = H4Lvertex.position().z();
								GlobalPoint Hver(pv_x, pv_y, pv_z);
								for(unsigned int k = 0; k < H4Lvertex.refittedTracks().size(); k++){
									TrajectoryStateClosestToPoint Htrajs = H4Lvertex.refittedTracks().at(k).trajectoryStateClosestToPoint(Hver);
									track_Hvtx_dxycorr.push_back(Htrajs.perigeeParameters().transverseImpactParameter());
									track_Hvtx_pt.push_back(H4Lvertex.refittedTracks().at(k).track().pt());
									track_Hvtx_eta.push_back(H4Lvertex.refittedTracks().at(k).track().eta());
									track_Hvtx_phi.push_back(H4Lvertex.refittedTracks().at(k).track().phi());
									track_Hvtx_pterr.push_back(H4Lvertex.refittedTracks().at(k).track().ptError());
								}
							}
							TransientVertex H4Lbeamspot = H4Lfitter.vertex(H4Lrefittedtracks, BS);
							if(H4Lbeamspot.isValid()&&H4Lbeamspot.hasRefittedTracks()&&H4Lbeamspot.refittedTracks().size()==4){
								track_Hbs_xcorr = H4Lbeamspot.position().x();
								track_Hbs_ycorr = H4Lbeamspot.position().y();
								track_Hbs_zcorr = H4Lbeamspot.position().z();
								GlobalPoint Hver(pv_x, pv_y, pv_z);
								for(unsigned int k = 0; k < H4Lbeamspot.refittedTracks().size(); k++){
									TrajectoryStateClosestToPoint Htrajs = H4Lbeamspot.refittedTracks().at(k).trajectoryStateClosestToPoint(Hver);
									track_Hbs_dxycorr.push_back(Htrajs.perigeeParameters().transverseImpactParameter());
									track_Hbs_pt.push_back(H4Lbeamspot.refittedTracks().at(k).track().pt());
									track_Hbs_eta.push_back(H4Lbeamspot.refittedTracks().at(k).track().eta());
									track_Hbs_phi.push_back(H4Lbeamspot.refittedTracks().at(k).track().phi());
									track_Hbs_pterr.push_back(H4Lbeamspot.refittedTracks().at(k).track().ptError());
								}
							}
							if(Checkvtx&&H4Lvertex.isValid()&&H4Lvertex.hasRefittedTracks()&&(H4Lvertex.refittedTracks().size()==4&&RecoFourMuEvent)){//only check H4mu events
								std::cout<<"!!!!!!!!!!!!!!!!!!"<<std::endl;
								std::cout<<"H4muevent oldphi: "<<lep_phi[lep_Hindex[0]]<<std::endl;
								std::cout<<"H4muevent newphi: "<<track_Hvtx_phi[0]<<std::endl;
							}
						}	
					}
				
                    H_noFSR_pt.push_back(HVecNoFSR.Pt());
                    H_noFSR_eta.push_back(HVecNoFSR.Eta());
                    H_noFSR_phi.push_back(HVecNoFSR.Phi());
                    H_noFSR_mass.push_back(HVecNoFSR.M());
                    
                    // check number of failing leptons
                    for(unsigned int i = 0; i <= 3; i++) {
                        if (!(abs(lep_id[lep_Hindex[i]])==11 && (lep_tightId[lep_Hindex[i]] && lep_RelIsoNoFSR[lep_Hindex[i]]<isoCutEl)) &&
                            !(abs(lep_id[lep_Hindex[i]])==13 && (lep_tightId[lep_Hindex[i]] && lep_RelIsoNoFSR[lep_Hindex[i]]<isoCutMu))){ nZXCRFailedLeptons++; }
                    }
                    if (verbose) cout << nZXCRFailedLeptons<<" failing leptons in higgs candidate"<<endl;
                    if (nZXCRFailedLeptons>0) { // at least one lepton has failed 
                        passedZ4lZXCRSelection = true;
                        if (Z2Vec.M() > mZ2Low && passedTrig) passedZXCRSelection = true;
                    } else { //  signal region candidate                    
                        passedZ4lSelection = true;
                        if(Z2Vec.M() > mZ2Low && passedTrig) passedFullSelection = true;
                    }
                    
                } // found higgs candidate 
                else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed higgs candidate"<<endl;}
                 
                
                //Set All the Variables for Saved Trees (after finding higgs candidate)
                if (verbose) cout<<"begin setting tree variables"<<endl;
                setTreeVariables(iEvent, iSetup, selectedMuons, selectedElectrons, recoMuons, recoElectrons, goodJets, goodJetQGTagger,goodJetaxis2, goodJetptD, goodJetmult, selectedMergedJets);
                if (verbose) cout<<"finshed setting tree variables"<<endl;
            
                
                // Comput Matrix Elelements After filling jets, Do Kinematic fit, add scale factors
                //if (foundHiggsCandidate || lep_pt.size()>=4) {
                if (foundHiggsCandidate) {
                    
                    if (foundHiggsCandidate) {
                        dataMCWeight = lep_dataMC[lep_Hindex[0]]*lep_dataMC[lep_Hindex[1]]*lep_dataMC[lep_Hindex[2]]*lep_dataMC[lep_Hindex[3]];
                    } else {
                        dataMCWeight = 1.0;
                    }
                    eventWeight = crossSection*pileupWeight*dataMCWeight*prefiringWeight;
                    
                        
                    dataMCWeight = lep_dataMC[lep_Hindex[0]]*lep_dataMC[lep_Hindex[1]]*lep_dataMC[lep_Hindex[2]]*lep_dataMC[lep_Hindex[3]];
                    eventWeight = genWeight*crossSection*pileupWeight*dataMCWeight*prefiringWeight;
                    
                    if (verbose) cout<<"Kin fitter begin with lep size "<<selectedLeptons.size()<<" fsr size "<<selectedFsrMap.size()<<endl;
                    
                    kinZfitter->Setup(selectedLeptons, selectedFsrMap);
                    mass4lErr = (float)kinZfitter->GetM4lErr();

                    if (doRefit) {
                        kinZfitter->KinRefitZ();
                        mass4lREFIT = (float)kinZfitter->GetRefitM4l();
                        mass4lErrREFIT = (float)kinZfitter->GetRefitM4lErrFullCov();
                        massZ1REFIT = (float)kinZfitter->GetRefitMZ1(); 
                        massZ2REFIT = (float)kinZfitter->GetRefitMZ2(); 
                    }

                    if (verbose) cout<<"mass4l "<<mass4l<<" mass4lREFIT "<<mass4lREFIT<<" massErr "<<mass4lErr<<" massErrREFIT "<<mass4lErrREFIT<<" massZ1REFIT "<<massZ1REFIT<<endl;
                }
                
           
                //if (doMela) {
                if (doMela && foundHiggsCandidate) {
                    
                    TLorentzVector Lep1, Lep2, Lep3, Lep4,  Jet1, Jet2;
                    if (foundHiggsCandidate) {
                        Lep1.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[0]],lepFSR_eta[lep_Hindex[0]],lepFSR_phi[lep_Hindex[0]],lepFSR_mass[lep_Hindex[0]]);
                        Lep2.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[1]],lepFSR_eta[lep_Hindex[1]],lepFSR_phi[lep_Hindex[1]],lepFSR_mass[lep_Hindex[1]]);
                        Lep3.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[2]],lepFSR_eta[lep_Hindex[2]],lepFSR_phi[lep_Hindex[2]],lepFSR_mass[lep_Hindex[2]]);
                        Lep4.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[3]],lepFSR_eta[lep_Hindex[3]],lepFSR_phi[lep_Hindex[3]],lepFSR_mass[lep_Hindex[3]]);
                    } else {
                        Lep1.SetPtEtaPhiM(lepFSR_pt[0],lepFSR_eta[0],lepFSR_phi[0],lepFSR_mass[0]);
                        Lep2.SetPtEtaPhiM(lepFSR_pt[1],lepFSR_eta[1],lepFSR_phi[1],lepFSR_mass[1]);
                        Lep3.SetPtEtaPhiM(lepFSR_pt[2],lepFSR_eta[2],lepFSR_phi[2],lepFSR_mass[2]);
                        Lep4.SetPtEtaPhiM(lepFSR_pt[3],lepFSR_eta[3],lepFSR_phi[3],lepFSR_mass[3]);
                    }
                    
                    SimpleParticleCollection_t daughters;
                    if (foundHiggsCandidate) {
                        daughters.push_back(SimpleParticle_t(lep_id[lep_Hindex[0]], Lep1));
                        daughters.push_back(SimpleParticle_t(lep_id[lep_Hindex[1]], Lep2));
                        daughters.push_back(SimpleParticle_t(lep_id[lep_Hindex[2]], Lep3));
                        daughters.push_back(SimpleParticle_t(lep_id[lep_Hindex[3]], Lep4));
                    } else {
                        daughters.push_back(SimpleParticle_t(lep_id[0], Lep1));
                        daughters.push_back(SimpleParticle_t(lep_id[1], Lep2));
                        daughters.push_back(SimpleParticle_t(lep_id[2], Lep3));
                        daughters.push_back(SimpleParticle_t(lep_id[3], Lep4));
                    }
                    
                    
                    SimpleParticleCollection_t associated;               
                    if (njets_pt30_eta4p7 > 0) {
                        Jet1.SetPtEtaPhiM(jet_pt[jet1index],jet_eta[jet1index],jet_phi[jet1index],jet_mass[jet1index]);
                        associated.push_back(SimpleParticle_t(0, Jet1));
                    }
                    if (njets_pt30_eta4p7 > 1) {
                        Jet2.SetPtEtaPhiM(jet_pt[jet2index],jet_eta[jet2index],jet_phi[jet2index],jet_mass[jet2index]);
                        associated.push_back(SimpleParticle_t(0, Jet2));
                    }
                
                    mela->setInputEvent(&daughters, &associated, 0, 0);
                    mela->setCurrentCandidateFromIndex(0);
                    
                    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
                    mela->computeP(me_0plus_JHU, true);
                    
                    mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
                    mela->computeP(p0minus_VAJHU, true);
                    
                    pg1g4_VAJHU=0.0;
                    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
                    (mela->selfDHggcoupl)[0][0][0]=1.;
                    (mela->selfDHzzcoupl)[0][0][0]=1.;
                    (mela->selfDHzzcoupl)[0][3][0]=1.;
                    mela->computeP(pg1g4_VAJHU, true);
                    pg1g4_VAJHU -= me_0plus_JHU+p0minus_VAJHU;
                    
                    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
                    mela->computeP(me_qqZZ_MCFM, true);
                    
                    mela->computeD_gg(TVar::MCFM, TVar::D_gg10, Dgg10_VAMCFM);
                    
                    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
                    mela->computePM4l(TVar::SMSyst_None, p0plus_m4l);
                    
                    mela->setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
                    mela->computePM4l(TVar::SMSyst_None, bkg_m4l);

                    
                    D_bkg_kin = me_0plus_JHU/(me_0plus_JHU+me_qqZZ_MCFM*helper.getDbkgkinConstant(idL1*idL2*idL3*idL3,mass4l)); 
                    D_bkg = me_0plus_JHU*p0plus_m4l/(me_0plus_JHU*p0plus_m4l+me_qqZZ_MCFM*bkg_m4l*helper.getDbkgConstant(idL1*idL2*idL3*idL4,mass4l)); // superMELA 
                    D_g4 = me_0plus_JHU/(me_0plus_JHU+pow(2.521, 2)*p0minus_VAJHU); // D_0-                
                    D_g1g4 = pg1g4_VAJHU*2.521/(me_0plus_JHU+pow(2.521, 2)*p0minus_VAJHU); // D_CP, 2.521 since g1=1 and g4=1 is used
                    
                    TUtil::computeAngles(cosThetaStar,cosTheta1,cosTheta2,Phi,Phi1, \
                                         Lep1, lep_id[lep_Hindex[0]], Lep2, lep_id[lep_Hindex[1]], \
                                         Lep3, lep_id[lep_Hindex[2]], Lep4, lep_id[lep_Hindex[3]]);
                    
                    if (njets_pt30_eta4p7>=2){

                        mela->setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        mela->computeProdDecP(p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal);
                        mela->getConstant(pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal);

                        mela->setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::Had_ZH_S);
                        mela->computeProdDecP(p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal);
                        mela->getConstant(pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal);

                        mela->setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::Had_WH_S);
                        mela->computeProdDecP(p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal);
                        mela->getConstant(pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal);

                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
                        mela->computeProdP(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
                        mela->getConstant(pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);

                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
                        mela->computeProdP(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal);
                        mela->getConstant(pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal);

                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
                        mela->computeDijetConvBW(p_HadWH_mavjj_JECNominal, false);

                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
                        mela->computeDijetConvBW(p_HadWH_mavjj_true_JECNominal, true);

                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
                        mela->computeProdP(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal);
                        mela->getConstant(pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal);

                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
                        mela->computeDijetConvBW(p_HadZH_mavjj_JECNominal, false);

                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
                        mela->computeDijetConvBW(p_HadZH_mavjj_true_JECNominal, true);

                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
                        mela->computeProdP(p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);
                        mela->getConstant(pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);

                        mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::JJVBF);
                        mela->computeProdDecP(p_JJVBF_BKG_MCFM_JECNominal);
                        mela->getConstant(pConst_JJVBF_BKG_MCFM_JECNominal);

                        mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::Had_ZH);
                        mela->computeProdDecP(p_HadZH_BKG_MCFM_JECNominal);
                        mela->getConstant(pConst_HadZH_BKG_MCFM_JECNominal);

                        mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::Had_WH);
                        mela->computeProdDecP(p_HadWH_BKG_MCFM_JECNominal);
                        mela->getConstant(pConst_HadWH_BKG_MCFM_JECNominal);

                        mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::JJQCD);
                        mela->computeProdDecP(p_JJQCD_BKG_MCFM_JECNominal);
                        mela->getConstant(pConst_JJQCD_BKG_MCFM_JECNominal);

                        
                        /*
                        std::cout<<"p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal "<<p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal<<std::endl;
                        std::cout<<"pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal "<<pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal<<std::endl;
                        std::cout<<"p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal "<<p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal<<std::endl;
                        std::cout<<"pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal "<<pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal<<std::endl;
                        std::cout<<"p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal "<<p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal<<std::endl;
                        std::cout<<"pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal "<<pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal<<std::endl;
                        std::cout<<"p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal "<<p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal<<std::endl;
                        std::cout<<"pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal "<<pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal<<std::endl;
                        std::cout<<"p_HadWH_SIG_ghw1_1_JHUGen_JECNominal "<<p_HadWH_SIG_ghw1_1_JHUGen_JECNominal<<std::endl;
                        std::cout<<"pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal "<<pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal<<std::endl;
                        std::cout<<"p_HadZH_SIG_ghz1_1_JHUGen_JECNominal "<<p_HadZH_SIG_ghz1_1_JHUGen_JECNominal<<std::endl;
                        std::cout<<"pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal "<<pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal<<std::endl;
                        std::cout<<"p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal "<<p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal<<std::endl;
                        std::cout<<"pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal "<<pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal<<std::endl;
                        std::cout<<"p_JJVBF_BKG_MCFM_JECNominal "<<p_JJVBF_BKG_MCFM_JECNominal<<std::endl;
                        std::cout<<"pConst_JJVBF_BKG_MCFM_JECNominal "<<pConst_JJVBF_BKG_MCFM_JECNominal<<std::endl;
                        std::cout<<"p_HadZH_BKG_MCFM_JECNominal "<<p_HadZH_BKG_MCFM_JECNominal<<std::endl;
                        std::cout<<"pConst_HadZH_BKG_MCFM_JECNominal "<<pConst_HadZH_BKG_MCFM_JECNominal<<std::endl;
                        std::cout<<"p_HadWH_BKG_MCFM_JECNominal "<<p_HadWH_BKG_MCFM_JECNominal<<std::endl;
                        std::cout<<"pConst_HadWH_BKG_MCFM_JECNominal "<<pConst_HadWH_BKG_MCFM_JECNominal<<std::endl;
                        std::cout<<"p_JJQCD_BKG_MCFM_JECNominal "<<p_JJQCD_BKG_MCFM_JECNominal<<std::endl;
                        std::cout<<"pConst_JJQCD_BKG_MCFM_JECNominal "<<pConst_JJQCD_BKG_MCFM_JECNominal<<std::endl;
                        */
                        

                        D_VBF = 1./(1.+ helper.getDVBF2jetsConstant(mass4l)*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
                        D_HadWH = 1./(1.+ helper.getDWHhConstant(mass4l)*(p_HadWH_mavjj_true_JECNominal*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_HadWH_mavjj_JECNominal*p_HadWH_SIG_ghw1_1_JHUGen_JECNominal));
                        D_HadZH =  1./(1.+ helper.getDZHhConstant(mass4l)*(p_HadZH_mavjj_true_JECNominal*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_HadZH_mavjj_JECNominal*p_HadZH_SIG_ghz1_1_JHUGen_JECNominal));

                        
                        //D_bkg_VBFdec
                        float DbkgVBFdecConstant = helper.getDbkgVBFdecConstant(idL1*idL2*idL3*idL3,mass4l);

                        float vbf = p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal/pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal;
                        float zh = p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal/pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal;
                        float wh = p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal/pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal;
                        float constA = 1./(1./pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal+1./pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal+1./pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal);

                        float vbs = p_JJVBF_BKG_MCFM_JECNominal/pConst_JJVBF_BKG_MCFM_JECNominal;
                        float zzz = p_HadZH_BKG_MCFM_JECNominal/pConst_HadZH_BKG_MCFM_JECNominal;
                        float wzz = p_HadWH_BKG_MCFM_JECNominal/pConst_HadWH_BKG_MCFM_JECNominal;
                        float qcdzz = p_JJQCD_BKG_MCFM_JECNominal/pConst_JJQCD_BKG_MCFM_JECNominal;
                        float constB = 1./(1./pConst_JJVBF_BKG_MCFM_JECNominal+1./pConst_HadZH_BKG_MCFM_JECNominal+1./pConst_HadWH_BKG_MCFM_JECNominal+1./pConst_JJQCD_BKG_MCFM_JECNominal);

                        const float scale_Pmjj_vb=1;
                        float scale_Pmjj_z = p_HadZH_mavjj_JECNominal/p_HadZH_mavjj_true_JECNominal;
                        float scale_Pmjj_w = p_HadWH_mavjj_JECNominal/p_HadWH_mavjj_true_JECNominal;

                        vbf *= scale_Pmjj_vb;
                        vbs *= scale_Pmjj_vb;

                        zh *= scale_Pmjj_z;
                        zzz *= scale_Pmjj_z;

                        wh *= scale_Pmjj_w;
                        wzz *= scale_Pmjj_w;


                        float PA = (vbf + zh + wh)*constA;
                        float PB = (vbs + zzz + wzz + qcdzz)*constB;

                        //std::cout<<"DbkgVBFdecConstant: "<<DbkgVBFdecConstant<<std::endl;
                        //std::cout<<"PA: "<<PA<<std::endl;
                        //std::cout<<"PB: "<<PB<<std::endl;

                        D_bkg_VBFdec =  PA/(PA+DbkgVBFdecConstant*PB);


                        // D_bkg_VHdec
                        float DbkgVHdecConstant = helper.getDbkgVHdecConstant(idL1*idL2*idL3*idL3,mass4l);

                        vbf = p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal/pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal;
                        zh = p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal/pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal;
                        wh = p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal/pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal;
                        constA = 1./(1./pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal+1./pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal+1./pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal);

                        vbs = p_JJVBF_BKG_MCFM_JECNominal/pConst_JJVBF_BKG_MCFM_JECNominal;
                        zzz = p_HadZH_BKG_MCFM_JECNominal/pConst_HadZH_BKG_MCFM_JECNominal;
                        wzz = p_HadWH_BKG_MCFM_JECNominal/pConst_HadWH_BKG_MCFM_JECNominal;
                        qcdzz = p_JJQCD_BKG_MCFM_JECNominal/pConst_JJQCD_BKG_MCFM_JECNominal;
                        constB = 1./(1./pConst_JJVBF_BKG_MCFM_JECNominal+1./pConst_HadZH_BKG_MCFM_JECNominal+1./pConst_HadWH_BKG_MCFM_JECNominal+1./pConst_JJQCD_BKG_MCFM_JECNominal);

                        //scale_Pmjj_vb=1;
                        scale_Pmjj_z = p_HadZH_mavjj_JECNominal/p_HadZH_mavjj_true_JECNominal;
                        scale_Pmjj_w = p_HadWH_mavjj_JECNominal/p_HadWH_mavjj_true_JECNominal;

                        vbf *= scale_Pmjj_vb;
                        vbs *= scale_Pmjj_vb;

                        zh *= scale_Pmjj_z;
                        zzz *= scale_Pmjj_z;

                        wh *= scale_Pmjj_w;
                        wzz *= scale_Pmjj_w;


                        PA = (vbf + zh + wh)*constA;
                        PB = (vbs + zzz + wzz + qcdzz)*constB;

                        D_bkg_VHdec = PA/(PA+DbkgVHdecConstant*PB);


                        /* OLD BUT WORKING
                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
                        mela->computeProdP(pvbf_VAJHU, true);
                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
                        mela->computeProdP(phjj_VAJHU, true);
                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
                        mela->computeProdP(pwh_hadronic_VAJHU, true);
                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
                        mela->computeProdP(pzh_hadronic_VAJHU, true);
                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
                        mela->computeDijetConvBW(p_HadWH_mavjj_JECNominal, false);
                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
                        mela->computeDijetConvBW(p_HadWH_mavjj_true_JECNominal, true);
                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
                        mela->computeDijetConvBW(p_HadZH_mavjj_JECNominal, false);
                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
                        mela->computeDijetConvBW(p_HadZH_mavjj_true_JECNominal, true);

                        D_VBF = pvbf_VAJHU/(pvbf_VAJHU+phjj_VAJHU*helper.getDVBF2jetsConstant(mass4l) ); // VBF(2j) vs. gg->H+2j
                        //D_VBF =  1./(1.+ helper.getDVBF2jetsConstant(mass4l)*pvbf_VAJHU+phjj_VAJHU/pvbf_VAJHU);
                        
                        //D_HadWH = pwh_hadronic_VAJHU/(pwh_hadronic_VAJHU+phjj_VAJHU*helper.getDWHhConstant(mass4l) ); // W(->2j)H vs. gg->H+2j
                        //D_HadZH = pzh_hadronic_VAJHU/(pzh_hadronic_VAJHU+phjj_VAJHU*helper.getDZHhConstant(mass4l) ); // Z(->2j)H vs. gg->H+2j

                        D_HadWH =  1./(1.+ helper.getDWHhConstant(mass4l)*(p_HadWH_mavjj_true_JECNominal*phjj_VAJHU)/(p_HadWH_mavjj_JECNominal*pwh_hadronic_VAJHU));
                        D_HadZH =  1./(1.+ helper.getDZHhConstant(mass4l)*(p_HadZH_mavjj_true_JECNominal*phjj_VAJHU)/(p_HadZH_mavjj_JECNominal*pzh_hadronic_VAJHU));
                        */

                    } else {
                        D_VBF = -1.0; D_HadWH = -1.0; D_HadZH = -1.0;  D_bkg_VBFdec=-1.0; D_bkg_VHdec=-1.0;
                    }
                    
                    if (njets_pt30_eta4p7==1) {
                        
                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
                        mela->computeProdP(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,true);
                        mela->getPAux(pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal);

                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JQCD);
                        mela->computeProdP(p_JQCD_SIG_ghv1_1_JHUGen_JECNominal,true);

                        D_VBF1j = 1./(1.+ helper.getDVBF1jetConstant(mass4l)*p_JQCD_SIG_ghv1_1_JHUGen_JECNominal/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal));

                        /*
                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
                        mela->computeProdP(pvbf_VAJHU, true); // Un-integrated ME
                        mela->getPAux(pAux_vbf_VAJHU); // = Integrated / un-integrated

                        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JQCD);
                        mela->computeProdP(phj_VAJHU, true);

                        D_VBF1j = pvbf_VAJHU*pAux_vbf_VAJHU/(pvbf_VAJHU*pAux_vbf_VAJHU+phj_VAJHU*helper.getDVBF1jetConstant(mass4l)); // VBF(1j) vs. gg->H+1j
                        */

                    } else {
                        D_VBF1j = -1.0;
                    }
                    
                    if (njets_pt30_eta4p7>=2) {                    
                        float jetqgl0 =jet_QGTagger[jet1index]; 
                        float jetqgl1 =jet_QGTagger[jet2index]; 
                        if(jetqgl0<0.){ // if the q/g tagger has the error value (-1.), take a random one instead
                            TRandom3 rand;
                            rand.SetSeed(abs(static_cast<int>(sin(jet_phi[jet1index])*100000)));
                            jetqgl0 = rand.Uniform();
                        }
                        if(jetqgl1<0.){ // if the q/g tagger has the error value (-1.), take a random one instead
                            TRandom3 rand;
                            rand.SetSeed(abs(static_cast<int>(sin(jet_phi[jet2index])*100000)));
                            jetqgl1 = rand.Uniform();
                        }
                        float jetPgOverPq0 = 1./jetqgl0- 1.;
                        float jetPgOverPq1 = 1./jetqgl1- 1.;
                        
                        D_VBF_QG = 1./(1.+ (1./D_VBF - 1.) * pow(jetPgOverPq0*jetPgOverPq1, 1./3.));
                        D_HadWH_QG = 1./(1.+ (1./D_HadWH - 1.) * pow(jetPgOverPq0*jetPgOverPq1, 1./3.));
                        D_HadZH_QG = 1./(1. + (1./D_HadZH - 1.) * pow(jetPgOverPq0*jetPgOverPq1, 1./3.));
                    } else {
                        D_VBF_QG = -1.0; D_HadWH_QG = -1.0; D_HadZH_QG = -1.0;
                    }
                    
                    if (njets_pt30_eta4p7==1) {
                        float jetqgl0 =jet_QGTagger[jet1index]; 
                        if(jetqgl0<0.){ // if the q/g tagger has the error value (-1.), take a random one instead
                            TRandom3 rand;
                            rand.SetSeed(abs(static_cast<int>(sin(jet_phi[jet1index])*100000)));
                            jetqgl0 = rand.Uniform();
                        }
                        float jetPgOverPq0 = 1./jetqgl0- 1.;
                        D_VBF1j_QG = 1/(1+ (1./D_VBF1j - 1.) * pow(jetPgOverPq0, 1./3.));
                    } else {
                        D_VBF1j_QG = -1.0;
                    }
                    
                    if (verbose) cout<<"D_bkg_kin: "<<D_bkg_kin<< ", D_bkg: " << D_bkg << ", Dgg: " << Dgg10_VAMCFM << " ,D0-: " << D_g4 << endl;               
                    if (verbose) cout<<"D_VBF : "<<D_VBF<< ", D_VBF1j : "<< D_VBF1j<<", WH: " << D_HadWH << ", ZH: " << D_HadZH <<endl;
                    if (verbose) cout<<"cosThetaStar: "<<cosThetaStar<<", cosTheta1: "<<cosTheta1<<", cosTheta2: "<<cosTheta2<<", Phi: "<<Phi<<" , Phi1: "<<Phi1<<endl;
                    
                    mela->resetInputEvent(); 
                    
                    if(njets_pt30_eta4p7>1){
                        TLorentzVector jet1, jet2;
                        jet1.SetPtEtaPhiM(jet_pt[jet1index],jet_eta[jet1index],jet_phi[jet1index],jet_mass[jet1index]);
                        jet2.SetPtEtaPhiM(jet_pt[jet2index],jet_eta[jet2index],jet_phi[jet2index],jet_mass[jet2index]);
                        TLorentzVector Dijet;
                        Dijet = jet1+jet2; 
                        DijetMass = Dijet.M();
                        DijetDEta = fabs(jet1.Eta()-jet2.Eta());
                        // OLD MORIOND --- FisherDiscrim = 0.09407*fabs(VBFDeltaEta) + 4.1581e-4*VBFDiJetMass;
                        DijetFisher = 0.18*fabs(DijetDEta) + 1.92e-4*DijetMass;
                    }
                    
                    // Double loop over jets, for V-jet tagging
                    for (int i=0; i<njets_pt30_eta4p7; i++) {
                        for (int j=i+1; j<njets_pt30_eta4p7; j++) {
                            if (i==j) continue;
                            TLorentzVector ijet, jjet;
                            ijet.SetPtEtaPhiM(jet_pt[jet_iscleanH4l[i]],jet_eta[jet_iscleanH4l[i]],jet_phi[jet_iscleanH4l[i]],jet_mass[jet_iscleanH4l[i]]);
                            jjet.SetPtEtaPhiM(jet_pt[jet_iscleanH4l[j]],jet_eta[jet_iscleanH4l[j]],jet_phi[jet_iscleanH4l[j]],jet_mass[jet_iscleanH4l[j]]);
                            if (ijet.Pt()<40.0 || abs(ijet.Eta())>2.4) continue;
                            if (jjet.Pt()<40.0 || abs(jjet.Eta())>2.4) continue;
                            TLorentzVector Dijet;
                            Dijet = ijet+jjet;
                            double mass = Dijet.M();
                            if (mass > 60 && mass < 120) nvjets_pt40_eta2p4++;            
                        }
                    }
                    
                    int sumplus=0; int summinus=0;
                    int sumflavor=0;
                    int noverlapping=0;
                    for(unsigned int i = 0; i < lep_pt.size(); i++) {

                        bool goodi=false;
                        if ((int)i==lep_Hindex[0] || (int)i==lep_Hindex[1] || (int)i==lep_Hindex[2] || (int)i==lep_Hindex[3]) {
                            goodi=true;
                        } else {
                            if ( (abs(lep_id[i])==11 && lep_tightId[i]==1 && lep_RelIsoNoFSR[i]<isoCutEl) || 
                                 (abs(lep_id[i])==13 && lep_tightId[i]==1 && lep_RelIsoNoFSR[i]<isoCutMu) ) {
                                goodi=true;
                            }
                        }
                        if (!goodi) continue;
                                                
                        nisoleptons++;  

                        if (lep_id[i]>0) sumplus++;
                        if (lep_id[i]<0) summinus++;
                        sumflavor+=lep_id[i];

                        TLorentzVector lepi;
                        lepi.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

                        float minDr=9999.0;
                        for(unsigned int j = 0; j < lep_pt.size(); j++) {
                        
                            if (i==j) continue;

                            bool goodj=false;
                            if ((int)j==lep_Hindex[0] || (int)j==lep_Hindex[1] || (int)j==lep_Hindex[2] || (int)j==lep_Hindex[3]) {
                                goodj=true;
                            } else {
                                if ( (abs(lep_id[j])==11 && lep_tightId[j]==1 && lep_RelIsoNoFSR[j]<isoCutEl) || 
                                     (abs(lep_id[j])==13 && lep_tightId[j]==1 && lep_RelIsoNoFSR[j]<isoCutMu) ) {
                                    goodj=true;
                                }
                            }
                            if (!goodj) continue;
                           
                            TLorentzVector lepj;
                            lepj.SetPtEtaPhiM(lep_pt[j],lep_eta[j],lep_phi[j],lep_mass[j]);
                            float thisdR = lepi.DeltaR(lepj);
                            if (thisdR<minDr) minDr=thisdR;
                        }
                        if (minDr<0.02) noverlapping+=1;
                    }
                    nisoleptons-=(noverlapping/2);

                    // Event Categories 
                    if (nisoleptons==4 && (((njets_pt30_eta4p7==2||njets_pt30_eta4p7==3)&&nbjets_pt30_eta4p7<2)||(njets_pt30_eta4p7>=4&&nbjets_pt30_eta4p7==0)) && D_VBF>0.46386) {EventCat=2;}
                    else if (nisoleptons==4 && ((njets_pt30_eta4p7==2||njets_pt30_eta4p7==3)||(njets_pt30_eta4p7>=4&&nbjets_pt30_eta4p7==0)) && (D_HadWH>0.88384 || D_HadZH>0.91315)) {EventCat=4;}
                    else if (njets_pt30_eta4p7<=3 && nbjets_pt30_eta4p7==0 && (nisoleptons==5 || (nisoleptons>=6&&sumplus>=3&&summinus>=3&&sumflavor==0))) {EventCat=3;}
                    else if (njets_pt30_eta4p7==0 && nisoleptons>=5) {EventCat=3;}
                    else if (njets_pt30_eta4p7>=4 && nbjets_pt30_eta4p7>0 && nisoleptons==4) {EventCat=5;}
                    else if (nisoleptons>=5) {EventCat=6;}
                    else if (nisoleptons==4 && njets_pt30_eta4p7<=1 && met>100.0) {EventCat=7;}
                    else if (nisoleptons==4 && njets_pt30_eta4p7==1 && D_VBF1j>0.37605) {EventCat=1;}
                    else {EventCat=0;}
                    
                }
                    
                // fill the vector<float>
				eleUncorr_pt_float.assign(eleUncorr_pt.begin(),eleUncorr_pt.end());
				eleUncorr_mass_float.assign(eleUncorr_mass.begin(),eleUncorr_mass.end());
				
				track_vtx_dxy_float.assign(track_vtx_dxy.begin(),track_vtx_dxy.end());
				track_bs_dxy_float.assign(track_bs_dxy.begin(),track_bs_dxy.end());
				
				track_vtx_dxycorr_float.assign(track_vtx_dxycorr.begin(),track_vtx_dxycorr.end());
				track_vtx_pt_float.assign(track_vtx_pt.begin(),track_vtx_pt.end());
				track_vtx_eta_float.assign(track_vtx_eta.begin(),track_vtx_eta.end());
				track_vtx_phi_float.assign(track_vtx_phi.begin(),track_vtx_phi.end());
				track_vtx_pterr_float.assign(track_vtx_pterr.begin(),track_vtx_pterr.end());
				
				track_bs_dxycorr_float.assign(track_bs_dxycorr.begin(),track_bs_dxycorr.end());
				track_bs_pt_float.assign(track_bs_pt.begin(),track_bs_pt.end());
				track_bs_eta_float.assign(track_bs_eta.begin(),track_bs_eta.end());
				track_bs_phi_float.assign(track_bs_phi.begin(),track_bs_phi.end());
				track_bs_pterr_float.assign(track_bs_pterr.begin(),track_bs_pterr.end());
				
//				track_vtx_ptcorr_float.assign(track_vtx_ptcorr.begin(),track_vtx_ptcorr.end());
//				track_vtx_pterrcorr_float.assign(track_vtx_pterrcorr.begin(),track_vtx_pterrcorr.end());
//				track_bs_ptcorr_float.assign(track_bs_ptcorr.begin(),track_bs_pterrcorr.end());
//				track_bs_pterrcorr_float.assign(track_bs_pterrcorr.begin(),track_bs_pterrcorr.end());
				
//				track_vtx_FSRpt_float.assign(track_vtx_FSRpt.begin(),track_vtx_FSRpt.end());
//				track_vtx_FSReta_float.assign(track_vtx_FSReta.begin(),track_vtx_FSReta.end());
//				track_vtx_FSRphi_float.assign(track_vtx_FSRphi.begin(),track_vtx_FSRphi.end());
//				track_vtx_FSRmass_float.assign(track_vtx_FSRmass.begin(),track_vtx_FSRmass.end());
				
//				track_bs_FSRpt_float.assign(track_bs_FSRpt.begin(),track_bs_FSRpt.end());
//				track_bs_FSReta_float.assign(track_bs_FSReta.begin(),track_bs_FSReta.end());
//				track_bs_FSRphi_float.assign(track_bs_FSRphi.begin(),track_bs_FSRphi.end());
//				track_bs_FSRmass_float.assign(track_bs_FSRmass.begin(),track_bs_FSRmass.end());
				
				
				
				track_Hvtx_dxycorr_float.assign(track_Hvtx_dxycorr.begin(),track_Hvtx_dxycorr.end());
				track_Hvtx_pt_float.assign(track_Hvtx_pt.begin(),track_Hvtx_pt.end());
				track_Hvtx_eta_float.assign(track_Hvtx_eta.begin(),track_Hvtx_eta.end());
				track_Hvtx_phi_float.assign(track_Hvtx_phi.begin(),track_Hvtx_phi.end());
				track_Hvtx_pterr_float.assign(track_Hvtx_pterr.begin(),track_Hvtx_pterr.end());
				
				track_Hbs_dxycorr_float.assign(track_Hbs_dxycorr.begin(),track_Hbs_dxycorr.end());
				track_Hbs_pt_float.assign(track_Hbs_pt.begin(),track_Hbs_pt.end());
				track_Hbs_eta_float.assign(track_Hbs_eta.begin(),track_Hbs_eta.end());
				track_Hbs_phi_float.assign(track_Hbs_phi.begin(),track_Hbs_phi.end());
				track_Hbs_pterr_float.assign(track_Hbs_pterr.begin(),track_Hbs_pterr.end());
				
                lep_pt_float.assign(lep_pt.begin(),lep_pt.end());
                lep_pterr_float.assign(lep_pterr.begin(),lep_pterr.end());
                lep_pterrold_float.assign(lep_pterrold.begin(),lep_pterrold.end());
                lep_p_float.assign(lep_p.begin(),lep_p.end());
                lep_ecalEnergy_float.assign(lep_ecalEnergy.begin(),lep_ecalEnergy.end());
                lep_eta_float.assign(lep_eta.begin(),lep_eta.end());
                lep_phi_float.assign(lep_phi.begin(),lep_phi.end());
                lep_mass_float.assign(lep_mass.begin(),lep_mass.end());
                lepFSR_pt_float.assign(lepFSR_pt.begin(),lepFSR_pt.end());
                lepFSR_eta_float.assign(lepFSR_eta.begin(),lepFSR_eta.end());
                lepFSR_phi_float.assign(lepFSR_phi.begin(),lepFSR_phi.end());
                lepFSR_mass_float.assign(lepFSR_mass.begin(),lepFSR_mass.end());
                tau_pt_float.assign(tau_pt.begin(),tau_pt.end());
                tau_eta_float.assign(tau_eta.begin(),tau_eta.end());
                tau_phi_float.assign(tau_phi.begin(),tau_phi.end());
                tau_mass_float.assign(tau_mass.begin(),tau_mass.end());
                pho_pt_float.assign(pho_pt.begin(),pho_pt.end());
                pho_eta_float.assign(pho_eta.begin(),pho_eta.end());
                pho_phi_float.assign(pho_phi.begin(),pho_phi.end());
                photonCutBasedIDLoose_float.assign(photonCutBasedIDLoose.begin(),photonCutBasedIDLoose.end());
                H_pt_float.assign(H_pt.begin(),H_pt.end());
                H_eta_float.assign(H_eta.begin(),H_eta.end());
                H_phi_float.assign(H_phi.begin(),H_phi.end());
                H_mass_float.assign(H_mass.begin(),H_mass.end());
                H_noFSR_pt_float.assign(H_noFSR_pt.begin(),H_noFSR_pt.end());
                H_noFSR_eta_float.assign(H_noFSR_eta.begin(),H_noFSR_eta.end());
                H_noFSR_phi_float.assign(H_noFSR_phi.begin(),H_noFSR_phi.end());
                H_noFSR_mass_float.assign(H_noFSR_mass.begin(),H_noFSR_mass.end());
                Z_pt_float.assign(Z_pt.begin(),Z_pt.end());
                Z_eta_float.assign(Z_eta.begin(),Z_eta.end());
                Z_phi_float.assign(Z_phi.begin(),Z_phi.end());
                Z_mass_float.assign(Z_mass.begin(),Z_mass.end());
                Z_noFSR_pt_float.assign(Z_noFSR_pt.begin(),Z_noFSR_pt.end());
                Z_noFSR_eta_float.assign(Z_noFSR_eta.begin(),Z_noFSR_eta.end());
                Z_noFSR_phi_float.assign(Z_noFSR_phi.begin(),Z_noFSR_phi.end());
                Z_noFSR_mass_float.assign(Z_noFSR_mass.begin(),Z_noFSR_mass.end());
                jet_pt_float.assign(jet_pt.begin(),jet_pt.end());
                jet_pt_raw_float.assign(jet_pt_raw.begin(),jet_pt_raw.end());
                jet_eta_float.assign(jet_eta.begin(),jet_eta.end());
                jet_phi_float.assign(jet_phi.begin(),jet_phi.end());
                jet_mass_float.assign(jet_mass.begin(),jet_mass.end());
                jet_jesup_pt_float.assign(jet_jesup_pt.begin(),jet_jesup_pt.end());
                jet_jesup_eta_float.assign(jet_jesup_eta.begin(),jet_jesup_eta.end());
                jet_jesup_phi_float.assign(jet_jesup_phi.begin(),jet_jesup_phi.end());
                jet_jesup_mass_float.assign(jet_jesup_mass.begin(),jet_jesup_mass.end());
                jet_jesdn_pt_float.assign(jet_jesdn_pt.begin(),jet_jesdn_pt.end());
                jet_jesdn_eta_float.assign(jet_jesdn_eta.begin(),jet_jesdn_eta.end());
                jet_jesdn_phi_float.assign(jet_jesdn_phi.begin(),jet_jesdn_phi.end());
                jet_jesdn_mass_float.assign(jet_jesdn_mass.begin(),jet_jesdn_mass.end());
                jet_jerup_pt_float.assign(jet_jerup_pt.begin(),jet_jerup_pt.end());
                jet_jerup_eta_float.assign(jet_jerup_eta.begin(),jet_jerup_eta.end());
                jet_jerup_phi_float.assign(jet_jerup_phi.begin(),jet_jerup_phi.end());
                jet_jerup_mass_float.assign(jet_jerup_mass.begin(),jet_jerup_mass.end());
                jet_jerdn_pt_float.assign(jet_jerdn_pt.begin(),jet_jerdn_pt.end());
                jet_jerdn_eta_float.assign(jet_jerdn_eta.begin(),jet_jerdn_eta.end());
                jet_jerdn_phi_float.assign(jet_jerdn_phi.begin(),jet_jerdn_phi.end());
                jet_jerdn_mass_float.assign(jet_jerdn_mass.begin(),jet_jerdn_mass.end());
                fsrPhotons_pt_float.assign(fsrPhotons_pt.begin(),fsrPhotons_pt.end());
                fsrPhotons_pterr_float.assign(fsrPhotons_pterr.begin(),fsrPhotons_pterr.end());
                fsrPhotons_eta_float.assign(fsrPhotons_eta.begin(),fsrPhotons_eta.end());
                fsrPhotons_phi_float.assign(fsrPhotons_phi.begin(),fsrPhotons_phi.end());
                fsrPhotons_mass_float.assign(fsrPhotons_mass.begin(),fsrPhotons_mass.end());
                
                if (!isMC) passedEventsTree_All->Fill();        
            } // 2 tight ID
            else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed  ntight ID"<<endl;}
        } //if 2 lepID
        else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed  nloose ID"<<endl;}
    } //primary vertex,notDuplicate
    else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed primary vertex"<<endl;}
    GENlep_pt_float.clear(); GENlep_pt_float.assign(GENlep_pt.begin(),GENlep_pt.end());
    GENlep_eta_float.clear(); GENlep_eta_float.assign(GENlep_eta.begin(),GENlep_eta.end());
    GENlep_phi_float.clear(); GENlep_phi_float.assign(GENlep_phi.begin(),GENlep_phi.end());
    GENlep_mass_float.clear(); GENlep_mass_float.assign(GENlep_mass.begin(),GENlep_mass.end());
    GENH_pt_float.clear(); GENH_pt_float.assign(GENH_pt.begin(),GENH_pt.end());
    GENH_eta_float.clear(); GENH_eta_float.assign(GENH_eta.begin(),GENH_eta.end());
    GENH_phi_float.clear(); GENH_phi_float.assign(GENH_phi.begin(),GENH_phi.end());
    GENH_mass_float.clear(); GENH_mass_float.assign(GENH_mass.begin(),GENH_mass.end());
    GENZ_pt_float.clear(); GENZ_pt_float.assign(GENZ_pt.begin(),GENZ_pt.end());
    GENZ_eta_float.clear(); GENZ_eta_float.assign(GENZ_eta.begin(),GENZ_eta.end());
    GENZ_phi_float.clear(); GENZ_phi_float.assign(GENZ_phi.begin(),GENZ_phi.end());
    GENZ_mass_float.clear(); GENZ_mass_float.assign(GENZ_mass.begin(),GENZ_mass.end());
    GENjet_pt_float.clear(); GENjet_pt_float.assign(GENjet_pt.begin(),GENjet_pt.end());
    GENjet_eta_float.clear(); GENjet_eta_float.assign(GENjet_eta.begin(),GENjet_eta.end());
    GENjet_phi_float.clear(); GENjet_phi_float.assign(GENjet_phi.begin(),GENjet_phi.end());
    GENjet_mass_float.clear(); GENjet_mass_float.assign(GENjet_mass.begin(),GENjet_mass.end());

    if (isMC) passedEventsTree_All->Fill();
    
    if (nEventsTotal==1000.0) passedEventsTree_All->OptimizeBaskets();

}



// ------------ method called once each job just before starting event loop  ------------
void 
UFHZZ4LAna::beginJob()
{
    using namespace edm;
    using namespace std;
    using namespace pat;

    bookPassedEventTree("passedEvents", passedEventsTree_All);

    firstEntry = true;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
UFHZZ4LAna::endJob() 
{
    histContainer_["NEVENTS"]->SetBinContent(1,nEventsTotal);
    histContainer_["NEVENTS"]->GetXaxis()->SetBinLabel(1,"N Events in Sample");
    histContainer_["SUMWEIGHTS"]->SetBinContent(1,sumWeightsTotal);
    histContainer_["SUMWEIGHTSPU"]->SetBinContent(1,sumWeightsTotalPU);
    histContainer_["SUMWEIGHTS"]->GetXaxis()->SetBinLabel(1,"sum Weights in Sample");
    histContainer_["SUMWEIGHTSPU"]->GetXaxis()->SetBinLabel(1,"sum Weights PU in Sample");
}

void
UFHZZ4LAna::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{

    //massErr.init(iSetup);
    if (isMC) {
        edm::Handle<LHERunInfoProduct> run;
        typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
        try {

            int pos=0;
            iRun.getByLabel( edm::InputTag("externalLHEProducer"), run );
            LHERunInfoProduct myLHERunInfoProduct = *(run.product());
            typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
            for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
                std::cout << iter->tag() << std::endl;
                std::vector<std::string> lines = iter->lines();
                for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
                    std::string pdfid=lines.at(iLine);
                    if (pdfid.substr(1,6)=="weight" && pdfid.substr(8,2)=="id") {
                        std::cout<<pdfid<<std::endl;
                        std::string pdf_weight_id = pdfid.substr(12,4);
                        int pdf_weightid=atoi(pdf_weight_id.c_str());
                        std::cout<<"parsed id: "<<pdf_weightid<<std::endl;
                        if (pdf_weightid==2001) {posNNPDF=int(pos);}
                        pos+=1;
                    }
                }
            }
        }
        catch(...) {
            std::cout<<"No LHERunInfoProduct"<<std::endl;
        }
    }

}


// ------------ method called when ending the processing of a run  ------------
void 
UFHZZ4LAna::endRun(const edm::Run& iRun, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
UFHZZ4LAna::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
UFHZZ4LAna::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg,edm::EventSetup const& eSetup)
{
    using namespace edm;
    using namespace std;
    // Keep track of all the events run over
    edm::Handle<MergeableCounter> numEventsCounter;
    lumiSeg.getByLabel("nEventsTotal", numEventsCounter);    
    if(numEventsCounter.isValid()) {
        std::cout<<"numEventsCounter->value "<<numEventsCounter->value<<endl;
        nEventsTotal += numEventsCounter->value;        
    }
}

// ============================ UF Functions ============================= //


//Find Z1,Z2, and Higgs candidate
//Pass good leptons for analysis as candMuons and candElectrons
//Pass empty vectors of pat leptons as selectedMuons and selectedElectrons
// these will be filled in the function and then useable for more analysis.
void
UFHZZ4LAna::findHiggsCandidate(std::vector< pat::Muon > &selectedMuons, std::vector< pat::Electron > &selectedElectrons,const edm::Event& iEvent )
{

    using namespace edm;
    using namespace pat;
    using namespace std;

    const double Zmass = 91.1876;

    unsigned int Nlep = lepFSR_pt.size();
    if (verbose) cout<<Nlep<<" leptons in total"<<endl;

    // First, make all Z candidates including any FSR photons
    int n_Zs=0;
    vector<int> Z_lepindex1;
    vector<int> Z_lepindex2;

    for(unsigned int i=0; i<Nlep; i++){
        //if (!(lep_tightId[i])) continue; // checking tight lepton ID
        //if (lep_RelIsoNoFSR[i]>((abs(lep_id[i])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
        for(unsigned int j=i+1; j<Nlep; j++){
            //if (!(lep_tightId[j])) continue; // checking tight lepton ID
            //if (lep_RelIsoNoFSR[j]>((abs(lep_id[j])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed

            // same flavor opposite charge
            if((lep_id[i]+lep_id[j])!=0) continue;

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
            lj.SetPtEtaPhiM(lep_pt[j],lep_eta[j],lep_phi[j],lep_mass[j]);

            TLorentzVector lifsr, ljfsr;
            lifsr.SetPtEtaPhiM(lepFSR_pt[i],lepFSR_eta[i],lepFSR_phi[i],lepFSR_mass[i]);
            ljfsr.SetPtEtaPhiM(lepFSR_pt[j],lepFSR_eta[j],lepFSR_phi[j],lepFSR_mass[j]);

            TLorentzVector liljfsr = lifsr+ljfsr;

            if (verbose) {
                cout<<"OSSF pair: i="<<i<<" id1="<<lep_id[i]<<" j="<<j<<" id2="<<lep_id[j]<<" pt1: "<<lifsr.Pt()<<" pt2: "<<ljfsr.Pt()<<" M: "<<liljfsr.M()<<endl;    
            }
            
            TLorentzVector Z, Z_noFSR;
            Z = lifsr+ljfsr;
            Z_noFSR = li+lj;
            
            if (verbose) cout<<"this Z mass: "<<Z.M()<<" mZ2Low: "<<mZ2Low<<endl;

            if (Z.M()>0.0) {
                n_Zs++;
                Z_pt.push_back(Z.Pt());
                Z_eta.push_back(Z.Eta());
                Z_phi.push_back(Z.Phi());
                Z_mass.push_back(Z.M());
                Z_noFSR_pt.push_back(Z_noFSR.Pt());
                Z_noFSR_eta.push_back(Z_noFSR.Eta());
                Z_noFSR_phi.push_back(Z_noFSR.Phi());
                Z_noFSR_mass.push_back(Z_noFSR.M());
                Z_lepindex1.push_back(i);
                Z_lepindex2.push_back(j);
                if (verbose) cout<<" add Z_lepindex1: "<<i<<" Z_lepindex2: "<<j<<endl;
            }

        } // lep i
    } // lep j
    
    if( (recoMuons.size() + recoElectrons.size()) < 4 ) return;
        
    if (verbose) cout<<"found four leptons"<<endl;     

    bool properLep_ID = false; int Nmm = 0; int Nmp = 0; int Nem = 0; int Nep = 0;
    for(unsigned int i =0; i<recoMuons.size(); i++) {
        if(recoMuons[i].pdgId()<0) Nmm = Nmm+1;
        if(recoMuons[i].pdgId()>0) Nmp = Nmp+1;
    }
    for(unsigned int i =0; i<recoElectrons.size(); i++) {
        if(recoElectrons[i].pdgId()<0) Nem = Nem+1;
        if(recoElectrons[i].pdgId()>0) Nep = Nep+1;
    }
    
    if(Nmm>=2 && Nmp>=2) properLep_ID = true; //4mu
    if(Nem>=2 && Nep>=2) properLep_ID = true; //4e
    if(Nmm>0 && Nmp>0 && Nem>0 && Nep>0) properLep_ID = true; //2e2mu
    
    // four proper charge flavor combination
    if(!properLep_ID) return;

    // Consider all ZZ candidates
    double minZ1DeltaM_SR=9999.9; double minZ1DeltaM_CR=99999.9;
    double maxZ2SumPt_SR=0.0; double maxZ2SumPt_CR=0.0;
    double max_D_bkg_kin_SR=0.0; double max_D_bkg_kin_CR=0.0;
    bool foundSRCandidate=false;

    for (int i=0; i<n_Zs; i++) {
        for (int j=i+1; j<n_Zs; j++) {
 
            int i1 = Z_lepindex1[i]; int i2 = Z_lepindex2[i];                            
            int j1 = Z_lepindex1[j]; int j2 = Z_lepindex2[j];                            

            if (i1 == j1 || i1 == j2 || i2 == j1 || i2 == j2) continue;

            TLorentzVector lep_i1, lep_i2, lep_j1, lep_j2;
            lep_i1.SetPtEtaPhiM(lepFSR_pt[i1],lepFSR_eta[i1],lepFSR_phi[i1],lepFSR_mass[i1]);
            lep_i2.SetPtEtaPhiM(lepFSR_pt[i2],lepFSR_eta[i2],lepFSR_phi[i2],lepFSR_mass[i2]);
            lep_j1.SetPtEtaPhiM(lepFSR_pt[j1],lepFSR_eta[j1],lepFSR_phi[j1],lepFSR_mass[j1]);
            lep_j2.SetPtEtaPhiM(lepFSR_pt[j2],lepFSR_eta[j2],lepFSR_phi[j2],lepFSR_mass[j2]);

            TLorentzVector lep_i1_nofsr, lep_i2_nofsr, lep_j1_nofsr, lep_j2_nofsr;
            lep_i1_nofsr.SetPtEtaPhiM(lep_pt[i1],lep_eta[i1],lep_phi[i1],lep_mass[i1]);
            lep_i2_nofsr.SetPtEtaPhiM(lep_pt[i2],lep_eta[i2],lep_phi[i2],lep_mass[i2]);
            lep_j1_nofsr.SetPtEtaPhiM(lep_pt[j1],lep_eta[j1],lep_phi[j1],lep_mass[j1]);
            lep_j2_nofsr.SetPtEtaPhiM(lep_pt[j2],lep_eta[j2],lep_phi[j2],lep_mass[j2]);

            TLorentzVector Zi, Zj;
            Zi.SetPtEtaPhiM(Z_pt[i],Z_eta[i],Z_phi[i],Z_mass[i]);
            Zj.SetPtEtaPhiM(Z_pt[j],Z_eta[j],Z_phi[j],Z_mass[j]);
            
            if (verbose) {cout<<"ZZ candidate Zi->M() "<<Zi.M()<<" Zj->M() "<<Zj.M()<<endl;}

            TLorentzVector Z1, Z2;
            int Z1index, Z2index;
            int Z1_lepindex[2] = {0,0};
            int Z2_lepindex[2] = {0,0};
            double Z1DeltaM, Z2SumPt;

            if (abs(Zi.M()-Zmass)<abs(Zj.M()-Zmass)) { 
                Z1index = i; Z2index = j;
                Z1 = Zi; Z2 = Zj;                 
                if (lep_i1.Pt()>lep_i2.Pt()) { Z1_lepindex[0] = i1;  Z1_lepindex[1] = i2; }
                else { Z1_lepindex[0] = i2;  Z1_lepindex[1] = i1; }                
                if (lep_j1.Pt()>lep_j2.Pt()) { Z2_lepindex[0] = j1;  Z2_lepindex[1] = j2; } 
                else { Z2_lepindex[0] = j2;  Z2_lepindex[1] = j1; }                
                Z1DeltaM = abs(Zi.M()-Zmass); 
                Z2SumPt = lep_j1_nofsr.Pt()+lep_j2_nofsr.Pt();
            }
            else { 
                Z1index = j; Z2index = i;
                Z1 = Zj; Z2 = Zi; 
                if (lep_j1.Pt()>lep_j2.Pt()) { Z1_lepindex[0] = j1;  Z1_lepindex[1] = j2; }
                else { Z1_lepindex[0] = j2;  Z1_lepindex[1] = j1; }
                if (lep_i1.Pt()>lep_i2.Pt()) { Z2_lepindex[0] = i1;  Z2_lepindex[1] = i2; }
                else { Z2_lepindex[0] = i2;  Z2_lepindex[1] = i1; }
                Z1DeltaM = abs(Zj.M()-Zmass); 
                Z2SumPt = lep_i1_nofsr.Pt()+lep_i2_nofsr.Pt();
            }

            // Check isolation cut (without FSR ) for Z1 leptons
            if (lep_RelIsoNoFSR[Z1_lepindex[0]]>((abs(lep_id[Z1_lepindex[0]])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
            if (lep_RelIsoNoFSR[Z1_lepindex[1]]>((abs(lep_id[Z1_lepindex[1]])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
            // Check tight ID cut for Z1 leptons
            if (!(lep_tightId[Z1_lepindex[0]])) continue; // checking tight lepton ID
            if (!(lep_tightId[Z1_lepindex[1]])) continue; // checking tight lepton ID 
           
            // Check Leading and Subleading pt Cut
            vector<double> allPt;
            allPt.push_back(lep_i1_nofsr.Pt()); allPt.push_back(lep_i2_nofsr.Pt());
            allPt.push_back(lep_j1_nofsr.Pt()); allPt.push_back(lep_j2_nofsr.Pt());
            std::sort(allPt.begin(), allPt.end());
            if (verbose) cout<<" leading pt: "<<allPt[3]<<" cut: "<<leadingPtCut<<" subleadingPt: "<<allPt[2]<<" cut: "<<subleadingPtCut<<endl;
            if (allPt[3]<leadingPtCut || allPt[2]<subleadingPtCut ) continue;
            
            // Check dR(li,lj)>0.02 for any i,j
            vector<double> alldR;
            alldR.push_back(deltaR(lep_i1_nofsr.Eta(),lep_i1_nofsr.Phi(),lep_i2_nofsr.Eta(),lep_i2_nofsr.Phi()));
            alldR.push_back(deltaR(lep_i1_nofsr.Eta(),lep_i1_nofsr.Phi(),lep_j1_nofsr.Eta(),lep_j1_nofsr.Phi()));
            alldR.push_back(deltaR(lep_i1_nofsr.Eta(),lep_i1_nofsr.Phi(),lep_j2_nofsr.Eta(),lep_j2_nofsr.Phi()));
            alldR.push_back(deltaR(lep_i2_nofsr.Eta(),lep_i2_nofsr.Phi(),lep_j1_nofsr.Eta(),lep_j1_nofsr.Phi()));
            alldR.push_back(deltaR(lep_i2_nofsr.Eta(),lep_i2_nofsr.Phi(),lep_j2_nofsr.Eta(),lep_j2_nofsr.Phi()));
            alldR.push_back(deltaR(lep_j1_nofsr.Eta(),lep_j1_nofsr.Phi(),lep_j2_nofsr.Eta(),lep_j2_nofsr.Phi()));            
            if (verbose) cout<<" minDr: "<<*min_element(alldR.begin(),alldR.end())<<endl;
            if (*min_element(alldR.begin(),alldR.end())<0.02) continue;

            // Check M(l+,l-)>4.0 GeV for any OS pair
            // Do not include FSR photons
            vector<double> allM;
            TLorentzVector i1i2;
            TLorentzVector _4l_temp;
            i1i2 = (lep_i1_nofsr)+(lep_i2_nofsr); allM.push_back(i1i2.M());
            TLorentzVector j1j2;
            j1j2 = (lep_j1_nofsr)+(lep_j2_nofsr); allM.push_back(j1j2.M());            
             _4l_temp = Z1 + Z2;

            if (lep_id[i1]*lep_id[j1]<0) {
                TLorentzVector i1j1;
                i1j1 = (lep_i1_nofsr)+(lep_j1_nofsr); allM.push_back(i1j1.M());
                TLorentzVector i2j2;
                i2j2 = (lep_i2_nofsr)+(lep_j2_nofsr); allM.push_back(i2j2.M());
            } else {
                TLorentzVector i1j2;
                i1j2 = (lep_i1_nofsr)+(lep_j2_nofsr); allM.push_back(i1j2.M());
                TLorentzVector i2j1;
                i2j1 = (lep_i2_nofsr)+(lep_j1_nofsr); allM.push_back(i2j1.M());
            }
            if (verbose) cout<<" min m(l+l-): "<<*min_element(allM.begin(),allM.end())<<endl;
            if (*min_element(allM.begin(),allM.end())<4.0) {passedQCDcut=false; continue;}

            // Check the "smart cut": !( |mZa-mZ| < |mZ1-mZ| && mZb<12)
            // only for 4mu or 4e ZZ candidates
            bool passSmartCut=true;
            if ( abs(lep_id[i1])==abs(lep_id[j1])) {
                TLorentzVector Za, Zb;
                if (lep_id[i1]==lep_id[j1]) {                  
                    Za = (lep_i1)+(lep_j2);
                    Zb = (lep_i2)+(lep_j1);                    
                } else {
                    Za = (lep_i1)+(lep_j1);
                    Zb = (lep_i2)+(lep_j2);
                }                
                if ( abs(Za.M()-Zmass)<abs(Zb.M()-Zmass) ) {
                    if (verbose) cout<<"abs(Za.M()-Zmass)-abs(Z1.M()-Zmass): "<<abs(Za.M()-Zmass)-abs(Z1.M()-Zmass)<<" Zb.M(): "<<Zb.M()<<endl;
                    if ( abs(Za.M()-Zmass)<abs(Z1.M()-Zmass) && Zb.M()<mZ2Low ) passSmartCut=false;
                }
                else {
                    if (verbose) cout<<"abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass): "<<abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass)<<" Za.M(): "<<Za.M()<<endl;
                    if ( abs(Zb.M()-Zmass)<abs(Z1.M()-Zmass) && Za.M()<mZ2Low ) passSmartCut=false;
                }

            }
            if (!passSmartCut) continue;
            
            if (verbose) cout<<" massZ1: "<<Z1.M()<<" massZ2: "<<Z2.M()<<endl;
            if ( (Z1.M() < mZ1Low) || (Z1.M() > mZ1High) || (Z2.M() < mZ2Low) || (Z2.M() > mZ2High) ) continue;
            if (verbose) cout<<" mass4l: "<<_4l_temp.M()<<endl;
            if ( _4l_temp.M() < m4lLowCut ) continue;

            // Signal region if Z2 leptons are both tight ID Iso
            bool signalRegion=true;
            if (lep_RelIsoNoFSR[Z2_lepindex[0]]>((abs(lep_id[Z2_lepindex[0]])==11) ? isoCutEl : isoCutMu)) signalRegion=false; // checking iso with FSR removed
            if (lep_RelIsoNoFSR[Z2_lepindex[1]]>((abs(lep_id[Z2_lepindex[1]])==11) ? isoCutEl : isoCutMu)) signalRegion=false; // checking iso with FSR removed
            if (!(lep_tightId[Z2_lepindex[0]])) signalRegion=false; // checking tight lepton ID
            if (!(lep_tightId[Z2_lepindex[1]])) signalRegion=false; // checking tight lepton ID          
            
            // Check if this candidate has the highest D_bkg_kin
            vector<TLorentzVector> P4s;
            P4s.clear();
            vector<int> tmpIDs;
            tmpIDs.clear();

            if (Z1_lepindex[0] == i1) {
                P4s.push_back(lep_i1); P4s.push_back(lep_i2);
                if (Z2_lepindex[0] == j1) {
                    P4s.push_back(lep_j1); P4s.push_back(lep_j2);
                } else {
                    P4s.push_back(lep_j2); P4s.push_back(lep_j1);
                }
            } else if (Z1_lepindex[0] == i2) {
                P4s.push_back(lep_i2); P4s.push_back(lep_i1);
                if (Z2_lepindex[0] == j1) {
                    P4s.push_back(lep_j1); P4s.push_back(lep_j2);
                } else {
                    P4s.push_back(lep_j2); P4s.push_back(lep_j1);
                }
            } else if (Z1_lepindex[0] == j1) {
                P4s.push_back(lep_j1); P4s.push_back(lep_j2);
                if (Z2_lepindex[0] == i1) {
                    P4s.push_back(lep_i1); P4s.push_back(lep_i2);
                } else {
                    P4s.push_back(lep_i2); P4s.push_back(lep_i1);
                }
            } else if (Z1_lepindex[0] == j2) {
                P4s.push_back(lep_j2); P4s.push_back(lep_j1);
                if (Z2_lepindex[0] == i1) {
                    P4s.push_back(lep_i1); P4s.push_back(lep_i2);
                } else {
                    P4s.push_back(lep_i2); P4s.push_back(lep_i1);
                }
            }

            tmpIDs.push_back(lep_id[Z1_lepindex[0]]); tmpIDs.push_back(lep_id[Z1_lepindex[1]]);
            tmpIDs.push_back(lep_id[Z2_lepindex[0]]); tmpIDs.push_back(lep_id[Z2_lepindex[1]]);

            SimpleParticleCollection_t daughters;
            daughters.push_back(SimpleParticle_t(tmpIDs[0],P4s[0]));
            daughters.push_back(SimpleParticle_t(tmpIDs[1],P4s[1]));
            daughters.push_back(SimpleParticle_t(tmpIDs[2],P4s[2]));
            daughters.push_back(SimpleParticle_t(tmpIDs[3],P4s[3]));

            SimpleParticleCollection_t associated;
	    float D_bkg_kin_tmp; 
	if(doMela){
            mela->setInputEvent(&daughters, &associated, 0, 0);
            mela->setCurrentCandidateFromIndex(0);

            float me_0plus_JHU_tmp, me_qqZZ_MCFM_tmp;
            mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
            mela->computeP(me_0plus_JHU_tmp, true);            
            mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
            mela->computeP(me_qqZZ_MCFM_tmp, true);
            D_bkg_kin_tmp = me_0plus_JHU_tmp / (me_0plus_JHU_tmp + me_qqZZ_MCFM_tmp);

            mela->resetInputEvent(); 
	 }
            if (verbose) cout<<"good ZZ candidate, D_bkg_kin: "<<D_bkg_kin_tmp<<" max D_bkg_kin SR: "<<max_D_bkg_kin_SR<<" max D_bkg_kin CR: "<<max_D_bkg_kin_CR<<endl;
            
            bool same4l=false;
            bool foundZ11=false; bool foundZ12=false; bool foundZ21=false; bool foundZ22=false;
            for(int l = 0; l < 4; l++){
                if (lep_Hindex[l]==Z1_lepindex[0]) foundZ11 = true;
                if (lep_Hindex[l]==Z1_lepindex[1]) foundZ12 = true;
                if (lep_Hindex[l]==Z2_lepindex[0]) foundZ21 = true;
                if (lep_Hindex[l]==Z2_lepindex[1]) foundZ22 = true;
            }
            same4l = (foundZ11 && foundZ12 && foundZ21 && foundZ22);

            if (signalRegion) { // Signal Region has priority
                
                if (!foundSRCandidate) same4l=false;                
                if ( (bestCandMela && ((!same4l && D_bkg_kin_tmp>max_D_bkg_kin_SR) || (same4l && Z1DeltaM<=minZ1DeltaM_SR))) 
                     || (!bestCandMela && Z1DeltaM<=minZ1DeltaM_SR) ) {                 
                //if ( (!same4l && D_bkg_kin_tmp>max_D_bkg_kin_SR) || (same4l && Z1DeltaM<=minZ1DeltaM_SR) ) {
                    
                    max_D_bkg_kin_SR = D_bkg_kin_tmp;
                    minZ1DeltaM_SR = Z1DeltaM;
                    
                    if (!bestCandMela && Z_Hindex[0]==Z1index && Z2SumPt<maxZ2SumPt_SR) continue;
                    
                    Z_Hindex[0] = Z1index;
                    lep_Hindex[0] = Z1_lepindex[0];
                    lep_Hindex[1] = Z1_lepindex[1];
                    
                    maxZ2SumPt_SR = Z2SumPt;
                    Z_Hindex[1] = Z2index;
                    lep_Hindex[2] = Z2_lepindex[0];
                    lep_Hindex[3] = Z2_lepindex[1];
                    
                    Z1Vec = Z1; Z2Vec = Z2; HVec = Z1+Z2;                   
                    massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); mass4l = HVec.M();
                    
                    if (verbose) cout<<" new best candidate SR: mass4l: "<<HVec.M()<<endl;
                    if (HVec.M()>m4lLowCut)  {  //m4lLowCut move forward
                        foundHiggsCandidate=true;                    
                        foundSRCandidate=true;
                    }
                }                

            } else if (!foundSRCandidate) { // Control regions get second priority

                if ( (bestCandMela && ((!same4l && D_bkg_kin_tmp>max_D_bkg_kin_CR) || (same4l && Z1DeltaM<=minZ1DeltaM_CR)))
                     || (!bestCandMela && Z1DeltaM<=minZ1DeltaM_CR) ) {                 
                //if ( (!same4l && D_bkg_kin_tmp>max_D_bkg_kin_CR) || (same4l && Z1DeltaM<=minZ1DeltaM_CR) ) {

                    max_D_bkg_kin_CR = D_bkg_kin_tmp;
                    minZ1DeltaM_CR = Z1DeltaM;
                    
                    if (!bestCandMela && Z_Hindex[0]==Z1index && Z2SumPt<maxZ2SumPt_CR) continue;

                    Z_Hindex[0] = Z1index;
                    lep_Hindex[0] = Z1_lepindex[0];
                    lep_Hindex[1] = Z1_lepindex[1];
                    
                    maxZ2SumPt_CR = Z2SumPt;
                    Z_Hindex[1] = Z2index;
                    lep_Hindex[2] = Z2_lepindex[0];
                    lep_Hindex[3] = Z2_lepindex[1];

                    Z1Vec = Z1; Z2Vec = Z2; HVec = Z1+Z2;                   
                    massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); mass4l = HVec.M();

                    if (verbose) cout<<" new best candidate CR: mass4l: "<<HVec.M()<<endl;
                    if (HVec.M()>m4lLowCut) foundHiggsCandidate=true;  //m4lLowCut move forward                  
                }
            }

            if (verbose) cout<<"Z_Hindex[0]: "<<Z_Hindex[0]<<" lep_Hindex[0]: "<<lep_Hindex[0]<<" lep_Hindex[1]: "<<lep_Hindex[1]
                             <<"Z_Hindex[1]: "<<Z_Hindex[1]<<" lep_Hindex[2]: "<<lep_Hindex[2]<<" lep_Hindex[3]: "<<lep_Hindex[3]<<endl;

        } // Zj
    } // Zi


    if(foundHiggsCandidate) {
        
        if (verbose) cout<<" lep_Hindex[0]: "<<lep_Hindex[0]<<" lep_Hindex[1]: "<<lep_Hindex[1]<<" lep_Hindex[2]: "<<lep_Hindex[2]<<" lep_Hindex[3]: "<<lep_Hindex[3]<<endl;

        if (verbose) cout<<" lep_id[lep_Hindex[0]]: "<<lep_id[lep_Hindex[0]]<<" lep_id[lep_Hindex[1]]: "<<lep_id[lep_Hindex[1]]
                         <<" lep_id[lep_Hindex[2]]: "<<lep_id[lep_Hindex[2]]<<" lep_id[lep_Hindex[3]]: "<<lep_id[lep_Hindex[3]]<<endl;

        if ( abs(lep_id[lep_Hindex[0]])==13 && abs(lep_id[lep_Hindex[2]])==13 ) {
            RecoFourMuEvent = true;
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[0]]]);
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[1]]]);
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[2]]]);
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[3]]]);
        }
        else if ( abs(lep_id[lep_Hindex[0]])==13 && abs(lep_id[lep_Hindex[2]])==11 ) {
            RecoTwoMuTwoEEvent = true;
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[0]]]);
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[1]]]);
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[2]]]);
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[3]]]);
        }
        else if ( abs(lep_id[lep_Hindex[0]])==11 && abs(lep_id[lep_Hindex[2]])==13 ) {
            RecoTwoETwoMuEvent = true;
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[0]]]);
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[1]]]);
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[2]]]);
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[3]]]);
        }
        else if ( abs(lep_id[lep_Hindex[0]])==11 && abs(lep_id[lep_Hindex[2]])==11 ) {
            RecoFourEEvent = true;
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[0]]]);
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[1]]]);
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[2]]]);
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[3]]]);
        }


    }

}

//Find Z + 1L candidate for fake rate study
void
UFHZZ4LAna::findZ1LCandidate(const edm::Event& iEvent )
{

    using namespace edm;
    using namespace pat;
    using namespace std;

    const double Zmass = 91.1876;

    unsigned int Nlep = lepFSR_pt.size();
    if (verbose) cout<<Nlep<<" leptons in total"<<endl;
    if( Nlep != 3 ) return;

    // First, make all Z candidates including any FSR photons
    int n_Zs=0;
    vector<int> Z_Z1L_lepindex1;
    vector<int> Z_Z1L_lepindex2;

    for(unsigned int i=0; i<Nlep; i++){
        for(unsigned int j=i+1; j<Nlep; j++){

            // same flavor opposite charge
            if((lep_id[i]+lep_id[j])!=0) continue;

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
            lj.SetPtEtaPhiM(lep_pt[j],lep_eta[j],lep_phi[j],lep_mass[j]);

            TLorentzVector lifsr, ljfsr;
            lifsr.SetPtEtaPhiM(lepFSR_pt[i],lepFSR_eta[i],lepFSR_phi[i],lepFSR_mass[i]);
            ljfsr.SetPtEtaPhiM(lepFSR_pt[j],lepFSR_eta[j],lepFSR_phi[j],lepFSR_mass[j]);

            TLorentzVector liljfsr = lifsr+ljfsr;

            if (verbose) {
                cout<<"OSSF pair: i="<<i<<" id1="<<lep_id[i]<<" j="<<j<<" id2="<<lep_id[j]<<" pt1: "<<lifsr.Pt()<<" pt2: "<<ljfsr.Pt()<<" M: "<<liljfsr.M()<<endl;
            }

            TLorentzVector Z, Z_noFSR;
            Z = lifsr+ljfsr;
            Z_noFSR = li+lj;

            if (verbose) cout<<"this Z mass: "<<Z.M()<<" mZ2Low: "<<mZ2Low<<endl;

            if (Z.M()>0.0) {
                n_Zs++;
                Z_Z1L_lepindex1.push_back(i);
                Z_Z1L_lepindex2.push_back(j);
                if (verbose) cout<<" add Z_lepindex1: "<<i<<" Z_lepindex2: "<<j<<endl;
            }

        } // lep i
    } // lep j

    bool properLep_ID = false; int Nmm = 0; int Nmp = 0; int Nem = 0; int Nep = 0;
    for(unsigned int i =0; i<recoMuons.size(); i++) {
        if(recoMuons[i].pdgId()<0) Nmm = Nmm+1;
        if(recoMuons[i].pdgId()>0) Nmp = Nmp+1;
    }
    for(unsigned int i =0; i<recoElectrons.size(); i++) {
        if(recoElectrons[i].pdgId()<0) Nem = Nem+1;
        if(recoElectrons[i].pdgId()>0) Nep = Nep+1;
    }

    if(Nmm>=1 && Nmp>=1) properLep_ID = true; //2mu + x
    if(Nem>=1 && Nep>=1) properLep_ID = true; //2e + x

    // proper charge flavor combination for Z + 1L
    if(!properLep_ID) return;

    if (verbose) cout<<"found three leptons"<<endl;

    // Consider all Z candidates
    double minZ1DeltaM=9999.9;
    for (int i=0; i<n_Zs; i++) {

        int i1 = Z_Z1L_lepindex1[i]; int i2 = Z_Z1L_lepindex2[i];
        int j1 = 3 - i1 - i2; // index of the third lepton (check if this works)

        TLorentzVector lep_i1, lep_i2, lep_j1;
        lep_i1.SetPtEtaPhiM(lepFSR_pt[i1],lepFSR_eta[i1],lepFSR_phi[i1],lepFSR_mass[i1]);
        lep_i2.SetPtEtaPhiM(lepFSR_pt[i2],lepFSR_eta[i2],lepFSR_phi[i2],lepFSR_mass[i2]);
        lep_j1.SetPtEtaPhiM(lepFSR_pt[j1],lepFSR_eta[j1],lepFSR_phi[j1],lepFSR_mass[j1]);

        TLorentzVector lep_i1_nofsr, lep_i2_nofsr, lep_j1_nofsr;
        lep_i1_nofsr.SetPtEtaPhiM(lep_pt[i1],lep_eta[i1],lep_phi[i1],lep_mass[i1]);
        lep_i2_nofsr.SetPtEtaPhiM(lep_pt[i2],lep_eta[i2],lep_phi[i2],lep_mass[i2]);
        lep_j1_nofsr.SetPtEtaPhiM(lep_pt[j1],lep_eta[j1],lep_phi[j1],lep_mass[j1]);

        TLorentzVector Zi;
        Zi = lep_i1+lep_i2;
        //Zi.SetPtEtaPhiM(Z_Z1L_pt[i],Z_Z1L_eta[i],Z_Z1L_phi[i],Z_Z1L_mass[i]);

        if (verbose) {cout<<"Z candidate Zi->M() "<<Zi.M()<<endl;}

        TLorentzVector Z1 = Zi;
        double Z1DeltaM = abs(Zi.M()-Zmass);
        int Z1_lepindex[2] = {0,0};
        if (lep_i1.Pt()>lep_i2.Pt()) { Z1_lepindex[0] = i1;  Z1_lepindex[1] = i2; }
        else { Z1_lepindex[0] = i2;  Z1_lepindex[1] = i1; }

        // Check Leading and Subleading pt Cut
        vector<double> allPt;
        allPt.push_back(lep_i1.Pt()); allPt.push_back(lep_i2.Pt());
        std::sort(allPt.begin(), allPt.end());
        if (verbose) cout<<" leading pt: "<<allPt[1]<<" cut: "<<leadingPtCut<<" subleadingPt: "<<allPt[0]<<" cut: "<<subleadingPtCut<<endl;
        if (allPt[1]<leadingPtCut || allPt[0]<subleadingPtCut ) continue;

        // Check dR(li,lj)>0.02 for any i,j
        vector<double> alldR;
        alldR.push_back(deltaR(lep_i1.Eta(),lep_i1.Phi(),lep_i2.Eta(),lep_i2.Phi()));
        alldR.push_back(deltaR(lep_i1.Eta(),lep_i1.Phi(),lep_j1.Eta(),lep_j1.Phi()));
        alldR.push_back(deltaR(lep_i2.Eta(),lep_i2.Phi(),lep_j1.Eta(),lep_j1.Phi()));
        if (verbose) cout<<" minDr: "<<*min_element(alldR.begin(),alldR.end())<<endl;
        if (*min_element(alldR.begin(),alldR.end())<0.02) continue;

        // Check M(l+,l-)>4.0 GeV for any OS pair
        // Do not include FSR photons
        vector<double> allM;
        TLorentzVector i1i2;
        i1i2 = (lep_i1_nofsr)+(lep_i2_nofsr); allM.push_back(i1i2.M());
        if (lep_id[i1]*lep_id[j1]<0) {
            TLorentzVector i1j1;
            i1j1 = (lep_i1_nofsr)+(lep_j1_nofsr); allM.push_back(i1j1.M());
        } else {
            TLorentzVector i2j1;
            i2j1 = (lep_i2_nofsr)+(lep_j1_nofsr); allM.push_back(i2j1.M());
        }
        if (verbose) cout<<" min m(l+l-): "<<*min_element(allM.begin(),allM.end())<<endl;
        if (*min_element(allM.begin(),allM.end())<4.0) {passedQCDcut=false; continue;}

        // Check isolation cut (without FSR ) for Z1 leptons
        if (lep_RelIsoNoFSR[Z1_lepindex[0]]>((abs(lep_id[Z1_lepindex[0]])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
        if (lep_RelIsoNoFSR[Z1_lepindex[1]]>((abs(lep_id[Z1_lepindex[1]])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
        // Check tight ID cut for Z1 leptons
        if (!(lep_tightId[Z1_lepindex[0]])) continue; // checking tight lepton ID
        if (!(lep_tightId[Z1_lepindex[1]])) continue; // checking tight lepton ID

        if ( (Z1.M() < mZ1Low) || (Z1.M() > mZ1High) ) continue;

        if (verbose) cout<<"good Z1L candidate, Z1DeltaM: "<<Z1DeltaM<<" minZ1DeltaM: "<<minZ1DeltaM<<endl;

        // Check if this candidate has the best Z1 and highest scalar sum of Z2 lepton pt

        if ( Z1DeltaM<=minZ1DeltaM ) {

            minZ1DeltaM = Z1DeltaM;

            TLorentzVector Z1L;
            Z1L = Z1+lep_j1;

            massZ1_Z1L = Z1.M();
            mass3l = Z1L.M();

            lep_Hindex[0] = Z1_lepindex[0];
            lep_Hindex[1] = Z1_lepindex[1];
            lep_Hindex[2] = j1;

            if (verbose) cout<<" new best Z1L candidate: massZ1: "<<massZ1<<" (mass3l: "<<mass3l<<")"<<endl;
            foundZ1LCandidate=true;

        }
    }   
}


void UFHZZ4LAna::bookPassedEventTree(TString treeName, TTree *tree)
{     

    using namespace edm;
    using namespace pat;
    using namespace std;


    // -------------------------                                                                                                                                                                        
    // RECO level information                                                                                                                                                                           
    // -------------------------                                                                                                                                                                        
    // Event variables
    tree->Branch("Run",&Run,"Run/l");
    tree->Branch("Event",&Event,"Event/l");
    tree->Branch("LumiSect",&LumiSect,"LumiSect/l");
    tree->Branch("nVtx",&nVtx,"nVtx/I");
    tree->Branch("nInt",&nInt,"nInt/I");
    tree->Branch("finalState",&finalState,"finalState/I");
    tree->Branch("triggersPassed",&triggersPassed);
    tree->Branch("passedTrig",&passedTrig,"passedTrig/O");
    tree->Branch("passedFullSelection",&passedFullSelection,"passedFullSelection/O");
    tree->Branch("passedZ4lSelection",&passedZ4lSelection,"passedZ4lSelection/O");
    tree->Branch("passedQCDcut",&passedQCDcut,"passedQCDcut/O");
    tree->Branch("passedZ1LSelection",&passedZ1LSelection,"passedZ1LSelection/O");
    tree->Branch("passedZ4lZ1LSelection",&passedZ4lZ1LSelection,"passedZ4lZ1LSelection/O");
    tree->Branch("passedZ4lZXCRSelection",&passedZ4lZXCRSelection,"passedZ4lZXCRSelection/O");
    tree->Branch("passedZXCRSelection",&passedZXCRSelection,"passedZXCRSelection/O");
    tree->Branch("nZXCRFailedLeptons",&nZXCRFailedLeptons,"nZXCRFailedLeptons/I");
    tree->Branch("genWeight",&genWeight,"genWeight/F");
    tree->Branch("k_ggZZ",&k_ggZZ,"k_ggZZ/F");
    tree->Branch("k_qqZZ_qcd_dPhi",&k_qqZZ_qcd_dPhi,"k_qqZZ_qcd_dPhi/F");
    tree->Branch("k_qqZZ_qcd_M",&k_qqZZ_qcd_M,"k_qqZZ_qcd_M/F");
    tree->Branch("k_qqZZ_qcd_Pt",&k_qqZZ_qcd_Pt,"k_qqZZ_qcd_Pt/F");
    tree->Branch("k_qqZZ_ewk",&k_qqZZ_ewk,"k_qqZZ_ewk/F");
    tree->Branch("qcdWeights",&qcdWeights);
    tree->Branch("nnloWeights",&nnloWeights);
    tree->Branch("pdfWeights",&pdfWeights);
    tree->Branch("pdfRMSup",&pdfRMSup,"pdfRMSup/F");
    tree->Branch("pdfRMSdown",&pdfRMSdown,"pdfRMSdown/F");
    tree->Branch("pdfENVup",&pdfENVup,"pdfENVup/F");
    tree->Branch("pdfENVdown",&pdfENVdown,"pdfENVdown/F");
    tree->Branch("pileupWeight",&pileupWeight,"pileupWeight/F");
    tree->Branch("pileupWeightUp",&pileupWeightUp,"pileupWeightUp/F");
    tree->Branch("pileupWeightDn",&pileupWeightDn,"pileupWeightDn/F");
    tree->Branch("dataMCWeight",&dataMCWeight,"dataMCWeight/F");
    tree->Branch("eventWeight",&eventWeight,"eventWeight/F");
    tree->Branch("prefiringWeight",&prefiringWeight,"prefiringWeight/F");
    tree->Branch("crossSection",&crossSection,"crossSection/F");

    // Lepton variables
    tree->Branch("lep_id",&lep_id);
    tree->Branch("lep_pt",&lep_pt_float);
    tree->Branch("lep_pterr",&lep_pterr_float);
    tree->Branch("lep_pterrold",&lep_pterrold_float);
    tree->Branch("lep_p",&lep_p_float);
    tree->Branch("lep_ecalEnergy",&lep_ecalEnergy_float);
    tree->Branch("lep_isEB",&lep_isEB);
    tree->Branch("lep_isEE",&lep_isEE);
    tree->Branch("lep_eta",&lep_eta_float);
    tree->Branch("lep_phi",&lep_phi_float);
    tree->Branch("lep_mass",&lep_mass_float);
    tree->Branch("lepFSR_pt",&lepFSR_pt_float);
    tree->Branch("lepFSR_eta",&lepFSR_eta_float);
    tree->Branch("lepFSR_phi",&lepFSR_phi_float);
    tree->Branch("lepFSR_mass",&lepFSR_mass_float);
    tree->Branch("lep_Hindex",&lep_Hindex,"lep_Hindex[4]/I");
    tree->Branch("lep_genindex",&lep_genindex);
    tree->Branch("lep_matchedR03_PdgId",&lep_matchedR03_PdgId);
    tree->Branch("lep_matchedR03_MomId",&lep_matchedR03_MomId);
    tree->Branch("lep_matchedR03_MomMomId",&lep_matchedR03_MomMomId);
    tree->Branch("lep_missingHits",&lep_missingHits);
    tree->Branch("lep_mva",&lep_mva);
    tree->Branch("lep_ecalDriven",&lep_ecalDriven);
    tree->Branch("lep_tightId",&lep_tightId);
    tree->Branch("lep_tightIdSUS",&lep_tightIdSUS);
    tree->Branch("lep_tightIdHiPt",&lep_tightIdHiPt);
    //tree->Branch("lep_tightId_old",&lep_tightId_old);
    tree->Branch("lep_Sip",&lep_Sip);
    tree->Branch("lep_IP",&lep_IP);
    tree->Branch("lep_isoNH",&lep_isoNH);
    tree->Branch("lep_isoCH",&lep_isoCH);
    tree->Branch("lep_isoPhot",&lep_isoPhot);
    tree->Branch("lep_isoPU",&lep_isoPU);
    tree->Branch("lep_isoPUcorr",&lep_isoPUcorr);
    tree->Branch("lep_RelIso",&lep_RelIso);
    tree->Branch("lep_RelIsoNoFSR",&lep_RelIsoNoFSR);
    tree->Branch("lep_MiniIso",&lep_MiniIso);
    tree->Branch("lep_ptRatio",&lep_ptRatio);
    tree->Branch("lep_ptRel",&lep_ptRel);
    tree->Branch("lep_filtersMatched",&lep_filtersMatched);
    tree->Branch("lep_dataMC",&lep_dataMC);
    tree->Branch("lep_dataMCErr",&lep_dataMCErr);
    tree->Branch("nisoleptons",&nisoleptons,"nisoleptons/I");
    tree->Branch("muRho",&muRho,"muRho/F");
    tree->Branch("elRho",&elRho,"elRho/F");
    tree->Branch("pTL1",&pTL1,"pTL1/F");
    tree->Branch("pTL2",&pTL2,"pTL2/F");
    tree->Branch("pTL3",&pTL3,"pTL3/F");
    tree->Branch("pTL4",&pTL4,"pTL4/F");
    tree->Branch("idL1",&idL1,"idL1/I");
    tree->Branch("idL2",&idL2,"idL2/I");
    tree->Branch("idL3",&idL3,"idL3/I");
    tree->Branch("idL4",&idL4,"idL4/I");
    tree->Branch("etaL1",&etaL1,"etaL1/F");
    tree->Branch("etaL2",&etaL2,"etaL2/F");
    tree->Branch("etaL3",&etaL3,"etaL3/F");
    tree->Branch("etaL4",&etaL4,"etaL4/F");
    tree->Branch("mL1",&mL1,"mL1/F");
    tree->Branch("mL2",&mL2,"mL2/F");
    tree->Branch("mL3",&mL3,"mL3/F");
    tree->Branch("mL4",&mL4,"mL4/F");
    tree->Branch("pTErrL1",&pTErrL1,"pTErrL1/F");
    tree->Branch("pTErrL2",&pTErrL2,"pTErrL2/F");
    tree->Branch("pTErrL3",&pTErrL3,"pTErrL3/F");
    tree->Branch("pTErrL4",&pTErrL4,"pTErrL4/F");
    tree->Branch("phiL1",&phiL1,"phiL1/F");
    tree->Branch("phiL2",&phiL2,"phiL2/F");
    tree->Branch("phiL3",&phiL3,"phiL3/F");
    tree->Branch("phiL4",&phiL4,"phiL4/F");
    tree->Branch("pTL1FSR",&pTL1FSR,"pTL1FSR/F");
    tree->Branch("pTL2FSR",&pTL2FSR,"pTL2FSR/F");
    tree->Branch("pTL3FSR",&pTL3FSR,"pTL3FSR/F");
    tree->Branch("pTL4FSR",&pTL4FSR,"pTL4FSR/F");
    tree->Branch("tau_id",&tau_id);
    tree->Branch("tau_pt",&tau_pt_float);
    tree->Branch("tau_eta",&tau_eta_float);
    tree->Branch("tau_phi",&tau_phi_float);
    tree->Branch("tau_mass",&tau_mass_float);
    tree->Branch("pho_pt",&pho_pt_float);
    tree->Branch("pho_eta",&pho_eta_float);
    tree->Branch("pho_phi",&pho_phi_float);
    tree->Branch("photonCutBasedIDLoose",&photonCutBasedIDLoose_float);

    //Higgs Candidate Variables
    tree->Branch("H_pt",&H_pt_float);
    tree->Branch("H_eta",&H_eta_float);
    tree->Branch("H_phi",&H_phi_float);
    tree->Branch("H_mass",&H_mass_float);
    tree->Branch("H_noFSR_pt",&H_noFSR_pt_float);
    tree->Branch("H_noFSR_eta",&H_noFSR_eta_float);
    tree->Branch("H_noFSR_phi",&H_noFSR_phi_float);
    tree->Branch("H_noFSR_mass",&H_noFSR_mass_float);
    tree->Branch("mass4l",&mass4l,"mass4l/F");
    tree->Branch("mass4l_noFSR",&mass4l_noFSR,"mass4l_noFSR/F");
    tree->Branch("mass4lErr",&mass4lErr,"mass4lErr/F");
    tree->Branch("mass4lREFIT",&mass4lREFIT,"mass4lREFIT/F");
    tree->Branch("mass4lErrREFIT",&mass4lErrREFIT,"mass4lErrREFIT/F");
    tree->Branch("massZ1REFIT",&massZ1REFIT,"massZ1REFIT/F");
    tree->Branch("massZ2REFIT",&massZ2REFIT,"massZ2REFIT/F");
    tree->Branch("mass4mu",&mass4mu,"mass4mu/F");
    tree->Branch("mass4e",&mass4e,"mass4e/F");
    tree->Branch("mass2e2mu",&mass2e2mu,"mass2e2mu/F");
    tree->Branch("pT4l",&pT4l,"pT4l/F");
    tree->Branch("eta4l",&eta4l,"eta4l/F");
    tree->Branch("phi4l",&phi4l,"phi4l/F");
    tree->Branch("rapidity4l",&rapidity4l,"rapidity4l/F");
    tree->Branch("cosTheta1",&cosTheta1,"cosTheta1/F");
    tree->Branch("cosTheta2",&cosTheta2,"cosTheta2/F");
    tree->Branch("cosThetaStar",&cosThetaStar,"cosThetaStar/F");
    tree->Branch("Phi",&Phi,"Phi/F");
    tree->Branch("Phi1",&Phi1,"Phi1/F");
    tree->Branch("mass3l",&mass3l,"mass3l/F");

    // Z candidate variables
    tree->Branch("Z_pt",&Z_pt_float);
    tree->Branch("Z_eta",&Z_eta_float);
    tree->Branch("Z_phi",&Z_phi_float);
    tree->Branch("Z_mass",&Z_mass_float);
    tree->Branch("Z_noFSR_pt",&Z_noFSR_pt_float);
    tree->Branch("Z_noFSR_eta",&Z_noFSR_eta_float);
    tree->Branch("Z_noFSR_phi",&Z_noFSR_phi_float);
    tree->Branch("Z_noFSR_mass",&Z_noFSR_mass_float);
    tree->Branch("Z_Hindex",&Z_Hindex,"Z_Hindex[2]/I");
    tree->Branch("massZ1",&massZ1,"massZ1/F");
    tree->Branch("massZ1_Z1L",&massZ1_Z1L,"massZ1_Z1L/F");
    tree->Branch("massZ2",&massZ2,"massZ2/F");  
    tree->Branch("pTZ1",&pTZ1,"pTZ1/F");
    tree->Branch("pTZ2",&pTZ2,"pTZ2/F");

    // MET
    tree->Branch("met",&met,"met/F");
    tree->Branch("met_phi",&met_phi,"met_phi/F");
    tree->Branch("met_jesup",&met_jesup,"met_jesup/F");
    tree->Branch("met_phi_jesup",&met_phi_jesup,"met_phi_jesup/F");
    tree->Branch("met_jesdn",&met_jesdn,"met_jesdn/F");
    tree->Branch("met_phi_jesdn",&met_phi_jesdn,"met_phi_jesdn/F");
    tree->Branch("met_uncenup",&met_uncenup,"met_uncenup/F");
    tree->Branch("met_phi_uncenup",&met_phi_uncenup,"met_phi_uncenup/F");
    tree->Branch("met_uncendn",&met_uncendn,"met_uncendn/F");
    tree->Branch("met_phi_uncendn",&met_phi_uncendn,"met_phi_uncendn/F");

    // Jets
    tree->Branch("jet_iscleanH4l",&jet_iscleanH4l);
    tree->Branch("jet1index",&jet1index,"jet1index/I");
    tree->Branch("jet2index",&jet2index,"jet2index/I");
    tree->Branch("jet_pt",&jet_pt_float);
    tree->Branch("jet_pt_raw",&jet_pt_raw_float);
    tree->Branch("jet_relpterr",&jet_relpterr);    
    tree->Branch("jet_eta",&jet_eta_float);
    tree->Branch("jet_phi",&jet_phi_float);
    tree->Branch("jet_phierr",&jet_phierr);
    tree->Branch("jet_bTagEffi",&jet_bTagEffi);
    tree->Branch("jet_cTagEffi",&jet_cTagEffi);
    tree->Branch("jet_udsgTagEffi",&jet_udsgTagEffi);
    tree->Branch("jet_mass",&jet_mass_float);    
    tree->Branch("jet_jesup_iscleanH4l",&jet_jesup_iscleanH4l);
    tree->Branch("jet_jesup_pt",&jet_jesup_pt_float);
    tree->Branch("jet_jesup_eta",&jet_jesup_eta_float);
    tree->Branch("jet_jesup_phi",&jet_jesup_phi_float);
    tree->Branch("jet_jesup_mass",&jet_jesup_mass_float);
    tree->Branch("jet_jesdn_iscleanH4l",&jet_jesdn_iscleanH4l);
    tree->Branch("jet_jesdn_pt",&jet_jesdn_pt_float);
    tree->Branch("jet_jesdn_eta",&jet_jesdn_eta_float);
    tree->Branch("jet_jesdn_phi",&jet_jesdn_phi_float);
    tree->Branch("jet_jesdn_mass",&jet_jesdn_mass_float);
    tree->Branch("jet_jerup_iscleanH4l",&jet_jerup_iscleanH4l);
    tree->Branch("jet_jerup_pt",&jet_jerup_pt_float);
    tree->Branch("jet_jerup_eta",&jet_jerup_eta_float);
    tree->Branch("jet_jerup_phi",&jet_jerup_phi_float);
    tree->Branch("jet_jerup_mass",&jet_jerup_mass_float);
    tree->Branch("jet_jerdn_iscleanH4l",&jet_jerdn_iscleanH4l);
    tree->Branch("jet_jerdn_pt",&jet_jerdn_pt_float);
    tree->Branch("jet_jerdn_eta",&jet_jerdn_eta_float);
    tree->Branch("jet_jerdn_phi",&jet_jerdn_phi_float);
    tree->Branch("jet_jerdn_mass",&jet_jerdn_mass_float);
    tree->Branch("jet_pumva",&jet_pumva);
    tree->Branch("jet_csvv2",&jet_csvv2);
    tree->Branch("jet_csvv2_",&jet_csvv2_);
    tree->Branch("jet_isbtag",&jet_isbtag);
    tree->Branch("jet_hadronFlavour",&jet_hadronFlavour);
    tree->Branch("jet_partonFlavour",&jet_partonFlavour);    
    tree->Branch("jet_QGTagger",&jet_QGTagger);
    tree->Branch("jet_QGTagger_jesup",&jet_QGTagger_jesup);
    tree->Branch("jet_QGTagger_jesdn",&jet_QGTagger_jesdn);
    tree->Branch("jet_axis2",&jet_axis2);
    tree->Branch("jet_ptD",&jet_ptD);
    tree->Branch("jet_mult",&jet_mult);
    tree->Branch("njets_pt30_eta4p7",&njets_pt30_eta4p7,"njets_pt30_eta4p7/I");
    tree->Branch("njets_pt30_eta4p7_jesup",&njets_pt30_eta4p7_jesup,"njets_pt30_eta4p7_jesup/I");
    tree->Branch("njets_pt30_eta4p7_jesdn",&njets_pt30_eta4p7_jesdn,"njets_pt30_eta4p7_jesdn/I");
    tree->Branch("njets_pt30_eta4p7_jerup",&njets_pt30_eta4p7_jerup,"njets_pt30_eta4p7_jerup/I");
    tree->Branch("njets_pt30_eta4p7_jerdn",&njets_pt30_eta4p7_jerdn,"njets_pt30_eta4p7_jerdn/I");
    tree->Branch("pt_leadingjet_pt30_eta4p7",&pt_leadingjet_pt30_eta4p7,"pt_leadingjet_pt30_eta4p7/F");
    tree->Branch("pt_leadingjet_pt30_eta4p7_jesup",&pt_leadingjet_pt30_eta4p7_jesup,"pt_leadingjet_pt30_eta4p7_jesup/F");
    tree->Branch("pt_leadingjet_pt30_eta4p7_jesdn",&pt_leadingjet_pt30_eta4p7_jesdn,"pt_leadingjet_pt30_eta4p7_jesdn/F");
    tree->Branch("pt_leadingjet_pt30_eta4p7_jerup",&pt_leadingjet_pt30_eta4p7_jerup,"pt_leadingjet_pt30_eta4p7_jerup/F");
    tree->Branch("pt_leadingjet_pt30_eta4p7_jerdn",&pt_leadingjet_pt30_eta4p7_jerdn,"pt_leadingjet_pt30_eta4p7_jerdn/F");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7",&absrapidity_leadingjet_pt30_eta4p7,"absrapidity_leadingjet_pt30_eta4p7/F");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jesup",&absrapidity_leadingjet_pt30_eta4p7_jesup,"absrapidity_leadingjet_pt30_eta4p7_jesup/F");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jesdn",&absrapidity_leadingjet_pt30_eta4p7_jesdn,"absrapidity_leadingjet_pt30_eta4p7_jesdn/F");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jerup",&absrapidity_leadingjet_pt30_eta4p7_jerup,"absrapidity_leadingjet_pt30_eta4p7_jerup/F");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jerdn",&absrapidity_leadingjet_pt30_eta4p7_jerdn,"absrapidity_leadingjet_pt30_eta4p7_jerdn/F");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7",&absdeltarapidity_hleadingjet_pt30_eta4p7,"absdeltarapidity_hleadingjet_pt30_eta4p7/F");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jesup",&absdeltarapidity_hleadingjet_pt30_eta4p7_jesup,"absdeltarapidity_hleadingjet_pt30_eta4p7_jesup/F");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn",&absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn,"absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn/F");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jerup",&absdeltarapidity_hleadingjet_pt30_eta4p7_jerup,"absdeltarapidity_hleadingjet_pt30_eta4p7_jerup/F");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn",&absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn,"absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn/F");
    tree->Branch("nbjets_pt30_eta4p7",&nbjets_pt30_eta4p7,"nbjets_pt30_eta4p7/I");
    tree->Branch("nvjets_pt40_eta2p4",&nvjets_pt40_eta2p4,"nvjets_pt40_eta2p4/I");
    tree->Branch("DijetMass",&DijetMass,"DijetMass/F");
    tree->Branch("DijetDEta",&DijetDEta,"DijetDEta/F");
    tree->Branch("DijetFisher",&DijetFisher,"DijetFisher/F");
    tree->Branch("njets_pt30_eta2p5",&njets_pt30_eta2p5,"njets_pt30_eta2p5/I");
    tree->Branch("njets_pt30_eta2p5_jesup",&njets_pt30_eta2p5_jesup,"njets_pt30_eta2p5_jesup/I");
    tree->Branch("njets_pt30_eta2p5_jesdn",&njets_pt30_eta2p5_jesdn,"njets_pt30_eta2p5_jesdn/I");
    tree->Branch("njets_pt30_eta2p5_jerup",&njets_pt30_eta2p5_jerup,"njets_pt30_eta2p5_jerup/I");
    tree->Branch("njets_pt30_eta2p5_jerdn",&njets_pt30_eta2p5_jerdn,"njets_pt30_eta2p5_jerdn/I");
    tree->Branch("pt_leadingjet_pt30_eta2p5",&pt_leadingjet_pt30_eta2p5,"pt_leadingjet_pt30_eta2p5/F");
    tree->Branch("pt_leadingjet_pt30_eta2p5_jesup",&pt_leadingjet_pt30_eta2p5_jesup,"pt_leadingjet_pt30_eta2p5_jesup/F");
    tree->Branch("pt_leadingjet_pt30_eta2p5_jesdn",&pt_leadingjet_pt30_eta2p5_jesdn,"pt_leadingjet_pt30_eta2p5_jesdn/F");
    tree->Branch("pt_leadingjet_pt30_eta2p5_jerup",&pt_leadingjet_pt30_eta2p5_jerup,"pt_leadingjet_pt30_eta2p5_jerup/F");
    tree->Branch("pt_leadingjet_pt30_eta2p5_jerdn",&pt_leadingjet_pt30_eta2p5_jerdn,"pt_leadingjet_pt30_eta2p5_jerdn/F");


    // merged jets
    tree->Branch("mergedjet_iscleanH4l",&mergedjet_iscleanH4l);
    tree->Branch("mergedjet_pt",&mergedjet_pt);
    tree->Branch("mergedjet_eta",&mergedjet_eta);
    tree->Branch("mergedjet_phi",&mergedjet_phi);
    tree->Branch("mergedjet_mass",&mergedjet_mass);    
    tree->Branch("mergedjet_tau1",&mergedjet_tau1);
    tree->Branch("mergedjet_tau2",&mergedjet_tau2);
    tree->Branch("mergedjet_btag",&mergedjet_btag);
    
    tree->Branch("mergedjet_L1",&mergedjet_L1);
    tree->Branch("mergedjet_softdropmass",&mergedjet_softdropmass);
    tree->Branch("mergedjet_prunedmass",&mergedjet_prunedmass);

    tree->Branch("mergedjet_nsubjet",&mergedjet_nsubjet);
    tree->Branch("mergedjet_subjet_pt",&mergedjet_subjet_pt);
    tree->Branch("mergedjet_subjet_eta",&mergedjet_subjet_eta);
    tree->Branch("mergedjet_subjet_phi",&mergedjet_subjet_phi);
    tree->Branch("mergedjet_subjet_mass",&mergedjet_subjet_mass);
    tree->Branch("mergedjet_subjet_btag",&mergedjet_subjet_btag);
    tree->Branch("mergedjet_subjet_partonFlavour",&mergedjet_subjet_partonFlavour);
    tree->Branch("mergedjet_subjet_hadronFlavour",&mergedjet_subjet_hadronFlavour);

    // FSR Photons
    tree->Branch("nFSRPhotons",&nFSRPhotons,"nFSRPhotons/I");
    tree->Branch("allfsrPhotons_dR",&allfsrPhotons_dR);
    tree->Branch("allfsrPhotons_iso",&allfsrPhotons_iso);
    tree->Branch("allfsrPhotons_pt",&allfsrPhotons_pt);
    tree->Branch("fsrPhotons_lepindex",&fsrPhotons_lepindex);
    tree->Branch("fsrPhotons_pt",&fsrPhotons_pt_float);
    tree->Branch("fsrPhotons_pterr",&fsrPhotons_pterr_float);
    tree->Branch("fsrPhotons_eta",&fsrPhotons_eta_float);
    tree->Branch("fsrPhotons_phi",&fsrPhotons_phi_float);
    tree->Branch("fsrPhotons_dR",&fsrPhotons_dR);
    tree->Branch("fsrPhotons_iso",&fsrPhotons_iso);

    // Z4l? FIXME
    tree->Branch("theta12",&theta12,"theta12/F"); 
    tree->Branch("theta13",&theta13,"theta13/F"); 
    tree->Branch("theta14",&theta14,"theta14/F");
    tree->Branch("minM3l",&minM3l,"minM3l/F"); 
    tree->Branch("Z4lmaxP",&Z4lmaxP,"Z4lmaxP/F"); 
    tree->Branch("minDeltR",&minDeltR,"minDeltR/F"); 
    tree->Branch("m3l_soft",&m3l_soft,"m3l_soft/F");
    tree->Branch("minMass2Lep",&minMass2Lep,"minMass2Lep/F"); 
    tree->Branch("maxMass2Lep",&maxMass2Lep,"maxMass2Lep/F");
    tree->Branch("thetaPhoton",&thetaPhoton,"thetaPhoton/F"); 
    tree->Branch("thetaPhotonZ",&thetaPhotonZ,"thetaPhotonZ/F");

    // Event Category
    tree->Branch("EventCat",&EventCat,"EventCat/I");

    // -------------------------                                                                                                                                                                        
    // GEN level information                                                                                                                                                                            
    // -------------------------                                                                                                                                                                        
    //Event variables
    tree->Branch("GENfinalState",&GENfinalState,"GENfinalState/I");
    tree->Branch("passedFiducialSelection",&passedFiducialSelection,"passedFiducialSelection/O");

    // lepton variables
    tree->Branch("GENlep_pt",&GENlep_pt_float);
    tree->Branch("GENlep_eta",&GENlep_eta_float);
    tree->Branch("GENlep_phi",&GENlep_phi_float);
    tree->Branch("GENlep_mass",&GENlep_mass_float);
    tree->Branch("GENlep_id",&GENlep_id);
    tree->Branch("GENlep_status",&GENlep_status);
    tree->Branch("GENlep_MomId",&GENlep_MomId);
    tree->Branch("GENlep_MomMomId",&GENlep_MomMomId);
    tree->Branch("GENlep_Hindex",&GENlep_Hindex,"GENlep_Hindex[4]/I");
    tree->Branch("GENlep_isoCH",&GENlep_isoCH);
    tree->Branch("GENlep_isoNH",&GENlep_isoNH);
    tree->Branch("GENlep_isoPhot",&GENlep_isoPhot);
    tree->Branch("GENlep_RelIso",&GENlep_RelIso);

    // Higgs candidate variables (calculated using selected gen leptons)
    tree->Branch("GENH_pt",&GENH_pt_float);
    tree->Branch("GENH_eta",&GENH_eta_float);
    tree->Branch("GENH_phi",&GENH_phi_float);
    tree->Branch("GENH_mass",&GENH_mass_float);
    tree->Branch("GENmass4l",&GENmass4l,"GENmass4l/F");
    tree->Branch("GENmass4mu",&GENmass4mu,"GENmass4mu/F");
    tree->Branch("GENmass4e",&GENmass4e,"GENmass4e/F");
    tree->Branch("GENmass2e2mu",&GENmass2e2mu,"GENmass2e2mu/F");
    tree->Branch("GENpT4l",&GENpT4l,"GENpT4l/F");
    tree->Branch("GENeta4l",&GENeta4l,"GENeta4l/F");
    tree->Branch("GENrapidity4l",&GENrapidity4l,"GENrapidity4l/F");
    tree->Branch("GENcosTheta1",&GENcosTheta1,"GENcosTheta1/F");
    tree->Branch("GENcosTheta2",&GENcosTheta2,"GENcosTheta2/F");
    tree->Branch("GENcosThetaStar",&GENcosThetaStar,"GENcosThetaStar/F");
    tree->Branch("GENPhi",&GENPhi,"GENPhi/F");
    tree->Branch("GENPhi1",&GENPhi1,"GENPhi1/F");
    tree->Branch("GENMH",&GENMH,"GENMH/F");

    // Z candidate variables
    tree->Branch("GENZ_pt",&GENZ_pt_float);
    tree->Branch("GENZ_eta",&GENZ_eta_float);
    tree->Branch("GENZ_phi",&GENZ_phi_float);
    tree->Branch("GENZ_mass",&GENZ_mass_float);
    tree->Branch("GENZ_DaughtersId",&GENZ_DaughtersId); 
    tree->Branch("GENZ_MomId",&GENZ_MomId);
    tree->Branch("GENmassZ1",&GENmassZ1,"GENmassZ1/F");
    tree->Branch("GENmassZ2",&GENmassZ2,"GENmassZ2/F");  
    tree->Branch("GENpTZ1",&GENpTZ1,"GENpTZ1/F");
    tree->Branch("GENpTZ2",&GENpTZ2,"GENpTZ2/F");
    tree->Branch("GENdPhiZZ",&GENdPhiZZ,"GENdPhiZZ/F");
    tree->Branch("GENmassZZ",&GENmassZZ,"GENmassZZ/F");
    tree->Branch("GENpTZZ",&GENpTZZ,"GENpTZZ/F");

    // Higgs variables directly from GEN particle
    tree->Branch("GENHmass",&GENHmass,"GENHmass/F");

    // STXS 
    tree->Branch("stage0cat",&stage0cat,"stage0cat/I");
    tree->Branch("stage1cat",&stage1cat,"stage1cat/I");
    tree->Branch("stage1p1cat",&stage1p1cat,"stage1p1cat/I");
    tree->Branch("stage1p2cat",&stage1p2cat,"stage1p2cat/I");
    // Fiducial Rivet
    tree->Branch("passedFiducialRivet",&passedFiducialRivet,"passedFiducialRivet/I");
    tree->Branch("GENpT4lRivet",&GENpT4lRivet,"GENpT4lRivet/F");
    tree->Branch("GENnjets_pt30_eta4p7Rivet",&GENnjets_pt30_eta4p7Rivet,"GENnjets_pt30_eta4p7Rivet/I");
    tree->Branch("GENpt_leadingjet_pt30_eta4p7Rivet",&GENpt_leadingjet_pt30_eta4p7Rivet,"GENpt_leadingjet_pt30_eta4p7Rivet/F");

    // Jets
    tree->Branch("GENjet_pt",&GENjet_pt_float);
    tree->Branch("GENjet_eta",&GENjet_eta_float);
    tree->Branch("GENjet_phi",&GENjet_phi_float);
    tree->Branch("GENjet_mass",&GENjet_mass_float);
    tree->Branch("GENnjets_pt30_eta4p7",&GENnjets_pt30_eta4p7,"GENnjets_pt30_eta4p7/I");
    tree->Branch("GENpt_leadingjet_pt30_eta4p7",&GENpt_leadingjet_pt30_eta4p7,"GENpt_leadingjet_pt30_eta4p7/F");
    tree->Branch("GENabsrapidity_leadingjet_pt30_eta4p7",&GENabsrapidity_leadingjet_pt30_eta4p7,"GENabsrapidity_leadingjet_pt30_eta4p7/F");
    tree->Branch("GENabsdeltarapidity_hleadingjet_pt30_eta4p7",&GENabsdeltarapidity_hleadingjet_pt30_eta4p7,"GENabsdeltarapidity_hleadingjet_pt30_eta4p7/F");
    tree->Branch("GENnjets_pt30_eta2p5",&GENnjets_pt30_eta2p5,"GENnjets_pt30_eta2p5/I");
    tree->Branch("GENpt_leadingjet_pt30_eta2p5",&GENpt_leadingjet_pt30_eta2p5,"GENpt_leadingjet_pt30_eta2p5/F");
    tree->Branch("lheNj",&lheNj,"lheNj/I");
    tree->Branch("lheNb",&lheNb,"lheNb/I");
    tree->Branch("nGenStatus2bHad",&nGenStatus2bHad,"nGenStatus2bHad/I");

    //ME
    tree->Branch("me_0plus_JHU", &me_0plus_JHU, "me_0plus_JHU/F");
    tree->Branch("me_qqZZ_MCFM", &me_qqZZ_MCFM, "me_qqZZ_MCFM/F");
    tree->Branch("p0plus_m4l", &p0plus_m4l, "p0plus_m4l/F");
    tree->Branch("bkg_m4l", &bkg_m4l, "bkg_m4l/F");
    tree->Branch("D_bkg_kin", &D_bkg_kin, "D_bkg_kin/F");
    tree->Branch("D_bkg", &D_bkg, "D_bkg/F");
    tree->Branch("Dgg10_VAMCFM", &Dgg10_VAMCFM, "Dgg10_VAMCFM/F");
    tree->Branch("D_g4", &D_g4, "D_g4/F");
    tree->Branch("D_g1g4", &D_g4, "D_g1g4/F");
    tree->Branch("D_VBF",&D_VBF,"D_VBF/F");
    tree->Branch("D_VBF1j",&D_VBF1j,"D_VBF1j/F");
    tree->Branch("D_HadWH",&D_HadWH,"D_HadWH/F");
    tree->Branch("D_HadZH",&D_HadZH,"D_HadZH/F");
    tree->Branch("D_bkg_VBFdec",&D_bkg_VBFdec,"D_bkg_VBFdec/F");
    tree->Branch("D_bkg_VHdec",&D_bkg_VHdec,"D_bkg_VHdec/F");
    tree->Branch("D_VBF_QG",&D_VBF_QG,"D_VBF_QG/F");
    tree->Branch("D_VBF1j_QG",&D_VBF1j_QG,"D_VBF1j_QG/F");
    tree->Branch("D_HadWH_QG",&D_HadWH_QG,"D_HadWH_QG/F");
    tree->Branch("D_HadZH_QG",&D_HadZH_QG,"D_HadZH_QG/F");
	
	//ele scale unc
	tree->Branch("eleUncorr_pt",&eleUncorr_pt_float);
	tree->Branch("eleUncorr_mass",&eleUncorr_mass_float);
	tree->Branch("muonUncorr_pt",&muonUncorr_pt);
	tree->Branch("muon_scale_sf",&muon_scale_sf);
	//vertex constriant
	tree->Branch("pv_x",&pv_x,"pv_x/D");
	tree->Branch("pv_y",&pv_y,"pv_y/D");
	tree->Branch("pv_z",&pv_z,"pv_z/D");
	tree->Branch("bs_x",&bs_x,"bs_x/D");
	tree->Branch("bs_y",&bs_y,"bs_y/D");
	tree->Branch("bs_z",&bs_z,"bs_z/D");
	tree->Branch("passLepindex",&passLepindex);
	tree->Branch("track_vtx_dxy",&track_vtx_dxy_float);
	tree->Branch("track_bs_dxy",&track_bs_dxy_float);
	tree->Branch("track_vtx_xcorr",&track_vtx_xcorr,"track_vtx_xcorr/D");
	tree->Branch("track_vtx_ycorr",&track_vtx_ycorr,"track_vtx_ycorr/D");
	tree->Branch("track_vtx_zcorr",&track_vtx_zcorr,"track_vtx_zcorr/D");
	tree->Branch("track_vtx_dxycorr",&track_vtx_dxycorr_float);
	tree->Branch("track_vtx_pt",&track_vtx_pt_float);
	tree->Branch("track_vtx_eta",&track_vtx_eta_float);
	tree->Branch("track_vtx_phi",&track_vtx_phi_float);
	tree->Branch("track_vtx_pterr",&track_vtx_pterr_float);
	tree->Branch("track_bs_xcorr",&track_bs_xcorr,"track_bs_xcorr/D");
	tree->Branch("track_bs_ycorr",&track_bs_ycorr,"track_bs_ycorr/D");
	tree->Branch("track_bs_zcorr",&track_bs_zcorr,"track_bs_zcorr/D");
	tree->Branch("track_bs_dxycorr",&track_bs_dxycorr_float);
	tree->Branch("track_bs_pt",&track_bs_pt_float);
	tree->Branch("track_bs_eta",&track_bs_eta_float);
	tree->Branch("track_bs_phi",&track_bs_phi_float);
	tree->Branch("track_bs_pterr",&track_bs_pterr_float);
//	tree->Branch("track_vtx_ptcorr",&track_vtx_ptcorr_float);
//	tree->Branch("track_vtx_pterrcorr",&track_vtx_pterrcorr_float);
//	tree->Branch("track_vtx_FSRpt",&track_vtx_FSRpt_float);
//	tree->Branch("track_vtx_FSReta",&track_vtx_FSReta_float);
//	tree->Branch("track_vtx_FSRphi",&track_vtx_FSRphi_float);
//	tree->Branch("track_vtx_FSRmass",&track_vtx_FSRmass_float);
//	tree->Branch("track_bs_ptcorr",&track_bs_ptcorr_float);
//	tree->Branch("track_bs_pterrcorr",&track_bs_pterrcorr_float);
//	tree->Branch("track_bs_FSRpt",&track_bs_FSRpt_float);
//	tree->Branch("track_bs_FSReta",&track_bs_FSReta_float);
//	tree->Branch("track_bs_FSRphi",&track_bs_FSRphi_float);
//	tree->Branch("track_bs_FSRmass",&track_bs_FSRmass_float);

	
	tree->Branch("track_Hvtx_xcorr",&track_Hvtx_xcorr,"track_Hvtx_xcorr/D");
	tree->Branch("track_Hvtx_ycorr",&track_Hvtx_ycorr,"track_Hvtx_ycorr/D");
	tree->Branch("track_Hvtx_zcorr",&track_Hvtx_zcorr,"track_Hvtx_zcorr/D");
	tree->Branch("track_Hbs_xcorr",&track_Hbs_xcorr,"track_Hbs_xcorr/D");
	tree->Branch("track_Hbs_ycorr",&track_Hbs_ycorr,"track_Hbs_ycorr/D");
	tree->Branch("track_Hvtx_pt",&track_Hvtx_pt_float);
	tree->Branch("track_Hvtx_eta",&track_Hvtx_eta_float);
	tree->Branch("track_Hvtx_phi",&track_Hvtx_phi_float);
	tree->Branch("track_Hvtx_dxycorr",&track_Hvtx_dxycorr_float);
	tree->Branch("track_Hvtx_pterr",&track_Hvtx_pterr_float);
	tree->Branch("track_Hbs_pt",&track_Hbs_pt_float);
	tree->Branch("track_Hbs_eta",&track_Hbs_eta_float);
	tree->Branch("track_Hbs_phi",&track_Hbs_phi_float);
	tree->Branch("track_Hbs_dxycorr",&track_Hbs_dxycorr_float);
	tree->Branch("track_Hbs_pterr",&track_Hbs_pterr_float);
	
	tree->Branch("mu_pt_roc",&mu_pt_roc);
	tree->Branch("mu_eta_roc",&mu_eta_roc);
	tree->Branch("mu_phi_roc",&mu_phi_roc);
	tree->Branch("mu_gen_pt_roc",&mu_gen_pt_roc);
	tree->Branch("mu_genparticle_roc",&mu_genparticle_roc);
	tree->Branch("mu_nl_roc",&mu_nl_roc);
	tree->Branch("mu_q_roc",&mu_q_roc);
	tree->Branch("mu_tracktype_roc",&mu_tracktype_roc);
	

}

void UFHZZ4LAna::setTreeVariables( const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                   std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons, 
                                   std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons, 
                                   std::vector<pat::Jet> goodJets, std::vector<float> goodJetQGTagger, 
                                   std::vector<float> goodJetaxis2, std::vector<float> goodJetptD, std::vector<int> goodJetmult,
                                   std::vector<pat::Jet> selectedMergedJets)
{

    using namespace edm;
    using namespace pat;
    using namespace std;

    // Jet Info
    double tempDeltaR = 999.0;
    for( unsigned int k = 0; k < goodJets.size(); k++) {

        if (verbose) cout<<"jet pt: "<<goodJets[k].pt()<<" eta: "<<goodJets[k].eta()<<" phi: "<<goodJets[k].phi()<<endl;

        bool isclean_H4l = true;

        // check overlap with tight ID isolated leptons OR higgs candidate leptons
        unsigned int Nleptons = lep_pt.size();
        for (unsigned int i=0; i<Nleptons; i++) {            
            bool passed_idiso=true;
            if (abs(lep_id[i])==13 && lep_RelIsoNoFSR[i]>isoCutMu) passed_idiso=false;
            if (abs(lep_id[i])==11 && lep_RelIsoNoFSR[i]>isoCutEl) passed_idiso=false;
            if (!(lep_tightId[i])) passed_idiso=false;
            bool cand_lep=false;
            for (int l=0; l<4; l++) {
                if ((int)i==lep_Hindex[l]) cand_lep=true;
            }
            if (!(passed_idiso || cand_lep)) continue;
            TLorentzVector thisLep;
            thisLep.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
            tempDeltaR=999.0;
            tempDeltaR=deltaR(goodJets[k].eta(),goodJets[k].phi(),thisLep.Eta(),thisLep.Phi());
            if (verbose) cout<<" jet DeltaR between Lep"<<i<<" : "<<tempDeltaR;
            if (tempDeltaR<0.4) {
                isclean_H4l = false;
            }
        }
        if (verbose) cout<<endl;

        // check overlap with fsr photons
        unsigned int N = fsrPhotons_pt.size();
        for(unsigned int i=0; i<N; i++) {
            
            // don't clean jet from fsr if the photon wasn't matched to tight Id and Isolated lepton
            if (!lep_tightId[fsrPhotons_lepindex[i]]) continue;
            double RelIsoNoFSR=lep_RelIsoNoFSR[fsrPhotons_lepindex[i]];
            if (RelIsoNoFSR>((abs(lep_id[fsrPhotons_lepindex[i]])==11) ? isoCutEl : isoCutMu)) continue;

            TLorentzVector pho;
            pho.SetPtEtaPhiM(fsrPhotons_pt[i],fsrPhotons_eta[i],fsrPhotons_phi[i],0.0);
            tempDeltaR=999.0;
            tempDeltaR=deltaR(goodJets[k].eta(),goodJets[k].phi(),pho.Eta(),pho.Phi());
            if (verbose) cout<<" jet DeltaR between pho"<<i<<" : "<<tempDeltaR;
            if (tempDeltaR<0.4) {
                isclean_H4l = false;
            }
        }
        if (N>0&&verbose) cout<<endl;        

        //JER from database
        JME::JetParameters parameters;
        parameters.setJetPt(goodJets[k].pt());
        parameters.setJetEta(goodJets[k].eta());
        parameters.setRho(muRho);
        float relpterr = resolution_pt.getResolution(parameters);
        float phierr = resolution_phi.getResolution(parameters);

        double jercorr = 1.0; double jercorrup = 1.0; double jercorrdn = 1.0;
        
        if (isMC && doJER) {
            JME::JetParameters sf_parameters = {{JME::Binning::JetPt, goodJets[k].pt()}, {JME::Binning::JetEta, goodJets[k].eta()}, {JME::Binning::Rho, muRho}};
            float factor = resolution_sf.getScaleFactor(sf_parameters);
            float factorup = resolution_sf.getScaleFactor(sf_parameters, Variation::UP);
            float factordn = resolution_sf.getScaleFactor(sf_parameters, Variation::DOWN);
            
            double pt_jer, pt_jerup, pt_jerdn;
            const reco::GenJet * genJet = goodJets[k].genJet();
            if (genJet && deltaR(goodJets[k].eta(),goodJets[k].phi(),genJet->eta(),genJet->phi())<0.2 
                && (abs(goodJets[k].pt()-genJet->pt())<3*relpterr*goodJets[k].pt())) {
                double gen_pt = genJet->pt();
                pt_jer = max(0.0,gen_pt+factor*(goodJets[k].pt()-gen_pt));
                pt_jerup = max(0.0,gen_pt+factorup*(goodJets[k].pt()-gen_pt));
                pt_jerdn = max(0.0,gen_pt+factordn*(goodJets[k].pt()-gen_pt));
            } else {
                TRandom3 rand;
                rand.SetSeed(abs(static_cast<int>(sin(goodJets[k].phi())*100000)));
                float smear = rand.Gaus(0,1.0);
                float sigma = sqrt(factor*factor-1.0)*relpterr*goodJets[k].pt();
                float sigmaup = sqrt(factorup*factorup-1.0)*relpterr*goodJets[k].pt();
                float sigmadn = sqrt(factordn*factordn-1.0)*relpterr*goodJets[k].pt();
                pt_jer = max(0.0,smear*sigma+goodJets[k].pt());
                pt_jerup = max(0.0,smear*sigmaup+goodJets[k].pt());
                pt_jerdn = max(0.0,smear*sigmadn+goodJets[k].pt());
            }
            
            jercorr = pt_jer/goodJets[k].pt();
            jercorrup = pt_jerup/goodJets[k].pt();
            jercorrdn = pt_jerdn/goodJets[k].pt();
        }

        TLorentzVector *jet_jer = new TLorentzVector(jercorr*goodJets[k].px(),jercorr*goodJets[k].py(),jercorr*goodJets[k].pz(),jercorr*goodJets[k].energy());
        TLorentzVector *jet_jerup = new TLorentzVector(jercorrup*goodJets[k].px(),jercorrup*goodJets[k].py(),jercorrup*goodJets[k].pz(),jercorrup*goodJets[k].energy());
        TLorentzVector *jet_jerdn = new TLorentzVector(jercorrdn*goodJets[k].px(),jercorrdn*goodJets[k].py(),jercorrdn*goodJets[k].pz(),jercorrdn*goodJets[k].energy());


        bool passPU_;
        if (doJEC && (year==2017 || year==2018)) {
            passPU_ = bool(goodJets[k].userInt("pileupJetIdUpdated:fullId") & (1 << 0));
        } else {
            passPU_ = bool(goodJets[k].userInt("pileupJetId:fullId") & (1 << 0));
        }
        if(!(passPU_ || !doPUJetID || jet_jer->Pt()>50)) continue;
        if(jet_jer->Pt()<10) continue;

        if (verbose) cout<<"Jet nominal: "<<goodJets[k].pt()<<" JER corrected: "<<jet_jer->Pt()<<" JER up: "<<jet_jerup->Pt()<<" JER dn: "<<jet_jerdn->Pt()<<" check Delta between jet and lep / pho: "<<isclean_H4l<<std::endl;

        jecunc->setJetPt(jet_jer->Pt());
        jecunc->setJetEta(goodJets[k].eta());
        double jecunc_up = 1.0+jecunc->getUncertainty(true);
        jecunc->setJetPt(jet_jer->Pt());
        jecunc->setJetEta(goodJets[k].eta());
        double jecunc_dn = 1.0-jecunc->getUncertainty(false);

        if (jet_jer->Pt() > 30.0 && fabs(goodJets[k].eta())<4.7) {
            if (isclean_H4l) { 
                njets_pt30_eta4p7++;
                jet_iscleanH4l.push_back((int)jet_pt.size());
                if (jet_jer->Pt() > pt_leadingjet_pt30_eta4p7) {
                    pt_leadingjet_pt30_eta4p7 = jet_jer->Pt();
                    absrapidity_leadingjet_pt30_eta4p7 = jet_jer->Rapidity(); //take abs later
                }
                if (jet_jer->Pt()>jet1pt) {
                    jet2pt=jet1pt; jet2index=jet1index;
                    jet1pt=jet_jer->Pt(); jet1index=(int)jet_pt.size();
                } else if (jet_jer->Pt()>jet2pt) {
                    jet2pt=jet_jer->Pt(); jet2index=(int)jet_pt.size();
                }
                if (fabs(goodJets[k].eta())<2.5) {
                    njets_pt30_eta2p5++;
                    if (jet_jer->Pt() > pt_leadingjet_pt30_eta2p5) {
                        pt_leadingjet_pt30_eta2p5 = jet_jer->Pt();
                    }
                }

            }
            jet_pt.push_back(jet_jer->Pt());
            jet_pt_raw.push_back(goodJets[k].correctedJet("Uncorrected").pt());///jet Pt without JEC applied
            jet_eta.push_back(jet_jer->Eta());
            jet_phi.push_back(jet_jer->Phi());
            jet_mass.push_back(jet_jer->M());
            if (doJEC && (year==2017 || year==2018)) {
                jet_pumva.push_back(goodJets[k].userFloat("pileupJetIdUpdated:fullDiscriminant"));
            } else {
                jet_pumva.push_back(goodJets[k].userFloat("pileupJetId:fullDiscriminant"));
            }
            
            //jet_csvv2.push_back(goodJets[k].bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll"));
            jet_csvv2.push_back(goodJets[k].bDiscriminator("pfDeepCSVJetTags:probb")+goodJets[k].bDiscriminator("pfDeepCSVJetTags:probbb"));
            jet_csvv2_.push_back(goodJets[k].bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll"));


            TRandom3 rand;
            rand.SetSeed(abs(static_cast<int>(sin(goodJets[k].phi())*100000)));
            float coin = rand.Uniform(1.); 

            double jet_scalefactor    = 1.0;
            
            jet_scalefactor    = reader->eval_auto_bounds(
                "central", 
                BTagEntry::FLAV_B, 
                goodJets[k].eta(), 
                goodJets[k].pt()
                ); 
            

            //if (goodJets[k].bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll")>BTagCut && coin>(1.0-jet_scalefactor)) {
            if ((goodJets[k].bDiscriminator("pfDeepCSVJetTags:probb")+goodJets[k].bDiscriminator("pfDeepCSVJetTags:probbb"))>BTagCut && coin>(1.0-jet_scalefactor)) {
                jet_isbtag.push_back(1);
            } else {
                jet_isbtag.push_back(0);
            }

            //if (goodJets[k].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>BTagCut && isclean_H4l) nbjets_pt30_eta4p7++;
            if ((goodJets[k].bDiscriminator("pfDeepCSVJetTags:probb")+goodJets[k].bDiscriminator("pfDeepCSVJetTags:probbb"))>BTagCut && isclean_H4l) nbjets_pt30_eta4p7++;
            if (isMC) {
                jet_hadronFlavour.push_back(goodJets[k].hadronFlavour());
                jet_partonFlavour.push_back(goodJets[k].partonFlavour());
            } else {
                jet_hadronFlavour.push_back(-1);
                jet_partonFlavour.push_back(-1);
            }
            jet_QGTagger.push_back(goodJetQGTagger[k]);
            jet_axis2.push_back(goodJetaxis2[k]);
            jet_ptD.push_back(goodJetptD[k]);
            jet_mult.push_back(goodJetmult[k]);
            jet_relpterr.push_back(relpterr);
            jet_phierr.push_back(phierr);
            ///jet_bTagEffi.push_back(helper.get_bTagEffi(goodJets[k].pt(), goodJets[k].eta(), hbTagEffi));
            jet_bTagEffi.push_back(helper.get_bTagEffi(jet_jer->Pt(), jet_jer->Eta(), hbTagEffi));
            jet_cTagEffi.push_back(helper.get_bTagEffi(jet_jer->Pt(), jet_jer->Eta(), hcTagEffi));
            jet_udsgTagEffi.push_back(helper.get_bTagEffi(jet_jer->Pt(), jet_jer->Eta(), hudsgTagEffi));
            //if (goodJets[k].bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll")>BTagCut && coin>(1.0-jet_scalefactor)) jet_bTagEffi.push_back(helper.get_bTagEffi(goodJets[k].pt(), goodJets[k].eta(), hbTagEffi));
            //else jet_bTagEffi.push_back(-1);
            //if (goodJets[k].bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll")>BTagCut && coin>(1.0-jet_scalefactor)) jet_bTagEffi_.push_back(helper.get_bTagEffi(jet_jer->Pt(), jet_jer->Eta(), hbTagEffi));
            //else jet_bTagEffi_.push_back(-1);
        }
        
        // JER up
        if (jet_jerup->Pt() > 30.0 && fabs(jet_jerup->Eta())<4.7) {
            if (isclean_H4l) {
                njets_pt30_eta4p7_jerup++;
                jet_jerup_iscleanH4l.push_back((int)jet_jerup_pt.size());
                if (jet_jerup->Pt() > pt_leadingjet_pt30_eta4p7_jerup) {
                    pt_leadingjet_pt30_eta4p7_jerup = jet_jerup->Pt();
                    absrapidity_leadingjet_pt30_eta4p7_jerup = jet_jerup->Rapidity(); //take abs later
                }
                if (fabs(jet_jerup->Eta())<2.5) {
                    njets_pt30_eta2p5_jerup++;
                    if (jet_jerup->Pt() > pt_leadingjet_pt30_eta2p5_jerup) {
                        pt_leadingjet_pt30_eta2p5_jerup = jet_jerup->Pt();
                    }
                }

            }
            jet_jerup_pt.push_back(jet_jerup->Pt());
            jet_jerup_eta.push_back(jet_jerup->Eta());
            jet_jerup_phi.push_back(jet_jerup->Phi());
            jet_jerup_mass.push_back(jet_jerup->M());            
        }

        // JER dn
        if (jet_jerdn->Pt() > 30.0 && fabs(jet_jerdn->Eta())<4.7) {
            if (isclean_H4l) {
                njets_pt30_eta4p7_jerdn++;
                jet_jerdn_iscleanH4l.push_back((int)jet_jerdn_pt.size());
                if (jet_jerdn->Pt() > pt_leadingjet_pt30_eta4p7_jerdn) {
                    pt_leadingjet_pt30_eta4p7_jerdn = jet_jerdn->Pt();
                    absrapidity_leadingjet_pt30_eta4p7_jerdn = jet_jerdn->Rapidity(); //take abs later
                }
                if (fabs(jet_jerdn->Eta())<2.5) {
                    njets_pt30_eta2p5_jerdn++;
                    if (jet_jerdn->Pt() > pt_leadingjet_pt30_eta2p5_jerdn) {
                        pt_leadingjet_pt30_eta2p5_jerdn = jet_jerdn->Pt();
                    }
                }

            }
            jet_jerdn_pt.push_back(jet_jerdn->Pt());
            jet_jerdn_eta.push_back(jet_jerdn->Eta());
            jet_jerdn_phi.push_back(jet_jerdn->Phi());
            jet_jerdn_mass.push_back(jet_jerdn->M());
        }

        double jetPx_jesup = jecunc_up * jet_jer->Px();
        double jetPy_jesup = jecunc_up * jet_jer->Py();
        double jetPz_jesup = jecunc_up * jet_jer->Pz();
        double jetE_jesup = sqrt(jetPx_jesup*jetPx_jesup + jetPy_jesup*jetPy_jesup + jetPz_jesup*jetPz_jesup + jet_jer->M()*jet_jer->M());
        TLorentzVector *jet_jesup = new TLorentzVector(jetPx_jesup,jetPy_jesup,jetPz_jesup,jetE_jesup);
            
        if (jet_jesup->Pt() > 30.0 && fabs(jet_jesup->Eta())<4.7) {
            if (isclean_H4l) {
                njets_pt30_eta4p7_jesup++;
                jet_jesup_iscleanH4l.push_back((int)jet_jesup_pt.size());
                if (jet_jesup->Pt() > pt_leadingjet_pt30_eta4p7_jesup) {
                    pt_leadingjet_pt30_eta4p7_jesup = jet_jesup->Pt();
                    absrapidity_leadingjet_pt30_eta4p7_jesup = jet_jesup->Rapidity(); //take abs later
                }                
                if (fabs(jet_jesup->Eta())<2.5) {
                    njets_pt30_eta2p5_jesup++;
                    if (jet_jesup->Pt() > pt_leadingjet_pt30_eta2p5_jesup) {
                        pt_leadingjet_pt30_eta2p5_jesup = jet_jesup->Pt();
                    }                
                }
             
            }
            TLorentzVector jet_jesup(jetPx_jesup,jetPy_jesup,jetPz_jesup,jetE_jesup);
            jet_jesup_pt.push_back(jet_jesup.Pt());
            jet_jesup_eta.push_back(jet_jesup.Eta());
            jet_jesup_phi.push_back(jet_jesup.Phi());
            jet_jesup_mass.push_back(jet_jesup.M());
            jet_QGTagger_jesup.push_back(goodJetQGTagger[k]);
        }

        double jetPx_jesdn = jecunc_dn * jet_jer->Px();
        double jetPy_jesdn = jecunc_dn * jet_jer->Py();
        double jetPz_jesdn = jecunc_dn * jet_jer->Pz();
        double jetE_jesdn = sqrt(jetPx_jesdn*jetPx_jesdn + jetPy_jesdn*jetPy_jesdn + jetPz_jesdn*jetPz_jesdn + jet_jer->M()*jet_jer->M());
        TLorentzVector *jet_jesdn = new TLorentzVector(jetPx_jesdn,jetPy_jesdn,jetPz_jesdn,jetE_jesdn);

        if (jet_jesdn->Pt() > 30.0 && fabs(jet_jesdn->Eta())<4.7) {
            if (isclean_H4l) {
                njets_pt30_eta4p7_jesdn++;
                jet_jesdn_iscleanH4l.push_back(jet_jesdn_pt.size());
                if (jet_jesdn->Pt() > pt_leadingjet_pt30_eta4p7_jesdn) {
                    pt_leadingjet_pt30_eta4p7_jesdn = jet_jesdn->Pt();
                    absrapidity_leadingjet_pt30_eta4p7_jesdn = jet_jesdn->Rapidity(); //take abs later
                }        
                if (fabs(jet_jesdn->Eta())<2.5) {
                    njets_pt30_eta2p5_jesdn++;
                    if (jet_jesdn->Pt() > pt_leadingjet_pt30_eta2p5_jesdn) {
                        pt_leadingjet_pt30_eta2p5_jesdn = jet_jesdn->Pt();
                    }                
                }            
            }
            TLorentzVector jet_jesdn(jetPx_jesdn,jetPy_jesdn,jetPz_jesdn,jetE_jesdn);
            jet_jesdn_pt.push_back(jet_jesdn.Pt());
            jet_jesdn_eta.push_back(jet_jesdn.Eta());
            jet_jesdn_phi.push_back(jet_jesdn.Phi());
            jet_jesdn_mass.push_back(jet_jesdn.M());
            jet_QGTagger_jesdn.push_back(goodJetQGTagger[k]);

        }

       
    } // loop over jets

    /*
    // merged jet
    for( unsigned int k = 0; k < selectedMergedJets.size(); k++) {
        
        double tempDeltaR = 999.0;
        bool isclean_H4l = true;
        
        unsigned int Nleptons = lep_pt.size();
        for (unsigned int i=0; i<Nleptons; i++) {
            
            if (abs(lep_id[i])==13 && lep_RelIsoNoFSR[i]>isoCutMu) continue;
            if (abs(lep_id[i])==11 && lep_RelIsoNoFSR[i]>isoCutEl) continue;
            if (!(lep_tightId[i])) continue;
            TLorentzVector thisLep;
            thisLep.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
            tempDeltaR=999.0;
            tempDeltaR=deltaR(selectedMergedJets[k].eta(),selectedMergedJets[k].phi(),thisLep.Eta(),thisLep.Phi());
            if (tempDeltaR<0.8) {
                isclean_H4l = false;
            }
        }
             
        // check overlap with fsr photons  
        unsigned int N = fsrPhotons_pt.size();
        for(unsigned int i=0; i<N; i++) {

            // don't clean jet from fsr if the photon wasn't matched to tight Id and Isolated lepton                                                                                                                                    
            if (!lep_tightId[fsrPhotons_lepindex[i]]) continue;
            double RelIsoNoFSR=lep_RelIsoNoFSR[fsrPhotons_lepindex[i]];
            if (RelIsoNoFSR>((abs(lep_id[fsrPhotons_lepindex[i]])==11) ? isoCutEl : isoCutMu)) continue;

            TLorentzVector pho;
            pho.SetPtEtaPhiM(fsrPhotons_pt[i],fsrPhotons_eta[i],fsrPhotons_phi[i],0.0);
            tempDeltaR=999.0;
            tempDeltaR=deltaR(selectedMergedJets[k].eta(),selectedMergedJets[k].phi(),pho.Eta(),pho.Phi());
            if (tempDeltaR<0.8) {
                isclean_H4l = false;
            }
       } 

        if (isclean_H4l) mergedjet_iscleanH4l.push_back((int)mergedjet_pt.size());
        mergedjet_pt.push_back((float)selectedMergedJets[k].pt());
        mergedjet_eta.push_back((float)selectedMergedJets[k].eta());             
        mergedjet_phi.push_back((float)selectedMergedJets[k].phi());
        mergedjet_mass.push_back((float)selectedMergedJets[k].mass());            
        mergedjet_L1.push_back((float)selectedMergedJets[k].jecFactor("L1FastJet")); // current JEC to L1

        mergedjet_softdropmass.push_back((float)selectedMergedJets[k].userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass"));
        mergedjet_prunedmass.push_back((float)selectedMergedJets[k].userFloat("ak8PFJetsCHSCorrPrunedMass"));
        mergedjet_tau1.push_back((float)selectedMergedJets[k].userFloat("NjettinessAK8Puppi:tau1") );
        mergedjet_tau2.push_back((float)selectedMergedJets[k].userFloat("NjettinessAK8Puppi:tau2") );
        mergedjet_btag.push_back((float)selectedMergedJets[k].bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") );

        if (verbose) cout<<"double btag: "<<selectedMergedJets[k].bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags")<<endl;
        
        auto wSubjets = selectedMergedJets[k].subjets("SoftDrop");
        int nsub = 0;
        vector<float> subjets_pt, subjets_eta, subjets_phi, subjets_mass, subjets_btag; 
        vector<int> subjets_hadronFlavour, subjets_partonFlavour;
        for ( auto const & iw : wSubjets ) {                        
            nsub = nsub + 1;            
            subjets_pt.push_back((float)iw->pt());
            subjets_eta.push_back((float)iw->eta());
            subjets_phi.push_back((float)iw->phi());
            subjets_mass.push_back((float)iw->mass());
            subjets_btag.push_back((float)iw->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
            subjets_hadronFlavour.push_back(iw->hadronFlavour());
            subjets_partonFlavour.push_back(iw->partonFlavour());
            if (verbose) cout<<"subjet parton "<<iw->partonFlavour()<<" hadron "<<iw->hadronFlavour()<<endl;

        }        

        mergedjet_nsubjet.push_back(nsub);
        mergedjet_subjet_pt.push_back(subjets_pt);
        mergedjet_subjet_eta.push_back(subjets_eta);
        mergedjet_subjet_phi.push_back(subjets_phi);
        mergedjet_subjet_mass.push_back(subjets_mass);
        mergedjet_subjet_btag.push_back(subjets_btag);        
        mergedjet_subjet_partonFlavour.push_back(subjets_partonFlavour);
        mergedjet_subjet_hadronFlavour.push_back(subjets_hadronFlavour);
            
    }
    */

    // Higgs Variables
    if( RecoFourMuEvent ){ finalState = 1;}
    if( RecoFourEEvent  ){ finalState = 2;}
    if( RecoTwoETwoMuEvent ){ finalState = 3;}
    if( RecoTwoMuTwoEEvent ){ finalState = 4;}

    H_pt.push_back(HVec.Pt());
    H_eta.push_back(HVec.Eta());
    H_phi.push_back(HVec.Phi());
    H_mass.push_back(HVec.M());
    mass4l = HVec.M();
    mass4l_noFSR = HVecNoFSR.M();

    if(RecoFourMuEvent){mass4mu = HVec.M();}
    else{mass4mu = -1;}
    if(RecoFourEEvent){ mass4e = HVec.M();}
    else{ mass4e = -1;}
    if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){ mass2e2mu = HVec.M(); }
    else{ mass2e2mu = -1;}

    pT4l = HVec.Pt(); eta4l = HVec.Eta(); rapidity4l = HVec.Rapidity(); phi4l = HVec.Phi();

    pTZ1 = Z1Vec.Pt(); pTZ2 = Z2Vec.Pt(); massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M();

    if (njets_pt30_eta4p7>0) absdeltarapidity_hleadingjet_pt30_eta4p7 = fabs(rapidity4l-absrapidity_leadingjet_pt30_eta4p7);
    if (njets_pt30_eta4p7_jesup>0) absdeltarapidity_hleadingjet_pt30_eta4p7_jesup = fabs(rapidity4l-absrapidity_leadingjet_pt30_eta4p7_jesup);
    if (njets_pt30_eta4p7_jesdn>0) absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn = fabs(rapidity4l-absrapidity_leadingjet_pt30_eta4p7_jesdn);
    if (njets_pt30_eta4p7_jerup>0) absdeltarapidity_hleadingjet_pt30_eta4p7_jerup = fabs(rapidity4l-absrapidity_leadingjet_pt30_eta4p7_jerup);
    if (njets_pt30_eta4p7_jerdn>0) absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn = fabs(rapidity4l-absrapidity_leadingjet_pt30_eta4p7_jerdn);
    
    if (njets_pt30_eta4p7>0) absrapidity_leadingjet_pt30_eta4p7 = fabs(absrapidity_leadingjet_pt30_eta4p7);
    if (njets_pt30_eta4p7_jesup>0) absrapidity_leadingjet_pt30_eta4p7_jesup = fabs(absrapidity_leadingjet_pt30_eta4p7_jesup);
    if (njets_pt30_eta4p7_jesdn>0) absrapidity_leadingjet_pt30_eta4p7_jesdn = fabs(absrapidity_leadingjet_pt30_eta4p7_jesdn);
    if (njets_pt30_eta4p7_jerup>0) absrapidity_leadingjet_pt30_eta4p7_jerup = fabs(absrapidity_leadingjet_pt30_eta4p7_jerup);
    if (njets_pt30_eta4p7_jerdn>0) absrapidity_leadingjet_pt30_eta4p7_jerdn = fabs(absrapidity_leadingjet_pt30_eta4p7_jerdn);

    if (foundHiggsCandidate) {

        TLorentzVector Lep1FSR, Lep2FSR, Lep3FSR, Lep4FSR;
        Lep1FSR.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[0]],lepFSR_eta[lep_Hindex[0]],lepFSR_phi[lep_Hindex[0]],lepFSR_mass[lep_Hindex[0]]);
        Lep2FSR.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[1]],lepFSR_eta[lep_Hindex[1]],lepFSR_phi[lep_Hindex[1]],lepFSR_mass[lep_Hindex[1]]);
        Lep3FSR.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[2]],lepFSR_eta[lep_Hindex[2]],lepFSR_phi[lep_Hindex[2]],lepFSR_mass[lep_Hindex[2]]);
        Lep4FSR.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[3]],lepFSR_eta[lep_Hindex[3]],lepFSR_phi[lep_Hindex[3]],lepFSR_mass[lep_Hindex[3]]);
        pTL1FSR = Lep1FSR.Pt(); pTL2FSR = Lep2FSR.Pt(); pTL3FSR = Lep3FSR.Pt(); pTL4FSR = Lep4FSR.Pt();

        TLorentzVector Lep1, Lep2, Lep3, Lep4;
        Lep1.SetPtEtaPhiM(lep_pt[lep_Hindex[0]],lep_eta[lep_Hindex[0]],lep_phi[lep_Hindex[0]],lep_mass[lep_Hindex[0]]);
        Lep2.SetPtEtaPhiM(lep_pt[lep_Hindex[1]],lep_eta[lep_Hindex[1]],lep_phi[lep_Hindex[1]],lep_mass[lep_Hindex[1]]);
        Lep3.SetPtEtaPhiM(lep_pt[lep_Hindex[2]],lep_eta[lep_Hindex[2]],lep_phi[lep_Hindex[2]],lep_mass[lep_Hindex[2]]);
        Lep4.SetPtEtaPhiM(lep_pt[lep_Hindex[3]],lep_eta[lep_Hindex[3]],lep_phi[lep_Hindex[3]],lep_mass[lep_Hindex[3]]);

        //idL1 = lep_id[lep_Hindex[0]]; pTL1 = Lep1.Pt(); etaL1 = Lep1.Eta();
        //idL2 = lep_id[lep_Hindex[1]]; pTL2 = Lep2.Pt(); etaL2 = Lep2.Eta();
        //idL3 = lep_id[lep_Hindex[2]]; pTL3 = Lep3.Pt(); etaL3 = Lep3.Eta();       
        //idL4 = lep_id[lep_Hindex[3]]; pTL4 = Lep4.Pt(); etaL4 = Lep4.Eta();
        idL1 = lep_id[lep_Hindex[0]]; pTL1 = Lep1.Pt(); etaL1 = Lep1.Eta(); pTErrL1 = lep_pterr[lep_Hindex[0]]; mL1 = lep_mass[lep_Hindex[0]]; phiL1 = lep_phi[lep_Hindex[0]];
        idL2 = lep_id[lep_Hindex[1]]; pTL2 = Lep2.Pt(); etaL2 = Lep2.Eta(); pTErrL2 = lep_pterr[lep_Hindex[1]]; mL2 = lep_mass[lep_Hindex[1]]; phiL2 = lep_phi[lep_Hindex[1]];
        idL3 = lep_id[lep_Hindex[2]]; pTL3 = Lep3.Pt(); etaL3 = Lep3.Eta(); pTErrL3 = lep_pterr[lep_Hindex[2]]; mL3 = lep_mass[lep_Hindex[2]]; phiL3 = lep_phi[lep_Hindex[2]];
        idL4 = lep_id[lep_Hindex[3]]; pTL4 = Lep4.Pt(); etaL4 = Lep4.Eta(); pTErrL4 = lep_pterr[lep_Hindex[3]]; mL4 = lep_mass[lep_Hindex[3]]; phiL4 = lep_phi[lep_Hindex[3]];

    }


}




void UFHZZ4LAna::setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles,
                                 edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles,
                                 edm::Handle<edm::View<reco::GenJet> > genJets)
{

    reco::GenParticleCollection::const_iterator genPart;
    int j = -1;
    int nGENLeptons=0;

    if (verbose) cout<<"begin looping on gen particles"<<endl;
    for(genPart = prunedgenParticles->begin(); genPart != prunedgenParticles->end(); genPart++) {
        j++;

        if (abs(genPart->pdgId())==11  || abs(genPart->pdgId())==13 || abs(genPart->pdgId())==15) {

            if (!(genPart->status()==1 || abs(genPart->pdgId())==15)) continue;
            if (!(genAna.MotherID(&prunedgenParticles->at(j))==23 || abs(genAna.MotherID(&prunedgenParticles->at(j)))==24) ) continue;
            
            nGENLeptons++;
            if (verbose) cout<<"found a gen lepton: id "<<genPart->pdgId()<<" pt: "<<genPart->pt()<<" eta: "<<genPart->eta()<<" status: "<<genPart->status()<<endl;

            // Collect FSR photons
            TLorentzVector lep_dressed;            
            lep_dressed.SetPtEtaPhiE(genPart->pt(),genPart->eta(),genPart->phi(),genPart->energy());
            set<int> gen_fsrset;
            for(size_t k=0; k<packedgenParticles->size();k++){
                if( (*packedgenParticles)[k].status() != 1) continue; // stable particles only
                if( (*packedgenParticles)[k].pdgId() != 22) continue; // only photons
                double this_dR_lgamma = deltaR(genPart->eta(), genPart->phi(), (*packedgenParticles)[k].eta(), (*packedgenParticles)[k].phi());
                bool idmatch=false;
                if ((*packedgenParticles)[k].mother(0)->pdgId()==genPart->pdgId() ) idmatch=true;
                const reco::Candidate * mother = (*packedgenParticles)[k].mother(0);
                for(size_t m=0;m<mother->numberOfMothers();m++) {
                    if ( (*packedgenParticles)[k].mother(m)->pdgId() == genPart->pdgId() ) idmatch=true;
                }
                if (!idmatch) continue;                
                if(this_dR_lgamma<((abs(genPart->pdgId())==11)?genIsoConeSizeEl:genIsoConeSizeMu)) {
                    gen_fsrset.insert(k);
                    TLorentzVector gamma;
                    gamma.SetPtEtaPhiE((*packedgenParticles)[k].pt(),(*packedgenParticles)[k].eta(),(*packedgenParticles)[k].phi(),(*packedgenParticles)[k].energy());
                    lep_dressed = lep_dressed+gamma;
                }
            } // Dressed leptons loop
            if (verbose) cout<<"gen lep pt "<<genPart->pt()<< " dressed pt: " << lep_dressed.Pt()<<endl;
  
            GENlep_id.push_back( genPart->pdgId() );
            GENlep_status.push_back(genPart->status());
            GENlep_pt.push_back( lep_dressed.Pt() );
            GENlep_eta.push_back( lep_dressed.Eta() );
            GENlep_phi.push_back( lep_dressed.Phi() );
            GENlep_mass.push_back( lep_dressed.M() );
            GENlep_MomId.push_back(genAna.MotherID(&prunedgenParticles->at(j)));
            GENlep_MomMomId.push_back(genAna.MotherMotherID(&prunedgenParticles->at(j)));
       
            TLorentzVector thisLep;
            thisLep.SetPtEtaPhiM(lep_dressed.Pt(),lep_dressed.Eta(),lep_dressed.Phi(),lep_dressed.M());
            // GEN iso calculation
            double this_GENiso=0.0;
            double this_GENneutraliso=0.0;
            double this_GENchargediso=0.0;
            if (verbose) cout<<"gen iso calculation"<<endl;
            for(size_t k=0; k<packedgenParticles->size();k++){
                if( (*packedgenParticles)[k].status() != 1 ) continue; // stable particles only         
                if (abs((*packedgenParticles)[k].pdgId())==12 || abs((*packedgenParticles)[k].pdgId())==14 || abs((*packedgenParticles)[k].pdgId())==16) continue; // exclude neutrinos
                if ((abs((*packedgenParticles)[k].pdgId())==11 || abs((*packedgenParticles)[k].pdgId())==13)) continue; // exclude leptons
                if (gen_fsrset.find(k)!=gen_fsrset.end()) continue; // exclude particles which were selected as fsr photons
                double this_dRvL = deltaR(thisLep.Eta(), thisLep.Phi(), (*packedgenParticles)[k].eta(), (*packedgenParticles)[k].phi());
                if(this_dRvL<((abs(genPart->pdgId())==11)?genIsoConeSizeEl:genIsoConeSizeMu)) {
                    if (verbose) cout<<"adding to geniso id: "<<(*packedgenParticles)[k].pdgId()<<" status: "<<(*packedgenParticles)[k].status()<<" pt: "<<(*packedgenParticles)[k].pt()<<" dR: "<<this_dRvL<<endl;
                    this_GENiso = this_GENiso + (*packedgenParticles)[k].pt();
                    if ((*packedgenParticles)[k].charge()==0) this_GENneutraliso = this_GENneutraliso + (*packedgenParticles)[k].pt();
                    if ((*packedgenParticles)[k].charge()!=0) this_GENchargediso = this_GENchargediso + (*packedgenParticles)[k].pt();
                }
            } // GEN iso loop
            this_GENiso = this_GENiso/thisLep.Pt();
            if (verbose) cout<<"gen lep pt: "<<thisLep.Pt()<<" rel iso: "<<this_GENiso<<endl;
            GENlep_RelIso.push_back(this_GENiso);
            // END GEN iso calculation

        } // leptons
        
        if (genPart->pdgId()==25) {
            GENMH=genPart->mass();
            GENH_pt.push_back(genPart->pt());
            GENH_eta.push_back(genPart->eta());
            GENH_phi.push_back(genPart->phi());
            GENH_mass.push_back(genPart->mass());
        }

        
        if (genPart->pdgId()==23 && (genPart->status()>=20 && genPart->status()<30) ) {
            const reco::Candidate *Zdau0=genPart->daughter(0);
            int ZdauId = fabs(Zdau0->pdgId());
            if (fabs(Zdau0->pdgId())==23) {
                int ndau = genPart->numberOfDaughters();
                for (int d=0; d<ndau; d++) {
                    const reco::Candidate *Zdau=genPart->daughter(d);
                    if (verbose) cout<<"ZDau "<<d<<" id "<<fabs(Zdau->pdgId())<<endl;
                    if (fabs(Zdau->pdgId())<17) {
                        ZdauId = fabs(Zdau->pdgId());
                        break;
                    }
                }
            }
            if (verbose) cout<<"GENZ status "<<genPart->status()<<" MomId: "<<genAna.MotherID(&prunedgenParticles->at(j))<< "DauId: "<<ZdauId<<endl;

            if (Zdau0) GENZ_DaughtersId.push_back(ZdauId);           
            GENZ_MomId.push_back(genAna.MotherID(&prunedgenParticles->at(j)));                
            GENZ_pt.push_back(genPart->pt());
            GENZ_eta.push_back(genPart->eta());
            GENZ_phi.push_back(genPart->phi());
            GENZ_mass.push_back(genPart->mass());
        }

        if (abs(genPart->pdgId())>500 && abs(genPart->pdgId())<600 && genPart->status()==2) {
            nGenStatus2bHad+=1;
        }

    }

    if (GENlep_pt.size()>=4) {
        
        unsigned int L1_nocuts=99; unsigned int L2_nocuts=99; unsigned int L3_nocuts=99; unsigned int L4_nocuts=99;
        bool passedFiducialSelectionNoCuts = mZ1_mZ2(L1_nocuts, L2_nocuts, L3_nocuts, L4_nocuts, false);
        if (passedFiducialSelectionNoCuts) {
            TLorentzVector Z1_1, Z1_2, Z2_1, Z2_2;
            Z1_1.SetPtEtaPhiM(GENlep_pt[L1_nocuts],GENlep_eta[L1_nocuts],GENlep_phi[L1_nocuts],GENlep_mass[L1_nocuts]);
            Z1_2.SetPtEtaPhiM(GENlep_pt[L2_nocuts],GENlep_eta[L2_nocuts],GENlep_phi[L2_nocuts],GENlep_mass[L2_nocuts]);
            Z2_1.SetPtEtaPhiM(GENlep_pt[L3_nocuts],GENlep_eta[L3_nocuts],GENlep_phi[L3_nocuts],GENlep_mass[L3_nocuts]);
            Z2_2.SetPtEtaPhiM(GENlep_pt[L4_nocuts],GENlep_eta[L4_nocuts],GENlep_phi[L4_nocuts],GENlep_mass[L4_nocuts]);         
            GENdPhiZZ = deltaPhi((Z1_1+Z1_2).Phi(),(Z2_1+Z2_2).Phi());            
            GENmassZZ = (Z1_1+Z1_2+Z2_1+Z2_2).M();            
            GENpTZZ = (Z1_1+Z1_2+Z2_1+Z2_2).Pt();            
            int genfs;
            if (abs(GENlep_id[L1_nocuts])==abs(GENlep_id[L3_nocuts])) genfs=1;
            else genfs=2;
            k_qqZZ_qcd_dPhi = helper.kfactor_qqZZ_qcd_dPhi(GENdPhiZZ,genfs);
            k_qqZZ_qcd_M = helper.kfactor_qqZZ_qcd_M(GENmassZZ,genfs);
            k_qqZZ_qcd_Pt = helper.kfactor_qqZZ_qcd_Pt(GENpTZZ,genfs);
            if (verbose) cout<<"qcd kfactor qqZZ: "<<k_qqZZ_qcd_M<<endl;
        }
    }    
    
   
    /////// DO THE FIDUCIAL VOLUME CALCULATION //////////////
    if (verbose) cout<<"begin fiducial volume calculation"<<endl;
    passedFiducialSelection=false;
    int nFiducialLeptons = 0;
    int nFiducialPtLead=0;
    int nFiducialPtSublead=0;
    
    for (unsigned int i=0; i<GENlep_id.size(); ++i) {    

        TLorentzVector thisLep;
        thisLep.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
        
        if ( ( (abs(GENlep_id[i]) == 13 && thisLep.Pt() > 5.0 && abs(thisLep.Eta()) < 2.4)
               || (abs(GENlep_id[i]) == 11 && thisLep.Pt() > 7.0 && abs(thisLep.Eta()) < 2.5) )
             && GENlep_RelIso[i]<((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu) ) {
            nFiducialLeptons++;
            if (verbose) cout<<nFiducialLeptons<<" fiducial leptons, id;"<<GENlep_id[i]<<" status: "<<GENlep_status[i]<<" pt: "<<thisLep.Pt()<<" eta: "<<thisLep.Eta()<<endl; 
            if (thisLep.Pt()>leadingPtCut) nFiducialPtLead++;
            if (thisLep.Pt()>subleadingPtCut) nFiducialPtSublead++;
        }                
    }
    
    if (nFiducialLeptons>=4 && nFiducialPtLead>=1 && nFiducialPtSublead>=2) {
                                
        // START FIDUCIAL EVENT TOPOLOGY CUTS
        unsigned int L1=99; unsigned int L2=99; unsigned int L3=99; unsigned int L4=99;
        GENmass4l = -1.0; GENmass4e = -1.0; GENmass4mu = -1.0; GENmass2e2mu = -1.0;
        GENmassZ1 = -1.0; GENmassZ2 = -1.0; GENpT4l = -1.0; GENeta4l = 999.; GENrapidity4l = 999.;
        //cout<<"Run: "<<Run<<" LumiSect: "<<LumiSect<<" Event: "<<Event<<endl;
        passedFiducialSelection = mZ1_mZ2(L1, L2, L3, L4, true);      
        if (verbose) cout<<"passedFiducialSelection? "<<passedFiducialSelection<<endl;
        
        GENlep_Hindex[0] = L1; GENlep_Hindex[1] = L2; GENlep_Hindex[2] = L3; GENlep_Hindex[3] = L4;

        if (passedFiducialSelection) {
            
            TLorentzVector LS3_Z1_1, LS3_Z1_2, LS3_Z2_1, LS3_Z2_2;
            LS3_Z1_1.SetPtEtaPhiM(GENlep_pt[L1],GENlep_eta[L1],GENlep_phi[L1],GENlep_mass[L1]);
            LS3_Z1_2.SetPtEtaPhiM(GENlep_pt[L2],GENlep_eta[L2],GENlep_phi[L2],GENlep_mass[L2]);
            LS3_Z2_1.SetPtEtaPhiM(GENlep_pt[L3],GENlep_eta[L3],GENlep_phi[L3],GENlep_mass[L3]);
            LS3_Z2_2.SetPtEtaPhiM(GENlep_pt[L4],GENlep_eta[L4],GENlep_phi[L4],GENlep_mass[L4]);
            
            GENmass4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).M();

            if (abs(GENlep_id[L1])==11 && abs(GENlep_id[L3])==11) {GENmass4e = GENmass4l;};
            if (abs(GENlep_id[L1])==13 && abs(GENlep_id[L3])==13) {GENmass4mu = GENmass4l;};
            if ( (abs(GENlep_id[L1])==11 || abs(GENlep_id[L1])==13) &&
                 (abs(GENlep_id[L3])==11 || abs(GENlep_id[L3])==13) &&
                 (abs(GENlep_id[L1])!=abs(GENlep_id[L3]) ) ) {GENmass2e2mu = GENmass4l;};
            GENpT4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).Pt();
            GENeta4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).Eta();
            GENrapidity4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).Rapidity();
            GENmassZ1 = (LS3_Z1_1+LS3_Z1_2).M();
            GENmassZ2 = (LS3_Z2_1+LS3_Z2_2).M();
                    
            int tmpIdL1,tmpIdL2,tmpIdL3,tmpIdL4;
            TLorentzVector GENL11P4, GENL12P4, GENL21P4, GENL22P4;
            if(GENlep_id[L1] < 0){ GENL11P4.SetPxPyPzE(LS3_Z1_1.Px(),LS3_Z1_1.Py(),LS3_Z1_1.Pz(),LS3_Z1_1.E()); tmpIdL1 = GENlep_id[L1];}
            else{ GENL11P4.SetPxPyPzE(LS3_Z1_2.Px(),LS3_Z1_2.Py(),LS3_Z1_2.Pz(),LS3_Z1_2.E()); tmpIdL1 = GENlep_id[L2];}
            if(GENlep_id[L2] > 0){ GENL12P4.SetPxPyPzE(LS3_Z1_2.Px(),LS3_Z1_2.Py(),LS3_Z1_2.Pz(),LS3_Z1_2.E()); tmpIdL2 = GENlep_id[L2];}
            else{ GENL12P4.SetPxPyPzE(LS3_Z1_1.Px(),LS3_Z1_1.Py(),LS3_Z1_1.Pz(),LS3_Z1_1.E()); tmpIdL2 = GENlep_id[L1];}
            if(GENlep_id[L3] < 0){ GENL21P4.SetPxPyPzE(LS3_Z2_1.Px(),LS3_Z2_1.Py(),LS3_Z2_1.Pz(),LS3_Z2_1.E()); tmpIdL3 = GENlep_id[L3];}
            else{ GENL21P4.SetPxPyPzE(LS3_Z2_2.Px(),LS3_Z2_2.Py(),LS3_Z2_2.Pz(),LS3_Z2_2.E()); tmpIdL3 = GENlep_id[L4];}
            if(GENlep_id[L4] > 0) { GENL22P4.SetPxPyPzE(LS3_Z2_2.Px(),LS3_Z2_2.Py(),LS3_Z2_2.Pz(),LS3_Z2_2.E()); tmpIdL4 = GENlep_id[L4];}
            else{ GENL22P4.SetPxPyPzE(LS3_Z2_1.Px(),LS3_Z2_1.Py(),LS3_Z2_1.Pz(),LS3_Z2_1.E()); tmpIdL4 = GENlep_id[L3];}
                    
            TUtil::computeAngles(GENcosThetaStar,GENcosTheta1,GENcosTheta2,GENPhi,GENPhi1, \
                                 GENL11P4, tmpIdL1, GENL12P4, tmpIdL2,  \
                                 GENL21P4, tmpIdL3, GENL22P4, tmpIdL4);
                                 

            
        }
        
        bool passedMassOS = true; bool passedElMuDeltaR = true; bool passedDeltaR = true;            
        unsigned int N=GENlep_pt.size();
        for(unsigned int i = 0; i<N; i++) {
            for(unsigned int j = i+1; j<N; j++) {
                        
                // only consider the leptons from Z1 and Z2
                if (!(i==L1 || i==L2 || i==L3 || i==L4)) continue; 
                if (!(j==L1 || j==L2 || j==L3 || j==L4)) continue;
                        
                TLorentzVector li, lj;
                li.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
                lj.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);

                TLorentzVector mll = li+lj;
                        
                if(GENlep_id[i]*GENlep_id[j]<0) {
                    if(mll.M()<=4) { passedMassOS = false; break; }
                }
                        
                if(abs(GENlep_id[i]) != abs(GENlep_id[j])) {
                    double deltaR = li.DeltaR(lj);
                    if(deltaR<=0.02) { passedElMuDeltaR = false; break; }
                }
                double deltaRll = li.DeltaR(lj);
                if(deltaRll<=0.02) { passedDeltaR = false; break; }
            }
        }
        
        if(passedMassOS==false || passedElMuDeltaR==false || passedDeltaR==false) passedFiducialSelection=false;
                
        if (verbose) cout<<"passedFiducialSelection after other cuts? "<<passedFiducialSelection<<endl;

        if (passedFiducialSelection) {

            // DO GEN JETS
            if (verbose) cout<<"begin filling gen jets"<<endl;    
            edm::View<reco::GenJet>::const_iterator genjet;            
            for(genjet = genJets->begin(); genjet != genJets->end(); genjet++) {
                                        
                double pt = genjet->pt();  double eta = genjet->eta();
                if (pt<30.0 || abs(eta)>4.7) continue;

                bool inDR_pt30_eta4p7 = false;
                unsigned int N=GENlep_pt.size();
                for(unsigned int i = 0; i<N; i++) {
                    //if (GENlep_status[i]!=1) continue;
                    if (!(abs(GENlep_id[i])==11 || abs(GENlep_id[i])==13)) continue;
                    TLorentzVector genlep;
                    genlep.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
                    double dR = deltaR(genlep.Eta(), genlep.Phi(), genjet->eta(),genjet->phi());                        
                    if(dR<0.4) {
                        inDR_pt30_eta4p7=true;
                    }                                
                }
                
                if (verbose) cout<<"check overlap of gen jet with gen leptons"<<endl;
                // count number of gen jets which no gen leptons are inside its cone
                if (!inDR_pt30_eta4p7) { 
                    GENnjets_pt30_eta4p7++;
                    GENjet_pt.push_back(genjet->pt());
                    GENjet_eta.push_back(genjet->eta());
                    GENjet_phi.push_back(genjet->phi());
                    GENjet_mass.push_back(genjet->mass());
                    if (pt>GENpt_leadingjet_pt30_eta4p7) {
                        GENpt_leadingjet_pt30_eta4p7=pt;
                        GENabsrapidity_leadingjet_pt30_eta4p7=genjet->rapidity(); //take abs later
                    }
                    if (abs(genjet->eta())<2.5) {
                        GENnjets_pt30_eta2p5++;
                        if (pt>GENpt_leadingjet_pt30_eta2p5) {
                            GENpt_leadingjet_pt30_eta2p5=pt;
                        }
                    }
                }

            }// loop over gen jets

            if (GENnjets_pt30_eta4p7>0) GENabsdeltarapidity_hleadingjet_pt30_eta4p7 = fabs(GENrapidity4l-GENabsrapidity_leadingjet_pt30_eta4p7);
            if (GENnjets_pt30_eta4p7>0) GENabsrapidity_leadingjet_pt30_eta4p7 = fabs(GENabsrapidity_leadingjet_pt30_eta4p7);
            
        } //passedFiducialSelection

    } // 4 fiducial leptons

}

bool UFHZZ4LAna::mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4, bool makeCuts)
{

    double offshell = 999.0; bool findZ1 = false; bool passZ1 = false;

    L1 = 0; L2 = 0;

    unsigned int N = GENlep_pt.size();

    for(unsigned int i=0; i<N; i++){
        for(unsigned int j=i+1; j<N; j++){


            if((GENlep_id[i]+GENlep_id[j])!=0) continue;

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
            lj.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);

            if (verbose) cout<<"gen lep i id: "<<GENlep_id[i]<<" pt: "<<li.Pt()<<" lep j id: "<<GENlep_id[j]<<" pt: "<<lj.Pt()<<endl;

            if (makeCuts) {
                if ( abs(GENlep_id[i]) == 13 && (li.Pt() < 5.0 || abs(li.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[i]) == 11 && (li.Pt() < 7.0 || abs(li.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[i]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;
                
                if ( abs(GENlep_id[j]) == 13 && (lj.Pt() < 5.0 || abs(lj.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[j]) == 11 && (lj.Pt() < 7.0 || abs(lj.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[j]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;                
            }

            TLorentzVector mll = li+lj;
            if (verbose) cout<<"gen mass ij: "<<mll.M()<<endl;

            if(abs(mll.M()-Zmass)<offshell){
                double mZ1 = mll.M();
                if (verbose) cout<<"foundZ1"<<endl;
                L1 = i; L2 = j; findZ1 = true; offshell = abs(mZ1-Zmass);          
            }
        }    
    }

    TLorentzVector l1, l2;
    l1.SetPtEtaPhiM(GENlep_pt[L1],GENlep_eta[L1],GENlep_phi[L1],GENlep_mass[L1]);
    l2.SetPtEtaPhiM(GENlep_pt[L2],GENlep_eta[L2],GENlep_phi[L2],GENlep_mass[L2]);
    TLorentzVector ml1l2 = l1+l2;

    if(ml1l2.M()>40 && ml1l2.M()<120 && findZ1) passZ1 = true;
    if (!makeCuts) passZ1 = true;

    double pTL34 = 0.0; bool findZ2 = false; 
    //bool m4lwindow=false; double window_lo=70.0; double window_hi=140.0;
   
    //cout<<"findZ2"<<endl;
    for(unsigned int i=0; i<N; i++){
        if(i==L1 || i==L2) continue; // can not be the lep from Z1
        for(unsigned int j=i+1; j<N; j++){
            if(j==L1 || j==L2) continue; // can not be the lep from Z1
            if((GENlep_id[i]+GENlep_id[j])!=0) continue;            

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
            lj.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);
            TLorentzVector Z2 = li+lj;

            if (makeCuts) {
                if ( abs(GENlep_id[i]) == 13 && (li.Pt() < 5.0 || abs(li.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[i]) == 11 && (li.Pt() < 7.0 || abs(li.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[i]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;
                
                if ( abs(GENlep_id[j]) == 13 && (lj.Pt() < 5.0 || abs(lj.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[j]) == 11 && (lj.Pt() < 7.0 || abs(lj.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[j]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;
            }

            if ( (li.Pt()+lj.Pt())>=pTL34 ) {
                double mZ2 = Z2.M();
                if (verbose) cout<<"GEN mZ2: "<<mZ2<<endl;
                if( (mZ2>12 && mZ2<120) || (!makeCuts) ) {
                    L3 = i; L4 = j; findZ2 = true; 
                    pTL34 = li.Pt()+lj.Pt();
                    if (verbose) cout<<"is the new GEN cand"<<endl;
                    //if (m4l>window_lo && m4l<window_hi) m4lwindow=true;
                } else {
                    // still assign L3 and L4 to this pair if we don't have a passing Z2 yet
                    if (findZ2 == false) {L3 = i; L4 = j;}
                    //cout<<"is not new GEN cand"<<endl;
                }
            }
            
        } // lj
    } // li

    if(passZ1 && findZ2) return true;
    else return false;
    
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
UFHZZ4LAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(UFHZZ4LAna);
