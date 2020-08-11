#ifndef ZZ4LAnalysisTree_h
#define ZZ4LAnalysisTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <map>
#include <utility>
#include <iterator>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "Math/QuantFuncMathCore.h"

#include "TSystem.h"
#include "TStyle.h"
#include "TPaveText.h"

#include "TPaveLabel.h"
#include "TLegend.h"

#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>
#include <fstream>
#include "TRandom3.h"

#include <algorithm>

#include "TColor.h"

using namespace std;  

////
// output tree
  
bool passedFullSelection, passedZ4lSelection;
bool passedZXCRSelection, passedZ4lZXCRSelection;
bool passedFiducialSelection;
int nZXCRFailedLeptons;
int finalState;

float GENMH;
float GENmass4l,GENmassZ1,GENmassZ2;

bool passTrig;
float pTL1, etaL1, pTErrL1, mL1, phiL1, pTREFITL1, pTErrREFITL1;
float pTL2, etaL2, pTErrL2, mL2, phiL2, pTREFITL2, pTErrREFITL2;
float pTL3, etaL3, pTErrL3, mL3, phiL3, pTREFITL3, pTErrREFITL3;
float pTL4, etaL4, pTErrL4, mL4, phiL4, pTREFITL4, pTErrREFITL4;

int ecalDrivenL1, ecalDrivenL2, ecalDrivenL3, ecalDrivenL4;

float pTGENL1, pTGENL2, pTGENL3, pTGENL4;

float correlation;

int idL1, idL2, idL3, idL4;

float mass4l, mass4l_noFSR, mass4lErr_old, mass4lErr;
float mass4l_up, mass4l_dn;
float mass4lREFIT, mass4lErrREFIT;
float massZ1REFIT, massZ2REFIT;
float mass4mu, mass4e, mass2e2mu;
float pT4l;
float massZ1, massZ1Err, massZ2;
float massZ1_noFSR, massZ2_noFSR;

float njets_pt30_eta4p7;
float pTj1, etaj1;
float pTj2, etaj2;

float D_bkg_kin;
double D_bkg;
double Dgg10_VAMCFM;
double D_g4;
double Djet_VAJHU;
double D_VBF1j_VAJHU;
double D_WHh_VAJHU;
double D_ZHh_VAJHU;
double D_VBF2j;
double D_VBF1j;
double D_WHh;
double D_ZHh;

float D_VBF;
float D_HadWH;
float D_HadZH;
float D_VBF_QG;
float D_VBF1j_QG;
float D_HadWH_QG;
float D_HadZH_QG;

float eventWeight;
int EventCat;

// input tree variables
std::string *triggersPassed;
ULong64_t Run, LumiSect, Event;
bool passedTrig;

float dataMCWeight, genWeight;
float k_qqZZ_qcd_M,k_qqZZ_ewk,k_ggZZ;
float sumW;
int nVtx, nInt; //nPV
int nFSRPhotons;

std::vector<float>* lep_mass;
std::vector<float> *lep_pt; std::vector<float> *lep_eta; std::vector<float> *lep_phi;
std::vector<float>* lepFSR_mass;
std::vector<float> *lepFSR_pt; std::vector<float> *lepFSR_eta; std::vector<float> *lepFSR_phi;

std::vector<float> *lep_genindex;
std::vector<float> *GENlep_pt;

std::vector<int> *lep_tightId;
std::vector<int> *lep_id;
int lep_Hindex[4];
std::vector<float> *lep_RelIso;
std::vector<float> *lep_RelIsoNoFSR;
std::vector<float> *lep_pterr;
std::vector<float> *lep_dataMC;
std::vector<int> *lep_ecalDriven;

std::vector<float> *jet_mass;
std::vector<float> *jet_pt; std::vector<float> *jet_eta; std::vector<float> *jet_phi;
std::vector<int> *jet_iscleanH4l;

std::vector<int> *fsrPhotons_lepindex;
std::vector<float> *fsrPhotons_pt; std::vector<float> *fsrPhotons_eta; std::vector<float> *fsrPhotons_phi;
std::vector<float> *fsrPhotons_pterr;

namespace ZZ4LAnalysisTree {

    void setAddresses(TTree* tree, TString filename){
        
        tree->SetBranchStatus("*",0);

        tree->SetBranchStatus("Run",1);
        tree->SetBranchStatus("LumiSect",1);
        tree->SetBranchStatus("Event",1);
        tree->SetBranchStatus("passedTrig",1);
        tree->SetBranchStatus("triggersPassed",1);
        tree->SetBranchStatus("lep_id",1);
        tree->SetBranchStatus("lep_tightId",1);
        tree->SetBranchStatus("lep_pt",1);
        tree->SetBranchStatus("lep_pterr",1);
        tree->SetBranchStatus("lep_eta",1);
        tree->SetBranchStatus("lep_phi",1);
        tree->SetBranchStatus("lep_mass",1);
        tree->SetBranchStatus("lep_RelIso",1);
        tree->SetBranchStatus("lep_RelIsoNoFSR",1);
        tree->SetBranchStatus("lepFSR_pt",1);
        tree->SetBranchStatus("lepFSR_eta",1);
        tree->SetBranchStatus("lepFSR_phi",1);
        tree->SetBranchStatus("lepFSR_mass",1);
        tree->SetBranchStatus("lep_genindex", 1);
        tree->SetBranchStatus("GENlep_pt", 1);
        
        tree->SetBranchStatus("jet_pt",1);
        tree->SetBranchStatus("jet_eta",1);
        tree->SetBranchStatus("jet_phi",1);
        tree->SetBranchStatus("jet_mass",1);
        //tree->SetBranchStatus("jet_iscleanH4l",1);

        tree->SetBranchStatus("fsrPhotons_pt",1);
        tree->SetBranchStatus("fsrPhotons_pterr",1);
        tree->SetBranchStatus("fsrPhotons_eta",1);
        tree->SetBranchStatus("fsrPhotons_phi",1);
        tree->SetBranchStatus("fsrPhotons_lepindex",1);

        tree->SetBranchAddress("Run",&Run);
        tree->SetBranchAddress("LumiSect",&LumiSect);
        tree->SetBranchAddress("Event",&Event);
        tree->SetBranchAddress("passedTrig",&passedTrig);
        tree->SetBranchAddress("triggersPassed",&triggersPassed);
        
        tree->SetBranchAddress("lep_tightId", &lep_tightId);
        tree->SetBranchAddress("lep_id", &lep_id);
        tree->SetBranchAddress("lep_pt",&lep_pt);
        tree->SetBranchAddress("lep_pterr",&lep_pterr);
        tree->SetBranchAddress("lep_eta",&lep_eta);
        tree->SetBranchAddress("lep_phi",&lep_phi);
        tree->SetBranchAddress("lep_mass",&lep_mass);
        tree->SetBranchAddress("lep_RelIso",&lep_RelIso);
        tree->SetBranchAddress("lep_RelIsoNoFSR",&lep_RelIsoNoFSR);
        tree->SetBranchAddress("lepFSR_pt",&lepFSR_pt);
        tree->SetBranchAddress("lepFSR_eta",&lepFSR_eta);
        tree->SetBranchAddress("lepFSR_phi",&lepFSR_phi);
        tree->SetBranchAddress("lepFSR_mass",&lepFSR_mass);

        tree->SetBranchAddress("jet_pt",&jet_pt);
        tree->SetBranchAddress("jet_eta",&jet_eta);
        tree->SetBranchAddress("jet_phi",&jet_phi);
        tree->SetBranchAddress("jet_mass",&jet_mass);
        //tree->SetBranchAddress("jet_iscleanH4l",&jet_iscleanH4l);

        tree->SetBranchAddress("fsrPhotons_pt",&fsrPhotons_pt);
        tree->SetBranchAddress("fsrPhotons_pterr",&fsrPhotons_pterr);
        tree->SetBranchAddress("fsrPhotons_eta",&fsrPhotons_eta);
        tree->SetBranchAddress("fsrPhotons_phi",&fsrPhotons_phi);
        tree->SetBranchAddress("fsrPhotons_lepindex",&fsrPhotons_lepindex);

        // Event Selection
        tree->SetBranchStatus("lep_Hindex",1);
        tree->SetBranchStatus("passedZ4lSelection",1);
        tree->SetBranchStatus("passedFullSelection",1);
        tree->SetBranchStatus("passedFiducialSelection",1);
        tree->SetBranchStatus("passedZXCRSelection",1);
        tree->SetBranchStatus("nZXCRFailedLeptons",1);
        tree->SetBranchStatus("finalState",1);
        tree->SetBranchStatus("GENMH",1);
        tree->SetBranchStatus("GENmass4l",1);
        tree->SetBranchStatus("GENmassZ1",1);
        tree->SetBranchStatus("GENmassZ2",1);
        tree->SetBranchStatus("dataMCWeight",1);
        tree->SetBranchStatus("k_qqZZ_qcd_M",1);
        tree->SetBranchStatus("k_qqZZ_ewk",1);
        tree->SetBranchStatus("k_ggZZ",1);
        tree->SetBranchStatus("D_bkg_kin",1);
        tree->SetBranchStatus("njets_pt30_eta4p7",1);
        tree->SetBranchStatus("nFSRPhotons",1);
        tree->SetBranchStatus("lep_ecalDriven",1);
        tree->SetBranchStatus("EventCat",1);
        tree->SetBranchStatus("eventWeight",1);
 
	tree->SetBranchAddress("lep_Hindex",lep_Hindex);
        tree->SetBranchAddress("passedZ4lSelection",&passedZ4lSelection);
        tree->SetBranchAddress("passedFullSelection",&passedFullSelection);
        tree->SetBranchAddress("passedFiducialSelection",&passedFiducialSelection);
        tree->SetBranchAddress("passedZXCRSelection",&passedZXCRSelection);
        tree->SetBranchAddress("nZXCRFailedLeptons",&nZXCRFailedLeptons);
        tree->SetBranchAddress("finalState",&finalState);
        tree->SetBranchAddress("GENMH",&GENMH);
        tree->SetBranchAddress("GENmass4l",&GENmass4l);
        tree->SetBranchAddress("GENmassZ1",&GENmassZ1);
        tree->SetBranchAddress("GENmassZ2",&GENmassZ2);
        tree->SetBranchAddress("dataMCWeight",&dataMCWeight);
        tree->SetBranchAddress("k_qqZZ_qcd_M",&k_qqZZ_qcd_M);
        tree->SetBranchAddress("k_qqZZ_ewk",&k_qqZZ_ewk);
        tree->SetBranchAddress("k_ggZZ",&k_ggZZ);
        tree->SetBranchAddress("D_bkg_kin",&D_bkg_kin);
        tree->SetBranchAddress("njets_pt30_eta4p7",&njets_pt30_eta4p7);
        tree->SetBranchAddress("lep_genindex", &lep_genindex);
        tree->SetBranchAddress("GENlep_pt", &GENlep_pt);
        tree->SetBranchAddress("nFSRPhotons", &nFSRPhotons);
        tree->SetBranchAddress("lep_ecalDriven", &lep_ecalDriven);
        tree->SetBranchAddress("EventCat",&EventCat);
        tree->SetBranchAddress("eventWeight",&eventWeight); 
        cout <<"filename "<<filename<<endl;
        cout<<"end loading the tree"<<endl;
        
    }

}

#endif
    
