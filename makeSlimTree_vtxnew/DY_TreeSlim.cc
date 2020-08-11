
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

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
//
#include <vector>
#include <fstream>
//
#include "TRandom3.h"
  
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooCBShape.h"
#include "RooMinuit.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"

#include "RooPlot.h"

// customized PDF
//#include "ZmassConstraintLinkDef.h"
//#include "include/HZZ2L2QRooPdfs.cc"
#include "RooClassFactory.h"
#ifndef __CINT__
#include "MyRooDoubleCBShape.h"
#endif

#include "RooSimultaneous.h"

#include <algorithm>

#include "KalmanMuonCalibrationsProducer/src/RoccoR.cc"

//#include "loader.C"
using namespace std;

/////////////////////

void setAddresses(TTree* tree);
void ReadTree(TTree* tree, TString fs, TTree* & newtree);

bool debug_;

TString filename;

// For HRes Reweighting

using namespace RooFit;
using namespace std;

// tree content

  std::vector<float>* Z_mass=0; 
  std::vector<float>* Z_noFSR_mass=0;
  std::vector<float>* Z_massErr=0;

  std::string *triggersPassed;
  ULong64_t Run, LumiSect, Event;

  std::vector<int> *lep_id=0;
  std::vector<float> *lep_Sip=0;
  std::vector<int> *lep_tightId=0;

  std::vector<float>* lep_mass=0;
  std::vector<float> *lep_pt=0;
  std::vector<float> *lep_eta=0;
  std::vector<float> *lep_phi=0;
  
  std::vector<float>* lepFSR_pt=0;
  std::vector<float>* lepFSR_eta=0;
  std::vector<float>* lepFSR_phi=0;
  std::vector<float>* lepFSR_mass=0;

  std::vector<int> *lep_genindex=0;
  std::vector<float> *lep_RelIso=0;
  std::vector<float> *lep_pterr=0;
  std::vector<float> *lep_pterrold=0;
//  std::vector<double> *lep_ecalTrkEnergyPostCorr=0;
//  std::vector<double> *lep_ecalTrkEnergyPreCorr=0;
  std::vector<float> *lep_dataMC=0;
  std::vector<float> *GENZ_mass=0;
  std::vector<float> *GENlep_pt=0;
  std::vector<float> *GENlep_eta=0;
  std::vector<float> *GENlep_phi=0;
  std::vector<float> *GENlep_mass=0;
  std::vector<float> *eleUncorr_pt=0;
  std::vector<float> *eleUncorr_mass=0;
  std::vector<float> *track_vtx_pt=0;
  std::vector<float> *track_vtx_eta=0;
  std::vector<float> *track_vtx_phi=0;
  std::vector<float> *track_vtx_pterr=0;
  std::vector<float> *track_bs_pt=0;
  std::vector<float> *track_bs_eta=0;
  std::vector<float> *track_bs_phi=0;
  std::vector<float> *track_bs_pterr=0;
  std::vector<int> *passZlepsIndex=0;
//  std::vector<int> *track_nl=0;
//  std::vector<int> *track_type=0;
//  std::vector<int> *track_gen_particle=0;
//  std::vector<float> *track_gen_pt=0;

  std::vector<double> *muon_uncorr_pt=0;
  std::vector<double>* muon_scale_sf=0;

  std::vector<double>* mu_pt_roc=0;
  std::vector<double>* mu_eta_roc=0;
  std::vector<double>* mu_phi_roc=0;
  std::vector<double>* mu_gen_pt_roc=0;
  std::vector<int>* mu_genparticle_roc=0;
  std::vector<int>* mu_nl_roc=0;
  std::vector<int>* mu_q_roc=0;
  std::vector<int>* mu_tracktype_roc=0;

  double massZ, massZErr, massZErrOld, pt1, pt2, eta1, eta2, RelIso1, RelIso2, Sip1, Sip2;
  double FSR_pt1, FSR_pt2, FSR_eta1, FSR_eta2, FSR_phi1, FSR_phi2, FSR_m1, FSR_m2;
  double m1,m2, phi1,phi2;
  double pterr1, pterr2;
  double pterr1old, pterr2old;  
  double genzm, GENmass2l;
  double weight;
  double weight1, weight2;
  
  double massZ_vtx=-999, massZErr_vtx=-999, massZ_bs=-999, massZErr_bs=-999;
  double pt1_vtx=-999, pt2_vtx=-999, eta1_vtx=-999, eta2_vtx=-999, phi1_vtx=-999, phi2_vtx=-999, m1_vtx=-999, m2_vtx=-999, pterr1_vtx=-999, pterr2_vtx=-999;
  double pt1_bs=-999, pt2_bs=-999, eta1_bs=-999, eta2_bs=-999, phi1_bs=-999, phi2_bs=-999, m1_bs=-999, m2_bs=-999, pterr1_bs=-999, pterr2_bs=-999;
  
  double massZ_vtxcorr=-999, massZErr_vtxcorr=-999, massZ_bscorr=-999, massZErr_bscorr=-999;
  double pt1_vtxcorr=-999, pt2_vtxcorr=-999, pt1_bscorr=-999, pt2_bscorr=-999; 
  double pterr1_vtxcorr=-999, pterr2_vtxcorr=-999, pterr1_bscorr=-999, pterr2_bscorr=-999;
  
  double pt1_eleuncorr=-999, pt2_eleuncorr=-999;
  double mass1_eleuncorr=-999, mass2_eleuncorr=-999;

  double pt1_mu_uncorr=-999, pt2_mu_uncorr=-999;
  double pt1_mu_corr=-999, pt2_mu_corr=-999;
  double pterr1_mu_corr=-999, pterr2_mu_corr=-999;
  double mass2mu_corr=-999, mass2muErr_corr=-0.999;
  double mass2mu_uncorr=-999, mass2muErr_uncorr=-999;
  double scale_sf1=-999, scale_sf2=-999;
  double scale_sf1_byhand=-999, scale_sf2_byhand=-999;

//  int lep1_nl=-999, lep2_nl=-999;
//  int lep1_genparticle=-999, lep2_genparticle=-999;
//  double lep1_gen_pt=-999, lep2_gen_pt=-999;
//  int lep1_type=-999, lep2_type=-999;
  
  double mu1_pt_roc=-999, mu2_pt_roc=-999;
  double mu1_eta_roc=-999, mu2_eta_roc=-999;
  double mu1_phi_roc=-999, mu2_phi_roc=-999;
  double mu1_gen_pt_roc=-999, mu2_gen_pt_roc=-999;
  int mu1_genparticle_roc=-999, mu2_genparticle_roc=-999;
  int mu1_nl_roc=-999, mu2_nl_roc=-999;
  int mu1_q_roc=-999, mu2_q_roc=-999;
  int mu1_tracktype_roc=-999, mu2_tracktype_roc=-999;

  int q1, q2;
  
  double genLep_pt1=-999, genLep_pt2=-999;
  double genLep_eta1=-999, genLep_eta2=-999;
  double genLep_phi1=-999, genLep_phi2=-999;

  int lep1_ecalDriven = -1, lep2_ecalDriven = -1;

  int nFSRPhotons;
  double Met;
  float_t met;
  bool passedTrig;
  bool passedFullSelection; 

  std::vector<int> *lep_ecalDriven;
  
  bool isMC;
  int year;

int main(int argc, char *argv[])
{    
     
  gStyle->SetOptStat(0000);
  //gStyle->SetOptTitle(0);
     
  using namespace std;
     
  if(argc != 7)  {
      cout<<argv[0]<<" filename "<<argv[1]<<" fs "<<argv[2]<< " indir " << argv[3] << " outdir" <<argv[4]<<endl;
      return -1;
    }

  gStyle->SetOptStat(0000);
  //gStyle->SetOptTitle(0);
  //gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libHiggsAnalysisCombinedLimit.so");
  //gSystem->Load("$CMSSW_BASE/lib/slc6_amd64_gcc493/pluginUFHZZAnalysisRun2UFHZZ4LAna.so");

  // control ALL txt format
  debug_ =   true;

  /////////////////////

  TString filename = argv[1];
  TString fs = argv[2];
  TString indir = argv[3];
  TString outdir = argv[4];
  TString isMC_s = argv[5];
  TString year_s = argv[6]; 
  
  if(isMC_s=="MC"){
	  isMC = 1;
  }else{
	  isMC = 0;
  }
  
  if(year_s=="2016"){
	  year = 2016;
  }else if(year_s=="2017"){
	  year = 2017;
  }else{
	  year = 2018;
  }
  
  if(fs!="2e" && fs!="2mu")
  cout<<"fs has to be 2e, or 2mu"<<endl;

  /////////////////////

  cout<<"fs is "<<fs<<endl;

  cout<<"read file"<<endl;

  TFile* infile = TFile::Open(indir+filename+".root");
  TTree* tree; 

  if(infile) tree = (TTree*) infile->Get("Ana/passedEvents");
  if(!tree) { cout<<"ERROR could not find the tree for "<<filename<<endl; return -1;}

  // read tree     
  TString name = filename+"_m"+fs;
  TFile* tmpFile =  new TFile(outdir + name+".root","RECREATE");

  TTree* newtree = new TTree("passedEvents","passedEvents");

  cout<<"start setting tree "<<endl;
  newtree->Branch("passedFullSelection",&passedFullSelection,"passedFullSelection/O");
  newtree->Branch("massZ",&massZ,"massZ/D");
  newtree->Branch("massZErr",&massZErr,"massZErr/D");
  newtree->Branch("massZErrOld",&massZErrOld,"massZErrOld/D");

  newtree->Branch("pt1",&pt1,"pt1/D");
  newtree->Branch("pt2",&pt2,"pt2/D");
  newtree->Branch("eta1",&eta1,"eta1/D");
  newtree->Branch("eta2",&eta2,"eta2/D");
  newtree->Branch("phi1",&phi1,"phi1/D");
  newtree->Branch("phi2",&phi2,"phi2/D");
  newtree->Branch("m1",&m1,"m1/D");
  newtree->Branch("m2",&m2,"m2/D");

  newtree->Branch("FSR_pt1",&FSR_pt1,"FSR_pt1/D");
  newtree->Branch("FSR_pt2",&FSR_pt2,"FSR_pt2/D");
  newtree->Branch("FSR_eta1",&FSR_eta1,"FSR_eta1/D");
  newtree->Branch("FSR_eta2",&FSR_eta2,"FSR_eta2/D");
  newtree->Branch("FSR_phi1",&FSR_phi1,"FSR_phi1/D");
  newtree->Branch("FSR_phi2",&FSR_phi2,"FSR_phi2/D");
  newtree->Branch("FSR_m1",&FSR_m1,"FSR_m1/D");
  newtree->Branch("FSR_m2",&FSR_m2,"FSR_m2/D");

  newtree->Branch("RelIso1",&RelIso1,"RelIso1/D");
  newtree->Branch("RelIso2",&RelIso2,"RelIso2/D");
  newtree->Branch("Sip1",&Sip1,"Sip1/D");
  newtree->Branch("Sip2",&Sip2,"Sip2/D");
  newtree->Branch("pterr1",&pterr1,"pterr1/D");
  newtree->Branch("pterr2",&pterr2,"pterr2/D");
  newtree->Branch("pterr1old",&pterr1old,"pterr1old/D");
  newtree->Branch("pterr2old",&pterr2old,"pterr2old/D");
  newtree->Branch("Met", &Met, "Met/D");
  newtree->Branch("weight",&weight,"weight/D");
  newtree->Branch("weight1",&weight1,"weight1/D");
  newtree->Branch("weight2",&weight2,"weight2/D");
  newtree->Branch("genzm",&genzm,"genzm/D");
  newtree->Branch("GENmass2l",&GENmass2l,"GENmass2l/D");
  newtree->Branch("genLep_pt1", &genLep_pt1, "genLep_pt1/D");
  newtree->Branch("genLep_pt2", &genLep_pt2, "genLep_pt2/D");
  newtree->Branch("genLep_eta1", &genLep_eta1, "genLep_eta1/D");
  newtree->Branch("genLep_eta2", &genLep_eta2, "genLep_eta2/D");
  newtree->Branch("genLep_phi1", &genLep_phi1, "genLep_phi1/D");
  newtree->Branch("genLep_phi2", &genLep_phi2, "genLep_phi2/D");

  newtree->Branch("nFSRPhotons", &nFSRPhotons, "nFSRPhotons/I");
  newtree->Branch("lep1_ecalDriven", &lep1_ecalDriven, "lep1_ecalDriven/I");
  newtree->Branch("lep2_ecalDriven", &lep2_ecalDriven, "lep2_ecalDriven/I");

  newtree->Branch("pt1_vtx", &pt1_vtx, "pt1_vtx/D");
  newtree->Branch("pt2_vtx", &pt2_vtx, "pt2_vtx/D");
  newtree->Branch("pterr1_vtx", &pterr1_vtx, "pterr1_vtx/D");
  newtree->Branch("pterr2_vtx", &pterr2_vtx, "pterr2_vtx/D");
  newtree->Branch("eta1_vtx", &eta1_vtx, "eta1_vtx/D");
  newtree->Branch("eta2_vtx", &eta2_vtx, "eta2_vtx/D");
  newtree->Branch("phi1_vtx", &phi1_vtx, "phi1_vtx/D");
  newtree->Branch("phi2_vtx", &phi2_vtx, "phi2_vtx/D");
  newtree->Branch("pt1_vtxcorr", &pt1_vtxcorr, "pt1_vtxcorr/D");
  newtree->Branch("pt2_vtxcorr", &pt2_vtxcorr, "pt2_vtxcorr/D");
  newtree->Branch("pterr1_vtxcorr", &pterr1_vtxcorr, "pterr1_vtxcorr/D");
  newtree->Branch("pterr2_vtxcorr", &pterr2_vtxcorr, "pterr2_vtxcorr/D");
  
  newtree->Branch("pt1_bs", &pt1_bs, "pt1_bs/D");
  newtree->Branch("pt2_bs", &pt2_bs, "pt2_bs/D");
  newtree->Branch("pterr1_bs", &pterr1_bs, "pterr1_bs/D");
  newtree->Branch("pterr2_bs", &pterr2_bs, "pterr2_bs/D");
  newtree->Branch("eta1_bs", &eta1_bs, "eta1_bs/D");
  newtree->Branch("eta2_bs", &eta2_bs, "eta2_bs/D");
  newtree->Branch("phi1_bs", &phi1_bs, "phi1_bs/D");
  newtree->Branch("phi2_bs", &phi2_bs, "phi2_bs/D");
  newtree->Branch("pt1_bscorr", &pt1_bscorr, "pt1_bscorr/D");
  newtree->Branch("pt2_bscorr", &pt2_bscorr, "pt2_bscorr/D");
  newtree->Branch("pterr1_bscorr", &pterr1_bscorr, "pterr1_bscorr/D");
  newtree->Branch("pterr2_bscorr", &pterr2_bscorr, "pterr2_bscorr/D");
  
  newtree->Branch("massZ_vtx", &massZ_vtx, "massZ_vtx/D");
  newtree->Branch("massZ_bs", &massZ_bs, "massZ_bs/D");
  newtree->Branch("massZ_vtxcorr", &massZ_vtxcorr, "massZ_vtxcorr/D");
  newtree->Branch("massZ_bscorr", &massZ_bscorr, "massZ_bscorr/D");
  
  newtree->Branch("massZErr_vtx", &massZErr_vtx, "massZErr_vtx/D");
  newtree->Branch("massZErr_bs", &massZErr_bs, "massZErr_bs/D");
  newtree->Branch("massZErr_vtxcorr", &massZErr_vtxcorr, "massZErr_vtxcorr/D");
  newtree->Branch("massZErr_bscorr", &massZErr_bscorr, "massZErr_bscorr/D");
  
  newtree->Branch("pt1_eleuncorr", &pt1_eleuncorr, "pt1_eleuncorr/D");
  newtree->Branch("pt2_eleuncorr", &pt2_eleuncorr, "pt2_eleuncorr/D");
  newtree->Branch("mass1_eleuncorr", &mass1_eleuncorr, "mass1_eleuncorr/D");
  newtree->Branch("mass2_eleuncorr", &mass2_eleuncorr, "mass2_eleuncorr/D");

  newtree->Branch("pt1_mu_uncorr",&pt1_mu_uncorr, "pt1_mu_uncorr/D");
  newtree->Branch("pt2_mu_uncorr",&pt2_mu_uncorr, "pt2_mu_uncorr/D");
  newtree->Branch("pt1_mu_corr", &pt1_mu_corr, "pt1_mu_corr/D");
  newtree->Branch("pt2_mu_corr", &pt2_mu_corr, "pt2_mu_corr/D");
  newtree->Branch("pterr1_mu_corr",&pterr1_mu_corr, "pterr1_mu_corr/D");
  newtree->Branch("pterr2_mu_corr",&pterr2_mu_corr, "pterr2_mu_corr/D");
  newtree->Branch("mass2mu_corr",&mass2mu_corr, "mass2mu_corr/D");
  newtree->Branch("mass2muErr_corr",&mass2muErr_corr, "mass2muErr_corr/D");
  newtree->Branch("mass2mu_uncorr",&mass2mu_uncorr,"mass2mu_uncorr/D");
  newtree->Branch("mass2muErr_uncorr",&mass2muErr_uncorr,"mass2muErr_uncorr/D");
  newtree->Branch("scale_sf1",&scale_sf1,"scale_sf1/D");
  newtree->Branch("scale_sf2",&scale_sf2,"scale_sf2/D");
  newtree->Branch("scale_sf1_byhand",&scale_sf1_byhand,"scale_sf1_byhand/D");
  newtree->Branch("scale_sf2_byhand",&scale_sf2_byhand,"scale_sf2_byhand/D");
  
//  newtree->Branch("lep1_nl",&lep1_nl, "lep1_nl/I");
//  newtree->Branch("lep2_nl",&lep2_nl, "lep2_nl/I");
//  newtree->Branch("lep1_type",&lep1_type, "lep1_type/I");
//  newtree->Branch("lep2_type",&lep2_type, "lep2_type/I");
//  newtree->Branch("lep1_genparticle",&lep1_genparticle,"lep1_genparticle/I");
//  newtree->Branch("lep2_genparticle",&lep2_genparticle,"lep2_genparticle/I");
//  newtree->Branch("lep1_gen_pt",&lep1_gen_pt,"lep1_gen_pt/D");
//  newtree->Branch("lep2_gen_pt",&lep2_gen_pt,"lep2_gen_pt/D");

  newtree->Branch("mu1_pt_roc",&mu1_pt_roc,"mu1_pt_roc/D");
  newtree->Branch("mu2_pt_roc",&mu2_pt_roc,"mu2_pt_roc/D");
  newtree->Branch("mu1_eta_roc",&mu1_eta_roc,"mu1_eta_roc/D");
  newtree->Branch("mu2_eta_roc",&mu2_eta_roc,"mu2_eta_roc/D");
  newtree->Branch("mu1_phi_roc",&mu1_phi_roc,"mu1_phi_roc/D");
  newtree->Branch("mu2_phi_roc",&mu2_phi_roc,"mu2_phi_roc/D");
  newtree->Branch("mu1_gen_pt_roc",&mu1_gen_pt_roc,"mu1_gen_pt_roc/D");
  newtree->Branch("mu2_gen_pt_roc",&mu2_gen_pt_roc,"mu2_gen_pt_roc/D");
  newtree->Branch("mu1_genparticle_roc",&mu1_genparticle_roc,"mu1_genparticle_roc/I");
  newtree->Branch("mu2_genparticle_roc",&mu2_genparticle_roc,"mu2_genparticle_roc/I");
  newtree->Branch("mu1_nl_roc",&mu1_nl_roc,"mu1_nl_roc/I");
  newtree->Branch("mu2_nl_roc",&mu2_nl_roc,"mu2_nl_roc/I");
  newtree->Branch("mu1_q_roc",&mu1_q_roc,"mu1_q_roc/I");
  newtree->Branch("mu2_q_roc",&mu2_q_roc,"mu2_q_roc/I");
  newtree->Branch("mu1_tracktyep_roc",&mu1_tracktype_roc,"mu1_tracktype_roc/I");
  newtree->Branch("mu2_tracktype_roc",&mu2_tracktype_roc,"mu2_tracktype_roc/I");

//  newtree->Branch("q1",&q1,"q1/I");
//  newtree->Branch("q2",&q2,"q2/I");

  cout<<"start reading tree "<<endl;

  ReadTree(tree, fs, newtree);

  cout<<"end reading tree"<<endl;

  tmpFile->cd();

  newtree->Write("passedEvents",TObject::kOverwrite);

  tmpFile->Write();
  tmpFile->Close(); 

  //delete infile; delete tmpFile;


}



void ReadTree(TTree* tree, TString fs, TTree* & newtree){
    
	std::string name;
	if(year==2016)name="KalmanMuonCalibrationsProducer/data/roccor.Run2.v3/RoccoR2016.txt";
	if(year==2017)name="KalmanMuonCalibrationsProducer/data/roccor.Run2.v3/RoccoR2017.txt";
	if(year==2018)name="KalmanMuonCalibrationsProducer/data/roccor.Run2.v3/RoccoR2018.txt";
	RoccoR* calibrator = new RoccoR(name);
	setAddresses(tree);
        for(int mcfmEvt_HZZ=0; mcfmEvt_HZZ < tree->GetEntries(); mcfmEvt_HZZ++) { //event loop
    	    tree->GetEntry(mcfmEvt_HZZ);
//          if(!passedTrig) continue;
	    if((*lep_id).size() < 2)continue;
	    vector<int> passLepIndex; passLepIndex.clear();
	    for(unsigned int il = 0; il < (*lep_id).size(); il++){
			if(!(*lep_tightId)[il])continue;
			if((*lep_RelIso)[il]>0.35)continue;
			passLepIndex.push_back(il);
	    }
            if(passLepIndex.size()!=2)continue;
	    unsigned int L1 = passLepIndex[0]; unsigned int L2 = passLepIndex[1];
	    int idL1 = (*lep_id)[L1]; int idL2 = (*lep_id)[L2];
	    if((idL1+idL2)!=0)continue;
	    if(fs=="2e"&&abs(idL1)!=11)continue;
	    if(fs=="2mu"&&abs(idL1)!=13)continue;

	    
	    mu1_pt_roc=(*mu_pt_roc)[L1]; mu2_pt_roc=(*mu_pt_roc)[L2];
	    mu1_eta_roc=(*mu_eta_roc)[L1]; mu2_eta_roc=(*mu_eta_roc)[L2];
	    mu1_phi_roc=(*mu_phi_roc)[L1]; mu2_phi_roc=(*mu_phi_roc)[L2];
	    mu1_gen_pt_roc=(*mu_gen_pt_roc)[L1]; mu2_gen_pt_roc=(*mu_gen_pt_roc)[L2];
	    mu1_genparticle_roc=(*mu_genparticle_roc)[L1]; mu2_genparticle_roc=(*mu_genparticle_roc)[L2];
	    mu1_nl_roc=(*mu_nl_roc)[L1]; mu2_nl_roc=(*mu_nl_roc)[L2];
	    mu1_q_roc=(*mu_q_roc)[L1]; mu2_q_roc=(*mu_q_roc)[L2];
	    mu1_tracktype_roc=(*mu_tracktype_roc)[L1]; mu2_tracktype_roc=(*mu_tracktype_roc)[L2];
	    
	    q1 = mu1_q_roc;
	    q2 = mu2_q_roc;

	    //for no vertex constraint			
	    weight1 = (*lep_dataMC)[L1];
	    weight2 = (*lep_dataMC)[L2];
            weight = (*lep_dataMC)[L1]*(*lep_dataMC)[L2];
  	 
            pt1 = (*lep_pt)[L1]; pt2 = (*lep_pt)[L2];
	    eta1 = (*lep_eta)[L1]; eta2 = (*lep_eta)[L2];
	    phi1 = (*lep_phi)[L1]; phi2 = (*lep_phi)[L2];
	    m1 = (*lep_mass)[L1]; m2 = (*lep_mass)[L2];

	    FSR_eta1 = (*lepFSR_eta)[L1]; FSR_eta2 = (*lepFSR_eta)[L2];
            FSR_phi1 = double((*lepFSR_phi)[L1]); FSR_m1 = double((*lepFSR_mass)[L1]);//old variable is *lep_phi, change it to lep_mass
            FSR_phi2 = double((*lepFSR_phi)[L2]); FSR_m2 = double((*lepFSR_mass)[L2]);
            FSR_pt1 = (*lepFSR_pt)[L1]; FSR_pt2 = (*lepFSR_pt)[L2];
	    
	    RelIso1 = (*lep_RelIso)[L1]; RelIso2 = (*lep_RelIso)[L2];
	    Sip1 = (*lep_Sip)[L1]; Sip2 = (*lep_Sip)[L2];
            scale_sf1 = (*muon_scale_sf)[L1]; scale_sf2 = (*muon_scale_sf)[L2];

	    TLorentzVector lep1(0,0,0,0);
	    TLorentzVector lep2(0,0,0,0);
	    lep1.SetPtEtaPhiM(double((*lepFSR_pt)[L1]),double((*lepFSR_eta)[L1]),double((*lepFSR_phi)[L1]),double((*lepFSR_mass)[L1]));
            lep2.SetPtEtaPhiM(double((*lepFSR_pt)[L2]),double((*lepFSR_eta)[L2]),double((*lepFSR_phi)[L2]),double((*lepFSR_mass)[L2]));
            massZ = (lep1+lep2).M();
            
	    pterr1 = double((*lep_pterr)[L1]); pterr2 = double((*lep_pterr)[L2]);
            pterr1old = double((*lep_pterrold)[L1]); pterr2old = double((*lep_pterrold)[L2]);
            TLorentzVector lep1p, lep2p;
            lep1p.SetPtEtaPhiM(double((*lepFSR_pt)[L1]+pterr1),double((*lepFSR_eta)[L1]),double((*lepFSR_phi)[L1]),double((*lepFSR_mass)[L1]));
            lep2p.SetPtEtaPhiM(double((*lepFSR_pt)[L2]+pterr2),double((*lepFSR_eta)[L2]),double((*lepFSR_phi)[L2]),double((*lepFSR_mass)[L2]));
            double dm1 = (lep1p+lep2).M()-(lep1+lep2).M();
            double dm2 = (lep1+lep2p).M()-(lep1+lep2).M(); 
            massZErr = TMath::Sqrt(dm1*dm1+dm2*dm2);
			
            lep1p.SetPtEtaPhiM(double((*lepFSR_pt)[L1]+pterr1old),double((*lepFSR_eta)[L1]),double((*lepFSR_phi)[L1]),double((*lepFSR_mass)[L1]));
            lep2p.SetPtEtaPhiM(double((*lepFSR_pt)[L2]+pterr2old),double((*lepFSR_eta)[L2]),double((*lepFSR_phi)[L2]),double((*lepFSR_mass)[L2]));
            dm1 = (lep1p+lep2).M()-(lep1+lep2).M();
            dm2 = (lep1+lep2p).M()-(lep1+lep2).M();
            massZErrOld = TMath::Sqrt(dm1*dm1+dm2*dm2);
			
	    //gen level information
            Met = met; 
	    genzm=0; GENmass2l=0;
            if(GENZ_mass->size()>0) genzm = (*GENZ_mass)[0];
            TLorentzVector GENlep1p, GENlep2p;
            if((*lep_genindex)[L1] >= 0 && (*lep_genindex)[L2] >= 0) {
            		genLep_pt1=(*GENlep_pt)[(*lep_genindex)[L1]]; genLep_pt2=(*GENlep_pt)[(*lep_genindex)[L2]];
            		genLep_eta1=(*GENlep_eta)[(*lep_genindex)[L1]]; genLep_eta2=(*GENlep_eta)[(*lep_genindex)[L2]];
            		genLep_phi1=(*GENlep_phi)[(*lep_genindex)[L1]]; genLep_phi2=(*GENlep_phi)[(*lep_genindex)[L2]];		
            		int genindex1 = (*lep_genindex)[L1];
            		int genindex2 = (*lep_genindex)[L2];
              		GENlep1p.SetPtEtaPhiM(double((*GENlep_pt)[genindex1]),double((*GENlep_eta)[genindex1]),double((*GENlep_phi)[genindex1]),double((*GENlep_mass)[genindex1]));
              		GENlep2p.SetPtEtaPhiM(double((*GENlep_pt)[genindex2]),double((*GENlep_eta)[genindex2]),double((*GENlep_phi)[genindex2]),double((*GENlep_mass)[genindex2]));
              		GENmass2l = (GENlep1p+GENlep2p).M();
            }
            if (lep_ecalDriven->size() > 0) {
               lep1_ecalDriven = (*lep_ecalDriven)[L1];
               lep2_ecalDriven = (*lep_ecalDriven)[L2];
            }
			
            //ele scale unc
	    if(fs=="2e"&&abs(idL1)==11){
			pt1_eleuncorr = (*eleUncorr_pt)[L1]; pt2_eleuncorr = (*eleUncorr_pt)[L2];
			mass1_eleuncorr = (*eleUncorr_mass)[L1]; mass2_eleuncorr = (*eleUncorr_mass)[L2];
	    }			
	    //roche corr for muons			
	    if(fs=="2mu"&&(*passZlepsIndex).size()==2&&(*track_vtx_pt).size()==2&&(*track_bs_pt).size()==2&&((*lep_id)[(*passZlepsIndex)[0]]+(*lep_id)[(*passZlepsIndex)[1]])==0 && (*muon_uncorr_pt).size()>=2 ){
				
				int vtx_idL1 = (*lep_id)[(*passZlepsIndex)[0]];
				int vtx_idL2 = (*lep_id)[(*passZlepsIndex)[1]];
				unsigned int vtx_L1 = (*passZlepsIndex)[0]; unsigned int vtx_L2 = (*passZlepsIndex)[1];//used to extract track_type, track_nl, track_gen_particle, and track_gen_pt, and muon uncorr pt
				
				pt1_vtx = (*track_vtx_pt)[0]; pt2_vtx = (*track_vtx_pt)[1];
				pterr1_vtx = (*track_vtx_pterr)[0]; pterr2_vtx = (*track_vtx_pterr)[1];
				eta1_vtx = (*track_vtx_eta)[0]; eta2_vtx = (*track_vtx_eta)[1];
				phi1_vtx = (*track_vtx_phi)[0]; phi2_vtx = (*track_vtx_phi)[1];
				
				pt1_bs = (*track_bs_pt)[0]; pt2_bs = (*track_bs_pt)[1];
				pterr1_bs = (*track_bs_pterr)[0]; pterr2_bs = (*track_bs_pterr)[1];
				eta1_bs = (*track_bs_eta)[0]; eta2_bs = (*track_bs_eta)[1];
				phi1_bs = (*track_bs_phi)[0]; phi2_bs = (*track_bs_phi)[1];
				
				pt1_mu_uncorr = (*muon_uncorr_pt)[vtx_L1];
				pt2_mu_uncorr = (*muon_uncorr_pt)[vtx_L2];

	//			if(vtx_idL1>0)q1 = 1;if(vtx_idL1<0)q1 = -1;
	//			if(vtx_idL2>0)q2 = 1;if(vtx_idL2<0)q1 = -1;
				
				for(int c = 0; c < 3; c++){
					std::vector<float> tmp_pt; std::vector<float> tmp_eta; std::vector<float> tmp_phi; std::vector<float> tmp_pterr;
					std::vector<int> tmp_nl; std::vector<int> tmp_genParticle; std::vector<float> tmp_type; std::vector<int> tmp_q; std::vector<float> tmp_gen_pt;
					tmp_pt.clear(); tmp_eta.clear(); tmp_phi.clear(); tmp_pterr.clear(); 
					tmp_nl.clear(); tmp_genParticle.clear(); tmp_type.clear(); tmp_q.clear(); tmp_gen_pt.clear();
				
					tmp_nl.push_back((*mu_nl_roc)[vtx_L1]); tmp_nl.push_back((*mu_nl_roc)[vtx_L2]);
					//lep1_nl = tmp_nl[0]; lep2_nl = tmp_nl[1];
					tmp_genParticle.push_back((*mu_genparticle_roc)[vtx_L1]); tmp_genParticle.push_back((*mu_genparticle_roc)[vtx_L2]);
					//lep1_genparticle = tmp_genParticle[0]; lep2_genparticle = tmp_genParticle[1];
					tmp_type.push_back((*mu_tracktype_roc)[vtx_L1]); tmp_type.push_back((*mu_tracktype_roc)[vtx_L2]);
					//lep1_type = tmp_type[0]; lep2_type = tmp_type[1];
					tmp_gen_pt.push_back((*mu_gen_pt_roc)[vtx_L1]); tmp_gen_pt.push_back((*mu_gen_pt_roc)[vtx_L2]);
					//lep1_gen_pt = tmp_gen_pt[0]; lep2_gen_pt = tmp_gen_pt[1];
					tmp_q.push_back(q1); tmp_q.push_back(q2);
				
					if(c==0){
						tmp_pt.push_back(pt1_vtx); tmp_pt.push_back(pt2_vtx);
						tmp_eta.push_back(eta1_vtx); tmp_eta.push_back(eta2_vtx);
						tmp_phi.push_back(phi1_vtx); tmp_phi.push_back(phi2_vtx);
						tmp_pterr.push_back(pterr1_vtx); tmp_pterr.push_back(pterr2_vtx);	
					}
					if(c==1){
						tmp_pt.push_back(pt1_bs); tmp_pt.push_back(pt2_bs);
						tmp_eta.push_back(eta1_bs); tmp_eta.push_back(eta2_bs);
						tmp_phi.push_back(phi1_bs); tmp_phi.push_back(phi2_bs);
						tmp_pterr.push_back(pterr1_bs); tmp_pterr.push_back(pterr2_bs);
					}
					if(c==2){
						tmp_pt.push_back(pt1_mu_uncorr); tmp_pt.push_back(pt2_mu_uncorr);
						tmp_eta.push_back(eta1); tmp_eta.push_back(eta2);
						tmp_phi.push_back(phi1); tmp_phi.push_back(phi2);
						tmp_pterr.push_back(pterr1old); tmp_pterr.push_back(pterr2old);
					}
					for(int i = 0; i < 2; i++){
						double oldpt = tmp_pt[i];
						double newpt = oldpt;
						double oldpterr = tmp_pterr[i];
						double newpterr = oldpterr;
						double scale_factor = 1.0;
						TRandom3 rand;
						rand.SetSeed(abs(static_cast<int>(sin(phi1_vtx))));
						double u1=rand.Uniform(1.0);
						if(tmp_type[i]==1 && oldpt < 200){
							if(isMC && tmp_nl[i] > 5){
								if(tmp_genParticle[i] != 0 && tmp_gen_pt[i] != 0){
									scale_factor = calibrator->kSpreadMC(tmp_q[i], oldpt, tmp_eta[i], tmp_phi[i], tmp_gen_pt[i]);
								}
								else{
									scale_factor = calibrator->kSmearMC(tmp_q[i], oldpt, tmp_eta[i], tmp_phi[i], tmp_nl[i], u1);
								}
								newpt = scale_factor*oldpt;
								newpterr = scale_factor*oldpterr;
							}
							else if(!isMC && tmp_nl[i] > 5){
								if(oldpt > 2.0 && fabs(tmp_eta[i]) < 2.4){
									scale_factor = calibrator->kScaleDT(tmp_q[i], oldpt, tmp_eta[i], tmp_phi[i]);
								}
								else{
									scale_factor = 1.0;
								}
								newpt = scale_factor*oldpt;
								newpterr = scale_factor*oldpterr;
							}
						}
						if(c==0 && i == 0){
							pt1_vtxcorr = newpt; pterr1_vtxcorr = newpterr;
						}
						else if(c==0 && i==1){
							pt2_vtxcorr = newpt; pterr1_vtxcorr = newpterr;
						}
						else if(c==1 && i==0){
							pt1_bscorr = newpt; pterr1_bscorr = newpterr;
						}
						else if(c==1 && i==1){
							pt2_bscorr = newpt; pterr2_bscorr = newpterr;
						}
						else if(c==2 && i==0){
							pt1_mu_corr = newpt; pterr1_mu_corr = newpterr; scale_sf1_byhand = scale_factor;
						}
						else{
							pt2_mu_corr = newpt; pterr2_mu_corr = newpterr; scale_sf2_byhand = scale_factor;
						}

					}//track_vtx_pt.size()==2
				}//vtx and bs 
				
				TLorentzVector mu1, mu2, mu1p, mu2p;
				
				mu1.SetPtEtaPhiM(pt1_vtx, eta1_vtx, phi1_vtx, 0.1057);
				mu2.SetPtEtaPhiM(pt2_vtx, eta2_vtx, phi2_vtx, 0.1057);
				massZ_vtx = (mu1+mu2).M();
				mu1p.SetPtEtaPhiM(pt1_vtx+pterr1_vtx, eta1_vtx, phi1_vtx, 0.1057);
				mu2p.SetPtEtaPhiM(pt2_vtx+pterr2_vtx, eta2_vtx, phi2_vtx, 0.1057);
				double dm1 = (mu1p+mu2).M()-massZ_vtx;
				double dm2 = (mu1+mu2p).M()-massZ_vtx;
				massZErr_vtx = sqrt(dm1*dm1+dm2*dm2);
				
				mu1.SetPtEtaPhiM(pt1_bs, eta1_bs, phi1_bs, 0.1057);
				mu2.SetPtEtaPhiM(pt2_bs, eta2_bs, phi2_bs, 0.1057);
				massZ_bs = (mu1+mu2).M();
				mu1p.SetPtEtaPhiM(pt1_bs+pterr1_bs, eta1_bs, phi1_bs, 0.1057);
				mu2p.SetPtEtaPhiM(pt2_bs+pterr2_bs, eta2_bs, phi2_bs, 0.1057);
				dm1 = (mu1p+mu2).M()-massZ_bs;
				dm2 = (mu1+mu2p).M()-massZ_bs;
				massZErr_bs = sqrt(dm1*dm1+dm2*dm2);
				
				mu1.SetPtEtaPhiM(pt1_vtxcorr, eta1_vtx, phi1_vtx, 0.1057);
				mu2.SetPtEtaPhiM(pt2_vtxcorr, eta2_vtx, phi2_vtx, 0.1057);
				massZ_vtxcorr = (mu1+mu2).M();
				mu1p.SetPtEtaPhiM(pt1_vtxcorr+pterr1_vtxcorr, eta1_vtx, phi1_vtx, 0.1057);
				mu2p.SetPtEtaPhiM(pt2_vtxcorr+pterr2_vtxcorr, eta2_vtx, phi2_vtx, 0.1057);
				dm1 = (mu1p+mu2).M()-massZ_vtxcorr;
				dm2 = (mu1+mu2p).M()-massZ_vtxcorr;
				massZErr_vtxcorr = sqrt(dm1*dm1+dm2*dm2);
				
				mu1.SetPtEtaPhiM(pt1_bscorr, eta1_bs, phi1_bs, 0.1057);
				mu2.SetPtEtaPhiM(pt2_bscorr, eta2_bs, phi2_bs, 0.1057);
				massZ_bscorr = (mu1+mu2).M();
				mu1p.SetPtEtaPhiM(pt1_bscorr+pterr1_bscorr, eta1_bs, phi1_bs, 0.1057);
				mu2p.SetPtEtaPhiM(pt2_bscorr+pterr2_bscorr, eta2_bs, phi2_bs, 0.1057);
				dm1 = (mu1p+mu2).M()-massZ_bscorr;
				dm2 = (mu1+mu2p).M()-massZ_bscorr;
				massZErr_bscorr = sqrt(dm1*dm1+dm2*dm2);

				mu1.SetPtEtaPhiM(pt1_mu_uncorr,eta1,phi1,0.1057);
				mu2.SetPtEtaPhiM(pt2_mu_uncorr,eta2,phi2,0.1057);
				mass2mu_uncorr = (mu1+mu2).M();
				mu1p.SetPtEtaPhiM(pt1_mu_uncorr+pterr1old,eta1,phi1,0.1057);
				mu2p.SetPtEtaPhiM(pt2_mu_uncorr+pterr2old,eta2,phi2,0.1057);
				dm1 = (mu1p+mu2).M();
				dm2 = (mu1+mu2p).M();
				mass2muErr_uncorr = sqrt(dm1*dm1+dm2*dm2);

				mu1.SetPtEtaPhiM(pt1_mu_corr,eta1,phi1,0.1057);
				mu2.SetPtEtaPhiM(pt2_mu_corr,eta2,phi2,0.1057);
				mass2mu_corr = (mu1+mu2).M();
				mu1p.SetPtEtaPhiM(pt1_mu_corr+pterr1_mu_corr,eta1,phi1,0.1057);
				mu2p.SetPtEtaPhiM(pt2_mu_corr+pterr2_mu_corr,eta2,phi2,0.1057);
				dm1 = (mu1p+mu2).M()-mass2mu_corr;
				dm2 = (mu1+mu2p).M()-mass2mu_corr;
				mass2muErr_corr = sqrt(dm1*dm1+dm2*dm2);
				
			}//rocher

		
         newtree->Fill();
        }//event loop

}

void setAddresses(TTree* tree){

    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("passedFullSelection",1);
    tree->SetBranchStatus("passedTrig",1);
    tree->SetBranchStatus("triggersPassed",1);
    tree->SetBranchStatus("lep_id",1);
    tree->SetBranchStatus("lep_tightId",1);
    tree->SetBranchStatus("lep_pt",1);
    tree->SetBranchStatus("lep_eta",1);
    tree->SetBranchStatus("lep_phi",1);
    tree->SetBranchStatus("lep_mass",1);

    tree->SetBranchStatus("lepFSR_pt",1);
    tree->SetBranchStatus("lepFSR_eta",1);
    tree->SetBranchStatus("lepFSR_phi",1);
    tree->SetBranchStatus("lepFSR_mass",1);
    
    tree->SetBranchStatus("lep_RelIso",1);
    tree->SetBranchStatus("lep_pterr",1);
    tree->SetBranchStatus("lep_pterrold",1);
    tree->SetBranchStatus("lep_Sip",1);
    tree->SetBranchStatus("lep_dataMC",1);
    tree->SetBranchStatus("lep_genindex",1);
//    tree->SetBranchStatus("lep_ecalTrkEnergyPostCorr",1);
//    tree->SetBranchStatus("lep_ecalTrkEnergyPreCorr",1);

    tree->SetBranchAddress("passedFullSelection",&passedFullSelection);
    tree->SetBranchAddress("passedTrig",&passedTrig);
    tree->SetBranchAddress("lep_tightId", &lep_tightId);
    tree->SetBranchAddress("lep_id", &lep_id);
    tree->SetBranchAddress("lep_pt",&lep_pt);
    tree->SetBranchAddress("lepFSR_pt",&lepFSR_pt);
    tree->SetBranchAddress("lep_eta",&lep_eta);
    tree->SetBranchAddress("lep_phi",&lep_phi);
    tree->SetBranchAddress("lep_mass",&lep_mass);
    tree->SetBranchAddress("lepFSR_eta",&lepFSR_eta);
    tree->SetBranchAddress("lepFSR_phi",&lepFSR_phi);
    tree->SetBranchAddress("lepFSR_mass",&lepFSR_mass);
    tree->SetBranchAddress("lep_RelIso",&lep_RelIso);
    tree->SetBranchAddress("lep_pterr",&lep_pterr);
    tree->SetBranchAddress("lep_pterrold",&lep_pterrold);
    tree->SetBranchAddress("lep_Sip", &lep_Sip); 
    tree->SetBranchAddress("lep_dataMC", &lep_dataMC); 
    tree->SetBranchAddress("lep_genindex", &lep_genindex);
//    tree->SetBranchAddress("lep_ecalTrkEnergyPostCorr",&lep_ecalTrkEnergyPostCorr);
//    tree->SetBranchAddress("lep_ecalTrkEnergyPreCorr",&lep_ecalTrkEnergyPreCorr);
    tree->SetBranchStatus("Run",1);
    tree->SetBranchStatus("LumiSect",1);
    tree->SetBranchStatus("Event",1);
    tree->SetBranchAddress("Run",&Run);
    tree->SetBranchAddress("LumiSect",&LumiSect);
    tree->SetBranchAddress("Event",&Event);
    tree->SetBranchStatus("met",1);
    tree->SetBranchAddress("met", &met);
    tree->SetBranchStatus("GENZ_mass",1);
    tree->SetBranchAddress("GENZ_mass", &GENZ_mass);

    tree->SetBranchStatus("GENlep_pt",1);
    tree->SetBranchAddress("GENlep_pt", &GENlep_pt);
    tree->SetBranchStatus("GENlep_eta",1);
    tree->SetBranchAddress("GENlep_eta", &GENlep_eta);
    tree->SetBranchStatus("GENlep_phi",1);
    tree->SetBranchAddress("GENlep_phi", &GENlep_phi);
    tree->SetBranchStatus("GENlep_mass",1);
    tree->SetBranchAddress("GENlep_mass", &GENlep_mass);

    tree->SetBranchStatus("nFSRPhotons",1);
    tree->SetBranchAddress("nFSRPhotons", &nFSRPhotons);

    tree->SetBranchStatus("lep_ecalDriven", 1);
    tree->SetBranchAddress("lep_ecalDriven", &lep_ecalDriven);

    tree->SetBranchStatus("track_vtx_pt", 1);
    tree->SetBranchAddress("track_vtx_pt", &track_vtx_pt);
    tree->SetBranchStatus("track_vtx_eta", 1);
    tree->SetBranchAddress("track_vtx_eta", &track_vtx_eta);
    tree->SetBranchStatus("track_vtx_phi", 1);
    tree->SetBranchAddress("track_vtx_phi", &track_vtx_phi);
    tree->SetBranchStatus("track_vtx_pterr", 1);
    tree->SetBranchAddress("track_vtx_pterr", &track_vtx_pterr);
    
//    tree->SetBranchStatus("track_nl", 1);
//    tree->SetBranchAddress("track_nl", &track_nl);
//    tree->SetBranchStatus("track_type", 1);
//    tree->SetBranchAddress("track_type", &track_type);
//    tree->SetBranchStatus("track_gen_particle", 1);
//    tree->SetBranchAddress("track_gen_particle", &track_gen_particle);
//	tree->SetBranchStatus("track_gen_pt",1);
//	tree->SetBranchAddress("track_gen_pt",&track_gen_pt);

    tree->SetBranchStatus("track_bs_pt", 1);
    tree->SetBranchStatus("track_bs_eta", 1);
    tree->SetBranchStatus("track_bs_phi", 1);
    tree->SetBranchStatus("track_bs_pterr", 1);
    tree->SetBranchAddress("track_bs_pt", &track_bs_pt);
    tree->SetBranchAddress("track_bs_eta", &track_bs_eta);
    tree->SetBranchAddress("track_bs_phi", &track_bs_phi);
    tree->SetBranchAddress("track_bs_pterr", &track_bs_pterr);
    tree->SetBranchStatus("passLepindex", 1);
    tree->SetBranchAddress("passLepindex", &passZlepsIndex);
	
	tree->SetBranchStatus("eleUncorr_pt", 1);
	tree->SetBranchStatus("eleUncorr_mass", 1);
	tree->SetBranchAddress("eleUncorr_pt", &eleUncorr_pt);
	tree->SetBranchAddress("eleUncorr_mass", &eleUncorr_mass);
 	
	tree->SetBranchStatus("muonUncorr_pt",1);
	tree->SetBranchAddress("muonUncorr_pt", &muon_uncorr_pt);

	tree->SetBranchStatus("muon_scale_sf",1);
	tree->SetBranchAddress("muon_scale_sf",&muon_scale_sf);

	tree->SetBranchStatus("mu_pt_roc",1);
	tree->SetBranchStatus("mu_eta_roc",1);
	tree->SetBranchStatus("mu_phi_roc",1);
	tree->SetBranchStatus("mu_gen_pt_roc",1);
	tree->SetBranchStatus("mu_genparticle_roc",1);
	tree->SetBranchStatus("mu_nl_roc",1);
	tree->SetBranchStatus("mu_q_roc",1);
	tree->SetBranchStatus("mu_tracktype_roc",1);

	tree->SetBranchAddress("mu_pt_roc",&mu_pt_roc);
	tree->SetBranchAddress("mu_eta_roc",&mu_eta_roc);
	tree->SetBranchAddress("mu_phi_roc",&mu_phi_roc);
	tree->SetBranchAddress("mu_gen_pt_roc",&mu_gen_pt_roc);
	tree->SetBranchAddress("mu_genparticle_roc",&mu_genparticle_roc);
	tree->SetBranchAddress("mu_nl_roc",&mu_nl_roc);
	tree->SetBranchAddress("mu_q_roc",&mu_q_roc);
	tree->SetBranchAddress("mu_tracktype_roc",&mu_tracktype_roc);

}


