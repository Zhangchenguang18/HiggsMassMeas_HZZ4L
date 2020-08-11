#include "ZZ4LConfig.h"
#include "ZZ4LAnalysisTree.h"
#include "LeptonEfficiency.h"
#include "PileupWeight.h"
#include "computeAngles.cc"
#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h"
#include "ZZMatrixElement/MELA/interface/Mela.h"
#include "KinZfitter/KinZfitter.h"

using namespace std;
using namespace MEMNames;

/////////////////////

void SetNewTree(TTree* newtree);
void ReadTree(TTree* tree, TTree* & newtree, TString filename);
double GetCorr(int id, double pt, double corrs[6]);
double GetAllCorr(int id, double pt, double eta);


TString filename;
bool debug;


MEMs* combinedMEM;
Mela* myMela;

KinZfitter *kinZfitter;
//EBE LUT
TFile* f_corr_e_1 = new TFile("KinZfitter/2017_MC_e1.root","READ");
TFile* f_corr_e_3 = new TFile("KinZfitter/2017_MC_e3.root","READ");
TFile* f_corr_mu = new TFile("KinZfitter/2017_MC_mu.root","READ");


TH2F* el_corr_1 = (TH2F*)f_corr_e_1->Get("e1");
TH2F* el_corr_3 = (TH2F*)f_corr_e_3->Get("e3");

TH2F* mu_corr = (TH2F*)f_corr_mu->Get("mu");//new LUT e1 for ecal electron, e3 for tracker electron, no e2

TAxis* x_elpTaxis_1 = el_corr_1->GetXaxis(); TAxis* y_eletaaxis_1 = el_corr_1->GetYaxis();
double maxPtEl_1 = x_elpTaxis_1->GetXmax(); double minPtEl_1 = x_elpTaxis_1->GetXmin();

TAxis* x_elpTaxis_3 = el_corr_3->GetXaxis(); TAxis* y_eletaaxis_3 = el_corr_3->GetYaxis();
double maxPtEl_3 = x_elpTaxis_3->GetXmax(); double minPtEl_3 = x_elpTaxis_3->GetXmin();

TAxis* x_mupTaxis = mu_corr->GetXaxis(); TAxis* y_muetaaxis = mu_corr->GetYaxis();
double maxPtMu = x_mupTaxis->GetXmax(); double minPtMu = x_mupTaxis->GetXmin();
//Energy Scale LUT
TFile* f = new TFile("EnergyScaleLUT/2017_Mu.root","READ");
TH2D* LUT_Mu = (TH2D*)f->Get("LUT");
TFile* f_ = new TFile("EnergyScaleLUT/2017_Ele.root","READ");
TH2D* LUT_Ele = (TH2D*)f_->Get("LUT");

int main(int argc, char *argv[])
{    
     
  debug = false;     

  if(argc > 6)  {
      cout<<argv[0]<<" filename "<<argv[1]<<" outfile "<<argv[2]<<" isData "<<argv[3]<<endl;
      return -1;
    }

  //gSystem->Load("$CMSSW_BASE/lib/slc6_amd64_gcc530/pluginUFHZZAnalysisRun2UFHZZ4LAna.so");  

  combinedMEM = new MEMs(13.0,125,"CTEQ6L",false);
  myMela = combinedMEM->m_MELA;

  /////////////////////

  filename = argv[1];
  TString outfilename = argv[2];
  
  if(atof(argv[3])>0) { isData = true; }
  else {isData = false; }

  if(atof(argv[4])>0) { 
      job=strtol(argv[4], NULL, 10); njobs=strtol(argv[5], NULL, 10);
      cout<<"job "<<job<<" of "<<njobs<<endl;
      outfilename+="_";
      outfilename+=argv[4];
  }

  kinZfitter = new KinZfitter(isData);

  TFile* infile = TFile::Open("root://cmsio5.rc.ufl.edu:1094/"+filename+".root");
  if(!infile) infile = new TFile(filename+".root");
  TTree* tree; 
  
  tree = (TTree*) infile->Get("Ana/passedEvents");
  if(!tree) tree = (TTree*) infile->Get("passedEvents");
  if(!tree) tree = (TTree*) infile->Get("selectedEvents");
  if(!tree) { cout<<"ERROR could not find the tree for "<<filename<<endl; return -1;}

  TH1F* sumWeights = (TH1F*) infile->Get("Ana/sumWeights");
  sumW = 1.0;
  if(sumWeights) sumW = sumWeights->GetBinContent(1); 

  cout<<"sumW is "<<sumW<<endl;

  // read tree     

  TString name = outfilename;
  TFile* tmpFile =  new TFile(name+".root","RECREATE");
  TTree* newtree = new TTree("passedEvents","passedEvents");

  if(debug)cout<<"start setting new tree "<<endl;

  SetNewTree(newtree);

  if(debug)cout<<"start reading tree "<<endl;

  ReadTree(tree, newtree, filename);

  if(debug)cout<<"end reading tree"<<endl;

  tmpFile->cd();

  newtree->Write("passedEvents",TObject::kOverwrite);

  tmpFile->Close(); 


}


void ReadTree(TTree* tree, TTree* & newtree, TString filename){
    
    ZZ4LAnalysisTree::setAddresses(tree, filename);

    float npass = 0.0;
    float sumpuweight = 0.0;

    if(debug) cout<<"start looping"<<endl;

    double firstevt=0; double lastevt=tree->GetEntries();
    if (job>0) {
        firstevt = tree->GetEntries()*(job-1)/njobs;
        lastevt = tree->GetEntries()*(job)/(njobs);
        cout << "firstevt: " << firstevt << ", lastevt: " << lastevt << endl;
    }


int count = 0;    
    for(int evt=0; evt < tree->GetEntries() ; evt++) { //event loop                              
   //if(evt>1000)break; 
        if (evt<firstevt) continue;
        if (evt>lastevt) continue;

        if(evt%1000==0) cout<<"Event "<<evt<<"/"<<tree->GetEntries()<<endl;                                     

        tree->GetEntry(evt);
//        if (Event != 273157367) continue;
//        if (Event != 163484) continue;
//
//count++;  if (count > 100) continue; 
        if (isData) {
            passTrig = passedTrig;
        } else {
            passTrig = true;
        }

        if((*lep_id).size()<4) continue;
        unsigned int Nlep = (*lep_id).size();
        if (debug) cout<<Nlep<<" leptons in total"<<endl;

        bool foundHiggsCandidate=false;

/*
        for(unsigned int i=0; i<Nlep; i++){

           if (abs((*lep_id)[i]) != 13 ) continue;
           if (abs((*lep_eta)[i]) < 0.9 )                                 (*lep_pt)[i] = (*lep_pt)[i]/(1-0.00070374);
           if (abs((*lep_eta)[i]) < 1.8 && abs((*lep_eta)[i]) > 0.9 ) (*lep_pt)[i] = (*lep_pt)[i]/(1-0.0015881);
           if (abs((*lep_eta)[i]) < 2.4 && abs((*lep_eta)[i]) > 1.8 ) (*lep_pt)[i] = (*lep_pt)[i]/(1-0.0029359);

           }
*/

        if (redoEventSelection) {
            // 2 OSSF Pairs
            bool properLep_ID = false; int Nmm = 0; int Nmp = 0; int Nem = 0; int Nep = 0;
            for(unsigned int i =0; i<Nlep; i++) {
                if((*lep_id)[i]==-13) Nmm = Nmm+1;
                if((*lep_id)[i]==13) Nmp = Nmp+1;
            }
            for(unsigned int i =0; i<Nlep; i++) {
                if((*lep_id)[i]==-11) Nem = Nem+1;
                if((*lep_id)[i]==11) Nep = Nep+1;
            }
            
            if(Nmm>=2 && Nmp>=2) properLep_ID = true; //4mu
            if(Nem>=2 && Nep>=2) properLep_ID = true; //4e
            if(Nmm>0 && Nmp>0 && Nem>0 && Nep>0) properLep_ID = true; //2e2mu
            
            if(!properLep_ID) continue;
            
            // First, make all Z candidates including any FSR photons
            const double Zmass = 91.1876;
            int n_Zs=0;
            vector<int> Z_lepindex1;
            vector<int> Z_lepindex2;
            vector<float> Z_pt, Z_eta, Z_phi, Z_mass;

            for(unsigned int i=0; i<Nlep; i++){
                for(unsigned int j=i+1; j<Nlep; j++){
                    // same flavor opposite charge
                    if(((*lep_id)[i]+(*lep_id)[j])!=0) continue;
                    
                    TLorentzVector li, lj;
                    li.SetPtEtaPhiM((*lep_pt)[i],(*lep_eta)[i],(*lep_phi)[i],(*lep_mass)[i]);
                    lj.SetPtEtaPhiM((*lep_pt)[j],(*lep_eta)[j],(*lep_phi)[j],(*lep_mass)[j]);
                    
                    TLorentzVector lifsr, ljfsr;
                    lifsr.SetPtEtaPhiM((*lepFSR_pt)[i],(*lepFSR_eta)[i],(*lepFSR_phi)[i],(*lepFSR_mass)[i]);
                    ljfsr.SetPtEtaPhiM((*lepFSR_pt)[j],(*lepFSR_eta)[j],(*lepFSR_phi)[j],(*lepFSR_mass)[j]);
                    
                    TLorentzVector liljfsr = lifsr+ljfsr;
                    
                    if (debug) {
                        cout<<"OSSF pair: i="<<i<<" id1="<<(*lep_id)[i]<<" j="<<j<<" id2="<<(*lep_id)[j]<<" pt1: "
                            <<lifsr.Pt()<<" pt2: "<<ljfsr.Pt()<<" M: "<<liljfsr.M()<<endl;    
                    }
                    
                    TLorentzVector Z, Z_noFSR;
                    Z = lifsr+ljfsr;
                    Z_noFSR = li+lj;
                    
                    if (debug) cout<<"this Z mass: "<<Z.M()<<" mZ2Low: "<<mZ2Low<<endl;
                    
                    if (Z.M()>0.0) {
                        n_Zs++;
                        Z_pt.push_back(Z.Pt());
                        Z_eta.push_back(Z.Eta());
                        Z_phi.push_back(Z.Phi());
                        Z_mass.push_back(Z.M());
                        Z_lepindex1.push_back(i);
                        Z_lepindex2.push_back(j);
                        if (debug) cout<<" add Z_lepindex1: "<<i<<" Z_lepindex2: "<<j<<endl;
                    }
                    
                } // lep i
            } // lep j
            
            // Consider all ZZ candidates
            TLorentzVector Z1Vec, Z2Vec, HVec;
            double minZ1DeltaM_SR=9999.9; double minZ1DeltaM_CR=99999.9;
            double max_D_bkg_kin_SR=0.0; double max_D_bkg_kin_CR=0.0;
            bool foundSRCandidate=false;
            
            passedZ4lSelection=false;
            passedFullSelection=false;
            
            vector<int> lep_Hindex, Z_Hindex;
            for (int i=0; i<4; i++) {
                if (i<2) Z_Hindex.push_back(-1);
                lep_Hindex.push_back(-1);
            }
            
            for (int i=0; i<n_Zs; i++) {
                for (int j=i+1; j<n_Zs; j++) {
                    
                    int i1 = Z_lepindex1[i]; int i2 = Z_lepindex2[i];                            
                    int j1 = Z_lepindex1[j]; int j2 = Z_lepindex2[j];                            
                    
                    if (i1 == j1 || i1 == j2 || i2 == j1 || i2 == j2) continue;
                    
                    TLorentzVector lep_i1, lep_i2, lep_j1, lep_j2;
                    lep_i1.SetPtEtaPhiM((*lepFSR_pt)[i1],(*lepFSR_eta)[i1],(*lepFSR_phi)[i1],(*lepFSR_mass)[i1]);
                    lep_i2.SetPtEtaPhiM((*lepFSR_pt)[i2],(*lepFSR_eta)[i2],(*lepFSR_phi)[i2],(*lepFSR_mass)[i2]);
                    lep_j1.SetPtEtaPhiM((*lepFSR_pt)[j1],(*lepFSR_eta)[j1],(*lepFSR_phi)[j1],(*lepFSR_mass)[j1]);
                    lep_j2.SetPtEtaPhiM((*lepFSR_pt)[j2],(*lepFSR_eta)[j2],(*lepFSR_phi)[j2],(*lepFSR_mass)[j2]);
                    
                    TLorentzVector lep_i1_nofsr, lep_i2_nofsr, lep_j1_nofsr, lep_j2_nofsr;
                    lep_i1_nofsr.SetPtEtaPhiM((*lep_pt)[i1],(*lep_eta)[i1],(*lep_phi)[i1],(*lep_mass)[i1]);
                    lep_i2_nofsr.SetPtEtaPhiM((*lep_pt)[i2],(*lep_eta)[i2],(*lep_phi)[i2],(*lep_mass)[i2]);
                    lep_j1_nofsr.SetPtEtaPhiM((*lep_pt)[j1],(*lep_eta)[j1],(*lep_phi)[j1],(*lep_mass)[j1]);
                    lep_j2_nofsr.SetPtEtaPhiM((*lep_pt)[j2],(*lep_eta)[j2],(*lep_phi)[j2],(*lep_mass)[j2]);
                    
                    TLorentzVector Zi, Zj;
                    Zi.SetPtEtaPhiM(Z_pt[i],Z_eta[i],Z_phi[i],Z_mass[i]);
                    Zj.SetPtEtaPhiM(Z_pt[j],Z_eta[j],Z_phi[j],Z_mass[j]);
                    
                    if (debug) {cout<<"ZZ candidate Zi->M() "<<Zi.M()<<" Zj->M() "<<Zj.M()<<endl;}
                    
                    TLorentzVector Z1, Z2;
                    int Z1index, Z2index;
                    int Z1_lepindex[2] = {0,0};
                    int Z2_lepindex[2] = {0,0};
                    double Z1DeltaM;
                    
                    if (abs(Zi.M()-Zmass)<abs(Zj.M()-Zmass)) { 
                        Z1index = i; Z2index = j;
                        Z1 = Zi; Z2 = Zj;                 
                        if (lep_i1.Pt()>lep_i2.Pt()) { Z1_lepindex[0] = i1;  Z1_lepindex[1] = i2; }
                        else { Z1_lepindex[0] = i2;  Z1_lepindex[1] = i1; }                
                        if (lep_j1.Pt()>lep_j2.Pt()) { Z2_lepindex[0] = j1;  Z2_lepindex[1] = j2; } 
                        else { Z2_lepindex[0] = j2;  Z2_lepindex[1] = j1; }                
                        Z1DeltaM = abs(Zi.M()-Zmass); 
                    }
                    else { 
                        Z1index = j; Z2index = i;
                        Z1 = Zj; Z2 = Zi; 
                        if (lep_j1.Pt()>lep_j2.Pt()) { Z1_lepindex[0] = j1;  Z1_lepindex[1] = j2; }
                        else { Z1_lepindex[0] = j2;  Z1_lepindex[1] = j1; }
                        if (lep_i1.Pt()>lep_i2.Pt()) { Z2_lepindex[0] = i1;  Z2_lepindex[1] = i2; }
                        else { Z2_lepindex[0] = i2;  Z2_lepindex[1] = i1; }
                        Z1DeltaM = abs(Zj.M()-Zmass); 
                    }
                    
                    // Check isolation cut (without FSR ) for Z1 leptons
                    if ((*lep_RelIsoNoFSR)[Z1_lepindex[0]]>((abs((*lep_id)[Z1_lepindex[0]])==11) ? isoCutEl : isoCutMu)) continue;
                    if ((*lep_RelIsoNoFSR)[Z1_lepindex[1]]>((abs((*lep_id)[Z1_lepindex[1]])==11) ? isoCutEl : isoCutMu)) continue;
                    // Check tight ID cut for Z1 leptons
                    if (!((*lep_tightId)[Z1_lepindex[0]])) continue;
                    if (!((*lep_tightId)[Z1_lepindex[1]])) continue;
                    
                    // Check Leading and Subleading pt Cut
                    vector<double> allPt;
                    allPt.push_back(lep_i1_nofsr.Pt()); allPt.push_back(lep_i2_nofsr.Pt());
                    allPt.push_back(lep_j1_nofsr.Pt()); allPt.push_back(lep_j2_nofsr.Pt());
                    std::sort(allPt.begin(), allPt.end());
                    if (debug) cout<<" leading pt: "<<allPt[3]<<" cut: "<<leadingPtCut
                                   <<" subleadingPt: "<<allPt[2]<<" cut: "<<subleadingPtCut<<endl;
                    if (allPt[3]<leadingPtCut || allPt[2]<subleadingPtCut ) continue;
                    
                    // Check dR(li,lj)>0.02 for any i,j
                    vector<double> alldR;
                    alldR.push_back(lep_i1_nofsr.DeltaR(lep_i2_nofsr));
                    alldR.push_back(lep_i1_nofsr.DeltaR(lep_j1_nofsr));
                    alldR.push_back(lep_i1_nofsr.DeltaR(lep_j2_nofsr));
                    alldR.push_back(lep_i2_nofsr.DeltaR(lep_j1_nofsr));
                    alldR.push_back(lep_i2_nofsr.DeltaR(lep_j2_nofsr));
                    alldR.push_back(lep_j1_nofsr.DeltaR(lep_j2_nofsr));
                    if (debug) cout<<" minDr: "<<*min_element(alldR.begin(),alldR.end())<<endl;
                    if (*min_element(alldR.begin(),alldR.end())<0.02) continue;
                    
                    // Check M(l+,l-)>4.0 GeV for any OS pair
                    // Do not include FSR photons
                    vector<double> allM;
                    TLorentzVector i1i2;
                    i1i2 = (lep_i1_nofsr)+(lep_i2_nofsr); allM.push_back(i1i2.M());
                    TLorentzVector j1j2;
                    j1j2 = (lep_j1_nofsr)+(lep_j2_nofsr); allM.push_back(j1j2.M());            
                    
                    if ((*lep_id)[i1]*(*lep_id)[j1]<0) {
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
                    if (debug) cout<<" min m(l+l-): "<<*min_element(allM.begin(),allM.end())<<endl;
                    if (*min_element(allM.begin(),allM.end())<4.0) { continue;}
                    
                    // Check the "smart cut": !( |mZa-mZ| < |mZ1-mZ| && mZb<12)
                    // only for 4mu or 4e ZZ candidates
                    bool passSmartCut=true;
                    if ( abs((*lep_id)[i1])==abs((*lep_id)[j1])) {
                        TLorentzVector Za, Zb;
                        if ((*lep_id)[i1]==(*lep_id)[j1]) {                  
                            Za = (lep_i1)+(lep_j2);
                            Zb = (lep_i2)+(lep_j1);                    
                        } else {
                            Za = (lep_i1)+(lep_j1);
                            Zb = (lep_i2)+(lep_j2);
                        }                
                        if ( abs(Za.M()-Zmass)<abs(Zb.M()-Zmass) ) {
                            if (debug) cout<<"abs(Za.M()-Zmass)-abs(Z1.M()-Zmass): "
                                           <<abs(Za.M()-Zmass)-abs(Z1.M()-Zmass)<<" Zb.M(): "<<Zb.M()<<endl;
                            if ( abs(Za.M()-Zmass)<abs(Z1.M()-Zmass) && Zb.M()<mZ2Low ) passSmartCut=false;
                        }
                        else {
                            if (debug) cout<<"abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass): "
                                           <<abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass)<<" Za.M(): "<<Za.M()<<endl;
                            if ( abs(Zb.M()-Zmass)<abs(Z1.M()-Zmass) && Za.M()<mZ2Low ) passSmartCut=false;
                        }
                        
                    }
                    if (!passSmartCut) continue;
                    
                    if (debug) cout<<" massZ1: "<<Z1.M()<<" massZ2: "<<Z2.M()<<endl;
                    if ( (Z1.M() < mZ1Low) || (Z1.M() > mZ1High) || (Z2.M() < mZ2Low) || (Z2.M() > mZ2High) ) continue;
                    if (debug) cout<<" pass Z mass cuts"<<endl;
                    
                    
                    // Signal region if Z2 leptons are both tight ID Iso
                    bool signalRegion=true;
                    if ((*lep_RelIsoNoFSR)[Z2_lepindex[0]]>((abs((*lep_id)[Z2_lepindex[0]])==11) ? isoCutEl : isoCutMu)) signalRegion=false;
                    if ((*lep_RelIsoNoFSR)[Z2_lepindex[1]]>((abs((*lep_id)[Z2_lepindex[1]])==11) ? isoCutEl : isoCutMu)) signalRegion=false;
                    if (!((*lep_tightId)[Z2_lepindex[0]])) signalRegion=false; // checking tight lepton ID
                    if (!((*lep_tightId)[Z2_lepindex[1]])) signalRegion=false; // checking tight lepton ID          
                    
                    if (debug) cout<<"signalRegion? "<<signalRegion<<endl;
                    
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
                    
                    tmpIDs.push_back((*lep_id)[Z1_lepindex[0]]); tmpIDs.push_back((*lep_id)[Z1_lepindex[1]]);
                    tmpIDs.push_back((*lep_id)[Z2_lepindex[0]]); tmpIDs.push_back((*lep_id)[Z2_lepindex[1]]);
                    
                    if (debug) cout<<"good ZZ candidate? "<<endl;
                    
                    double me_0plus_JHU_tmp, me_qqZZ_MCFM_tmp;
                    combinedMEM->computeME(MEMNames::kSMHiggs, MEMNames::kJHUGen, P4s, tmpIDs, me_0plus_JHU_tmp); // higgs, vector algebra, JHUgen
                    combinedMEM->computeME(MEMNames::kqqZZ, MEMNames::kMCFM, P4s, tmpIDs, me_qqZZ_MCFM_tmp); // background, vector algebra, MCFM 
                    double D_bkg_kin_tmp = me_0plus_JHU_tmp / (me_0plus_JHU_tmp + me_qqZZ_MCFM_tmp);
                    
                    if (debug) cout<<"good ZZ candidate, D_bkg_kin: "<<D_bkg_kin_tmp<<" max D_bkg_kin SR: "
                                   <<max_D_bkg_kin_SR<<" max D_bkg_kin CR: "<<max_D_bkg_kin_CR<<endl;
                    
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
                        
                        if ( (!same4l && D_bkg_kin_tmp>max_D_bkg_kin_SR) || (same4l && Z1DeltaM<=minZ1DeltaM_SR)) {                 
                            
                            max_D_bkg_kin_SR = D_bkg_kin_tmp;
                            minZ1DeltaM_SR = Z1DeltaM;
                            
                            Z_Hindex[0] = Z1index;
                            lep_Hindex[0] = Z1_lepindex[0];
                            lep_Hindex[1] = Z1_lepindex[1];
                            
                            Z_Hindex[1] = Z2index;
                            lep_Hindex[2] = Z2_lepindex[0];
                            lep_Hindex[3] = Z2_lepindex[1];
                            
                            Z1Vec = Z1; Z2Vec = Z2; HVec = Z1+Z2;                   
                            massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); mass4l = HVec.M();
                            
                            if (debug) cout<<" new best candidate SR: mass4l: "<<HVec.M()<<endl;
                            if (HVec.M()>m4lLowCut)  {
                                foundHiggsCandidate=true;                    
                                foundSRCandidate=true;
                            }
                        }
                    } else if (!foundSRCandidate) { // Control regions get second priority
                        
                        //if ( D_bkg_kin_tmp>max_D_bkg_kin_CR ) {
                        if ( (!same4l && D_bkg_kin_tmp>max_D_bkg_kin_CR) || (same4l && Z1DeltaM<=minZ1DeltaM_CR)) {                 
                            
                            max_D_bkg_kin_CR = D_bkg_kin_tmp;
                            minZ1DeltaM_CR = Z1DeltaM;
                            
                            Z_Hindex[0] = Z1index;
                            lep_Hindex[0] = Z1_lepindex[0];
                            lep_Hindex[1] = Z1_lepindex[1];
                            
                            Z_Hindex[1] = Z2index;
                            lep_Hindex[2] = Z2_lepindex[0];
                            lep_Hindex[3] = Z2_lepindex[1];
                        
                            Z1Vec = Z1; Z2Vec = Z2; HVec = Z1+Z2;                   
                            massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); mass4l = HVec.M(); pT4l = HVec.Pt();
                            
                            if (debug) cout<<" new best candidate CR: mass4l: "<<HVec.M()<<endl;
                            if (HVec.M()>m4lLowCut) foundHiggsCandidate=true;                    
                        }
                    }
                    
                    if (debug) cout<<"Z_Hindex[0]: "<<Z_Hindex[0]<<" lep_Hindex[0]: "<<lep_Hindex[0]<<" lep_Hindex[1]: "<<lep_Hindex[1]
                                   <<"Z_Hindex[1]: "<<Z_Hindex[1]<<" lep_Hindex[2]: "<<lep_Hindex[2]<<" lep_Hindex[3]: "<<lep_Hindex[3]<<endl;
                    
                } // Zj
            } // Zi
        }
                
        int nZXCRFailedLeptons=0;
        if ( (redoEventSelection&&foundHiggsCandidate) || passedFullSelection || passedZXCRSelection  || passedZ4lSelection) {

/*** correct muon pt scale ***/
/*
            for(unsigned int i=0; i<4; i++){
               int index = lep_Hindex[i];
               if (abs((*lep_id)[index]) != 13 ) continue;
               if (abs((*lep_eta)[index]) < 0.9 )                                 (*lepFSR_pt)[index] = (*lepFSR_pt)[index]/(1-0.000562032);
               if (abs((*lep_eta)[index]) < 1.8 && abs((*lep_eta)[index]) > 0.9 ) (*lepFSR_pt)[index] = (*lepFSR_pt)[index]/(1-0.00110609);
               if (abs((*lep_eta)[index]) < 2.4 && abs((*lep_eta)[index]) > 1.8 ) (*lepFSR_pt)[index] = (*lepFSR_pt)[index]/(1-0.00177786);

               }

            for(unsigned int i=0; i<4; i++){
               int index = lep_Hindex[i];
               if (abs((*lep_id)[index]) != 13 ) continue;

               if (abs((*lep_eta)[index]) < 0.9 )                                 (*lep_pt)[index] = (*lep_pt)[index]/(1-0.000562032);
               if (abs((*lep_eta)[index]) < 1.8 && abs((*lep_eta)[index]) > 0.9 ) (*lep_pt)[index] = (*lep_pt)[index]/(1-0.00110609);
               if (abs((*lep_eta)[index]) < 2.4 && abs((*lep_eta)[index]) > 1.8 ) (*lep_pt)[index] = (*lep_pt)[index]/(1-0.00177786);

               }

*/
            TLorentzVector Lep1, Lep2, Lep3, Lep4,  Jet1, Jet2;            
            TLorentzVector nullFourVector(0, 0, 0, 0);                 
            Lep1.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[0]],(*lepFSR_eta)[lep_Hindex[0]],(*lepFSR_phi)[lep_Hindex[0]],(*lepFSR_mass)[lep_Hindex[0]]);
            Lep2.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[1]],(*lepFSR_eta)[lep_Hindex[1]],(*lepFSR_phi)[lep_Hindex[1]],(*lepFSR_mass)[lep_Hindex[1]]);
            Lep3.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[2]],(*lepFSR_eta)[lep_Hindex[2]],(*lepFSR_phi)[lep_Hindex[2]],(*lepFSR_mass)[lep_Hindex[2]]);
            Lep4.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[3]],(*lepFSR_eta)[lep_Hindex[3]],(*lepFSR_phi)[lep_Hindex[3]],(*lepFSR_mass)[lep_Hindex[3]]);


            if (redoEventSelection) {                
                for(unsigned int i = 0; i <= 3; i++) {
                    if (!(abs((*lep_id)[lep_Hindex[i]])==11 && ((*lep_tightId)[lep_Hindex[i]] && (*lep_RelIsoNoFSR)[lep_Hindex[i]]<isoCutEl)) &&
                        !(abs((*lep_id)[lep_Hindex[i]])==13 && ((*lep_tightId)[lep_Hindex[i]] && (*lep_RelIsoNoFSR)[lep_Hindex[i]]<isoCutMu))){ 
                        nZXCRFailedLeptons++; 
                    }
                }
                
                if (debug) cout << nZXCRFailedLeptons<<" failing leptons in higgs candidate"<<endl;
                
                if (nZXCRFailedLeptons>0) { // at least one lepton has failed 
                    passedZ4lZXCRSelection = true;
                    if ((Lep1+Lep2).M() > mZ2Low && passedTrig) passedZXCRSelection = true;
                } else { //  signal region candidate                    
                    passedZ4lSelection = true;
                    if((Lep1+Lep2).M() > mZ2Low && passedTrig) passedFullSelection = true;
                }
                if (passedZ4lSelection) npass+=1;
                
                /*
                njets_pt30_eta4p7=0;
                for( unsigned int k = 0; k < (*jet_iscleanH4l).size(); k++) {
                    if ((*jet_pt)[(*jet_iscleanH4l)[k]]>30.0 && abs((*jet_eta)[(*jet_iscleanH4l)[k]])<4.7) {
                        njets_pt30_eta4p7+=1;
                    }
                }
                */

                if (debug) cout<<njets_pt30_eta4p7<<" jets"<<endl;
            }

            //cout<<evt<<" Run: "<<Run<<" LumiSect "<<LumiSect<<" Event "<<Event<<endl;
            
            idL1 = (*lep_id)[lep_Hindex[0]]; pTL1 = Lep1.Pt(); etaL1 = Lep1.Eta(); pTErrL1 = (*lep_pterr)[lep_Hindex[0]]; mL1 = (*lep_mass)[lep_Hindex[0]];
            idL2 = (*lep_id)[lep_Hindex[1]]; pTL2 = Lep2.Pt(); etaL2 = Lep2.Eta(); pTErrL2 = (*lep_pterr)[lep_Hindex[1]]; mL2 = (*lep_mass)[lep_Hindex[1]];
            idL3 = (*lep_id)[lep_Hindex[2]]; pTL3 = Lep3.Pt(); etaL3 = Lep3.Eta(); pTErrL3 = (*lep_pterr)[lep_Hindex[2]]; mL3 = (*lep_mass)[lep_Hindex[2]];      
            idL4 = (*lep_id)[lep_Hindex[3]]; pTL4 = Lep4.Pt(); etaL4 = Lep4.Eta(); pTErrL4 = (*lep_pterr)[lep_Hindex[3]]; mL4 = (*lep_mass)[lep_Hindex[3]];

            phiL1 = (*lep_phi)[lep_Hindex[0]];
            phiL2 = (*lep_phi)[lep_Hindex[1]];
            phiL3 = (*lep_phi)[lep_Hindex[2]];
            phiL4 = (*lep_phi)[lep_Hindex[3]];

            int genIndex1 = (*lep_genindex)[lep_Hindex[0]];
            int genIndex2 = (*lep_genindex)[lep_Hindex[1]];
            int genIndex3 = (*lep_genindex)[lep_Hindex[2]];
            int genIndex4 = (*lep_genindex)[lep_Hindex[3]];

            pTGENL1 = -999; pTGENL2 = -999; pTGENL3 = -999; pTGENL4 = -999;

            if (genIndex1 >= 0 && genIndex2 >= 0 && genIndex3 >= 0 && genIndex4 >= 0) {

               pTGENL1 = (*GENlep_pt)[genIndex1];
               pTGENL2 = (*GENlep_pt)[genIndex2];
               pTGENL3 = (*GENlep_pt)[genIndex3];
               pTGENL4 = (*GENlep_pt)[genIndex4];

               } 

            vector<TLorentzVector> P4s; vector<int> tmpIDs;             
            P4s.push_back(Lep1); P4s.push_back(Lep2);
            P4s.push_back(Lep3); P4s.push_back(Lep4);
            tmpIDs.push_back(idL1); tmpIDs.push_back(idL2);
            tmpIDs.push_back(idL3); tmpIDs.push_back(idL4);

            TLorentzVector higgs_undec = Lep1+Lep2+Lep3+Lep4;            

            // scale up and dn lep energy scale
            double corr1 = GetAllCorr(idL1, pTL1, etaL1);
            double corr2 = GetAllCorr(idL2, pTL2, etaL2);
            double corr3 = GetAllCorr(idL3, pTL3, etaL3);
            double corr4 = GetAllCorr(idL4, pTL4, etaL4);

            TLorentzVector higgs_undec_up, higgs_undec_dn, higgs_nofsr;
            TLorentzVector Lep1_up, Lep2_up, Lep3_up, Lep4_up;
            TLorentzVector Lep1_dn, Lep2_dn, Lep3_dn, Lep4_dn;
//cout << corr1 << ", " << corr2 << ", " << corr3 << ", " << corr4 << endl;
            Lep1_up.SetPtEtaPhiM(pTL1*(1+corr1),etaL1,phiL1,mL1);
            Lep2_up.SetPtEtaPhiM(pTL2*(1+corr2),etaL2,phiL2,mL2);
            Lep3_up.SetPtEtaPhiM(pTL3*(1+corr3),etaL3,phiL3,mL3);
            Lep4_up.SetPtEtaPhiM(pTL4*(1+corr4),etaL4,phiL4,mL4);

            Lep1_dn.SetPtEtaPhiM(pTL1*(1-corr1),etaL1,phiL1,mL1);
            Lep2_dn.SetPtEtaPhiM(pTL2*(1-corr2),etaL2,phiL2,mL2);
            Lep3_dn.SetPtEtaPhiM(pTL3*(1-corr3),etaL3,phiL3,mL3);
            Lep4_dn.SetPtEtaPhiM(pTL4*(1-corr4),etaL4,phiL4,mL4);

            higgs_undec_up = Lep1_up + Lep2_up + Lep3_up + Lep4_up;
            higgs_undec_dn = Lep1_dn + Lep2_dn + Lep3_dn + Lep4_dn;

            mass4l_up = higgs_undec_up.M();
            mass4l_dn = higgs_undec_dn.M();



            TLorentzVector LepNoFSR1, LepNoFSR2, LepNoFSR3, LepNoFSR4;
            LepNoFSR1.SetPtEtaPhiM((*lep_pt)[lep_Hindex[0]],(*lep_eta)[lep_Hindex[0]],(*lep_phi)[lep_Hindex[0]],(*lep_mass)[lep_Hindex[0]]);
            LepNoFSR2.SetPtEtaPhiM((*lep_pt)[lep_Hindex[1]],(*lep_eta)[lep_Hindex[1]],(*lep_phi)[lep_Hindex[1]],(*lep_mass)[lep_Hindex[1]]);
            LepNoFSR3.SetPtEtaPhiM((*lep_pt)[lep_Hindex[2]],(*lep_eta)[lep_Hindex[2]],(*lep_phi)[lep_Hindex[2]],(*lep_mass)[lep_Hindex[2]]);
            LepNoFSR4.SetPtEtaPhiM((*lep_pt)[lep_Hindex[3]],(*lep_eta)[lep_Hindex[3]],(*lep_phi)[lep_Hindex[3]],(*lep_mass)[lep_Hindex[3]]);
            higgs_nofsr = LepNoFSR1+LepNoFSR2+LepNoFSR3+LepNoFSR4;
         
            mass4l_noFSR = higgs_nofsr.M();
            massZ1_noFSR = (LepNoFSR1+LepNoFSR2).M();
            massZ2_noFSR = (LepNoFSR3+LepNoFSR4).M();


            TLorentzVector LepNoFSR1p, LepNoFSR2p;
            LepNoFSR1p.SetPtEtaPhiM((*lep_pt)[lep_Hindex[0]]+pTErrL1, etaL1, phiL1, mL1);
            LepNoFSR2p.SetPtEtaPhiM((*lep_pt)[lep_Hindex[1]]+pTErrL2, etaL2, phiL2, mL2);
            double dm1 = (LepNoFSR1p+LepNoFSR2).M() - (LepNoFSR1+LepNoFSR2).M();
            double dm2 = (LepNoFSR1+LepNoFSR2p).M() - (LepNoFSR1+LepNoFSR2).M();
            massZ1Err = TMath::Sqrt(dm1*dm1+dm2*dm2);
// FIX
            vector<TLorentzVector> P4sNoFSR; 
//            if (LepNoFSR1.Pt() > LepNoFSR2.Pt()) {
               P4sNoFSR.push_back(LepNoFSR1); P4sNoFSR.push_back(LepNoFSR2);
//               } else {
//                      P4sNoFSR.push_back(LepNoFSR2); P4sNoFSR.push_back(LepNoFSR1);
//                      }
 
            P4sNoFSR.push_back(LepNoFSR3); P4sNoFSR.push_back(LepNoFSR4);
            
            vector<double> lep_pterr_correction;
            if (!isData) {
                for (int i=0; i<4; i++) {
                    if (abs((*lep_id)[lep_Hindex[i]])==13) {
                        int xbin = x_mupTaxis->FindBin((*lep_pt)[lep_Hindex[i]]); 
                        int ybin = y_muetaaxis->FindBin(fabs((*lep_eta)[lep_Hindex[i]]));
                        if((*lep_pt)[lep_Hindex[i]]>minPtMu && (*lep_pt)[lep_Hindex[i]]<maxPtMu ){  
                            lep_pterr_correction.push_back(mu_corr->GetBinContent(xbin,ybin));  
                        } else {
                            lep_pterr_correction.push_back(1.0);
                        }
                    }
                    if (abs((*lep_id)[lep_Hindex[i]])==11) {

                       double pT_e = (*lep_pt)[lep_Hindex[i]];
                       double pTErr_e = (*lep_pterr)[lep_Hindex[i]];
                       double eta_e = (*lep_eta)[lep_Hindex[i]];

                       if ((*lep_ecalDriven)[lep_Hindex[i]]) {

                          int xbin = x_elpTaxis_1->FindBin(fabs(eta_e));
			  int ybin = y_eletaaxis_1->FindBin(pTErr_e/pT_e);
			  if(pT_e<100 && pT_e>7){
				  lep_pterr_correction.push_back(el_corr_1->GetBinContent(xbin,ybin));
				  }else{
				  lep_pterr_correction.push_back(1.0);
				  }
			  
                          //ecal dirven
                          } else {

                                 int xbin = x_elpTaxis_3->FindBin(pT_e);
                                 int ybin = y_eletaaxis_3->FindBin(fabs(eta_e));
                                 if(pT_e > minPtEl_3 && pT_e < maxPtEl_3 ){
                                   lep_pterr_correction.push_back(el_corr_3->GetBinContent(xbin,ybin));
                                   } else {
                                          lep_pterr_correction.push_back(1.0);
                                          }
                                 }

/*
                        int xbin = x_elpTaxis->FindBin((*lep_pt)[lep_Hindex[i]]); 
                        int ybin = y_eletaaxis->FindBin(fabs((*lep_eta)[lep_Hindex[i]]));
                        if((*lep_pt)[lep_Hindex[i]]>minPtEl && (*lep_pt)[lep_Hindex[i]]<maxPtEl ){  
                            lep_pterr_correction.push_back(el_corr->GetBinContent(xbin,ybin));  
                        } else {
                            lep_pterr_correction.push_back(1.0);
                        }

                    // overwrite correction if election is non-EcalDriven, hardcoded
                    if ((*lep_ecalDriven)[lep_Hindex[i]] == 0) {

                       if (fabs((*lep_eta)[lep_Hindex[i]]) < 1.44) lep_pterr_correction[i] = 2.44838;
                       if (fabs((*lep_eta)[lep_Hindex[i]]) > 1.44 && fabs((*lep_eta)[lep_Hindex[i]]) > 1.6) lep_pterr_correction[i] = 4;
                       if (fabs((*lep_eta)[lep_Hindex[i]]) > 1.6 && fabs((*lep_eta)[lep_Hindex[i]]) > 2) lep_pterr_correction[i] = 2.55443;
                       if (fabs((*lep_eta)[lep_Hindex[i]]) > 2 && fabs((*lep_eta)[lep_Hindex[i]]) > 2.5) lep_pterr_correction[i] = 1.95906;

                       }
*/

                    }
                }
            }
            else {
                for (int i=0; i<4; i++) lep_pterr_correction.push_back(1.0);
            }

            ecalDrivenL1 = (*lep_ecalDriven)[lep_Hindex[0]];
            ecalDrivenL2 = (*lep_ecalDriven)[lep_Hindex[1]];
            ecalDrivenL3 = (*lep_ecalDriven)[lep_Hindex[2]];
            ecalDrivenL4 = (*lep_ecalDriven)[lep_Hindex[3]];


            vector<float> ptErrs;
// FIX
//            if (LepNoFSR1.Pt() > LepNoFSR2.Pt()) {
            ptErrs.push_back(lep_pterr_correction[0]*(*lep_pterr)[lep_Hindex[0]]); ptErrs.push_back(lep_pterr_correction[1]*(*lep_pterr)[lep_Hindex[1]]);
//            } else {
//                   ptErrs.push_back(lep_pterr_correction[1]*(*lep_pterr)[lep_Hindex[1]]); ptErrs.push_back(lep_pterr_correction[0]*(*lep_pterr)[lep_Hindex[0]]);
//                   }
            ptErrs.push_back(lep_pterr_correction[2]*(*lep_pterr)[lep_Hindex[2]]); ptErrs.push_back(lep_pterr_correction[3]*(*lep_pterr)[lep_Hindex[3]]);

            //cout<<"idL1 "<<idL1<<" pTL1 "<<LepNoFSR1.Pt()<<" pTErrL1: "<<lep_pterr_correction[0]*(*lep_pterr)[lep_Hindex[0]]<<endl;
            //cout<<"idL2 "<<idL2<<" pTL2 "<<LepNoFSR2.Pt()<<" pTErrL2: "<<lep_pterr_correction[1]*(*lep_pterr)[lep_Hindex[1]]<<endl;
            //cout<<"idL3 "<<idL3<<" pTL3 "<<LepNoFSR3.Pt()<<" pTErrL3: "<<lep_pterr_correction[2]*(*lep_pterr)[lep_Hindex[2]]<<endl;
            //cout<<"idL4 "<<idL4<<" pTL4 "<<LepNoFSR4.Pt()<<" pTErrL4: "<<lep_pterr_correction[3]*(*lep_pterr)[lep_Hindex[3]]<<endl;
            
            map<unsigned int, TLorentzVector> selectedFsrMap;
            vector<float> phoPtErrs;
            for(unsigned int i = 0; i<4; i++) {
                selectedFsrMap[i] = nullFourVector;
                phoPtErrs.push_back(0.0);
            }
            for(unsigned int i = 0; i<4; i++) {
                int index = lep_Hindex[i];
                for (unsigned int j = 0; j<(*fsrPhotons_pt).size();j++) {
                    if ( (*fsrPhotons_lepindex)[j]==index ) {
                        TLorentzVector pho;
                        pho.SetPtEtaPhiM((*fsrPhotons_pt)[j],(*fsrPhotons_eta)[j],(*fsrPhotons_phi)[j],0.0);
                        selectedFsrMap[i] = pho;
                        phoPtErrs[i] = (*fsrPhotons_pterr)[j];
                        //cout<<"fsrpho pt: "<<pho.Pt()<<" index: "<<index+1<<endl;
                    }
                }
            }

// FIX
/*
            TLorentzVector tmpLV; double tmpPhoErr;
            tmpLV = selectedFsrMap[1]; tmpPhoErr = phoPtErrs[1];
            selectedFsrMap[1] = selectedFsrMap[0]; phoPtErrs[1] = phoPtErrs[0];
            selectedFsrMap[0] = tmpLV; phoPtErrs[0] = tmpPhoErr;
*/            
            //cout<<" MZ1: "<<(Lep1+Lep2).M()<<" MZ2: "<<(Lep3+Lep4).M()<<endl;

//if (pTErrL1/pTL1 < 0.025 && pTErrL2/pTL2 < 0.025){// && passedFullSelection && higgs_undec.M() > 105 && higgs_undec.M() < 140 ) {

          kinZfitter->Setup(P4sNoFSR, tmpIDs, ptErrs, selectedFsrMap, phoPtErrs);

            kinZfitter->KinRefitZ();
            mass4lREFIT = (float)kinZfitter->GetRefitM4l();
            mass4lErrREFIT = (float)kinZfitter->GetRefitM4lErrFullCov();
            mass4lErr = (float)kinZfitter->GetM4lErr();
            massZ1REFIT = (float)kinZfitter->GetRefitMZ1(); 
            massZ2REFIT = (float)kinZfitter->GetRefitMZ2(); 
            correlation = (float)kinZfitter->GetCorrelation();



            vector<double> refittedLepErrors = kinZfitter->GetRefitLepErrors();
            pTErrREFITL1 = refittedLepErrors[0];
            pTErrREFITL2 = refittedLepErrors[1];
            pTErrREFITL3 = refittedLepErrors[2];
            pTErrREFITL4 = refittedLepErrors[3];

            vector<TLorentzVector> refittedLepPts = kinZfitter->GetRefitLepPts();
            pTREFITL1 = refittedLepPts[0].Pt();
            pTREFITL2 = refittedLepPts[1].Pt();
            pTREFITL3 = refittedLepPts[2].Pt();
            pTREFITL4 = refittedLepPts[3].Pt();




//}
            
            massZ1 = (Lep1+Lep2).M(); massZ2 = (Lep3+Lep4).M(); 
            mass4l = higgs_undec.M(); pT4l = higgs_undec.Pt();


/*
            if (pTErrL1/pTL1 > 0.03 || pTErrL2/pTL2 > 0.03) {

               mass4lErrREFIT = mass4lErr;
               mass4lREFIT = mass4l;

            }


          
            if (finalState == 1 &&  (pTErrREFITL1 > 0.5 || pTErrREFITL2 > 0.5) ){//&& massZ1REFIT > 91 && massZ1REFIT < 92) {

               cout << "mass4lErrREFIT: " << mass4lErrREFIT << ", mass4lREFIT: " << mass4lREFIT << ", mass4lErrREFIT/mass4lREFIT: " << mass4lErrREFIT/mass4lREFIT << endl;
               cout << "mass4lErr: " << mass4lErr << ", mass4l: " << mass4l << ", mass4lErr/mass4l: " << mass4lErr/mass4l << endl;
               cout << "massZ1: " << massZ1 << ", massZ1REFIT: " << massZ1REFIT << ", Event: " << Event << endl;
               cout << "pTL1: " << pTL1 << ", pTREFITL1: " << pTREFITL1 << endl;
               cout << "pTL2: " << pTL2 << ", pTREFITL2: " << pTREFITL2 << endl;
               cout << "correlation: " << correlation << endl;

               }
*/


            if (abs(idL1)==11 && abs(idL3)==11) {mass4e=mass4l; mass4mu=-1.0; mass2e2mu=-1.0;}
            else if (abs(idL3)==13 && abs(idL3)==13) {mass4e=-1.0; mass4mu=mass4l; mass2e2mu=-1.0; }
            else if (abs(idL1)!=abs(idL3)) {mass4e=-1.0; mass4mu=-1.0; mass2e2mu=mass4l;}

            //cout<<"mass4l "<<mass4l<<" mass4lErr "<<mass4lErr<<" mass4lREFIT "<<mass4lREFIT<<" mass4lErrREFIT "<<mass4lErrREFIT<<endl;
                                
            /*
            if (njets_pt30_eta4p7 > 0) {
                Jet1.SetPtEtaPhiM((*jet_pt)[(*jet_iscleanH4l)[0]],(*jet_eta)[(*jet_iscleanH4l)[0]],(*jet_phi)[(*jet_iscleanH4l)[0]],(*jet_mass)[(*jet_iscleanH4l)[0]]);
            }
            if (njets_pt30_eta4p7 > 1) {
                Jet2.SetPtEtaPhiM((*jet_pt)[(*jet_iscleanH4l)[1]],(*jet_eta)[(*jet_iscleanH4l)[1]],(*jet_phi)[(*jet_iscleanH4l)[1]],(*jet_mass)[(*jet_iscleanH4l)[1]]);
            }

            vector<TLorentzVector> partPprod; vector<int> partIdprod;
            partPprod.push_back(Lep1); partPprod.push_back(Lep2);
            partPprod.push_back(Lep3); partPprod.push_back(Lep4);
            partPprod.push_back(njets_pt30_eta4p7 > 0 ? Jet1 : nullFourVector);
            partPprod.push_back(njets_pt30_eta4p7 > 1 ? Jet2 : nullFourVector);
            
            partIdprod.push_back(idL1); partIdprod.push_back(idL2);
            partIdprod.push_back(idL3); partIdprod.push_back(idL4);
            partIdprod.push_back(0); partIdprod.push_back(0);
            
            if (debug) cout<<"higgs_undec.M() "<<higgs_undec.M()<<endl;

            double p0plus_m4l, bkg_m4l, me_0plus_JHU, me_qqZZ_MCFM, p0minus_VAJHU;
            combinedMEM->computePm4l(P4s,tmpIDs,MEMNames::kNone, p0plus_m4l, bkg_m4l);
            combinedMEM->computeME(MEMNames::kSMHiggs, MEMNames::kJHUGen, P4s, tmpIDs, me_0plus_JHU); // higgs, vector algebra, JHUgen
            combinedMEM->computeME(MEMNames::kqqZZ, MEMNames::kMCFM, P4s, tmpIDs, me_qqZZ_MCFM); // background, vector algebra, MCFM
            combinedMEM->computeME(MEMNames::k0minus, MEMNames::kJHUGen, P4s, tmpIDs, p0minus_VAJHU); // Calculation of PS (0-, fa3=1) gg->H->4l JHUGen ME
            combinedMEM->computeME(MEMNames::kggHZZ_10, MEMNames::kMCFM, P4s, tmpIDs, Dgg10_VAMCFM); // Direct calculation of Dgg (D^kin for off-shell) from MCFM MEs
                        
            D_bkg_kin = me_0plus_JHU / (me_0plus_JHU + me_qqZZ_MCFM);
            D_bkg = me_0plus_JHU * p0plus_m4l / (me_0plus_JHU * p0plus_m4l + me_qqZZ_MCFM * bkg_m4l); // superMELA 
            D_g4 = me_0plus_JHU / ( me_0plus_JHU + p0minus_VAJHU ); // D_0-                
            //mela::computeAngles(P4s[0], tmpIDs[0], P4s[1], tmpIDs[1], P4s[2], tmpIDs[2], P4s[3], tmpIDs[3], cosThetaStar,cosTheta1,cosTheta2,Phi,Phi1);
            */


            //cout<<"idL1: "<<idL1<<" passedFullSelection "<<passedFullSelection<<" passedZ4lSelection: "<<passedZ4lSelection<<endl;
            
            if(debug) cout<<"fill tree"<<endl;
            if(debug) cout<<endl;
            newtree->Fill();
            
        }
        
    }
    cout<<"npass: "<<npass<<endl;
    
}

void SetNewTree(TTree* newtree){

    newtree->Branch("Run",&Run,"Run/l");
    newtree->Branch("Event",&Event,"Event/l");
    newtree->Branch("LumiSect",&LumiSect,"LumiSect/l");
    newtree->Branch("nVtx",&nVtx,"nVtx/I");
    newtree->Branch("passedTrig",&passedTrig,"passedTrig/O");
    newtree->Branch("passedFullSelection",&passedFullSelection,"passedFullSelection/O");
    newtree->Branch("passedFiducialSelection",&passedFiducialSelection,"passedFiducialSelection/O");
    newtree->Branch("passedZ4lSelection",&passedZ4lSelection,"passedZ4lSelection/O");
    newtree->Branch("passedZXCRSelection",&passedZXCRSelection,"passedZXCRSelection/O");
    newtree->Branch("nZXCRFailedLeptons",&nZXCRFailedLeptons,"nZXCRFailedLeptons/I");
    newtree->Branch("finalState",&finalState,"finalState/I");    
    newtree->Branch("GENMH",&GENMH,"GENMH/F");
    newtree->Branch("GENmass4l",&GENmass4l,"GENmass4l/F");
    newtree->Branch("GENmassZ1",&GENmassZ1,"GENmassZ1/F");
    newtree->Branch("GENmassZ2",&GENmassZ2,"GENmassZ2/F");

    newtree->Branch("dataMCWeight",&dataMCWeight,"dataMCWeight/F"); 
    newtree->Branch("k_qqZZ_qcd_M",&k_qqZZ_qcd_M,"k_qqZZ_qcd_M/F");
    newtree->Branch("k_qqZZ_ewk",&k_qqZZ_ewk,"k_qqZZ_ewk/F");
    newtree->Branch("k_ggZZ",&k_ggZZ,"k_ggZZ/F");

    newtree->Branch("pTGENL1",&pTGENL1,"pTGENL1/F");
    newtree->Branch("pTGENL2",&pTGENL2,"pTGENL2/F");
    newtree->Branch("pTGENL3",&pTGENL3,"pTGENL3/F");
    newtree->Branch("pTGENL4",&pTGENL4,"pTGENL4/F");

    newtree->Branch("pTL1",&pTL1,"pTL1/F");
    newtree->Branch("pTL2",&pTL2,"pTL2/F");
    newtree->Branch("pTL3",&pTL3,"pTL3/F");
    newtree->Branch("pTL4",&pTL4,"pTL4/F");
    newtree->Branch("idL1",&idL1,"idL1/I");
    newtree->Branch("idL2",&idL2,"idL2/I");
    newtree->Branch("idL3",&idL3,"idL3/I");
    newtree->Branch("idL4",&idL4,"idL4/I");
    newtree->Branch("etaL1",&etaL1,"etaL1/F");
    newtree->Branch("etaL2",&etaL2,"etaL2/F");
    newtree->Branch("etaL3",&etaL3,"etaL3/F");
    newtree->Branch("etaL4",&etaL4,"etaL4/F");
    newtree->Branch("pTErrL1",&pTErrL1,"pTErrL1/F");
    newtree->Branch("pTErrL2",&pTErrL2,"pTErrL2/F");
    newtree->Branch("pTErrL3",&pTErrL3,"pTErrL3/F");
    newtree->Branch("pTErrL4",&pTErrL4,"pTErrL4/F");
    newtree->Branch("pTErrREFITL1",&pTErrREFITL1,"pTErrREFITL1/F");
    newtree->Branch("pTErrREFITL2",&pTErrREFITL2,"pTErrREFITL2/F");
    newtree->Branch("pTErrREFITL3",&pTErrREFITL3,"pTErrREFITL3/F");
    newtree->Branch("pTErrREFITL4",&pTErrREFITL4,"pTErrREFITL4/F");
    newtree->Branch("pTREFITL1",&pTREFITL1,"pTREFITL1/F");
    newtree->Branch("pTREFITL2",&pTREFITL2,"pTREFITL2/F");
    newtree->Branch("pTREFITL3",&pTREFITL3,"pTREFITL3/F");
    newtree->Branch("pTREFITL4",&pTREFITL4,"pTREFITL4/F");
    newtree->Branch("mL1",&mL1,"mL1/F");
    newtree->Branch("mL2",&mL2,"mL2/F");
    newtree->Branch("mL3",&mL3,"mL3/F");
    newtree->Branch("mL4",&mL4,"mL4/F");
    newtree->Branch("phiL1",&phiL1,"phiL1/F");
    newtree->Branch("phiL2",&phiL2,"phiL2/F");
    newtree->Branch("phiL3",&phiL3,"phiL3/F");
    newtree->Branch("phiL4",&phiL4,"phiL4/F");
    newtree->Branch("correlation",&correlation,"correlation/F");

    newtree->Branch("nFSRPhotons", &nFSRPhotons, "nFSRPhotons/I");

    newtree->Branch("mass4l",&mass4l,"mass4l/F");
    newtree->Branch("mass4l_noFSR",&mass4l_noFSR,"mass4l_noFSR/F");
    newtree->Branch("massZ1_noFSR",&massZ1_noFSR,"massZ1_noFSR/F");
    newtree->Branch("massZ2_noFSR",&massZ2_noFSR,"massZ2_noFSR/F");

    newtree->Branch("mass4l_up",&mass4l_up,"mass4l_up/F");
    newtree->Branch("mass4l_dn",&mass4l_dn,"mass4l_dn/F");
    newtree->Branch("mass4lErr",&mass4lErr,"mass4lErr/F");
    newtree->Branch("mass4lErr_old",&mass4lErr_old,"mass4lErr_old/F");
    newtree->Branch("mass4lREFIT",&mass4lREFIT,"mass4lREFIT/F");
    newtree->Branch("mass4lErrREFIT",&mass4lErrREFIT,"mass4lErrREFIT/F");
    newtree->Branch("massZ1REFIT",&massZ1REFIT,"massZ1REFIT/F");
    newtree->Branch("mass4mu",&mass4mu,"mass4mu/F");
    newtree->Branch("mass4e",&mass4e,"mass4e/F");
    newtree->Branch("mass2e2mu",&mass2e2mu,"mass2e2mu/F");
    newtree->Branch("pT4l",&pT4l,"pT4l/F");
    newtree->Branch("massZ1",&massZ1,"massZ1/F");
    newtree->Branch("massZ1Err",&massZ1Err,"massZ1Err/F");
    newtree->Branch("massZ2",&massZ2,"massZ2/F"); 
    newtree->Branch("njets_pt30_eta4p7",&njets_pt30_eta4p7,"njets_pt30_eta4p7/F");

    newtree->Branch("pTj1",&pTj1,"pTj1/D");
    newtree->Branch("etaj1",&etaj1,"etaj1/D");
    newtree->Branch("pTj2",&pTj2,"pTj2/D");
    newtree->Branch("etaj2",&etaj2,"etaj2/D");

    newtree->Branch("D_bkg_kin", &D_bkg_kin, "D_bkg_kin/F");
    newtree->Branch("D_bkg", &D_bkg, "D_bkg/D");
    newtree->Branch("Dgg10_VAMCFM", &Dgg10_VAMCFM, "Dgg10_VAMCFM/D");
    newtree->Branch("D_g4", &D_g4, "D_g4/D");
    newtree->Branch("Djet_VAJHU", &Djet_VAJHU, "Djet_VAJHU/D");
    newtree->Branch("D_VBF1j_VAJHU",&D_VBF1j_VAJHU,"D_VBF1j_VAJHU/D");
    newtree->Branch("D_WHh_VAJHU",&D_WHh_VAJHU,"D_WHh_VAJHU/D");
    newtree->Branch("D_ZHh_VAJHU",&D_ZHh_VAJHU,"D_ZHh_VAJHU/D");
    newtree->Branch("D_VBF2j",&D_VBF2j,"D_VBF2j/D");
    newtree->Branch("D_VBF1j",&D_VBF1j,"D_VBF1j/D");
    newtree->Branch("D_WHh",&D_WHh,"D_WHh/D");
    newtree->Branch("D_ZHh",&D_ZHh,"D_ZHh/D");

    newtree->Branch("ecalDrivenL1",&ecalDrivenL1,"ecalDrivenL1/I");
    newtree->Branch("ecalDrivenL2",&ecalDrivenL1,"ecalDrivenL2/I");
    newtree->Branch("ecalDrivenL3",&ecalDrivenL1,"ecalDrivenL3/I");
    newtree->Branch("ecalDrivenL4",&ecalDrivenL1,"ecalDrivenL4/I");

    newtree->Branch("EventCat",&EventCat,"EventCat/I");//added by myslef for ZH WH lineshape
    newtree->Branch("eventWeight",&eventWeight,"eventWeight/F");//add by myself for event yeild

}

double GetAllCorr(int id, double pt, double eta) 
{

       double corr = 1;
	if(abs(id)==13){
                  int xbin = LUT_Mu->GetXaxis()->FindBin(pt);
                  int ybin = LUT_Mu->GetYaxis()->FindBin(abs(eta));
                  corr = LUT_Mu->GetBinContent(xbin,ybin);
                  corr = abs(corr);
          }
	if(abs(id)==11){
		int xbin = LUT_Ele->GetXaxis()->FindBin(pt);
		int ybin = LUT_Ele->GetYaxis()->FindBin(abs(eta));
		corr = LUT_Ele->GetBinContent(xbin,ybin);
		corr = abs(corr);
	}


//       double corr_mu_etaRange1[6] = {6.10521e-03, 9.99495e-03, -0.0204909, -0.0138003, -3.06942e-03, -0.0216152};
//       double corr_mu_etaRange2[6] = {-0.0330147, -3.65458e-03, -0.010035, 0.0127915, 0.0348044, 9.42388e-03};
//       double corr_mu_etaRange3[6] = {-0.0273004, -0.0237043, -2.7429e-03, 3.31006e-03, 5.56494e-03, -0.0104418};

//       double corr_mu_etaRange1[12] = {-0.0038, 0.0423, -0.019, -0.0635, -0.05, -0.0344, 0.0079, 0.033, 0.0256, 0.0003, 0.0179, 0.0182};
//       double corr_mu_etaRange2[12] = {-0.0638, 0.0428, -0.0468, -0.0401, -0.0172, -0.0611, -0.0089, 0.0975, 0.0151, -0.0112, 0.086, -0.0339};
//       double corr_mu_etaRange3[12] = {-0.1002, -0.1054, -0.0294, -0.026, 0.0027, -0.0106, -0.0314, 0.0415, 0.0566, 0.0789, -0.0755, -0.0691};

//       double corr_e_etaRange1[6] = {0.0101802, 1.33612e-03, -0.0248805, 8.97247e-03, 0.0261508, -0.0153908};
//       double corr_e_etaRange2[6] = {0.229518, 0.0277794, -0.0136254, 0.01628, 0.0631584, 0.0146108};
//       double corr_e_etaRange3[6] = {-1.32521, 9.26452e-04, -0.0183938, 0.0559622, 0.0785587, 0.133244};

//       double eta1_mu = 0.9;
//       double eta2_mu = 1.4;
//       double eta1_e = 0.9;
//       double eta2_e = 1.4;

//       if (abs(id) == 13 && abs(eta) < eta1_mu ) corr *= GetCorr(id, pt, corr_mu_etaRange1);
//       if (abs(id) == 13 && abs(eta) > eta1_mu && abs(eta) < eta2_mu) corr *= GetCorr(id, pt, corr_mu_etaRange2);
//       if (abs(id) == 13 && abs(eta) > eta2_mu ) corr *= GetCorr(id, pt, corr_mu_etaRange3);

//       if (abs(id) == 11 && abs(eta) < eta1_e ) corr *= GetCorr(id, pt, corr_e_etaRange1);
//       if (abs(id) == 11 && abs(eta) > eta1_e && abs(eta) < eta2_e) corr *= GetCorr(id, pt, corr_e_etaRange2);
//       if (abs(id) == 11 && abs(eta) > eta2_e ) corr *= GetCorr(id, pt, corr_e_etaRange3);

       return corr;
}



double GetCorr(int id, double pt, double corrs[6])
{

  double corr = 1;

  if (abs(id)==13) {


     if (pt > 5 && pt < 20) corr*= abs(corrs[0]);
     if (pt > 20 && pt < 30) corr*= abs(corrs[1]);
     if (pt > 30 && pt < 40) corr*= abs(corrs[2]);
     if (pt > 40 && pt < 50) corr*= abs(corrs[3]);
     if (pt > 50 && pt < 60) corr*= abs(corrs[4]);
     if (pt > 60 && pt < 100) corr*= abs(corrs[5]);
/*
     if (pt > 5 && pt < 15) corr*= abs(corrs[0]);
     if (pt > 15 && pt < 20) corr*= abs(corrs[1]);
     if (pt > 20 && pt < 25) corr*= abs(corrs[2]);
     if (pt > 25 && pt < 30) corr*= abs(corrs[3]);
     if (pt > 30 && pt < 35) corr*= abs(corrs[4]);
     if (pt > 35 && pt < 40) corr*= abs(corrs[5]);
     if (pt > 40 && pt < 45) corr*= abs(corrs[6]);
     if (pt > 45 && pt < 50) corr*= abs(corrs[7]);
     if (pt > 50 && pt < 55) corr*= abs(corrs[8]);
     if (pt > 55 && pt < 60) corr*= abs(corrs[9]);
     if (pt > 60 && pt < 65) corr*= abs(corrs[10]);
     if (pt > 65 && pt < 100) corr*= abs(corrs[11]);
*/
//     corr = 0;

     }


  if (abs(id)==11) {

     if (/*pt > 5 &&*/ pt < 20) corr*= abs(corrs[0]);
     if (pt > 20 && pt < 30) corr*= abs(corrs[1]);
     if (pt > 30 && pt < 40) corr*= abs(corrs[2]);
     if (pt > 40 && pt < 50) corr*= abs(corrs[3]);
     if (pt > 50 && pt < 60) corr*= abs(corrs[4]);
     if (pt > 60 && pt < 100) corr*= abs(corrs[5]);

//     corr = 0;

     }

  corr = corr/100.0;

  return corr;

}
