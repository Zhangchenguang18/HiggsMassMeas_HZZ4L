TString year = "2016";

//set variable
TString variable="relative_mass4lErr";

//set channel
TString fs = "4mu";
Int_t fs1_, fs2_;

//setup refit reco
TString corr="RECO";

void GetHist(TTree* t,TString isData,TString corr, TString process);
void SetAddress(TTree* t,TString isData);

bool passedFullSelection;
Int_t finalState;
Float_t mass4l,mass4lErr,mass4lErr_old,eventWeight,dataMCWeight, mass4lErrREFIT, mass4lREFIT;
Float_t pTL1,pTL2,pTL3,pTL4,etaL1,etaL2,etaL3,etaL4,phiL1,phiL2,phiL3,phiL4,mL1,mL2,mL3,mL4,pTErrL1,pTErrL2,pTErrL3,pTErrL4;
Int_t idL1,idL2,idL3,idL4,ecalDrivenL1,ecalDrivenL2,ecalDrivenL3,ecalDrivenL4;

void mass4lErrTemp(){
	
	if(fs=="4mu"){fs1_=1;fs2_=1;}
	if(fs=="4e"){fs1_=2;fs2_=2;}
	if(fs=="2e2mu"||fs=="2mu2e"){fs1_=3;fs2_=4;}
		
	TFile* f_qqzz = new TFile("/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/slimmed_"+year+"qqZZ.root");
	TTree* t_qqzz = (TTree*)f_qqzz->Get("passedEvents");
	SetAddress(t_qqzz,"MC");
	GetHist(t_qqzz,"MC",corr,"qqZZ");

	TFile* f_ggzz = new TFile("/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/"+year+"ggZZ.root");
	TTree* t_ggzz = (TTree*)f_ggzz->Get("passedEvents");
	SetAddress(t_ggzz,"MC");
	GetHist(t_ggzz,"MC",corr,"ggZZ");

	TFile* f_sig = new TFile("/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/"+year+"GGH_125.root");
	TTree* t_sig = (TTree*)f_sig->Get("passedEvents");
	SetAddress(t_sig,"MC");
	GetHist(t_sig,"MC",corr,"signal");

//	TFile* f_zx = new TFile();
//	TTree* t_zx = (TTree*)f_zx->Get("passedEvents");
//	SetAddress(t_zx,"MC");
//	TH1F* h_zx = GetHist(t_zx,"MC",corr,"pdfErrZX");
	
	
}


void GetHist(TTree* t,TString isData, TString corr, TString process){
	
	TH1F* h = new TH1F("h_Dm","h_Dm",200,0,0.1);
	h->SetTitle("");
	h->Sumw2();
	for(Int_t i=0; i<t->GetEntries(); i++){
		t->GetEntry(i);
		if(passedFullSelection==1&&mass4l>105&&mass4l<140&&(finalState==fs1_||finalState==fs2_)){
	
				if(corr=="RECO"){//miss the uncorr mass4lErr when slimming, re-calculate them
					
					if(isData=="MC"){
						h->Fill(mass4lErr/mass4l,dataMCWeight);
					}//cout<<mass4lErr<<"  "<<mass4lErr_corr<<endl;
					if(isData=="Data")h->Fill(mass4lErr/mass4l);
					
				}		
				if(corr=="REFIT"){
					if(isData=="MC"){
						h->Fill(mass4lErrREFIT/mass4lREFIT,dataMCWeight);
					}
					if(isData=="Data")h->Fill(mass4lErrREFIT/mass4lREFIT);
				}
		}
	}
	h->Smooth();
	h->Smooth();
	h->Scale(1/h->Integral());
	
	TString fs_name;
        if(fs1_==1)fs_name="4mu";
        if(fs1_==2)fs_name="4e";
        if(fs1_>2)fs_name="2e2mu";

	TCanvas c("c","c",1400,1000);
	c.cd();
	h->SetStats(0);
	h->Draw("hist");
	c.SaveAs("/home/chenguan/public_html/EBETemplate/"+year+"/"+corr+"/Dm_"+process+"_"+fs_name+".png");
	
	TFile* name = new TFile("/home/chenguan/public_html/EBETemplate/"+year+"/"+corr+"/Dm_"+process+"_"+fs_name+".root","RECREATE");
        name->cd();
        h->Write();
        name->Close();

}


void SetAddress(TTree* t,TString isData){
	
	t->SetBranchAddress("mass4l",&mass4l);
	t->SetBranchAddress("mass4lErr",&mass4lErr);
	t->SetBranchAddress("mass4lErr_old",&mass4lErr_old);
	t->SetBranchAddress("passedFullSelection",&passedFullSelection);
	t->SetBranchAddress("finalState",&finalState);
	t->SetBranchAddress("pTL1",&pTL1);
	t->SetBranchAddress("pTL2",&pTL2);
	t->SetBranchAddress("pTL3",&pTL3);
	t->SetBranchAddress("pTL4",&pTL4);
	t->SetBranchAddress("etaL1",&etaL1);
	t->SetBranchAddress("etaL2",&etaL2);
	t->SetBranchAddress("etaL3",&etaL3);
	t->SetBranchAddress("etaL4",&etaL4);
	t->SetBranchAddress("phiL1",&phiL1);
	t->SetBranchAddress("phiL2",&phiL2);
	t->SetBranchAddress("phiL3",&phiL3);
	t->SetBranchAddress("phiL4",&phiL4);
	t->SetBranchAddress("mL1",&mL1);
	t->SetBranchAddress("mL2",&mL2);
	t->SetBranchAddress("mL3",&mL3);
	t->SetBranchAddress("mL4",&mL4);
	t->SetBranchAddress("pTErrL1",&pTErrL1);
	t->SetBranchAddress("pTErrL2",&pTErrL2);
	t->SetBranchAddress("pTErrL3",&pTErrL3);
	t->SetBranchAddress("pTErrL4",&pTErrL4);
	t->SetBranchAddress("ecalDrivenL1",&ecalDrivenL1);
	t->SetBranchAddress("ecalDrivenL2",&ecalDrivenL2);
	t->SetBranchAddress("ecalDrivenL3",&ecalDrivenL3);
	t->SetBranchAddress("ecalDrivenL4",&ecalDrivenL4);
	t->SetBranchAddress("eventWeight",&eventWeight);
	t->SetBranchAddress("dataMCWeight",&dataMCWeight);
	t->SetBranchAddress("mass4lErrREFIT",&mass4lErrREFIT);
	t->SetBranchAddress("mass4lREFIT",&mass4lREFIT);

}

