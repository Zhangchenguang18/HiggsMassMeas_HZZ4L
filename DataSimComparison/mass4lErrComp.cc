//setup mass region
TString massregion = "high";
Double_t massregion_low, massregion_high;

//
Double_t Lumi;

//setup hist
Int_t bins=30;
Double_t low=0;
Double_t high=0.06;
Double_t ptlow;

TString year = "2016";

//set variable
TString variable="relative_mass4lErr";

//set channel
TString fs = "4mu";
Int_t fs1_, fs2_;

//setup corr or uncorr
TString corr="Corr";

//
Double_t qqzz_=0;
Double_t ggzz_4e=0;
Double_t ggzz_4mu=0;
Double_t ggzz_2e2mu=0;

TH1F* GetHist(TTree* t,TString isData,TString corr, TString process);
void SetAddress(TTree* t,TString isData);
void Compare(TH1F* data, TH1F* qqzz, TH1F* ggzz);
TH1F* GetRatio(TH1F* data, TH1F* qqzz, TH1F* ggzz);
void mass4lErrComp_(TString year, TString fs, TString corr, TString massregion);

bool passedFullSelection;
Int_t finalState;
Float_t mass4l,mass4lErr,mass4lREFIT,mass4lErrREFIT,mass4lErr_old,eventWeight,dataMCWeight,k_ggZZ,k_qqZZ_ewk,k_qqZZ_qcd_M,k;
Float_t pTL1,pTL2,pTL3,pTL4,etaL1,etaL2,etaL3,etaL4,phiL1,phiL2,phiL3,phiL4,mL1,mL2,mL3,mL4,pTErrL1,pTErrL2,pTErrL3,pTErrL4;
Int_t idL1,idL2,idL3,idL4,ecalDrivenL1,ecalDrivenL2,ecalDrivenL3,ecalDrivenL4;

void mass4lErrComp(){
	for(Int_t i=0; i<3; i++){
		if(i==0)year="2016";
		if(i==1)year="2017";
		if(i==2)year="2018";
		for(Int_t j=0; j<3; j++){
			if(j==0)fs="4mu";
			if(j==1)fs="4e";
			if(j==2)fs="2e2mu";
			for(Int_t k=0; k<2; k++){
				if(k==0)corr="Corr";
				if(k==1)corr="Uncorr";
				for(Int_t m=0; m<2; m++){
					if(m==0)massregion="high";
					if(m==1)massregion="low";
					mass4lErrComp_(year, fs, corr, massregion);
				}
			}
		}
	}
}


void mass4lErrComp_(TString year, TString fs, TString corr, TString massregion){
	
	if(massregion=="low"){massregion_low=80;massregion_high=100;}
	if(massregion=="high"){massregion_low=180;massregion_high=200;}
	if(massregion=="full"){massregion_low=70;massregion_high=999999;}

	if(fs=="4mu"){fs1_=1;fs2_=1;}
	if(fs=="4e"){fs1_=2;fs2_=2;}
	if(fs=="2e2mu"||fs=="2mu2e"){fs1_=3;fs2_=4;}
		
	if(year=="2016")Lumi=35.92;
	if(year=="2017")Lumi=41.53;
	if(year=="2018")Lumi=59.68;

	TFile* f_qqzz = new TFile();
	f_qqzz = new TFile("/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/slimmed_"+year+"qqZZ.root");
	TTree* t_qqzz = (TTree*)f_qqzz->Get("passedEvents");
	SetAddress(t_qqzz,"MC");
	TH1F* h_qqzz = GetHist(t_qqzz,"MC",corr,"qqzz");

	TFile* f_ggzz = new TFile();
	if(fs=="2e2mu")f_ggzz = new TFile("/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/"+year+"ggZZ_2e2mu.root");
	if(fs=="4e")f_ggzz = new TFile("/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/"+year+"ggZZ_4e.root");
	if(fs=="4mu")f_ggzz = new TFile("/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/"+year+"ggZZ_4mu.root");
	TTree* t_ggzz = (TTree*)f_ggzz->Get("passedEvents");
	SetAddress(t_ggzz,"MC");
	TH1F* h_ggzz = GetHist(t_ggzz,"MC",corr,"ggzz_"+fs);


	TFile* f_data = new TFile();
	f_data=new TFile("/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/"+year+"_noDuplicates.root");
	TTree* t_data = (TTree*)f_data->Get("passedEvents");
	SetAddress(t_data,"Data");
	TH1F* h_data = GetHist(t_data,"Data",corr,"0");
	
	Compare(h_data,h_qqzz,h_ggzz);
	
}


TH1F* GetHist(TTree* t,TString isData, TString corr, TString process){
	
	TH1F* h = new TH1F("h","h",bins,low,high);
	for(Int_t i=0; i<t->GetEntries(); i++){
		t->GetEntry(i);
		if(passedFullSelection==1&&mass4l>massregion_low&&mass4l<massregion_high&&(finalState==fs1_||finalState==fs2_)){
				if(process=="qqzz")k=k_qqZZ_ewk*k_qqZZ_qcd_M;
				if(process=="ggzz_4mu"||process=="ggzz_4e"||process=="gzz_2e2mu")k=k_ggZZ;
				if(corr=="Uncorr"){//miss the uncorr mass4lErr when slimming, re-calculate them
					
					if(isData=="MC"){
						h->Fill(mass4lErr/mass4l,eventWeight);
						if(process=="qqzz")qqzz_ = qqzz_ + eventWeight*k;
						if(process=="ggzz_4e")ggzz_4e = ggzz_4e + eventWeight*k;
						if(process=="ggzz_4mu")ggzz_4mu = ggzz_4mu + eventWeight*k;
						if(process=="ggzz_2e2mu")ggzz_2e2mu = ggzz_2e2mu + eventWeight*k;
					}//cout<<mass4lErr<<"  "<<mass4lErr_corr<<endl;
					if(isData=="Data")h->Fill(mass4lErr/mass4l);
					
				}		
				if(corr=="Corr"){
					if(isData=="MC"){
						h->Fill(mass4lErrREFIT/mass4lREFIT,eventWeight);
						if(process=="qqzz")qqzz_ = qqzz_ + eventWeight*k;
					        if(process=="ggzz_4e")ggzz_4e = ggzz_4e + eventWeight*k;
					        if(process=="ggzz_4mu")ggzz_4mu = ggzz_4mu + eventWeight*k;       
						if(process=="ggzz_2e2mu")ggzz_2e2mu = ggzz_2e2mu + eventWeight*k;
					}
					if(isData=="Data")h->Fill(mass4lErrREFIT/mass4lREFIT);
				}
			}
		}
		h->Sumw2();
		return h;

}


TH1F* GetRatio(TH1F* data, TH1F* qqzz, TH1F* ggzz){

	qqzz->Add(ggzz);
	qqzz->Add(data,-1);
	qqzz->Divide(data);
	return qqzz;

}

void Normalize(TH1F* qqzz, TH1F* ggzz){
	
	Double_t N_qqzz, N_ggzz;
	
//	TFile* f_qqzz = TFile::Open("root://cmsio5.rc.ufl.edu//store/user/t2/users/ferrico/Full_RunII/qqZZ/ZZTo4L_powheg_"+year+".root");
//	TTree* t_qqzz = (TTree*)f_qqzz->Get("Ana/passedEvents");// dont know if Filippo finishe all crab jobs, have to produce NTuple by myself
//	N_qqzz = t_qqzz->GetEntries();
	
//	TFile* f_ggzz;
//	if(fs=="4mu")f_ggzz = TFile::Open("root://cmsio5.rc.ufl.edu//store/user/t2/users/ferrico/Full_RunII/ggZZ/GluGluToContinToZZTo4mu_"+year+".root");
//	if(fs=="4e")f_ggzz = TFile::Open("root://cmsio5.rc.ufl.edu//store/user/t2/users/ferrico/Full_RunII/ggZZ/GluGluToContinToZZTo4e_"+year+".root");
//	if(fs=="2e2mu")f_ggzz = TFile::Open("root://cmsio5.rc.ufl.edu//store/user/t2/users/ferrico/Full_RunII/ggZZ/GluGluToContinToZZTo2e2mu_"+year+".root");
//	TTree* t_ggzz = (TTree*)f_ggzz->Get("Ana/passedEvents");
//	N_ggzz = t_ggzz->GetEntries();
//cout<<N_ggzz<<"  "<<N_qqzz<<endl;
	
	if(fs=="4e"){
                 if(year=="2016")N_ggzz=965000;
                 if(year=="2017")N_ggzz=976928;
                 if(year=="2018")N_ggzz=913200;
        }
        if(fs=="4mu"){
                 if(year=="2016")N_ggzz=995200;
                 if(year=="2017")N_ggzz=909000;
                 if(year=="2018")N_ggzz=982548;
        }
        if(fs=="2e2mu"){
                 if(year=="2016")N_ggzz=1469600;
                 if(year=="2017")N_ggzz=500000;
                 if(year=="2018")N_ggzz=494000;
        }
        if(year=="2016")N_qqzz=77998181+18470796;
        if(year=="2017")N_qqzz=50756359+36007127;
        if(year=="2018")N_qqzz=19089600;

cout<<"qqzz_"+fs<<": "<<qqzz_*(Lumi*1000/N_qqzz)<<endl;
cout<<"ggzz_4e: "<<ggzz_4e*(Lumi*1000/N_ggzz)<<endl;
cout<<"ggzz_4mu: "<<ggzz_4mu*(Lumi*1000/N_ggzz)<<endl;
cout<<"ggzz_2e2mu: "<<ggzz_2e2mu*(Lumi*1000/N_ggzz)<<endl;

	qqzz->Scale(Lumi*1000/N_qqzz);
	ggzz->Scale(Lumi*1000/N_ggzz);
//	delete f_qqzz; delete f_ggzz;

}


void Compare(TH1F* data, TH1F* qqzz, TH1F* ggzz){
	
	Normalize(qqzz, ggzz);
	
	TH1F* data_ = (TH1F*)data->Clone();
	TH1F* qqzz_ = (TH1F*)qqzz->Clone();	
	TH1F* ggzz_ = (TH1F*)ggzz->Clone();
	
	TH1F* ratio = GetRatio(data_,qqzz_,ggzz_);	
	
	Float_t binwidth = (high-low)/bins;
	TString binwidth_s = to_string(binwidth).substr(0,5);
	
	THStack* hs = new THStack("hs","");
	ggzz->SetFillColor(7);
	ggzz->SetLineColor(7);
	hs->Add(ggzz);
	qqzz->SetFillColor(3);
	qqzz->SetLineColor(3);
	hs->Add(qqzz);
	
	TH1F* h = new TH1F("h","",bins,low,high);
	h->SetStats(0);
	hs->SetHistogram(h);
	hs->GetHistogram()->GetXaxis()->SetLabelSize(0);
	hs->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
	hs->GetHistogram()->GetYaxis()->SetTitle("Events/ "+binwidth_s+"GeV");
	
	TCanvas c("c","c",1000,1000);
	c.SetTopMargin(0);
	c.SetBottomMargin(0);
	c.cd();
	TPad c_1("c_1","c_1",0,0.3,1,0.95);
	c_1.Draw();
	c_1.cd();
	c_1.SetTopMargin(0.08);
        c_1.SetBottomMargin(0.02);
	
	hs->Draw("hist");

	data->Draw("same");
	
	Double_t max = qqzz->GetMaximum();
	if(max<data->GetMaximum())max = data->GetMaximum();
	hs->SetMaximum(1.2*max);

        data->SetMarkerColor(kBlue);
        data->SetMarkerStyle(34);
	
	TLegend* legend = new TLegend(0.7,0.7,0.8,0.8);
        legend->AddEntry(data,"Data","P");
	legend->AddEntry(qqzz,"qqZZ","F");
	legend->AddEntry(ggzz,"ggZZ","F");
        legend->SetTextSize(0.04);
        legend->SetLineWidth(0);
        legend->SetFillColor(0);
        legend->SetBorderSize(0);
	legend->Draw("same");
	
	c.cd();
	TPad c_2("c_2","c_2",0.0,0.1,1,0.3);
	c_2.Draw();
	c_2.cd();
	c_2.SetTopMargin(0.05);
	c_2.SetBottomMargin(0.3);
	
	ratio->SetTitle("");
	ratio->SetStats(0);
	ratio->SetMarkerColor(kBlue);
	ratio->SetMarkerStyle(34);
	ratio->SetMaximum(1);
	ratio->SetMinimum(-1);
	ratio->GetYaxis()->SetTitle("(MC-Data)/Data");
	ratio->GetXaxis()->SetTitle("mass4lErr/mass4l");

        ratio->GetXaxis()->SetTitleFont(43);
        ratio->GetXaxis()->SetLabelFont(43);
        ratio->GetXaxis()->SetLabelSize(25);
        ratio->GetXaxis()->SetTitleSize(25);
 
        ratio->GetYaxis()->SetTitleFont(43);
        ratio->GetYaxis()->SetLabelFont(43);
        ratio->GetYaxis()->SetLabelSize(16);
        ratio->GetYaxis()->SetTitleSize(25);

        ratio->GetXaxis()->SetTitleOffset(5.);
        ratio->GetYaxis()->SetTitleOffset(1.5);

	TLine* line = new TLine(low,0,high,0);
        line->SetLineColor(2);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(2);
	TLine* line1 = new TLine(low,0.2,high,0.2);
	line1->SetLineColor(2);
	line1->SetLineStyle(kDashed);
	line1->SetLineWidth(2);
	TLine* line2 = new TLine(low,-0.2,high,-0.2);
	line2->SetLineColor(2);
	line2->SetLineStyle(kDashed);
	line2->SetLineWidth(2);
	ratio->Draw("p");
	line->Draw();
	line1->Draw("same");
	line2->Draw("same");
	
	TString region;
	if(massregion=="low")region="_80_100";
	if(massregion=="high")region="_180_200";
	if(massregion=="full")region="_70_99999";
	c.SaveAs("/home/chenguan/public_html/DataVSMCin4LEvents/"+year+"/"+fs+"_"+variable+region+"_"+corr+".png");
	c.SaveAs("/home/chenguan/public_html/DataVSMCin4LEvents/"+year+"/"+fs+"_"+variable+region+"_"+corr+".pdf");

}

void SetAddress(TTree* t,TString isData){
	
	t->SetBranchAddress("mass4l",&mass4l);
	t->SetBranchAddress("mass4lErr",&mass4lErr);
	t->SetBranchAddress("mass4lErrREFIT",&mass4lErrREFIT);
	t->SetBranchAddress("mass4lREFIT",&mass4lREFIT);
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
	t->SetBranchAddress("k_ggZZ",&k_ggZZ);
        t->SetBranchAddress("k_qqZZ_qcd_M",&k_qqZZ_qcd_M);
        t->SetBranchAddress("k_qqZZ_ewk",&k_qqZZ_ewk);

}

