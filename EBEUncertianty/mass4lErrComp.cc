//setup hist
Int_t bins=50;
Double_t low=0;
Double_t high=7;
Double_t ptlow;

TString year = "2018";

//set variable
TString variable="massZErr";

//set channel
TString fs = "2e";

//
TString corr="Corr";

TH1D* GetHist(TTree* t,TString isData,TString corr);
void SetAddress(TTree* t,TString isData);
void Compare(TH1D* data, TH1D* mc);
TH1D* GetRatio(TH1D* data, TH1D* mc);

Double_t weight,massZ, massZErr, m1, m2, phi1, phi2, pt1, pt2, pterr1, pterr2, eta1, eta2, RelIso1, RelIso2, GENmass2l;
Int_t lep1_ecalDriven,lep2_ecalDriven;

void mass2lErrComp(){
	
	if(fs=="2mu")ptlow=5;
	if(fs=="2e")ptlow=7;	

	TFile* f_mc = new TFile();
	if(fs=="2mu"&&year=="2016")f_mc=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2016_m2mu.root");
	if(fs=="2e"&&year=="2016")f_mc=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2016_m2e.root");
	if(fs=="2mu"&&year=="2017")f_mc=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2017_m2mu.root");
	if(fs=="2e"&&year=="2017")f_mc=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2017_m2e.root");
	if(fs=="2mu"&&year=="2018")f_mc=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2018_m2mu.root");
	if(fs=="2e"&&year=="2018")f_mc=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2018_m2e.root");
	TTree* t_mc = (TTree*)f_mc->Get("passedEvents");
	SetAddress(t_mc,"MC");
	TH1D* h_mc = GetHist(t_mc,"MC",corr);

	TFile* f_data = new TFile();
	if(fs=="2mu"&&year=="2016")f_data=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2016/Muon/SingleMuon_m2mu.root");
	if(fs=="2e"&&year=="2016")f_data=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2016/Electron/SingleElectron_m2e.root");
	if(fs=="2mu"&&year=="2017")f_data=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2017/Muon/SingleMuon_m2mu.root");
	if(fs=="2e"&&year=="2017")f_data=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2017/Electron/SingleElectron_m2e.root");
	if(fs=="2mu"&&year=="2018")f_data=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2018/Muon/SingleMuon_m2mu.root");
	if(fs=="2e"&&year=="2018")f_data=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2018/Electron/EGamma_m2e.root");
	TTree* t_data = (TTree*)f_data->Get("passedEvents");
	SetAddress(t_data,"Data");
	TH1D* h_data = GetHist(t_data,"Data",corr);
	
	Compare(h_data,h_mc);
	
}


TH1D* GetHist(TTree* t,TString isData, TString corr){
	
	TFile* f_mu = new TFile("/raid/raid9/chenguan/Mass/leptonpTErrCorrector_2/LUT_finalbinning_backup/"+year+"_"+isData+"_mu.root");
	TFile* f_e1 = new TFile("/raid/raid9/chenguan/Mass/leptonpTErrCorrector_2/LUT_finalbinning_backup/"+year+"_"+isData+"_e1.root");
	TFile* f_e3 = new TFile("/raid/raid9/chenguan/Mass/leptonpTErrCorrector_2/LUT_finalbinning_backup/"+year+"_"+isData+"_e3.root");
	TH2D* h_mu = (TH2D*)f_mu->Get("mu");
	TH2D* h_e1 = (TH2D*)f_e1->Get("e1");
	TH2D* h_e3 = (TH2D*)f_e3->Get("e3");
	
	TH1D* h = new TH1D("h","h",bins,low,high);
	for(Int_t i=0; i<t->GetEntries(); i++){
		t->GetEntry(i);
		if(GENmass2l==0&&isData==0)continue;
		if(massZ<100&&massZ>80&&massZErr>0.2&&massZErr<7.2&&pt1>ptlow&&pt1<100&&pt2>ptlow&&pt2<100&&abs(eta1)<2.4&&abs(eta2)<2.4&&RelIso1<0.35&&RelIso2<0.35){
	//		if(massZErr/massZ>0&&massZErr/massZ<0.007){
			if(corr=="Corr"){
				Double_t lambda1=1;
				Double_t lambda2=1;
				if(fs=="2mu"){
					Int_t xbin1 = h_mu->GetXaxis()->FindBin(pt1);
					Int_t ybin1 = h_mu->GetYaxis()->FindBin(abs(eta1));
					Int_t xbin2 = h_mu->GetXaxis()->FindBin(pt2);
					Int_t ybin2 = h_mu->GetYaxis()->FindBin(abs(eta2));
					lambda1 = h_mu->GetBinContent(xbin1,ybin1);
					lambda2 = h_mu->GetBinContent(xbin2,ybin2);
				}
				if(fs=="2e"){
					if(lep1_ecalDriven==1){
						Int_t xbin1 = h_e1->GetXaxis()->FindBin(abs(eta1));
						Int_t ybin1 = h_e1->GetYaxis()->FindBin(pterr1/pt1);
						lambda1 = h_e1->GetBinContent(xbin1,ybin1);
					}
					if(lep1_ecalDriven==0){
						Int_t xbin1 = h_e3->GetXaxis()->FindBin(pt1);
						Int_t ybin1 = h_e3->GetYaxis()->FindBin(abs(eta1));
						lambda1 = h_e3->GetBinContent(xbin1,ybin1);
					}
					if(lep2_ecalDriven==1){
						Int_t xbin2 = h_e1->GetXaxis()->FindBin(abs(eta2));
						Int_t ybin2 = h_e1->GetYaxis()->FindBin(pterr2/pt2);
						lambda2 = h_e1->GetBinContent(xbin2,ybin2);
					}
					if(lep2_ecalDriven==0){
						Int_t xbin2 = h_e3->GetXaxis()->FindBin(pt2);
						Int_t ybin2 = h_e3->GetYaxis()->FindBin(abs(eta2));
						lambda2 = h_e3->GetBinContent(xbin2,ybin2);
					}
				}
				Double_t pt1err_corr = pterr1*lambda1;
				Double_t pt2err_corr = pterr2*lambda2;
				TLorentzVector lep1, lep2;
				lep1.SetPtEtaPhiM(pt1,eta1,phi1,m1);
				lep2.SetPtEtaPhiM(pt2,eta2,phi2,m2);		
				TLorentzVector lep1p, lep2p;
				lep1p.SetPtEtaPhiM(pt1+pt1err_corr,eta1,phi1,m1);
				lep2p.SetPtEtaPhiM(pt2+pt2err_corr,eta2,phi2,m2);
				Double_t dm1corr = (lep1p+lep2).M()-(lep1+lep2).M();
				Double_t dm2corr = (lep1+lep2p).M()-(lep1+lep2).M();
				Double_t massZErr_corr = TMath::Sqrt(dm1corr*dm1corr+dm2corr*dm2corr);
				h->Fill(massZErr_corr,weight);
			}		
			if(corr=="Uncorr")h->Fill(massZErr,weight);
		}
	}
	h->Sumw2();
	h->Scale(1/h->Integral());
	return h;

}


TH1D* GetRatio(TH1D* data, TH1D* mc){

	mc->Add(data,-1);
	mc->Divide(data);
	return mc;

}


void Compare(TH1D* data, TH1D* mc ){
	
	TH1D* data_ = (TH1D*)data->Clone();
	TH1D* mc_ = (TH1D*)mc->Clone();	
	TH1D* ratio = GetRatio(data_,mc_);	
	
		
	TCanvas c("c","c",1000,1000);
	c.SetTopMargin(0);
	c.SetBottomMargin(0);
	c.cd();
	TPad c_1("c_1","c_1",0,0.3,1,0.95);
	c_1.Draw();
	c_1.cd();
	c_1.SetTopMargin(0.08);
        c_1.SetBottomMargin(0.02);
	
	mc->Draw("hist");	
	mc->SetTitle("");
	mc->GetYaxis()->SetTitleSize(25);
	mc->GetYaxis()->SetTitleFont(43);
	mc->GetYaxis()->SetTitleOffset(2.);
	mc->GetYaxis()->SetTitle("Events");
	mc->GetYaxis()->SetLabelFont(43);
	mc->GetYaxis()->SetLabelSize(25);
	mc->GetXaxis()->SetLabelSize(0);
	mc->SetStats(0);
	data->Draw("same");
	
	Double_t max = mc->GetMaximum();
	if(max<data->GetMaximum())max = data->GetMaximum();
	mc->SetMaximum(1.1*max);

	mc->SetFillColor(kGreen);
	mc->SetLineColor(kGreen);
        data->SetMarkerColor(kBlue);
        data->SetMarkerStyle(34);
	
	TLegend* legend = new TLegend(0.7,0.7,0.8,0.8);
        legend->AddEntry(data,"Data","P");
	legend->AddEntry(mc,"MC","F");
        legend->SetTextSize(0.06);
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
	ratio->SetMaximum(0.5);
	ratio->SetMinimum(-0.5);
	ratio->GetYaxis()->SetTitle("(MC-Data)/Data");
	ratio->GetXaxis()->SetTitle(variable);

        ratio->GetXaxis()->SetTitleFont(43);
        ratio->GetXaxis()->SetLabelFont(43);
        ratio->GetXaxis()->SetLabelSize(25);
        ratio->GetXaxis()->SetTitleSize(25);
 
        ratio->GetYaxis()->SetTitleFont(43);
        ratio->GetYaxis()->SetLabelFont(43);
        ratio->GetYaxis()->SetLabelSize(16);
        ratio->GetYaxis()->SetTitleSize(25);

        ratio->GetXaxis()->SetTitleOffset(5.);
        ratio->GetYaxis()->SetTitleOffset(2.);

	TLine* line = new TLine(low,0,high,0);
        line->SetLineColor(2);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(2);
	
	ratio->Draw("p");
	line->Draw("same");
	c.SaveAs("/home/chenguan/public_html/DataVSMCin2LEvents/"+year+"/"+fs+"_"+variable+"_"+corr+".png");
	c.SaveAs("/home/chenguan/public_html/DataVSMCin2LEvents/"+year+"/"+fs+"_"+variable+"_"+corr+".pdf");

}

void SetAddress(TTree* t,TString isData){
	t->SetBranchAddress("massZ",&massZ);
	t->SetBranchAddress("massZErr",&massZErr);
	t->SetBranchAddress("pT1",&pt1);
	t->SetBranchAddress("pT2",&pt2);
	t->SetBranchAddress("pterr1",&pterr1);
	t->SetBranchAddress("pterr2",&pterr2);
	t->SetBranchAddress("phi1",&phi1);
	t->SetBranchAddress("phi2",&phi2);
	t->SetBranchAddress("m1",&m1);
	t->SetBranchAddress("m2",&m2);
	t->SetBranchAddress("eta1",&eta1);
	t->SetBranchAddress("eta2",&eta2);
	t->SetBranchAddress("RelIso1",&RelIso1);
	t->SetBranchAddress("RelIso2",&RelIso2);
	t->SetBranchAddress("lep1_ecalDriven",&lep1_ecalDriven);
	t->SetBranchAddress("lep2_ecalDriven",&lep2_ecalDriven);
	if(isData=="MC")t->SetBranchAddress("GENmass2l",&GENmass2l);
	t->SetBranchAddress("weight",&weight);

}

