

TString year="2016";


void doPlots(Double_t p1, Double_t p2, bool isMu, bool isECAL, bool ispT, bool isData);

Double_t up, dn;
TString name, sub_name1, sub_name2;
up = 1.0;
dn = 0.0;
sub_name1 = to_string(dn);
sub_name2 = to_string(up);
name = sub_name1+sub_name2;
       
	
void Super_do2DMap(){
//         doPlots(-2.4,2.4,1,0,0,1);//the last three para are "isMu isECAL ispT and isData" respectively
// 	 doPlots(-2.4,2.4,1,0,0,0);
		 doPlots(-2.4,2.4,0,1,0,0);
		 doPlots(-2.4,2.4,0,1,0,1);
/*		 doPlots(-2.4,2.4,0,0,0);
		 doPlots(0,0.9,1,0,1);
		 doPlots(0.9,1.8,1,0,1);
		 doPlots(1.8,2.4,1,0,1);
		 doPlots(0,0.8,0,1,1);
		 doPlots(0.8,1,0,1,1);
		 doPlots(1,1.44,0,1,1);
		 doPlots(1.44,1.57,0,1,1);
		 doPlots(1.57,2,0,1,1);
		 doPlots(2,2.5,0,1,1);
*/
		 }





void doPlots(Double_t p1, Double_t p2, bool isMu, bool isECAL, bool ispT, bool isData){
//	Int_t counter=0;

	Double_t pT_low;
	if(isMu==1)pT_low=5;
	if(isMu==0)pT_low=7;
	
	char p1c[10];
	char p2c[10];
	sprintf(p1c,"%1.2f",p1);
	sprintf(p2c,"%1.2f",p2);
	TString p1s;
	TString p2s;
	p1s=p1c;
	p2s=p2c;
	TString tag;
	TString isdata;
	if(isData==1)isdata="Data";
	if(isData==0)isdata="Sim";
	if(isMu==1&&isECAL==0)tag="_Muon_";
	if(isMu==0&&isECAL==1)tag="_Electron_ECAL_";
	if(isMu==0&&isECAL==0)tag="_Electron_TRACK_";
	
	TString filename = year+"_"+isdata+tag+p1s+"_"+p2s;

	 TH2D *h=new TH2D();
         if(ispT==0)h=new TH2D("h", "", 200, -2.4, 2.4, 200, 0, 0.1);
         if(ispT==1)h=new TH2D("h", "", 200, 0, 100, 200, 0, 0.1);


	TFile *f=new TFile();
	
	
	if(year=="2016"&&isMu==1&&isData==0)f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2016_m2mu.root");
	if(year=="2017"&&isMu==1&&isData==0)f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2017_m2mu.root");
	if(year=="2018"&&isMu==1&&isData==0)f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2018_m2mu.root");
	if(year=="2016"&&isMu==0&&isData==0)f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2016_m2e.root");
	if(year=="2017"&&isMu==0&&isData==0)f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2017_m2e.root");
	if(year=="2018"&&isMu==0&&isData==0)f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2018_m2e.root");
	if(year=="2016"&&isMu==0&&isData==1)f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2016/Electron/SingleElectron_m2e.root");
	if(year=="2016"&&isMu==1&&isData==1)f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2016/Muon/SingleMuon_m2mu.root");
	//if(year=="hualin"&&isMu==0)f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/DYJetsToLL_M-50_kalman_v4_m2e.root");	
	
	TTree *t=(TTree*)f->Get("passedEvents");
	
	Double_t genLep_pt1, genLep_pt2;
	Double_t pterr1,pterr2;
	Int_t id1,id2;
	Double_t pT1,pT2;
	Double_t eta1,eta2;	
	Double_t weight;
	Double_t massZ;
	Double_t massZErr;
	Int_t lep1_ecalDriven,lep2_ecalDriven;
	Double_t RelIso1, RelIso2;
	t->SetBranchAddress("RelIso1",&RelIso1);	
	t->SetBranchAddress("RelIso2",&RelIso2);
	t->SetBranchAddress("genLep_pt1",&genLep_pt1);
	t->SetBranchAddress("genLep_pt2",&genLep_pt2);
	t->SetBranchAddress("pterr1",&pterr1);
	t->SetBranchAddress("pterr2",&pterr2);
//	t->SetBranchAddress("id1",&id1);
//	t->SetBranchAddress("id2",&id2);
	t->SetBranchAddress("pT1",&pT1);
	t->SetBranchAddress("pT2",&pT2);
	t->SetBranchAddress("eta1",&eta1);
	t->SetBranchAddress("eta2",&eta2);
	t->SetBranchAddress("weight",&weight);
	t->SetBranchAddress("massZ",&massZ);
	t->SetBranchAddress("massZErr",&massZErr);
	t->SetBranchAddress("lep1_ecalDriven",&lep1_ecalDriven);
	t->SetBranchAddress("lep2_ecalDriven",&lep2_ecalDriven);
	

	Int_t sum = (Int_t)t->GetEntries();
	
	for(Int_t i=0; i<sum; i++){

		t->GetEntry(i);
		
		if(RelIso1>0.35||RelIso2>0.35)continue;

		if(massZErr<7.2&&massZErr>0.2&&pT1<100&&pT1>pT_low&&eta1<p2&&eta1>p1&&lep1_ecalDriven==isECAL&&massZErr/massZ<up&&dn<massZErr/massZ){
			if(ispT==0)h->Fill(eta1,pterr1/pT1);
			if(ispT==1)h->Fill(pT1,pterr1/pT1);
//			counter=counter+1;
//		
	}
		if(massZErr>0.2&&massZErr<7.2&&pT2<100&&pT2>pT_low&&eta2<p2&&eta2>p1&&lep2_ecalDriven==isECAL&&massZErr/massZ<up&&dn<massZErr/massZ){
			if(ispT==0)h->Fill(eta2,pterr2/pT2);
			if(ispT==1)h->Fill(pT2,pterr2/pT2);
//			counter=counter+1;
		}

	}
	
	if(ispT==0)h->GetXaxis()->SetTitle("eta");
	if(ispT==1)h->GetXaxis()->SetTitle("pT");
        h->GetYaxis()->SetTitle("PtErr/Pt");
	h->Scale(1/h->Integral());
	h->SetStats(kFALSE);
		
	TCanvas *c=new TCanvas("c","c",1000,800);
	c->cd();
	h->Draw("COLZ");
	c->SetLeftMargin(0.15);
	c->SetRightMargin(0.15);
	c->SaveAs("/home/chenguan/public_html/2Dmaps/"+year+"/"+filename+"_"+name+".png");
	
	h->Print();
	//cout<<sum<<endl;
//	cout<<counter<<endl;	
//	cout<<p1<<"****"<<p2<<endl;
 delete c; delete h; delete t; delete f; 


}

