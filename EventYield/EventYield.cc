void EventYield_(TString process);
void SetTreeAddress(TTree* t);
Double_t* ReadTree(TTree* t, Int_t fs1, Int_t fs2, Int_t masspoint, TString process);
void MakeResult(vector<Double_t*> points_4mu_v, vector<Double_t*> points_4e_v, vector<Double_t*> points_2e2mu_v, TString process);
Int_t GetNtot(TString filename);

Float_t mass4l, eventWeight, pTL1, pTL2, pTL3, pTL4, etaL1, etaL2, etaL3, etaL4, k_ggZZ, k_qqZZ_qcd_M, k_qqZZ_ewk;
Int_t finalState, EventCat;
bool passedFullSelection;
Double_t Lumi;

TString year = "2018";


void EventYield(){
	
//	EventYield_("GGH");
//	EventYield_("VBF");
	EventYield_("WplusH");
//	EventYield_("WminusH");
//	EventYield_("ZH");
//	EventYield_("ttH");

}

void EventYield_(TString process){
	
	gStyle->SetTitleOffset(2,"Y");	
	
	TString inputpath="/raid/raid9/chenguan/input/AfterSync_200424/";

	Int_t fs1, fs2;
	
	TString str_fs;
	
	Int_t masspoint[5] = {120, 124, 125, 126, 130};

	Double_t* points = new Double_t[2];
	vector<Double_t*> points_4mu;
	points_4mu.clear();
	vector<Double_t*> points_4e;
	points_4e.clear();
	vector<Double_t*> points_2e2mu;
	points_2e2mu.clear();	

	for(Int_t k=0; k<5; k++){               //masspoint
				
		TString mh = to_string(masspoint[k]);
		TString NTfilename;
		NTfilename = inputpath+year+process+"_"+mh+".root";
		if(process=="ttH"&&k<4)NTfilename = inputpath+year+process+"_"+mh+"_v14_v2.root";
		TFile* f = new TFile(NTfilename);
		TTree* t = (TTree*)f->Get("Ana/passedEvents");
		SetTreeAddress(t);
		
		points = ReadTree(t, 1, 1, masspoint[k], process);
		points_4mu.push_back(points);
		
		points = ReadTree(t, 2, 2, masspoint[k], process);
		points_4e.push_back(points);
		
		points = ReadTree(t, 3, 4, masspoint[k], process);
		points_2e2mu.push_back(points);
		
		delete t;
					
					
	}
	
	MakeResult(points_4mu, points_4e, points_2e2mu, process);
	
	delete[] points;
			
			
		
}

void MakeResult(vector<Double_t*> points_4mu_v, vector<Double_t*> points_4e_v, vector<Double_t*> points_2e2mu_v, TString process){
	
	if(year=="2016")Lumi=35.92;
	if(year=="2017")Lumi=41.53;
	if(year=="2018")Lumi=59.68;

	Double_t* x = new Double_t[5];
	x[0] = 120;
	x[1] = 124;
	x[2] = 125;
	x[3] = 126;
	x[4] = 130;
	
	Double_t* points_4mu = new Double_t[5];
	Double_t* points_4e = new Double_t[5];
	Double_t* points_2e2mu = new Double_t[5];
	
	Double_t* Err_points_4mu = new Double_t[5];
	Double_t* Err_points_4e = new Double_t[5];
	Double_t* Err_points_2e2mu = new Double_t[5];
	

		
		for(Int_t j=0; j<5; j++){
			
			points_4mu[j]=points_4mu_v[j][0];
			Err_points_4mu[j]=points_4mu_v[j][1];
			
			points_4e[j]=points_4e_v[j][0];
			Err_points_4e[j]=points_4e_v[j][1];
			
			points_2e2mu[j]=points_2e2mu_v[j][0];
			Err_points_2e2mu[j]=points_2e2mu_v[j][1];
			
		}


	
	TGraph* g_4mu = new TGraphErrors(5, x, points_4mu, 0, Err_points_4mu);
	TGraph* g_4e = new TGraphErrors(5, x, points_4e, 0, Err_points_4e);
	TGraph* g_2e2mu = new TGraphErrors(5, x, points_2e2mu, 0, Err_points_2e2mu);
	
	TF1* f1 = new TF1("f1","[0]+[1]*x+[2]*x*x",118,130);
	TF1* f2 = new TF1("f2","[0]+[1]*x+[2]*x*x",118,130);
	TF1* f3 = new TF1("f3","[0]+[1]*x+[2]*x*x",118,130);
	
/*	f1->SetParameter(2,0);
	f2->SetParameter(2,0);
	f3->SetParameter(2,0);

	f1->SetParLimits(2,0,0);
	f2->SetParLimits(2,0,0);
	f3->SetParLimits(2,0,0);	
*/	
	TFitResultPtr r1 = g_4mu->Fit(f1,"S E");
	TFitResultPtr r2 = g_4e->Fit(f2,"S E");
	TFitResultPtr r3 = g_2e2mu->Fit(f3,"S E");
	
	f1 = g_4mu->GetFunction("f1");
	f2 = g_4e->GetFunction("f2");
	f3 = g_2e2mu->GetFunction("f3");
	
	f1->SetLineColor(7);
	f1->SetLineWidth(2);
	f2->SetLineColor(7);
	f2->SetLineWidth(2);
	f3->SetLineColor(7);
	f3->SetLineWidth(2);
	
	Double_t p1_0 = r1->Value(0);
	Double_t p1_1 = r1->Value(1);
	Double_t p1_2 = r1->Value(2);
	Double_t chi2_1 = r1->Chi2();
		
	Double_t p2_0 = r2->Value(0);
	Double_t p2_1 = r2->Value(1);
	Double_t p2_2 = r2->Value(2);
	Double_t chi2_2 = r2->Chi2();

	Double_t p3_0 = r3->Value(0);
	Double_t p3_1 = r3->Value(1);
	Double_t p3_2 = r3->Value(2);
	Double_t chi2_3 = r3->Chi2();

	TString str_chi2_1 = to_string(chi2_1/3);
	TString str_chi2_2 = to_string(chi2_2/3);
	TString str_chi2_3 = to_string(chi2_3/3);
	
	
	TString lumi_s = to_string(Lumi).substr(0,5);
	g_4mu->SetTitle("");
	g_4mu->GetXaxis()->SetTitle("m_{H}");
	g_4mu->GetYaxis()->SetTitle("expected yield in "+lumi_s+" fb^{-1}");
	g_4mu->GetHistogram()->SetMaximum(points_2e2mu[4]+0.05*points_2e2mu[4]);
	g_4mu->GetHistogram()->SetMinimum(0);
	g_4mu->GetXaxis()->SetLimits(115,135);
	g_4mu->SetMarkerStyle(20);
	g_4mu->SetMarkerColor(4);
	g_4mu->SetMarkerSize(2);
	g_4e->SetMarkerStyle(22);
	g_4e->SetMarkerColor(4);
	g_4e->SetMarkerSize(2);
	g_2e2mu->SetMarkerStyle(21);
	g_2e2mu->SetMarkerColor(4);
	g_2e2mu->SetMarkerSize(2);
	
	TLegend *legend=new TLegend(0.15,0.75,0.3,0.85);
	legend->AddEntry(g_4mu,"4mu","P");
	legend->AddEntry(g_4e,"4e","P");
	legend->AddEntry(g_2e2mu,"2e2mu","P");
	legend->SetTextSize(0.03);
	legend->SetLineWidth(0);
	legend->SetFillColor(0);
	legend->SetBorderSize();
	

	TCanvas* c = new TCanvas("c","c", 1000, 1000);
	c->cd();
	c->SetLeftMargin(0.14);
	g_4mu->Draw("ap");
	g_4e->Draw("p");
	g_2e2mu->Draw("p");	
	legend->Draw("same");
	TLatex *latex=new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.05);
    latex->SetTextFont(10);
    latex->SetTextAlign(23);
  //       latex->DrawLatex(0.7, 0.26, "4mu chi2/ndf: "+str_chi2_1);
  //       latex->DrawLatex(0.7, 0.22, "4e chi2/ndf: "+str_chi2_2);
  //       latex->DrawLatex(0.7, 0.18, "2e2mu chi2/ndf: "+str_chi2_3);

	c->SaveAs("/home/chenguan/public_html/EventYield/"+year+"/"+process+".png");
	
	ofstream in;
	in.open("/home/chenguan/public_html/EventYield/"+year+"/"+process+".txt");
	in<<"4mu"<<" "<<p1_0<<" "<<p1_1<<" "<<p1_2<<endl;
	in<<"4e"<<" "<<p2_0<<" "<<p2_1<<" "<<p2_2<<endl;
	in<<"2e2mu"<<" "<<p3_0<<" "<<p3_1<<" "<<p3_2<<endl;
	
	cout<<"4mu "<<"("<<p1_0<<")"<<"+("<<p1_1<<")*mH+("<<p1_2<<")*mH^2"<<"="<<(p1_0)+(p1_1*125)+(p1_2*125*125)<<endl;
	cout<<"4e "<<"("<<p2_0<<")+("<<p2_1<<")*mH+("<<p2_2<<")*mH^2"<<"="<<(p2_0)+(p2_1*125)+(p2_2*125*125)<<endl;
	cout<<"2e2mu "<<"("<<p3_0<<")+("<<p3_1<<")*mH+("<<p3_2<<")*mH^2"<<"="<<(p3_0)+(p3_1*125)+(p3_2*125*125)<<endl;
	
	delete c;
	delete legend;
	delete[] x;
	delete f1;
	delete f2;
	delete f3;
	delete g_4mu;
	delete g_4e;
	delete g_2e2mu;
	delete[] Err_points_2e2mu;
	delete[] Err_points_4e;
	delete[] Err_points_4mu;
	delete latex;	
	
}

Double_t* ReadTree(TTree* t, Int_t fs1, Int_t fs2, Int_t masspoint, TString process){
	
	if(year=="2016")Lumi=35.92;
    if(year=="2017")Lumi=41.53;
    if(year=="2018")Lumi=59.68;
	
	Double_t xs;
	if(process=="GGH"){
		if(masspoint==120)xs=0.008663/0.00792;
		if(masspoint==124)xs=0.012327/0.01122;
		if(masspoint==125)xs=0.013335/0.01212;
		if(masspoint==126)xs=0.014372/0.01308;
		if(masspoint==130)xs=0.018686/0.01706;
	}
	if(process=="VBF"){
		if(masspoint==120)xs=0.000653/0.0006573;
		if(masspoint==124)xs=0.000954/0.0009606;
		if(masspoint==125)xs=0.001038/0.001034;
		if(masspoint==126)xs=0.001126/0.001133;
		if(masspoint==130)xs=0.0015/0.001506;
	}
	if(process=="WminusH"){
		if(masspoint==120)xs=0.000101/0.0001015;
		if(masspoint==124)xs=0.000137/0.0001379;
		if(masspoint==125)xs=0.000146/0.0001471;
		if(masspoint==126)xs=0.000156/0.000157;
		if(masspoint==130)xs=0.000193/0.0001944;
	}
	if(process=="WplusH"){
		if(masspoint==120)xs=0.000159/0.0001606;
		if(masspoint==124)xs=0.000215/0.000219;
		if(masspoint==125)xs=0.000231/0.0002339;
		if(masspoint==126)xs=0.000246/0.0002497;
		if(masspoint==130)xs=0.000306/0.0003103;
	}
	if(process=="ZH"){
		if(masspoint==120)xs=0.000448/0.0004467;
		if(masspoint==124)xs=0.000615/0.0006145;
		if(masspoint==125)xs=0.000662/0.0006569;
		if(masspoint==126)xs=0.000708/0.0007026;
		if(masspoint==130)xs=0.000892/0.0008778;
	}
	if(process=="ttH"){
		if(masspoint==120)xs=0.000262/0.00021383;
		if(masspoint==124)xs=0.000363/0.00034867;
		if(masspoint==125)xs=0.00039/0.00038991;
		if(masspoint==126)xs=0.000418/0.00041842;
		if(masspoint==130)xs=0.000529/0.00046059;
	}

	Double_t Yeild_w = 0;//w means weight
	Double_t Yeild_e = 0;//e means event //the Yeild_e and weight_ave just for poisson error
	Double_t weight_ave = 0;
	
	Int_t Ntot = 0;
	Int_t NLoop = t->GetEntries();
	Ntot = NLoop;

//cout<<NLoop<<endl;	
	for(Int_t i=0; i<NLoop; i++){
		
		t->GetEntry(i);
		
		if(mass4l<140&&mass4l>105&&passedFullSelection==1&&(finalState==fs1||finalState==fs2)){
			
			Yeild_w=Yeild_w+eventWeight*xs;

			Yeild_e=Yeild_e+1;
		
		}
	
	}

	weight_ave = Yeild_w/Yeild_e;
	Yeild_e = sqrt(Yeild_e)*weight_ave*Lumi*1000/Ntot;
	Yeild_w = (Yeild_w*Lumi*1000)/Ntot;
		
	Double_t* results = new Double_t[2];
	
	results[0]=Yeild_w;
	results[1]=Yeild_e;
//    cout<<Yeild_e<<"  "<<Yeild_w<<endl;
	return results;

	
}

void SetTreeAddress(TTree* t){
	t->SetBranchStatus("*",0);
	t->SetBranchStatus("EventCat",1);
	t->SetBranchStatus("mass4l",1);
	t->SetBranchStatus("passedFullSelection",1);
	t->SetBranchStatus("finalState",1);
	t->SetBranchStatus("eventWeight",1);
	t->SetBranchStatus("k_ggZZ",1);
	t->SetBranchStatus("k_qqZZ_ewk",1);
	t->SetBranchStatus("k_qqZZ_qcd_M",1);
	
	t->SetBranchAddress("EventCat",&EventCat);
	t->SetBranchAddress("mass4l",&mass4l);
	t->SetBranchAddress("passedFullSelection",&passedFullSelection);
	t->SetBranchAddress("finalState",&finalState);
	t->SetBranchAddress("eventWeight",&eventWeight);
//	t->SetBranchAddress("etaL1",&etaL1);
//	t->SetBranchAddress("etaL2",&etaL2);
//	t->SetBranchAddress("etaL3",&etaL3);
//	t->SetBranchAddress("etaL4",&etaL4);
//	t->SetBranchAddress("pTL1",&pTL1);
//	t->SetBranchAddress("pTL2",&pTL2);
//	t->SetBranchAddress("pTL3",&pTL3);
//	t->SetBranchAddress("pTL4",&pTL4);
	t->SetBranchAddress("k_qqZZ_ewk",&k_qqZZ_ewk);
	t->SetBranchAddress("k_qqZZ_qcd_M",&k_qqZZ_ewk);
	t->SetBranchAddress("k_ggZZ",&k_ggZZ);



}
