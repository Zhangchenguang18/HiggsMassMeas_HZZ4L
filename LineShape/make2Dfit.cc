TString year = "2018";
TString type = "RECO";
void SetTreeAddress(TTree* t);
RooDataHist* MakeDataSet(TTree* t, Int_t fs1, Int_t fs2);
Float_t mass4l_tmp, mass4lErr_tmp, mass4lREFIT_tmp, mass4lErrREFIT_tmp;
bool passedFullSelection;
Int_t finalState;
Double_t low, high, a1_1D, n1_1D, a2_1D, n2_1D, a1_2D, n1_2D, a2_2D, n2_2D;
void make2Dfit(){

	TFile* f = new TFile("/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/"+year+"GGH_125.root");
	TTree* t = (TTree*)f->Get("passedEvents");
	SetTreeAddress(t);
	Int_t fs1, fs2;
	TString fs_s;
	for(Int_t i=0; i<3; i++){
		if(i==0){fs1=1;fs2=1;fs_s="4mu";low=0.006;high=0.025;}
		if(i==1){fs1=2;fs2=2;fs_s="4e";low=0.006;high=0.05;}
		if(i==2){fs1=3;fs2=4;fs_s="2e2mu";low=0.006;high=0.04;}
		
		RooDataHist* dataset = MakeDataSet(t,fs1,fs2);
		RooRealVar mH("mH","mH",125);
		mH.setConstant(kTRUE);
		RooRealVar mass4l("mass4l","mass4l",105,140);
		RooRealVar mass4lErr("mass4lErr","mass4lErr",low,high);
		RooArgSet argset(mass4lErr);

		RooRealVar mean("mean","mean",125,120,130);
		RooRealVar sigma("sigma","sigma",1,0,10);
		RooRealVar a1("a1","a1",1,0,10);
		RooRealVar n1("n1","n1",1,0,10);
		RooRealVar a2("a2","a2",1,0,10);
		RooRealVar n2("n2","n2",1,0,10);
		RooDoubleCB DCB("DCB","DCB",mass4l,mean,sigma,a1,n1,a2,n2);
		DCB.Print();
		
		RooRealVar mean2d("mean2d","mean2d",125,120,130);
		RooFormulaVar sigma2d("sigma2d","sigma2d","@0*@1",RooArgList(mass4lErr,mH));
		RooRealVar a12d("a12d","a12d",1,0,10);
		RooRealVar n12d("n12d","n12d",1,0,10);
		RooRealVar a22d("a22d","a22d",1,0,10);
		RooRealVar n22d("n22d","n22d",1,0,10);
		RooDoubleCB DCB2d("DCB2d","DCB2d",mass4l,mean2d,sigma2d,a12d,n12d,a22d,n22d);

//		TFile* f = new TFile("/home/chenguan/public_html/EBETemplate/"+year+"/"+type+"/"+"Dm_signal_"+fs_s+".root");
//		TH1F* h = (TH1F*)f->Get("h_Dm");
//		RooDataHist* hist = new RooDataHist("hist","hist",argset,h);
//		RooHistPdf* pdf = new RooHistPdf("pdf","pdf",argset,*hist);
//		//pdf->Print();
//		RooProdPdf* model = new RooProdPdf("model","model",RooArgSet(*pdf),Conditional(RooArgSet(DCB2d),RooArgSet(mass4l)));
//		cout<<"model print!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
//		model->Print();
		
		DCB.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1));
		a1_1D = a1.getVal();
		n1_1D = n1.getVal();
		a2_1D = a2.getVal();
		n2_1D = n2.getVal();
		RooPlot* frame = mass4l.frame(Bins(100));
		frame->SetTitle("");
		dataset->plotOn(frame,MarkerColor(kBlack));
		DCB.plotOn(frame,LineColor(kBlue),LineWidth(6));
		DCB.paramOn(frame,Layout(0.1,0.4,0.9),FillColor(kBlue));

		//mass4l dimentional
		DCB2d.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1),ConditionalObservables(mass4lErr));
		DCB2d.paramOn(frame,Layout(0.1,0.2,0.5));
		a1_2D = a12d.getVal();
		n1_2D = n12d.getVal();
		a2_2D = a22d.getVal();
		n2_2D = n22d.getVal();
		DCB2d.plotOn(frame,LineColor(kRed),ProjWData(mass4lErr,*dataset));
		TCanvas c("c","c",1400,1000);
		c.cd();
		frame->Draw();
		c.SaveAs("/home/chenguan/public_html/LineShape/"+year+"/2DFits_/"+fs_s+"_"+type+".png");

		//mass4lErr dimentional
		RooPlot* frame1 = mass4lErr.frame(Bins(100));
		dataset->plotOn(frame1,MarkerColor(kBlue));
		DCB2d.plotOn(frame1,LineColor(kBlue));
		c.Clear();
		c.cd();
		frame1->Draw();
		c.SaveAs("/home/chenguan/public_html/LineShape/"+year+"/2DFits_/"+fs_s+"_"+type+"masserr.png");	

		//text infor		
		ofstream out;
		out.open("/home/chenguan/public_html/LineShape/"+year+"/2DFits_/"+fs_s+"_"+type+".txt");
		out<<a1_2D<<endl;
		out<<n1_2D<<endl;
		out<<a2_2D<<endl;
		out<<n2_2D<<endl;
		out.close();

	}


}

RooDataHist* MakeDataSet(TTree* t, Int_t fs1, Int_t fs2){
	
	RooRealVar mass4l("mass4l","mass4l",105,140);
	RooRealVar mass4lErr("mass4lErr","mass4lErr",low,high);
	RooArgSet argset(mass4l,mass4lErr);
	TH2F* h = new TH2F("h","h",200,105,140,100,0,0.1);
	for(Int_t i=0; i<t->GetEntries(); i++){
		t->GetEntry(i);
		if(type=="RECO"&&mass4l_tmp<140&&mass4l_tmp>105&&passedFullSelection==1&&(finalState==fs1||finalState==fs2)){
			h->Fill(mass4l_tmp,mass4lErr_tmp/mass4l_tmp);
		}
		if(type=="REFIT"&&mass4lREFIT_tmp<140&&mass4lREFIT_tmp>105&&passedFullSelection==1&&(finalState==fs1||finalState==fs2)){
			h->Fill(mass4lREFIT_tmp,mass4lErrREFIT_tmp/mass4lREFIT_tmp);
		}
	}
	RooDataHist* dataset_h = new RooDataHist("hist","hist",argset,Import(*h));
	cout<<"print 2D dataset !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	dataset_h->Print();
	return dataset_h;

}

void SetTreeAddress(TTree* t){
	t->SetBranchAddress("passedFullSelection",&passedFullSelection);
	t->SetBranchAddress("finalState",&finalState);
	t->SetBranchAddress("mass4l",&mass4l_tmp);
	t->SetBranchAddress("mass4lErr",&mass4lErr_tmp);
	t->SetBranchAddress("mass4lREFIT",&mass4lREFIT_tmp);
	t->SetBranchAddress("mass4lErrREFIT",&mass4lErrREFIT_tmp);
}
