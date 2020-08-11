using namespace RooFit;

bool passedFullSelection;
Int_t finalState, EventCat;
Float_t mass4l, mass4lErr, dataMCWeight, mass4lREFIT, pTL1, pTL2, pTL3, pTL4, etaL1, etaL2, etaL3, etaL4;

void setTreeaddress(TTree* t);
RooDataSet* doDataset(TString inputpath, Int_t fs1, Int_t fs2, TString production);
void ZZLineShape(TString inputpath, Int_t fs1, Int_t fs2, TString production);



TString year = "2016";
bool useREFIT = 0;//based on HIG16041 datacard, always useREFIT = 0

TString plotpath = "/home/chenguan/public_html/LineShape/"+year+"/ZZShape/";



void ZZShape(){
	
	
	TString inputpath="/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/";//for 2016
	
	
	for(Int_t i=0; i<3; i++){
		Int_t fs1, fs2;
		if(i==0){fs1=1; fs2=1;}
		if(i==1){fs1=2; fs2=2;}
		if(i==2){fs1=3; fs2=4;}
		
		for(Int_t j=0; j<2; j++){
			TString production;
			if(j==0)production="qqZZ";
			if(j==1)production="ggZZ";
			ZZLineShape(inputpath, fs1, fs2, production);
		}
	}
			


}

void ZZLineShape(TString inputpath, Int_t fs1, Int_t fs2, TString production){
	
	TString str_fs;
	if(fs1==1)str_fs = "4mu";
	if(fs1==2)str_fs = "4e";
	if(fs1>2)str_fs = "2e2mu";
	
	RooDataSet* dataset = doDataset(inputpath, fs1, fs2, production);
	
	RooRealVar mass4l_d("mass4l_d","mass4l_d",105,140);
	RooRealVar c1_init("c1_init","c1_init",0,3);
	RooRealVar c2_init("c2_init","c2_init",0,3);
	RooRealVar c3_init("c3_init","c3_init",0,3);
	
	RooBernstein ZZ_init("ZZ_init","ZZ_init",mass4l_d,RooArgList(c1_init,c2_init,c3_init));
	
	ZZ_init.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1),Timer(kTRUE));
	
	Double_t c1_, c2_, c3_, c1Err, c2Err, c3Err;
	c1_=c1_init.getVal(); c2_=c2_init.getVal(); c3_=c3_init.getVal();
	c1Err=c1_init.getError(); c2Err=c2_init.getError(); c3Err=c3_init.getError();
	
	RooRealVar c1("c1","c1",c1_,c1_-c1Err,c1_+c1Err);
	RooRealVar c2("c2","c2",c2_,c2_-c2Err,c2_+c2Err);
	RooRealVar c3("c3","c3",c3_,c3_-c3Err,c3_+c3Err);
	
	RooBernstein ZZ("ZZ","ZZ",mass4l_d,RooArgList(c1,c2,c3));
	ZZ.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1),Timer(kTRUE));
	
	RooPlot* frame = mass4l_d.frame(Bins(50));
	frame->SetTitle("");
	frame->GetXaxis()->SetTitle("mass4l(GeV)");
	frame->GetYaxis()->SetTitleOffset(1.55);
	dataset->plotOn(frame,MarkerStyle(24));
	ZZ.plotOn(frame,LineColor(kBlue));
	//ZZ.paramOn(frame,Layout(0.1,0.1,0.9));
	Double_t chisquare=frame->chiSquare(3);
	TLatex *latex=new TLatex();
	latex->SetNDC();
	latex->SetTextSize(0.04);
	latex->SetTextFont(42);
	latex->SetTextAlign(13);
	char chi2[20];
	sprintf(chi2,"%s%1.4f","#chi^{2}/DOF=",chisquare);
	TString str_c1 = to_string(c1.getVal());
	TString str_c2 = to_string(c2.getVal());
	TString str_c3 = to_string(c3.getVal());
	TString str_c1Err = to_string(c1.getError());
	TString str_c2Err = to_string(c2.getError());
	TString str_c3Err = to_string(c3.getError());

	TCanvas c("c","c",1400,1000);
	c.cd();
	frame->Draw();
	
	latex->DrawLatex(0.1, 0.54, chi2);
	latex->DrawLatex(0.1, 0.44, "c1 = "+str_c1+" +/- "+str_c1Err);
	latex->DrawLatex(0.1, 0.34, "c2 = "+str_c2+" +/- "+str_c2Err);
	latex->DrawLatex(0.1, 0.24, "c3 = "+str_c3+" +/- "+str_c3Err);
 	
	cout<<"**************************************************"<<endl;	
	cout<<c1.getVal()<<" "<<c2.getVal()<<" "<<c3.getVal()<<endl;
	
	TString REFIT;
	if(useREFIT)REFIT="_REFIT";
	if(!useREFIT)REFIT="";
	c.SaveAs(plotpath+production+"_"+str_fs+REFIT+".png");
	
	ofstream fout;
	fout.open(plotpath+production+"_"+str_fs+REFIT+".txt");
	if(str_fs=="2e2mu"){
		str_c1 = "'chebPol1_3' : "+str_c1;
		str_c2 = "'chebPol2_3' : "+str_c2;
		str_c3 = "'chebPol3_3' : "+str_c3;
	}
	if(str_fs=="4e"){
		str_c1 = "'chebPol1_2' : "+str_c1;
		str_c2 = "'chebPol2_2' : "+str_c2;
		str_c3 = "'chebPol3_2' : "+str_c3;
	}
	if(str_fs=="4mu"){
		str_c1 = "'chebPol1_1' : "+str_c1;
		str_c2 = "'chebPol2_1' : "+str_c2;
		str_c3 = "'chebPol3_1' : "+str_c3;
	}
	fout<<str_c1<<endl;
	fout<<str_c2<<endl;
	fout<<str_c3<<endl;
	fout.close();
	
	
	delete latex;
	delete dataset;
	delete frame;
	
}


RooDataSet* doDataset(TString inputpath, Int_t fs1, Int_t fs2, TString production){
	
	Int_t sum;

	TString str_fs;
	if(fs1==1)str_fs = "4mu";
	if(fs1==2)str_fs = "4e";
	if(fs1>2)str_fs = "2e2mu";
	
	TString filename;
	if(production=="qqZZ")filename="slimmed_"+year+"qqZZ.root";
	if(production=="ggZZ")filename=year+"ggZZ_"+str_fs+".root";

	RooRealVar *dataMCWeight_d = new RooRealVar("dataMCWeight_d","dataMCWeight_d",0.0000001,1000);
	RooRealVar *mass4l_d = new RooRealVar("mass4l_d","mass4l_d",105,140); 
	RooDataSet *dataset4l_d = new RooDataSet("dataset4l_d","dataset4l_d",RooArgSet(*mass4l_d,*dataMCWeight_d));
	
	TFile* f = new TFile(inputpath+filename);
	TTree* t = (TTree*)f->Get("passedEvents");
	setTreeaddress(t);
	sum=t->GetEntries();
	for(Int_t i=0; i<sum; i++){
		t->GetEntry(i);
		if(passedFullSelection==1&&(finalState==fs1||finalState==fs2)&&mass4l>105&&mass4l<140){
			if(useREFIT)*mass4l_d=mass4lREFIT;
			if(!useREFIT)*mass4l_d=mass4l;
			*dataMCWeight_d=dataMCWeight;
			dataset4l_d->add(RooArgSet(*mass4l_d,*dataMCWeight_d));
		}
	}
	RooDataSet *dataset_w=new RooDataSet(dataset4l_d->GetName(),dataset4l_d->GetTitle(),dataset4l_d,*dataset4l_d->get(),"1","dataMCWeight_d");
	delete t;
	delete dataMCWeight_d;
	delete mass4l_d;
	delete dataset4l_d;
	
	
	return dataset_w;
		
	
	
	
}

void setTreeaddress(TTree* t){
	
	t->SetBranchAddress("mass4l",&mass4l);
	t->SetBranchAddress("mass4lREFIT",&mass4lREFIT);
	t->SetBranchAddress("passedFullSelection",&passedFullSelection);
	t->SetBranchAddress("finalState",&finalState);
	t->SetBranchAddress("dataMCWeight",&dataMCWeight);
	t->SetBranchAddress("pTL1",&pTL1);
	t->SetBranchAddress("pTL2",&pTL2);
	t->SetBranchAddress("pTL3",&pTL3);
	t->SetBranchAddress("pTL4",&pTL4);
	t->SetBranchAddress("etaL1",&etaL1);
	t->SetBranchAddress("etaL2",&etaL2);
	t->SetBranchAddress("etaL3",&etaL3);
	t->SetBranchAddress("etaL4",&etaL4);
	t->SetBranchAddress("EventCat",&EventCat);

}
