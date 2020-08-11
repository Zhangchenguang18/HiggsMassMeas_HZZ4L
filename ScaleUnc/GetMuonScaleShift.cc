#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "TFile.h"
#include "TLine.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TH1F.h"
using namespace RooFit;
Double_t massZ, massZErr, eta1, eta2, pT1, pT2, pterr1, pterr2, m1, m2, phi1, phi2, RelIso1, RelIso2, weight; Int_t lep1_ecalDriven, lep2_ecalDriven;
TString fs = "muon";

TString year = "2018";

//set up range
Double_t low = 80;
Double_t high = 100;

RooDataHist* MakeDataSet(Double_t* ptcut, Double_t* etacut, bool isData,Int_t seed, TString Sample);
//void DataMCComp(RooDataHist* mc, RooDataHist* data, TString name);
Double_t* GetMean(RooDataHist* dataset, TString name, bool isData, bool fixtail, TString Sample);
void ClosurePlot(vector<vector<Double_t*>> &Z, vector<vector<Double_t*>> &Jpsi);
void SetAddress(TTree* t);
void MakeLUT(Double_t* sc1, Double_t* sc2, Double_t* sc3);
void GetScaleShift_(TString Sample, vector<vector<Double_t*>> &results);

void GetMuonScaleShift(){

	vector<vector<Double_t*>> Z; Z.clear();
	GetScaleShift_("Z", Z);
	vector<vector<Double_t*>> JPsi; JPsi.clear();
//        GetScaleShift_("JPsi", JPsi);
	
	ClosurePlot(Z, JPsi);

}


void GetScaleShift_(TString Sample, vector<vector<Double_t*>> &results){	
	
	Double_t mass;
	Int_t pt_edges, eta_edges;
	if(Sample=="Z"){low=80;high=100;mass=91.19;pt_edges=7;}
//	if(Sample=="Upsilon"){low=8;high=10;mass=9.46;pt_edges=5;}
	if(Sample=="JPsi"){low=2.95;high=3.25;mass=3.09;pt_edges=4;}

	Double_t etabin[4] = {0,0.9,1.4,2.4};
	
	Double_t* ptbin = new Double_t[pt_edges];
	if(fs=="electron")ptbin[0]=7;
	if(fs=="muon")ptbin[0]=5;

	if(Sample=="Z"){ptbin[1]=20;ptbin[2]=30;ptbin[3]=40;ptbin[4]=50;ptbin[5]=60;ptbin[6]=100;}
	if(Sample=="JPsi"){ptbin[1]=10;ptbin[2]=15;ptbin[3]=30;}
//	if(Sample=="Upsilon"){}

	Double_t* ptcut = new Double_t[2];
	Double_t* etacut = new Double_t[2];
	Double_t* tmp_mc = new Double_t[2];
	Double_t* tmp_data = new Double_t[2];

	for(Int_t i=0; i<3; i++){//eta loop 
		vector<Double_t*> result;result.clear();
		for(Int_t j=0; j<pt_edges-1; j++){//pt loop
			ptcut[0] = ptbin[j];ptcut[1] = ptbin[j+1];
			etacut[0] = etabin[i];etacut[1] = etabin[i+1];
			
			RooDataHist* mc = MakeDataSet(ptcut, etacut, 0, 5*j+i, Sample);
			RooDataHist* data = MakeDataSet(ptcut, etacut, 1, 7*j+i, Sample);
			TString name = "pt_"+to_string(ptcut[0]).substr(0,3)+"_"+to_string(ptcut[1]).substr(0,3)+"_eta_"+to_string(etacut[0]).substr(0,3)+"_"+to_string(etacut[1]).substr(0,3);
			//DataMCComp(mc, data, name);
			tmp_mc = GetMean(mc, name, 0, 1, Sample);
			tmp_data = GetMean(data, name ,1, 1, Sample);
			
			Double_t* point = new Double_t[4];
			point[0] = (tmp_data[0]-tmp_mc[0])/mass;
			point[1] = sqrt(tmp_mc[1]*tmp_mc[1]+tmp_data[1]*tmp_data[1])/mass;
			point[2] = 0.5*abs(ptcut[1]+ptcut[0]);
			point[3] = 0.5*abs(ptcut[1]-ptcut[0]);
			result.push_back(point);
			delete mc; delete data;
		}
		results.push_back(result);
	}

}

/*
void DataMCComp(RooDataHist* mc, RooDataHist* data, TString name){
	RooRealVar massZ("massZ","massZ",low,high);
	RooHistPdf* pdf_mc = new RooHistPdf("pdf_mc","pdf_mc",massZ,*mc);
	RooHistPdf* pdf_data = new RooHistPdf("pdf_data","pdf_data",massZ,*data);
	RooPlot* frame = massZ.frame(Bins(100));
	frame->SetTitle("");
	pdf_mc->plotOn(frame,LineColor(kBlue));
	pdf_data->plotOn(frame,LineColor(kRed));
	TCanvas c("c","c",1400,1000);
	c.cd();
	frame->Draw();
	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+fs+"/"+plots+"/Comparison_"+name+".png");
	delete pdf_mc; delete pdf_data; delete frame;
}
*/

void ClosurePlot(vector<vector<Double_t*>> &Z, vector<vector<Double_t*>> &JPsi){
	
	gStyle->SetTitleOffset(2,"Y");	
	
	Int_t Zptbins = Z.at(1).size();
	Double_t* Zsc1 = new Double_t[Zptbins];Double_t* Zsc1Err = new Double_t[Zptbins];
	Double_t* Zsc2 = new Double_t[Zptbins];Double_t* Zsc2Err = new Double_t[Zptbins];
	Double_t* Zsc3 = new Double_t[Zptbins];Double_t* Zsc3Err = new Double_t[Zptbins];
	Double_t* Zx = new Double_t[Zptbins];Double_t* Zx_error = new Double_t[Zptbins];

	for(Int_t i=0; i<Zptbins; i++){
		Zsc1[i] = Z.at(0).at(i)[0];
		Zsc1Err[i] = Z.at(0).at(i)[1];
		Zsc2[i] = Z.at(1).at(i)[0];
		Zsc2Err[i] = Z.at(1).at(i)[1];
		Zsc3[i] = Z.at(2).at(i)[0];
		Zsc3Err[i] = Z.at(2).at(i)[1];
		Zx[i] = Z.at(0).at(i)[2];
		Zx_error[i] = Z.at(0).at(i)[3];
	}

//	Int_t JPsiptbins = JPsi.at(1).size();
//	Double_t* JPsisc1 = new Double_t[JPsiptbins];Double_t* JPsisc1Err = new Double_t[JPsiptbins];
//	Double_t* JPsisc2 = new Double_t[JPsiptbins];Double_t* JPsisc2Err = new Double_t[JPsiptbins];
//	Double_t* JPsisc3 = new Double_t[JPsiptbins];Double_t* JPsisc3Err = new Double_t[JPsiptbins];
//	Double_t* JPsix = new Double_t[JPsiptbins];Double_t* JPsix_err = new Double_t[JPsiptbins];
//	for(Int_t i=0; i<JPsiptbins; i++){
//		JPsisc1[i] = JPsi.at(0).at(i)[0];
//		JPsisc1Err[i] = JPsi.at(0).at(i)[1];
//		JPsisc2[i] = JPsi.at(1).at(i)[0];
//		JPsisc2Err[i] = JPsi.at(1).at(i)[1];
//		JPsisc3[i] = JPsi.at(2).at(i)[0];
//		JPsisc3Err[i] = JPsi.at(2).at(i)[1];
//		JPsix[i] = JPsi.at(0).at(i)[2];
//		JPsix_err[i] = JPsi.at(0).at(i)[3];
//	}
	Double_t ptmax = Zx[Zptbins-1]+Zx_error[Zptbins-1];

	MakeLUT(Zsc1,Zsc2,Zsc3);

	TGraph* g1 = new TGraphErrors(Zptbins,Zx,Zsc1,Zx_error,Zsc1Err);
	TGraph* g2 = new TGraphErrors(Zptbins,Zx,Zsc2,Zx_error,Zsc2Err);
	TGraph* g3 = new TGraphErrors(Zptbins,Zx,Zsc3,Zx_error,Zsc3Err);
//	TGraph* g4 = new TGraphErrors(JPsiptbins,JPsix,JPsisc1,JPsix_err,JPsisc1Err);
//	TGraph* g5 = new TGraphErrors(JPsiptbins,JPsix,JPsisc2,JPsix_err,JPsisc2Err);
//	TGraph* g6 = new TGraphErrors(JPsiptbins,JPsix,JPsisc3,JPsix_err,JPsisc3Err);

	g1->SetMarkerStyle(25);
	g1->SetMarkerColor(kBlue);
	g1->SetTitle("");
	g1->GetXaxis()->SetTitle(fs+" p_{T} (GeV)");
	g1->GetYaxis()->SetTitle("(m_{data}-m_{MC})/m_{PDG}");
	g1->GetXaxis()->SetLimits(0,ptmax);
	g1->GetHistogram()->SetMaximum(0.004);
	g1->GetHistogram()->SetMinimum(-0.004);

	g2->SetMarkerStyle(25);
	g2->SetMarkerColor(kBlack);
	g3->SetMarkerStyle(25);
	g3->SetMarkerColor(kRed);
//	g4->SetMarkerStyle(26);
//	g4->SetMarkerColor(kBlue);
//	g5->SetMarkerStyle(26);
//	g5->SetMarkerColor(kBlack);
//	g6->SetMarkerStyle(26);
//	g6->SetMarkerColor(kRed);

	TLine* l1 = new TLine(0,0,ptmax,0);
	l1->SetLineStyle(kDashed);
	TLine* l2 = new TLine(0,0.001,ptmax,0.001);
	l2->SetLineStyle(kDashed);
	TLine* l3 = new TLine(0,-0.001,ptmax,-0.001);
	l3->SetLineStyle(kDashed);
	TLine* l4 = new TLine(0,-0.0004,ptmax,-0.0004);
	l4->SetLineStyle(kDashed);
	TLine* l5 = new TLine(0,0.0004,ptmax,0.0004);
	l5->SetLineStyle(kDashed);

	TLegend* legend = new TLegend(0.7,0.15,0.75,0.35);
        legend->AddEntry(g1, "Z |#eta| 0.0-0.9", "P");
        legend->AddEntry(g2, "Z |#eta| 0.9-1.4", "P");
        legend->AddEntry(g3, "Z |#eta| 1.4-2.4", "P");
//	legend->AddEntry(g4, "J#Psi |#eta| 0.0-0.9", "P");
//	legend->AddEntry(g5, "J#Psi |#eta| 0.9-1.4", "P");
//	legend->AddEntry(g6, "J#Psi |#eta| 1.4-2.4", "P");
      legend->SetTextSize(0.03);
        legend->SetLineWidth(0);
        legend->SetFillColor(0);
        legend->SetBorderSize();
	
	TCanvas* c = new TCanvas("c","c",1000,1000);
	c->cd();
	c->SetLeftMargin(0.14);	
	g1->Draw("ap");
	g2->Draw("p");
	g3->Draw("p");
//	g4->Draw("p");
//	g5->Draw("p");
//	g6->Draw("p");
	legend->Draw("same");
	l1->Draw("same");
	l2->Draw("same");
	l3->Draw("same");
	l4->Draw("same");
	l5->Draw("same");

	c->SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/Closure_test_"+fs+".png");
//	c->SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/Closure_test_"+fs+".pdf");
	delete c; delete g1; delete g2; delete g3; delete legend; delete l1; delete l2; delete l3;

}

Double_t* GetMean(RooDataHist* dataset, TString name, bool isData, bool fixtail, TString Sample){	
	gStyle->SetTitleOffset(2,"Y");

	Double_t width, mean_, mean_min, mean_max, sigma_, sigma_min, sigma_max, n1_, n1_min, n1_max, n2_, n2_min, n2_max;
	if(Sample=="Z"){ 
		if(isData==1){width=2.49;}
		if(isData==0){width=2.44;}
		mean_=0; mean_min=-5; mean_max=5; 
		sigma_=1; sigma_min=0, sigma_max=10;
		n1_=5; n1_min=0; n1_max=60;
		n2_=5; n2_min=0; n2_max=60;
	}
	if(Sample=="JPsi"){
		mean_=3.09; mean_min=3; mean_max=3.15;
		sigma_=0.02; sigma_min=0; sigma_max=0.1;
		n1_=5; n1_min=0; n1_max=60;
		n2_=5; n2_min=0; n2_max=60;
	}
	if(Sample=="Upsilon"){
		mean_=9.46; mean_min=9; mean_max=10;
		sigma_=0; sigma_min=0; sigma_max=0;
		n1_=0; n1_min=0; n1_max=0;
		n2_=0; n2_min=0; n2_max=0;
	}

	RooRealVar massZ("massZ","massZ",low,high);
	RooRealVar PDGmZ("PDGmZ","PDGmZ",91.19);
	RooRealVar PDGwZ("PDGwZ","PDGwZ",width);
	PDGmZ.setConstant(kTRUE);
	PDGwZ.setConstant(kTRUE);
	RooBreitWigner PDGBW("PDGBW","PDGBW",massZ,PDGmZ,PDGwZ);
	
	RooRealVar mean("mean","mean",mean_,mean_min,mean_max);
	RooRealVar sigma("sigma","sigma",sigma_,sigma_min,sigma_max);
	RooRealVar alpha1("alpha1","alpha1",1,0,10);
	RooRealVar n1("n1","n1",n1_,n1_min,n1_max);
	RooRealVar alpha2("alpha2","alpha2",1,0,10);
	RooRealVar n2("n2","n2",n2_,n2_min,n2_max);
	RooDoubleCB DCB("DCB","DCB",massZ,mean,sigma,alpha1,n1,alpha2,n2);
	
	RooFFTConvPdf BW_DCB("BW_DCB","BW_DCB",massZ,PDGBW,DCB);
	RooRealVar tau("tau","tau",0,-1,1);
	RooExponential bkg("bkg","bkg",massZ,tau);
	RooRealVar fsig("fsig","fsig",0.9,0.5,1);
	RooAddPdf Model("Model","Model",BW_DCB,bkg,fsig);
	
	RooAddPdf Model2("Model2","Model2",DCB,bkg,fsig);//no conv for JPsi and Upsilon. Model and Model2 share same params

	if(Sample=="JPsi"){
		Model2.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1));
	}
	if(Sample=="Upsilon"){
		Model2.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1));
	}	
	if(Sample=="Z"){
		Model.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1));
	}
	
//	cout<<"1st mean "<<mean.getVal()<<" "<<mean.getError()<<endl;
//	cout<<"1st n1 "<<n1.getVal()<<endl;	
		
	if(fixtail==1){
		alpha1.setConstant(1);
		n1.setConstant(1);
		alpha2.setConstant(1);
		n2.setConstant(1);
		tau.setConstant(1);
		fsig.setConstant(1);
		mean.setConstant(0);
		sigma.setConstant(0);
		if(Sample!="Z")Model2.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1));
		if(Sample=="Z")Model.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1));
//		cout<<"2ed mean "<<mean.getVal()<<" "<<mean.getError()<<endl;
//		cout<<"2ed n1 "<<n1.getVal()<<endl;
	}
	
	RooPlot* frame = massZ.frame(Bins(100));
	frame->SetTitle("");
	if(fs=="muon")frame->SetXTitle("m_{#mu#mu}");
	dataset->plotOn(frame);
	if(Sample=="Z")Model.plotOn(frame,LineColor(kBlue),LineWidth(2));
	if(Sample!="Z")Model2.plotOn(frame,LineColor(kBlue),LineWidth(2));
	TString mean_s, mean_error_s,sigma_s, alpha1_s, n1_s, alpha2_s, n2_s, chi2_s;
	mean_s = to_string(mean.getVal()); sigma_s = to_string(sigma.getVal());
	alpha1_s = to_string(alpha1.getVal()); n1_s = to_string(n1.getVal());
	alpha2_s = to_string(alpha2.getVal()); n2_s = to_string(n2.getVal());
	mean_error_s = to_string(mean.getError());
	Double_t chi2;
	if(fixtail==0)chi2 = frame->chiSquare(8);
	if(fixtail==1)chi2 = frame->chiSquare(2);
	chi2_s = to_string(chi2);
	TCanvas c("c","c",1400,1000);
	c.SetLeftMargin(0.14);
	frame->Draw();
	TLatex *latex=new TLatex();
	latex->SetNDC();
        latex->SetTextSize(0.05);
        latex->SetTextFont(42);
        latex->SetTextAlign(23);
	latex->DrawLatex(0.7,0.8,"#chi^{2}/DOF="+chi2_s);
	latex->DrawLatex(0.7,0.7,"mean="+mean_s+"+/-"+mean_error_s);
	latex->DrawLatex(0.7,0.6,"sigma="+sigma_s);
	latex->DrawLatex(0.7,0.5,"alpha1="+alpha1_s);
	latex->DrawLatex(0.7,0.4,"n1="+n1_s);
	latex->DrawLatex(0.7,0.3,"alpha2="+alpha2_s);
	latex->DrawLatex(0.7,0.2,"n2="+n2_s);
	TString subname;
	if(isData==0)subname = "MC_";if(isData==1)subname = "Data_";
	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+Sample+"/"+fs+"/"+subname+name+".png");
//	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+Sample+"/"+fs+"/"+subname+name+".pdf");
	delete frame; delete latex;
	Double_t* results = new Double_t[2];
	results[0] = mean.getVal();results[1] = mean.getError();
	return results;	

}

RooDataHist* MakeDataSet(Double_t* ptcut, Double_t* etacut, bool isData, Int_t seed, TString Sample){
	
	RooRealVar massZ_("massZ","massZ",low,high);
	massZ_.setBins(300);
	RooArgSet argset(massZ_);
	//RooDataSet* dataset = new RooDataSet("dataset","dataset",argset);
	TFile *f=new TFile();
	
/*	//HIG-16-041
	if(isData==0&&year=="2016"&&fs=="muon")f = new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/DY_2016MC_v3_20170312_afterApproval/DYJetsToLL_M-50_kalman_v4_m2mu.root");	
	if(isData==1&&year=="2016"&&fs=="muon")f = new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/Data_2016_v1_20170312_afterApproval/DoubleLepton_m2mu.root");
	if(isData==0&&year=="2016"&&fs=="electron")f = new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/DY_2016MC_v3_20170312_afterApproval/DYJetsToLL_M-50_kalman_v4_m2e.root");
	if(isData==1&&year=="2016"&&fs=="electron")f = new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/Data_2016_v1_20170312_afterApproval/DoubleLepton_m2e.root");
	*/
	
	//Run2 Legacy
	if(Sample=="Z"&&isData==0&&year=="2016"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2016DY_m2mu.root");
	if(Sample=="Z"&&isData==0&&year=="2017"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2017DY_m2mu.root");
	if(Sample=="Z"&&isData==0&&year=="2018"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2018DY_m2mu.root");
//	if(Sample=="Z"&&isData==0&&year=="2016"&&fs=="electron")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2016_m2e.root");
//	if(Sample=="Z"&&isData==0&&year=="2017"&&fs=="electron")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2017_m2e.root");
//	if(Sample=="Z"&&isData==0&&year=="2018"&&fs=="electron")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2018_m2e.root");
	
	if(Sample=="Z"&&isData==1&&year=="2016"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/SingleMuon2016_m2mu.root");
	if(Sample=="Z"&&isData==1&&year=="2017"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/SingleMuon2017_m2mu.root");
	if(Sample=="Z"&&isData==1&&year=="2018"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/SingleMuon2018_m2mu.root");
//	if(Sample=="Z"&&isData==1&&year=="2016"&&fs=="electron")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2016/Electron/SingleElectron_m2e.root");
//	if(Sample=="Z"&&isData==1&&year=="2017"&&fs=="electron")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2017/Electron/SingleElectron_m2e.root");
//	if(Sample=="Z"&&isData==1&&year=="2018"&&fs=="electron")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2018/Electron/EGamma_m2e.root");
	
	//JPsi
//	if(Sample=="JPsi"&&isData==0&&year=="2016"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/JPsiandUpsilon/JPsiMuMu_2016_m2mu.root");
//	if(Sample=="JPsi"&&isData==0&&year=="2017"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/JPsiandUpsilon/JPsiMuMu_2017_m2mu.root");
//	if(Sample=="JPsi"&&isData==0&&year=="2018"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/JPsiandUpsilon/JPsiMuMu_2018_m2mu.root");
//	
//	if(Sample=="JPsi"&&isData==1&&year=="2016"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2016/Muon/SingleMuon_m2mu.root");
//	if(Sample=="JPsi"&&isData==1&&year=="2017"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2017/Muon/SingleMuon_m2mu.root");
//	if(Sample=="JPsi"&&isData==1&&year=="2018"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2018/Muon/SingleMuon_m2mu.root");

	//Upsilon
	if(Sample=="Upsilon"&&isData==0&&year=="2016"&&fs=="muon")f = new TFile("");
	if(Sample=="Upsilon"&&isData==0&&year=="2017"&&fs=="muon")f = new TFile("");
	if(Sample=="Upsilon"&&isData==0&&year=="2018"&&fs=="muon")f = new TFile("");

	if(Sample=="Upsilon"&&isData==1&&year=="2016"&&fs=="muon")f = new TFile("");
	if(Sample=="Upsilon"&&isData==1&&year=="2017"&&fs=="muon")f = new TFile("");
	if(Sample=="Upsilon"&&isData==1&&year=="2018"&&fs=="muon")f = new TFile("");

	TTree* t = (TTree*)f->Get("passedEvents");
	Int_t sum = t->GetEntries();
	SetAddress(t);
	TH1D* h = new TH1D("h","h",300,low,high);
	TRandom3 rand;
	rand.SetSeed(12345679);
	Double_t randcut;
	Double_t factor1, factor2;
	for(Int_t i=0; i<sum; i++){
		t->GetEntry(i);
		randcut=rand.Gaus(0,1);
		if(massZ<high&&massZ>low&&RelIso1<0.35&&RelIso2<0.35){
			if((pT1>ptcut[0]&&pT1<ptcut[1]&&eta1<etacut[1]&&eta1>etacut[0])||(eta2<etacut[1]&&eta2>etacut[0]&&pT2>ptcut[0]&&pT2<ptcut[1])){
				h->Fill(massZ,weight);
			}
		}		
	}
	RooDataHist* dataset = new RooDataHist("dataset","dataset",argset,Import(*h));
	f->Close();
	delete f; 
	return dataset;

}
				
void SetAddress(TTree* t){
	t->SetBranchStatus("*",0);
	t->SetBranchStatus("massZ",1);
	t->SetBranchStatus("massZErr",1);
	t->SetBranchStatus("eta1",1);
	t->SetBranchStatus("eta2",1);
	t->SetBranchStatus("pT1",1);
	t->SetBranchStatus("pT2",1);
	t->SetBranchStatus("pterr1",1);
	t->SetBranchStatus("pterr2",1);
	t->SetBranchStatus("m1",1);
	t->SetBranchStatus("m2",1);
	t->SetBranchStatus("phi1",1);
	t->SetBranchStatus("phi2",1);
	t->SetBranchStatus("lep1_ecalDriven",1);
	t->SetBranchStatus("lep2_ecalDriven",1);
	t->SetBranchStatus("RelIso1",1);
	t->SetBranchStatus("RelIso2",1);
	t->SetBranchStatus("weight",1);

	t->SetBranchAddress("massZ",&massZ);
	t->SetBranchAddress("massZErr",&massZErr);
	t->SetBranchAddress("eta1",&eta1);
	t->SetBranchAddress("eta2",&eta2);
	t->SetBranchAddress("pT1",&pT1);
	t->SetBranchAddress("pT2",&pT2);
	t->SetBranchAddress("pterr1",&pterr1);
	t->SetBranchAddress("pterr2",&pterr2);
	t->SetBranchAddress("m1",&m1);
	t->SetBranchAddress("m2",&m2);
	t->SetBranchAddress("phi1",&phi1);
	t->SetBranchAddress("phi2",&phi2);
	t->SetBranchAddress("lep1_ecalDriven",&lep1_ecalDriven);
	t->SetBranchAddress("lep2_ecalDriven",&lep2_ecalDriven);
	t->SetBranchAddress("RelIso1",&RelIso1);
	t->SetBranchAddress("RelIso2",&RelIso2);
	t->SetBranchAddress("weight",&weight);

}		

void MakeLUT(Double_t* sc1, Double_t* sc2, Double_t* sc3){
	

	Double_t etabin[4] = {0,0.9,1.4,2.4};
	Double_t ptbin[7] = {5,20,30,40,50,60,100};
	TH2D* LUT = new TH2D("LUT","LUT",6,ptbin,3,etabin);
	for(Int_t j=0; j<6; j++){
		LUT->SetBinContent(j+1,1,sc1[j]);
		LUT->SetBinContent(j+1,2,sc2[j]);
		LUT->SetBinContent(j+1,3,sc3[j]);
	}
	TCanvas c("c","c",1400,1000);
	c.cd();
        LUT->Draw("TEXT");
        c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/LUT_"+fs+"_.png");
//	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/LUT_"+fs+"_.pdf");
        TFile* f = new TFile("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/LUT_"+fs+"_.root","RECREATE");
        f->cd();
        LUT->Write();
        f->Close();

}
