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

TString year = "2018";
bool useREFIT = 1;
TString fs_="4mu";

bool passedFullSelection;
Float_t mass4l;
Float_t mass4lREFIT;
Int_t finalState, EventCat;
Float_t dataMCWeight, pTL1, pTL2, pTL3, pTL4, etaL1, etaL2, etaL3, etaL4;

Double_t* doFirstfit(TString inputpath, Int_t fs1, Int_t fs2, TString mh);
RooDataSet* doDataset(TString inputpath, Int_t fs1, Int_t fs2, TString mh);
Double_t* doSimultaneousfit(TString inputpath, Int_t fs1, Int_t fs2, Double_t* Parameters);
void setTreeaddress(TTree* t);
void makeplots(TString inputpath, Int_t fs1, Int_t fs2, Double_t* Shape);

void DCB_Parametrization(){
	
	TString inputpath="/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/";
	
	
		Int_t fs1, fs2;
		if(fs_=="4mu"){fs1=1; fs2=1;}
		if(fs_=="4e"){fs1=2; fs2=2;}
		if(fs_=="2e2mu"){fs1=3; fs2=4;}
		Double_t* Parameters = doFirstfit(inputpath, fs1, fs2, "125");
		Double_t* Shape = doSimultaneousfit(inputpath, fs1, fs2, Parameters);//doFirstfit and doSimultaneousfit finish all jobs need to be done, makeplots just like do closure test 
		makeplots(inputpath, fs1, fs2, Shape);

	

}

void makeplots(TString inputpath, Int_t fs1, Int_t fs2, Double_t* Shape){//画随着质量点变化的DCB，本来也可以画Hualin 的结果，被注释掉了。
	
	TString str_fs;
	if(fs1==1)str_fs = "4mu";
	if(fs1==2)str_fs = "4e";
	if(fs1>2)str_fs = "2e2mu";
	
	TString str_inter_mean = to_string(Shape[0]);
	TString str_inter_sigma = to_string(Shape[1]);
	TString str_inter_alpha1 = to_string(Shape[2]);
	TString str_inter_n1 = to_string(Shape[3]);
	TString str_inter_alpha2 = to_string(Shape[4]);
	TString str_inter_n2 = to_string(Shape[5]);
	TString str_slope_mean = to_string(Shape[6]);
	TString str_slope_sigma = to_string(Shape[7]);
	TString str_slope_alpha1 = to_string(Shape[8]);
	TString str_slope_n1 = to_string(Shape[9]);
	TString str_slope_alpha2 = to_string(Shape[10]);
	TString str_slope_n2 = to_string(Shape[11]);
	
	TString mean_form = str_inter_mean+"+"+"("+str_slope_mean+")"+"*(@0-125)";
	TString sigma_form = str_inter_sigma+"+"+"("+str_slope_sigma+")"+"*(@0-125)";
	TString alpha1_form = str_inter_alpha1+"+"+"("+str_slope_alpha1+")"+"*(@0-125)";
	TString n1_form = str_inter_n1+"+"+"("+str_slope_n1+")"+"*(@0-125)";
	TString alpha2_form = str_inter_alpha2+"+"+"("+str_slope_alpha2+")"+"*(@0-125)";
	TString n2_form = str_inter_n2+"+"+"("+str_slope_n2+")"+"*(@0-125)";
	
	RooRealVar mass4l_d("mass4_d","mass4l_d",105,140);
	RooPlot* frame = mass4l_d.frame();
	frame->SetTitle("");
	
	Double_t MH;
	TString mH;
	Double_t MH_tmp[11]={120,121,122,123,124,125,126,127,128,129,130};
	
	for(Int_t i=0; i<=10; i++){
		MH = MH_tmp[i];
		mH = to_string((int)(MH));
	
		RooRealVar mass("mass","mass",120,130);
		mass.setVal(MH);
		mass.setConstant(kTRUE);
	 
	
		RooFormulaVar meanDCB_v("meanDCB_v",mean_form,RooArgList(mass));
		RooFormulaVar sigmaDCB_v("sigmaDCB_v",sigma_form,RooArgList(mass));
		RooFormulaVar alphaDCB_v("alphaDCB_v",alpha1_form,RooArgList(mass));
		RooFormulaVar nDCB_v("nDCB_v",n1_form,RooArgList(mass));
		RooFormulaVar alpha2_v("alpha2_v",alpha2_form,RooArgList(mass));
		RooFormulaVar n2_v("n2_v",n2_form,RooArgList(mass));
	
		RooRealVar meanDCB("meanDCB","meanDCB",meanDCB_v.getVal(),120,130);TString init_red_mean = to_string(meanDCB.getVal());
		RooRealVar sigmaDCB("sigmaDCB","sigmaDCB",sigmaDCB_v.getVal(),0,10);TString init_red_sigma = to_string(sigmaDCB.getVal());
		RooRealVar alphaDCB("alphaDCB","alphaDCB",alphaDCB_v.getVal(),0,10);
		RooRealVar nDCB("nDCB","nDCB",nDCB_v.getVal(),0,10);
		RooRealVar alpha2("alpha2","alpha2",alpha2_v.getVal(),0,10);
		RooRealVar n2("n2","n2",n2_v.getVal(),0,10);
		RooDoubleCB DCB("DCB","DCB",mass4l_d,meanDCB,sigmaDCB,alphaDCB,nDCB,alpha2,n2);

		DCB.plotOn(frame,LineColor(kBlue),LineWidth(2));


	}
	
	TCanvas ccc("c","c",1400,1000);
	frame->GetXaxis()->SetTitle("mass4l(GeV)");
	frame->GetYaxis()->SetTitleOffset(1.75);
	frame->Draw();
	
	TString REFIT;
	if(useREFIT)REFIT="_REFIT";
	if(!useREFIT)REFIT="";
	ccc.SaveAs("/home/chenguan/public_html/LineShape/"+year+"/ggHShape/Parametrization_"+str_fs+REFIT+".png");
	
	ofstream in;
	in.open("/home/chenguan/public_html/LineShape/"+year+"/ggHShape/Parametrization_"+str_fs+REFIT+".txt");
	if(str_fs=="4e"){
		mean_form = "'mean_2' : '"+mean_form+"'";
		sigma_form = "'sigma_2' : '"+sigma_form+"'";
		alpha1_form = "'alpha_2' : '"+alpha1_form+"'";
		n1_form = "'n_2' : '"+n1_form+"'";
		alpha2_form = "'alpha2_2' : '"+alpha2_form+"'";
		n2_form = "'n2_2' : '"+n2_form+"'";
	}
	if(str_fs=="4mu"){
	        mean_form = "'mean_1' : '"+mean_form+"'";
	        sigma_form = "'sigma_1' : '"+sigma_form+"'";
	        alpha1_form = "'alpha_1' : '"+alpha1_form+"'";
	        n1_form = "'n_1' : '"+n1_form+"'";
	        alpha2_form = "'alpha2_1' : '"+alpha2_form+"'";
	        n2_form = "'n2_1' : '"+n2_form+"'";
	}
	if(str_fs=="2e2mu"){
	        mean_form = "'mean_3' : '"+mean_form+"'";
	        sigma_form = "'sigma_3' : '"+sigma_form+"'";
	        alpha1_form = "'alpha_3' : '"+alpha1_form+"'";
	        n1_form = "'n_3' : '"+n1_form+"'";
	        alpha2_form = "'alpha2_3' : '"+alpha2_form+"'";
	        n2_form = "'n2_3' : '"+n2_form+"'";
	}

	in<<mean_form<<endl;
	in<<sigma_form<<endl;
	in<<alpha1_form<<endl;
	in<<n1_form<<endl;
	in<<alpha2_form<<endl;
	in<<n2_form<<endl;
	in.close();
	
}


Double_t* doFirstfit(TString inputpath, Int_t fs1,Int_t fs2, TString mh){
	
	
	TString str_fs;
	if(fs1==1)str_fs = "4mu";
	if(fs1==2)str_fs = "4e";
	if(fs1>2)str_fs = "2e2mu";
	
	RooRealVar mass4l_d("mass4l_d","mass4l_d", 105,140);
	RooDataSet* dataset = doDataset(inputpath, fs1, fs2, mh);
	
	RooRealVar meanDCB("meanDCB","meanDCB",125,120,130);
	RooRealVar sigmaDCB("sigmaDCB","sigmaDCB",1,0,10);
	RooRealVar alphaDCB("alphaDCB","alphaDCB",1,0,10);
	RooRealVar nDCB("nDCB","nDCB",1,0,10);
    RooRealVar alpha2("alpha2","alpha2",1,0,10);
	RooRealVar n2("n2","n2",1,0,10);
	RooDoubleCB DCB("DCB","DCB",mass4l_d,meanDCB,sigmaDCB,alphaDCB,nDCB,alpha2,n2);
	
	 
	RooFitResult* r = DCB.fitTo(*dataset,PrintLevel(-1),SumW2Error(kTRUE),Timer(kTRUE),Save(kTRUE));
	RooPlot *frame=mass4l_d.frame(Bins(100));
    frame->SetTitle(str_fs+"_MH"+mh+".png");/////////////////////////////////////////////////////////////////////
    dataset->plotOn(frame);
    DCB.plotOn(frame);
    DCB.paramOn(frame,Layout(0.1,0.4,0.9));
	
	Double_t chisquare=frame->chiSquare(6);
	TLatex *latex=new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.05);
    latex->SetTextFont(42);
    latex->SetTextAlign(23);
    char chi2[20];
    sprintf(chi2,"%s%1.4f","#chi^{2}/DOF=",chisquare);
   TString REFIT;
   if(useREFIT)REFIT="_REFIT";
   if(!useREFIT)REFIT="";
	TCanvas* c = new TCanvas("c","c",1400,1000);
	c->cd();
	frame->Draw();
	latex->DrawLatex(0.7, 0.8, chi2);
	c->SaveAs("/home/chenguan/public_html/LineShape/"+year+"/Separate_Fits/"+str_fs+"_MH"+mh+REFIT+".png");
	
	//plotCorrelationMatrix
	TH2* corr = new TH2D();
	corr = r->correlationHist();
	corr->SetStats(kFALSE);
	TCanvas* cc = new TCanvas("cc","cc",1400,1000);
	cc->cd();
	corr->Draw("colztext");
	//corr->Draw("text");
	cc->SaveAs("/home/chenguan/public_html/LineShape/"+year+"/Separate_Fits/"+str_fs+"_CorrelationMatrix_"+mh+REFIT+".png");
	
	
	Double_t* results = new Double_t[12];
	
	results[0] = meanDCB.getVal();
	results[1] = sigmaDCB.getVal();
	results[2] = alphaDCB.getVal();
	results[3] = nDCB.getVal();
	results[4] = alpha2.getVal();
	results[5] = n2.getVal();
	
	results[6] = meanDCB.getError();
	results[7] = sigmaDCB.getError();
	results[8] = alphaDCB.getError();
	results[9] = nDCB.getError();
	results[10] = alpha2.getError();
	results[11] = n2.getError();

	for(Int_t i=0; i<6; i++){
		cout<<results[i]<<"  "<<results[i+6]<<endl;
	}
	
//	delete dataset;
//	delete r;
//	delete frame;
//	delete latex;
//	delete c;
//	delete cc;
//	delete corr;
	
	return results;
}

RooDataSet* doDataset(TString inputpath, Int_t fs1, Int_t fs2, TString mh){
	
	Int_t sum;
	//TString str_mh = to_string(mh); ///////////////////////////////////
	
	RooRealVar *dataMCWeight_d = new RooRealVar("dataMCWeight_d","dataMCWeight_d",0.0000001,1000);
	RooRealVar *mass4l_d = new RooRealVar("mass4l_d","mass4l_d",105,140); 
//	RooRealVar *mass4lREFIT_d = new RooRealVar("mass4lREFIT_d","mass4lREFIT_d",105,140);
	RooDataSet *dataset4l_d = new RooDataSet("dataset4l_d","dataset4l_d",RooArgSet(*mass4l_d,*dataMCWeight_d));
	
	TFile* f = new TFile(inputpath+year+"GGH_"+mh+".root");
	TTree* t = (TTree*)f->Get("passedEvents");
	setTreeaddress(t);
	sum=t->GetEntries();
	for(Int_t i=0; i<sum; i++){
		t->GetEntry(i);
		if(passedFullSelection==1&&(finalState==fs1||finalState==fs2)&&mass4l<140&&mass4l>105&&pTL1<100&&pTL2<100&&pTL3<100&&pTL4<100&&abs(etaL1)<2.4&&abs(etaL2)<2.4&&abs(etaL3)<2.4&&abs(etaL4)<2.4&&EventCat==0){
			if(useREFIT)*mass4l_d=mass4lREFIT;
			if(!useREFIT)*mass4l_d=mass4l;
			*dataMCWeight_d=dataMCWeight;
			dataset4l_d->add(RooArgSet(*mass4l_d,*dataMCWeight_d));
		}
	}
	RooDataSet *dataset_w=new RooDataSet(dataset4l_d->GetName(),dataset4l_d->GetTitle(),dataset4l_d,*dataset4l_d->get(),"1","dataMCWeight_d");
//	delete t;
//	delete dataMCWeight_d;
//	delete mass4l_d;
//	delete dataset4l_d;
	
	
	return dataset_w;
		
	
	
	
}

Double_t* doSimultaneousfit(TString inputpath, Int_t fs1, Int_t fs2, Double_t* Parameters){//830-950做同时拟合 955-1092对每一个DCB作图
	
	TString str_fs;
	if(fs1==1)str_fs = "4mu";
	if(fs1==2)str_fs = "4e";
	if(fs1>2)str_fs = "2e2mu";
	Double_t* chi2 = new Double_t[5];
	//make dataset
	RooRealVar mass4l_d("mass4l_d","mass4l_d",105,140);
	
	RooDataSet* dataset_120 = doDataset(inputpath, fs1, fs2, "120");
	RooDataSet* dataset_124 = doDataset(inputpath, fs1, fs2, "124");
	RooDataSet* dataset_125 = doDataset(inputpath, fs1, fs2, "125");
	RooDataSet* dataset_126 = doDataset(inputpath, fs1, fs2, "126");
	RooDataSet* dataset_130 = doDataset(inputpath, fs1, fs2, "130");
	
	RooCategory samples("samples","samples");
	
	samples.defineType("mh_120");
	samples.defineType("mh_124");
	samples.defineType("mh_125");
	samples.defineType("mh_126");
	samples.defineType("mh_130");
	
	RooDataSet comb_dataset("comb_dataset","comb_dataset",RooArgSet(mass4l_d),Index(samples),Import("mh_120",*dataset_120),Import("mh_124",*dataset_124),Import("mh_125",*dataset_125),Import("mh_126",*dataset_126),Import("mh_130",*dataset_130));
	
	
	RooRealVar inter_meanDCB("inter_meanDCB","inter_meanDCB",Parameters[0],Parameters[0]-Parameters[6],Parameters[0]+Parameters[6]);
	RooRealVar inter_sigmaDCB("inter_sigmaDCB","inter_sigmaDCB",Parameters[1],Parameters[1]-Parameters[7],Parameters[1]+Parameters[7]);
	RooRealVar inter_alphaDCB("inter_alphaDCB","inter_alphaDCB",Parameters[2],Parameters[2]-Parameters[8],Parameters[2]+Parameters[8]);
	RooRealVar inter_nDCB("inter_nDCB","inter_nDCB",Parameters[3],Parameters[3]-Parameters[9],Parameters[3]+Parameters[9]);
	RooRealVar inter_alpha2("inter_alpha2","inter_alpha2",Parameters[4],Parameters[4]-Parameters[10],Parameters[4]+Parameters[10]);
	RooRealVar inter_n2("inter_n2","inter_n2",Parameters[5],Parameters[5]-Parameters[11],Parameters[5]+Parameters[11]);
	
	
	
	RooRealVar MH_120("MH_120","MH_120",120);
	//MH_120.setConstant(kTRUE);
	RooRealVar MH_124("MH_124","MH_124",124);
	//MH_124.setConstant(kTRUE);
	RooRealVar MH_125("MH_125","MH_125",125);
	//MH_125.setConstant(kTRUE);
	RooRealVar MH_126("MH_126","MH_126",126);
	//MH_126.setConstant(kTRUE);
	RooRealVar MH_130("MH_130","MH_130",130);
	//MH_130.setConstant(kTRUE);
	
	RooRealVar slope_meanDCB("slope_meanDCB","slope_meanDCB",1,0.9,1.1);
	RooRealVar slope_sigmaDCB("slope_sigmaDCB","slope_sigmaDCB",0,-0.5,0.5);
	RooRealVar slope_alphaDCB("slope_alphaDCB","slope_alphaDCB",0,-0.5,0.5);
	RooRealVar slope_nDCB("slope_nDCB","slope_nDCB",0,-10,10);
	RooRealVar slope_alpha2("slope_alpha2","slope_alpha2",0,-0.5,0.5);
	RooRealVar slope_n2("slope_n2","slope_n2",0,-0.5,0.5);
	
	
	//freeze or float parameters accroding to correlation matrix  
	inter_meanDCB.setConstant(kFALSE);
	inter_sigmaDCB.setConstant(kFALSE);
	inter_alphaDCB.setConstant(kTRUE);
	inter_nDCB.setConstant(kTRUE);
	inter_alpha2.setConstant(kTRUE);
	inter_n2.setConstant(kFALSE);
	
	slope_meanDCB.setConstant(kFALSE);
	slope_sigmaDCB.setConstant(kFALSE);
	slope_alphaDCB.setConstant(kFALSE);
	slope_alpha2.setConstant(kFALSE);
	slope_nDCB.setConstant(kTRUE);
	slope_n2.setConstant(kTRUE);
	///////////////////////////////////////////////
	
	RooFormulaVar meanDCB_120("meanDCB_120","@0+@1*(@2-125)",RooArgList(inter_meanDCB,slope_meanDCB,MH_120));
	RooFormulaVar sigmaDCB_120("sigmaDCB_120","@0+@1*(@2-125)",RooArgList(inter_sigmaDCB,slope_sigmaDCB,MH_120));
	RooFormulaVar alphaDCB_120("alphaDCB_120","@0+@1*(@2-125)",RooArgList(inter_alphaDCB,slope_alphaDCB,MH_120));
	RooFormulaVar nDCB_120("nDCB_120","@0+@1*(@2-125)",RooArgList(inter_nDCB,slope_nDCB,MH_120));
	RooFormulaVar alpha2_120("alpha2_120","@0+@1*(@2-125)",RooArgList(inter_alpha2,slope_alpha2,MH_120));
	RooFormulaVar n2_120("n2_120","@0+@1*(@2-125)",RooArgList(inter_n2,slope_n2,MH_120));
	
	RooFormulaVar meanDCB_124("meanDCB_124","@0+@1*(@2-125)",RooArgList(inter_meanDCB,slope_meanDCB,MH_124));
	RooFormulaVar sigmaDCB_124("sigmaDCB_124","@0+@1*(@2-125)",RooArgList(inter_sigmaDCB,slope_sigmaDCB,MH_124));
	RooFormulaVar alphaDCB_124("alphaDCB_124","@0+@1*(@2-125)",RooArgList(inter_alphaDCB,slope_alphaDCB,MH_124));
	RooFormulaVar nDCB_124("nDCB_124","@0+@1*(@2-125)",RooArgList(inter_nDCB,slope_nDCB,MH_124));
	RooFormulaVar alpha2_124("alpha2_124","@0+@1*(@2-125)",RooArgList(inter_alpha2,slope_alpha2,MH_124));
	RooFormulaVar n2_124("n2_124","@0+@1*(@2-125)",RooArgList(inter_n2,slope_n2,MH_124));
	
	RooFormulaVar meanDCB_125("meanDCB_125","@0+@1*(@2-125)",RooArgList(inter_meanDCB,slope_meanDCB,MH_125));
	RooFormulaVar sigmaDCB_125("sigmaDCB_125","@0+@1*(@2-125)",RooArgList(inter_sigmaDCB,slope_sigmaDCB,MH_125));
	RooFormulaVar alphaDCB_125("alphaDCB_125","@0+@1*(@2-125)",RooArgList(inter_alphaDCB,slope_alphaDCB,MH_125));
	RooFormulaVar nDCB_125("nDCB_125","@0+@1*(@2-125)",RooArgList(inter_nDCB,slope_nDCB,MH_125));
	RooFormulaVar alpha2_125("alpha2_125","@0+@1*(@2-125)",RooArgList(inter_alpha2,slope_alpha2,MH_125));
	RooFormulaVar n2_125("n2_125","@0+@1*(@2-125)",RooArgList(inter_n2,slope_n2,MH_125));
	
	RooFormulaVar meanDCB_126("meanDCB_126","@0+@1*(@2-125)",RooArgList(inter_meanDCB,slope_meanDCB,MH_126));
	RooFormulaVar sigmaDCB_126("sigmaDCB_126","@0+@1*(@2-125)",RooArgList(inter_sigmaDCB,slope_sigmaDCB,MH_126));
	RooFormulaVar alphaDCB_126("alphaDCB_126","@0+@1*(@2-125)",RooArgList(inter_alphaDCB,slope_alphaDCB,MH_126));
	RooFormulaVar nDCB_126("nDCB_126","@0+@1*(@2-125)",RooArgList(inter_nDCB,slope_nDCB,MH_126));
	RooFormulaVar alpha2_126("alpha2_126","@0+@1*(@2-125)",RooArgList(inter_alpha2,slope_alpha2,MH_126));
	RooFormulaVar n2_126("n2_126","@0+@1*(@2-125)",RooArgList(inter_n2,slope_n2,MH_126));
	
	RooFormulaVar meanDCB_130("meanDCB_130","@0+@1*(@2-125)",RooArgList(inter_meanDCB,slope_meanDCB,MH_130));
	RooFormulaVar sigmaDCB_130("sigmaDCB_130","@0+@1*(@2-125)",RooArgList(inter_sigmaDCB,slope_sigmaDCB,MH_130));
	RooFormulaVar alphaDCB_130("alphaDCB_130","@0+@1*(@2-125)",RooArgList(inter_alphaDCB,slope_alphaDCB,MH_130));
	RooFormulaVar nDCB_130("nDCB_130","@0+@1*(@2-125)",RooArgList(inter_nDCB,slope_nDCB,MH_130));
	RooFormulaVar alpha2_130("alpha2_130","@0+@1*(@2-125)",RooArgList(inter_alpha2,slope_alpha2,MH_130));
	RooFormulaVar n2_130("n2_130","@0+@1*(@2-125)",RooArgList(inter_n2,slope_n2,MH_130));
	
	RooDoubleCB DCB_120("DCB_120","DCB_120",mass4l_d,meanDCB_120,sigmaDCB_120,alphaDCB_120,nDCB_120,alpha2_120,n2_120);
	RooDoubleCB DCB_124("DCB_124","DCB_124",mass4l_d,meanDCB_124,sigmaDCB_124,alphaDCB_124,nDCB_124,alpha2_124,n2_124);
	RooDoubleCB DCB_125("DCB_125","DCB_125",mass4l_d,meanDCB_125,sigmaDCB_125,alphaDCB_125,nDCB_125,alpha2_125,n2_125);
	RooDoubleCB DCB_126("DCB_126","DCB_126",mass4l_d,meanDCB_126,sigmaDCB_126,alphaDCB_126,nDCB_126,alpha2_126,n2_126);
	RooDoubleCB DCB_130("DCB_130","DCB_130",mass4l_d,meanDCB_130,sigmaDCB_130,alphaDCB_130,nDCB_130,alpha2_130,n2_130);
	
	
	
	RooSimultaneous comb_DCB("comb_DCB","comb_DCB",samples);
	comb_DCB.addPdf(DCB_120,"mh_120");
	comb_DCB.addPdf(DCB_124,"mh_124");
	comb_DCB.addPdf(DCB_125,"mh_125");
	comb_DCB.addPdf(DCB_126,"mh_126");
	comb_DCB.addPdf(DCB_130,"mh_130");
	
	comb_DCB.fitTo(comb_dataset,PrintLevel(-1),SumW2Error(kTRUE),Timer(kTRUE));
    
	 
	
	//make plots
	for(Int_t i=0; i<6; i++){
		TString str_mh;
		TString str_tag;
		if(i==0){str_tag="120";str_mh="mh_120";}
		if(i==1){str_tag="124";str_mh="mh_124";}
		if(i==2){str_tag="125";str_mh="mh_125";}
		if(i==3){str_tag="126";str_mh="mh_126";}
		if(i==4){str_tag="130";str_mh="mh_130";}
		
		
		RooPlot *frame=mass4l_d.frame(Bins(100));
		frame->SetTitle(str_fs+"_"+str_mh);
		comb_dataset.plotOn(frame,Cut("samples==samples::"+str_mh));
		comb_DCB.plotOn(frame,Slice(samples,str_mh),ProjWData(RooArgSet(samples),comb_dataset));
		//comb_DCB.paramOn(frame);
		Double_t chisquare=frame->chiSquare(6);
		chi2[i] = chisquare;
		TLatex *latex=new TLatex();
		latex->SetNDC();
		latex->SetTextSize(0.05);
		latex->SetTextFont(42);
		latex->SetTextAlign(23);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////a large block to make text about parameters untill 1204
		Double_t slope_mean_v = slope_meanDCB.getVal();
		Double_t slope_sigma_v = slope_sigmaDCB.getVal();
		Double_t slope_alpha1_v = slope_alphaDCB.getVal();
		Double_t slope_n1_v = slope_nDCB.getVal();
		Double_t slope_alpha2_v = slope_alpha2.getVal();
		Double_t slope_n2_v = slope_n2.getVal();
		
		Double_t meanDCB_120_v = meanDCB_120.getVal(); TString str_meanDCB_120 = to_string(meanDCB_120_v);
		Double_t meanDCB_124_v = meanDCB_124.getVal(); TString str_meanDCB_124 = to_string(meanDCB_124_v);
		Double_t meanDCB_125_v = meanDCB_125.getVal(); TString str_meanDCB_125 = to_string(meanDCB_125_v);
		Double_t meanDCB_126_v = meanDCB_126.getVal(); TString str_meanDCB_126 = to_string(meanDCB_126_v);
		Double_t meanDCB_130_v = meanDCB_130.getVal(); TString str_meanDCB_130 = to_string(meanDCB_130_v);
		
		Double_t sigmaDCB_120_v = sigmaDCB_120.getVal(); TString str_sigma_120 = to_string(sigmaDCB_120_v);
		Double_t sigmaDCB_124_v = sigmaDCB_124.getVal(); TString str_sigma_124 = to_string(sigmaDCB_124_v);
		Double_t sigmaDCB_125_v = sigmaDCB_125.getVal(); TString str_sigma_125 = to_string(sigmaDCB_125_v);
		Double_t sigmaDCB_126_v = sigmaDCB_126.getVal(); TString str_sigma_126 = to_string(sigmaDCB_126_v);
		Double_t sigmaDCB_130_v = sigmaDCB_130.getVal(); TString str_sigma_130 = to_string(sigmaDCB_130_v);
		
		Double_t alphaDCB_120_v = alphaDCB_120.getVal(); TString str_alphaDCB_120 = to_string(alphaDCB_120_v);
		Double_t alphaDCB_124_v = alphaDCB_124.getVal(); TString str_alphaDCB_124 = to_string(alphaDCB_124_v);
		Double_t alphaDCB_125_v = alphaDCB_125.getVal(); TString str_alphaDCB_125 = to_string(alphaDCB_125_v);
		Double_t alphaDCB_130_v = alphaDCB_130.getVal(); TString str_alphaDCB_130 = to_string(alphaDCB_130_v);
		Double_t alphaDCB_126_v = alphaDCB_126.getVal(); TString str_alphaDCB_126 = to_string(alphaDCB_126_v);
		
		Double_t nDCB_120_v = nDCB_120.getVal(); TString str_nDCB_120 = to_string(nDCB_120_v);
		Double_t nDCB_124_v = nDCB_124.getVal(); TString str_nDCB_124 = to_string(nDCB_124_v);
		Double_t nDCB_125_v = nDCB_125.getVal(); TString str_nDCB_125 = to_string(nDCB_125_v);
		Double_t nDCB_126_v = nDCB_126.getVal(); TString str_nDCB_126 = to_string(nDCB_126_v);
		Double_t nDCB_130_v = nDCB_130.getVal(); TString str_nDCB_130 = to_string(nDCB_130_v);
		
		Double_t alpha2_120_v = alpha2_120.getVal(); TString str_alpha2_120 = to_string(alpha2_120_v);
		Double_t alpha2_124_v = alpha2_124.getVal(); TString str_alpha2_124 = to_string(alpha2_124_v);
		Double_t alpha2_125_v = alpha2_125.getVal(); TString str_alpha2_125 = to_string(alpha2_125_v);
		Double_t alpha2_126_v = alpha2_126.getVal(); TString str_alpha2_126 = to_string(alpha2_126_v);
		Double_t alpha2_130_v = alpha2_130.getVal(); TString str_alpha2_130 = to_string(alpha2_130_v);
		
		Double_t n2_120_v = n2_120.getVal(); TString str_n2_120 = to_string(n2_120_v);
		Double_t n2_124_v = n2_124.getVal(); TString str_n2_124 = to_string(n2_124_v);
		Double_t n2_125_v = n2_125.getVal(); TString str_n2_125 = to_string(n2_125_v);
		Double_t n2_126_v = n2_126.getVal(); TString str_n2_126 = to_string(n2_126_v);
		Double_t n2_130_v = n2_130.getVal(); TString str_n2_130 = to_string(n2_130_v);
		
		TString str_slope_mean = to_string(slope_mean_v);
		TString str_slope_sigma = to_string(slope_sigma_v);
		TString str_slope_alpha1 = to_string(slope_alpha1_v);
		TString str_slope_n1 = to_string(slope_n1_v);
		TString str_slope_alpha2 = to_string(slope_alpha2_v);
		TString str_slope_n2 = to_string(slope_n2_v);
		TString str_chi2 = to_string(chisquare);
		cout<<"*************!!!!!!!!!!!!!!!!**************"<<chisquare<<endl;
		cout<<"*************!!!!!!!!!!!!!!!!**************"<<str_chi2<<endl;
		TCanvas c("c","c",1400,1000);
		c.cd();
		frame->Draw();
		if(i==0){
		latex->DrawLatex(0.7, 0.8, "#chi^{2}/DOF="+str_chi2);
		latex->DrawLatex(0.7, 0.7, "meanDCB="+str_meanDCB_120);
		latex->DrawLatex(0.7, 0.6, "sigam="+str_sigma_120);
		latex->DrawLatex(0.7, 0.5, "alpha1="+str_alphaDCB_120);
		latex->DrawLatex(0.7, 0.4, "n1="+str_nDCB_120);
		latex->DrawLatex(0.7, 0.3, "alpha2="+str_alpha2_120);
		latex->DrawLatex(0.7, 0.2, "n2="+str_n2_120);
		}
		if(i==1){
		latex->DrawLatex(0.7, 0.8, "#chi^{2}/DOF="+str_chi2);
		latex->DrawLatex(0.7, 0.7, "meanDCB="+str_meanDCB_124);
		latex->DrawLatex(0.7, 0.6, "sigam="+str_sigma_124);
		latex->DrawLatex(0.7, 0.5, "alpha1="+str_alphaDCB_124);
		latex->DrawLatex(0.7, 0.4, "n1="+str_nDCB_124);
		latex->DrawLatex(0.7, 0.3, "alpha2="+str_alpha2_124);
		latex->DrawLatex(0.7, 0.2, "n2="+str_n2_124);
		}
		if(i==2){
		latex->DrawLatex(0.7, 0.8, "#chi^{2}/DOF="+str_chi2);
		latex->DrawLatex(0.7, 0.7, "meanDCB="+str_meanDCB_125);
		latex->DrawLatex(0.7, 0.6, "sigam="+str_sigma_125);
		latex->DrawLatex(0.7, 0.5, "alpha1="+str_alphaDCB_125);
		latex->DrawLatex(0.7, 0.4, "n1="+str_nDCB_125);
		latex->DrawLatex(0.7, 0.3, "alpha2="+str_alpha2_125);
		latex->DrawLatex(0.7, 0.2, "n2="+str_n2_125);
		}
		if(i==3){
		latex->DrawLatex(0.7, 0.8, "#chi^{2}/DOF="+str_chi2);
		latex->DrawLatex(0.7, 0.7, "meanDCB="+str_meanDCB_126);
		latex->DrawLatex(0.7, 0.6, "sigam="+str_sigma_126);
		latex->DrawLatex(0.7, 0.5, "alpha1="+str_alphaDCB_126);
		latex->DrawLatex(0.7, 0.4, "n1="+str_nDCB_126);
		latex->DrawLatex(0.7, 0.3, "alpha2="+str_alpha2_126);
		latex->DrawLatex(0.7, 0.2, "n2="+str_n2_126);
		}
		if(i==4){
		latex->DrawLatex(0.7, 0.8, "#chi^{2}/DOF="+str_chi2);
		latex->DrawLatex(0.7, 0.7, "meanDCB="+str_meanDCB_130);
		latex->DrawLatex(0.7, 0.6, "sigam="+str_sigma_130);
		latex->DrawLatex(0.7, 0.5, "alpha1="+str_alphaDCB_130);
		latex->DrawLatex(0.7, 0.4, "n1="+str_nDCB_130);
		latex->DrawLatex(0.7, 0.3, "alpha2="+str_alpha2_130);
		latex->DrawLatex(0.7, 0.2, "n2="+str_n2_130);
		}
		if(i==5){
		latex->DrawLatex(0.7, 0.8, "#chi^{2}/DOF="+str_chi2);
		latex->DrawLatex(0.7, 0.7, "slope_mean="+str_slope_mean);
		latex->DrawLatex(0.7, 0.6, "slope_sigma="+str_slope_sigma);
		latex->DrawLatex(0.7, 0.5, "sloep_alpha1="+str_slope_alpha1);
		latex->DrawLatex(0.7, 0.4, "slope_n1="+str_slope_n1);
		latex->DrawLatex(0.7, 0.3, "slope_alpha2="+str_slope_alpha2);
		latex->DrawLatex(0.7, 0.2, "slope_n2="+str_slope_n2);
		}//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		TString REFIT;
		if(useREFIT)REFIT="_REFIT";
		if(!useREFIT)REFIT="";
		c.SaveAs("/home/chenguan/public_html/LineShape/"+year+"/ggHShape/"+str_fs+"_"+str_mh+REFIT+".png");
		
//		delete frame;
//		delete c;
//		delete latex;
	}

		
		cout<<"slope_mean = "<<slope_meanDCB.getVal()<<endl;
		cout<<"slope_sigma = "<<slope_sigmaDCB.getVal()<<endl;
		cout<<"slope_alpha1 = "<<slope_alphaDCB.getVal()<<endl;
		cout<<"slope_nDCB1 = "<<slope_nDCB.getVal()<<endl;
		cout<<"slope_alpha2 = "<<slope_alpha2.getVal()<<endl;
		cout<<"slope_n2 = "<<slope_n2.getVal()<<endl;
		
		Double_t* Shape = new Double_t[12];	
		for(Int_t i=0; i<6; i++){
			Shape[i] = Parameters[i];
		}
		Shape[6] = slope_meanDCB.getVal(); Shape[7] = slope_sigmaDCB.getVal();
		Shape[8] = slope_alphaDCB.getVal(); Shape[9] = slope_nDCB.getVal();
		Shape[10] = slope_alpha2.getVal(); Shape[11] = slope_n2.getVal();
		
		
//		delete dataset_120;
//		delete dataset_124;
//		delete dataset_125;
//		delete dataset_126;
//		delete dataset_130;
		//delete chi2;
		return Shape;
	
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

