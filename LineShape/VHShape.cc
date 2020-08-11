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
bool useREFIT = 0;//based on HIG16041 datacard this kind of signal always is useREFIT=0


bool passedFullSelection;
Float_t mass4l;
Float_t mass4lREFIT;
Int_t finalState;
Float_t dataMCWeight;
Int_t EventCat;
Float_t pTL1, pTL2, pTL3, pTL4, etaL1, etaL2, etaL3, etaL4;

void setTreeaddress(TTree* t);
RooDataSet* doDataset(TString inputpath, Int_t fs1, Int_t fs2, TString production);
void doPlots(Int_t fs1, Int_t fs2, RooDataSet* dataset, TString production);


void VHShape(){
	
	TString inputpath="/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/";
	
	TString production = "WH";
	
	for(Int_t i=0; i<3; i++){
		Int_t fs1, fs2;
		if(i==0){fs1=1; fs2=1;}
		if(i==1){fs1=2; fs2=2;}
		if(i==2){fs1=3; fs2=4;}
		RooDataSet* dataset = doDataset(inputpath, fs1, fs2, production);
		doPlots(fs1, fs2, dataset, production);
	}

		





}

void doPlots(Int_t fs1, Int_t fs2, RooDataSet* dataset, TString production){
	
	TString str_fs;
	if(fs1==1)str_fs = "4mu";
	if(fs1==2)str_fs = "4e";
	if(fs1>2)str_fs = "2e2mu";
	
	TString REFIT;
	if(useREFIT)REFIT="_REFIT";
	if(!useREFIT)REFIT="";
	
	TString path = "/home/chenguan/public_html/LineShape/"+year+"/ggHShape/Parametrization_"+str_fs+REFIT+".txt";
	string mean_form, sigma_form, alpha1_form, n1_form, alpha2_form, n2_form,tmp;
	ifstream in;
	in.open(path);
	for(Int_t i=0; i<6; i++){
		if(i==0)getline(in,mean_form);
		if(i==1)getline(in,sigma_form);
		if(i==2)getline(in,alpha1_form);
		if(i==3)getline(in,n1_form);
		if(i==4)getline(in,alpha2_form);
		if(i==5)getline(in,n2_form);
	}
	string mean_tstr(mean_form);
	string sigma_tstr(sigma_form);
	string alpha1_tstr(alpha1_form);
	string n1_tstr(n1_form);
	string alpha2_tstr(alpha2_form);
	string n2_tstr(n2_form);
	
	string mean_tstr_ = mean_tstr.substr(0,mean_tstr.length()-1);
	string sigma_tstr_ = sigma_tstr.substr(0,sigma_tstr.length()-1);
	string alpha1_tstr_ = alpha1_tstr.substr(0,alpha1_tstr.length()-1);
	string n1_tstr_ = n1_tstr.substr(0,n1_tstr.length()-1);
	string alpha2_tstr_ = alpha2_tstr.substr(0,alpha2_tstr.length()-1);
	string n2_tstr_ = n2_tstr.substr(0,n2_tstr.length()-1);

	string mean_str = mean_tstr_.substr(12);
	string sigma_str = sigma_tstr_.substr(13);
	string alpha1_str = alpha1_tstr_.substr(13);
	string n1_str = n1_tstr_.substr(9);
	string alpha2_str = alpha2_tstr_.substr(14);
	string n2_str = n2_tstr_.substr(10);

	TString mean_str_(mean_str);
	TString sigma_str_(sigma_str);
	TString alpha1_str_(alpha1_str);
	TString n1_str_(n1_str);
	TString alpha2_str_(alpha2_str);
	TString n2_str_(n2_str);

/*
	cout<<mean_str<<endl;
	cout<<sigma_str<<endl;
	cout<<alpha1_str<<endl;
	cout<<n1_str<<endl;
	cout<<alpha2_str<<endl;
	cout<<n2_str<<endl;
*/
	RooRealVar mass4l_d("mass4l_d","mass4l_d",105,140);
	RooRealVar mass("mass","mass",120,130);
	mass.setVal(125);
	mass.setConstant(kTRUE);
	
	RooFormulaVar meanDCB_v("meanDCB_v",mean_str_,RooArgList(mass));
	RooFormulaVar sigmaDCB_v("sigmaDCB_v",sigma_str_,RooArgList(mass));
	RooFormulaVar alphaDCB_v("alphaDCB_v",alpha1_str_,RooArgList(mass));
	RooFormulaVar nDCB_v("nDCB_v",n1_str_,RooArgList(mass));
	RooFormulaVar alpha2_v("alpha2_v",alpha2_str_,RooArgList(mass));
	RooFormulaVar n2_v("n2_v",n2_str_,RooArgList(mass));
	
	RooRealVar meanDCB("meanDCB","meanDCB",meanDCB_v.getVal());
	RooRealVar sigmaDCB("sigmaDCB","sigmaDCB",sigmaDCB_v.getVal());
	RooRealVar alphaDCB("alphaDCB","alphaDCB",alphaDCB_v.getVal());
	RooRealVar nDCB("nDCB","nDCB",nDCB_v.getVal());
	RooRealVar alpha2("alpha2","alpha2",alpha2_v.getVal());
	RooRealVar n2("n2","n2",n2_v.getVal());
	TString str_mean = to_string(meanDCB.getVal());
	TString str_sigma = to_string(sigmaDCB.getVal());
	
	RooDoubleCB DCB("DCB","DCB",mass4l_d,meanDCB,sigmaDCB,alphaDCB,nDCB,alpha2,n2);
	
	RooRealVar p0_v("p0_v","p0_v",130,150);
	RooRealVar p1_v("p1_v","p1_v",10,20);
	RooLandau landau("landau", "landau", mass4l_d, p0_v, p1_v);
	RooRealVar frac_v("frac_v","frac_v",0,1);
	RooAddPdf Model("Model", "Model", DCB, landau, frac_v);
	
	Model.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1),Timer(kTRUE));
	
	Double_t p0_, p1_, frac_, p0Err, p1Err, fracErr;
	p0_=p0_v.getVal(); p1_= p1_v.getVal(); frac_=frac_v.getVal(); p0Err=p0_v.getError(); p1Err=p1_v.getError(); fracErr=frac_v.getError();
	
//	cout<<p0_<<"  "<<p1_<<" "<<frac_<<endl;
//	cout<<p0Err<<" "<<p1Err<<" "<<fracErr<<endl;
	
	Double_t p0_dn = p0_-p0Err;
	Double_t p0_up = p0_+p0Err;
	Double_t p1_dn = p1_-p1Err;
	Double_t p1_up = p1_+p1Err;
	Double_t frac_dn = frac_-fracErr;
	Double_t frac_up = frac_+fracErr;
	
	RooRealVar p0("p0","p0",p0_,p0_dn,p0_up);
	RooRealVar p1("p1","p1",p1_,p1_dn,p1_up);
	RooRealVar frac("frac","frac",frac_,frac_dn,frac_up);
	RooLandau landau_("landau_","landau_",mass4l_d, p0,p1);
	RooAddPdf Model_("Model_","Model_", DCB, landau_, frac);
	
	Model_.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1),Timer(kTRUE));

	RooPlot* frame = mass4l_d.frame(Bins(70));
	frame->GetXaxis()->SetTitle("mass4l(GeV)");
	frame->GetYaxis()->SetTitleOffset(1.55);
	frame->SetTitle("");
	dataset->plotOn(frame,MarkerStyle(24));
	DCB.plotOn(frame,LineColor(kGreen));
	landau_.plotOn(frame,LineColor(7));
	Model_.plotOn(frame,LineColor(28));
	Model_.paramOn(frame,Layout(0.1,0.1,0.9));
	
	Double_t chisquare=frame->chiSquare(3);
	TLatex *latex=new TLatex();
	latex->SetNDC();
	latex->SetTextSize(0.04);
	latex->SetTextFont(42);
	latex->SetTextAlign(13);
	char chi2[20];
	sprintf(chi2,"%s%1.4f","#chi^{2}/DOF=",chisquare);
	
	TCanvas c("c","c",1400,1000);
	c.cd();
	frame->Draw();
	latex->DrawLatex(0.1, 0.54, chi2);
	latex->DrawLatex(0.1, 0.7, "mean: "+str_mean);
	latex->DrawLatex(0.1, 0.62, "sigma: "+str_sigma);
	
	string str_p0, str_p1, str_frac;
	str_p0 = to_string(p0.getVal());
	str_p1 = to_string(p1.getVal());
	str_frac = to_string(frac.getVal());
	
	c.SaveAs("/home/chenguan/public_html/LineShape/"+year+"/VHShape/"+production+"_"+str_fs+REFIT+".png");
	
	ofstream out;
	out.open("/home/chenguan/public_html/LineShape/"+year+"/VHShape/"+production+"_"+str_fs+REFIT+".txt");
	if(production=="WH"){
		if(str_fs == "4e"){
			str_p0 = "'WHp0_2' : '"+str_p0+"'";
			str_p1 = "'WHp1_2' : '"+str_p1+"'";
			str_frac = "'WHfrac_2' : '"+str_frac+"'";
		}
		if(str_fs == "4mu"){
			str_p0 = "'WHp0_1' : '"+str_p0+"'";
			str_p1 = "'WHp1_1' : '"+str_p1+"'";
			str_frac = "'WHfrac_1' : '"+str_frac+"'";
		}
		if(str_fs == "2e2mu"){
			str_p0 = "'WHp0_3' : '"+str_p0+"'";
			str_p1 = "'WHp1_3' : '"+str_p1+"'";
			str_frac = "'WHfrac_3' : '"+str_frac+"'";
		}
	}
	if(production=="ZH"){
		if(str_fs == "4e"){
		        str_p0 = "'ZHp0_2' : '"+str_p0+"'";
		        str_p1 = "'ZHp1_2' : '"+str_p1+"'";
		        str_frac = "'ZHfrac_2' : '"+str_frac+"'";
		}
		if(str_fs == "4mu"){
		        str_p0 = "'ZHp0_1' : '"+str_p0+"'";
		        str_p1 = "'ZHp1_1' : '"+str_p1+"'";
		        str_frac = "'ZHfrac_1' : '"+str_frac+"'";
		}
		if(str_fs == "2e2mu"){
		        str_p0 = "'ZHp0_3' : '"+str_p0+"'";
		        str_p1 = "'ZHp1_3' : '"+str_p1+"'";
		        str_frac = "'ZHfrac_3' : '"+str_frac+"'";
		}
	}
	out<<str_p0<<endl;
	out<<str_p1<<endl;
	out<<str_frac<<endl;
	out.close();

	cout<<p0.getVal()<<" "<<p1.getVal()<<" "<<frac.getVal()<<endl;

}



RooDataSet* doDataset(TString inputpath, Int_t fs1, Int_t fs2, TString production){
	
	Int_t sum;
	TString filename;
	if(production=="ZH")filename="ZH_125";
	if(production=="WH")filename="WH_125";
	
	RooRealVar *dataMCWeight_d = new RooRealVar("dataMCWeight_d","dataMCWeight_d",0.0000001,1000);
	RooRealVar *mass4l_d = new RooRealVar("mass4l_d","mass4l_d",105,140); 
	RooDataSet *dataset4l_d = new RooDataSet("dataset4l_d","dataset4l_d",RooArgSet(*mass4l_d,*dataMCWeight_d));
	
	TFile* f = new TFile(inputpath+year+filename+".root");
	TTree* t = (TTree*)f->Get("passedEvents");
	setTreeaddress(t);
	sum=t->GetEntries();
	for(Int_t i=0; i<sum; i++){
		t->GetEntry(i);
		if(passedFullSelection==1&&(finalState==fs1||finalState==fs2)&&mass4l<140&&mass4l>105&&pTL1<100&&pTL2<100&&pTL3<100&&pTL4<100&&abs(etaL1)<2.4&&abs(etaL2)<2.4&&abs(etaL3)<2.4&&abs(etaL4)<2.4){
			if(!useREFIT)*mass4l_d=mass4l;
			if(useREFIT)*mass4l_d=mass4lREFIT;
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
	t->SetBranchAddress("EventCat",&EventCat);
	t->SetBranchAddress("pTL1",&pTL1);
	t->SetBranchAddress("pTL2",&pTL2);
	t->SetBranchAddress("pTL3",&pTL3);
	t->SetBranchAddress("pTL4",&pTL4);
	t->SetBranchAddress("etaL1",&etaL1);
	t->SetBranchAddress("etaL2",&etaL2);
	t->SetBranchAddress("etaL3",&etaL3);
	t->SetBranchAddress("etaL4",&etaL4);

}
