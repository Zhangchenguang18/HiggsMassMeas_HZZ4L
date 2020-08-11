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
using namespace RooFit ;


Double_t* getLambda3(Double_t p1, Double_t p2, Double_t* Lambda_ex, Double_t Lambda2, bool dopara, int interation);//for trecker electron 
Double_t* getLambda2(Double_t p1, Double_t p2, Double_t q1, Double_t q2, Double_t* Lambda_ex, Double_t Lambda2, bool dopara, int interation);//for ECAL electron low statistic 
Double_t getLambda1(Double_t p1, Double_t p2, Double_t q1, Double_t q2);//for ECAL electron high statistic
Double_t getLambda4(Double_t p1, Double_t p2);//for muon 

//TString dataname="MC/2016/";bool isData=false;TString year="2016";
//TString dataname="Data/2016/";bool isData=true;TString year="2016";
//TString dataname="MC/2017/";bool isData=false;TString year="2017";
//TString dataname="Data/2017/";bool isData=true;TString year="2017";
//TString dataname="MC/2018/";bool isData=false;TString year="2018";
TString dataname="Data/2018/";bool isData=true;TString year="2018";


TString plotpathname="/home/chenguan/public_html/getLambda/";

double Zwidth;


double *mean_temp=new double(); 
double *alpha_temp=new double(); 
double *n_temp=new double(); 
double *tau_temp=new double(); 
double *fsig_temp=new double();

Double_t getLambda1(Double_t p1, Double_t p2, Double_t q1, Double_t q2){//for some bins with high statistic just need eta1/2 and pt1/2 electron

if(isData==true)Zwidth=2.49;
if(isData==false)Zwidth=2.44;

//transfer double to string
char p1c[10];
char p2c[10];
sprintf(p1c,"%1.2f",p1);
sprintf(p2c,"%1.2f",p2);
TString p1s;
TString p2s;
p1s=p1c;
p2s=p2c;
char q1c[10];
char q2c[10];
sprintf(q1c,"%1.3f",q1);
sprintf(q2c,"%1.3f",q2);
TString q1s;
TString q2s;
q1s=q1c;
q2s=q2c;

//make string for all filename, title and path 
//title 
TString getparatitle="getPara_eta_"+p1s+"_"+p2s+"_"+"pterr"+"_"+q1s+"_"+q2s;
TString getlambdatitle="getLambda_eta_"+p1s+"_"+p2s+"_"+"pterr"+"_"+q1s+"_"+q2s;
//name
TString getparaname=getparatitle+".png";
TString getlambdaname=getlambdatitle;
//path 
TString getparafullpath=plotpathname+"ecalDriven/"+dataname+getparaname;
TString getlambdafullpath=plotpathname+"ecalDriven/"+dataname+getlambdaname;
//information name and path 
TString inforname=getlambdatitle+".txt";
TString inforfullpath="/raid/raid9/chenguan/Mass/leptonpTErrCorrector/information/ecalDriven/"+dataname+inforname;

//get roodatahist of mass with cut.
TFile *f=new TFile();



/*
if(isData==false&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/DY_2017/DYJetsToLL_M50_m2e.root");//2017MC for e.  
if(isData==true&&year=="2016")f=new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/Data_2016_v1_20170223/DoubleLepton_m2e.root");//2016Data for e, using in AN-16-442.
if(isData==false&&year=="2016")f=new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/DY_2016MC_v1_20170222/DYJetsToLL_M-50_kalman_v4_m2e.root");//2016MC for e, from doAllClosure.py. when do correction on mu, i guess this is the input using in an-16-442. 
if(isData==true&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/DoubleEG_Run2017-17Nov2017_m2e.root");//2017Data for e.
if(isData==false&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/DY_2018_m2e.root");
if(isData==true&&year=="2018")f=new TFile();

//for Legacy 
if(isData==false&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2017_m2e.root");
if(isData==true&&year=="2016")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2016/Electron/SingleElectron_m2e.root");
//if(isData==false&&year=="2016")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2016_m2e.root");
if(isData==false&&year=="2016")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/forDataSimComparison_relIso_lessthan0_1/2016_m2e.root");
if(isData==true&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2017/Electron/SingleElectron_m2e.root");
if(isData==false&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2018_m2e.root");
//if(isData==true&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2018/Electron/EGamma_m2e.root");
if(isData==true&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2018/Muon/Electron_NoScale_m2e.root");
*/

if(isData==0&&year=="2016")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2016DY_m2e.root");
if(isData==1&&year=="2016")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/crab_SingleElectron_2016_vtx_m2e.root");
if(isData==0&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2017DY_m2e.root");
if(isData==1&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/crab_SingleElectron_2017_vtx_m2e.root");
if(isData==0&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2018DY_m2e.root");
if(isData==1&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/crab_SingleElectron_2018_vtx_m2e.root");

TTree *t=(TTree*)f->Get("passedEvents");
Int_t num;

RooRealVar *massZ_d=new RooRealVar("massZ","massZ",60,120);
RooRealVar *massZErr_d=new RooRealVar("massZErr","massZErr",0.2,7.2);
RooRealVar *weight_d=new RooRealVar("weight","weight",0.00001,100);

RooArgSet *argset=new RooArgSet(*massZ_d,*massZErr_d,*weight_d);

RooDataSet *dataset=new RooDataSet("dataset","dataset",*argset);
massZ_d->setBins(1000,"cache");
massZErr_d->setBins(100,"cache");


Double_t massZ1;
Double_t massZErr1;

Double_t weight1;

Double_t pterr1;
Double_t pterr2;

Double_t pterr_rel1;
Double_t pterr_rel2;

Double_t eta1;
Double_t eta2;

Double_t pT1;
Double_t pT2;

Int_t lep1_ecalDriven;
Int_t lep2_ecalDriven;
Double_t RelIso1, RelIso2;

Double_t GENmass2l;
t->SetBranchAddress("GENmass2l",&GENmass2l);
t->SetBranchAddress("massZ",&massZ1);
t->SetBranchAddress("massZErr",&massZErr1);
t->SetBranchAddress("weight",&weight1);
t->SetBranchAddress("eta1",&eta1);
t->SetBranchAddress("eta2",&eta2);
t->SetBranchAddress("pterr1",&pterr1);
t->SetBranchAddress("pterr2",&pterr2);
t->SetBranchAddress("pT1",&pT1);
t->SetBranchAddress("pT2",&pT2);
t->SetBranchAddress("lep1_ecalDriven",&lep1_ecalDriven);
t->SetBranchAddress("lep2_ecalDriven",&lep2_ecalDriven);
t->SetBranchAddress("RelIso1",&RelIso1);
t->SetBranchAddress("RelIso2",&RelIso2);

Int_t index=0;
Int_t num1=(Int_t)t->GetEntries();
for(Int_t i=0; i<num1; i++){
t->GetEntry(i);
if(RelIso1>0.35||RelIso2>0.35)continue;
if(isData==0&&GENmass2l==0)continue;
pterr_rel1=pterr1/pT1;
pterr_rel2=pterr2/pT2;

if(0.2<massZErr1&&massZErr1<7.2&&60<massZ1&&massZ1<120&&lep1_ecalDriven==1&&lep2_ecalDriven==1&&7<pT1&&pT1<100&&7<pT2&&pT2<100){
       	
        if(p1<abs(eta1)&&abs(eta1)<p2&&p1<abs(eta2)&&abs(eta2)<p2&&q1<pterr_rel1&&pterr_rel1<q2&&q1<pterr_rel2&&pterr_rel2<q2){
        
		
        index=index+1;
        
		massZ_d->setVal(massZ1);
		massZErr_d->setVal(massZErr1);
		weight_d->setVal(weight1);
         		
		dataset->add(*argset);	
        }
        
   }
}

RooDataSet *dataset_w=new RooDataSet(dataset->GetName(),dataset->GetTitle(),dataset,*dataset->get(),"1","weight");

RooDataHist *dataset_binned=dataset_w->binnedClone();




num=dataset->numEntries();


//make model to get para.
RooRealVar massZ_P("massZ","massZ",60,120);
RooRealVar massZErr_P("massZErr","massZErr",0.2,7.2);
RooRealVar PDGmZ("PDGmZ","PDGmZ",91.19);
RooRealVar PDGwZ("PDGwZ","PDGwZ",Zwidth);
PDGmZ.setConstant(kTRUE);
PDGwZ.setConstant(kTRUE);
RooBreitWigner PDGBW("PDGBW","PDGBW",massZ_P,PDGmZ,PDGwZ);

//make crystalball
RooRealVar mean("mean","mean", 0, -5, 5);                                                        
RooRealVar alpha("alpha","alpha", 1, 0, 10);
RooRealVar n("n","n", 5, 0, 60);
RooRealVar sigma("sigma", "sigma", 1, 0, 10);//first fit to get parameters of model and then instead sigma of lambda*masserror.
RooCBShape CB("CB","CB", massZ_P, mean, sigma, alpha, n);
RooFFTConvPdf CW("CW","CW",massZ_P,PDGBW,CB);
RooRealVar tau("tau","tau",  -1, 1);
RooExponential bkg("bkg","bkg", massZ_P, tau);
RooRealVar fsig("fsig","signal fraction", 0.7, 0.5, 1);
RooAddPdf model("model","model", CW, bkg, fsig);
//

model.fitTo(*dataset_binned,SumW2Error(kTRUE),PrintLevel(-1),Timer(kTRUE));
RooPlot *frame=massZ_P.frame(Bins(300));
frame->SetTitle("");
dataset_w->plotOn(frame);
model.plotOn(frame,LineColor(2),LineWidth(1));
model.paramOn(frame,Layout(0.1,0.4,0.9));
dataset_w->statOn(frame,Layout(0.1,0.4,0.5));
TCanvas *c=new TCanvas("c","c",1400,1000);
c->cd();
frame->Draw();

Double_t chisquare=frame->chiSquare(6);
TLatex *latex=new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.05);
    latex->SetTextFont(42);
    latex->SetTextAlign(23);
    char chi2[20];
    sprintf(chi2,"%s%1.4f","#chi^{2}/DOF=",chisquare);
    latex->DrawLatex(0.7, 0.8, chi2);


//c->SaveAs(getparafullpath);

//make model to get lambda
RooRealVar massZ_L("massZ","massZ",60,120);
RooRealVar massZErr_L("massZErr","massZErr",0.2,7.2);
RooRealVar PDGmZ1("PDGmZ","PDGmZ",91.19);
RooRealVar PDGwZ1("PDGwZ","PDGwZ",Zwidth);
PDGmZ1.setConstant(kTRUE);
PDGwZ1.setConstant(kTRUE);
RooBreitWigner PDGBW1("PDGBW","PDGBW",massZ_L,PDGmZ1,PDGwZ1);
RooRealVar mean1("mean1","mean1",mean.getValV());
RooRealVar alpha1("alpha1","alpha",alpha.getValV());
RooRealVar n1("n1","n1",n.getValV());
RooRealVar lambda("lambda","lambda", 1, 0.7, 5);


RooFormulaVar sigma1("sigma1","@1*@0", RooArgList(lambda,massZErr_L));
RooCBShape CB1("CB1","CB1", massZ_L, mean1,sigma1,alpha1,n1);
RooFFTConvPdf CW1("CW1","CW1",massZ_L,PDGBW1,CB1);
RooRealVar tau1("tau1","tau1",tau.getValV());
RooExponential bkg1("bkg1","bkg1", massZ_L, tau1);
RooRealVar fsig1("fsig1","signal fraction", fsig.getValV());
RooAddPdf model1("model1","model1", CW1, bkg1, fsig1);
//


model1.fitTo(*dataset_binned,ConditionalObservables(massZErr_L),SumW2Error(kTRUE),PrintLevel(-9999),Timer(kTRUE));//just need the name of variables are same.
RooPlot *frame2=massZ_L.frame(Bins(300));
frame2->SetTitle("");
dataset_binned->plotOn(frame2);
model1.plotOn(frame2,ProjWData(massZErr_L,*dataset_binned),LineColor(2),LineWidth(1));
//model1.paramOn(frame2);
//model.paramOn(frame2,Layout(0.1,0.4,0.9));
//dataset_w->statOn(frame2,Layout(0.1,0.4,0.5));
TCanvas *c1=new TCanvas("c1","c1",1400,1000);
c1->cd();
c1->SetLeftMargin(0.14);
frame2->Draw();
Double_t chisquare1=frame2->chiSquare(1);
TLatex *latex1=new TLatex();
    latex1->SetNDC();
    latex1->SetTextSize(0.05);
    latex1->SetTextFont(42);
    latex1->SetTextAlign(23);
    char chi21[20];
    sprintf(chi21,"%s%1.4f","#chi^{2}/DOF=",chisquare1);
	TString mean_s, alpha_s, n_s, fsig_s, tau_s, sigma_s, lambda_s, num_s;
	mean_s=to_string(mean1.getVal());
	sigma_s=to_string(sigma1.getVal());
	alpha_s=to_string(alpha1.getVal());
	n_s=to_string(n1.getVal());
	fsig_s=to_string(fsig1.getVal());
	tau_s=to_string(tau1.getVal());
	lambda_s=to_string(lambda.getVal());
	num_s=to_string(num);
    latex1->DrawLatex(0.7, 0.9, chi21);
	latex1->DrawLatex(0.7, 0.8, "lambda="+lambda_s);
	latex1->DrawLatex(0.7, 0.7, "mean="+mean_s);
	latex1->DrawLatex(0.7, 0.6, "sigma="+sigma_s);
	latex1->DrawLatex(0.7, 0.5, "alpha="+alpha_s);
	latex1->DrawLatex(0.7, 0.4, "n="+n_s);
	latex1->DrawLatex(0.7, 0.3, "tau="+tau_s);
	latex1->DrawLatex(0.7, 0.2, "Entries="+num_s);




c1->SaveAs(getlambdafullpath+".png");
c1->SaveAs(getlambdafullpath+".pdf");
Double_t results=lambda.getValV();


//ofstream fout;
//fout.open(inforfullpath);
//fout<<lambda.getValV();

cout<<"lambda="<<lambda.getValV()<<endl;

dataset_w->Print();
cout<<index<<endl;
cout<<num<<endl;
//cout<<chisquare<<endl;
//delete f;//if this pointer is deleted, the run will error!! why??
delete t;
delete massZ_d;
delete massZErr_d;
delete weight_d;
delete argset;
delete dataset;
delete dataset_binned;
delete dataset_w;
delete frame;
delete frame2;
delete c;
delete latex;
delete c1;
delete latex1;

return results;
}


Double_t* getLambda2(Double_t p1, Double_t p2, Double_t q1, Double_t q2, Double_t* Lambda_ex, Double_t Lambda2, bool dopara, int v){//**************************************************for electron with low statistic bins  


if(isData==true)Zwidth=2.49;
if(isData==false)Zwidth=2.44;


//transfer double and int to string
char p1c[10];
char p2c[10];
sprintf(p1c,"%1.2f",p1);
sprintf(p2c,"%1.2f",p2);
TString p1s;
TString p2s;
p1s=p1c;
p2s=p2c;
char q1c[10];
char q2c[10];
sprintf(q1c,"%1.3f",q1);
sprintf(q2c,"%1.3f",q2);
TString q1s;
TString q2s;
q1s=q1c;
q2s=q2c;
char interation_c[5];
sprintf(interation_c,"%d",v);
TString interation_s;
interation_s=interation_c;
//make string for all filename, title and path 
//title 
TString getparatitle="getPara_eta_"+p1s+"_"+p2s+"_"+"pterr"+"_"+q1s+"_"+q2s;
TString getlambdatitle="getLambda_eta_"+p1s+"_"+p2s+"_"+"pterr"+"_"+q1s+"_"+q2s+"_"+"v"+interation_s;
//name
TString getparaname=getparatitle+".png";
TString getlambdaname=getlambdatitle;
//path 
TString getparafullpath=plotpathname+"ecalDriven/"+dataname+getparaname;
TString getlambdafullpath=plotpathname+"ecalDriven/"+dataname+getlambdaname;
//information name and path 
TString inforname=getlambdatitle+".txt";
TString inforfullpath="/raid/raid9/chenguan/Mass/leptonpTErrCorrector/information/ecalDriven/"+dataname+inforname;


//get roodatahist of mass with cut.
TFile *f=new TFile();
/*
if(isData==false&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/DY_2017/DYJetsToLL_M50_m2e.root");//2017MC for e. 
if(isData==false&&year=="2016")f=new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/DY_2016MC_v1_20170222/DYJetsToLL_M-50_kalman_v4_m2e.root");//2016MC for e, from doAllClosure.py. when do correction on mu, i guess this is the input using in an-16-442. 
if(isData==true&&year=="2016")f=new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/Data_2016_v1_20170223/DoubleLepton_m2e.root");//16Data
if(isData==true&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/DoubleEG_Run2017-17Nov2017_m2e.root");//2017Data for e.
if(isData==false&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/DY_2018_m2e.root");

if(isData==false&&year=="2017")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2017_m2e.root");
//if(isData==false&&year=="2016")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2016_m2e.root");
if(isData==false&&year=="2016")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/forDataSimComparison_relIso_lessthan0_1/2016_m2e.root");
if(isData==false&&year=="2018")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2018_m2e.root");
if(isData==1&&year=="2016")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2016/Electron/SingleElectron_m2e.root");
if(isData==1&&year=="2017")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2017/Electron/SingleElectron_m2e.root");
//if(isData==1&&year=="2018")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2018/Electron/EGamma_m2e.root");
if(isData==true&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2018/Muon/Electron_NoScale_m2e.root");
*/

if(isData==0&&year=="2016")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2016DY_m2e.root");
if(isData==1&&year=="2016")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/crab_SingleElectron_2016_vtx_m2e.root");
if(isData==0&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2017DY_m2e.root");
if(isData==1&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/crab_SingleElectron_2017_vtx_m2e.root");
if(isData==0&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2018DY_m2e.root");
if(isData==1&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/crab_SingleElectron_2018_vtx_m2e.root");

TTree *t=(TTree*)f->Get("passedEvents");
Int_t num;

RooRealVar *massZ_d=new RooRealVar("massZ","massZ",60,120);
RooRealVar *massZErr_d=new RooRealVar("massZErr","massZErr",0.2,7.2);
RooRealVar *weight_d=new RooRealVar("weight","weight",0.00001,100);

RooArgSet *argset=new RooArgSet(*massZ_d,*massZErr_d,*weight_d);

RooDataSet *dataset=new RooDataSet("dataset","dataset",*argset);

massZ_d->setBins(1000,"cache");
massZErr_d->setBins(100,"cache");

Double_t massZ1;
Double_t massZErr1;

Double_t weight1;

Double_t pterr1;
Double_t pterr2;

Double_t pterr1_rel;
Double_t pterr2_rel;

Double_t eta1;
Double_t eta2;

Double_t pT1;
Double_t pT2;

Double_t phi1;
Double_t phi2;

Double_t m1;
Double_t m2;

Int_t lep1_ecalDriven;
Int_t lep2_ecalDriven;
Double_t RelIso1, RelIso2;
Double_t GENmass2l;
t->SetBranchAddress("GENmass2l",&GENmass2l);
t->SetBranchAddress("massZ",&massZ1);
t->SetBranchAddress("massZErr",&massZErr1);
t->SetBranchAddress("weight",&weight1);
t->SetBranchAddress("eta1",&eta1);
t->SetBranchAddress("eta2",&eta2);
t->SetBranchAddress("pterr1",&pterr1);
t->SetBranchAddress("pterr2",&pterr2);
t->SetBranchAddress("pT1",&pT1);
t->SetBranchAddress("pT2",&pT2);
t->SetBranchAddress("phi1",&phi1);
t->SetBranchAddress("phi2",&phi2);
t->SetBranchAddress("m1",&m1);
t->SetBranchAddress("m2",&m2);
t->SetBranchAddress("lep1_ecalDriven",&lep1_ecalDriven);
t->SetBranchAddress("lep2_ecalDriven",&lep2_ecalDriven);
t->SetBranchAddress("RelIso1",&RelIso1);
t->SetBranchAddress("RelIso2",&RelIso2);

Int_t index=0;
Int_t num1=(Int_t)t->GetEntries();
for(Int_t i=0; i<num1; i++){
t->GetEntry(i);

if(RelIso1>0.35||RelIso2>0.35)continue;
if(isData==0&&GENmass2l==0)continue;
pterr1_rel=pterr1/pT1;
pterr2_rel=pterr2/pT2;

if(0.2<massZErr1&&massZErr1<7.2&&60<massZ1&&massZ1<120&&lep1_ecalDriven==1&&lep2_ecalDriven==1&&7<pT1&&pT1<100&&7<pT2&&pT2<100){
       	
        if((0<abs(eta1)&&abs(eta1)<1.0&&pterr1_rel<0.025)&&q1<pterr2_rel&&pterr2_rel<q2&&p1<abs(eta2)&&abs(eta2)<p2){//set an electron at first bin and second electron at other bin.
        
		
		Double_t Lambda1;
			if(pterr1_rel<0.01&&abs(eta1)<0.8)Lambda1=Lambda_ex[0];
			if(pterr1_rel<0.015&&pterr1_rel>0.01&&abs(eta1)<0.8)Lambda1=Lambda_ex[1];
			if(pterr1_rel<0.025&&pterr1_rel>0.015&&abs(eta1)<0.8)Lambda1=Lambda_ex[2];
			if(pterr1_rel<0.025&&abs(eta1)>0.8&&abs(eta1)<1.0)Lambda1=Lambda_ex[3];
//			if(pterr1_rel<0.025&&abs(eta1)>1.0&&abs(eta1)<1.44)Lambda1=Lambda_ex[4];
//			if(pterr1_rel<0.06&&pterr1_rel>0.025&&abs(eta1)>1.0&&abs(eta1)<1.44)Lambda1=Lambda_ex[5];
			
		//updata massZErr
		TLorentzVector lep1 = TLorentzVector(0,0,0,0);
             TLorentzVector lep2 = TLorentzVector(0,0,0,0);
             lep1.SetPtEtaPhiM(pT1,eta1,phi1,m1);
             lep2.SetPtEtaPhiM(pT2,eta2,phi2,m2);

             TLorentzVector lep1p = TLorentzVector(0,0,0,0);
             TLorentzVector lep2p = TLorentzVector(0,0,0,0);
             lep1p.SetPtEtaPhiM( (pT1)+(pterr1)*Lambda1, eta1, phi1, m1);
             lep2p.SetPtEtaPhiM( (pT2)+(pterr2)*Lambda2, eta2, phi2, m2);

             double dm1 = (lep1p+lep2).M()-(lep1+lep2).M();
             double dm2 = (lep1+lep2p).M()-(lep1+lep2).M();

             double new_massZErr1 = TMath::Sqrt(dm1*dm1+dm2*dm2);
		
		
		index=index+1;
        massZ_d->setVal(massZ1);
		massZErr_d->setVal(new_massZErr1);
		weight_d->setVal(weight1);
        dataset->add(*argset);	
        }
        else if((0<abs(eta2)&&abs(eta2)<1.0&&pterr2_rel<0.025)&&q1<pterr1_rel&&pterr1_rel<q2&&p1<abs(eta1)&&abs(eta1)<p2){
		
		Double_t Lambda1;
			if(pterr2_rel<0.01&&abs(eta2)<0.8)Lambda1=Lambda_ex[0];
			if(pterr2_rel<0.015&&pterr2_rel>0.01&&abs(eta2)<0.8)Lambda1=Lambda_ex[1];
			if(pterr2_rel<0.025&&pterr2_rel>0.015&&abs(eta2)<0.8)Lambda1=Lambda_ex[2];
			if(pterr2_rel<0.025&&abs(eta2)>0.8&&abs(eta2)<1.0)Lambda1=Lambda_ex[3];
//			if(pterr2_rel<0.025&&abs(eta2)>1.0&&abs(eta2)<1.44)Lambda1=Lambda_ex[4];
//			if(pterr2_rel<0.06&&pterr2_rel>0.025&&abs(eta2)>1.0&&abs(eta2)<1.44)Lambda1=Lambda_ex[5];
		
		//updata massZErr
		TLorentzVector lep1 = TLorentzVector(0,0,0,0);
             TLorentzVector lep2 = TLorentzVector(0,0,0,0);
             lep1.SetPtEtaPhiM(pT1,eta1,phi1,m1);
             lep2.SetPtEtaPhiM(pT2,eta2,phi2,m2);

             TLorentzVector lep1p = TLorentzVector(0,0,0,0);
             TLorentzVector lep2p = TLorentzVector(0,0,0,0);
             lep1p.SetPtEtaPhiM( (pT1)+(pterr1)*Lambda2, eta1, phi1, m1);
             lep2p.SetPtEtaPhiM( (pT2)+(pterr2)*Lambda1, eta2, phi2, m2);

             double dm1 = (lep1p+lep2).M()-(lep1+lep2).M();
             double dm2 = (lep1+lep2p).M()-(lep1+lep2).M();

             double new_massZErr1 = TMath::Sqrt(dm1*dm1+dm2*dm2);
		
		
        index=index+1;
        massZ_d->setVal(massZ1);
		massZErr_d->setVal(new_massZErr1);
		weight_d->setVal(weight1);
        dataset->add(*argset);	
		}
   }
   
}

RooDataSet *dataset_w=new RooDataSet(dataset->GetName(),dataset->GetTitle(),dataset,*dataset->get(),"1","weight");

RooDataHist *dataset_binned=dataset_w->binnedClone();


num=dataset->numEntries();







//make model get para massZ_P and massZErr_P
RooRealVar massZ_P("massZ","massZ",60,120);
RooRealVar massZErr_P("massZErr","massZErr",0.2,7.2);
RooRealVar PDGmZ("PDGmZ","PDGmZ",91.19);
RooRealVar PDGwZ("PDGwZ","PDGwZ",Zwidth);
PDGmZ.setConstant(kTRUE);
PDGwZ.setConstant(kTRUE);
RooBreitWigner PDGBW("PDGBW","PDGBW",massZ_P,PDGmZ,PDGwZ);

//make crystalball
RooRealVar mean("mean","mean", 0, -5, 5);                                                        
RooRealVar alpha("alpha","alpha", 1, 0, 10);
RooRealVar n("n","n", 5, 0, 60);
RooRealVar sigma("sigma", "sigma", 1, 0, 10);
RooCBShape CB("CB","CB", massZ_P, mean, sigma, alpha, n);
RooFFTConvPdf CW("CW","CW",massZ_P,PDGBW,CB);
RooRealVar tau("tau","tau",  -1, 1);
RooExponential bkg("bkg","bkg", massZ_P, tau);
RooRealVar fsig("fsig","signal fraction", 0.7, 0.5, 1);
RooAddPdf model("model","model", CW, bkg, fsig);



if(dopara){


model.fitTo(*dataset_binned,SumW2Error(kTRUE),PrintLevel(-1),Timer(kTRUE));
RooPlot *frame=massZ_P.frame(Bins(300));
frame->SetTitle("");
dataset_binned->plotOn(frame);
model.plotOn(frame,LineColor(2),LineWidth(1));
model.paramOn(frame,Layout(0.1,0.4,0.9));
dataset_w->statOn(frame,Layout(0.1,0.4,0.5));

TCanvas *c=new TCanvas("c","c",1400,1000);
c->cd();
frame->Draw();

Double_t chisquare=frame->chiSquare(6);
TLatex *latex=new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.05);
    latex->SetTextFont(42);
    latex->SetTextAlign(23);
    char chi2[20];
    sprintf(chi2,"%s%1.4f","#chi^{2}/DOF=",chisquare);
    latex->DrawLatex(0.7, 0.8, chi2);


//c->SaveAs(getparafullpath);


*mean_temp=mean.getValV();
*alpha_temp=alpha.getValV();
*n_temp=n.getValV();
*tau_temp=tau.getValV();
*fsig_temp=fsig.getValV();
delete c; delete frame; delete latex;
}


cout<<*mean_temp<<endl;
cout<<*alpha_temp<<endl;
cout<<*n_temp<<endl;
cout<<*tau_temp<<endl;
cout<<*fsig_temp<<endl;




//make model to get lambda massZ_L and massZErr_L.
RooRealVar massZ_L("massZ","massZ",60,120);
RooRealVar massZErr_L("massZErr","massZErr",0.2,7.2);
RooRealVar PDGmZ1("PDGmZ","PDGmZ",91.19);
RooRealVar PDGwZ1("PDGwZ","PDGwZ",Zwidth);
PDGmZ1.setConstant(kTRUE);
PDGwZ1.setConstant(kTRUE);
RooBreitWigner PDGBW1("PDGBW","PDGBW",massZ_P,PDGmZ1,PDGwZ1);
RooRealVar mean1("mean1","mean1",*mean_temp);
RooRealVar alpha1("alpha1","alpha",*alpha_temp);
RooRealVar n1("n1","n1",*n_temp);
RooRealVar lambda("lambda","lambda", 1, 0, 10);


RooFormulaVar sigma1("sigma1","@1*@0", RooArgList(lambda,massZErr_L));
RooCBShape CB1("CB1","CB1", massZ_L, mean1,sigma1,alpha1,n1);
RooFFTConvPdf CW1("CW1","CW1",massZ_L,PDGBW1,CB1);
RooRealVar tau1("tau1","tau1",*tau_temp);
RooExponential bkg1("bkg1","bkg1", massZ_L, tau1);
RooRealVar fsig1("fsig1","signal fraction", *fsig_temp);
RooAddPdf model1("model1","model1", CW1, bkg1, fsig1);


model1.fitTo(*dataset_binned,ConditionalObservables(massZErr_L),SumW2Error(kTRUE),PrintLevel(-9999),Timer(kTRUE));
RooPlot *frame2=massZ_L.frame(Bins(300));
frame2->SetTitle("");
dataset_binned->plotOn(frame2);
model1.plotOn(frame2,ProjWData(massZErr_L,*dataset_binned),LineColor(2),LineWidth(1));
//model1.paramOn(frame2,Layout(0.1,0.4,0.9));
//dataset_w->statOn(frame2,Layout(0.1,0.4,0.5));
TCanvas *c1=new TCanvas("c1","c1",1400,1000);
c1->cd();
c1->SetLeftMargin(0.14);
frame2->Draw();
Double_t chisquare1=frame2->chiSquare(1);
TLatex *latex1=new TLatex();
    latex1->SetNDC();
    latex1->SetTextSize(0.05);
    latex1->SetTextFont(42);
    latex1->SetTextAlign(23);
    char chi21[20];
    sprintf(chi21,"%s%1.4f","#chi^{2}/DOF=",chisquare1);
	TString mean_s, alpha_s, n_s, fsig_s, tau_s, sigma_s, lambda_s, num_s;
	mean_s=to_string(mean1.getVal());
	sigma_s=to_string(sigma1.getVal());
	alpha_s=to_string(alpha1.getVal());
	n_s=to_string(n1.getVal());
	fsig_s=to_string(fsig1.getVal());
	tau_s=to_string(tau1.getVal());
	lambda_s=to_string(lambda.getVal());
	num_s=to_string(num);
    latex1->DrawLatex(0.7, 0.9, chi21);
	latex1->DrawLatex(0.7, 0.8, "lambda="+lambda_s);
	latex1->DrawLatex(0.7, 0.7, "mean="+mean_s);
	latex1->DrawLatex(0.7, 0.6, "sigma="+sigma_s);
	latex1->DrawLatex(0.7, 0.5, "alpha="+alpha_s);
	latex1->DrawLatex(0.7, 0.4, "n="+n_s);
	latex1->DrawLatex(0.7, 0.3, "tau="+tau_s);
	latex1->DrawLatex(0.7, 0.2, "Entries="+num_s);
  
c1->SaveAs(getlambdafullpath+".png");
c1->SaveAs(getlambdafullpath+".pdf");


Double_t result=lambda.getValV();

//ofstream fout;
//fout.open(inforfullpath);
//fout<<Lambda2*result;

cout<<"lambda="<<lambda.getValV()<<endl;
cout<<"Lambda2="<<Lambda2*result<<endl;//here Lambda2 will not return, just used to print out 
dataset_w->Print();
cout<<index<<endl;
cout<<num<<endl;

Double_t *results=new Double_t[2];
results[0]=result;
results[1]=Lambda2*result;

delete t;
delete massZ_d;
delete massZErr_d;
delete weight_d;
delete argset;
delete dataset;
delete dataset_binned;
delete dataset_w;
delete frame2;
delete c1;
delete latex1;

return results;

}


Double_t* getLambda3(Double_t p1, Double_t p2, Double_t* Lambda_ex, Double_t Lambda2, bool dopara, int v){//********************************************************************************************for tracker electron.


if(isData==true)Zwidth=2.49;
if(isData==false)Zwidth=2.44;

//transfer double and int to string
char p1c[10];
char p2c[10];
sprintf(p1c,"%1.2f",p1);
sprintf(p2c,"%1.2f",p2);
TString p1s;
TString p2s;
p1s=p1c;
p2s=p2c;

char interation_c[5];
sprintf(interation_c,"%d",v);
TString interation_s;
interation_s=interation_c;
//make string for all filename, title and path 
//title 
TString getparatitle="getPara_eta_"+p1s+"_"+p2s;
TString getlambdatitle="getLambda_eta_"+p1s+"_"+p2s+"_"+"v"+interation_s;
//name
TString getparaname=getparatitle+".png";
TString getlambdaname=getlambdatitle;
//path 
TString getparafullpath=plotpathname+"trackerDriven/"+dataname+getparaname;
TString getlambdafullpath=plotpathname+"trackerDriven/"+dataname+getlambdaname;
//information name and path 
TString inforname=getlambdatitle+".txt";
TString inforfullpath="/raid/raid9/chenguan/Mass/leptonpTErrCorrector/information/trackerDriven/"+dataname+inforname;




//get roodatahist of mass with cut.
TFile *f=new TFile();
/*
if(isData==false&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/DY_2017/DYJetsToLL_M50_m2e.root");//2017MC for e. 
if(isData==false&&year=="2016")f=new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/DY_2016MC_v1_20170222/DYJetsToLL_M-50_kalman_v4_m2e.root");//2016MC for e, from doAllClosure.py. when do correction on mu, i guess this is the input using in an-16-442. 
if(isData==true&&year=="2016")f=new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/Data_2016_v1_20170223/DoubleLepton_m2e.root");//16Data
if(isData==true&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/DoubleEG_Run2017-17Nov2017_m2e.root");//2017Data for e.
if(isData==false&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/DY_2018_m2e.root");
if(isData==true&&year=="2018")f=new TFile();

if(isData==false&&year=="2017")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2017_m2e.root");
//if(isData==false&&year=="2016")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2016_m2e.root");
if(isData==false&&year=="2016")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/forDataSimComparison_relIso_lessthan0_1/2016_m2e.root");
if(isData==false&&year=="2018")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2018_m2e.root");
if(isData==1&&year=="2016")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2016/Electron/SingleElectron_m2e.root");
if(isData==1&&year=="2017")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2017/Electron/SingleElectron_m2e.root");
//if(isData==1&&year=="2018")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2018/Electron/EGamma_m2e.root");
if(isData==true&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2018/Muon/Electron_NoScale_m2e.root");
*/

if(isData==0&&year=="2016")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2016DY_m2e.root");
if(isData==1&&year=="2016")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/crab_SingleElectron_2016_vtx_m2e.root");
if(isData==0&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2017DY_m2e.root");
if(isData==1&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/crab_SingleElectron_2017_vtx_m2e.root");
if(isData==0&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2018DY_m2e.root");
if(isData==1&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/crab_SingleElectron_2018_vtx_m2e.root");


TTree *t=(TTree*)f->Get("passedEvents");
Int_t num;

RooRealVar *massZ_d=new RooRealVar("massZ","massZ",60,120);
RooRealVar *massZErr_d=new RooRealVar("massZErr","massZErr",0.2,7.2);
RooRealVar *weight_d=new RooRealVar("weight","weight",0.00001,100);

RooArgSet *argset=new RooArgSet(*massZ_d,*massZErr_d,*weight_d);

RooDataSet *dataset=new RooDataSet("dataset","dataset",*argset);

massZ_d->setBins(1000,"cache");
massZErr_d->setBins(100,"cache");

Double_t massZ1;
Double_t massZErr1;

Double_t weight1;

Double_t pterr1;
Double_t pterr2;

Double_t pterr_rel1;
Double_t pterr_rel2;

Double_t eta1;
Double_t eta2;

Double_t pT1;
Double_t pT2;

Double_t phi1;
Double_t phi2;

Double_t m1;
Double_t m2;

Int_t lep1_ecalDriven;
Int_t lep2_ecalDriven;
Double_t RelIso1, RelIso2;
Double_t GENmass2l;
t->SetBranchAddress("GENmass2l",&GENmass2l);
t->SetBranchAddress("RelIso1",&RelIso1);
t->SetBranchAddress("RelIso2",&RelIso2);

t->SetBranchAddress("massZ",&massZ1);
t->SetBranchAddress("massZErr",&massZErr1);
t->SetBranchAddress("weight",&weight1);
t->SetBranchAddress("eta1",&eta1);
t->SetBranchAddress("eta2",&eta2);
t->SetBranchAddress("pterr1",&pterr1);
t->SetBranchAddress("pterr2",&pterr2);
t->SetBranchAddress("pT1",&pT1);
t->SetBranchAddress("pT2",&pT2);
t->SetBranchAddress("phi1",&phi1);
t->SetBranchAddress("phi2",&phi2);
t->SetBranchAddress("m1",&m1);
t->SetBranchAddress("m2",&m2);
t->SetBranchAddress("lep1_ecalDriven",&lep1_ecalDriven);
t->SetBranchAddress("lep2_ecalDriven",&lep2_ecalDriven);

Int_t index=0;
Int_t num1=(Int_t)t->GetEntries();
for(Int_t i=0; i<num1; i++){
t->GetEntry(i);

if(RelIso1>0.35||RelIso2>0.35)continue;
if(isData==0&&GENmass2l==0)continue;

Double_t pterr1_rel=pterr1/pT1;
Double_t pterr2_rel=pterr2/pT2;

if(0.2<massZErr1&&massZErr1<7.2&&60<massZ1&&massZ1<120&&7<pT1&&pT1<100&&7<pT2&&pT2<100){
       	
        if((0<abs(eta1)&&abs(eta1)<1.0&&lep1_ecalDriven==1&&pterr1_rel<0.025)&&p1<abs(eta2)&&abs(eta2)<p2&&lep2_ecalDriven==0){
        
		//updata massZErr
		
		
			Double_t Lambda1;
			if(pterr1_rel<0.01&&abs(eta1)<0.8)Lambda1=Lambda_ex[0];
			if(pterr1_rel<0.015&&pterr1_rel>0.01&&abs(eta1)<0.8)Lambda1=Lambda_ex[1];
			if(pterr1_rel<0.025&&pterr1_rel>0.015&&abs(eta1)<0.8)Lambda1=Lambda_ex[2];
			if(pterr1_rel<0.025&&abs(eta1)>0.8&&abs(eta1)<1.0)Lambda1=Lambda_ex[3];
//			if(pterr1_rel<0.025&&abs(eta1)>1.0&&abs(eta1)<1.44)Lambda1=Lambda_ex[4];
//			if(pterr1_rel<0.0&&pterr1_rel>0.023&&abs(eta1)>1.57&&abs(eta1)<2.5)Lambda1=Lambda_ex[5];
			
			 TLorentzVector lep1 = TLorentzVector(0,0,0,0);
             TLorentzVector lep2 = TLorentzVector(0,0,0,0);
             lep1.SetPtEtaPhiM(pT1,eta1,phi1,m1);
             lep2.SetPtEtaPhiM(pT2,eta2,phi2,m2);

             TLorentzVector lep1p = TLorentzVector(0,0,0,0);
             TLorentzVector lep2p = TLorentzVector(0,0,0,0);
             lep1p.SetPtEtaPhiM( (pT1)+(pterr1)*Lambda1, eta1, phi1, m1);
             lep2p.SetPtEtaPhiM( (pT2)+(pterr2)*Lambda2, eta2, phi2, m2);

             double dm1 = (lep1p+lep2).M()-(lep1+lep2).M();
             double dm2 = (lep1+lep2p).M()-(lep1+lep2).M();

             double new_massZErr1 = TMath::Sqrt(dm1*dm1+dm2*dm2);
		
		
			index=index+1;
			massZ_d->setVal(massZ1);
			massZErr_d->setVal(new_massZErr1);
			weight_d->setVal(weight1);
			dataset->add(*argset);
			
        }
        else if((0<abs(eta2)&&abs(eta2)<1.0&&lep2_ecalDriven==1&&pterr2_rel<0.025)&&p1<abs(eta1)&&abs(eta1)<p2&&lep1_ecalDriven==0){
		
		//updata massZErr
		
			Double_t Lambda1;
			if(pterr2_rel<0.01&&abs(eta2)<0.8)Lambda1=Lambda_ex[0];
			if(pterr2_rel<0.015&&pterr2_rel>0.01&&abs(eta2)<0.8)Lambda1=Lambda_ex[1];
			if(pterr2_rel<0.025&&pterr2_rel>0.015&&abs(eta2)<0.8)Lambda1=Lambda_ex[2];
			if(pterr2_rel<0.025&&abs(eta2)>0.8&&abs(eta2)<1.0)Lambda1=Lambda_ex[3];
//			if(pterr2_rel<0.025&&abs(eta2)>1.0&&abs(eta2)<1.44)Lambda1=Lambda_ex[4];
//			if(pterr2_rel<0.06&&pterr2_rel>0.025&&abs(eta2)>1.0&&abs(eta2)<1.44)Lambda1=Lambda_ex[5];
			
			TLorentzVector lep1 = TLorentzVector(0,0,0,0);
             TLorentzVector lep2 = TLorentzVector(0,0,0,0);
             lep1.SetPtEtaPhiM(pT1,eta1,phi1,m1);
             lep2.SetPtEtaPhiM(pT2,eta2,phi2,m2);

             TLorentzVector lep1p = TLorentzVector(0,0,0,0);
             TLorentzVector lep2p = TLorentzVector(0,0,0,0);
             lep1p.SetPtEtaPhiM( (pT1)+(pterr1)*Lambda2, eta1, phi1, m1);
             lep2p.SetPtEtaPhiM( (pT2)+(pterr2)*Lambda1, eta2, phi2, m2);

             double dm1 = (lep1p+lep2).M()-(lep1+lep2).M();
             double dm2 = (lep1+lep2p).M()-(lep1+lep2).M();

             double new_massZErr1 = TMath::Sqrt(dm1*dm1+dm2*dm2);
		
		
			index=index+1;
			massZ_d->setVal(massZ1);
			massZErr_d->setVal(new_massZErr1);
			weight_d->setVal(weight1);
			dataset->add(*argset);	
			
		}
   }
   
}

RooDataSet *dataset_w=new RooDataSet(dataset->GetName(),dataset->GetTitle(),dataset,*dataset->get(),"1","weight");
RooDataHist *dataset_binned=dataset_w->binnedClone();
num=dataset_w->numEntries();

//make model get para massZ_P and massZErr_P
RooRealVar massZ_P("massZ","massZ",60,120);
RooRealVar massZErr_P("massZErr","massZErr",0.2,7.2);
RooRealVar PDGmZ("PDGmZ","PDGmZ",91.19);
RooRealVar PDGwZ("PDGwZ","PDGwZ",Zwidth);
PDGmZ.setConstant(kTRUE);
PDGwZ.setConstant(kTRUE);
RooBreitWigner PDGBW("PDGBW","PDGBW",massZ_P,PDGmZ,PDGwZ);

//make crystalball
RooRealVar mean("mean","mean", 0, -7, 7);                                                        
RooRealVar alpha("alpha","alpha", 1, 0, 10);
RooRealVar n("n","n", 5, 0, 60);
RooRealVar sigma("sigma", "sigma", 1, 0, 10);
RooCBShape CB("CB","CB", massZ_P, mean, sigma, alpha, n);
RooFFTConvPdf CW("CW","CW",massZ_P,PDGBW,CB);
RooRealVar tau("tau","tau",  -1, 1);
RooExponential bkg("bkg","bkg", massZ_P, tau);
RooRealVar fsig("fsig","signal fraction", 0.7, 0.5, 1);
RooAddPdf model("model","model", CW, bkg, fsig);



if(dopara){


if(num>5000)model.fitTo(*dataset_binned,SumW2Error(kTRUE),PrintLevel(-1),Timer(kTRUE));
if(num<5000)model.fitTo(*dataset_w,SumW2Error(kTRUE),PrintLevel(-1),Timer(kTRUE));
RooPlot *frame=massZ_P.frame(Bins(300));
frame->SetTitle("");
dataset_binned->plotOn(frame);
model.plotOn(frame,LineColor(2),LineWidth(1));
model.paramOn(frame,Layout(0.1,0.4,0.9));
dataset_w->statOn(frame,Layout(0.1,0.4,0.5));

TCanvas *c=new TCanvas("c","c",1400,1000);
c->cd();
frame->Draw();
Double_t chisquare=frame->chiSquare(6);
TLatex *latex=new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.05);
    latex->SetTextFont(42);
    latex->SetTextAlign(23);
    char chi2[20];
    sprintf(chi2,"%s%1.4f","#chi^{2}/DOF=",chisquare);
    latex->DrawLatex(0.7, 0.8, chi2);
//c->SaveAs(getparafullpath);


*mean_temp=mean.getValV();
*alpha_temp=alpha.getValV();
*n_temp=n.getValV();
*tau_temp=tau.getValV();
*fsig_temp=fsig.getValV();
delete c; delete frame; delete latex;
}


cout<<*mean_temp<<endl;
cout<<*alpha_temp<<endl;
cout<<*n_temp<<endl;
cout<<*tau_temp<<endl;
cout<<*fsig_temp<<endl;




//make model to get lambda massZ_L and massZErr_L.
RooRealVar massZ_L("massZ","massZ",60,120);
RooRealVar massZErr_L("massZErr","massZErr",0.2,7.2);
RooRealVar PDGmZ1("PDGmZ","PDGmZ",91.19);
RooRealVar PDGwZ1("PDGwZ","PDGwZ",Zwidth);
PDGmZ1.setConstant(kTRUE);
PDGwZ1.setConstant(kTRUE);
RooBreitWigner PDGBW1("PDGBW","PDGBW",massZ_P,PDGmZ1,PDGwZ1);
RooRealVar mean1("mean1","mean1",*mean_temp);
RooRealVar alpha1("alpha1","alpha",*alpha_temp);
RooRealVar n1("n1","n1",*n_temp);
RooRealVar lambda("lambda","lambda", 1, 0, 10);


RooFormulaVar sigma1("sigma1","@1*@0", RooArgList(lambda,massZErr_L));
RooCBShape CB1("CB1","CB1", massZ_L, mean1,sigma1,alpha1,n1);
RooFFTConvPdf CW1("CW1","CW1",massZ_L,PDGBW1,CB1);
RooRealVar tau1("tau1","tau1",*tau_temp);
RooExponential bkg1("bkg1","bkg1", massZ_L, tau1);
RooRealVar fsig1("fsig1","signal fraction", *fsig_temp);
RooAddPdf model1("model1","model1", CW1, bkg1, fsig1);
	

if(num>5000)model1.fitTo(*dataset_binned,ConditionalObservables(massZErr_L),SumW2Error(kTRUE),PrintLevel(-9999),Timer(kTRUE));
if(num<5000)model1.fitTo(*dataset_w,ConditionalObservables(massZErr_L),SumW2Error(kTRUE),PrintLevel(-9999),Timer(kTRUE));
RooPlot *frame2=massZ_L.frame(Bins(300));
frame2->SetTitle("");
dataset_binned->plotOn(frame2);
model1.plotOn(frame2,ProjWData(massZErr_L,*dataset_binned),LineColor(2),LineWidth(1));
//model1.paramOn(frame2,Layout(0.1,0.4,0.9));
//dataset_w->statOn(frame2,Layout(0.1,0.4,0.5));
TCanvas *c1=new TCanvas("c1","c1",1400,1000);
c1->cd();
c1->SetLeftMargin(0.14);
frame2->Draw();

Double_t chisquare1=frame2->chiSquare(1);
TLatex *latex1=new TLatex();
    latex1->SetNDC();
    latex1->SetTextSize(0.05);
    latex1->SetTextFont(42);
    latex1->SetTextAlign(23);
    char chi21[20];
    sprintf(chi21,"%s%1.4f","#chi^{2}/DOF=",chisquare1);
	TString mean_s, alpha_s, n_s, fsig_s, tau_s, sigma_s, lambda_s, num_s;
	mean_s=to_string(mean1.getVal());
	sigma_s=to_string(sigma1.getVal());
	alpha_s=to_string(alpha1.getVal());
	n_s=to_string(n1.getVal());
	fsig_s=to_string(fsig1.getVal());
	tau_s=to_string(tau1.getVal());
	lambda_s=to_string(lambda.getVal());
	num_s=to_string(num);
    latex1->DrawLatex(0.7, 0.9, chi21);
	latex1->DrawLatex(0.7, 0.8, "lambda="+lambda_s);
	latex1->DrawLatex(0.7, 0.7, "mean="+mean_s);
	latex1->DrawLatex(0.7, 0.6, "sigma="+sigma_s);
	latex1->DrawLatex(0.7, 0.5, "alpha="+alpha_s);
	latex1->DrawLatex(0.7, 0.4, "n="+n_s);
	latex1->DrawLatex(0.7, 0.3, "tau="+tau_s);
	latex1->DrawLatex(0.7, 0.2, "Entries="+num_s);
    
c1->SaveAs(getlambdafullpath+".png");
c1->SaveAs(getlambdafullpath+".pdf");


Double_t result=lambda.getValV();

//ofstream fout;
//fout.open(inforfullpath);
//fout<<Lambda2*result;

cout<<"lambda="<<lambda.getValV()<<endl;
cout<<"Lambda2="<<Lambda2*result<<endl;
dataset_w->Print();
cout<<index<<endl;
cout<<num<<endl;

Double_t *results=new Double_t[2];
results[0]=result;
results[1]=Lambda2*result;

delete t;
delete massZ_d;
delete massZErr_d;
delete weight_d;
delete argset;
delete dataset;
delete dataset_binned;
delete dataset_w;
delete frame2;
delete c1;
delete latex1;


return  results;


}


Double_t getLambda4(Double_t p1, Double_t p2){//********************************************************************************************************for muon 
	
	
if(isData==true)Zwidth=2.49;
if(isData==false)Zwidth=2.44;
	
//transfer double to string
char p1c[10];
char p2c[10];
sprintf(p1c,"%1.2f",p1);
sprintf(p2c,"%1.2f",p2);
TString p1s;
TString p2s;
p1s=p1c;
p2s=p2c;
//make string for all filename, title and path 
//title 
TString getparatitle="getPara_eta_"+p1s+"_"+p2s;
TString getlambdatitle="getLambda_eta_"+p1s+"_"+p2s;
//name
TString getparaname=getparatitle+".png";
TString getlambdaname=getlambdatitle;
//path 
TString getparafullpath=plotpathname+"muon/"+dataname+getparaname;
TString getlambdafullpath=plotpathname+"muon/"+dataname+getlambdaname;
//information name and path 
TString inforname=getlambdatitle+".txt";
TString inforfullpath="/raid/raid9/chenguan/Mass/leptonpTErrCorrector/information/muon/"+dataname+inforname;


TFile *f=new TFile();
/*
if(isData==true&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/DoubleMuon_Run2017-17Nov2017_m2mu.root");//17Data
if(isData==false&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/DY_2017/DYJetsToLL_M50_m2mu.root");//17MC
if(isData==true&&year=="2016")f=new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/Data_2016_v1_20170223/DoubleLepton_m2mu.root");//16Data   //this input comes from huanlin's code doAllClosure_mZ.py. the result are 1.25468, 1.09411, 1.20506.
if(isData==false&&year=="2016")f=new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/DY_2016MC_v1_20170222/DYJetsToLL_M-50_kalman_v4_m2mu.root");//16MC //this input comes from code doAllClosure_mZ.py. //the result are 1.27743, 1.28987, 1.13031.   ！！严重怀疑这个就是AN-16-442使用的input ！！
if(isData==false&&year=="2016")f=new TFile("/raid/raid9/chenguan/input/DY_2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_m2mu.root");//this line aims to see if this sample is same with the one used by hualin for 2016
if(isData==false&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/DY_2018_m2mu.root");
if(isData==true&&year=="2018")f=new TFile();

if(isData==false&&year=="2017")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2017_m2mu.root");
//if(isData==false&&year=="2016")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2016_m2mu.root");
if(isData==false&&year=="2016")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/forDataSimComparison_relIso_lessthan0_1/2016_m2mu.root");
if(isData==false&&year=="2018")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2018_m2mu.root");
if(isData==true&&year=="2016")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2016/Muon/SingleMuon_m2mu.root");
if(isData==true&&year=="2017")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2017/Muon/SingleMuon_m2mu.root");
if(isData==true&&year=="2018")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2018/Muon/SingleMuon_m2mu.root");
*/
if(isData==0&&year=="2016")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2016DY_m2mu.root");
if(isData==1&&year=="2016")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/SingleMuon2016_m2mu.root");
if(isData==0&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2017DY_m2mu.root");
if(isData==1&&year=="2017")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/SingleMuon2017_m2mu.root");
if(isData==0&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/2018DY_m2mu.root");
if(isData==1&&year=="2018")f=new TFile("/raid/raid9/chenguan/input/Vertex_AfterSync_200424/SingleMuon2018_m2mu.root");


TTree *t=(TTree*)f->Get("passedEvents");
Int_t num;

RooRealVar *massZ_d=new RooRealVar("massZ","massZ",60,120);
RooRealVar *massZErr_d=new RooRealVar("massZErr","massZErr",0.2,7.2);
RooRealVar *weight_d=new RooRealVar("weight","weight",0.00001,100);

RooArgSet *argset=new RooArgSet(*massZ_d,*massZErr_d,*weight_d);

RooDataSet *dataset=new RooDataSet("dataset","dataset",*argset);

massZ_d->setBins(1000,"cache");
massZErr_d->setBins(100,"cache");



Double_t massZ1;
Double_t massZErr1;
Double_t weight1;
Double_t eta1;
Double_t eta2;
Double_t pt1;
Double_t pt2;
Double_t RelIso1, RelIso2;
Double_t GENmass2l;
t->SetBranchAddress("GENmass2l",&GENmass2l);
t->SetBranchAddress("RelIso1",&RelIso1);
t->SetBranchAddress("RelIso2",&RelIso2);

t->SetBranchAddress("weight",&weight1);
t->SetBranchAddress("massZ",&massZ1);
t->SetBranchAddress("massZErr",&massZErr1);
t->SetBranchAddress("eta1",&eta1);
t->SetBranchAddress("eta2",&eta2);
t->SetBranchAddress("pT1",&pt1);
t->SetBranchAddress("pT2",&pt2);
Int_t index=0;
Int_t num1=(Int_t)t->GetEntries();
for(Int_t i=0; i<num1; i++){
t->GetEntry(i);

if(RelIso1>0.35||RelIso2>0.35)continue;
if(isData==0&&GENmass2l==0)continue;
if(0.2<massZErr1&&massZErr1<7.2&&60<massZ1&&massZ1<120&&5<pt1&&pt1<100&&5<pt2&&pt2<100){//check the cut condition.
       	
        if(p1<fabs(eta1)&&fabs(eta1)<p2&&p1<fabs(eta2)&&fabs(eta2)<p2){
        
	    
        index=index+1;
        
		massZ_d->setVal(massZ1);
		massZErr_d->setVal(massZErr1);
		weight_d->setVal(weight1);
        
	
		dataset->add(*argset);	
		
        }
		else{continue;}
        
   }
   else{continue;}
}


RooDataSet *dataset_w=new RooDataSet(dataset->GetName(),dataset->GetTitle(),dataset,*dataset->get(),"1","weight");




RooDataHist *dataset_binned=dataset_w->binnedClone();

num=dataset->numEntries();




//make model to getPara

RooRealVar massZ_P("massZ","massZ",60,120);
RooRealVar massZErr_P("massZErr","massZErr",0.2,7.2);
RooRealVar PDGmZ("PDGmZ","PDGmZ",91.19);
RooRealVar PDGwZ("PDGwZ","PDGwZ",Zwidth);
PDGmZ.setConstant(kTRUE);
PDGwZ.setConstant(kTRUE);
RooBreitWigner PDGBW("PDGBW","PDGBW",massZ_P,PDGmZ,PDGwZ);

//make crystalball
RooRealVar mean("mean","mean", 0, -5, 5);                                                        
RooRealVar alpha("alpha","alpha", 1, 0, 10);
RooRealVar n("n","n", 5, 0, 60);
RooRealVar sigma("sigma", "sigma", 1, 0, 10);
RooCBShape CB("CB","CB", massZ_P, mean, sigma, alpha, n);
RooFFTConvPdf CW("CW","CW",massZ_P,PDGBW,CB);
RooRealVar tau("tau","tau",  -1, 1);
RooExponential bkg("bkg","bkg", massZ_P, tau);
RooRealVar fsig("fsig","signal fraction", 0.7, 0.5, 1);
RooAddPdf model("model","model", CW, bkg, fsig);

model.fitTo(*dataset_binned,SumW2Error(kTRUE),PrintLevel(-1),Timer(kTRUE));
RooPlot *frame=massZ_P.frame(Bins(300));
frame->SetTitle("");//
dataset_w->plotOn(frame);
model.plotOn(frame,LineColor(2),LineWidth(1));
//model.plotOn(frame,Components("bkg"),LineStyle(kDashed));
model.paramOn(frame,Layout(0.1,0.4,0.9));
dataset_w->statOn(frame,Layout(0.1,0.4,0.5));
TCanvas *c=new TCanvas("c","c",1400,1000);
c->cd();
frame->Draw();
Double_t chisquare=frame->chiSquare(6);
TLatex *latex=new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.05);
    latex->SetTextFont(42);
    latex->SetTextAlign(23);
    char chi2[20];
    sprintf(chi2,"%s%1.4f","#chi^{2}/DOF=",chisquare);
    latex->DrawLatex(0.7, 0.8, chi2);
//c->SaveAs(getparafullpath);

//make model to get lambda
RooRealVar massZ_L("massZ","massZ",60,120);
RooRealVar massZErr_L("massZErr","massZErr",0.2,7.2);

RooRealVar PDGmZ1("PDGmZ","PDGmZ",91.19);
RooRealVar PDGwZ1("PDGwZ","PDGwZ",Zwidth);
PDGmZ1.setConstant(kTRUE);
PDGwZ1.setConstant(kTRUE);

RooBreitWigner PDGBW1("PDGBW","PDGBW",massZ_L,PDGmZ1,PDGwZ1);
RooRealVar mean1("mean1","mean1",mean.getValV());
RooRealVar alpha1("alpha1","alpha",alpha.getValV());
RooRealVar n1("n1","n1",n.getValV());
RooRealVar lambda("lambda","lambda", 1, 0.7, 5);


RooFormulaVar sigma1("sigma1","@1*@0", RooArgList(lambda,massZErr_L));
RooCBShape CB1("CB1","CB1", massZ_L, mean1,sigma1,alpha1,n1);
RooFFTConvPdf CW1("CW1","CW1",massZ_L,PDGBW1,CB1);
RooRealVar tau1("tau1","tau1",tau.getValV());
RooExponential bkg1("bkg1","bkg1", massZ_L, tau1);
RooRealVar fsig1("fsig1","signal fraction", fsig.getValV());
RooAddPdf model1("model1","model1", CW1, bkg1, fsig1);


model1.fitTo(*dataset_binned,ConditionalObservables(massZErr_L),SumW2Error(kTRUE),PrintLevel(-9999),Timer(kTRUE));
RooPlot *frame2=massZ_L.frame(Bins(300));
frame2->SetTitle("");//
dataset_binned->plotOn(frame2);
model1.plotOn(frame2,ProjWData(massZErr_L,*dataset_binned),LineColor(2),LineWidth(1));
//model1.plotOn(frame2,Components("bkg1"),LineStyle(kDashed));
//model1.paramOn(frame2,Layout(0.1,0.4,0.9));
TCanvas *c1=new TCanvas("c1","c1",1400,1000);
c1->cd();
c1->SetLeftMargin(0.14);
frame2->Draw();
Double_t chisquare1=frame2->chiSquare(1);
TLatex *latex1=new TLatex();
    latex1->SetNDC();
    latex1->SetTextSize(0.05);
    latex1->SetTextFont(42);
    latex1->SetTextAlign(23);
    char chi21[20];
    sprintf(chi21,"%s%1.4f","#chi^{2}/DOF=",chisquare1);
	TString mean_s, alpha_s, n_s, fsig_s, tau_s, sigma_s, lambda_s, num_s;
	mean_s=to_string(mean.getVal());
	sigma_s=to_string(sigma1.getVal());
	alpha_s=to_string(alpha.getVal());
	n_s=to_string(n.getVal());
	fsig_s=to_string(fsig.getVal());
	tau_s=to_string(tau.getVal());
	lambda_s=to_string(lambda.getVal());
	num_s=to_string(num);
    latex1->DrawLatex(0.7, 0.9, chi21);
	latex1->DrawLatex(0.7, 0.8, "lambda="+lambda_s);
	latex1->DrawLatex(0.7, 0.7, "mean="+mean_s);
	latex1->DrawLatex(0.7, 0.6, "sigma="+sigma_s);
	latex1->DrawLatex(0.7, 0.5, "alpha="+alpha_s);
	latex1->DrawLatex(0.7, 0.4, "n="+n_s);
	latex1->DrawLatex(0.7, 0.3, "tau="+tau_s);
	latex1->DrawLatex(0.7, 0.2, "Entries="+num_s);
	
c1->SaveAs(getlambdafullpath+".png");
c1->SaveAs(getlambdafullpath+".pdf");

Double_t results=lambda.getValV();

//ofstream fout;
//fout.open(inforfullpath);
//fout<<lambda.getValV();

cout<<"lambda="<<lambda.getValV()<<endl;

dataset_w->Print();
dataset_binned->Print();
cout<<index<<endl;
cout<<num<<endl;

//return results;
delete t;
delete massZ_d;
delete massZErr_d;
delete weight_d;
delete argset;
delete dataset;
delete dataset_binned;
delete dataset_w;
delete frame;
delete frame2;
delete c;
delete latex;
delete c1;
delete latex1;

return results;
}










void getAllLambda(){

Double_t lambda_1,lambda_2;
Double_t lambda_3,lambda_4,lambda_5,lambda_6,lambda_7;
Double_t lambda_8,lambda_9;//
Double_t lambda_10,lambda_11,lambda_12,lambda_13;//for tracherDriven electron
Double_t lambda_14,lambda_15,lambda_16;//for muon

Double_t lambda_ex1, lambda_ex2, lambda_ex3, lambda_ex4, lambda_ex5;

Double_t *Lp;


int v=1;
Double_t Lambda2=1;


lambda_14=getLambda4(0,0.9);
lambda_15=getLambda4(0.9,1.8);
lambda_16=getLambda4(1.8,2.4);

lambda_1=getLambda1(0,0.8,0.0,0.01);
lambda_ex1=getLambda1(0,0.8,0.01,0.015);
lambda_ex2=getLambda1(0,0.8,0.015,0.025);


lambda_2=getLambda1(0.8,1.0,0,0.025);


lambda_3=getLambda1(1.0,2.5,0,0.02);
lambda_4=getLambda1(1.0,2.5,0.02,0.03);
lambda_5=getLambda1(1.0,2.5,0.03,0.04);
lambda_6=getLambda1(1.0,2.5,0.04,0.06);


Double_t* Lambda1 = new Double_t[4];
Lambda1[0] = lambda_1;
Lambda1[1] = lambda_ex1;
Lambda1[2] = lambda_ex2;
Lambda1[3] = lambda_2;

v=1;
Lambda2=1;
Lp=getLambda2(0.0,1.0,0.025,1,Lambda1,Lambda2,true,0);
for(Int_t interation=0;interation<5;interation++){
	if(abs(Lp[0]-1)>0.01){
		Lambda2=Lambda2*Lp[0];
	    Lp=getLambda2(0.0,1.0,0.025,1,Lambda1,Lambda2,false,v);
	    v=v+1;
		lambda_8=Lp[1];
	}
	else{break;}
}

Lambda2=1;v=1;
Lp=getLambda2(1.0,2.5,0.06,1,Lambda1,Lambda2,true,0);
for(Int_t interation=0;interation<5;interation++){
	if(abs(Lp[0]-1)>0.01){
		Lambda2=Lambda2*Lp[0];
	    Lp=getLambda2(1.0,2.5,0.06,1,Lambda1,Lambda2,false,v);
		v=v+1;
		lambda_9=Lp[1];
	}
	else{break;}
}

/////////////////////////////////////////////do trackerDriven electrons.

//Double_t* Lambda1_ = new Double_t[4];
//Lambda1_[0] = lambda_1;
//Lambda1_[1] = lambda_ex1;
//Lambda1_[2] = lambda_ex2;
//Lambda1_[3] = lambda_8;

Lambda2=1; v=1;
Lp=getLambda3(0,1.44,Lambda1,Lambda2,true,0);
for(Int_t interation=0;interation<5;interation++){
	if(abs(Lp[0]-1)>0.01){
		Lambda2=Lambda2*Lp[0];
	    Lp=getLambda3(0,1.44,Lambda1,Lambda2,false,v);
	    v=v+1;
		lambda_10=Lp[1];
	}
	else{break;}
} 

v=1;
Lambda2=1;
Lp=getLambda3(1.44,1.6,Lambda1,Lambda2,true,0);
for(Int_t interation=0;interation<5;interation++){
	if(abs(Lp[0]-1)>0.01){
		Lambda2=Lambda2*Lp[0];
	    Lp=getLambda3(1.44,1.6,Lambda1,Lambda2,false,v);
	    v=v+1;
		lambda_11=Lp[1];
	}
	else{break;}
} 

v=1;
Lambda2=1;
Lp=getLambda3(1.6,2,Lambda1,Lambda2,true,0);
for(Int_t interation=0;interation<5;interation++){
	if(abs(Lp[0]-1)>0.01){
		Lambda2=Lambda2*Lp[0];
	    Lp=getLambda3(1.6,2,Lambda1,Lambda2,false,v);
	    v=v+1;
		lambda_12=Lp[1];
	}
	else{break;}
} 

v=1;
Lambda2=1;
Lp=getLambda3(2,2.5,Lambda1,Lambda2,true,0);
for(Int_t interation=0;interation<5;interation++){
	if(abs(Lp[0]-1)>0.01){
		Lambda2=Lambda2*Lp[0];
	    Lp=getLambda3(2,2.5,Lambda1,Lambda2,false,v);
	    v=v+1;
		lambda_13=Lp[1];
	}
	else{break;}
} 


/////////make LUT 
	const Int_t mupt=1;  const Int_t ept=1;
	const Int_t mueta=3; const Int_t e3eta=4;
	
	Double_t muptEdges[mupt+1]={5,100}; Double_t eptEdges[ept+1]={7,100};
	Double_t muetaEdges[mueta+1]={0,0.9,1.8,2.4};
	Double_t e3etaEdges[e3eta+1]={0,1.44,1.6,2,2.5};//for trackerDriven electron 
	
	TH2F *mu=new TH2F("mu","mu",mupt,muptEdges,mueta,muetaEdges);
	TH2F *e3=new TH2F("e3","e3",ept,eptEdges,e3eta,e3etaEdges);

	mu->SetBinContent(1,1,lambda_14);
	mu->SetBinContent(1,2,lambda_15);
	mu->SetBinContent(1,3,lambda_16);
	
	e3->SetBinContent(1,1,lambda_10);
	e3->SetBinContent(1,2,lambda_11);
	e3->SetBinContent(1,3,lambda_12);
	e3->SetBinContent(1,4,lambda_13);
	
	//new LUT for ecaldriven electron
	const Int_t ecalDriveneta=3; const Int_t ecalDrivenpterr=8;
	Double_t ecalDrivenetaEdges[ecalDriveneta+1]={0, 0.8, 1.0, 2.5};
	Double_t ecalDrivenpterrEdges[ecalDrivenpterr+1]={0, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.06, 1};
	TH2F* e1 = new TH2F("e1","e1", ecalDriveneta, ecalDrivenetaEdges, ecalDrivenpterr, ecalDrivenpterrEdges);
	
	e1->SetBinContent(1,1,lambda_1);
	e1->SetBinContent(1,2,lambda_ex1);
	e1->SetBinContent(1,3,lambda_ex2);
	e1->SetBinContent(1,4,lambda_ex2);
	e1->SetBinContent(1,5,lambda_8);
	e1->SetBinContent(1,6,lambda_8);
	e1->SetBinContent(1,7,lambda_8);
	e1->SetBinContent(1,8,lambda_8);

	e1->SetBinContent(2,1,lambda_2);
	e1->SetBinContent(2,2,lambda_2);
	e1->SetBinContent(2,3,lambda_2);
	e1->SetBinContent(2,4,lambda_2);
	e1->SetBinContent(2,5,lambda_8);
	e1->SetBinContent(2,6,lambda_8);
	e1->SetBinContent(2,7,lambda_8);
	e1->SetBinContent(2,8,lambda_8);
	
	e1->SetBinContent(3,1,lambda_3);
	e1->SetBinContent(3,2,lambda_3);
	e1->SetBinContent(3,3,lambda_3);
	e1->SetBinContent(3,4,lambda_4);
	e1->SetBinContent(3,5,lambda_4);
	e1->SetBinContent(3,6,lambda_5);
	e1->SetBinContent(3,7,lambda_6);
	e1->SetBinContent(3,8,lambda_9);
	
	TString name;
	if(isData==false&&year=="2017")name="2017_MC";
	if(isData==true&&year=="2017")name="2017_Data";
	if(isData==false&&year=="2016")name="2016_MC";
	if(isData==true&&year=="2016")name="2016_Data";
	if(isData==false&&year=="2018")name="2018_MC";
	if(isData==true&&year=="2018")name="2018_Data";

	TCanvas *c1=new TCanvas("c1","c1",800,800);//for muon
	TCanvas *c2=new TCanvas("c2","c2",800,800);//for ecalDriven electron
	TCanvas *c4=new TCanvas("c3","c3",800,800);//for tracker driven electron
		
	c1->cd ();
	mu->Draw("TEXT");
	c1->SaveAs(plotpathname+"/"+"muon"+"/"+dataname+"/"+name+"_mu.png");
	c1->SaveAs(plotpathname+"/"+"muon"+"/"+dataname+"/"+name+"_mu.pdf");
	
	c2->cd();
	e1->Draw("TEXT");
	c2->SaveAs(plotpathname+"/"+"ecalDriven"+"/"+dataname+"/"+name+"_e1.png");
	c2->SaveAs(plotpathname+"/"+"ecalDriven"+"/"+dataname+"/"+name+"_e1.pdf");

	c4->cd();
	e3->Draw("TEXT");
	c4->SaveAs(plotpathname+"/"+"trackerDriven"+"/"+dataname+"/"+name+"_e3.png");
	c4->SaveAs(plotpathname+"/"+"trackerDriven"+"/"+dataname+"/"+name+"_e3.pdf");

	TFile *fmu=new TFile(plotpathname+"/"+"muon"+"/"+dataname+"/"+name+"_mu.root","RECREATE");
	fmu->cd();
	mu->Write();
	fmu->Close();
	
	TFile *fe1=new TFile(plotpathname+"/"+"ecalDriven"+"/"+dataname+"/"+name+"_e1.root","RECREATE");
	fe1->cd();
	e1->Write();
	fe1->Close();
	
	TFile *fe3=new TFile(plotpathname+"/"+"trackerDriven"+"/"+dataname+"/"+name+"_e3.root","RECREATE");
	fe3->cd();
	e3->Write();
	fe3->Close();

}
