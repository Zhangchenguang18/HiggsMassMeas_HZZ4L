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
TString year="2018";
TString type="pt_non-closure";
void doplot(TString fs,TString year);
//float BacktoPre(float pt, float eta, float phi, float m, float e_precorr);
//float GetCorr(float pt, float eta, int id);
//TFile* f_mu = new TFile("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/LUT_muon_.root");
//TH2D* LUT_mu = (TH2D*)f_mu->Get("LUT");
//TFile* f_ele = new TFile("/home/chenguan/public_html/Syst_uncertainty/ElectronScaleUncPlanB/"+year+"/ScaleShift/LUT.root");
//TH2D* LUT_ele = (TH2D*)f_ele->Get("LUT");
void ScaleSmearUnc(){
	        
//			TString year="2016";

			doplot("4e", year);//alpha1 alpha2 n1 n2 sigma
//			doplot("4mu", year);
//			doplot("2e2mu", year);
			
			
			// 0.76, 1.41, 7.32, 8.62, 2.0, 0.77, 1.41, 7.25, 9.89, 2.1, 0.76, 1.39, 7.47, 10.15, 2.08 ); 
			
}

void doplot(TString fs, TString year){
		
	Int_t fs1,fs2;
	if(fs=="4e"){fs1=2;fs2=2;}
	if(fs=="4mu"){fs1=1;fs2=1;}
	if(fs=="2e2mu"){fs1=3;fs2=4;}
		
	bool passedFullSelection;
	Int_t finalState, idL1, idL2, idL3, idL4;
	Int_t lep_Hindex[4];
	vector<double>* ele_scaleup=0; vector<double>* ele_scaledn=0;
	vector<double>* ele_smearup=0; vector<double>* ele_smeardn=0;
	vector<double>* ele_smearRhoup=0; vector<double>* ele_smearRhodn=0;
	vector<double>* ele_smearPhiup=0; vector<double>* ele_smearPhidn=0;
	vector<double>* ele_ecaltrkpostcorr=0; vector<double>* ele_ecaltrkprecorr=0;
	vector<double>* ele_scaleupecalonly=0; vector<double>* ele_scalednecalonly=0;
	vector<float>* lep_pt=0; vector<float>* lep_eta=0; vector<float>* lep_phi=0; vector<float>* lep_mass=0; vector<float>* lep_ecalEnergy=0;
	vector<double>* ele_ecalprecorr=0; vector<double>* ele_ecalpostcorr=0;

	RooRealVar *mass4lno=new RooRealVar("mass4lno","mass4lno",105,140); 
	RooDataSet *dataset4lno=new RooDataSet("dataset4lno","dataset4lno",RooArgSet(*mass4lno));
	
	RooRealVar *mass4lup=new RooRealVar("mass4lup","mass4lup",105,140); 
	RooDataSet *dataset4lup=new RooDataSet("dataset4lup","dataset4lup",RooArgSet(*mass4lup));
	
	RooRealVar *mass4ldn=new RooRealVar("mass4ldn","mass4ldn",105,140); 
	RooDataSet *dataset4ldn=new RooDataSet("dataset4ldn","dataset4ldn",RooArgSet(*mass4ldn));
	
	TFile* f=new TFile();
	f = new TFile("/raid/raid9/chenguan/input/EleScaleUnc/"+year+"GGH.root");
	TTree *t=(TTree*)f->Get("Ana/passedEvents");
	int sum=t->GetEntries();
	
	t->SetBranchAddress("finalState",&finalState);
	t->SetBranchAddress("passedFullSelection",&passedFullSelection);
	t->SetBranchAddress("lep_Hindex",&lep_Hindex);
	t->SetBranchAddress("lep_pt",&lep_pt);
	t->SetBranchAddress("lep_eta",&lep_eta);
	t->SetBranchAddress("lep_phi",&lep_phi);
	t->SetBranchAddress("lep_mass",&lep_mass);
	
	t->SetBranchAddress("ele_ecaltrkEnergypostcorr",&ele_ecaltrkpostcorr);
	t->SetBranchAddress("ele_ecaltrkEnergyprecorr",&ele_ecaltrkprecorr);
	if(type=="scale"){	
		t->SetBranchAddress("ele_scaleup",&ele_scaleup);
		t->SetBranchAddress("ele_scaledn",&ele_scaledn);
	}
//	t->SetBranchAddress("lep_ecalEnergy",&ele_ecaltrkpostcorr);
//	t->SetBranchAddress("ele_scaleupecalonly",&ele_scaleup);
//	t->SetBranchAddress("ele_scalednecalonly",&ele_scaledn);
	if(type=="smearing"){
		t->SetBranchAddress("ele_smearup",&ele_scaleup);
		t->SetBranchAddress("ele_smeardn",&ele_scaledn);
	}	
//	t->SetBranchAddress("ele_smearRhoup",&ele_scaleup);
//	t->SetBranchAddress("ele_smearRhodn",&ele_scaledn);
//	t->SetBranchAddress("ele_smearPhiup",&ele_scaleup);
//	t->SetBranchAddress("ele_smearPhidn",&ele_scaledn);
	
	for(Int_t i=0; i<t->GetEntries(); i++){
		t->GetEntry(i);
//		cout<<i<<endl;
		if(passedFullSelection==1&&finalState==2){
//	cout<<"11111"<<endl;		
//	cout<<lep_Hindex[0]<<"  "<<lep_Hindex[1]<<"  "<<lep_Hindex[2]<<"  "<<lep_Hindex[3]<<endl;
			float pt1 = (*lep_pt)[lep_Hindex[0]]; float pt2 = (*lep_pt)[lep_Hindex[1]]; float pt3 = (*lep_pt)[lep_Hindex[2]]; float pt4 = (*lep_pt)[lep_Hindex[3]];
//cout<<"22222"<<endl;
			float eta1 = (*lep_eta)[lep_Hindex[0]]; float eta2 = (*lep_eta)[lep_Hindex[1]]; float eta3 = (*lep_eta)[lep_Hindex[2]]; float eta4 = (*lep_eta)[lep_Hindex[3]];
			float phi1 = (*lep_phi)[lep_Hindex[0]]; float phi2 = (*lep_phi)[lep_Hindex[1]]; float phi3 = (*lep_phi)[lep_Hindex[2]]; float phi4 = (*lep_phi)[lep_Hindex[3]];
			float m1 = (*lep_mass)[lep_Hindex[0]]; float m2 = (*lep_mass)[lep_Hindex[1]]; float m3 = (*lep_mass)[lep_Hindex[2]]; float m4 = (*lep_mass)[lep_Hindex[3]];
			double corr1 = (*ele_ecaltrkpostcorr)[lep_Hindex[0]]/(*ele_ecaltrkprecorr)[lep_Hindex[0]];
			double corr2 = (*ele_ecaltrkpostcorr)[lep_Hindex[1]]/(*ele_ecaltrkprecorr)[lep_Hindex[1]];
			double corr3 = (*ele_ecaltrkpostcorr)[lep_Hindex[2]]/(*ele_ecaltrkprecorr)[lep_Hindex[2]];
			double corr4 = (*ele_ecaltrkpostcorr)[lep_Hindex[3]]/(*ele_ecaltrkprecorr)[lep_Hindex[3]];
			TLorentzVector lep1_no, lep2_no, lep3_no, lep4_no;
			lep1_no.SetPtEtaPhiM(pt1*corr1,eta1,phi1,m1*corr1);
			lep2_no.SetPtEtaPhiM(pt2*corr2,eta2,phi2,m2*corr2);
			lep3_no.SetPtEtaPhiM(pt3*corr3,eta3,phi3,m3*corr3);
			lep4_no.SetPtEtaPhiM(pt4*corr4,eta4,phi4,m4*corr4);
			double mass4l_no = (lep1_no + lep2_no + lep3_no + lep4_no).M();
			//for scale up
			corr1 = (*ele_scaleup)[lep_Hindex[0]]/(*ele_ecaltrkprecorr)[lep_Hindex[0]];
			corr2 = (*ele_scaleup)[lep_Hindex[1]]/(*ele_ecaltrkprecorr)[lep_Hindex[1]];
			corr3 = (*ele_scaleup)[lep_Hindex[2]]/(*ele_ecaltrkprecorr)[lep_Hindex[2]];
			corr4 = (*ele_scaleup)[lep_Hindex[3]]/(*ele_ecaltrkprecorr)[lep_Hindex[3]];
			TLorentzVector lep1_scaleup, lep2_scaleup, lep3_scaleup, lep4_scaleup;
			lep1_scaleup.SetPtEtaPhiM(pt1*corr1,eta1,phi1,m1*corr1);
			lep2_scaleup.SetPtEtaPhiM(pt2*corr2,eta2,phi2,m2*corr2);
			lep3_scaleup.SetPtEtaPhiM(pt3*corr3,eta3,phi3,m3*corr3);
			lep4_scaleup.SetPtEtaPhiM(pt4*corr4,eta4,phi4,m4*corr4);
			double mass4l_scaleup = (lep1_scaleup+lep2_scaleup+lep3_scaleup+lep4_scaleup).M();
			
			//for scale dn
			corr1 = (*ele_scaledn)[lep_Hindex[0]]/(*ele_ecaltrkprecorr)[lep_Hindex[0]];
			corr2 = (*ele_scaledn)[lep_Hindex[1]]/(*ele_ecaltrkprecorr)[lep_Hindex[1]];
			corr3 = (*ele_scaledn)[lep_Hindex[2]]/(*ele_ecaltrkprecorr)[lep_Hindex[2]];
			corr4 = (*ele_scaledn)[lep_Hindex[3]]/(*ele_ecaltrkprecorr)[lep_Hindex[3]];
			TLorentzVector lep1_scaledn, lep2_scaledn, lep3_scaledn, lep4_scaledn;
			lep1_scaledn.SetPtEtaPhiM(pt1*corr1,eta1,phi1,m1*corr1);
			lep2_scaledn.SetPtEtaPhiM(pt2*corr2,eta2,phi2,m2*corr2);
			lep3_scaledn.SetPtEtaPhiM(pt3*corr3,eta3,phi3,m3*corr3);
			lep4_scaledn.SetPtEtaPhiM(pt4*corr4,eta4,phi4,m4*corr4);
			double mass4l_scaledn = (lep1_scaledn+lep2_scaledn+lep3_scaledn+lep4_scaledn).M();
			//
			
			
			
			
			if(mass4l_no<140&&mass4l_no>105){//&&pTL1>7&&pTL2>7&&pTL3>7&&pTL4>7){
				*mass4lno=mass4l_no;
				dataset4lno->add(RooArgSet(*mass4lno));
			}
			if(mass4l_scaleup<140&&mass4l_scaleup>105){//&&pTL1_up<100&&pTL2_up<100&&pTL3_up<100&&pTL4_up<100){//&&pTL1_up>7&&pTL2_up>7&&pTL3_up>7&&pTL4_up>7){
				*mass4lup=mass4l_scaleup;
				dataset4lup->add(RooArgSet(*mass4lup));
			}		
			if(mass4l_scaledn<140&&mass4l_scaledn>105){//&&pTL1_dn<100&&pTL2_dn<100&&pTL3_dn<100&&pTL4_dn<100){//&&pTL1_dn<7&&pTL2_dn<7&&pTL3_dn<7&&pTL4_dn<7){
				*mass4ldn=mass4l_scaledn;
				dataset4ldn->add(RooArgSet(*mass4ldn));
			}
			
		}
		
	}
		//make model
	
	RooRealVar meanDCB1("meanDCB","meanDCB",125,120,130);
	RooRealVar sigmaDCB1("sigmaDCB","sigmaDCB",1,0,10);
	RooRealVar alphaDCB1("alphaDCB","alphaDCB",1,0,10);
	RooRealVar nDCB1("nDCB","nDCB",1,0,10);
        RooRealVar alpha21("alpha2","alpha2",1,0,10);
	RooRealVar n21("n2","n2",1,0,10);
	RooDoubleCB DCBno("DCB","DCB",*mass4lno,meanDCB1,sigmaDCB1,alphaDCB1,nDCB1,alpha21,n21);
				
	DCBno.fitTo(*dataset4lno,PrintLevel(-1),SumW2Error(kTRUE),Timer(kTRUE));
    RooPlot *frame1=mass4lno->frame(Bins(100));
    frame1->SetTitle(fs+"_normal");/////////////////////////////////////////////////////////////////////
    dataset4lno->plotOn(frame1);
    DCBno.plotOn(frame1);
    DCBno.paramOn(frame1,Layout(0.1,0.4,0.9));
	
	Double_t chisquare1=frame1->chiSquare(6);
	TLatex *latex1=new TLatex();
    latex1->SetNDC();
    latex1->SetTextSize(0.05);
    latex1->SetTextFont(42);
    latex1->SetTextAlign(23);
    char chi21[20];
    sprintf(chi21,"%s%1.4f","#chi^{2}/DOF=",chisquare1);
   
	TCanvas* c = new TCanvas("c","c",1400,1000);
	c->cd();
	frame1->Draw();
	latex1->DrawLatex(0.7, 0.8, chi21);
	c->SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+type+"unc_"+fs+"_Normal.png");
	
//	RooRealVar shift_up("shift","shift",0,-1,1);
//	RooFormulaVar meanDCB2("meanDCB2","@1+@0",RooArgList(meanDCB,shift_up);	
	RooRealVar meanDCB2("meanDCB","meanDCB",125,120,130);
//	RooFormulaVar meanDCB_up("meanDCB","@1+@0",RooArgList(meanDCB2,shift_up)); 
	RooRealVar sigmaDCB2("sigmaDCB","sigmaDCB",1,0,10);//sigmaDCB1.getVal());
	RooRealVar alphaDCB2("alphaDCB","alphaDCB",1,0,10);//alphaDCB1.getVal());
	RooRealVar nDCB2("nDCB","nDCB",1,0,10);//nDCB1.getVal());
    RooRealVar alpha22("alpha2","alpha2",1,0,10);//alpha21.getVal());
	RooRealVar n22("n2","n2",1,0,10);//n21.getVal());	
	RooDoubleCB DCBup("DCB","DCB",*mass4lup,meanDCB2,sigmaDCB2,alphaDCB2,nDCB2,alpha22,n22);
	
	DCBup.fitTo(*dataset4lup,PrintLevel(-1),SumW2Error(kTRUE),Timer(kTRUE));
    RooPlot *frame2=mass4lup->frame(Bins(100));
    frame2->SetTitle(fs+"_up");/////////////////////////////////////////////////////////////////////
    dataset4lup->plotOn(frame2);
    DCBup.plotOn(frame2);
    DCBup.paramOn(frame2,Layout(0.1,0.4,0.9));
	
	Double_t chisquare2=frame2->chiSquare(6);
	TLatex *latex2=new TLatex();
    latex2->SetNDC();
    latex2->SetTextSize(0.05);
    latex2->SetTextFont(42);
    latex2->SetTextAlign(23);
    char chi22[20];
    sprintf(chi22,"%s%1.4f","#chi^{2}/DOF=",chisquare2);
	
	c->Clear();
	c->cd();
	frame2->Draw();
	latex2->DrawLatex(0.7,0.8,chi22);
	c->SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+type+"unc_"+fs+"_Up.png");

	
//	RooRealVar shift_dn("shift","shift",0,1,-1);
	RooRealVar meanDCB3("meanDCB","meanDCB",125,120,130);
//	RooFormulaVar meanDCB_dn("meanDCB","@1+@0",RooArgList(meanDCB3,shift_dn));
	RooRealVar sigmaDCB3("sigmaDCB","sigmaDCB",1,0,10);//sigmaDCB1.getVal());
	RooRealVar alphaDCB3("alphaDCB","alphaDCB",1,0,10);//alphaDCB1.getVal());
	RooRealVar nDCB3("nDCB","nDCB",1,0,10);//nDCB1.getVal());
    RooRealVar alpha23("alpha2","alpha2",1,0,10);//alpha21.getVal());
	RooRealVar n23("n2","n2",1,0,10);//n21.getVal());
	RooDoubleCB DCBdn("DCB","DCB",*mass4ldn,meanDCB3,sigmaDCB3,alphaDCB3,nDCB3,alpha23,n23);
	
	DCBdn.fitTo(*dataset4ldn,PrintLevel(-1),SumW2Error(kTRUE),Timer(kTRUE));
    RooPlot *frame3=mass4ldn->frame(Bins(100));
    frame3->SetTitle(fs+"_dn");/////////////////////////////////////////////////////////////////////
    dataset4ldn->plotOn(frame3);
    DCBdn.plotOn(frame3);
    DCBdn.paramOn(frame3,Layout(0.1,0.4,0.9));
	
	Double_t chisquare3=frame3->chiSquare(6);
	TLatex *latex3=new TLatex();
    latex3->SetNDC();
    latex3->SetTextSize(0.05);
    latex3->SetTextFont(42);
    latex3->SetTextAlign(23);
    char chi23[20];
    sprintf(chi23,"%s%1.4f","#chi^{2}/DOF=",chisquare3);
	
	c->Clear();
	c->cd();
	frame3->Draw();
	latex3->DrawLatex(0.7,0.8,chi23);
	c->SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+type+"unc_"+fs+"_Down.png");
	
	
//	TString inforpath="/raid/raid9/chenguan/Mass/SystUceratinty/Uncertainty/"+year+"/"+fs+"/";
	
//	ofstream fout1;
//	fout1.open(inforpath+"normal_noabs.txt");
//	fout1<<meanDCB1.getValV();
//	
//	ofstream fout2;
//	fout2.open(inforpath+"up_noabs.txt");
//	fout2<<meanDCB2.getValV();
//	
//	ofstream fout3;
//	fout3.open(inforpath+"down_noabs.txt");
//	fout3<<meanDCB3.getValV();
	
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!result!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
cout<<"                              "<<"Scale"<<endl;
cout<<"(normal-down)/125 = "<<((meanDCB1.getVal()-meanDCB3.getVal())/125)*100<<" %"<<endl;
cout<<" "<<endl;
cout<<"(up-normal)/125 = "<<((meanDCB2.getVal()-meanDCB1.getVal())/125)*100<<" %"<<endl;
cout<<" "<<endl;
//cout<<"                             "<<"Resolution"<<endl;
//cout<<"normal-down     "<<(sigmaDCB1.getVal()-sigmaDCB2.getVal())/sigmaDCB1.getVal()<<endl;
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;




	delete mass4lno;
	delete dataset4lno;
	delete mass4lup;
	delete dataset4lup;
	delete mass4ldn;
	delete dataset4ldn;
	
	delete t;
	delete frame1;
	delete frame2;
	delete frame3;
	delete latex1;
	delete latex2;
	delete latex3;
	delete c;
}
/*	
float GetCorr(float pt, float eta, int id){
	
	float corr = 0;
//	TFile* f_mu = new TFile("/home/chenguan/public_html/Sys_uncertainty/LepEnergyScale/"+year+"/ScaleShift/LUT_muon_.root");
//	TH2D* LUT_mu = (TH2D*)f_mu->Get("LUT");
//	TFile* f_ele = new TFile("/home/chenguan/public_html/Sys_uncertainty/ElectronScaleUncPlanB/"+year+"/ScaleShift/LUT.root");
//	TH2D* LUT_ele = (TH2D*)f_ele->Get("LUT");
	if(abs(id)==13){
		int xbin = LUT_mu->GetXaxis()->FindBin(pt);
		int ybin = LUT_mu->GetYaxis()->FindBin(abs(eta));
		corr = LUT_mu->GetBinContent(xbin,ybin);
		corr = abs(corr);
	}
	if(abs(id)==11){
		int xbin = LUT_ele->GetXaxis()->FindBin(pt);
		int ybin = LUT_ele->GetYaxis()->FindBin(abs(eta));
		corr = LUT_ele->GetBinContent(xbin,ybin);
		corr = abs(corr);
	}
	return corr;

}
	
float BacktoPre(float pt, float eta, float phi, float m, float e_precorr){
	TLorentzVector lep;
	lep.SetPtEtaPhiM(pt,eta,phi,m);
	float pt_precorr = pt*(e_precorr/lep.Energy());
	return pt_precorr;
}
	
	
	
*/	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
