Double_t* CheckEventYield_(TString process);
TString year="2016";

void CheckEventYield(){
	
	TString process="bbH";
//	Double_t* bbH = CheckEventYield_(process);
	Double_t* bbH = new Double_t[4];bbH[0]=0;bbH[1]=0;bbH[2]=0;bbH[3]=0;	

	process = "ttH";
	Double_t* ttH = CheckEventYield_(process);

	process = "WminusH";
	Double_t* WminusH = CheckEventYield_(process);
	
	process = "WplusH";
	Double_t* WplusH = CheckEventYield_(process);

	process = "ZH";
	Double_t* ZH = CheckEventYield_(process);

	process = "VBF";
	Double_t* VBF = CheckEventYield_(process);

	process = "ggH";
	Double_t* ggH = CheckEventYield_(process);

	Double_t fourE = bbH[0] + ttH[0] + ZH[0] + WplusH[0] + WminusH[0] + VBF[0] + ggH[0];
	Double_t fourMu = bbH[1] + ttH[1] + ZH[1] + WplusH[1] + WminusH[1] + VBF[1] + ggH[1];
	Double_t twoEtwoMu = bbH[2] + ttH[2] + ZH[2] + WplusH[2] + WminusH[2] + VBF[2] + ggH[2];
	Double_t sum = bbH[3] + ttH[3] + ZH[3] + WplusH[3] + WminusH[3] + VBF[3] + ggH[3];
	
	cout<<"*******************************************************************"<<endl;
	cout<<"Sum: "<<sum<<endl;
	cout<<"4e: "<<fourE<<endl;
	cout<<"4mu: "<<fourMu<<endl;
	cout<<"2e2mu: "<<twoEtwoMu<<endl;

	process = "qqZZ";
	Double_t* qqzz = CheckEventYield_(process);

	process = "ggZZ_4e";
	Double_t* ggZZ_4e = CheckEventYield_(process);

	process = "ggZZ_4mu";
	Double_t* ggZZ_4mu = CheckEventYield_(process);

	process = "ggZZ_2e2mu";
	Double_t* ggZZ_2e2mu = CheckEventYield_(process);

	process = "ggZZ_4tau";
	Double_t* ggZZ_4tau = CheckEventYield_(process);

	process = "ggZZ_2e2tau";
	Double_t* ggZZ_2e2tau = CheckEventYield_(process);

	process = "ggZZ_2mu2tau";
	Double_t* ggZZ_2mu2tau = CheckEventYield_(process);

	
	cout<<"ggZZ_Sum: "<<ggZZ_4e[3]+ggZZ_4mu[3]+ggZZ_4tau[3]+ggZZ_2e2mu[3]+ggZZ_2e2tau[3]+ggZZ_2mu2tau[3]<<endl;
	cout<<"ggZZ_4e: "<<ggZZ_4e[0]+ggZZ_4mu[0]+ggZZ_4tau[0]+ggZZ_2e2mu[0]+ggZZ_2e2tau[0]+ggZZ_2mu2tau[0]<<endl;
	cout<<"ggZZ_4mu: "<<ggZZ_4e[1]+ggZZ_4mu[1]+ggZZ_4tau[1]+ggZZ_2e2mu[1]+ggZZ_2e2tau[1]+ggZZ_2mu2tau[1]<<endl;
	cout<<"ggZZ_2e2mu: "<<ggZZ_4e[2]+ggZZ_4mu[2]+ggZZ_4tau[2]+ggZZ_2e2mu[2]+ggZZ_2e2tau[2]+ggZZ_2mu2tau[2]<<endl;

}

Double_t* CheckEventYield_(TString process){
	
	Float_t mass4l, eventWeight, genWeight, dataMCWeight, pileupWeight,crossSection,k_ggZZ,k_qqZZ_ewk,k_qqZZ_qcd_M,k;
	bool passedFullSelection;
	Int_t finalState;
	Double_t Lumi;

	if(year=="2016")Lumi=35.92;
        if(year=="2017")Lumi=41.53;
        if(year=="2018")Lumi=59.68;

	Double_t xs;
	
	TFile* f = new TFile();
	if(process=="ggH"){f = new TFile("/raid/raid9/chenguan/input/AfterSync_200424/"+year+"GGH_125.root");xs=0.013335/0.01212;}
	if(process=="VBF"){f = new TFile("/raid/raid9/chenguan/input/AfterSync_200424/"+year+"VBF_125.root");xs=0.001038/0.0010339;}
	if(process=="WplusH"){f = new TFile("/raid/raid9/chenguan/input/AfterSync_200424/"+year+"WplusH_125.root");xs=0.000231/0.0002339;}
	if(process=="WminusH"){f = new TFile("/raid/raid9/chenguan/input/AfterSync_200424/"+year+"WminusH_125.root");xs=0.000146/0.0001471;}
	if(process=="ZH"){f = new TFile("/raid/raid9/chenguan/input/AfterSync_200424/"+year+"ZH_125.root");xs=0.000662/0.0006568;}
	if(process=="ttH"){f = new TFile("/raid/raid9/chenguan/input/AfterSync_200424/"+year+"ttH_125.root");xs=0.00039/0.0003899;}
//	if(process=="bbH"){f = new TFile("");}
	
	if(process=="qqZZ"){f = new TFile("/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/slimmed_"+year+"qqZZ.root");xs=1;}
	if(process=="ggZZ_4e"){f = new TFile("/raid/raid9/chenguan/input/AfterSync_200424/"+year+"ggZZ_4e.root");xs=1;}
	if(process=="ggZZ_4mu"){f = new TFile("/raid/raid9/chenguan/input/AfterSync_200424/"+year+"ggZZ_4mu.root");xs=1;}
	if(process=="ggZZ_2e2mu"){f = new TFile("/raid/raid9/chenguan/input/AfterSync_200424/"+year+"ggZZ_2e2mu.root");xs=1;}	
	if(process=="ggZZ_4tau"){f = new TFile("/raid/raid9/chenguan/input/AfterSync_200424/"+year+"ggZZ_4tau.root");xs=1;}
	if(process=="ggZZ_2mu2tau"){f = new TFile("/raid/raid9/chenguan/input/AfterSync_200424/"+year+"ggZZ_2mu2tau.root");xs=1;}
	if(process=="ggZZ_2e2tau"){f = new TFile("/raid/raid9/chenguan/input/AfterSync_200424/"+year+"ggZZ_2e2tau.root");xs=1;}

	TTree* t = (TTree*)f->Get("Ana/passedEvents");
	if(!t)t = (TTree*)f->Get("passedEvents");
	
	t->SetBranchStatus("*",0);
//	t->SetBranchStatus("crossSection",1);
	t->SetBranchStatus("mass4l",1);
	t->SetBranchStatus("passedFullSelection",1);
	t->SetBranchStatus("finalState",1);
	t->SetBranchStatus("eventWeight",1);
	t->SetBranchStatus("k_ggZZ",1);
	t->SetBranchStatus("k_qqZZ_qcd_M",1);
	t->SetBranchStatus("k_qqZZ_ewk",1);

//	t->SetBranchAddress("crossSection",&crossSection);
	t->SetBranchAddress("mass4l",&mass4l);
	t->SetBranchAddress("passedFullSelection",&passedFullSelection);
	t->SetBranchAddress("finalState",&finalState);
	t->SetBranchAddress("eventWeight",&eventWeight);
	t->SetBranchAddress("k_ggZZ",&k_ggZZ);
	t->SetBranchAddress("k_qqZZ_qcd_M",&k_qqZZ_qcd_M);
	t->SetBranchAddress("k_qqZZ_ewk",&k_qqZZ_ewk);

	Double_t e = 0;
	Double_t mu = 0;
	Double_t emu = 0;
	Double_t sum = 0;
	Double_t result1 = 0;
	Double_t result2 = 0;
	Double_t result3 = 0;
	unsigned long Nloop;
	unsigned long Ntot;
	Nloop = t->GetEntries();
	Ntot = t->GetEntries();
	
	if(process=="qqZZ"&&year=="2016")Ntot=77998181+18470796;
	if(process=="qqZZ"&&year=="2017")Ntot=50756359+36007127;
	if(process=="qqZZ"&&year=="2018")Ntot=19089600;
	for(Int_t i=0; i<Nloop; i++){
		t->GetEntry(i);
		if(process=="ggH"||process=="VBF"||"ZH"||"WplusH"||"WminusH"||"ttH")k=1;
		if(process=="ggZZ_4e"||process=="ggZZ_4mu"||process=="ggZZ_4tau"||process=="ggZZ_2e2mu"||process=="ggZZ_2e2tau"||process=="ggZZ_2mu2tau")k=k_ggZZ;
		if(process=="qqZZ")k=k_qqZZ_ewk*k_qqZZ_qcd_M;		

		if(mass4l>80&&mass4l<100&&passedFullSelection==1){
			if(finalState==1){mu = mu + eventWeight*k;}
			if(finalState==2){e = e + eventWeight*k;}
			if(finalState==3||finalState==4){emu = emu + eventWeight*k;}
		}
	}
	delete t;
	result1 = (e*xs*Lumi*1000)/Ntot;
	result2 = (mu*xs*Lumi*1000)/Ntot;
	result3 = (emu*xs*Lumi*1000)/Ntot;
	sum = result1 + result2 + result3;	
	cout<<"*******************************************************************************************************"<<endl;
	cout<<process<<" "<<"4e: "<<result1<<" 4mu: "<<result2<<" 2e2mu: "<<result3<<" tot: "<<sum<<endl;
	Double_t* results = new Double_t[4];
	results[0] = result1;
	results[1] = result2;
	results[2] = result3;
	results[3] = sum;
	return results;

}


