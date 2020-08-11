#ifndef PileupWeight_h
#define PileupWeight_h

TFile *f_pilup  = new TFile("include/pileup_MC_80x_2016B_Run_271036-274443_mb_69735.root","READ");
TH1D *puweight_dtmc = (TH1D*) f_pilup->Get("puweight_dtmc");

double puweight(int nPV){

 return puweight_dtmc->GetBinContent(nPV);

}

#endif
