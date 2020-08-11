void PrintElectronLUT(){
	
	TFile* f = new TFile("/home/chenguan/public_html/getLambda/ecalDriven/Data/2018/2018_Data_e1.root");
	TH2D* h = (TH2D*)f->Get("e1");
//	for(Int_t i=1; i<4; i++){
		for(Int_t j=8; j>0; j--){
			cout<<h->GetBinContent(1,j)<<"  "<<h->GetBinContent(2,j)<<"  "<<h->GetBinContent(3,j)<<endl;
		}
//	}
	h->SetStats(0);
//	h->GetYaxis()->SetLimits(0,0.1);
	TCanvas c("c","c",1000,1000);
	c.cd();
	c.SetLogy();
	h->Draw("text");
//	c.SaveAs("/home/chenguan/public_html/TESTPLOTS/llllllll.png");

}
