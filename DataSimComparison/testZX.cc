void testZX(){
	TFile* f = new TFile("/raid/raid9/chenguan/input/ZX.root");
	TTree* t = (TTree*)f->Get("passedEvents");
	Float_t mass4l;
	bool passedFullSelection;
	t->SetBranchAddress("mass4l",&mass4l);
	t->SetBranchAddress("passedFullSelection",&passedFullSelection);
	TH1F* h = new TH1F("h","h",100,70,170);
	for(Int_t i=0; i<t->GetEntries(); i++){
		t->GetEntry(i);
		if(passedFullSelection){
			h->Fill(mass4l);
		}
	}
	TCanvas c("c","c",1400,1000);
	c.cd();
	h->Draw();
	c.SaveAs("/home/chenguan/public_html/TESTPLOTS/ZX.png");
}

