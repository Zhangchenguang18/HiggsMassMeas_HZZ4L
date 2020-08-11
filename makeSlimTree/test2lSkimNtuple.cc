void test2lSkimNtuple(){
	std::vector<int> *lep_id;
	std::vector<float> *lep_pt;
	Double_t mass4l;
	TFile* f = TFile::Open("root://cmsio5.rc.ufl.edu//store/user/t2/users/ferrico/Full_RunII/madgraph/2016.root");
	TTree* t = (TTree*)f->Get("Ana/pessedEvents");
	t->SetBranchAddress("lep_id", &lep_id);
	t->SetBranchAddress("lep_pt",&lep_pt);
	t->SetBranchAddress("mass4l",&mass4l);
	Int_t lepptsize=0;
	Int_t lepidsize=0;

	for(Int_t i=0; i<t->GetEntries(); i++){
		if((*lep_pt).size()!=4)continue;
		cout<<mass4l<<endl;
	}
	
	
}
