void CheckEBE(){
	
	TFile* f_old = new TFile("/raid/raid9/chenguan/Mass/doMeasurement/CreateDataCards_Run2_testIssue_20200622/templates2D_refit/Dm_signal_4mu.root");
	TFile* f_new = new TFile("/raid/raid9/chenguan/Mass/doMeasurement/CreateDataCards_Run2_testIssue_20200622/templates2016REFIT/Dm_signal_4mu.root");

	TH1F* sig_old = (TH1F*)f_old->Get("h_Dm");
	sig_old->SetLineColor(kBlue);
//	TH1F* qqzz_reco = (TH1F*)f_reco->Get("");
//	TH1F* ggzz_reco = (TH1F*)f_reco->Get("");

	TH1F* sig_new = (TH1F*)f_new->Get("h_Dm");
	sig_new->SetLineColor(kRed);
//	TH1F* qqzz_refit = (TH1F*)f_refit->Get("");
//	TH1F* ggzz_refit = (TH1F*)f_refit->Get("");

	TCanvas c("c","c",1400,1000);
	c.cd();
	sig_old->Draw("hist");
	sig_new->Draw("histsame");


	TLatex *latex=new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.05);
        latex->SetTextFont(42);
        latex->SetTextAlign(13);
        latex->SetTextColor(kRed);
        latex->DrawLatex(0.7,0.7,"new");
        latex->SetTextColor(kBlue);
        latex->DrawLatex(0.7,0.6,"old");
	
	c.SaveAs("/home/chenguan/public_html/TESTPLOTS/sig_refit_new_old.png");

/*	c.Clear();
	sig_refit->Draw("hist");
	c.SaveAs("/home/chenguan/public_html/TESTPLOTS/sig_refit.png");
	c.Clear();
	sig_reco->Add(sig_refit,-1);
	sig_reco->Draw("hist");
	c.SaveAs("/home/chenguan/public_html/TESTPLOTS/sig_reco-refit.png");
	c.Clear();
	c.cd();
*/	

}
