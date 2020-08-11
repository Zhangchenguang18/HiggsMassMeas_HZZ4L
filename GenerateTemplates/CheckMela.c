void CheckMela(){
	TString fs = "2e2mu";
	TString Process = "Dbackground_ZX_"+fs;
	TFile* f_reco = new TFile("/raid/raid8/mhl/HZZ4L_Run2_post2016ICHEP/HiggsMass_HZZ4L/packages/doMeasurement/CreateDatacards_Moriond17_unblindv5_updateZX_sys_20170307/templates2D/"+Process+".root");
	TFile* f_refit = new TFile("/raid/raid8/mhl/HZZ4L_Run2_post2016ICHEP/HiggsMass_HZZ4L/packages/doMeasurement/CreateDatacards_Moriond17_unblindv5_updateZX_sys_20170307/templates2D_refit/"+Process+".root");
	
	TH2F* h1_reco = (TH2F*)f_reco->Get("h_mzzD");
	TH2F* h2_reco = (TH2F*)f_reco->Get("h_mzzD_up");
	TH2F* h3_reco = (TH2F*)f_reco->Get("h_mzzD_dn");

	TH2F* h1_refit = (TH2F*)f_refit->Get("h_mzzD");
	TH2F* h2_refit = (TH2F*)f_refit->Get("h_mzzD_up");
	TH2F* h3_refit = (TH2F*)f_refit->Get("h_mzzD_dn");
	
	double nor = h1_reco->Integral();
	cout<<"!!!!!!!!!!!!!"<<nor<<endl;	
	TCanvas c("c","c",1400,1000);
/*	
	c.cd();
	if(h1_reco){
		h1_reco->Draw("colz");
		c.SaveAs("/home/chenguan/public_html/TESTPLOTS/"+Process+"_reco_no.png");
		c.Clear();
	}
	
	c.cd();
	if(h2_reco){
		h2_reco->Draw("colz");
		c.SaveAs("/home/chenguan/public_html/TESTPLOTS/"+Process+"_reco_up.png");
		c.Clear();
	}
	
	c.cd();
	if(h3_reco){
		h3_reco->Draw("colz");
		c.SaveAs("/home/chenguan/public_html/TESTPLOTS/"+Process+"_reco_dn.png");
		c.Clear();
	}

	c.cd();
	if(h1_refit){
		h1_refit->Draw("colz");
		c.SaveAs("/home/chenguan/public_html/TESTPLOTS/"+Process+"_refit_no.png");
		c.Clear();
	}

	c.cd();
	if(h2_refit){
		h2_refit->Draw("colz");
		c.SaveAs("/home/chenguan/public_html/TESTPLOTS/"+Process+"_refit_up.png");
		c.Clear();
	}

	c.cd();
	if(h3_refit){
		h3_refit->Draw("colz");
		c.SaveAs("/home/chenguan/public_html/TESTPLOTS/"+Process+"_refit_dn.png");
		c.Clear();
	}
	
	if(h1_reco&&h2_reco){
		h1_reco->Add(h2_reco,-1);
		c.cd();
		h1_reco->Draw("colz");
		c.SaveAs("/home/chenguan/public_html/TESTPLOTS/"+Process+"_reco_no-up.png");
		c.Clear();
	}
	if(h2_reco&&h3_reco){
		h2_reco->Add(h3_reco,-1);
		c.cd();
		h2_reco->Draw("colz");
		c.SaveAs("/home/chenguan/public_html/TESTPLOTS/"+Process+"_reco_up-dn.png");
		c.Clear();

		}
		*/
//	if(h1_reco&&h3_reco){
///		h3_reco->Add(h3_reco,-1);
//		c.cd();
//		h3_reco->Draw("colz");
//		c.SaveAs("/home/chenguan/public_html/TESTPLOTS/"+fs+"/"+Process+"_reco_dn.png");
//		c.Clear();
//	}
/*	h1_reco->Add(h2_reco,-1);
	c.cd();
	h1_reco->Draw("colz");
	c.SaveAs("/home/chenguan/public_html/TESTPLOTS/"+fs+"/"+Process+"_reco_no-up.png");
	c.Clear();

	h1_refit->Add(h2_refit,-1);
	c.cd();
	h1_refit->Draw("colz");
	c.SaveAs("/home/chenguan/public_html/TESTPLOTS/"+fs+"/"+Process+"_refit_no-up.png");
	c.Clear();
*/

}
