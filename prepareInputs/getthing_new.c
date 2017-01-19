void getthing_new(){
//	TFile *br_zz = new TFile("../AnalyticMassShape_100_200.root","read");
//	TCanvas *c2=(TCanvas*)br_zz->Get("c1");
//	TGraph *gbr_2e2mu = (TGraph*)c2->GetPrimitive("pdf_Proj[m1,m2]_Norm[m4l]");
//
//	TFile *fnew = new TFile("br_heshy.root","update");
//	fnew->cd();
//	gbr_2e2mu->Write("br_heshy_100_400");
//	fnew->Close();
//	return;
/*
 TChain *tBR_ZZ=new TChain("SelectedTree");
//  //tBR_ZZ->Add("vbf_xsec.root");
  tBR_ZZ->Add("br_bsm.root");
  tBR_ZZ->Draw("xsec:ZZMass","","pl");
  TGraph *lo= (TGraph*) gROOT->FindObject("Graph")->Clone();  
	lo->Sort();
	lo->Draw("al");
	TFile *fxsec = new TFile("width_bsm.root","recreate");
	lo->Write("width_zz");
	fxsec->Close();
	return;
*/	
	TFile *fnew = new TFile("br_heshy.root","read");
	TGraph *gbr_2e2mu = (TGraph*)fnew->Get("br_heshy");
	TGraph *gbr_2e2mu_low= (TGraph*)fnew->Get("br_heshy_100_400");

	TFile *fnew = new TFile("width_bsm.root","read");
	TGraph *gbr_2e2mu_bsm= (TGraph*)fnew->Get("width_2e2mu");

	TFile *BR_ZZ=new TFile("width_ZZ.root","read");
	TH1F *hbr_2e2mu = (TH1F*)BR_ZZ->Get("br_2e2mu");	
	TH1F *hbr_4e= (TH1F*)BR_ZZ->Get("br_4e");	

	for (int i=0;i<gbr_2e2mu->GetN();i++) gbr_2e2mu->GetY()[i] /= gbr_2e2mu->GetX()[i]*3.730471584038695e-09*2.;//*0.9955476009830814;
	hbr_2e2mu->Scale(413.6654661864746/2.);
	hbr_4e->Scale(413.6654661864746/2.);
	cout << gbr_2e2mu->Eval(500)<<endl;
	cout << gbr_2e2mu_bsm->Eval(500)<<endl;
	int bin450=hbr_2e2mu->FindBin(500.);
	cout << hbr_2e2mu->GetBinContent(bin450)<<endl;

	gbr_2e2mu_bsm->Draw("al");
	gbr_2e2mu->Draw("lsame");
	gbr_2e2mu_bsm->SetLineColor(2);
	hbr_2e2mu->Draw("same");
//	return;

	//double yr = hbr_2e2mu->GetBinContent(bin450);
	//hbr_2e2mu->Scale(1./1553.900734131863);
	//hbr_2e2mu->SetLineColor(2);
	//hbr_2e2mu->Draw("same");
	//cout << yr<<endl;
	//return;
//	gbr_2e2mu->Draw("apl");
	cout <<gbr_2e2mu->Eval(400)<<endl;
//	cout<<gbr_2e2mu_low->Eval(400)<<endl;
//	for (int i=0;i<gbr_2e2mu_low->GetN();i++) gbr_2e2mu_low->GetY()[i] /= 3891.5750465418528;


	//gbr_2e2mu_low->SetLineColor(2);
	//gbr_2e2mu_low->Draw("plsame");
//	return;

////	TFile *fnew_200= new TFile("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/Fit/AnalyticMassShape_100_200.root","read");
////	 TCanvas *c2=(TCanvas*)fnew_200->Get("c1");
////	 TGraph *gbr_2e2mu = (TGraph*)c2->GetPrimitive("pdf_Proj[m1,m2]_Norm[m4l]");
////	TFile *fnew = new TFile("br_heshy.root","update");
////	fnew->cd();
////	gbr_2e2mu->Write("br_heshy_low");
////	fnew->Close();



	TFile *flo=new TFile("VBFSM_new.root","read");
  	TGraph *lo= (TGraph*)(flo->Get("lo"));
	lo->Sort();

	TFile *fxsec = new TFile("xsec.root","read");
	fxsec->cd();
  	TGraph *vbf= (TGraph*)(fxsec->Get("vbf"));
  	TGraph *wh= (TGraph*)(fxsec->Get("wh"));
  	TGraph *zh= (TGraph*)(fxsec->Get("zh"));




	TCanvas *c2= new TCanvas("c2","",800,600);

	
//	TGraph *gbr_4e = new TGraph ();
	TH1F *br_4enew= new TH1F("br_4enew","",5800,100.25,3000.25);
	TH1F *br_2e2munew= new TH1F("br_2e2munew","",5800,100.25,3000.25);
	TGraph *gbr_2e2munew = new TGraph(5800); 
	TGraph *gbr_4enew = new TGraph(5800); 
	double xi[5995];
	double yi[5995];
	int j=0;
	for(int i =0; i< 5800;i++){
		double ma = br_4enew->GetBinCenter(i+1); 
//		if(ma>183 && ma<189)
//			continue;
		double wid2e2mu;
		int binn = hbr_2e2mu->FindBin(ma);
		double f2e2mu = hbr_2e2mu->GetBinContent(binn);
		double f4e = hbr_4e->GetBinContent(binn);
		if(ma<180)
			wid2e2mu = f2e2mu; 
		else if(ma<500)
			wid2e2mu = gbr_2e2mu_bsm->Eval(ma);
		else
			wid2e2mu = gbr_2e2mu->Eval(ma);
		if(ma>1000){
//		continue;
			f2e2mu=1.;
			f4e=1./2.;
		}
//		double pro =lo->Eval(ma);
		double pro=vbf->Eval(ma)*0.955;
		if(ma<2000)	
		pro += (wh->Eval(ma)*0.654+zh->Eval(ma)*0.669);
		if(ma==125)
		cout <<"mH125: "<< pro<<endl;

//		double pro =vbf->Eval(ma);

//		if(ma>350&&ma<550)
//		cout << pro<< "\t"<<vbf->Eval(ma)<<endl;
		

//		br_2e2munew->Fill(ma,wid2e2mu);
		gbr_4enew->SetPoint(gbr_4enew->GetN(), ma,wid2e2mu*f4e/f2e2mu*pro);
		gbr_2e2munew->SetPoint(gbr_2e2munew->GetN(), ma,wid2e2mu*pro);
//		gbr_2e2munew->SetPoint(gbr_2e2munew->GetN(), ma,f2e2mu*pro);
//		gbr_4enew->SetPoint(gbr_4enew->GetN(), ma,f4e*pro);
//		gbr_4enew->SetPoint(gbr_4enew->GetN(), ma,wid2e2mu*f4e/(f2e2mu+2*f4e)*pro);
		xi[j]=ma;
		//yi[j]=wid2e2mu*pro;
		//yi[j]=wid2e2mu*pro*f2e2mu/(f2e2mu+2*f4e);
		yi[j]=wid2e2mu*pro;
		j++;
			
		}

		cout<<j<<endl;
		gbr_2e2munew->SetLineColor(8);
		gbr_2e2munew->SetLineWidth(2);
		gbr_2e2munew->Draw("al");
//		gbr_2e2munew->GetXaxis()->SetRangeUser(120,180);
//		gbr_2e2munew->GetYaxis()->SetRangeUser(1.E-10,8.E-8);

//		TSpline3* sp_2e2mu= new TSpline3("sp_2e2mu",xi,yi,j);
//		TSpline3* sp_4e= new TSpline3("sp_4e",gbr_4enew);
//		sp_2e2mu->SetLineColor(1);
//		sp_2e2mu->SetLineWidth(3);
//		sp_2e2mu->Draw("lsame");
////		gbr_2e2mu->Draw("samel");
		gbr_4enew->SetLineColor(2);
		gbr_4enew->Draw("samel");

//		fnew->cd();
//		gbr_2e2munew->Write("width_2e2mu");
//		gbr_4enew->Write("width_4e");
//		fnew->Close();
		TFile *width_new = new TFile("width_new.root","recreate");
		gbr_2e2munew->Write("br_2e2mu");
		gbr_4enew->Write("br_4e");
		//sp_2e2mu->Write("spline_2e2mu");
		return;
/*
	for(int i =0; i< 1500;i++){
		double ma = br_4enew->GetBinCenter(i+1); 
		int binn = hbr_2e2mu->FindBin(ma);
		double f2e2mu = hbr_2e2mu->GetBinContent(binn);
		double f4e = hbr_4e->GetBinContent(binn);
		double wid4e = f4e/f2e2mu*wid2e2mu;
		br_4enew->Fill(ma,wid4e);
		//cout << ma << "\t"<<wid2e2mu <<"\t"<< f2e2mu<<endl;
		//br_4enew->Fill(ma,f2e2mu/wid2e2mu/1500.);
	}
	//gbr_2e2mu->Draw("al");
	
	br_4enew->SetLineColor(2);
	br_4enew->Draw("l");
	return;

  TGraph *gr0 = new TGraph(br_4enew); 
	gr0->Draw("Al");
	gr0->GetXaxis()->SetRangeUser(100,200);
	gPad->Update();
	gbr_2e2mu->Draw("lsame");
	return;
//	double scal = 0.000705877;
//	hbr_2e2mu->Scale(scal);
//	hbr_2e2mu->Draw("same");
//	cout << gbr_2e2mu->Eval(450)/yr<<endl;
*/
}
