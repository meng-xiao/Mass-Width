/** Things to remember later
 * mzz range defined by ggZZ workspace, create 2 versions for onshell and offshell
 * give different ggzz RooKeysPdf name for 2e2mu and 4e
 *
 ****/
using namespace std;
using namespace RooFit;

void dosomething(TString chan ="2e2mu",bool cate_vbf =false,bool onshell=false){

//	double lumi = 9.235;
	double lumi = 12.9;
	double ggzz_xsec = 2.827;
	double x_xsec = 1.1584*1.19*1.13; 

	if(onshell){
		if(chan=="2e2mu"){
		ggzz_xsec = 0.2184;
		x_xsec = 0.99794*1.19*1.13; 
		}
		else{
			ggzz_xsec = 0.1092;
			x_xsec = 0.540*1.19*1.13;
		}
	}
		else{
			if(chan!="2e2mu"){
			ggzz_xsec=1.414;
			x_xsec = 0.5989*1.19*1.13; 
			}
		}

	gStyle->SetOptStat(0);
	double dcbPara_2e2mu_2nd[6][11]={
		1000,2000,0.741753, 0.00110575, -3.38327e-07, 1.15192, 0.000285416, 7.18395e-08, -1.39387, 0.00283121, -5.64608e-07, 
		1000,2000,1.75226, 0.000169345, -1.25396e-07, 1.92894, -0.000184015, 5.1284e-08, 1.28478, 0.000460145, -1.09756e-07, 
		1000,2000,0.032571, -0.00340172, 9.33694e-07, -1.82225, 0.000307922, -9.21127e-07, 8.73742, -0.0102517, 1.71879e-06, 
		1000,2000,3.38794, -0.00216507, 7.89938e-07, 2.41936, -0.00022791, -1.78642e-07, 6.53698, -0.00434553, 8.50763e-07, 
		1000,2000,2.862, 0.00144161, 5.2692e-07, 0.880902, 0.00540381, -1.45418e-06, 2.86996, 0.00341475, -9.56913e-07, 
		1000,2000,0.564531, 0.00639651, 6.34132e-06, -1.98043, 0.0114864, 3.79636e-06, -22.3317, 0.0318377, -1.29146e-06
	};
	double dcbPara_4mu_2nd[6][11]=
	{
		1000,2000,1.2122, 0.000173393, 5.12448e-07, -0.27018, 0.00313815, -9.69932e-07, 2.15928, 0.000708693, -3.62567e-07, 
		1000,2000,1.84735, 5.4914e-05, -3.83834e-08, 1.33726, 0.00107509, -5.48473e-07, 3.46976, -0.00105741, -1.53484e-08, 
		1000,2000,0.0923286, -0.00227067, -3.25584e-07, 5.44305, -0.0129721, 5.02514e-06, -58.1031, 0.050574, -1.08614e-05, 
		1000,2000,1.85848, 0.00123506, -1.40928e-06, 3.96505, -0.00297808, 6.97288e-07, 5.83482, -0.00484785, 1.16473e-06, 
		1000,2000,2.78568, 0.00166347, -1.48862e-07, 7.89859, -0.00856235, 4.96405e-06, -30.2704, 0.0296066, -4.5782e-06, 
		1000,2000,-0.0175753, 0.00735663, 1.1293e-05, -9.43701, 0.0261955, 1.87352e-06, -35.9749, 0.0527334, -4.76095e-06 
	};
	double dcbPara_4e_2nd[6][11]=
	{
		1000,2000,0.917418, -0.000280427, 2.75712e-07, 0.492322, 0.000569765, -1.49384e-07, 1.9491, -0.000887013, 2.14811e-07, 
		1000,2000,1.88487, -0.000477813, 3.69221e-07, 1.47359, 0.000344747, -4.2059e-08, 2.37719, -0.000558853, 1.83841e-07, 
		1000,2000,-0.622834, -0.00150395, 7.52586e-07, -3.02999, 0.00331036, -1.65457e-06, 6.31176, -0.00603139, 6.80867e-07, 
		1000,2000,5.8324, -0.00801147, 4.30593e-06, 1.04451, 0.00156431, -4.81962e-07, 2.41795, 0.00019087, -1.38602e-07, 
		1000,2000,2.82558, 0.00162296, -7.92314e-07, 2.88448, 0.00150516, -7.33414e-07, 5.11166, -0.00072202, -1.76619e-07, 
		1000,2000,1.55048, 0.00397066, -3.49938e-07, 3.16794, 0.00073574, 1.26752e-06, -1.96982, 0.0058735, -1.69175e-08 
	};

	double eff[3][11]={
		//10% before
		//    4e
//		-4.413834E+00,4.616797E+00,-4.159802E+02,3.640200E+02, 1.981811E+00, 1.056963E-03,-3.426846E-07,1.598360E-01, 1.600000E+02,5.811111E+01,2.579534E-11,
//		//    4mu
//		-4.434596E+00, 4.608567E+00, -1.103128E+02, 1.314750E+02, 2.738879E+00, 2.937410E-03, -1.450741E-06, 1.554714E-02,1.600001E+02, 2.815378E+01,2.282859E-10,
//		//    2e2mu
//		-4.430343E+00,4.609577E+00, -7.804443E+01, 1.183081E+02, 2.305252E+00, 2.203730E-03, -9.924391E-07, 3.939525E-11, 1.776605E+02, 5.012696E+01,1.400622E-10

//		4e
//		 -4.43896,4.60067,-88.9451,124.26,2.00739,0.00166935,-7.03222e-07,6.37113e-13,208.193,12.1879,8.99317e-11,
//		 //4mu
//		 -4.44961,4.5883,-493.54,364.615,3.30893,0.00282812,-1.34105e-06,0.133779,160.001,55.7017,2.01351e-10,
//		 //2e2mu
//		 -4.43738,4.60043,-380.597,314.768,2.48908,0.00177297,-7.67946e-07,0.127317,160,49.2503,1.02808e-10

////		ggH->0p_4e_param_0:
//		-4.4663,4.5754,-124.303,136.744,2.82937,0.0025389,-1.07711e-06,0.0157579,172.842,24.8505,1.38559e-10,
////		ggH->0p_4mu_param_0:
//		-4.46301,4.57444,-479.667,341.556,4.33312,0.00347696,-1.5298e-06,0.107783,160,49.898,2.11056e-10,
////		ggH->0p_2e2mu_param_0:
//		-4.45,4.58669,-398.384,316.466,2.99933,0.00206663,-8.44269e-07,0.124652,160,48.1899,1.05772e-10

//		ggH->0p_4e_param_0:
		-4.4781,4.56478,-410.689,302.175,3.54126,0.00226419,-9.08486e-07,0.0861973,160,45.6037,1.07622e-10,
	//	ggH->0p_4mu_param_0:
		-4.44313,4.59816,-2313.39,1129.31,3.13553,0.00269624,-1.20507e-06,-9.71235,-105.225,79.0966,1.66983e-10,
	//	ggH->0p_2e2mu_param_0:
		-4.45465,4.5858,-313.199,295.457,3.09284,0.00192647,-7.78952e-07,0.36715,105.723,81.4614,9.69396e-11
	};

	double dcbPara_2nd[6][11];
	double effsig[11];

	if (chan=="4e")   {
		for (int i=0;i<6;i++){for(int j=0;j<11;j++){dcbPara_2nd[i][j]= dcbPara_4e_2nd[i][j];}}
		for (int i=0;i<11;i++){effsig[i]=eff[0][i];}

	}
	if (chan=="4mu")  {
		for (int i=0;i<6;i++){for(int j=0;j<11;j++){dcbPara_2nd[i][j]= dcbPara_4mu_2nd[i][j];}}
		for (int i=0;i<11;i++){effsig[i]=eff[1][i];}
	}
	if (chan=="2e2mu") {
		for (int i=0;i<6;i++){for(int j=0;j<11;j++){dcbPara_2nd[i][j]= dcbPara_2e2mu_2nd[i][j];}}
		for (int i=0;i<11;i++){effsig[i]=eff[2][i];}
	}


	double lowarr[2]={100.5,100.5};
	double higharr[2]={1500.5,150.5};
	const int nbinsarr[2]={1500,500};

	double recolowarr[2]={104,105};
	double recohigharr[2]={1604.,140.};
	const int reconbinsarr[2]={750,100};

	const double low= lowarr[onshell];
	const double high=higharr[onshell];
	const int nbins= nbinsarr[onshell]; 

	const double low_reco=recolowarr[onshell];
	const double high_reco=recohigharr[onshell];
	const int nbins_reco=reconbinsarr[onshell];

	cout << low<<"\t"<<high<<endl;
	cout << low_reco<<"\t"<<high_reco<<endl;

 TString ap = "";
   if(onshell)
		     ap="_onshell";
	TFile *fpdfbkg = new TFile("pdfs"+ap+".root");
	RooWorkspace *wbkg =( RooWorkspace*)fpdfbkg->Get ("w"); 

	//RooRealVar* mzz = new RooRealVar("ZZMass","M_{ZZ} (GeV)",125,low,high);
	RooRealVar* mzz = wbkg->var("ZZMass"); 
	RooRealVar* mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",125,low_reco,high_reco);
	RooRealVar* mdiff= new RooRealVar("mdiff","M_{ZZ} (GeV)",125,low_reco,high_reco);

	RooRealVar *r= new RooRealVar("r","signal strength",1.,0.0001,1000);

	//	mreco->setBins(nbins_reco);
	RooRealVar* mean = new RooRealVar("mean_pole","mean_pole",125,100,180);
	RooRealVar* sigma= new RooRealVar("sigma_pole","sigma_pole",0.00418,0.,10.);

	RooConstVar* mean_125= new RooConstVar("mean_125","mean_125",125);
	RooConstVar* sigma_125= new RooConstVar("sigma_125","sigma_125",0.00407);

	RooPlot* frame= mreco->frame(low_reco,high_reco) ;
	//	RooPlot* frame= mreco->frame(150,250) ;
	RooPlot* frame_mzz= mzz->frame(Title("Z mass")) ;
	RooPlot* frame_width= sigma->frame(Title("width")) ;
	RooPlot* frame_mean= mean->frame(Title("mean")) ;

	TFile *flo=new TFile("ggh_input_spline.root","read");
	//TFile *flo=new TFile("xsec_ggzz4l_13TeV_4e.root","read");
//	TFile *flo=new TFile("width_new.root","read");
	TString chn = "";
	if(chan!="2e2mu")
		chn="_4e";
//	TGraph *lo=(TGraph*) flo->Get("gr_"+chn);
	//TGraph *lo=(TGraph*) flo->Get("br_"+chn);
	TSpline3 *lo=(TSpline3*) flo->Get("sp_xsec_statnom"+chn);

	SplinePdf par2_int ("par2_int"+chan+Form("%d",cate_vbf),"",*mzz,*mean,*sigma,*lo);
	RooRealVar m_gauss ("m_gauss","",125);
	RooRealVar w_gauss ("w_gauss","",0.004);
	RooGaussian gauss("gauss","",*mzz,m_gauss,w_gauss);

	TString pdfn = "2e2mu";
	if(chan!="2e2mu")
		pdfn = "4e";
	RooKeysPdf *pdfbkg = wbkg->pdf("pdfbkg_"+pdfn);

	//	RooKeysPdf *pdfbkg = (RooKeysPdf*)fpdfbkg->Get("pdfbkg");

	//	TChain *ggzz = new TChain ("SelectedTree");
	//	ggzz->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/rootfiles/ggZZ_Bkg_xcheck.root");
	//	ggzz->Draw("ZZMass>>hggzz","ZZMass>100");
	//	TTree *cuttree = ggzz->CloneTree(0);
	//
	////	    for(int i =0;i<100000;i++){
	//	for(int i =0;i<2000;i++){
	//		ggzz->GetEntry(i);
	//		cuttree->Fill();
	//	}
	//	RooDataSet bkgdata ("bkgdata","",cuttree,*mzz);
	//	RooKeysPdf *pdfbkg=new RooKeysPdf("pdfbkg","",*mzz,bkgdata,RooKeysPdf::MirrorLeft);


	RooConstVar *ggzznorm= new RooConstVar("ggzznorm"+chan+Form("%d",cate_vbf),"",lumi*ggzz_xsec);
	RooExtendPdf pdf_ggzz("pdf_ggzz"+chan+Form("%d",cate_vbf),"pdf_ggzz"+chan+Form("%d",cate_vbf),*pdfbkg,*ggzznorm);

	RooConstVar *xnorm= new RooConstVar("xnorm"+chan+Form("%d",cate_vbf),"",lumi*x_xsec);
	RooExtendPdf pdf_x("pdf_x"+chan+Form("%d",cate_vbf),"pdf_x"+chan+Form("%d",cate_vbf),par2_int,*xnorm);


	TFile *fphase_noweight=new TFile("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/Fit/whatthefuck/fphase_ggH.root");
	TGraph *phase_sin = fphase_noweight->Get("sinspline");
	TGraph *phase_cos = fphase_noweight->Get("cosspline");

//	TFile *fkfactor = new TFile("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/Fit/whatthefuck/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
//		//TSpline3* ggZZ_kf [9];//= (TSpline3*)fkfactor->Get("sp_kfactor_Nominal");
//		TSpline3* ggZZ_kf = (TSpline3*)fkfactor->Get("sp_kfactor_Nominal");
 TString strSystTitle[5] ={
  "Nominal",
  "qcd_dn",
  "qcd_up",
  "pdf_dn",
  "pdf_up"
};
	TFile *fkfactor = new TFile("kfactor.root");

/**** calculate the xsec after kfactor ****/
/*
    double sum=0;
    double sum_ori=0;
    
    mean->setVal(125);
    sigma->setVal(0.004);
    
    for (int i =0;i<100;i++){
      double ma = 125-10*0.004+0.2*0.004*i; 
      mzz->setVal(ma);
      double  kfa = ggZZ_kf->Eval(ma);
      double  pdfval = pdf_ggzz.getVal(*mzz);
      double rev = pdfval*kfa*0.2*0.004;
      sum+=rev;
      sum_ori+=rev/kfa;
    }
    cout << sum<< "\t"<< sum_ori<<endl;
    for(int i =1;i<20;i++){
      double ma = 125-10*0.004-i; 
      mzz->setVal(ma);
      double  kfa = ggZZ_kf->Eval(ma);
      double  pdfval = pdf_ggzz.getVal(*mzz);
      double rev = pdfval*kfa;
      sum+=rev;
      sum_ori+=rev/kfa;
    }
    cout << sum<< "\t"<< sum_ori<<endl;
    for(int i =1;i<1000;i++){
//    for(int i =1;i<16;i++){
      double ma = 125+10*0.004+i; 
      mzz->setVal(ma);
      double  kfa = ggZZ_kf->Eval(ma);
      double  pdfval = pdf_ggzz.getVal(*mzz);
      double rev = pdfval*kfa;
      sum+=rev;
      sum_ori+=rev/kfa;
			if(ma>140&&ma<140.5)
    	cout << sum<< "\t"<< sum_ori<<endl;
    }
    	cout << sum<< "\t"<< sum_ori<<endl;
return;
*/

//	TFile *fbkge = new TFile ("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/Fit/whatthefuck/bkg_eff.root");
//	TGraph *eff_bkg =  (TGraph*)fbkge->Get("bkgeff_"+chan);

	TFile *fbkge = new TFile ("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/Reduce/CMSSW_7_1_5/src/160726/bkg_reg_eff.root");
	TGraph *eff_bkg =  (TGraph*)fbkge->Get("ggZZ_reg_"+chan);


	TGraph *phase_sin_fix= new TGraph(nbins*2);
	TGraph *phase_cos_fix= new TGraph(nbins*2);

	RooWorkspace w("w");
TGraph* ggZZ_kf [5];

	mean->setRange(100,1500);
	sigma->setRange(0.00005,100.);
	mean->setVal(125);
	sigma->setVal(0.004165);

	TString formu_2nd=" (@0<@1)*(@3+@0*@4+@0*@0*@5 ) + ( @0>=@1 && @0<@2)*(@6+@0*@7+@0*@0*@8) + (@0>=@2)*(@9+@0*@10+@0*@0*@11)";	

	RooArgList formuList_a1;
	RooArgList formuList_a2;
	RooArgList formuList_mean;
	RooArgList formuList_n1;
	RooArgList formuList_n2;
	RooArgList formuList_sigma;
	formuList_a1.add(*mzz);
	formuList_a2.add(*mzz);
	formuList_mean.add(*mzz);
	formuList_n1.add(*mzz);
	formuList_n2.add(*mzz);
	formuList_sigma.add(*mzz);

	RooConstVar* a1_p0_0_2nd[11] ;
	RooConstVar* a2_p0_0_2nd[11] ;
	RooConstVar* mean_p0_0_2nd[11] ;
	RooConstVar* n1_p0_0_2nd[11] ;
	RooConstVar* n2_p0_0_2nd[11] ;
	RooConstVar* sigma_p0_0_2nd[11] ;
	for(int i =0; i<11;i++){
		a1_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_a1_p0_0_2nd",chan.Data(),i),Form("%s_%d_a1_p0_0_2nd",chan.Data(),i),dcbPara_2nd[0][i]);
		a2_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_a2_p0_0_2nd",chan.Data(),i),Form("%s_%d_a2_p0_0_2nd",chan.Data(),i),dcbPara_2nd[1][i]);
		mean_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_mean_p0_0_2nd",chan.Data(),i),Form("%s_%d_mean_p0_0_2nd",chan.Data(),i),dcbPara_2nd[2][i]);
		n1_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_n1_p0_0_2nd",chan.Data(),i),Form("%s_%d_n1_p0_0_2nd",chan.Data(),i),dcbPara_2nd[3][i]);
		n2_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_n2_p0_0_2nd",chan.Data(),i),Form("%s_%d_n2_p0_0_2nd",chan.Data(),i),dcbPara_2nd[4][i]);
		sigma_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_sigma_p0_0_2nd",chan.Data(),i),Form("%s_%d_sigma_p0_0_2nd",chan.Data(),i),dcbPara_2nd[5][i]);

		formuList_a1.add(*a1_p0_0_2nd[i]);
		formuList_a2.add(*a2_p0_0_2nd[i]);
		formuList_mean.add(*mean_p0_0_2nd[i]);
		formuList_n1.add(*n1_p0_0_2nd[i]);
		formuList_n2.add(*n2_p0_0_2nd[i]);
		formuList_sigma.add(*sigma_p0_0_2nd[i]);
	}

	RooFormulaVar* a1_p0_2nd= new RooFormulaVar("a1_p0_2nd"+chan,"a1_p0_2nd"+chan,formu_2nd,formuList_a1);
	RooFormulaVar* a2_p0_2nd= new RooFormulaVar("a2_p0_2nd"+chan,"a2_p0_2nd"+chan,formu_2nd,formuList_a2);
	RooFormulaVar* mean_p0_2nd= new RooFormulaVar("mean_p0_2nd"+chan,"mean_p0_2nd"+chan,formu_2nd,formuList_mean);
	RooFormulaVar* n1_p0_2nd= new RooFormulaVar("n1_p0_2nd"+chan,"n1_p0_2nd"+chan,formu_2nd,formuList_n1);
	RooFormulaVar* n2_p0_2nd= new RooFormulaVar("n2_p0_2nd"+chan,"n2_p0_2nd"+chan,formu_2nd,formuList_n2);
	RooFormulaVar* sigma_p0_2nd= new RooFormulaVar("sigma_p0_2nd"+chan,"sigma_p0_2nd"+chan,formu_2nd,formuList_sigma);


	//	RooFormulaVar *mzz_shift = new RooFormulaVar("mzz_shift","","@0+@1",RooArgList(*mean_p0_2nd,*mzz));

//    RooFormulaVar *n2_p0_up = new RooFormulaVar("n2_p0_up","","@0+0.2*@0",*n2_p0_2nd);
//    RooFormulaVar *n2_p0_dn = new RooFormulaVar("n2_p0_dn","","@0-0.2*@0",*n2_p0_2nd);

    RooFormulaVar *sigma_p0_up = new RooFormulaVar("sigma_p0_up","","@0+0.2*@0",*sigma_p0_2nd);
    RooFormulaVar *sigma_p0_dn = new RooFormulaVar("sigma_p0_dn","","@0-0.2*@0",*sigma_p0_2nd);
	RooDoubleCB dcrReso("dcrReso"+chan,"Double Crystal ball ",*mreco,*mzz,*mean_p0_2nd,*sigma_p0_2nd,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);
	RooDoubleCB dcrReso_up("dcrReso"+chan+"_up","dcb up",*mreco,*mzz,*mean_p0_2nd,*sigma_p0_up,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);
	RooDoubleCB dcrReso_dn("dcrReso"+chan+"_dn","dcb up",*mreco,*mzz,*mean_p0_2nd,*sigma_p0_dn,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);


//	Width_conv *convpdf_spline[5];
//for(int f = 0; f<1;f++){
int f=0;
	TString sysname = strSystTitle[f];
	ggZZ_kf[f] =(TGraph*)fkfactor->Get("sp_kfactor_"+sysname);
	TGraph *effxkf_sig= new TGraph(nbins*2);
	TGraph *effxkf_bkg= new TGraph(nbins*2);

	for(int i =0;i<nbins*2;i++){
		double cva = low+ i*(high-low)/double(nbins)/2.;
		double effval_sig = (effsig[0]+effsig[1]*TMath::Erf( (cva-effsig[2])/effsig[3] ))*(effsig[4]+effsig[5]*cva+effsig[6]*cva*cva+effsig[10]*cva*cva*cva)+effsig[7]*TMath::Gaus(cva,effsig[8],effsig[9]);
//		double effcate = 2.572353e-02;
		double effcate = 0.0387026; 
		if(!cate_vbf)
			effcate = 1-effcate;

		double va_bkg= eff_bkg->Eval(cva)*ggZZ_kf[f]->Eval(cva)*effcate; 	
		double va_sig= effval_sig*ggZZ_kf[f]->Eval(cva)*effcate; 	
//		if(chan=="2e2mu"){
//			va_bkg*=0.95;
//			va_sig*=0.95;
//		}
		effxkf_sig->SetPoint(effxkf_sig->GetN(),cva,va_sig);
		effxkf_bkg->SetPoint(effxkf_bkg->GetN(),cva,va_bkg);
		phase_sin_fix->SetPoint(phase_sin_fix->GetN(),cva,phase_sin->Eval(cva)/2.*1.76);
		phase_cos_fix->SetPoint(phase_cos_fix->GetN(),cva,phase_cos->Eval(cva)/2.*1.76);
	}
	effxkf_sig->SetName("gghsigeffxkf"+chan+Form("%d",cate_vbf)+strSystTitle[f]);
	effxkf_bkg->SetName("gghbkgeffxkf"+chan+Form("%d",cate_vbf)+strSystTitle[f]);
	phase_sin_fix->SetName("gghsin"+chan+Form("%d",cate_vbf));
	phase_cos_fix->SetName("gghcos"+chan+Form("%d",cate_vbf));
	//phase_cos_fix->Draw("al");
	//phase_sin_fix->Draw("lsame");
	//return;
	RooAbsReal *final_integral;

	if(onshell){
		convpdf_spline=new Width_conv("ggH"+strSystTitle[f], "ggH"+strSystTitle[f],*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso),*phase_cos_fix, *phase_sin_fix, *effxkf_sig,*effxkf_bkg); 
    convpdf_spline->SetNameTitle("ggH","ggH");
		Width_conv convpdf_spline_up("ggH_Res"+chan+"Up", "ggH"+chan+"Up",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso_up),*phase_cos_fix, *phase_sin_fix, *effxkf_sig,*effxkf_bkg); 
		Width_conv convpdf_spline_dn("ggH_Res"+chan+"Down", "ggH"+chan+"Down",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso_dn),*phase_cos_fix, *phase_sin_fix, *effxkf_sig,*effxkf_bkg); 
		w.import(*convpdf_spline,RecycleConflictNodes());
		w.import(convpdf_spline_up,RecycleConflictNodes());
		w.import(convpdf_spline_dn,RecycleConflictNodes());
		final_integral = convpdf_spline->createIntegral(*mreco);
	}
	else{
		Width_conv_offshell *convpdf_spline=new Width_conv_offshell("bggH"+strSystTitle[f], "bggH"+strSystTitle[f],*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso),*phase_cos_fix, *phase_sin_fix, *effxkf_sig,*effxkf_bkg); 
//    convpdf_spline->SetNameTitle("ggH","ggH");
//		Width_conv_offshell convpdf_spline_up("ggH_Res"+chan+"Up", "ggH"+chan+"Up",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso_up),*phase_cos_fix, *phase_sin_fix, *effxkf_sig,*effxkf_bkg); 
//		Width_conv_offshell convpdf_spline_dn("ggH_Res"+chan+"Down", "ggH"+chan+"Down",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso_dn),*phase_cos_fix, *phase_sin_fix, *effxkf_sig,*effxkf_bkg); 
//		w.import(*convpdf_spline,RecycleConflictNodes());
//		w.import(convpdf_spline_up,RecycleConflictNodes());
//		w.import(convpdf_spline_dn,RecycleConflictNodes());
		final_integral = convpdf_spline->createIntegral(*mreco);
	}
//else{
//      strSystTitle[f].ReplaceAll("dn","Down");
//      strSystTitle[f].ReplaceAll("up","Up");
//      convpdf_spline[f]->SetNameTitle("ggH_"+strSystTitle[f],"ggH_"+strSystTitle[f]);
//}

	mean->setVal(125);

	//pdf_ggzz.plotOn(frame_mzz);
	pdf_x.plotOn(frame_mzz);

//	sigma->setVal(5.);
//	convpdf_spline.plotOn(frame);
//	sigma->setVal(1.);
//	convpdf_spline.plotOn(frame,LineColor(2));
	sigma->setVal(0.004);
//	convpdf_spline[f]->plotOn(frame,LineColor(f+1));
	//frame->Draw();
////frame_mzz->Draw();

	r->setVal(0);
	double bexp =  final_integral->getVal();

	RooConstVar *bkg_integral= new RooConstVar("bkg_integral"+chan+Form("%d",cate_vbf)+strSystTitle[f],"",bexp);

	mean->setVal(125);
	//ROOT::Math::Interpolator inter(200, ROOT::Math::Interpolation::kCSPLINE);
	sigma->setRange(0.00005,100.);

	TH2F *hint ; 
	TH2F *hsig ;
	if(!onshell){	
		//hint= new TH2F("hint","",101,99.75,150.25,101,-0.0005,0.1005);
		//hsig= new TH2F("hsig","",101,99.75,150.25,101,-0.0005,0.1005);
		hint= new TH2F("hint","",101,119.95,130.05,101,-0.0005,0.1005);
		hsig= new TH2F("hsig","",101,119.95,130.05,101,-0.0005,0.1005);
	//	hint= new TH2F("hint","",11,100.25,150.25,11,-0.0005,0.1005);
	//	hsig= new TH2F("hsig","",11,100.25,150.25,11,-0.0005,0.1005);
	}
	else{
		hint= new TH2F("hint","",101,119.95,130.05,101,-0.025,5.025);
		hsig= new TH2F("hsig","",101,119.95,130.05,101,-0.025,5.025);
		//hint= new TH2F("hint","",101,122.475,127.525,101,-0.025,5.025);
		//hsig= new TH2F("hsig","",101,122.475,127.525,101,-0.025,5.025);
	//	hint= new TH2F("hint","",11,122.475,127.525,11,-0.025,5.025);
	//	hsig= new TH2F("hsig","",11,122.475,127.525,11,-0.025,5.025);
	}

	for(int i = 0;i <101;i++){
		if(i%10==0)
			cout<<i<<endl;
		for(int j = 0;j <101;j++){
			double mv=hint->GetXaxis()->GetBinCenter(i+1);
			double sv=hint->GetYaxis()->GetBinCenter(j+1);
			sigma->setVal(sv);
			mean->setVal(mv);
			//sigma->setVal(0.001*(i+1));
			//mean->setVal(100+0.5*(j+1));
			r->setVal(1);
			double sbi =  final_integral->getVal();
			r->setVal(4);
			double sbi2 =  final_integral->getVal();
			double sexp = ((sbi2-sbi*2)+bexp)/2.;
			double iexp = sbi -sexp -bexp; 
			//			cout << sigma->getVal()<<"\t"<< mean->getVal()<<endl;
			//	float integral = final_integral->getVal();
			hint->Fill(mean->getVal(),sigma->getVal(),iexp);
			hsig->Fill(mean->getVal(),sigma->getVal(),sexp);
		}
	}

	RooDataHist* hinthist= new RooDataHist("hinthist"+chan+Form("%d",cate_vbf)+strSystTitle[f],"hinthist"+chan+Form("%d",cate_vbf)+strSystTitle[f],RooArgSet(*mean,*sigma),hint);
	RooHistFunc *hintfunc = new RooHistFunc("hintfunc"+chan+Form("%d",cate_vbf)+strSystTitle[f],"",RooArgSet(*mean,*sigma),*hinthist);
	Width_integral inter_intergral ("int_integral"+chan+Form("%d",cate_vbf)+strSystTitle[f],"",*mean,*sigma,RooArgList(*hintfunc));

	RooDataHist* hsighist= new RooDataHist("hsighist"+chan+Form("%d",cate_vbf)+strSystTitle[f],"hsighist"+chan+Form("%d",cate_vbf)+strSystTitle[f],RooArgSet(*mean,*sigma),hsig);
	RooHistFunc *hsigfunc = new RooHistFunc("hsigfunc"+chan+Form("%d",cate_vbf)+strSystTitle[f],"",RooArgSet(*mean,*sigma),*hsighist);
	Width_integral sig_intergral ("sig_integral"+chan+Form("%d",cate_vbf)+strSystTitle[f],"",*mean,*sigma,RooArgList(*hsigfunc));
	RooFormulaVar overall_intergral("ggH_norm"+strSystTitle[f],"","@0*@2+ @1 + sqrt(@2)*@3",RooArgList(sig_intergral, *bkg_integral, *r,inter_intergral));
	if(f==0)
    overall_intergral->SetNameTitle("ggH_norm","ggH_norm");

		Width_conv_integral *convpdf_spline_integral=new Width_conv_integral("ggH", "ggH",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso,*hsigfunc,*hintfunc),*phase_cos_fix, *phase_sin_fix, *effxkf_sig,*effxkf_bkg); 
		Width_conv_integral convpdf_spline_integral_up("ggH_Res"+chan+"Up", "ggH"+chan+"Up",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso_up,*hsigfunc,*hintfunc),*phase_cos_fix, *phase_sin_fix, *effxkf_sig,*effxkf_bkg); 
		Width_conv_integral convpdf_spline_integral_dn("ggH_Res"+chan+"Down", "ggH"+chan+"Down",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso_dn,*hsigfunc,*hintfunc),*phase_cos_fix, *phase_sin_fix, *effxkf_sig,*effxkf_bkg); 
		w.import(*convpdf_spline_integral,RecycleConflictNodes());
		w.import(convpdf_spline_integral_up,RecycleConflictNodes());
		w.import(convpdf_spline_integral_dn,RecycleConflictNodes());

	mean->setVal(125);
	sigma->setVal(0.004);
	cout<< sig_intergral.getVal()<<"\t"<< bkg_integral->getVal()<<"\t"<<inter_intergral.getVal()<<endl;

	overall_intergral.plotOn(frame_width);
//	frame_width->Draw();

	mzz->setConstant(0);
	mean->setConstant(0);
	sigma->setConstant(0);
	mreco->setConstant(0);

	mreco->setRange(low_reco,high_reco);
	mean->setVal(125);
	sigma->setVal(0.004);
	r->setVal(1);

	w.import(overall_intergral,RecycleConflictNodes());

	//frame->Draw();
	//return;
	mreco->Print("v");
	cout << mreco->getBins();

	TString filename = "workspace125/hzz4l_"+chan+Form("%dS_13TeV.input_func.root",cate_vbf);
	if(onshell)
		filename = "workspace125_onshell/hzz4l_"+chan+Form("%dS_13TeV.input_func.root",cate_vbf);
	TFile *fwor=new TFile(filename,"recreate");	
	fwor->cd();
	w.Write();
	fwor->Close();
	return;
}

void clean(TString chan="4e",bool cat_vbf=0, bool onshell=0){
	gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
	gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
	gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");
	dosomething(chan,cat_vbf,onshell);
//	dosomething("2e2mu",0,0);
//	dosomething("4e",0,0);
//	dosomething("4mu",0,0);
//	dosomething("2e2mu",1,0);
//	dosomething("4e",1,0);
//	dosomething("4mu",1,0);

}
