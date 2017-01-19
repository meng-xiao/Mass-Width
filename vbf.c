/** Things to remember later
 * mzz range defined by ggZZ workspace, create 2 versions for onshell and offshell
 * give different ggzz RooKeysPdf name for 2e2mu and 4e
 *
 ****/
using namespace std;
using namespace RooFit;

void dosomething(TString chan ="2e2mu",bool cate_vbf =false,bool onshell=false){

	double lumi = 12.9;
	double ggzz_xsec = 0.488;
	double x_xsec = 0.312*1.108; 

	if(onshell){
		if(chan=="2e2mu"){
		ggzz_xsec = 0.1009;
		x_xsec = 0.275*1.108;; 
		}
		else{
			ggzz_xsec = 0.0473;
			x_xsec = 0.149*1.108;;
		}
	}
		else{
			if(chan!="2e2mu"){
			ggzz_xsec=0.239;
			x_xsec = 0.167*1.108;; 
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
//// 4e
//     -4.409236E+00, 4.628614E+00,-8.559665E+01,1.271156E+02,1.593683E+00,1.585147E-03,-6.754524E-07,1.789210E-02,1.636138E+02,4.799723E+01,8.811266E-11,
////    4mu
//      -4.362745E+00,4.659031E+00,-9.937283E+01,1.333036E+02,1.685288E+00,1.736186E-03,-8.699309E-07,2.834451E-02,1.637809E+02,3.856606E+01, 1.387049E-10,
////    2e2mu
//     -4.410208E+00, 4.628288E+00, -8.355707E+01, 1.227452E+02, 1.952608E+00, 1.877126E-03, -8.679733E-07, 1.838955E-02, 1.717942E+02, 2.085859E+01, 1.263515E-10
//4e
//-4.41139,4.62286,-46.1825,100.856,1.53391,0.0013975,-6.15721e-07,0.0147732,192.54,13.5543,8.38385e-11,
////4mu
//-4.42781,4.61174,-77.899,111.718,2.43934,0.00247607,-1.24247e-06,0.0218535,184.11,10,1.97244e-10,
////2e2mu
//-4.39108,4.64497,-153.433,172.265,1.53766,0.00137698,-6.3804e-07,0.0474389,160,37.524,9.24412e-11
//
////VBF->0p_4e_param_0:
//-4.44635,4.59326,-59.0953,102.861,2.17093,0.00183632,-7.35685e-07,0.0187474,188.212,13.7126,8.75885e-11,
////VBF->0p_4mu_param_0:
//-4.39564,4.64072,-244.111,226.617,1.9519,0.00173683,-7.75353e-07,0.0829745,160,46.8172,1.06192e-10,
////VBF->0p_2e2mu_param_0:
//-4.40817,4.62868,-139.383,159.416,1.76159,0.00160184,-7.16129e-07,0.0416353,160,39.5239,9.86506e-11

//VBF->0p_4e_param_0:
-4.47849,4.5637,-80.8308,107.137,3.38191,0.00298955,-1.33756e-06,0.0179383,188.843,13.8356,1.83965e-10,
//VBF->0p_4mu_param_0:
-4.44681,4.59543,-87.5324,114.192,3.22478,0.00318981,-1.53157e-06,0.0216856,182.532,11.7483,2.26166e-10,
//VBF->0p_2e2mu_param_0:
-4.43582,4.60461,-196.311,191.783,2.22227,0.00204376,-9.41486e-07,0.0718321,148.027,50.0609,1.35394e-10
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
	const int nbinsarr[2]={1500,100};

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
	TFile *fpdfbkg = new TFile("vbfpdfs"+ap+".root");
	RooWorkspace *wbkg =( RooWorkspace*)fpdfbkg->Get ("w"); 

	//RooRealVar* mzz = new RooRealVar("ZZMass","M_{ZZ} (GeV)",125,low,high);
	RooRealVar* mzz = wbkg->var("ZZMass"); 
	RooRealVar* mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",125,low_reco,high_reco);
	RooRealVar* mdiff= new RooRealVar("mdiff","M_{ZZ} (GeV)",125,low_reco,high_reco);

	RooRealVar *r= new RooRealVar("r","signal strength",1.,0.000,1000);
	RooRealVar *rvbf_ggh= new RooRealVar("rvbf_ggh","rvb_ggh",1.,0.000,1000);
	RooFormulaVar *rvbf= new RooFormulaVar("rvbf","@0*@1",RooArgSet(*r,*rvbf_ggh));

	RooRealVar* mean = new RooRealVar("mean_pole","mean_pole",125,100,180);
	RooRealVar* sigma= new RooRealVar("sigma_pole","sigma_pole",0.00418,0.,10.);

	RooConstVar* mean_125= new RooConstVar("mean_125","mean_125",125);
	RooConstVar* sigma_125= new RooConstVar("sigma_125","sigma_125",0.00407);

	RooPlot* frame= mreco->frame(low_reco,high_reco) ;
	//	RooPlot* frame= mreco->frame(150,250) ;
	RooPlot* frame_mzz= mzz->frame(Title("Z mass")) ;
	RooPlot* frame_width= sigma->frame(Title("width")) ;
	RooPlot* frame_mean= mean->frame(Title("mean")) ;

	TFile *flo=new TFile("width_new_spline.root","read");
	TString chn = "2e2mu";
	if(chan!="2e2mu")
		chn="4e";
	TSpline3 *lo=(TSpline3*) flo->Get("br_"+chn);

	SplinePdf par2_int ("vbfpar2_int"+chan+Form("%d",cate_vbf),"",*mzz,*mean,*sigma,*lo);
	RooRealVar m_gauss ("m_gauss","",125);
	RooRealVar w_gauss ("w_gauss","",0.004);
	RooGaussian gauss("gauss","",*mzz,m_gauss,w_gauss);

	TString pdfn = "2e2mu";
	if(chan!="2e2mu")
		pdfn = "4e";
	RooKeysPdf *pdfbkg = wbkg->pdf("vbfpdfbkg_"+pdfn);


	RooConstVar *ggzznorm= new RooConstVar("vbfbkgnorm"+chan+Form("%d",cate_vbf),"",lumi*ggzz_xsec);
	RooExtendPdf pdf_ggzz("vbfpdf_bkg"+chan+Form("%d",cate_vbf),"vbfpdf_bkg"+chan+Form("%d",cate_vbf),*pdfbkg,*ggzznorm);
//	pdf_ggzz.plotOn(frame_mzz);
//	frame_mzz->Draw();
//	return;

	RooConstVar *xnorm= new RooConstVar("vbfxnorm"+chan+Form("%d",cate_vbf),"",lumi*x_xsec);
	RooExtendPdf pdf_x("vbfpdf_x"+chan+Form("%d",cate_vbf),"vbfpdf_x"+chan+Form("%d",cate_vbf),par2_int,*xnorm);


  TFile *fphase_noweight=new TFile("fgraph_vbf_phase.root");
  TGraph *cosfunc = (TGraph*)fphase_noweight->Get("cos");
  TGraph *sinfunc = (TGraph*)fphase_noweight->Get("sin");


//	TFile *fbkge = new TFile ("VBFbkg_eff.root");
//	TGraph *eff_bkg =  (TGraph*)fbkge->Get("bkgeff_"+chan);
	
	TFile *fbkge = new TFile ("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/Reduce/CMSSW_7_1_5/src/160726/bkg_reg_eff.root");
	TGraph *eff_bkg =  (TGraph*)fbkge->Get("VBF_reg_"+chan);


	TGraph *effxkf_sig= new TGraph(nbins);
	TGraph *effxkf_bkg= new TGraph(nbins);

  TFile *fxsec=new TFile("xsec.root");
  TGraph *vbfxs= (TGraph*)fxsec->Get("vbf");
  TGraph *whxs= (TGraph*)fxsec->Get("wh");
  TGraph *zhxs= (TGraph*)fxsec->Get("zh");


	for(int i =0;i<nbins;i++){
		double cva = low+ i*(high-low)/double(nbins);
		double effval_sig = (effsig[0]+effsig[1]*TMath::Erf( (cva-effsig[2])/effsig[3] ))*(effsig[4]+effsig[5]*cva+effsig[6]*cva*cva+effsig[10]*cva*cva*cva)+effsig[7]*TMath::Gaus(cva,effsig[8],effsig[9]);
    double m4l = cva;
    //if (m4l > 900.) m4l = 900.;
    //double effcate = 3.844411e-01 + 2.755289e-04*m4l + -5.098494e-07*m4l*m4l + 2.491489e-10*m4l*m4l*m4l;
		double effcate_sig = 0.476383 - 4.94232e-05*m4l;
		double effcate_bkg = (m4l < 200 ? (3.849250e-01 + 2.677220e-03*(m4l-200)) : (m4l < 300 ? 3.849250e-01+1.520690e-04*(m4l-200) + -1.379680e-07*(m4l-200)*(m4l-200) : +1.426950e+02/(m4l+3.971120e+02) - 1.426950e+02/(300+3.971120e+02) + 3.849250e-01+1.520690e-04*100 + -1.379680e-07*100*100 )) ;
		double vbffrac=1;
		if(cva<2000)
				vbffrac = vbfxs->Eval(cva)*0.955/(whxs->Eval(cva)*0.654+zhxs->Eval(cva)*0.669+vbfxs->Eval(cva)*0.955);
		if(cate_vbf){
			effcate_sig= effcate_sig*vbffrac;
		}
   	else{ 
       effcate_sig= (1-effcate_sig)*vbffrac + (1-vbffrac)*1./0.7;
			 effcate_bkg = 1-effcate_bkg;
		}
		double va_bkg= eff_bkg->Eval(cva)*effcate_bkg; 	
		double va_sig= effval_sig*effcate_sig; 	
//		if(chan=="2e2mu"){
//			va_bkg*=0.95;
//			va_sig*=0.95;
//		}
		effxkf_sig->SetPoint(effxkf_sig->GetN(),cva,va_sig);
		effxkf_bkg->SetPoint(effxkf_bkg->GetN(),cva,va_bkg);
	}
	effxkf_sig->SetName("vbfsigeffxkf"+chan+Form("%d",cate_vbf));
	effxkf_bkg->SetName("vbfbkgeffxkf"+chan+Form("%d",cate_vbf));

	//effxkf_sig->Draw("al");
	//effxkf_bkg->Draw("lsame");
	//return;

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

	//a1_p0_2nd->plotOn(frame_mzz);
	//a2_p0_2nd->plotOn(frame_mzz,LineColor(6));
	//n1_p0_2nd->plotOn(frame_mzz,LineColor(kCyan));
	//n2_p0_2nd->plotOn(frame_mzz,LineColor(kOrange));
	//mean_p0_2nd->plotOn(frame_mzz,LineColor(2));
	//sigma_p0_2nd->plotOn(frame_mzz,LineColor(3));
	//frame_mzz->Draw();
	//return;

	//	RooFormulaVar *mzz_shift = new RooFormulaVar("mzz_shift","","@0+@1",RooArgList(*mean_p0_2nd,*mzz));


    //RooFormulaVar *n2_p0_up = new RooFormulaVar("n2_p0_up","","@0+0.2*@0",*n2_p0_2nd);
    //RooFormulaVar *n2_p0_dn = new RooFormulaVar("n2_p0_dn","","@0-0.2*@0",*n2_p0_2nd);
    RooFormulaVar *sigma_p0_up = new RooFormulaVar("sigma_p0_up","","@0+0.2*@0",*sigma_p0_2nd);
    RooFormulaVar *sigma_p0_dn = new RooFormulaVar("sigma_p0_dn","","@0-0.2*@0",*sigma_p0_2nd);

	RooWorkspace w("w");
	RooDoubleCB dcrReso("dcrReso"+chan,"Double Crystal ball"+chan,*mreco,*mzz,*mean_p0_2nd,*sigma_p0_2nd,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);
	RooDoubleCB dcrReso_up("dcrReso"+chan+"_up","dcb up"+chan,*mreco,*mzz,*mean_p0_2nd,*sigma_p0_up,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);
	RooDoubleCB dcrReso_dn("dcrReso"+chan+"_dn","dcb up"+chan,*mreco,*mzz,*mean_p0_2nd,*sigma_p0_dn,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);
	RooAbsReal *final_integral; 

if(onshell){
	Width_conv convpdf_spline("qqH", "qqH",*mreco, *mean, *sigma, *rvbf, RooArgList(pdf_x, pdf_ggzz,dcrReso),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
//	convpdf_spline.plotOn(frame);
	w.import(convpdf_spline,RecycleConflictNodes());
	Width_conv convpdf_spline_up("qqH_Res"+chan+"Up", "qqH_Res"+chan+"Up",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso_up),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
	Width_conv convpdf_spline_dn("qqH_Res"+chan+"Down", "qqH_Res"+chan+"Down",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso_dn),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
	w.import(convpdf_spline_up,RecycleConflictNodes());
	w.import(convpdf_spline_dn,RecycleConflictNodes());
	final_integral = convpdf_spline.createIntegral(*mreco);
}
else{
	Width_conv_offshell convpdf_spline("bqqH", "bqqH",*mreco, *mean, *sigma, *rvbf, RooArgList(pdf_x, pdf_ggzz,dcrReso),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
//	convpdf_spline.plotOn(frame);
//	w.import(convpdf_spline,RecycleConflictNodes());
//	Width_conv_offshell convpdf_spline_up("qqH_Res"+chan+"Up", "qqH_Res"+chan+"Up",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso_up),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
//	Width_conv_offshell convpdf_spline_dn("qqH_Res"+chan+"Down", "qqH_Res"+chan+"Down",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso_dn),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
	final_integral = convpdf_spline.createIntegral(*mreco);
}

	mean->setVal(125);
	mreco->setVal(126);

	//pdf_ggzz.plotOn(frame_mzz);
	//pdf_x.plotOn(frame_mzz);

//	sigma->setVal(5.);
//	convpdf_spline.plotOn(frame);
//	sigma->setVal(1.);
//	convpdf_spline.plotOn(frame,LineColor(2));
//	sigma->setVal(0.004);
//	convpdf_spline.plotOn(frame,LineColor(8));
//	frame->Draw();
//	//frame_mzz->Draw();
//	return;

//	RooAbsReal *final_integral = convpdf_spline.createIntegral(*mreco);
	r->setVal(0);
	double bexp =  final_integral->getVal();
	cout << bexp<<endl;

	RooConstVar *bkg_integral= new RooConstVar("bkg_integral"+chan+Form("%d",cate_vbf),"",bexp);

	mean->setVal(125);
	//ROOT::Math::Interpolator inter(200, ROOT::Math::Interpolation::kCSPLINE);
	sigma->setRange(0.00005,100.);

	TH2F *hint ; 
	TH2F *hsig ;
	if(!onshell){	
		hint= new TH2F("hint","",101,119.95,130.05,101,-0.0005,0.1005);
		hsig= new TH2F("hsig","",101,119.95,130.05,101,-0.0005,0.1005);
	}
	else{
		hint= new TH2F("hint","",101,119.95,130.05,101,-0.025,5.025);
		hsig= new TH2F("hsig","",101,119.95,130.05,101,-0.025,5.025);
		//hint= new TH2F("hint","",101,122.475,127.525,101,-0.025,5.025);
		//hsig= new TH2F("hsig","",101,122.475,127.525,101,-0.025,5.025);
	}
	//		double xi[201]; 
	//		double yi[201] ;

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
			//cout<< sv<<"\t"<<mv<< "\t"<<sexp<< "\t"<<iexp<<endl;
			//			cout << sigma->getVal()<<"\t"<< mean->getVal()<<endl;
			//	float integral = final_integral->getVal();
			hint->Fill(mean->getVal(),sigma->getVal(),iexp);
			hsig->Fill(mean->getVal(),sigma->getVal(),sexp);
			//if(sigma->getVal()>1)
			//cout << sexp << endl;
			////		xi[i] = sigma->getVal();
			////		yi[i] = integral; 
		}
	}
	RooDataHist* hinthist= new RooDataHist("vbfhinthist"+chan+Form("%d",cate_vbf),"vbfhinthist"+chan+Form("%d",cate_vbf),RooArgSet(*mean,*sigma),hint);
	RooHistFunc *hintfunc = new RooHistFunc("vbfhintfunc"+chan+Form("%d",cate_vbf),"",RooArgSet(*mean,*sigma),*hinthist);
	Width_integral inter_intergral ("vbfint_integral"+chan+Form("%d",cate_vbf),"",*mean,*sigma,RooArgList(*hintfunc));

	RooDataHist* hsighist= new RooDataHist("vbfhsighist"+chan+Form("%d",cate_vbf),"hsighist"+chan+Form("%d",cate_vbf),RooArgSet(*mean,*sigma),hsig);
	RooHistFunc *hsigfunc = new RooHistFunc("vbfhsigfunc"+chan+Form("%d",cate_vbf),"",RooArgSet(*mean,*sigma),*hsighist);
	Width_integral sig_intergral ("vbfsig_integral"+chan+Form("%d",cate_vbf),"",*mean,*sigma,RooArgList(*hsigfunc));
	RooFormulaVar overall_intergral("qqH_norm","","@0*@2+ @1 + sqrt(@2)*@3",RooArgList(sig_intergral, *bkg_integral, *rvbf,inter_intergral));
	mean->setVal(125);
	sigma->setVal(0.004);
	cout<< sig_intergral.getVal()<<"\t"<< bkg_integral->getVal()<<"\t"<<inter_intergral.getVal()<<endl;

	Width_conv_integral convpdf_spline_integral("qqH", "qqH",*mreco, *mean, *sigma, *rvbf, RooArgList(pdf_x, pdf_ggzz,dcrReso,*hsigfunc,*hintfunc),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
	Width_conv_integral convpdf_spline_integral_up("qqH_Res"+chan+"Up", "qqH_Res"+chan+"Up",*mreco, *mean, *sigma, *rvbf, RooArgList(pdf_x, pdf_ggzz,dcrReso_up,*hsigfunc,*hintfunc),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
	Width_conv_integral convpdf_spline_integral_dn("qqH_Res"+chan+"Down", "qqH_Res"+chan+"Down",*mreco, *mean, *sigma, *rvbf, RooArgList(pdf_x, pdf_ggzz,dcrReso_dn,*hsigfunc,*hintfunc),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 

	w.import(convpdf_spline_integral,RecycleConflictNodes());
	w.import(convpdf_spline_integral_up,RecycleConflictNodes());
	w.import(convpdf_spline_integral_dn,RecycleConflictNodes());
	overall_intergral.plotOn(frame_width);
	frame_width->Draw();

	mzz->setConstant(0);
	mean->setConstant(0);
	sigma->setConstant(0);
	mreco->setConstant(0);

	mreco->setRange(low_reco,high_reco);
	mean->setVal(125);
	sigma->setVal(0.004);
	r->setVal(1);

	convpdf_spline.plotOn(frame,LineColor(2));
	frame->Draw();

	w.import(overall_intergral,RecycleConflictNodes());

	mreco->Print("v");
	cout << mreco->getBins();

	TString filename = "workspace125/hzz4l_"+chan+Form("%dS_13TeV.input_func_vbf.root",cate_vbf);
	if(onshell)
		filename = "workspace125_onshell/hzz4l_"+chan+Form("%dS_13TeV.input_func_vbf.root",cate_vbf);
	TFile *f=new TFile(filename,"recreate");	
	f->cd();
	w.Write();
	f->Close();
	return;
}

void vbf(TString chan="4e",bool cat_vbf=1, bool onshell=0){
	gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
	gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
	gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");
	dosomething(chan,cat_vbf,onshell);
}
