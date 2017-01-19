
//#include "HZZ2L2QRooPdfs.cc+"
//#include "HZZ4L_RooHighmass.cc+"
//#include "RooHighmass_conv.cc+"

using namespace RooFit;
using namespace RooStats;
void readWorkspace(){
	gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
	 gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
	 gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");
//	TFile *f=new TFile("hzz4l_13TeV.root");
	TFile *f=new TFile("testConv.root");
//	TFile *f=new TFile("workspaceWithAsimov_gen.root");

//	TFile *f=new TFile("hzz4l_13TeV_onshell.root");
	//TFile *f=new TFile("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/Fit/whatthefuck/testConv.root");
	RooWorkspace *w= (RooWorkspace*) f->Get("w");
	w->Print();
	w->exportToCint("wo");
	cout <<wo::mean_pole->getVal()<<endl;
	cout<<wo::sigma_pole->getVal()<<endl;;
	//wo::CMS_channel->Print("v");
	//wo::model_s->Print("v");
	RooPlot *frame = wo::ZZMass2->frame();
	RooPlot *frame_sig= wo::r->frame(0,5);
	RooPlot *frame_mean= wo::mean_pole->frame();
	RooPlot *frame_sigma= wo::sigma_pole->frame();
	RooRealVar mreco= wo::ZZMass2;
	return;

	//TFile *ftoy = new TFile ("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/Fit/whatthefuck/toys_0.5.root");
	//TFile *ftoy = new TFile ("higgsCombinescan_5.MultiDimFit.mH120.123456.root");
//	TFile *ftoy = new TFile ("higgsCombinescan.GenerateOnly.mH120.123456.root");
//	TFile *ftoy = new TFile ("asimov_onshell.root");
//	RooDataSet *toys = (RooDataSet*) ftoy->Get("asimov");
	RooDataSet *toys = (RooDataSet*) w->data("toys/asimov");
	RooCategory cat("CMS_channel","CMS_channel");
	cat.defineType("ch1",0);
//	RooDataSet *toys_cate= new RooDataSet("asimov_cate","asimov_cate",RooArgSet(mreco,cat),Index(cat),Import("ch1",*toys));
//	TFile *ftoy_cate= new TFile ("asimov_cate.root","recreate");
//	ftoy_cate->cd();
//	toys_cate->Write();
		RooBinning tbins(104,1604) ;
		tbins.addUniform(500,104,154) ;
		tbins.addUniform(400,154,554) ;
		tbins.addUniform(525,554,1604) ;

//	toys->plotOn(frame,Binning(tbins));
//toys->plotOn(frame,Binning(75));
//toys->plotOn(frame,Binning(1500));
//cout << toys->sumEntries()<<endl;
//	RooDataSet *toys = (RooDataSet*) ftoy->Get("smear_0.5");
//	RooDataSet *toys_01 = (RooDataSet*) ftoy->Get("smear_0.1");
//	RooDataSet *toys = (RooDataSet*) ftoy->Get("toys/toy_asimov");
//	RooDataSet *toys = wo::toys/smear_0.5;
//	toys->plotOn(frame);
//	toys_01->plotOn(frame,Binning(1500),MarkerColor(4));
//
	//TFile *ftoy = new TFile("asimov.root");
	////RooDataSet *toy = w->data("toys/asimov_cate");
	//RooDataSet *toy = ftoy->Get("asimov");
	//toy->plotOn(frame,Binning(750));

	wo::mean_pole.setVal(125);
	wo::sigma_pole.setVal(0.004);
	wo::r.setVal(1);
	wo::pdf_binch1.plotOn(frame);
//	frame->Draw();
//	return;
	//wo::r.setVal(0);
	//cout<< wo::n_exp_final_binch1_proc_ggH.getVal()<<endl;	
	//wo::r.setVal(1);
	//cout<< wo::n_exp_final_binch1_proc_ggH.getVal()<<endl;	
	//wo::r.setVal(2);
	//cout<< wo::n_exp_final_binch1_proc_ggH.getVal()<<endl;	
//	wo::n_exp_final_binch1_proc_ggH.plotOn(frame_mean);
//	//wo::hsigfunc.plotOn(frame_mean);
//	//wo::hsigfunc.plotOn(frame_sigma);
//	//wo::mean_pole.setVal(143.5);	
//	//cout<<wo::hsigfunc.getVal()<<endl;
//	//wo::mean_pole.setVal(142.5);	
//	//cout<<wo::hsigfunc.getVal()<<endl;
//
//////	frame_sigma->Draw();
//	frame_mean->Draw();
	//return;
////
	
	double bkgexp = wo::n_exp_binch1_proc_bkg_qqzz.getVal();
	double sigexp = wo::n_exp_final_binch1_proc_ggH.getVal();
	////cout<<bkgexp	<<endl;
	//cout<<sigexp<<endl;
	//wo::r.setVal(0);
	//cout<<	wo::n_exp_final_binch1_proc_ggH.getVal()<<endl;
	//return;
    RooRealVar weight ("weight","",1.,0.,2000.);
    RooDataSet *asimov= new RooDataSet("asimov","",RooArgSet(wo::ZZMass2,weight,cat),"weight"); 
//    TH1F *toyhs = (TH1F*)wo::pdf_binch1.createHistogram("dataset",wo::ZZMass2,Binning(15000,105.5,1605.5));
    TH1F *toyhs = (TH1F*)wo::pdf_binch1.createHistogram("dataset",wo::ZZMass2,Binning(tbins));
		//cout << toyhs->Integral()<<"\t"<<toyhs->Integral("width")<<endl;
		//return;
   double ow = (bkgexp+sigexp)/toyhs->Integral("width");
   //double ow = (bkgexp+sigexp);
//    for(int i =0;i<1500;i++){
    for(int i =0;i<1425;i++){
//    for(int i =0;i<350;i++){
			int bc = i+1;
			double ww = 1.;
			if(i<500)
//				bc=i+1;
				ww=0.1;
			else if(i>=500&&i<900){
//				bc = 495+(i-500)*10+1;
				ww =1.;
			}
			else{ 
			//	bc = 4450+(i-900)*100+1;
				ww =2.;
			}

      wo::ZZMass2.setVal(toyhs->GetBinCenter(bc));
      weight.setVal(ow*ww*toyhs->GetBinContent(bc));

//     wo::ZZMass2.setVal(105+2.*i);
		 
    // weight.setVal(ow*wo::pdf_binch1.getVal(wo::ZZMass2));
     //wo::ZZMass2.setVal(105+2*i);
    // wo::ZZMass2.setVal(105+0.1*i);
 //    weight.setVal(sigexp*wo::shapeSig_ggH_ch1.getVal(wo::ZZMass2)+bkgexp*wo::shapeBkg_bkg_qqzz_ch1.getVal(wo::ZZMass2));
		   cat.setLabel("ch1") ;

     asimov->add(RooArgSet(wo::ZZMass2,weight,cat),weight.getVal());
     cout << wo::ZZMass2.getVal()<<"\t"<< weight.getVal()<<endl;
    }
//		cout << asimov->sumEntries()<<endl;
//   cout<<(bkgexp+sigexp)<<endl;
		asimov->plotOn(frame,Binning(750));
		//asimov->plotOn(frame,Binning(tbins),DataError(RooAbsData::None));
//		asimov->plotOn(frame,DataError(RooAbsData::None));
//	asimov->plotOn(frame,Binning(35));
	wo::pdf_binch1.plotOn(frame);
//	wo::pdf_binch1.fitTo(*asimov);
	frame->Draw();
	//return;
	//wo::pdf_binch1->plotOn(frame);
	//wo::sigma_pole->setVal(0.01);
	//wo::pdf_binch1->plotOn(frame,LineColor(2));
	TFile *fasi = new TFile("asimov.root","recreate");
	//TFile *fasi = new TFile("asimov_onshell.root","recreate");
	fasi->cd();
	asimov->Write();
	fasi->Close();
	return;
	////RooSimultaneous *model_s= w->pdf("model_s");
	////RooAbsPdf *pdf = model_s->getPdf("a1_1");
  //RooRealVar *mreco= w->var("mreco");
	//RooDataSet *toy=wo::pdf_bina1_1->generate(*mreco,1000);
	//toy->plotOn(frame);
	//frame->Draw();
	//wo::pdf_bina1_1->fitTo(*toy);

	/*
//  w->exportToCint("wo");
  RooRealVar *mzz= w->var("ZZMass");
  RooRealVar *signorm= w->var("signorm");
  RooAbsReal *ggH_norm= w->function("ggH_norm");
	RooPlot *frame = mreco->frame();
//	wo::bsi_Int[ZZMass]->plotOn(frame);
//	RooRealVar *norm = (RooRealVar*)w->var("bsi_Int"); 
	RooAbsPdf *pdf= (RooAbsPdf*)w->pdf("ggH"); 
	RooAbsPdf *newconv= (RooAbsPdf*)w->pdf("newconv"); 
//	RooAbsPdf *pdf= (RooAbsPdf*)w->pdf("conv_brut"); 
	RooAbsPdf *pdf_gen= (RooAbsPdf*)w->pdf("bsi_hist"); 
	RooAbsPdf *dcrReso_2nd= (RooAbsPdf*)w->pdf("dcrReso_2nd"); 

	RooAbsReal *mean_p0_2nd = w->function("mean_p0_2nd");
	RooAbsReal *sigma_p0_2nd=w->function("sigma_p0_2nd");
	RooAbsReal *a1_p0_2nd =w->function("a1_p0_2nd");
	RooAbsReal *n1_p0_2nd=w->function("n1_p0_2nd");
	RooAbsReal *a2_p0_2nd=w->function("a2_p0_2nd");
	RooAbsReal *n2_p0_2nd=w->function("n2_p0_2nd");
	a1_p0_2nd->Print();
	pdf->defaultIntegratorConfig()->setEpsRel(1e-7) ;
	pdf->defaultIntegratorConfig()->setEpsAbs(1e-7) ;
//		pdf->defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
	RooPlot *frame_norm= signorm->frame();

	RooArgList *cachepdf =new RooArgList("cachepdf");
	RooArgSet *cachevar =new RooArgSet("cachevar");

	RooRealSumPdf *cond_pdf[100];
	RooRealVar *weight[100];
	RooRealVar *weight_1[100];
//	RooHighmass_conv *newconv = new RooHighmass_conv("newconv","",*mreco,RooArgList(*pdf_gen,*dcrReso_2nd));
//	mreco->setVal(200);
//	cout<<	newconv->getVal()<<endl;;
//	mreco->setVal(250);
//	cout<<newconv->getVal()<<endl;
	for(int i=0;i<1;i++){
		signorm->setVal(i*1.);
		double norm= ggH_norm->getVal();
		cout << norm<<endl;
		pdf->plotOn(frame,LineColor(i+1),Normalization(norm,RooAbsReal::NumEvent));
  	newconv->plotOn(frame,LineColor(i+1),Normalization(norm,RooAbsReal::NumEvent));	
	}
	frame->Draw();
//	gPad->Print("norm_check.png");
//	gPad->Print("norm_check.pdf");
//	gPad->Print("norm_check.root");

	RooArgSet *reso_arg = dcrReso_2nd->getParameters(*mreco);		
		//reso_arg->Print();

	for(int j=0;j<100;j++){
		double mzzval=j*5+200.;
		mzz->setVal(mzzval);
		double gen_weight = pdf_gen->getVal(); 
		weight[j]=new RooRealVar (Form("w%d",j),"",gen_weight);
		weight_1[j]=new RooRealVar (Form("w_1%d",j),"",1);
		cachevar->add(*weight[j]);
		RooRealVar *mean_p0=new RooRealVar("mean_p0","",0);
		RooRealVar *sigma_p0=new RooRealVar("sigmap0","",0);
		RooRealVar *a1_p0=new RooRealVar("a1_p0","",0);
		RooRealVar *n1_p0=new RooRealVar("n1_p0","",0);
		RooRealVar *a2_p0=new RooRealVar("a2_p0","",0);
		RooRealVar *n2_p0=new RooRealVar("n2_p0","",0);
		mean_p0->setVal( mean_p0_2nd->getVal());
		sigma_p0->setVal( sigma_p0_2nd->getVal());
		a1_p0->setVal( a1_p0_2nd->getVal());
		a2_p0->setVal( a2_p0_2nd->getVal());
		n2_p0->setVal( n2_p0_2nd->getVal());
		n1_p0->setVal( n1_p0_2nd->getVal());
		RooDoubleCB dcbtmp("dcbtmp","",*mreco,*mzz,*mean_p0,*sigma_p0,*a1_p0,*n1_p0,*a2_p0,*n2_p0);
		RooAbsPdf *tmppdf = (RooAbsPdf*)dcbtmp.clone(Form("pd%d",j));
//		cond_pdf[j]= new RooRealSumPdf(Form("cond_pdf%d",j), "",dcbtmp, *weight[j]); 
		cachepdf->add(*tmppdf);
		tmppdf->plotOn(frame,LineColor(4),Normalization(gen_weight,RooAbsReal::NumEvent));
	}
	RooRealSumPdf hand("hand","",*cachepdf,*cachevar);
	hand.plotOn(frame,LineColor(4));

	frame->Draw();
	//ggH_norm->plotOn(frame_norm,Range(0.,10.));
	//TCanvas *c2 = new TCanvas("c2","",800,600);
	//frame_norm->Draw();
	
	*/
}
