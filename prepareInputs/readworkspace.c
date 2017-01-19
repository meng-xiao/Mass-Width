using namespace RooFit;
using namespace RooStats;
void readworkspace(){
	gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
	 gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
	 gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");

TFile *f=new TFile("batch_offshell/hzz4l_all_13TeV.root");
//TFile *f=new TFile("vbfpdfs.root");
//TFile *f=new TFile("pdfs.root");
RooWorkspace *w=(RooWorkspace*)f->Get("w");
w->Print("v");
return;
w->exportToCint("wo");
//wo::ggH.Print("v");

//return;
  RooRealVar *mreco= w->var("mreco");
  RooPlot *frame = mreco->frame();
  RooPlot *frame_mean= wo::mean_pole->frame();//(122.5,127.5);
  RooPlot *frame_sigma= wo::sigma_pole->frame(0.,0.1);//(0,5.);
  wo::mean_pole.setRange(100,1500);
  wo::sigma_pole.setRange(0,5);
  RooPlot *frame_mzz= wo::ZZMass->frame();

//	wo::vbfpdfbkg_2e2mu.plotOn(frame_mzz);
//	wo::vbfpdfbkg_4e.plotOn(frame_mzz,LineColor(2));
//	wo::pdfbkg_2e2mu.plotOn(frame_mzz);
//	wo::pdfbkg_4e.plotOn(frame_mzz,LineColor(2));
//	frame_mzz->Draw();
//	return:

  wo::mean_pole.setVal(125);
  wo::sigma_pole.setVal(0.1);
//  wo::sigma_pole.setVal(0.004);
  wo::r.setVal(1);
  //wo::ggH.plotOn(frame,Normalization(wo::wo::ggH_norm.getVal(),RooAbsReal::NumEvent));
//////  wo::ggH_Res4muUp.plotOn(frame);
//////  wo::ggH_Res4muDown.plotOn(frame);
  //wo::sigma_pole.setVal(0.004);
  //wo::ggH.plotOn(frame,Normalization(wo::wo::ggH_norm.getVal(),RooAbsReal::NumEvent),LineColor(2));
	//cout << wo::ggH_norm.getVal()<<endl;
	//frame->Draw();
	//return;

  //wo::ggH_qcd_Up.plotOn(frame,LineColor(2));
  //wo::ggH_qcd_Down.plotOn(frame,LineColor(8));
  //wo::ggH_pdf_Up.plotOn(frame,LineColor(2));
  //wo::ggH_pdf_Down.plotOn(frame,LineColor(8));
  //wo::ggH_Res2e2muUp.plotOn(frame,LineColor(2));
  //wo::ggH_Res2e2muDown.plotOn(frame,LineColor(8));
  //frame->Draw();
	//return;

//	wo::hinthist2e2mu0pdf_Up.plotOn(frame_mean,LineColor(2));
//	wo::hinthist2e2mu0pdf_Down.plotOn(frame_mean,LineColor(8));
  wo::ggH_norm.plotOn(frame_mean);
  //wo::ggH_normpdf_Up.plotOn(frame_mean,LineColor(2));
  //wo::ggH_normpdf_Down.plotOn(frame_mean,LineColor(8));
  //wo::ggH_normqcd_Up.plotOn(frame_mean,LineColor(2));
  //wo::ggH_normqcd_Down.plotOn(frame_mean,LineColor(8));
  //wo::sig_integral.plotOn(frame_mean);
  //wo::int_integral.plotOn(frame_mean);

  wo::ggH_norm.plotOn(frame_sigma);
  //wo::ggH_normpdf_Up.plotOn(frame_sigma,LineColor(2));
  //wo::ggH_normqcd_Up.plotOn(frame_sigma,LineColor(2));
  //wo::ggH_normqcd_Down.plotOn(frame_sigma,LineColor(8));
  //wo::ggH_normpdf_Down.plotOn(frame_sigma,LineColor(8));
//  //wo::sig_integral.plotOn(frame_sigma);
//  //wo::int_integral.plotOn(frame_sigma);
  frame_sigma->Draw();
  TCanvas *c2=new TCanvas("c2","",800,600);
  frame_mean->Draw();
 
//wo::ggH.plotOn(frame,LineColor(2),Normalization(wo::overallIntegral.getVal(),RooAbsReal::NumEvent));
return;

wo::sigma_pole.setVal(0.0418);
	TH1F *toys = (TH1F*)wo::ggH.createHistogram("dataset",*mreco,Binning(1500));
	toys->Scale(50./toys->Integral());
	RooDataHist* toyhist= new RooDataHist("toyhist","toyhist",RooArgSet(*mreco),toys);
	toyhist->plotOn(frame,DataError(RooAbsData::None));
wo::sigma_pole.setVal(0.00418);
wo::ggH.plotOn(frame,LineColor(2));
wo::sigma_pole.setVal(0.0418);
wo::ggH.plotOn(frame,LineColor(4));
//wo::r.setVal(10);
//wo::ggH.plotOn(frame,LineColor(8));
frame->Draw("");

}
