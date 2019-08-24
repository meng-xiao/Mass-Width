
using namespace RooFit;
using namespace RooStats;
void readCombineWorkspace(){
	gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
	 gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
	 gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");
	TFile *f=new TFile("workspace125_onshell_simple/hzz4l_2e2mu0S_13TeV.input_func.root");

	RooWorkspace *w= (RooWorkspace*) f->Get("w");
	w->Print();
	RooRealVar *mreco = (RooRealVar*)w->var("mreco");
	RooRealVar *mean_pole = (RooRealVar*)w->var("mean_pole");
	RooRealVar *sigma_pole = (RooRealVar*)w->var("sigma_pole");
	RooFormulaVar *ggH_norm= (RooFormulaVar*)w->function("ggH_norm");
	RooRealVar *r= (RooRealVar*)w->var("r");
	RooAbsPdf *pdf = (RooAbsPdf*)w->pdf("ggH");
	RooPlot *frame = mreco->frame();

	mean_pole->setVal(125);
	sigma_pole->setVal(1);
	r->setVal(1);

	pdf->plotOn(frame,Normalization(ggH_norm->getVal(),RooAbsReal::NumEvent));

	mean_pole->setVal(125);
	sigma_pole->setVal(0.5);
	r->setVal(1);
	pdf->plotOn(frame,Normalization(ggH_norm->getVal(),RooAbsReal::NumEvent));

	mean_pole->setVal(125);
	sigma_pole->setVal(0.5);
	r->setVal(0);
	pdf->plotOn(frame,Normalization(ggH_norm->getVal(),RooAbsReal::NumEvent));
	frame->Draw();
	gPad->Print("sigshape.png");
	return;

}
