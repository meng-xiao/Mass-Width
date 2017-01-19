//#include "RooKeysPdf.h"

using namespace RooFit;
void convpdf(){
	gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
	gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
	gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");
  	//const double low=100.5;
    //const double high=1500.5;
		//const int nbins=1400;
		bool onshell = 1;

  	const double low=100.5;
    const double high=150.5;
		const int nbins=500;
	RooRealVar* mzz = new RooRealVar("ZZMass","M_{ZZ} (GeV)",125,low,high);
	RooPlot* frame_mzz= mzz->frame(Title("Z mass")) ;

	TChain *ggzz = new TChain ("SelectedTree");
	ggzz->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/rootfiles/ggZZ_Bkg_xcheck.root");
	TChain *ggzz_4e= new TChain ("SelectedTree");
	ggzz_4e->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/rootfiles/ggZZ_4e.root");

	//float ZZMass;
	 TTree *cuttree = ggzz->CloneTree(0);
	 TTree *cuttree_4e = ggzz_4e->CloneTree(0);

//	TTree *cuttree =new TTree("SelectedTree","SelectedTree");
//	TBranch *zzmass = cuttree->Branch("ZZMass",&ZZMass,"ZZMass/F");

	//for(int i =0;i<100000;i++){
	for(int i =0;i<900000;i++){
		ggzz->GetEntry(i);
		ggzz_4e->GetEntry(i);
		cuttree->Fill();
		cuttree_4e->Fill();
	}

	double rho = 0.9;
	if(onshell)
		rho = 4.;
	RooDataSet bkgdata ("bkgdata","",cuttree,*mzz);
	RooKeysPdf pdf_bkg("pdfbkg_2e2mu","",*mzz,bkgdata,RooKeysPdf::MirrorBoth,rho);

	RooDataSet bkgdata_4e ("bkgdata_4e","",cuttree_4e,*mzz);
	RooKeysPdf pdf_bkg_4e("pdfbkg_4e","",*mzz,bkgdata_4e,RooKeysPdf::MirrorBoth,rho);

	bkgdata.plotOn(frame_mzz,Binning(50));
	pdf_bkg.plotOn(frame_mzz);
	bkgdata_4e.plotOn(frame_mzz,Binning(50),MarkerColor(2));
	pdf_bkg_4e.plotOn(frame_mzz,LineColor(2));
	frame_mzz->Draw();
	RooWorkspace w("w");
	w.import(pdf_bkg,RecycleConflictNodes());
	w.import(pdf_bkg_4e,RecycleConflictNodes());
	TFile *f=new TFile("pdfs_onshell.root","recreate");	
	f->cd();
//	pdf_bkg.Write();
	w.Write();
	f->Close();
	return;
}
