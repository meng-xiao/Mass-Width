//#include "RooKeysPdf.h"

using namespace RooFit;
void convpdf_vbf(){
	gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
	gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
	gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");
//  	const double low=100.5;
//    const double high=1500.5;
//	const int nbins=1400;
	bool onshell = 1;
  	const double low=100.5;
    const double high=151.5;
		const int nbins=500;

	RooRealVar* mzz = new RooRealVar("ZZMass","M_{ZZ} (GeV)",125,low,high);
	RooPlot* frame_mzz= mzz->frame(Title("Z mass")) ;

	TChain *ggzz = new TChain ("SelectedTree");
	ggzz->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/rootfiles/vbf/phantom_bkg.root");
	TChain *ggzz_4e= new TChain ("SelectedTree");
	ggzz_4e->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/rootfiles/vbf/phantom_bkg_4e.root");

   float ZZMass;
   float ZZMass_4e;

    ggzz->SetBranchAddress("GenHMass",&ZZMass);
    ggzz_4e->SetBranchAddress("GenHMass",&ZZMass_4e);


	//float ZZMass;
//	 TTree *cuttree = ggzz->CloneTree(0);
//	 TTree *cuttree_4e = ggzz_4e->CloneTree(0);

	TTree *cuttree =new TTree("SelectedTree","SelectedTree");
	TBranch *zzmass = cuttree->Branch("ZZMass",&ZZMass,"ZZMass/F");

	TTree *cuttree_4e=new TTree("SelectedTree","SelectedTree");
	TBranch *zzmass_4e= cuttree_4e->Branch("ZZMass",&ZZMass_4e,"ZZMass/F");

//	for(int i =0;i<200000;i++){
	for(int i =0;i<500000;i++){
		ggzz->GetEntry(i);
		ggzz_4e->GetEntry(i);
		if(ZZMass>=low&&ZZMass<=high)
		cuttree->Fill();
		if(ZZMass_4e>=low&&ZZMass_4e<=high)
		cuttree_4e->Fill();
	}
	double rho =0.6;
	if(onshell)
		rho=3.;
	RooDataSet bkgdata ("bkgdata","",cuttree,*mzz);
	RooKeysPdf pdf_bkg("vbfpdfbkg_2e2mu","",*mzz,bkgdata,RooKeysPdf::MirrorBoth,rho);

	RooDataSet bkgdata_4e ("bkgdata_4e","",cuttree_4e,*mzz);
	RooKeysPdf pdf_bkg_4e("vbfpdfbkg_4e","",*mzz,bkgdata_4e,RooKeysPdf::MirrorBoth,rho);

	bkgdata.plotOn(frame_mzz,Binning(50));
	pdf_bkg.plotOn(frame_mzz);
	bkgdata_4e.plotOn(frame_mzz,Binning(50),MarkerColor(2));
	pdf_bkg_4e.plotOn(frame_mzz,LineColor(2));
	frame_mzz->Draw();
	RooWorkspace w("w");
	w.import(pdf_bkg,RecycleConflictNodes());
	w.import(pdf_bkg_4e,RecycleConflictNodes());
	TFile *f=new TFile("vbfpdfs_onshell.root","recreate");	
	f->cd();
//	pdf_bkg.Write();
	w.Write();
	f->Close();
	return;
}
