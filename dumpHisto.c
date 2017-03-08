using namespace RooFit;
using namespace RooStats;

void dumpHisto(){
	gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
	gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
	gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");

//	gROOT->ProcessLine(".x tdrstyle.cc");
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadTopMargin(0.05);	


	TFile *f=new TFile("/afs/cern.ch/user/w/wahung/work/public/CMSSW_7_1_5/src/HiggsAnalysis/CombinedLimit/Mass-Width/workspace125_onshell/hzz4l_2e2mu0S_13TeV.input_func.root");
	RooWorkspace *w= (RooWorkspace*) f->Get("w");

	w->exportToCint("wo");


	const int nbins=140;
	TH1F *hsig[10] ;


	TFile *fnew = new TFile("histo.root","recreate");

	for(int i=0;i<10;i++){
		float mass = 122+0.5*i;
		wo::mean_pole.setVal(mass);
		wo::sigma_pole.setVal(0.004);
		hsig[i]= new TH1F(Form("hsig_%.0f",mass*10),"",nbins,105,140);

		wo::r.setVal(1);

		hsig[i]=(TH1F*)wo::ggH.createHistogram(hsig[i]->GetName(),wo::mreco,Binning(nbins));
		double expect= wo::ggH_norm->getVal();
		hsig[i]->Scale(expect/hsig[i]->Integral());

		fnew->cd();
		hsig[i]->Write(hsig[i]->GetName());
		hsig[i]->SetLineColor(1);
		hsig[i]->SetLineWidth(2);
		if(i==0)
			hsig[i]->Draw();
		else
			hsig[i]->Draw("same");
	}
	gPad->Print("mass.png");
	fnew->Close();
}
