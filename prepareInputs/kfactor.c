#include <algorithm>    // std::max

using namespace std;
void kfactor(){
	TFile *fkfactor = new TFile("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/Fit/whatthefuck/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
	TSpline3* ggZZ_kf[9];
  TString strSystTitle[9] ={
  "Nominal",
  "PDFScaleDn",
  "PDFScaleUp",
  "QCDScaleDn",
  "QCDScaleUp",
  "AsDn",
  "AsUp",
  "PDFReplicaDn",
  "PDFReplicaUp"
  };
 
for(int i =0;i<9;i++){
	ggZZ_kf[i]=(TSpline3*)fkfactor->Get("sp_kfactor_"+strSystTitle[i]);
}

	double xi[5801];	
	double yi_pdf_up[5801];	
	double yi_qcd_up[5801];	
	double yi_pdf_dn[5801];	
	double yi_qcd_dn[5801];	
	double yi[5801];	

	for(int i =0;i<5801;i++){
	double ma = 100.+0.5*i;
	double pdfrep_up = ggZZ_kf[8]->Eval(ma);
	double pdfrep_dn = ggZZ_kf[7]->Eval(ma);

	double pdfscale_up = ggZZ_kf[2]->Eval(ma);
	double pdfscale_dn = ggZZ_kf[1]->Eval(ma);

	double as_up = ggZZ_kf[6]->Eval(ma);
	double as_dn = ggZZ_kf[5]->Eval(ma);

	double qcdscale_up = ggZZ_kf[4]->Eval(ma);
	double qcdscale_dn = ggZZ_kf[3]->Eval(ma);

	double nominal= ggZZ_kf[0]->Eval(ma);

	double variation_pdf_up = sqrt( pow( max(pdfrep_up,pdfrep_dn)-nominal, 2 ) + pow( max(as_up,as_dn)-nominal,2) );
	double variation_pdf_dn = sqrt( pow( min(pdfrep_up,pdfrep_dn)-nominal, 2 ) + pow( min(as_up,as_dn)-nominal,2) );
	double variation_qcd_up = sqrt( pow( max(pdfscale_up,pdfscale_dn)-nominal, 2 ) + pow( max(qcdscale_up,qcdscale_dn)-nominal,2) );
	double variation_qcd_dn = sqrt( pow( min(pdfscale_up,pdfscale_dn)-nominal, 2 ) + pow( min(qcdscale_up,qcdscale_dn)-nominal,2) );
	xi[i]=ma;
	yi[i]=nominal;
	yi_pdf_up[i]=nominal+variation_pdf_up;
	yi_pdf_dn[i]=nominal-variation_pdf_dn;
	yi_qcd_up[i]=nominal+variation_qcd_up;
	yi_qcd_dn[i]=nominal-variation_qcd_dn;
	}
TGraph *sp_kfactor_qcd_dn = new TGraph (5801,xi,yi_qcd_dn);
TGraph *sp_kfactor_qcd_up = new TGraph (5801,xi,yi_qcd_up);
TGraph *sp_kfactor_pdf_dn = new TGraph (5801,xi,yi_pdf_dn);
TGraph *sp_kfactor_pdf_up = new TGraph (5801,xi,yi_pdf_up);
TGraph *sp_kfactor_nominal = new TGraph (5801,xi,yi);
TFile *fkf = new TFile("kfactor.root","recreate");
fkf->cd();
sp_kfactor_nominal->Write("sp_kfactor_Nominal");
sp_kfactor_qcd_up->Write("sp_kfactor_qcd_up");
sp_kfactor_qcd_dn->Write("sp_kfactor_qcd_dn");
sp_kfactor_pdf_up->Write("sp_kfactor_pdf_up");
sp_kfactor_pdf_dn->Write("sp_kfactor_pdf_dn");
fkf->Close();
}
