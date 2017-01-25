
#include "external_cConstants.h"
//float getDVBF2jetsConstant(float ZZMass){
//  float par[9]={
//    1.876,
//    -55.488,
//    403.32,
//    0.3906,
//    80.8,
//    27.7,
//    -0.06,
//    54.97,
//    309.96
//  };
//  float kappa =
//    pow(1.-atan((ZZMass-par[1])/par[2])*2./TMath::Pi(), par[0])
//    + par[3]*exp(-pow((ZZMass-par[4])/par[5], 2))
//    + par[6]*exp(-pow((ZZMass-par[7])/par[8], 2));
//  float constant = kappa/(1.-kappa);
//  return constant;
//}

void data(){

		
		TString treename [3]={"ZZTree/candTree","ZZTreelooseEle/candTree","ZZTreetle/candTree"};
		TString newtreename[3]={"","_rse","_tle"};
		TFile* fnew = new TFile("data.root","recreate");
		TTree *tnew = new TTree("SelectedTree","SelectedTree");
		for(int t =0;t<1;t++){
	  TChain *tqqzz= new TChain(treename[t]);
//	  tqqzz->Add("root://lxcms03//data3/Higgs/160725/AllData/ZZ4lAnalysis.root");
	  tqqzz->Add("root://lxcms03//data3/Higgs/160729_complete/AllData/ZZ4lAnalysis.root");
		float ZZPt,ZZMass;
		vector<float> *LepPt=new vector<float>;
		short Z1Flav,Z2Flav;
		short nCleanedJetsPt30;
		float pvbf_VAJHU_old;
		float phjj_VAJHU_old;
		float bkg_VAMCFM,p0plus_VAJHU;
		short ZZsel;
		float TLE_dR_Z;
		vector<short> *LepLepId=0;
		tqqzz->SetBranchAddress("pvbf_VAJHU_highestPTJets",&pvbf_VAJHU_old);
		tqqzz->SetBranchAddress("phjj_VAJHU_highestPTJets",&phjj_VAJHU_old);
		tqqzz->SetBranchAddress("p0plus_VAJHU",&p0plus_VAJHU);
		tqqzz->SetBranchAddress("bkg_VAMCFM",&bkg_VAMCFM);
		tqqzz->SetBranchAddress("ZZPt",&ZZPt);
		tqqzz->SetBranchAddress("ZZMass",&ZZMass);
		tqqzz->SetBranchAddress("Z1Flav",&Z1Flav);
		tqqzz->SetBranchAddress("Z2Flav",&Z2Flav);
		tqqzz->SetBranchAddress("ZZsel",&ZZsel);
		tqqzz->SetBranchAddress("TLE_dR_Z",&TLE_dR_Z);
		tqqzz->SetBranchAddress("LepLepId",&LepLepId);
		tqqzz->SetBranchAddress("LepPt",&LepPt);
		tqqzz->SetBranchAddress("nCleanedJetsPt30",&nCleanedJetsPt30);
		int chan;
		int vbfcate;
		float dbkg_kin;
		tnew->Branch("mreco",&ZZMass,"mreco/F");
		tnew->Branch("dbkg_kin",&dbkg_kin,"dbkg_kin/F");
		tnew->Branch("chan",&chan,"chan/I");
		tnew->Branch("vbfcate",&vbfcate,"vbfcate/I");
		for(int i=0;i<tqqzz->GetEntries();i++){
			tqqzz->GetEntry(i);
			if(ZZsel!=120 && t>0)
				continue;
			if(t>0)
				if(ZZMass<300)
					continue;

			if(t!=2){
			if(abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)==121){
				chan=2;
			}
			else if (abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)!=121){
				chan=1;
			}
			else{
				chan=3;
			}
			}
			bool patle = true;
			if(t==2){
							for (int k=0 ; k<LepLepId->size() ; k++){ 
										if (abs(LepLepId->at(k))==22 && LepPt->at(k)<=30) 
											patle = false;
							}
							 if(TLE_dR_Z<=1.6)
								 patle = false;
							if(abs(Z1Flav*Z2Flav)==29282)
								chan =2;
							else if(abs(Z1Flav*Z2Flav)==40898)
								chan =3;
			}
			if(!patle)
				continue;	
			double c =getDVBF2jetsConstant(ZZMass);
			float vbfMela = 1/(1+ c*phjj_VAJHU_old/pvbf_VAJHU_old);
//			float vbfMela = pvbf_VAJHU_old / ( phjj_VAJHU_old + pvbf_VAJHU_old );
//			if(vbfMela>0.5  && nCleanedJetsPt30>=2)
			if(t==1)
				 vbfcate =2;
			else{
				 if(vbfMela> (1.043-460./(ZZMass+634.)) && nCleanedJetsPt30>=2)
					vbfcate=1;
					else
						vbfcate=0;
			}

			short ZZFlav = Z1Flav*Z2Flav;
			dbkg_kin = p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM*getDbkgkinConstant(ZZFlav, ZZMass));
			tnew->Fill();
		}
		tqqzz->SetLineColor(2);
		tqqzz->SetMarkerColor(2);
		}
		fnew->cd();
		tnew->Write();
		fnew->Close();
}
