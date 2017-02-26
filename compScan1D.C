
#include <map>
const int nmax=2;
int co[14]={1,2,2,8,kViolet,kCyan, kMagenta,kBlue-10, 3,kYellow,93};
//bool sigma=true;
//bool onshell=false;
TString tex[nmax]={
	"Observed",
	"Expected" 
};

void dodo( bool sigma=false, bool onshell=false){
	//	bool pro = true;
	gROOT->ProcessLine(".x tdrstyle.cc");
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadTopMargin(0.05);

	TString fna = "offshell";
	if(onshell)
		fna = "onshell";
	TString rootfname[nmax]={
		"higgsCombineobsscan_2d_"+fna+"_merge.MultiDimFit.mH120.root",
		"higgsCombinescan_2d_"+fna+"_merge.MultiDimFit.mH120.root"
	};
	gROOT->ProcessLine(".x tdrstyle.cc");
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadTopMargin(0.05);
	TCanvas *c1=new TCanvas("can1","CANVAS-SCAN1D",800,800);
	TCanvas *c2=new TCanvas("can2","CANVAS-SCAN1D",800,800);
	TLegend *leg = new TLegend(0.24,0.70,0.4,0.83);
	leg->SetFillColor(0);
	leg->SetLineColor(0);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.04);


	int start = 0;

	for (int j=0;j<2;j++){
		TChain *t=new TChain("limit");
		t->Add( rootfname[j]);

		c2->cd();
		float offset=0;
		float offset_125=0;
		float CMS_zz4l_fg4,CMS_zz4l_fg2,deltaNLL;
		t->SetBranchAddress("mean_pole",&CMS_zz4l_fg4);
		t->SetBranchAddress("sigma_pole",&CMS_zz4l_fg2);
		t->SetBranchAddress("deltaNLL",&deltaNLL);

		float lowest = 55;
		float lowest_125= 55;


		std::map <float, float> a4;
		for (int entry=0;entry<t->GetEntries();entry++){
			t->GetEntry(entry);
			if(deltaNLL==0)
				continue;
			float content=2*deltaNLL;

			if(!sigma){
				if(a4.count(CMS_zz4l_fg4)){
					if(content<a4[CMS_zz4l_fg4])
						a4[CMS_zz4l_fg4]= content;
				}
				else{
					a4[CMS_zz4l_fg4]= content;
				}
			}
			else{
				if(a4.count(CMS_zz4l_fg2)){
					if(content<a4[CMS_zz4l_fg2])
						a4[CMS_zz4l_fg2]= content;
				}
				else{
					a4[CMS_zz4l_fg2]= content;
				}
			}

			if (deltaNLL<lowest){
				lowest = deltaNLL;
				offset = deltaNLL;
			}
			if(CMS_zz4l_fg4>125.08&&CMS_zz4l_fg4<125.11){
				if(deltaNLL<lowest_125){
					offset_125= deltaNLL;
					lowest_125= deltaNLL;
				}
			}

		}


		cout<<offset<<endl;
		cout<<offset_125<<endl;

		int point=0;
		TGraph *gr0 ;
		//	if(pro)
		gr0=new TGraph(a4.size());
		for(std::map<float,float>::iterator it = a4.begin(); it != a4.end(); it++) {
//			if(it->first>4.9&&sigma)
//				continue;
			if(!onshell&&j==1&&sigma)
				gr0->SetPoint(point,it->first,(it->second-offset*2)); 
			else
				gr0->SetPoint(point,it->first,it->second-offset*2); 

			//		cout << it->first<<" "<<it->second<<endl;
			point++;
		}

		//		offset = 0.3585195;
		//		if(onshell)
		//			offset=0;
		TString offstring = Form("%f",offset);
		TString offstring_125= Form("%f",offset_125);


		//	if(!pro){
		if(sigma) 
			if(j!=0 &&!onshell)
				t->Draw("2*(deltaNLL-"+offstring_125+"):sigma_pole", "deltaNLL&&(mean_pole<125.11&&mean_pole>125.08)","PL");
			else
				t->Draw("2*(deltaNLL-"+offstring_125+"):sigma_pole", "deltaNLL&&(mean_pole<125.11&&mean_pole>125.08)","PL");
		else{
			if(!onshell)
				t->Draw("2*(deltaNLL-"+offstring+"):mean_pole", "deltaNLL&&(sigma_pole>0.01&&sigma_pole<0.011)","PL");
			else
				t->Draw("2*(deltaNLL-"+offstring+"):mean_pole", "deltaNLL&&(sigma_pole>0.&&sigma_pole<0.1)","PL");
		}
		TGraph *gr1 = (TGraph*) gROOT->FindObject("Graph")->Clone();
		//	}

		gr0->Sort();
		gr1->Sort();
		{
			double sigmav[5] = {0.,1.,3.84,1.,3.84};
			double x1p[5]={-0.0001};
			double sigma1=500;
			double x=0;
			while(x<5){
				double cu=gr0->Eval(x);
				for(int k=0;k<3;k++){
					//    if (fabs(cu-sigmav[k])<0.03){
					if (fabs(cu-sigmav[k])<0.0001){
						cout<<  "Expected: "<<std::fixed<<sigmav[k]<< " "<<x<< " "<< cu<<endl;
					}
				}
				x+=0.0001;
				}
			}
			cout.precision(4);
			double y;
			double ypre;
			double ypre_pre;
			double yaft;
			double yaft_aft;
			double x;
			for(int p= 2;p<gr0->GetN()-1;p++){
				gr0->GetPoint(p+1,x,yaft);	
				gr0->GetPoint(p-1,x,ypre);	
				gr0->GetPoint(p,x,y);	
				if(y!=0&&y<ypre&&y<yaft){
					gr0->RemovePoint(p);
					p--;
				}
			}
//			gr0->Draw("psame");


			gr0->SetName("Exp1D");
			gr0->SetLineWidth(2);
			gr1->SetLineWidth(2);
			gr1->SetLineColor(1);
			gr0->SetTitle("");
			gr0->GetXaxis()->SetLabelSize(0.04);
			gr0->GetYaxis()->SetLabelSize(0.04);
			gr0->GetXaxis()->SetTitleSize(0.06);
			gr0->GetYaxis()->SetTitleSize(0.06);
			gr0->GetXaxis()->SetLabelFont(42);
			gr0->GetYaxis()->SetLabelFont(42);
			gr0->GetXaxis()->SetTitleFont(42);
			gr0->GetYaxis()->SetTitleFont(42);
			c1->cd();
			//	gr0->SetLineColor(co[j]);
			if(j!=0){
				gr0->SetLineStyle(2);
				gr1->SetLineStyle(2);
			}
			TString			    title = "";
			title="m_{H} [GeV]";
			gr0->GetYaxis()->SetRangeUser(0.,16);
			gr1->GetYaxis()->SetRangeUser(0.,16);
			//gr0->GetXaxis()->SetRangeUser(122.5,127.5);
			gr0->GetXaxis()->SetRangeUser(123,127);
			if(sigma){
				title="#Gamma_{H} [GeV]";
			}

			gr0->GetXaxis()->SetTitle(title);
			gr0->GetYaxis()->SetTitle("-2 #Delta lnL");
			gr0->GetXaxis()->SetLabelSize(0.04);
			gr0->GetYaxis()->SetLabelSize(0.04);
			gr1->GetXaxis()->SetTitle(title);
			gr1->GetYaxis()->SetTitle("-2 #Delta lnL");
			gr1->GetXaxis()->SetLabelSize(0.04);
			gr1->GetYaxis()->SetLabelSize(0.04);
			double max = gr0->GetYaxis()->GetXmax();

			//gr0->GetXaxis()->SetRangeUser(-0.98,0.98);
			//  gr0->GetXaxis()->CenterTitle();
			//  gr0->GetYaxis()->CenterTitle();
			if(j==0){
				//gr0->Draw("AL");
				//gr0->SetLineColor(1);
				gr1->Draw("AL");
			}
			if(j!=0){
				//	gr0->Draw("same");
				gr1->Draw("same");
			}
			if(sigma){
				if(!onshell){
					gr0->GetXaxis()->SetRangeUser(0.00,0.096);
					gr1->GetXaxis()->SetRangeUser(0.00,0.096);
				}
				else{
					gr0->GetXaxis()->SetRangeUser(0.00,4.8);
					gr1->GetXaxis()->SetRangeUser(0.00,4.8);
				}
			}
			else
				gr0->GetXaxis()->SetRangeUser(123,127);
			gPad->Update();

			cout.precision(4);
			leg->AddEntry(gr1, tex[j],"l");
			//leg->AddEntry(gr1,  tex[j]+", m_{H}=125.09 GeV","l");

		}
		leg->Draw();

		TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
		pt->SetBorderSize(0);
		pt->SetTextAlign(12);
		pt->SetFillStyle(0);
		pt->SetTextFont(42);
		pt->SetTextSize(0.04);
		TText *text = pt->AddText(0.02,0.45,"#font[61]{CMS}");
		text->SetTextSize(0.044);
			text = pt->AddText(0.14, 0.45, "#it{Preliminary}");
			text->SetTextSize(0.0315);
			text = pt->AddText(0.72,0.45,"#font[42]{12.9 fb^{-1} (13 TeV)}");
			text->SetTextSize(0.0315);
		pt->Draw();

		TPaveText *pts = new TPaveText(0.1577181,0.8562937,0.9580537,0.8947552,"brNDC");
		pts->SetBorderSize(0);
		pts->SetTextAlign(12);
		pts->SetFillStyle(0);
		pts->SetTextFont(42);
		pts->SetTextSize(0.04);
		if(onshell){
			TText *texts = pts->AddText(0.05,0.45,"105 GeV < m_{4l}< 140 GeV");
			texts->SetTextSize(0.044);
		}
		else{
			TText *texts = pts->AddText(0.05,0.45,"105 GeV < m_{4l}< 1600 GeV");
			texts->SetTextSize(0.044);
		}
		pts->Draw();
		//text = pt->AddText(0.2,0.6,Form("#sqrt{s} = 7 TeV, L = %.1f fb^{-1}  #sqrt{s} = 8 TeV, L = %.1f fb^{-1}",lumi7TeV,lumi8TeV));
		//  text = pt->AddText(0.2,0.6,Form("#sqrt{s} = 8 TeV, L = %.1f fb^{-1}, 2e2#mu",lumi8TeV));
		//pt->Draw();  

		TPaveText *oneSig = new TPaveText(0.175,0.18,0.275,0.21,"NDC");
		if(j==0)
			oneSig = new TPaveText(0.175,0.25,0.275,0.28,"NDC");

		oneSig->SetFillColor(0);
		oneSig->SetTextFont(42);
		oneSig->SetTextColor(1);
		oneSig->SetBorderSize(0);
		//oneSig->AddText("1#sigma"); 
		oneSig->AddText("68\% CL"); 
		//	if(version!=3)
		oneSig->Draw();

		TPaveText *twoSig = new TPaveText(0.175,0.35,0.275,0.38,"NDC");
		if(j==0)
			twoSig=new TPaveText(0.175,0.44,0.275,0.48,"NDC");
		twoSig->SetFillColor(0);
		twoSig->SetTextFont(42);
		twoSig->SetTextColor(1);
		twoSig->SetBorderSize(0);
		twoSig->AddText("95\% CL");
		twoSig->Draw();

		TLine *l1=new TLine();
		l1->SetLineStyle(9);
		l1->SetLineWidth(2);
		//l1->SetLineColor(kRed);
		if(sigma){
			if(!onshell)
				l1->DrawLine(0,1,0.096,1.0);
			else
				l1->DrawLine(0,1,4.8,1.0);
		}
		else
			l1->DrawLine(122.5,1,127.5,1.0);

		l1->Draw("same");
		TLine *l2=new TLine();
		l2->SetLineStyle(9);
		l2->SetLineWidth(2);
		if(sigma){
			if(!onshell)
				l2->DrawLine(0,3.84,0.096,3.84);
			else
				l2->DrawLine(0,3.84,4.8,3.84);
		}
		else
			l2->DrawLine(122.5,3.84,127.5,3.84);
		l2->Draw("same");

		TString app = "obs_"+fna+"_mean";
		if(sigma)
			app = "obs_"+fna+"_sigma_125.09";

		c1->SaveAs(app+".png");
		c1->SaveAs(app+".pdf");
		c1->SaveAs(app+".eps");
		c1->SaveAs(app+".root");
		return;


	}
	void compScan1D(){
		//dodo(0,1);
		//dodo(0,0);
		//
		dodo(1,1);
	//	dodo(1,0);
	}
