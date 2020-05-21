#include <TH1D.h>

#include <iostream>

using namespace std;

void MakeMyPlotPretty(TH1D* histo){

		int FontStyle = 132;

		histo->SetLineWidth(3);
		histo->GetXaxis()->CenterTitle();
		histo->GetXaxis()->SetTitleFont(FontStyle);
		histo->GetXaxis()->SetLabelFont(FontStyle);
		histo->GetXaxis()->SetTitleSize(0.06);
		histo->GetXaxis()->SetLabelSize(0.04);
		histo->GetXaxis()->SetTitleOffset(0.7);
		histo->GetXaxis()->SetNdivisions(5);

		histo->GetYaxis()->CenterTitle();
		histo->GetYaxis()->SetTitleFont(FontStyle);
		histo->GetYaxis()->SetTitleSize(0.12);
		histo->GetYaxis()->SetLabelFont(FontStyle);
		histo->GetYaxis()->SetNdivisions(6);
		histo->GetYaxis()->SetTitleOffset(0.8);
		histo->GetYaxis()->SetTitleSize(0.06);
		histo->GetYaxis()->SetLabelSize(0.04);
};
