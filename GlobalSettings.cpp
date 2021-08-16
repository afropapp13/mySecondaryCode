#include <TStyle.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>

#include "ubana/myClasses/Constants.h"

// Aug 16 2021: Testing new authentication method

void GlobalSettings() {

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	
	gStyle->SetPalette(55); 
	gStyle->SetNumberContours(NCont); 
	gStyle->SetTitleSize(TextSize,"t"); 
	gStyle->SetTitleFont(FontStyle,"t");

	gStyle->SetOptStat(0);

	gROOT->ForceStyle();

};
