#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

using namespace std;

// -------------------------------------------------------------------------------------------------------------------------------------

float round(float var,int acc = 0) 
{ 
    float value = (int)(var * TMath::Power(10.,acc) + .5); 
    return (float)value / TMath::Power(10.,acc); 
} 

// -------------------------------------------------------------------------------------------------------------------------------------

double IntegratedXSec(TH1D* h) {

	int NBinsX = h->GetXaxis()->GetNbins();

	double IntegratedXSec = 0;

	for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

		double BinWidth = h->GetBinWidth(WhichXBin+1);
		double GenieBinEntry = h->GetBinContent(WhichXBin+1);

		IntegratedXSec += GenieBinEntry * BinWidth;

	}

	return IntegratedXSec;

}

// -------------------------------------------------------------------------------------------------------------------------------------

double IntegratedXSecError(TH1D* h) {
/*
	int NBinsX = h->GetXaxis()->GetNbins();

	double IntegratedXSecErrorSquared = 0;

	for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

		double BinWidth = h->GetBinWidth(WhichXBin+1);
		double GenieBinError = h->GetBinError(WhichXBin+1);

		IntegratedXSecErrorSquared += TMath::Power(GenieBinError,2.) * BinWidth;

	}

	double IntegratedXSecError = TMath::Sqrt(IntegratedXSecErrorSquared);

	return IntegratedXSecError;
*/

	int NBinsX = h->GetXaxis()->GetNbins();

	double IntegratedXSecError = 0;

	for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

		double BinWidth = h->GetBinWidth(WhichXBin+1);
		double GenieBinError = h->GetBinError(WhichXBin+1);

		IntegratedXSecError += GenieBinError * BinWidth;

	}

	return IntegratedXSecError;

}

// -------------------------------------------------------------------------------------------------------------------------------------

double Chi2(TH1D* h1,TH1D* h2) {

	int NBinsX = h1->GetXaxis()->GetNbins();

	double chi2 = 0;

	for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

		double h1Entry = h1->GetBinContent(WhichXBin+1);
		double h1Error = h1->GetBinError(WhichXBin+1);
		double h2Entry = h2->GetBinContent(WhichXBin+1);
		double h2Error = h2->GetBinError(WhichXBin+1);

		double num = TMath::Power(h1Entry - h2Entry,2.);
		double den = TMath::Power(h1Error,2.) + TMath::Power(h2Error,2.);
		chi2 += (num / den); 

	}

	return chi2;

}

// -------------------------------------------------------------------------------------------------------------------------------------

TString ToStringInt(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

TString ToString(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

