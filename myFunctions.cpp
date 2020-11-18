#include "TMath.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TMatrix.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

using namespace std;

// -------------------------------------------------------------------------------------------------------------------------------------

void Reweight2D(TH2D* h, double SF) {

	int NBinsX = h->GetXaxis()->GetNbins();
	int NBinsY = h->GetYaxis()->GetNbins();

	for (int i = 0; i < NBinsX; i++) {

		for (int j = 0; j < NBinsX; j++) {

			double CurrentEntry = h->GetBinContent(i+1,j+1);
			double NewEntry = CurrentEntry * SF / ( h->GetXaxis()->GetBinWidth(i+1) * h->GetYaxis()->GetBinWidth(j+1) );

			double CurrentError = h->GetBinError(i+1,j+1);
			double NewError = CurrentError * SF / ( h->GetXaxis()->GetBinWidth(i+1) * h->GetYaxis()->GetBinWidth(j+1) );

			h->SetBinContent(i+1,j+1,NewEntry); 
//			h->SetBinError(i+1,j+1,NewError); 
			h->SetBinError(i+1,j+1,0.000001); 

		}

	}

}

// -------------------------------------------------------------------------------------------------------------------------------------

void Reweight(TH1D* h, double SF) {

	int NBins = h->GetXaxis()->GetNbins();

	for (int i = 0; i < NBins; i++) {

		double CurrentEntry = h->GetBinContent(i+1);
		double NewEntry = CurrentEntry * SF / h->GetBinWidth(i+1);

		double CurrentError = h->GetBinError(i+1);
		double NewError = CurrentError * SF / h->GetBinWidth(i+1);

		h->SetBinContent(i+1,NewEntry); 
//		h->SetBinError(i+1,NewError); 
		h->SetBinError(i+1,0.000001); 

	}

}

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

	int NBinsX = h->GetXaxis()->GetNbins();

	double IntegratedXSecErrorSquared = 0;

	for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

		double BinWidth = h->GetBinWidth(WhichXBin+1);
		double BinError = h->GetBinError(WhichXBin+1);

		IntegratedXSecErrorSquared += TMath::Power(BinError,2.) * BinWidth;

	}

	double IntegratedXSecError = TMath::Sqrt(IntegratedXSecErrorSquared);

	return IntegratedXSecError;

/*
	int NBinsX = h->GetXaxis()->GetNbins();

	double IntegratedXSecError = 0;

	for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

		double BinWidth = h->GetBinWidth(WhichXBin+1);
		double BinError = h->GetBinError(WhichXBin+1);

		IntegratedXSecError += BinError * BinWidth;

	}

	return IntegratedXSecError;
*/
}

// -------------------------------------------------------------------------------------------------------------------------------------

double Chi2(TH1D* h1,TH1D* h2, int LowBin = -1, int HighBin = -1) {

	int NBinsX = h1->GetXaxis()->GetNbins();

	double chi2 = 0;
	
	if (LowBin == -1) { LowBin = 0; }
	if (HighBin == -1) { HighBin = NBinsX; }	

	for (int WhichXBin = LowBin; WhichXBin < HighBin; WhichXBin++) {

		double h1Entry = h1->GetBinContent(WhichXBin+1);
		double h1Error = h1->GetBinError(WhichXBin+1);
		double h2Entry = h2->GetBinContent(WhichXBin+1);
		double h2Error = h2->GetBinError(WhichXBin+1);

		double num = TMath::Power(h1Entry - h2Entry,2.);
		double den = TMath::Power(h1Error,2.) + TMath::Power(h2Error,2.);
		if (den != 0) { chi2 += (num / den); }

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

// -------------------------------------------------------------------------------------------------------------------------------------

TString ToString(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

// -------------------------------------------------------------------------------------------------------------------------------------

double computeMean(std::vector<double> numbers) {

	if(numbers.empty()) return 0;

	double total = 0;
	for (int number = 0; number < (int)(numbers.size()); number ++) {
		total += numbers[number];
	}

	double average = total / numbers.size();
	return average;
}

// -------------------------------------------------------------------------------------------------------------------------------------

double computeStd(double mean, std::vector<double> numbers) {

	float result = 0;
	for (int number = 0; number < (int)(numbers.size()); number ++) {
		result += (numbers[number] - mean)*(numbers[number] - mean);
	}

	return sqrt(result / (numbers.size() - 1));
}

// --------------------------------------------------------------------------------------------------------------------------------------------

TH1D* ForwardFold(TH1D* True, TH1D* Reco, TH2D* MigrationMatrix) {

	TH1D* ForwardFoldEfficiency = (TH1D*)(Reco->Clone("ForwardFoldEfficiency"));

	int XBins = True->GetXaxis()->GetNbins();
	int YBins = True->GetYaxis()->GetNbins();

	if (XBins != YBins) { std::cout << "Not symmetric matrix" << std::endl; }


	for (int WhichXBin = 0; WhichXBin < XBins; WhichXBin++) {

		double Num = 0;
		double Den = 0;

		for (int WhichYBin = 0; WhichYBin < YBins; WhichYBin++) {

			double RecoInBin = Reco->GetBinContent(WhichYBin + 1); 
			double TrueInBin = True->GetBinContent(WhichYBin + 1);
			double MigrationInBin = MigrationMatrix->GetBinContent(YBins - WhichYBin,WhichXBin + 1);

			Num +=  MigrationInBin * RecoInBin;
			Den +=  MigrationInBin * TrueInBin;

	
		}

		double FFefficiency = 0.;
		if (Den > 0) { FFefficiency = Num / Den; }
		ForwardFoldEfficiency->SetBinContent(WhichXBin+1,FFefficiency);

	}

	return 	ForwardFoldEfficiency;

}

// --------------------------------------------------------------------------------------------------------------------------------------------


