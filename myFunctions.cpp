#include "TMath.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TMatrixD.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

double FindTwoDimHistoMaxValue(TH2D* h){

	int NXBins = h->GetXaxis()->GetNbins();
	int NYBins = h->GetYaxis()->GetNbins();	
	double StartHistoMax = -9999.;
	double HistoMax = StartHistoMax;	

	for (int xbin = 1; xbin <= NXBins; xbin++) {

		for (int ybin = 1; ybin <= NYBins; ybin++) {		

			double LocalMax = h->GetBinContent(xbin,ybin);
			if (LocalMax > HistoMax && !isinf(LocalMax) ) { HistoMax = LocalMax; }

		}

	}

	if (HistoMax == StartHistoMax) { cout << "HistoMax = " << HistoMax << endl; }
	return HistoMax;

}

//----------------------------------------//

double FindTwoDimHistoMinValue(TH2D* h){

	int NXBins = h->GetXaxis()->GetNbins();
	int NYBins = h->GetYaxis()->GetNbins();	
	double StartHistoMin = 999999.;
	double HistoMin = StartHistoMin;	

	for (int xbin = 1; xbin<= NXBins; xbin++) {

		for (int ybin = 1; ybin<= NYBins; ybin++) {		

			double LocalMin = h->GetBinContent(xbin,ybin);
			if (LocalMin < HistoMin) { HistoMin = LocalMin; }

		}

	}

	if (HistoMin == StartHistoMin) { cout << "HistoMin = " << HistoMin << endl; }
	return HistoMin;

}

//----------------------------------------//

double FindOneDimHistoMaxValue(TH1D* h){

	int NBins = h->GetXaxis()->GetNbins();
	double HistoMax = -9999.;	

	for (int ibin = 1; ibin<= NBins; ibin++) {

		double LocalMax = h->GetBinContent(ibin);
		if (LocalMax > HistoMax) { HistoMax = LocalMax; }

	}

	return HistoMax;

}

//----------------------------------------//

double FindOneDimHistoMinValue(TH1D* h){

	int NBins = h->GetXaxis()->GetNbins();
	double HistoMin = 999999.;	

	for (int ibin = 1; ibin<= NBins; ibin++) {

		double LocalMin = h->GetBinContent(ibin);
		if (LocalMin < HistoMin) { HistoMin = LocalMin; }

	}

	return HistoMin;

}

//----------------------------------------//

TString to_string_with_precision(double a_value, const int n = 3)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return TString(out.str());
}

//----------------------------------------//

void CalcChiSquared(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval) {

	// Clone them so we can scale them 

	TH1D* h_model_clone = (TH1D*)h_model->Clone();
	TH1D* h_data_clone  = (TH1D*)h_data->Clone();
	TH2D* h_cov_clone   = (TH2D*)cov->Clone();
	int NBins = h_cov_clone->GetNbinsX();

	// Getting covariance matrix in TMatrix form

	TMatrixD cov_m;
	cov_m.Clear();
	cov_m.ResizeTo(NBins,NBins);

	// loop over rows

	for (int i = 0; i < NBins; i++) {			

		// loop over columns

		for (int j = 0; j < NBins; j++) {

			cov_m[i][j] = h_cov_clone->GetBinContent(i+1, j+1);
 
		}
	
	}

	TMatrixD copy_cov_m = cov_m;

	// Inverting the covariance matrix
	TMatrixD inverse_cov_m = cov_m.Invert();

	// Calculating the chi2 = Summation_ij{ (x_i - mu_j)*E_ij^(-1)*(x_j - mu_j)  }
	// x = data, mu = model, E^(-1) = inverted covariance matrix 

	chi = 0.;
	
	for (int i = 0; i < NBins; i++) {

		//double XWidth = h_data_clone->GetBinWidth(i+1);

		for (int j = 0; j < NBins; j++) {

			//double YWidth = h_data_clone->GetBinWidth(i+1);

			double diffi = h_data_clone->GetBinContent(i+1) - h_model_clone->GetBinContent(i+1);
			double diffj = h_data_clone->GetBinContent(j+1) - h_model_clone->GetBinContent(j+1);
			double LocalChi = diffi * inverse_cov_m[i][j] * diffj; 
			chi += LocalChi;
		}

	}

	ndof = h_data_clone->GetNbinsX();
	pval = TMath::Prob(chi, ndof);

	delete h_model_clone;
	delete h_data_clone;
	delete h_cov_clone;

}

// -------------------------------------------------------------------------------------------------------------------------------------

double PeLEE_ReturnBeamOnRunPOT(TString Run) {

	double DataPOT = -99.;

	if (Run == "Run1") { DataPOT = Fulltor860_wcut_Run1 ; }
	if (Run == "Run2") { DataPOT = Fulltor860_wcut_Run2 ; }
	if (Run == "Run3") { DataPOT = Fulltor860_wcut_Run3 ; }
	if (Run == "Run4") { DataPOT = Fulltor860_wcut_Run4 ; }
	if (Run == "Run4a") { DataPOT = Fulltor860_wcut_Run4a ; }
	if (Run == "Run4aRutgers") { DataPOT = Fulltor860_wcut_Run4aRutgers ; }		
	if (Run == "Run5") { DataPOT = Fulltor860_wcut_Run5 ; }
	if (Run == "Combined") { DataPOT = Fulltor860_wcut_Combined ; }

	return DataPOT;

}

// -------------------------------------------------------------------------------------------------------------------------------------

double ReturnBeamOnRunPOT(TString Run) {

	double DataPOT = -99.;

	if (Run == "Run1") { DataPOT = tor860_wcut_Run1 ; }
	if (Run == "Run2") { DataPOT = tor860_wcut_Run2 ; }
	if (Run == "Run3") { DataPOT = tor860_wcut_Run3 ; }
	if (Run == "Run4") { DataPOT = tor860_wcut_Run4 ; }
	if (Run == "Run5") { DataPOT = tor860_wcut_Run5 ; }

	return DataPOT;

}

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

double round(double var,int acc = 0) 
{ 
    double value = (int)(var * TMath::Power(10.,acc) + .5); 
    return (double)value / TMath::Power(10.,acc); 
} 

// -------------------------------------------------------------------------------------------------------------------------------------

double IntegratedXSec(TH1D* h) {

	int NBinsX = h->GetXaxis()->GetNbins();

	double IntegratedXSec = 0;

	for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

		double BinWidth = h->GetBinWidth(WhichXBin+1);
		double BinEntry = h->GetBinContent(WhichXBin+1);

		IntegratedXSec += BinEntry * BinWidth;

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

		IntegratedXSecErrorSquared += TMath::Power(BinError* BinWidth,2.);

	}

	double IntegratedXSecError = TMath::Sqrt(IntegratedXSecErrorSquared);

	return IntegratedXSecError;

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

TString ToStringFloat(float num) {

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

TH1D* SmearTrueToReco(TH1D* True, TFile* MigrationMatrixFile, TString PlotName, TString identifier) {

	TH2D* MigrationMatrix = (TH2D*)MigrationMatrixFile->Get("CC1pReco"+PlotName+"2D");

	TH1D* Reco = (TH1D*)(True->Clone("Reco"+PlotName+identifier));

	int XBins = MigrationMatrix->GetXaxis()->GetNbins();
	int YBins = MigrationMatrix->GetYaxis()->GetNbins();

	if (XBins != YBins) { std::cout << "Not symmetric matrix" << std::endl; }

	for (int WhichXBin = 0; WhichXBin < XBins; WhichXBin++) {

		double Entry = 0.;

		for (int WhichYBin = 0; WhichYBin < YBins; WhichYBin++) {

			double TrueInBin = True->GetBinContent(WhichYBin + 1);
			double MigrationInBin = MigrationMatrix->GetBinContent(WhichYBin + 1,WhichXBin + 1);

			Entry +=  MigrationInBin * TrueInBin;
	
		}

		// Bin entry in reco space

		Reco->SetBinContent(WhichXBin+1,Entry);

	}

	return 	Reco;

}

// --------------------------------------------------------------------------------------------------------------------------------------------

TH1D* ForwardFold(TH1D* True, TH1D* Reco, TH2D* MigrationMatrix) {

	// References (Marco & Lu)
	//https://microboone-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=30139&filename=NCE_MCC9_Internal_Note__v1_1%20_temp.pdf&version=1
	//https://microboone-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=14937&filename=numu_cc_inclusive_internal_note_v3.1.pdf&version=8

//	TH1D* ForwardFoldEfficiency = (TH1D*)(Reco->Clone("ForwardFoldEfficiency"));
	TH1D* ForwardFoldEfficiency = (TH1D*)(Reco->Clone());

	int XBins = MigrationMatrix->GetXaxis()->GetNbins();
	int YBins = MigrationMatrix->GetYaxis()->GetNbins();

	if (XBins != YBins) { std::cout << "Not symmetric matrix" << std::endl; }

	for (int WhichXBin = 0; WhichXBin < XBins; WhichXBin++) {

		double Num = 0.;
		double Den = 0.;

		double NumErrSquared = 0.;
		double DenErrSquared = 0.;

		for (int WhichYBin = 0; WhichYBin < YBins; WhichYBin++) {

			double RecoInBin = Reco->GetBinContent(WhichYBin + 1); 
			double TrueInBin = True->GetBinContent(WhichYBin + 1);
			double MigrationInBin = MigrationMatrix->GetBinContent(WhichYBin + 1,WhichXBin + 1);

			double ErrorRecoInBin = Reco->GetBinError(WhichYBin + 1); 
			double ErrorTrueInBin = True->GetBinError(WhichYBin + 1);
			double ErrorMigrationInBin = MigrationMatrix->GetBinError(WhichYBin + 1,WhichXBin + 1);

			Num +=  MigrationInBin * RecoInBin;
			Den +=  MigrationInBin * TrueInBin;

			NumErrSquared += TMath::Power(ErrorMigrationInBin * RecoInBin,2.) + TMath::Power(MigrationInBin * ErrorRecoInBin,2.);
			DenErrSquared += TMath::Power(ErrorMigrationInBin * TrueInBin,2.) + TMath::Power(MigrationInBin * ErrorTrueInBin,2.);

	
		}

		double NumErr = TMath::Sqrt(NumErrSquared);
		double DenErr = TMath::Sqrt(NumErrSquared);

		// FFEfficiency

		double FFefficiency = 0.;
		if (Den > 0) { FFefficiency = Num / Den; }
		ForwardFoldEfficiency->SetBinContent(WhichXBin+1,FFefficiency);

		// FFEfficiency Error
		double ErrorFFefficiency = 0.;
		if (Den > 0) { ErrorFFefficiency = FFefficiency * TMath::Sqrt( TMath::Power(NumErr/Num,2.) + TMath::Power(DenErr/Den,2.) ); }
		ForwardFoldEfficiency->SetBinError(WhichXBin+1,ErrorFFefficiency);

	}

	return 	ForwardFoldEfficiency;

}

// --------------------------------------------------------------------------------------------------------------------------------------------


