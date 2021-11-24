#include "TMath.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TMatrix.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

using namespace std;

//----------------------------------------//

double round(double var,int acc = 0) {

    double value = (int)(var * TMath::Power(10.,acc) + .5); 
    return (double)value / TMath::Power(10.,acc); 
} 

//----------------------------------------//

void CalcChiSquared(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval) {

	// Clone them so we can scale them 

	TH1D* h_model_clone = (TH1D*)h_model->Clone();
	TH1D* h_data_clone  = (TH1D*)h_data->Clone();
	TH2D* h_cov_clone   = (TH2D*)cov->Clone();
	int NBins = h_cov_clone->GetNbinsX();

	double ScaleFactor = 1.;

	// Getting covariance matrix in TMatrix form

	TMatrixD cov_m;
	cov_m.Clear();
	cov_m.ResizeTo(NBins,NBins);

	// loop over rows

	for (int i = 0; i < NBins; i++) {

		//double BinWidth = h_data_clone->GetBinWidth(i+1);

		//double MCEntry = h_model_clone->GetBinContent(i+1);
		//double DataEntry = h_data_clone->GetBinContent(i+1);

		//h_model_clone->SetBinContent(i+1,MCEntry*BinWidth);
		//h_data_clone->SetBinContent(i+1,DataEntry*BinWidth);				

		// loop over columns

		for (int j = 0; j < NBins; j++) {

			cov_m[i][j] = h_cov_clone->GetBinContent(i+1, j+1) * ScaleFactor; // Scale by ScaleFactor otherwise infinities

//			cov_m[i][j] = 1.;
//			if (i != j) { cov_m[i][j] = 0.; }
 
		}
	
	}

	TMatrixD copy_cov_m = cov_m;

	// Inverting the covariance matrix
	//cov_m.SetTol(1.e-23);
	TMatrixD inverse_cov_m = cov_m.Invert();
	//inverse_cov_m.SetTol(1.e-23);

	// Calculating the chi2 = Summation_ij{ (x_i - mu_j)*E_ij^(-1)*(x_j - mu_j)  }
	// x = data, mu = model, E^(-1) = inverted covariance matrix 

	chi = 0.;
	
	for (int i = 0; i < NBins; i++) {

		//double XWidth = h_data_clone->GetBinWidth(i+1);

		for (int j = 0; j < NBins; j++) {

			//double YWidth = h_data_clone->GetBinWidth(i+1);

			double diffi = h_data_clone->GetBinContent(i+1) - h_model_clone->GetBinContent(i+1);
			double diffj = h_data_clone->GetBinContent(j+1) - h_model_clone->GetBinContent(j+1);
			double LocalChi = diffi * inverse_cov_m[i][j] * diffj * ScaleFactor; 
			chi += LocalChi;
//if (i == j) { cout << "i = " << i << " j = " << j << " diffi = " << diffi << " diffj = " << diffj << "  inv = " << inverse_cov_m[i][j] << "  reg m = " << copy_cov_m[i][j] << " chi = " << chi << " LocalChi = " << LocalChi << endl; }
//cout << "i = " << i << " j = " << j << " diffi = " << diffi << " diffj = " << diffj << "  inv = " << inverse_cov_m[i][j] << "  reg m = " << copy_cov_m[i][j] << " chi = " << chi << " LocalChi = " << LocalChi << endl;

		}

	}

//cout << endl;

//TCanvas* Canvas = new TCanvas("canvas","canvas",205,34,1024,768);
//copy_cov_m.Draw("coltz text");

//TCanvas* InvCanvas = new TCanvas("invcanvas","invcanvas",205,34,1024,768);
//inverse_cov_m.Draw("coltz text");

//TCanvas* ProdCanvas = new TCanvas("Prodcanvas","Prodcanvas",205,34,1024,768);
//TMatrixD product = inverse_cov_m * copy_cov_m;
//product.Draw("coltz text");

	ndof = h_data_clone->GetNbinsX();
	pval = TMath::Prob(chi, ndof);

	//std::cout << "Chi2/dof: " << chi << "/" << h_data_clone->GetNbinsX() << " = " << chi/double(ndof) <<  std::endl;
	//std::cout << "p-value: " <<  pval << "\n" <<  std::endl;

	delete h_model_clone;
	delete h_data_clone;
	delete h_cov_clone;

}

// --------------------------------------------------------------------------------------------------------------------------------------------

TH1D* Multiply(TH1D* True, TH2D* SmearMatrix) {

	TH1D* TrueClone = (TH1D*)(True->Clone());

	int XBins = SmearMatrix->GetXaxis()->GetNbins();
	int YBins = SmearMatrix->GetYaxis()->GetNbins();

	if (XBins != YBins) { std::cout << "Not symmetric matrix" << std::endl; }

	TVectorD signal(XBins);
	TMatrixD response(XBins,XBins);

	H2V(TrueClone, signal);
	H2M(SmearMatrix, response, kTRUE);

	TVectorD RecoSpace = response * signal;
	V2H(RecoSpace, TrueClone);	

	return TrueClone;

}

// -------------------------------------------------------------------------------------------------------------------------------------

double Chi2Func(TH1D* h1,TH1D* h2, int LowBin = -1, int HighBin = -1) {

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
		double den = TMath::Sqrt( TMath::Power(h1Error,2.) + TMath::Power(h2Error,2.) );
		if (den != 0) { chi2 += (num / den); }

	}

	return chi2;

}

// -------------------------------------------------------------------------------------------------------------------------------------

TString to_string_with_precision(double a_value, const int n = 3)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return TString(out.str());
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

		IntegratedXSecErrorSquared += TMath::Power(BinError * BinWidth,2.);

	}

	double IntegratedXSecError = TMath::Sqrt(IntegratedXSecErrorSquared);

	return IntegratedXSecError;

}

// -------------------------------------------------------------------------------------------------------------------------------------

TString ToString(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

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
