#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

// This code is designed to calculate the Bethe-Heitler bremsstralung yield using recattering flux calculations.

// The code uses atomic units

//Define Constants
double Z = 54;
double c = 137;
double m = 1;
double e = 1;
double hbar = 1;

// Short and Long flux array length
const int fluxlength = 50;

double extern longflux[fluxlength][2]; // extern is so I can forward declare the array
double extern shortflux[fluxlength][2];

double brem_per_Egamma(double, double);
void sum_flux(vector<vector<double> >, int*);


// Begin main //
int main() {
	vector<vector<double> > flux;
	int count = 0;
	sum_flux(flux, &count);
	cout << count << endl;
	cout << flux[0][0] << endl;
	system("pause");
	return 0;
}
// End main //


// Define subroutines // 
double brem_per_Egamma(double Er, double Eg) {
	double E0 = Er + m * c * c;
	double Ef = E0 - Eg;
	double pf = sqrt(pow(Ef, 2.0) - pow(m * pow(c, 2.0), 2.0)); // "pf" and "p0" is defined as "momentum x c" in the Bethe-Heitler expression
	double p0 = sqrt(pow(E0, 2.0) - pow(m * pow(c, 2.0), 2.0));
	double mu = m * pow(c, 2.0);
	double eps0 = 2 * log((E0 + p0) / mu);
	double eps = 2 * log((Ef + pf) / mu);
	double L = 2 * log((E0 * Ef + p0 * pf - pow(mu, 2.0)) / (mu * Eg));
	double alpha = pow(e, 2.0) / (hbar * c);
	// Bethe-Heitler equation
	double BH = alpha * pow(Z, 2.0) * pow(pow(e, 2.0) / mu, 2.0) * (pf / p0) * (1 / Eg) * (
		(4.0 / 3.0) - (2 * Ef * E0 * ((pow(pf, 2.0) + pow(p0, 2.0)) / pow(p0 * pf, 2.0))) +
		pow(mu, 2.0) * ((eps0 * Ef / pow(p0, 3.0)) + (eps * E0 / pow(pf, 3.0)) - (eps * eps0 / (p0 * pf))) +
		(((8.0 / 3.0) * (Ef * E0 / (pf * p0))) + (pow(Eg, 2.0) / pow(p0 * pf, 3.0)) * (pow(E0 * Ef, 2.0) + pow(pf * p0, 2.0))) * L +
		L * (pow(mu, 2.0) * Eg / (2 * p0 * pf)) * (((E0 * Ef + pow(p0, 2.0)) / pow(p0, 3.0)) * eps0 - ((E0 * Ef + pow(pf, 2.0)) / pow(pf, 3.0)) * eps +
			(2 * Eg * E0 * Ef / pow(pf * p0, 2.0))));
	return BH;
}


void sum_flux(vector<vector<double>> flux, int *count) {
	for (int i=1;i<=fluxlength;i++) {
		if ((longflux[i - 1][1] + shortflux[i - 1][1]) != double(0.0)) {
			*count = *count + 1;
			flux.push_back({ longflux[i - 1][0], longflux[i - 1][1] + shortflux[i - 1][1] });
		}
	}
}

//Define Long Flux array
double longflux[fluxlength][2]{
{0.01033,0},
{0.01296,0},
{0.01626,0},
{0.02039,0},
{0.02557,0},
{0.03208,0},
{0.04024,0},
{0.05047,0},
{0.06331,0},
{0.07941,0},
{0.0996,1.28898E-10},
{0.12493,0},
{0.15671,0},
{0.19656,0},
{0.24655,0},
{0.30925,0},
{0.3879,0},
{0.48656,0},
{0.6103,0},
{0.76551,0},
{0.9602,0},
{1.2044,0},
{1.5107,0},
{1.89491,0},
{2.37683,0},
{2.98131,0},
{3.73952,0},
{4.69057,0},
{5.88349,0},
{7.37979,0},
{9.25664,0},
{11.6108,0},
{14.5637,0},
{18.2676,0},
{22.9134,0},
{28.7409,0},
{36.0503,0},
{45.2187,0},
{56.7188,0},
{71.1437,6.74325E-11},
{89.2372,0},
{111.932,0},
{140.399,5.99144E-11},
{176.106,0},
{220.893,6.19163E-11},
{277.072,6.73793E-11},
{347.537,8.58706E-11},
{435.924,1.34232E-10},
{546.789,3.29621E-10},
{685.85,7.12869E-9}
};

//Define Short Flux array
double shortflux[fluxlength][2]{
{0.01033,1.45482E-298},
{0.01296,0},
{0.01626,8.23573E-196},
{0.02039,0},
{0.02557,0},
{0.03208,0},
{0.04024,0},
{0.05047,8.9211E-144},
{0.06331,1.06472E-113},
{0.07941,0},
{0.0996,0},
{0.12493,0},
{0.15671,0},
{0.19656,3.62986E-79},
{0.24655,0},
{0.30925,0},
{0.3879,0},
{0.48656,8.7493E-60},
{0.6103,0},
{0.76551,7.23917E-53},
{0.9602,0},
{1.2044,1.69748E-47},
{1.5107,1.15667E-42},
{1.89491,2.18111E-39},
{2.37683,5.82378E-36},
{2.98131,1.49036E-33},
{3.73952,6.51133E-31},
{4.69057,3.64022E-29},
{5.88349,4.56013E-27},
{7.37979,8.96791E-26},
{9.25664,4.69186E-24},
{11.6108,4.37353E-23},
{14.5637,3.87545E-21},
{18.2676,7.33798E-20},
{22.9134,1.83485E-18},
{28.7409,1.62842E-17},
{36.0503,1.6233E-16},
{45.2187,1.98582E-15},
{56.7188,1.93467E-14},
{71.1437,1.19788E-13},
{89.2372,6.06612E-13},
{111.932,3.54333E-12},
{140.399,1.98764E-11},
{176.106,8.29703E-11},
{220.893,2.70027E-10},
{277.072,8.28815E-10},
{347.537,2.17526E-9},
{435.924,5.08796E-9},
{546.789,9.04089E-9},
{685.85,2.37257E-8}
};
