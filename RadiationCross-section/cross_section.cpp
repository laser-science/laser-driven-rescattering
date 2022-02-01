#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

// This code is designed to calculate the Bethe-Heitler bremsstralung yield using recattering flux calculations.

// The code uses atomic units

//Define Constants
const double IP = 5.27414; // IP of Ar+6 in hartree
const double Z = 18;
const double c = 137;
const double m = 1;
const double e = 1;
const double hbar = 1;
// Short and Long flux array length
const int fluxlength = 55;


// Forward declare flux arrays
double extern longflux[fluxlength][2];
double extern shortflux[fluxlength][2];

// Forward declare subroutines
double brem_per_Egamma(double, double);
double recombination(double);
void sum_flux(vector<vector<double> >&, int&);
void integrate_flux(vector<vector<double> >, vector<vector<double> >&, int);
void integrate_flux_2(vector<vector<double> >, vector<vector<double> >&, int);
void recomb_flux(vector<vector<double> >, vector<vector<double> >&, int);


/************* Begin main ***************/
int main() {
	ofstream outfile,outfile2;
	outfile.open("bremyield.dat");
	outfile2.open("recombyield.dat");
	vector<vector<double> > recomb;
	vector<vector<double> > radyield;
	vector<vector<double> > flux;
	int count = 0;
	sum_flux(flux, count);
	integrate_flux(flux, radyield, count);
	recomb_flux(flux, recomb, count);
	for (int i = 1; i <= count; i++) {
		outfile << radyield[i - 1][0] << " " << radyield[i - 1][1] << endl;
		outfile2 << recomb[i - 1][0] << " " << recomb[i - 1][1] << endl;
	}
	outfile.close();
	outfile2.close();
	return 0;
}
/************* End main ***************/


/********************************** Define subroutines *******************************************/ 
// brem_per_Egamma gets the Bethe-Heitler bremsstrahlung cross section 
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
	// This conditional checks to see if the photon energy exceeds the rescatter kinetic energy. 
	// In this case, the cross-section is not a number. As the photon energy approaches the kinetic energy, 
	// the cross-section goes to zero.
	if (isnan(BH)) {
		return 1e-300; // Return approximate machine zero
	}
	else {
		return BH;
	}
}


// This function calculates the recombination cross section (currently set for Ar+7 -> Ar+6)
double recombination(double Er) {
	// Define fit parameters for Argon +7
	double E0, sig0, ya, P, yw, y0, y1;
	E0 = 3.884; // eV
	sig0 = 3.295e1; //Mb
	ya = 7.082e2;
	P = 4.645;
	yw = 0.0;
	y0 = 0.0;
	y1 = 0.0;
	// Define fit
	double x, y, F, sigPI, Egamma;
	Egamma = Er + IP;
	x = (Egamma*27.2 / E0) - y0; // convert Egamma to eV
	y = sqrt(x * x + y1 * y1);
	F = (pow(x - 1, 2.0) + yw * yw) * pow(y, 0.5 * P - 5.5) * pow(1 + sqrt(y / ya), -P);
	sigPI = sig0 * F; // Mb
	// Define recombination using Milne relation Ar+7 -> Ar+6
	double sigRC;
	sigRC = sigPI * (2.0 / 1.0) * Egamma * Egamma / (2 * m * c * c * Er); // Mb
	sigRC = sigRC * 1e-22 / pow(5.291e-11, 2.0); // Conversion from Mb -> m^2 -> atomic units (bohr^2)
	return sigRC;
}


// This function takes the sum of the long and short flux while removing the empty bins
void sum_flux(vector<vector<double> >& flux, int& count) {
	for (int i = 1; i <= fluxlength; i++) {
		if ((longflux[i - 1][1] + shortflux[i-1][1]) != double(0.0)) {
			count++;
			flux.push_back({ longflux[i - 1][0], longflux[i - 1][1] + shortflux[i - 1][1] });
		}
	}
}


// This function integrates the rescattering flux to get the photon resolved radiation yield
void integrate_flux(vector<vector<double> > flux, vector<vector<double> >& rad, int count) {
	// Trapezoid rule is used to approximate integral (useful because of the variable step size for the flux)
	for (int i = 1; i <= count; i++) { // Photon energy loop
		double RadYieldSum = 0.0;
		double photonEnergy = flux[i - 1][0];
		for (int j = 1; j < count; j++) { // Rescatter integral

			RadYieldSum += (flux[j][0] - flux[j - 1][0]) * 0.5 * (brem_per_Egamma(flux[j][0], photonEnergy) * flux[j][1] + brem_per_Egamma(flux[j - 1][0], photonEnergy) * flux[j - 1][1]);
		}
		rad.push_back({ flux[i - 1][0], photonEnergy * RadYieldSum });
	}
}


// This function integrates the rescattering flux to get the kinetic energy resolved radiation yield
void integrate_flux_2(vector<vector<double> > flux, vector<vector<double> >& rad, int count) {
	// Trapezoid rule is used to approximate integral (useful because of the variable step size for the flux)
	for (int i = 1; i <= count; i++) { // Rescatter energy loop
		double RadYieldSum = 0.0;
		double RescatterEnergy = flux[i - 1][0];
		for (int j = 1; j < count; j++) { // Photon integral

			RadYieldSum += (flux[j][0] - flux[j - 1][0]) * 0.5 * (flux[j][0]*brem_per_Egamma(RescatterEnergy, flux[j][0]) + flux[j - 1][0]*brem_per_Egamma(RescatterEnergy, flux[j-1][0]));
		}
		rad.push_back({ flux[i - 1][0], flux[i-1][1] * RadYieldSum});
	}
}


// This function calculates the recombination radition yield (currently set for Ar+7 -> Ar+6)
void recomb_flux(vector<vector<double> > flux, vector<vector<double> >& rad, int count) {
	double RadYield, Egamma, Er;
	for (int i = 1; i <= count; i++) {
		Er = flux[i - 1][0];
		Egamma = Er + IP;
		RadYield = Egamma * recombination(Er) * flux[i - 1][1];
		rad.push_back({ Egamma,RadYield });
	}
}



//Define Long Flux array
double longflux[fluxlength][2]{
{0.00148,0},
{0.00185,0},
{0.00232,0},
{0.0029,0},
{0.00363,0},
{0.00455,0},
{0.00569,0},
{0.00713,0},
{0.00893,0},
{0.01118,0},
{0.014,0},
{0.01753,0},
{0.02195,0},
{0.02748,0},
{0.03441,0},
{0.04308,0},
{0.05394,0},
{0.06754,0},
{0.08457,2.72482E-8},
{0.1059,0},
{0.13259,0},
{0.16602,0},
{0.20788,0},
{0.26029,0},
{0.32592,0},
{0.40809,0},
{0.51097,0},
{0.6398,0},
{0.80111,0},
{1.00308,0},
{1.25598,0},
{1.57264,0},
{1.96913,0},
{2.46559,0},
{3.08721,0},
{3.86556,0},
{4.84015,0},
{6.06045,0},
{7.58841,2.67037E-8},
{9.50161,0},
{11.8972,3.34605E-8},
{14.8967,0},
{18.6524,2.80488E-8},
{23.3551,3.14407E-8},
{29.2434,2.82452E-8},
{36.6163,3.07189E-8},
{45.848,3.27143E-8},
{57.4072,3.19017E-8},
{71.8807,3.52946E-8},
{90.0033,3.70406E-8},
{112.695,4.28243E-8},
{141.108,4.92506E-8},
{176.684,6.53759E-8},
{221.23,1.14729E-7},
{277.006,9.13722E-7}
};

//Define Short Flux array
double shortflux[fluxlength][2]{
{0.00148,0},
{0.00185,0},
{0.00232,0},
{0.0029,2.28482E-255},
{0.00363,0},
{0.00455,0},
{0.00569,1.69457E-174},
{0.00713,0},
{0.00893,2.03179E-140},
{0.01118,0},
{0.014,2.88341E-117},
{0.01753,0},
{0.02195,6.58941E-96},
{0.02748,0},
{0.03441,4.05782E-84},
{0.04308,0},
{0.05394,6.74116E-75},
{0.06754,0},
{0.08457,2.69073E-65},
{0.1059,2.33833E-59},
{0.13259,1.91824E-54},
{0.16602,0},
{0.20788,4.95469E-49},
{0.26029,1.97504E-45},
{0.32592,2.18811E-42},
{0.40809,5.9988E-39},
{0.51097,1.55552E-36},
{0.6398,3.2974E-32},
{0.80111,5.46232E-30},
{1.00308,4.84794E-28},
{1.25598,4.03296E-26},
{1.57264,6.42225E-25},
{1.96913,1.75527E-23},
{2.46559,1.91118E-21},
{3.08721,2.48168E-20},
{3.86556,7.99962E-19},
{4.84015,1.11043E-17},
{6.06045,1.14646E-16},
{7.58841,1.67552E-15},
{9.50161,1.01197E-14},
{11.8972,1.00108E-13},
{14.8967,5.47977E-13},
{18.6524,2.94791E-12},
{23.3551,1.50138E-11},
{29.2434,6.85901E-11},
{36.6163,1.6747E-10},
{45.848,5.85686E-10},
{57.4072,1.5979E-9},
{71.8807,4.13527E-9},
{90.0033,1.09032E-8},
{112.695,2.36865E-8},
{141.108,5.55828E-8},
{176.684,9.65616E-8},
{221.23,2.05951E-7},
{277.006,2.52087E-6}
};
