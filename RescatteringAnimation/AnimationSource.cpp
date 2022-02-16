/*
* Rescattering code for animation figure
*/
#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include<random>
#include<string>
#include"rksuite.h"

using namespace std;

// Initialize global random number generator engine
default_random_engine generator;

// Constants //

// charge of electron in atomic units
const double e = 1.0;
// rest mass of electron in atomic units
const double me = 1.0;
// hbar in atomic units
const double hb = 1.0;
// speed of light in atomic units
const double c = 137.03545;
// Pi
const double pi = 4.0 * atan(1.0);


// number of equations of motion corresponding to force and momentum (also momentum and displacement) in three dimensional Cartesian coordinates
const int neq = 6;


// wavelength of drive laser
const double wavelength_nm = 409; // Input wavelength here in nanometers
const double wavelength = wavelength_nm * 18.897; // Conversion to atomic units (do not change)
const double freq = (2.0 * pi * c) / wavelength;
const double period = wavelength / c;


// intensity of laser
const double intensity = 2.2e18; // W/cm^2
const double int_au = 6.43640931e15; // atomic units to W/cm^2 conversion factor (do not change)
const double Emag = sqrt((8.0 * pi * (intensity / int_au)) / c); // conversion to atomic units field


/*Define atomic number and ionization potential*/
const int Z = 36; // atomic number of desired species calculated
const double IP = 4.625; // hartree


// electrons sampled per wavepacket 
const int nSample = 100;
// Number of phase steps 
const int phaseSteps = 10;


// Declare subroutines //
void PhaseIterator(); 
void wavepacket(double, double, bool, vector<vector<double>>&, vector<vector<double>>&);
void propagator(double, double, double[neq], double[neq]);
void Derivs(double, double[neq], double[neq]);
void GetEmField(double, double* , double&);
void GetInitial(double, double, int, vector<double>*,
	vector<double>*, vector<double>*, vector<double>*, vector<double>*, vector<double>*);


/********************* Main function ***************************/
/////////////////////////////////////////////////////////////////
int main() {
	PhaseIterator();
	return 0;
}
//////////////////////////////////////////////////////////////////
/****************************************************************/


// Define subroutines //


void PhaseIterator() {
	vector<vector<vector<double>>> frame1(phaseSteps), frame2(phaseSteps);
	ofstream outfile;
	string tmp;
	double time_stepsize = period / double(phaseSteps);
	double t = 0.0;
	// Begin Iterating over phase steps
	for (int i = 1; i <= phaseSteps; i++) {
		t += time_stepsize;
		// open image data file
		tmp = "frame" + to_string(i) + ".dat";
		const char* name = tmp.c_str();
		outfile.open(name);
		// clear previous frame data
		frame2.clear();
		frame2.resize(phaseSteps);
		// Iterate over all wavepackets (create a new wave packet per frame)
		for (int j = 1; j <= i; j++) {	
			if (j == i) {
				wavepacket(t, time_stepsize, true, frame1[j - 1], frame2[j - 1]);
			}
			else {
				wavepacket(t, time_stepsize, false, frame1[j - 1], frame2[j - 1]);
				//cout << i << endl;
				//cout << frame2[j - 1].size() << endl;
				//system("pause");
			}
		}
		// Output to data file
		// Loop over all samples per wavepacket (rows)
		for (int j = 1; j <= nSample; j++) {
			// Loop over all wavepackets (columns)
			for (int k = 1; k <= i; k++) {
				// output (x,z) coordinate
				//cout << frame2.size() << endl;
				//cout << frame2[k - 1].size() << endl;
				//cout << frame2[k - 1][j - 1].size() << endl;
				//system("pause");
				outfile << frame2[k - 1][j - 1][3] << " " << frame2[k - 1][j - 1][5] << " ";
			}
			outfile << endl;
		}
		// Set the current frame to the initial frame next phase step
		frame1 = frame2;
		// Close current frame data file
		outfile.close();
	} // end phase step loop
}


// This calculates every electron in the wave packet
void wavepacket(double tstart, double tdelta, bool firststep, vector<vector<double>>& wavepacket_initial, vector<vector<double>>& wavepacket_final) {
	// Initialize wave packet if just ionized
	if (firststep) {
		vector<double> X, Y, Z, Px, Py, Pz;
		GetInitial(tstart, IP, nSample, &X, &Y, &Z, &Px, &Py, &Pz);
		for (int i = 1; i <= nSample; i++) {
			wavepacket_final.push_back({ Px[i -1 ],Py[i - 1],Pz[i - 1],X[i - 1],Y[i - 1],Z[i - 1] });
		}
	}
	else {
		for (int i = 1; i <= nSample; i++) {
			double initial[neq]{ wavepacket_initial[i - 1][0],wavepacket_initial[i - 1][1],wavepacket_initial[i - 1][2],
			wavepacket_initial[i - 1][3], wavepacket_initial[i - 1][4], wavepacket_initial[i - 1][5] };
			double final[neq] = { 0 };
			propagator(tstart, tdelta, initial, final);
			wavepacket_final.push_back({ final[0],final[1],final[2],final[3],final[4],final[5] });
		}
	}
}


void GetInitial(double t, double ip, int nSample, vector<double>* iniX, vector<double>* iniY,
	vector<double>* iniZ, vector<double>* iniPx, vector<double>* iniPy, vector<double>* iniPz) {
	if ((t == 0.0) || (t == period / 2.0) || (t == period)) {
		for (int i = 0; i < nSample; i++) {
			iniZ->push_back(0);
			iniY->push_back(0);
			iniPz->push_back(0); // This is a random sample from the momentum distribution in z
			iniPy->push_back(0); // This is a random sample from the momentum distribution in y
			iniX->push_back(0);
			iniPx->push_back(0);
		}
	}
	else {
		/*calculate EM field*/
		double e_cpn[neq], eff;
		GetEmField(t, e_cpn, eff);

		/*calculate the spatial uncertainty width from the ionized electron tranverse momentum spectrum */
		double yz_width = sqrt(hb / (2 * e * sqrt(2.0 * me)) * sqrt(ip) / abs(e_cpn[0]));

		/*calculate the spatial uncertainty width from the ionized electron in the direction of polarization momentum spectrum */
		double gamma_k = freq * sqrt(2 * me * ip) / (e * abs(e_cpn[0]));
		double x_width = sqrt(hb * pow(gamma_k, 3.0) / (12 * me * freq));

		/*declare normal distribution generator*/
		normal_distribution<double> zy_position(0.0, yz_width);
		normal_distribution<double> zy_momentum(0.0, hb * 0.5 / yz_width);
		normal_distribution<double> x_position(0.0, x_width);
		normal_distribution<double> x_momentum(0.0, hb * 0.5 / x_width);

		/*declare spatial and momentum variables*/
		double delx, delz, dely, delPx, delPz, delPy;

		for (int i = 0; i < nSample; i++) {
			delz = zy_position(generator); // This is a random sample from the z distribution
			dely = zy_position(generator); // This is a random sample from the y distribution
			iniZ->push_back(delz);
			iniY->push_back(dely);

			delPz = zy_momentum(generator);
			delPy = zy_momentum(generator);
			iniPz->push_back(delPz); // This is a random sample from the momentum distribution in z
			iniPy->push_back(delPy); // This is a random sample from the momentum distribution in y

			delx = x_position(generator); // This is a random sample from the x distribution
			iniX->push_back(delx);

			delPx = x_momentum(generator);
			iniPx->push_back(delPx);
		}
	}
	/*end of function*/
}


// This propogates an electron in the laser field for a given time step
void propagator(double tstart, double tdelta, double yinitial[neq], double ygot[neq]) {
	RKSUITE rksuite;
	double ypgot[neq], ymax[neq];
	double tgot;
	double tol = 1e-6; // =1e-6 for method 2
	double thres[neq] = { 1e-10,1e-10,1e-10,1e-10,1e-10,1e-10 }; //1e-10 for method 2
	double hstart = 0.0; //Gives control of choosing initial integration point
	int method = 2; //Method 2 is most efficient, medium-error
	char jobtype = 'U';
	char* jobpointer = &(jobtype); // RKsuite only takes pointers to characters.
	int uflag; //For error flagging
	bool mesage = false;
	bool errass = false;
	// Setup RKsuite
	rksuite.setup(neq, tstart, yinitial, tstart + tdelta, tol, thres, method, jobpointer, errass, hstart, mesage);
	// Solve for tdelta
	rksuite.ut(Derivs, tstart + tdelta, tgot, ygot, ypgot, ymax, uflag);
}


void GetEmField(double t, double* e_cpn, double& eff) {
	/*The EM wave in it's default condition provided is a homogeneous plane wave E(t). If it is desired
	 * to make it a traveling plane wave E(z,t), simply uncomment the "-wn*zz" term.*/

	 //double wn = 2*pi/wavelength;
	eff = Emag * sin(freq * t);
	/*calculate EM components*/
	e_cpn[0] = eff; // Electric field is directed in the +x-axis
	e_cpn[1] = 0.0;
	e_cpn[2] = 0.0;
	e_cpn[3] = 0.0;
	e_cpn[4] = e_cpn[0] / c; // Magnetic field is directed in the +y-axis and is equal to E/c
	e_cpn[4] = 0.0;
	e_cpn[5] = 0.0;
	/*end of subroutine*/
}


// this subroutine defines the equation of motion
void Derivs(double tgot, double ygot[neq], double ypgot[neq]) {
	/*declare EM field variable*/
	double e_cpn[neq];
	double eff;

	/*get EM field components*/
	GetEmField(tgot, e_cpn, eff);

	//calculate gamma factor.
	double gamma = sqrt(pow(me * c, 2.0) + pow(ygot[0], 2.0) + pow(ygot[1], 2.0) + pow(ygot[2], 2.0)) / (me * c);

	ypgot[3] = ygot[0] / (me * gamma); // dx/dt = vx = px/(m0*gamma)
	ypgot[4] = ygot[1] / (me * gamma); // dy/dt = vy = py/(m0*gamma)
	ypgot[5] = ygot[2] / (me * gamma); // dz/dt = vz = pz/(m0*gamma)
	ypgot[0] = -e * e_cpn[0] - e * (ygot[1] * e_cpn[5] - ygot[2] * e_cpn[4]) / (me * gamma); // dp_x/dt = f_x  /* The velocity in the Lorentz force has been substituted for the relativistic momentum */
	ypgot[1] = -e * e_cpn[1] - e * (ygot[2] * e_cpn[3] - ygot[0] * e_cpn[5]) / (me * gamma); // dp_y/dt = f_y
	ypgot[2] = -e * e_cpn[2] - e * (ygot[0] * e_cpn[4] - ygot[1] * e_cpn[3]) / (me * gamma); // dp_z/dt = f_z

/*end of subroutine*/
}


