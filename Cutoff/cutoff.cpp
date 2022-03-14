/****************************************************************************/
/*
 * cutoff.cpp
 *
 * cutoff.cpp is intended to calculate the rescattering dynamics of a photo-ionized electron
 * wavepacket in a linearly polarized laser field via a combination of analytical quantum
 * tunneling and relativistic monte-carlo trajectory ensembles. Comments are provided to explain
 * the process and physics step by step.
 *
 * Author: Sui Luo.
 * Revisions by Evan Jones
 * Department of Physics and Astronomy, University of Delaware
 * Last Updated: 1:00 12/19/2021
*/
/****************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <vector>


// rksuite
#include "rksuite.h"

// distribution generation
#include "nr3.h"
#include "ran.h"
#include "deviates.h"

using namespace std;

/****************************************************************************/
// CONSTANTS AND USER INPUT
/****************************************************************************/
// charge of electron in atomic units
const double e = 1.0;

// rest mass of electron in atomic units
const double me = 1.0;

// hbar in atomic units
const double hb = 1.0;

// speed of light in atomic units
const double c = 137.03545;

// Pi
const double pi = 4.0*atan(1.0);

// atomic units to W/cm^2 conversion factor
const double int_au = 6.43640931e15;

// wavelength of drive laser
const double wavelength_nm = 720; // Input wavelength here in nanometers
const double wavelength = wavelength_nm*18.897; // Conversion to atomic units (do not change)

// frequency of drive laser
const double freq = (2.0*pi*c)/wavelength;

// period of drive laser
const double period = wavelength/c;

// number of equations of motion corresponding to force and momentum (also momentum and displacement) in three dimensional Cartesian coordinates
const int neq = 6;

/*Define atomic number, ion charges, ionization potentials, angular quantum number (l), and desired laser intensities (W/cm^2 input)*/
const int Z = 36; // atomic number of desired species calculated

const int ionNum[Z] = { // all possible ion charges for desired species from +1 to +Z
		1,
		2,
		3,
		4,
		5,
		6,
		7,
		8,
		9,
		10,
		11,
		12,
		13,
		14,
		15,
		16,
		17,
		18
};

const double ip[Z] = { // ionization potential for each charge of desired species; taken from NIST (https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html)
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	4.625
	/*0.5794,
	1.0158,
	1.49761,
	2.19044,
	2.75147,
	3.35625,
	4.5739,
	5.27414,
	15.53676,
	17.63824,
	19.86765,
	22.75735,
	25.20221,
	27.76213,
	31.45221,
	33.76379,
	151.49506,
	162.72878*/
};

const int l[Z] = { // angular quantum number for each charge of desired species
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0
	/*1,
	1,
	1,
	1,
	1,
	1,
	0,
	0,
	1,
	1,
	1,
	1,
	1,
	1,
	0,
	0,
	0,
	0*/
};

const double intensityVect[Z] = { /*Desired laser intensity for each charge state in W/cm^2. Set to zero if the corresponding ion charge is not desired.
									Typically full chamber saturation intensities are used, which can be calculated theoretically using ADK rate
									population in a laser pulse (see LaserIonizationYield code in LaserScience repository).*/
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		2.2e17
		/*6.3e14,
		1.5e15,
		3.5e15,
		7.9e15,
		1.1e16,
		1.7e16,
		4e16,
		5.8e16,
		2e18,
		2.8e18,
		3.9e18,
		5.5e18,
		7.3e18,
		1e19,
		1.4e19,
		2.2e19,
		2.8e21,
		4e21*/
};

/*These are random seeds associated with the random number generators utilized in the code.
  Can take any integer value but array should not exceed size Z*/
const int seedVect[Z] =
{
		0,
		1,
		2,
		3,
		4,
		5,
		6,
		7,
		8,
		9,
		10,
		11,
		12,
		13,
		14,
		15,
		16,
		17
};

/* Number of monte-carlo generated trajectories sampled to form the wave packet */
int nSample = 50;

/*Time steps for initial ionization phase (birth phase) of laser cycle*/
int BirthSteps = 4000;

/*These are the initial time steps used to find the x-axis. Leave this at around 100.*/
int TrajSteps = 200;

/*Precision of bisection algorithm when solving for xwant. If the algorithm gets x within the absolute value of this number,
 * it is considered to be equal to xwant. */
double bisection_precision = 0.0005;
double xwant = 0.0;

/****************************************************************************/
// END CONSTANTS AND USER INPUT
/****************************************************************************/


/****************************************************************************/
/*global*/
/****************************************************************************/
/*input intensity*/
double int_W_per_cmSqr;
/*intensity in atomic unit*/
double intensity;
/*EM peak amplitude from laser in atomic units*/
double eAmpMax;
/*initial phase of pulse*/
double iniphase;
/*job number*/
int n_job;

/****************************************************************************/
/*subroutine and function declare*/
/****************************************************************************/
double GetADK(const double,
	const int,
	const int,
	const double[],
	const double);

void Derivs(double,
	double[],
	double[]);

void GetEmField(const double,
	const double,
	double*,
	double&);

void GetInitial(const double,
	const double,
	const double,
	const int,
	const int,
	vector<double>*,
	vector<double>*,
	vector<double>*,
	vector<double>*);

int start(int);

string getStringFromNumber(const int);

/****************************************************************************/
/*main program*/
/****************************************************************************/
int main(){
	for(int i = 1; i <= Z; i++){
		if(intensityVect[i-1]!=0){
			start(i);
			cout << "Ion number " << ionNum[i-1] << " finished" << endl;
		}
		else
		{
			cout << "Ion number " << ionNum[i-1] << " had intensity 0" << endl;
		}
	}
	return 0;
}

int start(int ion){

	/*record running time*/
	clock_t elapsed1,elapsed2;
	elapsed1 = clock();

	/*declare population evolution variables*/
	double rate_adk;

	/*declare trajectory integration variables*/
	double y_start[neq], x_pre, x_now, rescatter_time, rescatter_kin, rescatter_Etot;

	/*initialize rksuite parameters*/

	////////////////////////////////////////////////////////////////////////////////////
	/*This code utilizes rksuite to solve for the electron trajectories in the laser field. It is a differential equation solver based on runge-kutta methods.
	Documentation on usage is provided here at "https://www.netlib.org/ode/rksuite/" */
	////////////////////////////////////////////////////////////////////////////////////

	RKSUITE rksuite;
	double y[neq],yp[neq],ymax[neq];
	double twant, tnow;
	double tol = 1e-6; // =1e-6 for method 2
	double thres[neq] = {1e-10,1e-10,1e-10,1e-10,1e-10,1e-10}; //1e-10 for method 2
	double hstart = 0.0; //Gives control of choosing initial integration point
	int method = 2; //Method 2 is most efficient, medium-error
	char jobtype = 'U';
	char* jobpointer = &(jobtype); // RKsuite only takes pointers to characters.
	int uflag; //For error flagging
	bool mesage = false;
	bool errass = false;

	n_job = ion - 1; // This converts the desired ion to the index "n_job", which is used to select the the atomic parameters from their respective arrays.

	/*generate initial phase of pulse, which is distinct from the birth phase. Should be set to zero by default.*/
	iniphase = double(0);

	/*set intensity and EM parameters*/
	int_W_per_cmSqr = intensityVect[n_job]; // selects intensity specific to desired ion
	intensity = int_W_per_cmSqr/int_au; // conversion to atomic units
	eAmpMax = sqrt((8.0*pi*intensity)/c); // electric field amplitude for laser in atomic units

	/*output declare*/
	/*This names the files based on ion index and ionization number*/
	string file_str_Log = "cutoff_output/Z=" + getStringFromNumber(Z) + "_log_" + getStringFromNumber(n_job) + "+" + getStringFromNumber(ionNum[n_job]) + "_" + getStringFromNumber(wavelength_nm) + "nm.txt";
	string file_str_Data = "cutoff_output/Z=" + getStringFromNumber(Z) + "_data_" + getStringFromNumber(n_job)  + "+" + getStringFromNumber(ionNum[n_job]) + "_" + getStringFromNumber(wavelength_nm) + "nm.dat";
	const char *file_char_Log = file_str_Log.c_str();
	const char *file_char_Data = file_str_Data.c_str();
	ofstream outLog(file_char_Log, ios::out);
	ofstream outData(file_char_Data, ios::out);
	if (!outLog) {
		cout << "Could not open directory 'cutoff_output'. Make sure you create the directory" << endl;
		cout << "'cutoff_output' and place it in same directory as the executable." << endl;
		system("pause");
		exit(1);
	};
	if (!outData) {
		cout << "Could not open directory 'cutoff_output'. Make sure you have created the directory" << endl;
		cout << "'cutoff_output' and place it in same directory as the executable." << endl;
		system("pause");
		exit(1);
	};

	/*set time domain parameters*/
	double t_delta = period/BirthSteps/4.0; /* Time step for birth phase. Divides the period of the drive laser by 4 (quarter cycle), then divides that time by the
												number of birth time steps. */
	double t_integ_delta = (period/TrajSteps); /* Time step for trajectory calculation. Divides the period of the drive laser by the number of trajectory time steps. */

	double omega = (2.0*pi*c)/wavelength; // angular frequency of the drive laser
	double a0 = e*eAmpMax/(me*c*omega); // a0 and gamma_R are from the Lorentz force deflection parameter found in "10.1103/PhysRevLett.118.093001"
	double gamma_R = (sqrt(2.0*ip[n_job]*me*pow(c,2.0))*pow(a0,3.0))/(16.0*hb*omega);

	/*declare time domain variable*/
	double t_start = period/4.0 - t_delta; /*This is the initial phase of the laser phase (birth phase) set to the period/4. The first time step is subtracted off
											so that the algorithm begins at the birth phase. */
	double t_final, t_integ_final;

	/*write log of project parameters*/
	outLog << ">>> wavelength(a.u.) = " << wavelength << endl;
	outLog << ">>> period(a.u.) = " << period << endl;
	outLog << ">>> freq(a.u.) = " << freq << endl;
	outLog << ">>> intensity(W/cm^2) = " << int_W_per_cmSqr << endl;
	outLog << ">>> intensity(a.u.) = " << intensity << endl;
	outLog << ">>> eAmpMax(a.u.) = " << eAmpMax << endl;
	outLog << ">>> t_delta(a.u.) = " << t_delta << endl;
	outLog << ">>> t_integ_delta(a.u.) = " << t_integ_delta << endl;
	outLog << ">>> range = " << (t_start + t_delta) / period * 2 * pi << " to " << (BirthSteps * t_delta) / period * 2 * pi << " rad" << endl;
	outLog << ">>> Ion Number: " << ion << endl;
	outLog << ">>> Gamma_r: " << gamma_R << endl;

	/*write title of rescatter data*/
	outData << "birth_phase" << " "	<< "adk_rate_dt" << " " << "return_phase" << " " << "ini_X" << " " << "ini_Y" << " " 
		<< "ini_Z" << " " << "res_X" << " " << "res_Y" << " " << "res_Z" << " " << "return_kin" << endl;

	/*Begin loop over every birth phase (outer loop)*/
	for(int pp = 0; pp < BirthSteps; pp++){

		/*current start*/
		t_start += t_delta;
		t_final = t_start + 1.1*period; /* This is the final integration time setting the duration of the trajectory integration to the period of the laser.*/ // TESTESTSTSETSTSET
		t_integ_final = t_final + t_integ_delta; /*The additional term added to t_final is to prevent sampling times outside
												   of the final time provided to rkuite.setup(), which causes an error
												   (see rksuite documentation https://www.netlib.org/ode/rksuite/ for more details)*/

		/*initialize trajectory setting the initial displacement and momentum to zero*/
		for (int i = 0; i < neq; i++){
			y_start[i] = 0.0;
		}

		/*get ADK rate*/
		rate_adk = GetADK(ip[n_job],ionNum[n_job],l[n_job],y_start,t_start);



		/*initialize vector to store initial position and momentum*/
		vector<double>* iniY = new vector<double>();
		vector<double>* iniZ = new vector<double>();
		vector<double>* iniPy = new vector<double>();
		vector<double>* iniPz = new vector<double>();

		/*get initial conditions for entire set of trajectories*/
		GetInitial(t_start,y_start[5],ip[n_job],nSample,seedVect[n_job],iniY,iniZ,iniPy,iniPz);

		/*Calculate all trajectories per birth phase (inner loop)*/
		for (int i = 0; i < nSample; i++){

			/*get initial conditions for momentum and position*/
			y_start[0] = 0.0;
			y_start[1] = (*iniPy)[i];
			y_start[2] = (*iniPz)[i];
			y_start[3] = 0.0;
			y_start[4] = (*iniY)[i];
			y_start[5] = (*iniZ)[i];

			/*send initial conditions, method, and precision requirements among other things to rksuite
									 (see rksuite documentation https://www.netlib.org/ode/rksuite/ for more details)*/
			rksuite.setup(neq,t_start,y_start,t_integ_final,tol,thres,method,jobpointer,errass,hstart,mesage);

			/*initialize twant to the current birth phase, where twant is the time rksuite will solve for*/
			twant = t_start;

			/*initialize x_pre which is the x position at the beginning of the integration*/
			x_pre = (y_start[3]-xwant);

			/*trajectory integration*/
			/*This while loop integrates until the break is reached or
			* the electron isn't rescattered by t_final
			*/
			while (twant <= t_final){
				/*time step increment*/
				twant += t_integ_delta;

				/*rksuite takes the equations of motion (Derivs) and twant, then returns tnow, y vector, yp (derivative of y) vector, ymax, and uflag.
				 * y is the solution to the differential equations provided in "Derivs".
				 * (see rksuite documentation https://www.netlib.org/ode/rksuite/ for more details)*/
				rksuite.ut(Derivs,twant,tnow,y,yp,ymax,uflag);

				/*conditional if rksuite has an error.*/
				if (uflag > 3){
					cout << "UT uflag > 3 happens: twant = " << twant << endl;
				}

				/*update x_now*/
				x_now = (y[3]-xwant); 

				/*check if rescatter happens by determining if the electron has crossed the x=0 plane*/
				if ((x_pre*x_now) < 0.0){

					/*Begin bisection algorithm*/

					double ta = tnow;
					double tb = tnow - t_integ_delta;
					double xa = (x_now);
					int counter = 0; // This counts the bisection algorithm iterations
					int maxloop = 100; // This is the maximum iterations allowed.
					while(abs(y[3]-xwant) > bisection_precision){ // Set minimum acceptable absolute value in the x-axis
						if(counter==maxloop){ // This conditional is for the event that there is floating point error in the bisection algorithm.
							cout << "A trajectory failed to converge on the x-axis. It will not be included in sample." << endl;
							twant = t_final + t_integ_delta; //Set twant and tnow so that current electron trajectory is skipped and not recorded
							tnow = period + t_integ_delta;
							break;
						}
						else{
							double tx = (ta+tb)/2.0; // Bisect time interval
							rksuite.setup(neq,t_start,y_start,t_integ_final,tol,thres,method,jobpointer,errass,hstart,mesage); // Reset RKsuite
							rksuite.ut(Derivs,tx,tnow,y,yp,ymax,uflag); // Solve for new time
							double xx = (y[3]-xwant); // calculate x at new time 
							if(xa*xx < 0){ // Check if xx crosses x-axis relative to xa
								tb = tx; // Set tb to new time
							}
							else{
								xa = xx;
								ta = tx;
							}
							counter++;
						}
					}

					/*End bisection algorithm*/


					rescatter_time = tnow-t_start;
					/*This will print out the data if the electron rescatters
					* Then, it breaks the trajectory loop
					*/
					if (rescatter_time <= 1.1*period){
						/* relativistic total energy */
						rescatter_Etot = sqrt( (pow(y[0],2.0)+pow(y[1],2.0)+pow(y[2],2.0))*pow(c,2.0) +  pow(me,2.0)*pow(c,4.0) );

						/* relativistic kinetic energy */
						rescatter_kin = rescatter_Etot - me*pow(c,2.0);

						outData << t_start/period*2.0*pi << " "
							<< rate_adk*t_delta<< " " // This output is the ADK rate times the birth phase time delta. This simplifies computation in the PostProcess code.
							<< rescatter_time/period*2.0*pi << " "
							<< y_start[3] << " "
							<< y_start[4] << " "
							<< y_start[5] << " "
							<< y[3] << " "
							<< y[4] << " "
							<< y[5] << " "
							<< rescatter_kin << endl;

						break;
					}
				}

				/*update x_pre*/
				x_pre = x_now;

				/*end of trajectory integration*/
			}

			// This prints ADK rate and Birth phase in the event rescattering does not occur
			if (twant > t_final) {
				outData << t_start / period * 2.0 * pi << " "
					<< rate_adk * t_delta << " " // This output is the ADK rate times the birth phase time delta. This simplifies computation in the PostProcess code.
					<< 0 << " "
					<< 0 << " "
					<< 0 << " "
					<< 0 << " "
					<< 0 << " "
					<< 0 << " "
					<< 0 << " "
					<< 0 << endl;
			}

			/*end of loop through all trajecotries*/
		}


		/*end of birth phase loop*/
	}

	/*here the elapsed code calculation time is printed to the log file*/
	elapsed2 = clock();
	float elapsed_diff = ((float)elapsed2-(float)elapsed1);
	outLog << "The code elapsed time is = "
		<< elapsed_diff/CLOCKS_PER_SEC << "seconds";

	/*close output*/
	outLog.close();
	outData.close();

	/*end of main program*/
	return 0;


}

/****************************************************************************/
/*function adk
  this function returns the adk rate*/
  /****************************************************************************/
double GetADK(const double ip,const int ionNumber,const int l,const double y[neq],const double t){
	/*calculate ADK rate and ionization probability.
	The following derivation follows page 25~27
	in David Neal Fittinghoff's Ph.D. Thesis Dec 1993, University of California.*/

	/*variable declare*/
	double rate_adk, epsilon, flm, nstar, Zeff, nmpower, lstar, c2nl, arg;
	double e_cpn[neq], efield;

	/*get EM components and amplitude*/
	GetEmField(t, y[5], e_cpn, efield); /* This function takes the time point within the period of the optical cycle and position along the z-axis
														and outputs the electromagnetic field vector (e_cpn) and electric field magnitude (efield) */

	/*check if EM == 0*/
	if (efield == 0.0){
		rate_adk = 0.0;
	}
	else{
		/*calculate factors*/
		flm = (2.0 * double(l) + 1.0);
		Zeff = ionNumber; // We set the effective nuclear charge Zeff to the ion number. This creates a sequence from +1 to +Z for charges +0 to +(Z-1) respectively
		nstar = Zeff / sqrt(2.0 * ip);
		nmpower = 2.0 * nstar - 1; // We assume m=0 for all charge states
		lstar = nstar - 1.0;
		epsilon = pow(2.0 * ip, 1.5);
		c2nl = pow(2.0, 2.0 * nstar) / (nstar * tgamma(nstar + lstar + 1.0) * tgamma(nstar - lstar));
		arg = -2.0 * epsilon / (3.0 * abs(efield)); // This is the argument for the exponential factor in the final expression

		/*calculate rate*/
		rate_adk = c2nl * ip * flm * pow(2.0 * epsilon / abs(efield), nmpower) * exp(arg); // outputs ADK rate in atomic units with the dimension "electrons/time"
	}
	return rate_adk;
	/*end of function*/
}

/****************************************************************************/
/*function derivs
  this subroutine defines the equation of motion*/
  /****************************************************************************/
void Derivs(double tgot,double ygot[neq],double ypgot[neq]){
	/*declare EM field variable*/
	double e_cpn[neq];
	double eff;

	/*get EM field components*/
	GetEmField(tgot, ygot[5], e_cpn, eff);

	//calculate gamma factor.
	double gamma = sqrt(pow(me*c,2.0)+pow(ygot[0],2.0)+pow(ygot[1],2.0)+ pow(ygot[2],2.0))/(me*c); /* This derivation of the Lorentz factor is a little unusual as it takes
																					  the magnitude of the momentum as opposed to the speed, but it is a useful
																					  formulation as it allows one to avoid calculating the change in velocity.
																					  Generally speaking, with relativistic dynamics it is often much easier
																					  to deal with the change in momentum as the relativistic acceleration
																					  equations are lengthy and difficult to work with. This form can be obtained
																					  by using the relativistic definition of energy and momentum, and the
																					  momentum-energy relation. */
	//define equation of motion - relativistic

	ypgot[3] = ygot[0]/(me*gamma); // dx/dt = vx = px/(m0*gamma)
	ypgot[4] = ygot[1]/(me*gamma); // dy/dt = vy = py/(m0*gamma)
	ypgot[5] = ygot[2]/(me*gamma); // dz/dt = vz = pz/(m0*gamma)
	ypgot[0] = -e*e_cpn[0]-e*(ygot[1]*e_cpn[5]-ygot[2]*e_cpn[4])/(me*gamma); // dp_x/dt = f_x  /* The velocity in the Lorentz force has been substituted for the relativistic momentum */
	ypgot[1] = -e*e_cpn[1]-e*(ygot[2]*e_cpn[3]-ygot[0]*e_cpn[5])/(me*gamma); // dp_y/dt = f_y
	ypgot[2] = -e*e_cpn[2]-e*(ygot[0]*e_cpn[4]-ygot[1]*e_cpn[3])/(me*gamma); // dp_z/dt = f_z

/*end of subroutine*/
}

/****************************************************************************/
/*function emfield
  this subroutine calculates the EM field*/
  /****************************************************************************/
void GetEmField(const double t,const double zz,double* e_cpn,double& eff){
	/*The EM wave in it's default condition provided is a homogeneous plane wave E(t). If it is desired
	 * to make it a traveling plane wave E(z,t), simply uncomment the "-wn*zz" term.*/

	 //double wn = 2*pi/wavelength;
	eff = eAmpMax*sin(freq*t/*-wn*zz*/);
	/*calculate EM components*/
	e_cpn[0] = eff; // Electric field is directed in the +x-axis
	e_cpn[1] = 0.0;
	e_cpn[2] = 0.0;
	e_cpn[3] = 0.0;
	e_cpn[4] = e_cpn[0]/c; // Magnetic field is directed in the +y-axis and is equal to E/c
	e_cpn[5] = 0.0;
	/*end of subroutine*/
}

/****************************************************************************/
/*function initial condition
  this function returns the initial position and momentum*/
  /****************************************************************************/
void GetInitial(const double t,const double zz,const double ip,const int nSample,const int seed,vector<double>* iniY,
	vector<double>* iniZ,vector<double>* iniPy,vector<double>* iniPz){

	/*calculate EM field*/
	double e_cpn[neq], eff;
	GetEmField(t, zz, e_cpn, eff);

	/*calculate the spatial uncertainty width from the ionized electron tranverse momentum spectrum */
	double yz_width = sqrt( (hb / (2 * e * sqrt(2.0 * me))) * sqrt(ip) / abs(e_cpn[0]) );

	/*declare normal distribution generator*/
	Normaldev ng1(0.0,yz_width,seed); // Spacial distribution
	Normaldev ng2(0.0,hb*0.5/yz_width,seed); // Momentum distribution

	/*declare random generator*/
	Ran myran(seed);

	/*declare spatial and momentum variables*/
	double delz,dely,delPz,delPy;

	for (int i = 0; i < nSample; i++){
		delz = ng1.dev(); // This is a random sample from the z distribution
		dely = ng1.dev(); // This is a random sample from the y distribution
		iniZ->push_back(delz);
		iniY->push_back(dely);

		delPz = ng2.dev();
		delPy = ng2.dev();
		iniPz->push_back(delPz); // This is a random sample from the momentum distribution in z
		iniPy->push_back(delPy); // This is a random sample from the momentum distribution in y
	}

	/*end of function*/
}

/****************************************************************************/
/*this function return string of input number*/
/****************************************************************************/
string getStringFromNumber(const int n){

	stringstream ss;
	ss << n;
	return ss.str();

	/*end of function*/
}
