/************************************************************
 * GetFlux.cpp
 *
 * GetFlux.cpp is designed to take calculated rescattering quantities from
 * cutoff.cpp and output the normalized rescattering flux per unit energy.
 * Comments are provided to explain the process and physics step by step.
 *
 * Author: Sui Luo.
 * Revisions by Evan Jones
 * Department of Physics and Astronomy, University of Delaware
 * Last Updated: 1:00 12/19/2021
 ************************************************************/


/*********************** Header ***********************/
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;

/***************** Declare Subroutines ****************/
void procMean(string, string, string, string,string);
void getMean(vector<double>, 
		vector<double>,
		vector<double>,
		vector<double>,
		vector<double>,
		vector<double>,
		double&, double&,
		double&, double&,
		double&, double&);
double getPsiStarPsi(double, double);

/******************** Main Program ********************/
int main()
{
	/*generate input file name&path*/
	string path_str;
	path_str="cutoff_output/"; /* Put here the path to the folder your input data is in */

	/*generate input file*/
	string file_str;
	for (int i = 7; i <= 7; i++) /* The range i runs through should be the range of indexes from your input files. *
	 	 	 	 	 	 	 	   * This is the first integer which appears in the file name after the element.    *
	 	 	 	 	 	 	 	   * If only one charge is desired, set the start i and final i to the file index.  */
	{
		stringstream ind,ion,wavelength,Z;
		ind << i; //ind will be the index value currently being called
		ion << i+1; //ion will be the ion currently being called - for us this was usually 1 greater than the index as we would run ions from 1+ to some maximum
		wavelength << 1280; //input the wavelength you're running
		Z << 36; //input the atomic number of the species calculated
		file_str = path_str + "Z=" + Z.str() + "_data_" + ind.str() + "+" + ion.str() + "_" +  wavelength.str() + "nm.dat";
		procMean(file_str,ind.str(),ion.str(),wavelength.str(),Z.str());
	}

	/*end of the program*/
	return 0;
}

/******************** Subroutines *********************/

/******************************************************/

void procMean(string file_str, string flag, string ion, string wavelength, string Z) {

	/*generate output file and title columns */
	string out_str = "getflux_output/Z=" + Z + "_out_" + flag + "+" + ion + "_" + wavelength + "nm.dat";
	const char *output_file_char = out_str.c_str();
	ofstream output(output_file_char, ios::out);
	output << "release_phase" << " "
			<< "return_phase" << " "
			<< "adk_rate_dt" << " "
			<< "ini_deflection" << " "
			<< "ini_width" << " "
			<< "Wavepacket_deflection" << " "
			<< "Wavepacket_width" << " "
			<< "Gamma_R" << " "
			<< "Return_kinetic" << " "
			<< "Return_psi*psi_ADK_dt" << " "
			<< "dF_R/dE" << endl;

	/*generate target file name with full path&name*/
	const char *input_file_char = file_str.c_str();

	/*open file*/
	ifstream input;
	input.open(input_file_char);

	/*check existance of file*/
	if(!input)
	{
		cout << "Error: file could not be opened >> " + file_str << endl;
		system("pause");
		exit(1);
	}

	/* This vector will serve as a container to hold the calculated results from the rescattering cutoff code.*/
	vector<vector<double> > datas;

	/*temperoral variable*/
	string t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;//title in each file
	double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9;//calculated result in each file

	/*read in titles*/
	input >> t0 >> t1 >> t2 >> t3 >> t4 >> t5 >> t6 >> t7 >> t8 >> t9; /* This reads in the titles for each column for each quantity output
																		  from the rescattering code and sets it to the title "t" variables.*/

	/*variables to record*/
	double tstart = 0.0;
	double rate = 0.0; // ADK_rate times dt
	vector<double> vt, vy, vz, vyi, vzi, vkin;
	bool firstRow = true; // This boolean is used to determine if the first row of the input quantities have already been recorded.


	while(input >> d0 >> d1 >> d2 >> d3 >> d4 >> d5 >> d6 >> d7 >> d8 >> d9) /* This while loop states that while there are quantities to read-in from
																				the file source, it will continue.*/
	{
		if (!firstRow && abs(d0-tstart) > 0.0)/* This conditional is to calculate quantities for a given birth phase and ensure that quantities corresponding to a
		 	 	 	 	 	 	 	 	 	 	 	 	     different birth phase does not contribute. */
		{
			double mean, std, meani, stdi, kin, rePhase;
			getMean(vt,vy,vz,vyi,vzi,vkin,mean,std,meani,stdi,kin,rePhase); /* This function calculates the average and standard for the given phase quantities obtained
			 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	   from the previous birth phase.*/
			/*get new record*/
			vector<double> tmpdata;
			tmpdata.push_back(tstart); // birth phase
			tmpdata.push_back(rePhase); // rescatter time in radians
			tmpdata.push_back(rate); // adk_rate times dt
			tmpdata.push_back(meani); // initial wavepacket deflection (should be close to zero)
			tmpdata.push_back(stdi); // initial wavepacket width
			tmpdata.push_back(mean); // rescattered wavepacket deflection
			tmpdata.push_back(std); // rescattered wavepacket width
			tmpdata.push_back(mean*mean/(2*std*std)); // This is the literal definition of Gamma_R from 10.1103/PhysRevA.74.033403
			tmpdata.push_back(kin); // return kinetic energy
			tmpdata.push_back(getPsiStarPsi(mean,std)*rate); /* rescattered wave function magnitude (psi*psi) calculated at the parent ion (x=y=z=0) times ADK_population times dt.
																this results in the unit electrons per bohr^2.*/
			datas.push_back(tmpdata); // places the vector of quantities into the vector "datas" to be used for writing to the output file

			/*reset*/
			vt.clear();
			vyi.clear();
			vzi.clear();
			vy.clear();
			vz.clear();
			vkin.clear();
		} /* end of if */

		tstart = d0;
		if (firstRow) /* check to see if first row of input quantities has been recorded */
		{
			if (tstart > 0.0)
			{
				firstRow = false;
			}
		}

		/* These variables and vectors record the calculated result for each row so that the average and standard deviation can be
		 * calculated in the getMean() function. */
		rate = d1; // ADK_rate times dt
		vt.push_back(d2); // rescatter time in radians
		vyi.push_back(d4); // initial y
		vzi.push_back(d5); // initial z
		vy.push_back(d7); // return y
		vz.push_back(d8); // return z
		vkin.push_back(d9); // return kinetic energy

	} /* end of while loop */


	/* close input */
	input.close();

	/*get rescattering flux per rescattering kinetic energy (dF_R/dE) and add it to each vector in datas*/
	datas[0].push_back(datas[0][9]/abs(datas[0][8])); // First entry is outside the loop to prevent selecting "datas[i][8]" with a negative integer
	for (int i = 1; i < int(datas.size()); i++) /* For loop takes the rescattered wave function magnitude (psi*psi) calculated at the parent ion (x=y=z=0)
												   times ADK_population times dt and divides it all by the change in rescattering kinetic energy dE. This quantity is
												   the rescattering flux per rescattering kinetic energy (dF_R/dE) in the form of "Eq. 2" from 10.1103/PhysRevLett.118.093001 .*/
	{
		datas[i].push_back(datas[i][9]/abs(datas[i][8]-datas[i-1][8]));
	}


	/*This vector will serve as a container to hold the flux and kinetic energy from the flux code.*/
	vector<vector<double> > longtraj;
	vector<vector<double> > shorttraj;

	double ADK_dt_sum = 0;
	double E_pre=0,E_now;
	/*output all calculated quantities per birth phase*/
	for (int i = 0; i < int(datas.size()); i++)
	{
		ADK_dt_sum += datas[i][2];
		E_now = datas[i][8];
		if(E_now - E_pre > 0.0){ // This conditional defines long trajectories. It assumes that the kinetic energy input is monotonic up to the maximum.
			longtraj.push_back({E_now,datas[i][10]});
			E_pre = E_now;
		}
		else{
			shorttraj.push_back({E_now,datas[i][10]});
		}
		/*This for loop prints to a separate outfile*/
		for (int j = 0; j < int(datas[i].size()); j++)
		{
			output << datas[i][j] << " ";
		}
		output << endl;
	}
	output.close();

	/*normalize and output dF_R/dE by long and short trajectories*/

	/*generate output file for long and short trajectories*/
	string normLong_out_str = "getflux_output/Z=" + Z + "_norm_Long" + flag + "+" + ion + "_" + wavelength + "nm.dat";
	string normShort_out_str = "getflux_output/Z=" + Z + "_norm_Short" + flag + "+" + ion + "_" + wavelength + "nm.dat";
	const char *normLong_output_file_char = normLong_out_str.c_str();
	const char *normShort_output_file_char = normShort_out_str.c_str();
    
	ofstream longfile, shortfile;
	longfile.open(normLong_output_file_char);
	shortfile.open(normShort_output_file_char);
	longfile << "E" << " " << "norm_dF_R/dE" << endl;
	shortfile << "E" << " " << "norm_dF_R/dE" << endl;
	/* This for loop writes the long and short trajectories to their respective files */
	for(int i = 0; i < int(longtraj.size()) + int(shorttraj.size()); i++){
		if (i < int(longtraj.size())) {
			// Check if output is numeric before continuing
			if (isnan(longtraj[i][1] / ADK_dt_sum)) {
				cout << "Warning: A non-numeric long trajectory flux was computed for ion charge " + ion + ". Substituing with zero." << endl;
				//system("pause");
				longfile << longtraj[i][0] << " " << 0 << endl;
			}
			else {
			longfile << longtraj[i][0] << " " << longtraj[i][1] / ADK_dt_sum << endl;
			}
		}
		else{
			// Check if output is numeric before continuing
			if (isnan(shorttraj[(int(shorttraj.size()) - 1) - (i - int(longtraj.size()))][1] / ADK_dt_sum)) {
				cout << "Warning: A non-numeric short trajectory flux was computed for ion charge " + ion + ". Substituing with zero." << endl;
				//system("pause");
				shortfile << shorttraj[(int(shorttraj.size()) - 1) - (i - int(longtraj.size()))][0] << " " << 0 << endl;
			}
			/* The lengthy argument calls each element of the shorttraj vector but in reverse so that the output is increasing in energy */
			else {
				shortfile << shorttraj[(int(shorttraj.size()) - 1) - (i - int(longtraj.size()))][0] << " "
					<< shorttraj[(int(shorttraj.size()) - 1) - (i - int(longtraj.size()))][1] / ADK_dt_sum << endl;
			}
		}
	}/*end for loop*/
	longfile.close();
	shortfile.close();

	/*Print to command prompt that ion is compelete */
	cout << "Ion charge " + ion + " complete" << endl;

	/*end of the subroutine*/
	return;

}

/******************************************************/
// finds the means and standard deviations
void getMean(vector<double> vt,
		vector<double> vy,
		vector<double> vz,
		vector<double> vyi,
		vector<double> vzi,
		vector<double> vkin,
		double &mean, double &std,
		double &meani, double &stdi,
		double &kin, double &rePhase) {

	/*mean time*/
	double mt = 0.0;
	for (int i = 0; i < int(vt.size()); i++)
	{
		mt += vt[i];
	}
	mt /= vt.size();
	rePhase = mt;

	/*mean return y*/
	double my = 0.0;
	for (int i = 0; i < int(vy.size()); i++)
	{
		my += vy[i];
	}
	my /= vy.size();

	/*mean return z*/
	double mz = 0.0;
	for (int i = 0; i < int(vz.size()); i++)
	{
		mz += vz[i];
	}
	mz /= vz.size();

	mean = sqrt(my*my + mz*mz);

	/*std return y*/
	double sqry = 0.0;
	for (int i = 0; i < int(vy.size()); i++)
	{
		sqry += (vy[i]-my)*(vy[i]-my);
	}

	/*std return z*/
	double sqrz = 0.0;
	for (int i = 0; i < int(vz.size()); i++)
	{
		sqrz += (vz[i]-mz)*(vz[i]-mz);
	}

	std = sqrt((sqry+sqrz)/vy.size());

	/*mean initial y*/
	double myi = 0.0;
	for (int i = 0; i < int(vyi.size()); i++)
	{
		myi += vyi[i];
	}
	myi /= vyi.size();

	/*mean initial z*/
	double mzi = 0.0;
	for (int i = 0; i < int(vzi.size()); i++)
	{
		mzi += vzi[i];
	}
	mzi /= vzi.size();

	meani = sqrt(myi*myi + mzi*mzi);

	/*std initial y*/
	double sqryi = 0.0;
	for (int i = 0; i < int(vyi.size()); i++)
	{
		sqryi += (vyi[i]-myi)*(vyi[i]-myi);
	}

	/*std initial z*/
	double sqrzi = 0.0;
	for (int i = 0; i < int(vzi.size()); i++)
	{
		sqrzi += (vzi[i]-mzi)*(vzi[i]-mzi);
	}

	stdi = sqrt((sqryi+sqrzi)/vyi.size());

	/*mean return kin*/
	double mkin = 0.0;
	for (int i = 0; i < int(vkin.size()); i++)
	{
		mkin += vkin[i];
	}
	mkin /= vkin.size();
	kin = mkin;
}

/******************************************************/
//just a comment for this one to explain what is does (E: Barry needs to look at this and decide what we should call it)
double getPsiStarPsi(double mean, double std) {
	double pi = 4.0*atan(1.0);
	double psiStarpsi;
	/*2D Gaussian wave function*/
	psiStarpsi = 1.0/(2.0*pi*std*std)*exp(-mean*mean/(2.0*std*std));

	return psiStarpsi;
}

