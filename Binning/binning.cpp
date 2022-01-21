/*
 * binning.cpp
 *
 *  Created on: Dec 23, 2021
 *      Author: Evan
 * Department of Physics and Astronomy, University of Delaware
 * Last Updated: 2:00 12/23/2021
 *
 * binning.cpp takes the calculated rescattering flux per unit energy from
 * GetFlux.cpp code and performs a bin-and-interpolate algorithm on it for
 * presentation purposes. Comments are provided to explain the process
 * step by step.
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

void LogBinInterpolate(string,string,string,string,string,double,int);
void LinearBin(string,string,string,string,string,double,int);

/******************** Main Program ********************/
int main()
{
	/*generate input file name&path*/
	string path_str;
	path_str="getflux_output/"; /* Put here the path to the folder your input data is in */

	/*generate input file*/
	string file_str;
	for (int i = 0; i <= 53; i++) /* The range i runs through should be the range of indexes from your input files. *
								   * This is the first integer which appears in the file name after the element.    *
								   * If only one charge is desired, set the start i and final i to the file index.  */
	{
		stringstream ind,ion,wavelength,Z;
		int Steps;
		double logStepFactor, linearStepSize;
		bool logbin;
		ind << i; //ind will be the index value currently being called
		ion << i+1; //ion will be the ion currently being called - for us this was usually 1 greater than the index as we would run ions from 1+ to some maximum
		Z << 54; //input the atomic number of the species calculated
		wavelength << 300; //input the wavelength you're running

		/*********************************************************************************************************/
		/* Use either of the following options to determine step size. If using one, set the others to 0. */
		/* Note that selecting a number of steps that is too small can decrease the precision and */
		/* smoothness of the interpolation algorithm.                                                            */
		/*********************************************************************************************************/
		logbin = true; // Select if log binning or linear binning is to be performed
		Steps = 0; // Select the number of bin steps desired per ion. This separates the minimum and maximum energy per ion by the given number of steps.
		/* OR if log binning */
		logStepFactor = 1.25; /*input the log step increase factor in binning desired. i.e., if you want energy bins separated by a factor of two, input "2".
		 	 	 	 	 	* This will be the same across all ions. Note that the factor must be greater than one. */
		/* OR if linear binning */
		linearStepSize = 10; /*input the linear step size in binning desired */
		/**********************************************************************************************************/
		if(logbin == true){
			LogBinInterpolate(path_str,ind.str(),ion.str(),wavelength.str(),Z.str(),logStepFactor,Steps);
		}
		else{
			LinearBin(path_str,ind.str(),ion.str(),wavelength.str(),Z.str(),linearStepSize,Steps);
		}

	}

	/*end of the program*/
	return 0;
}

/******************** Subroutines *********************/

/******************************************************/

// This subroutine calculates log-spaced bins and averages the flux within them. Then, it performs a linear interpolation 
// algorithm between the bins.
void LogBinInterpolate(string path_str, string flag, string ion, string wavelength, string Z, double LogStepFactor_input, int LogSteps_input){
	/*temperoral variable*/
	string t0, t1;//title in each file
	double d0, d1;//calculated result in each file

	/*generate target file name with full path&name*/
	string file_str_long = path_str + "Z=" + Z + "_norm_Long" + flag + "+" + ion + "_" +  wavelength + "nm.dat";
	string file_str_short = path_str + "Z=" + Z + "_norm_Short" + flag + "+" + ion + "_" +  wavelength + "nm.dat";
	const char *input_file_char_long = file_str_long.c_str();
	const char *input_file_char_short = file_str_short.c_str();

	/*open file*/
	ifstream input_long;
	ifstream input_short;
	input_long.open(input_file_char_long);
	input_short.open(input_file_char_short);


	/*check existance of files*/
	if(!input_long)
	{
		cout << "Error: file could not be opened >> " + file_str_long << endl;
		system("pause");
		exit(1);
	}
	if(!input_short)
	{
		cout << "Error: file could not be opened >> " + file_str_short << endl;
		system("pause");
		exit(1);
	}


	/*This vector will serve as a container to hold the flux and kinetic energy from the flux code.*/
	vector<vector<double> > longtraj;
	vector<vector<double> > shorttraj;

	/*read in titles*/
	input_long >> t0 >> t1; /* This reads in the titles for each column for each quantity output
 	 	 	 	 	 	 	 from the rescattering code and sets it to the title "t" variables.*/
	input_short >> t0 >> t1;


	/*This takes calculated result from long trajectory file. The maximum and minimum kinetic energy from both long and short
	 * trajectories (not individually, but together) is recorded. */
	double E_min=1e150, E_max=0;
	while(input_long >> d0 >> d1)/* Read in long trajectories*/
	{
		if(d0<E_min){
			E_min = d0;
		}
		if(d0>E_max){
			E_max = d0;
		}
		longtraj.push_back({d0,d1});
	} /*end file read loop*/
	input_long.close();
	while(input_short >> d0 >> d1)/* Read in short trajectories*/
	{
		if(d0<E_min){
			E_min = d0;
		}
		if(d0>E_max){
			E_max = d0;
		}
		shorttraj.push_back({d0,d1});
	}/*end file read loop*/
	input_short.close();


	/**************************************************************************/
	// Begin log binning algorithm
	/**************************************************************************/
	/* Define log bins */
	vector<vector<double> > longbin;
	vector<vector<double> > shortbin;
	int i,s_counter,s_counter_pre=0,l_counter,l_counter_pre=0,k;
	double etmp1,etmp2,emean,shorttmp_sum,longtmp_sum; // These are temporary variables
	double kappa = log(E_max/E_min); // This is the constant that determines the rate the exponentially sized bins grow


	/* Check if log steps or step size are specified */
	int LogSteps;
	if((LogStepFactor_input > 0)&&(LogSteps_input > 0)){
		cout << "Error: Both inputs specifying log steps size are non-zero. Please select one of the two options provided " << endl;
		cout << "in the code while setting the other to zero." << endl;
		system("pause");
		exit(1);
	}
	else if(LogStepFactor_input > 0){
		LogSteps = kappa/log(LogStepFactor_input);
	}
	else if(LogSteps_input > 0){
		LogSteps = LogSteps_input;
	}
	else{
		cout << "Error: Both inputs specifying log steps size are zero. Please select one of the two options provided " << endl;
		cout << "in the code while setting the other to zero." << endl;
		system("pause");
		exit(1);
	}
	/* End of check if log steps or step size are specified */

// This for loop iterates over all bins and fills them with the average flux per bin
	for(i=0; i<=LogSteps; i++){

		/* energy bounds */
		etmp1 = E_min*exp(kappa*(double(i)-0.5)/double(LogSteps)); // left bound of energy bin
		etmp2 = E_min*exp(kappa*(double(i)+0.5)/double(LogSteps)); // right bound of energy bin

		/* short steps */
		k=0; // this integer counts the energies in each bin
		shorttmp_sum = 0.0; // this takes the sum of the flux in each bin
		s_counter = s_counter_pre;
		if(s_counter < int(shorttraj.size())){ // Check to see if there are more energies to check from input file
			while(shorttraj[s_counter][0] < etmp2){ // Check to see if energy exceeds right bound of energy bin
				shorttmp_sum += shorttraj[s_counter][1];
				k++;
				s_counter++; // iterate to the next energy from input file
				if(s_counter >= int(shorttraj.size())){ // This conditional breaks the loop if the index exceeds the array size
					break;
				}
			}
			s_counter_pre = s_counter; // Save the last index to be checked next bin
			if(k!=0){
				emean = E_min*exp(kappa*(double(i))/double(LogSteps)); // center of bin in log representation
				shortbin.push_back({emean,shorttmp_sum/double(k)}); // calculate the average flux within short bin
			}
			else{
				emean = E_min*exp(kappa*(double(i))/double(LogSteps)); // center of bin in log representation
				shortbin.push_back({emean,0.0}); // give a flux of zero to avoid dividing by zero (with k)
			}
		}

		/* long steps */
		k=0; // this integer counts the energies in each bin
		longtmp_sum = 0.0; // this takes the sum of the flux in each bin
		l_counter = l_counter_pre;
		if(l_counter < int(longtraj.size())){ // Check to see if there are more energies to check from input file
			while(longtraj[l_counter][0] < etmp2){ // Check to see if energy exceeds right bound of energy bin
				longtmp_sum += longtraj[l_counter][1];
				k++;
				l_counter++; // iterate to the next energy from input file
				if(l_counter >= int(longtraj.size())){
					break;
				}
			}
			l_counter_pre = l_counter; // Save the last index to be checked next bin
			if(k!=0){
				emean = E_min*exp(kappa*(double(i))/double(LogSteps)); // center of bin in log representation
				longbin.push_back({emean,longtmp_sum/double(k)}); // calculate the average flux within long bin
			}
			else{
				emean = E_min*exp(kappa*(double(i))/double(LogSteps)); // center of bin in log representation
				longbin.push_back({emean,0.0}); // give a flux of zero to avoid dividing by zero (with k)
			}
		}
	}// End for loop


	/*generate bin output file */
	string binnedLong_out_str =  "binning_output/Z=" + Z + "_binned_Long" + flag + "+" + ion + "_" + wavelength + "nm.dat";
	string binnedShort_out_str =  "binning_output/Z=" + Z + "_binned_Short" + flag + "+" + ion + "_" + wavelength + "nm.dat";
	const char *binnedLong_output_file_char = binnedLong_out_str.c_str();
	const char *binnedShort_output_file_char = binnedShort_out_str.c_str();

	ofstream longfile, shortfile;
	longfile.open(binnedLong_output_file_char);
	shortfile.open(binnedShort_output_file_char);

	/* create title columns */
	longfile << "E" << " " << "dF_R/dE" << endl;
	shortfile << "E" << " " << "dF_R/dE" << endl;

	/* output binned calculation */
	for(i=0;i<int(shortbin.size());i++){ // This could also be longbin.size(). Post-binning they are the same size.
		longfile << longbin[i][0] << " " << longbin[i][1] << endl;
		shortfile << shortbin[i][0] << " " << shortbin[i][1] << endl;
	}
	longfile.close();
	shortfile.close();


	/**************************************************************************/
	// Begin linear interpolation algorithm for log-log plot
	/**************************************************************************/
	// This algorithm is designed to fill the empty bins created in the
	// previous binning algorithm.
	/**************************************************************************/
	/*
	// short trajectories
	int jprev = 0; // jprev and j are used to denote filled bins with machine-zero between them
	int j = 0;
	i=0; // This is the primary iterator to find non-machine-zero bins
	bool firstflux = true;

	// This while loop fills the empty bins corresponding to the short trajectories 
	while(i<int(shortbin.size())){ // This checks if the primary iterator exceeds shortbin size
		if(shortbin[i][1] > 1e-300){ // find a flux value greater than machine zero
			if(firstflux){ // This conditional determines if this is the first non-zero flux
				jprev = i; // set jprev to known filled bin
				j=jprev+1; // iterate j one step more than jprev
				if(j<int(shortbin.size())){ // check if j index is larger than shortbin size
					while(shortbin[j][1] <= 1e-300){ // Find the next non-machine-zero flux from i
						if(j<int(shortbin.size())){ // check if j index is larger than shortbin size (within loop)
							j++;
						}
						else{
							j=j-1; // if j index is larger than shortbin size, go back one step (within loop)
						}
					}
				}
				else{
					j=j-1; // if j index is larger than shortbin size, go back one step
				}
				// calculate log-log point 1 
				double etmp1 = log(shortbin[jprev][0]);
				double ytmp1 = log(shortbin[jprev][1]);
				// calculate log-log point 2 
				double etmp2 = log(shortbin[j][0]);
				double ytmp2 = log(shortbin[j][1]);
				// calculate rate of change 
				double m = (ytmp2-ytmp1)/(etmp2-etmp1);
				// calculate constant term 
				double b = ytmp2 - m*etmp2;
				// calculate intermediary points between i and j. If there aren't any, than this loop is skipped
				for(int k=jprev+1; k<j; k++){
					shortbin[k][1] = exp(b + m*log(shortbin[k][0]));
				}
				i=j+1; // iterate i one step more than j
				jprev=j; // set jprev to be used after first flux is found
				firstflux = false; // set firstflux to false so that this conditional is skipped after the first non-machine-zero flux is found
			}
			// If the first flux has already been found, algorithm will skip the conditional above
			else{
				if(shortbin[i-1][1] > 1e-300){ // Check to see if the previous bin is already filled
					jprev=i; // set jprev to current i
					i++;
				}
				else{ // This is in the event a non-machine-zero bin is found with zeros between it and the last known filled bin (jprev)
					j=i; // set j to current i.
					// calculate log point 1 
					double etmp1 = log(shortbin[jprev][0]);
					double ytmp1 = log(shortbin[jprev][1]);
					// calculate log point 2 
					double etmp2 = log(shortbin[j][0]);
					double ytmp2 = log(shortbin[j][1]);
					// calculate rate of change 
					double m = (ytmp2-ytmp1)/(etmp2-etmp1);
					// calculate constant term 
					double b = ytmp2 - m*etmp2;
					// calculate intermediary points between jprev and j
					for(int k=jprev+1; k<j; k++){
						shortbin[k][1] = exp(b + m*log(shortbin[k][0]));
					}
					i=j+1; // iterate i one step more than j
					jprev=j; // set jprev to j
				}
			}
		}
		else{
			i++; // algorithm will continue iterating to find next non-machine-zero bin
		}
	} // end short trajectories while loop


	// long trajectories
	jprev = 0; // jprev and j are used to denote filled bins with machine-zero between them
	j=0;
	i=0; // This is the primary iterator to find non-machine-zero bins
	firstflux = true;

	// This while loop fills the empty bins corresponding to the short trajectories 
	while(i<int(longbin.size())){ // This checks if the primary iterator exceeds shortbin size
		if(longbin[i][1] > 1e-300){ // find a flux value greater than zero
			// This conditional determines if this is the first non-zero flux
			if(firstflux){ // This conditional determines if this is the first non-zero flux
				jprev = i; // set jprev to known filled bin
				j=jprev+1; // iterate j one step more than jprev
				if(j<int(longbin.size())){ // check if j index is larger than shortbin size
					while(longbin[j][1] <= 1e-300){ // Find the next non-machine-zero flux
						if(j<int(longbin.size())){ // check if j index is larger than shortbin size (within loop)
							j++;
						}
						else{
							j=j-1; // if j index is larger than shortbin size, go back one step (within loop)
						}
					}
				}
				else{
					j=j-1; // if j index is larger than shortbin size, go back one step
				}
				// calculate log point 1 
				double etmp1 = log(longbin[jprev][0]);
				double ytmp1 = log(longbin[jprev][1]);
				// calculate log point 2 
				double etmp2 = log(longbin[j][0]);
				double ytmp2 = log(longbin[j][1]);
				// calculate rate of change 
				double m = (ytmp2-ytmp1)/(etmp2-etmp1);
				// calculate constant term 
				double b = ytmp2 - m*etmp2;
				// calculate intermediary points between i and j
				for(int k=jprev+1; k<j; k++){
					longbin[k][1] = exp(b + m*log(longbin[k][0]));
				}
				i=j+1; // iterate i one step more than j
				jprev=j; // set jprev
				firstflux = false; // set firstflux to false so that this conditional is skipped
			}
			// If the first flux has already been found, algorithm will skip the conditional above
			else{
				if(longbin[i-1][1] > 1e-300){ // Check to see if the previous bin is already filled
					jprev=i; // set jprev to current i
					i++;
				}
				else{ // This is in the event a non-machine-zero bin is found with zeros between it and the last known filled bin (jprev)
					j=i; // set j to current i.
					// calculate log point 1 
					double etmp1 = log(longbin[jprev][0]);
					double ytmp1 = log(longbin[jprev][1]);
					// calculate log point 2 
					double etmp2 = log(longbin[j][0]);
					double ytmp2 = log(longbin[j][1]);
					// calculate rate of change 
					double m = (ytmp2-ytmp1)/(etmp2-etmp1);
					// calculate constant term 
					double b = ytmp2 - m*etmp2;
					// calculate intermediary points between jprev and j
					for(int k=jprev+1; k<j; k++){
						longbin[k][1] = exp(b + m*log(longbin[k][0]));
					}
					i=j+1; // iterate i one step more than j
					jprev=j; // set jprev to j
				}
			}
		}
		else{
			i++; // algorithm will continue iterating to find next non-machine-zero bin
		}
	} // end long trajectories while loop


	// generate interpolate output file 
	string interpolateLong_out_str = "binning_output/Z=" + Z + "_interpolate_Long" + flag + "+" + ion + "_" + wavelength + "nm.dat";
	string interpolateShort_out_str = "binning_output/Z=" + Z + "_interpolate_Short" + flag + "+" + ion + "_" + wavelength + "nm.dat";
	const char *interpolateLong_output_file_char = interpolateLong_out_str.c_str();
	const char *interpolateShort_output_file_char = interpolateShort_out_str.c_str();

	longfile.open(interpolateLong_output_file_char);
	shortfile.open(interpolateShort_output_file_char);

	// create title columns 
	longfile << "E" << " " << "dF_R/dE" << endl;
	shortfile << "E" << " " << "dF_R/dE" << endl;

	// output interpolated calculation 
	for(i=0;i<int(shortbin.size());i++){ // Developers note: This could also be longbin.size(). Post-binning they are the same size.
		longfile << longbin[i][0] << " " << longbin[i][1] << endl;
		shortfile << shortbin[i][0] << " " << shortbin[i][1] << endl;
	}
	longfile.close();
	shortfile.close();
	*/

	/*Print to command prompt that ion is compelete */
	cout << "Ion charge " + ion + " complete" << endl;
}


// This subroutine calculates log-spaced bins and averages the flux within them.
void LinearBin(string path_str, string flag, string ion, string wavelength, string Z, double LinStepSize_input, int LinSteps_input){
	/*temperoral variable*/
	string t0, t1;//title in each file
	double d0, d1;//calculated result in each file

	/*generate target file name with full path&name*/
	string file_str_long = path_str + "Z=" + Z + "_norm_Long" + flag + "+" + ion + "_" +  wavelength + "nm.dat";
	string file_str_short = path_str + "Z=" + Z + "_norm_Short" + flag + "+" + ion + "_" +  wavelength + "nm.dat";
	const char *input_file_char_long = file_str_long.c_str();
	const char *input_file_char_short = file_str_short.c_str();

	/*open file*/
	ifstream input_long;
	ifstream input_short;
	input_long.open(input_file_char_long);
	input_short.open(input_file_char_short);


	/*check existance of files*/
	if(!input_long)
	{
		cout << "Error: file could not be opened >> " + file_str_long << endl;
		system("pause");
		exit(1);
	}
	if(!input_short)
	{
		cout << "Error: file could not be opened >> " + file_str_short << endl;
		system("pause");
		exit(1);
	}


	/*This vector will serve as a container to hold the flux and kinetic energy from the flux code.*/
	vector<vector<double> > longtraj;
	vector<vector<double> > shorttraj;

	/*read in titles*/
	input_long >> t0 >> t1; /* This reads in the titles for each column for each quantity output
 	 	 	 	 	 	 	 from the rescattering code and sets it to the title "t" variables.*/
	input_short >> t0 >> t1;


	/*This takes calculated result from long trajectory file. The maximum and minimum kinetic energy from both long and short
	 * trajectories (not individually, but together) is recorded. */
	double E_min=1e150, E_max=0;
	while(input_long >> d0 >> d1)/* Read in long trajectories*/
	{
		if(d0<E_min){
			E_min = d0;
		}
		if(d0>E_max){
			E_max = d0;
		}
		longtraj.push_back({d0,d1});
	} /*end file read loop*/
	input_long.close();
	while(input_short >> d0 >> d1)/* Read in short trajectories*/
	{
		if(d0<E_min){
			E_min = d0;
		}
		if(d0>E_max){
			E_max = d0;
		}
		shorttraj.push_back({d0,d1});
	}/*end file read loop*/
	input_short.close();


	/**************************************************************************/
	// Begin linear binning algorithm
	/**************************************************************************/
	/* Define linear bins */
	vector<vector<double> > longbin;
	vector<vector<double> > shortbin;
	int i,s_counter,s_counter_pre=0,l_counter,l_counter_pre=0,k;
	double etmp1,etmp2,emean,shorttmp_sum,longtmp_sum; // These are temporary variables


	/* Check if linear steps or step size are specified */
	int LinSteps;
	if((LinStepSize_input > 0)&&(LinSteps_input > 0)){
		cout << "Error: Both inputs specifying step size are non-zero. Please select one of the two options provided " << endl;
		cout << "in the code while setting the other to zero." << endl;
		system("pause");
		exit(1);
	}
	else if(LinStepSize_input > 0){
		LinSteps = int((E_max - E_min)/LinStepSize_input);
	}
	else if(LinSteps_input > 0){
		LinSteps = LinSteps_input;
	}
	else{
		cout << "Error: Both inputs specifying step size are zero. Please select one of the two options provided " << endl;
		cout << "in the code while setting the other to zero." << endl;
		system("pause");
		exit(1);
	}
	/* End of check if linear steps or step size are specified */


	for(i=0; i<=LinSteps; i++){

		/* energy bounds */
		etmp1 = E_min + (double(i)-0.5)*LinStepSize_input; // left bound of energy bin
		etmp2 = E_min + (double(i)+0.5)*LinStepSize_input; // right bound of energy bin

		/* short steps */
		k=0; // this integer counts the energies in each bin
		shorttmp_sum = 0.0; // this takes the sum of the flux in each bin
		s_counter = s_counter_pre;
		if(s_counter < int(shorttraj.size())){ // Check to see if there are more energies to check from input file
			while(shorttraj[s_counter][0] < etmp2){ // Check to see if energy exceeds right bound of energy bin
				shorttmp_sum += shorttraj[s_counter][1];
				k++;
				s_counter++; // iterate to the next energy from input file
				if(s_counter >= int(shorttraj.size())){ // This conditional breaks the loop if the index exceeds the array size
					break;
				}
			}
			s_counter_pre = s_counter; // Save the last index to be checked next bin
			if(k!=0){
				emean = E_min + (double(i))*LinStepSize_input; // center of bin
				shortbin.push_back({emean,shorttmp_sum/double(k)}); // calculate the average flux within short bin
			}
			else{
				emean = E_min + (double(i))*LinStepSize_input; // center of bin
				shortbin.push_back({emean,0.0}); // give a flux of zero to avoid dividing by zero (with k)
			}
		}

		/* long steps */
		k=0; // this integer counts the energies in each bin
		longtmp_sum = 0.0; // this takes the sum of the flux in each bin
		l_counter = l_counter_pre;
		if(l_counter < int(longtraj.size())){ // Check to see if there are more energies to check from input file
			while(longtraj[l_counter][0] < etmp2){ // Check to see if energy exceeds right bound of energy bin
				longtmp_sum += longtraj[l_counter][1];
				k++;
				l_counter++; // iterate to the next energy from input file
				if(l_counter >= int(longtraj.size())){
					break;
				}
			}
			l_counter_pre = l_counter; // Save the last index to be checked next bin
			if(k!=0){
				emean = E_min + (double(i))*LinStepSize_input; // center of bin
				longbin.push_back({emean,longtmp_sum/double(k)}); // calculate the average flux within long bin
			}
			else{
				emean = E_min + (double(i))*LinStepSize_input; // center of bin
				longbin.push_back({emean,0.0}); // give a flux of zero to avoid dividing by zero (with k)
			}
		}
	}// End for loop


	/*generate bin output file */
	string binnedLong_out_str =  "binning_output/Z=" + Z + "_linear_binned_Long" + flag + "+" + ion + "_" + wavelength + "nm.dat";
	string binnedShort_out_str =  "binning_output/Z=" + Z + "_linear_binned_Short" + flag + "+" + ion + "_" + wavelength + "nm.dat";
	const char *binnedLong_output_file_char = binnedLong_out_str.c_str();
	const char *binnedShort_output_file_char = binnedShort_out_str.c_str();

	ofstream longfile, shortfile;
	longfile.open(binnedLong_output_file_char);
	shortfile.open(binnedShort_output_file_char);

	/* create title columns */
	longfile << "E" << " " << "dF_R/dE" << endl;
	shortfile << "E" << " " << "dF_R/dE" << endl;

	/* output binned calculation */
	for(i=0;i<int(shortbin.size());i++){ // This could also be longbin.size(). Post-binning they are the same size.
		longfile << longbin[i][0] << " " << longbin[i][1] << endl;
		shortfile << shortbin[i][0] << " " << shortbin[i][1] << endl;
	}
	longfile.close();
	shortfile.close();

	/*Print to command prompt that ion is compelete */
	cout << "Ion charge " + ion + " complete" << endl;
}




