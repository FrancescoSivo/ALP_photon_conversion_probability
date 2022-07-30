//program to generate a skymaps and axion-photon conversion probability
#include <iostream>
#include <vector>
#include <math.h>
#include <chrono>
#include <algorithm>
#include <stdio.h>
#include <cstring>
#include "ALP_photon_conversion.h"
#include "B_field.h"
#include "matrix.h"
using namespace std;

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

int main(){
    //! Main function
    //@return 0
    //TO DO: add more parameters choices
	// Initialize the parameters
    long double basetime = 15.0, g, ma, mai, omega, d, Deltaa, ma_f, ga_f, omega_f;
	int divmappa, energymeasure, single_multi, maxprobscale, minprobscale, scalesetting, type, Bringmul, omegaif, N_gauleg, energymeasure_f, sampling_type, N1, count = 0, plot_if, adornot;
	unsigned int multi, cpucount;
	float input;
    string filename;
    FILE* file;
	// Taking inputs from file
	if((file=fopen("ic.txt","r"))==NULL)
			cout<<"Error opening the file"<<endl;
	fscanf(file, "%d", &single_multi);
	if(single_multi == 0){
		fscanf(file, "%f", &input);
		g = (long double) input;
		g *= pow(10,-11);
		fscanf(file, "%f", &input);
		ma = (long double) input;
		mai = ma + 1.0;
		fscanf(file, "%f", &input);
		omega = (long double) input;
		fscanf(file, "%d", &energymeasure);
		omega*=pow(10,3*(energymeasure-2));
		fscanf(file, "%f", &input);
		d = (long double) input;
		fscanf(file, "%d", &divmappa);
		fscanf(file, "%d", &type);
		Bset(type);
		fscanf(file, "%d", &Bringmul);
		if(type>2)
			Bringmul = 0;
		fscanf(file, "%d", &multi);
		fscanf(file, "%d", &N);
		fscanf(file, "%f", &input);
		long double obsx = (long double) input;
		fscanf(file, "%f", &input);
		long double obsy = (long double) input;
		fscanf(file, "%f", &input);
		long double obsz = (long double) input;
		obs = {obsx,obsy,obsz};
		fscanf(file, "%d", &omegaif);
        omegaifset(omegaif);
		fscanf(file, "%d", &ABSIF);
		Deltaa=-c*0.5*pow(10,-15)*pow(ma,2.0)/omega;
		fscanf(file, "%d", &scalesetting);
		fscanf(file, "%f", &input);
		float minprob = input;
		fscanf(file, "%d", &minprobscale);
		fscanf(file, "%f", &input);
		float maxprob =  input;
		fscanf(file, "%d", &maxprobscale);
		fscanf(file, "%d", &cpucount);
		fscanf(file, "%d", &N_gauleg);
		//scan from file the string filename
		fscanf(file, "%s", filename.c_str());
		fscanf(file, "%d", &plot_if);
		fscanf(file, "%d", &adornot);
		fclose(file);
		if(type == 11)
			customBfielevaluator(filename);
		if(ABSIF){
			if(!gauleg(-1.0,1.0,N_gauleg))
				return 0;
			dngamma_dE_evaluator();
			Gammaabs_value(omega, N_gauleg);
		}
		cout<<"Skymap calculation in progress..." <<endl;
		long double li = -180.0;				//<------- initial latitude
		long double lf = 180.0;					//<------- final latitude
		long double bi = -90.0;					//<------- initial longitude
		long double bf = 90.0;					//<------- final longitude
		long double dl = (lf-li)/divmappa;		//<------- lattice step for latitude
		long double db = (bf-bi)/divmappa;		//<------- lattice step for longitude
		int divmappa1 = divmappa + 1;
		N1 = pow(divmappa1,2);
		// Skymap calculation
		for(int i = 0; i<N1; i++){
			if(adornot)
				cout << li + (i/divmappa1)*dl<<"   "<<bi + (i%divmappa1)*db<<"   "<<Prob(li + (i/divmappa1)*dl, bi + (i%divmappa1)*db, d, g, omega, Deltaa, mai, ABSIF)<<endl;
			else
				cout << li + (i/divmappa1)*dl<<"   "<<bi + (i%divmappa1)*db<<"   "<<Prob_nonadaptive(li + (i/divmappa1)*dl, bi + (i%divmappa1)*db, d, g, omega, Deltaa, mai, ABSIF)<<endl;
		}
	}
	else{
		fscanf(file, "%f", &input);
		g = (long double) input;
		g *= pow(10,-11);
		fscanf(file, "%f", &input);
		ma = (long double) input;
		mai = ma + 1.0;
		fscanf(file, "%f", &input);
		omega = (long double) input;
		fscanf(file, "%d", &energymeasure);
		omega*=pow(10,3*(energymeasure-2));
		fscanf(file, "%f", &input);
		d = (long double) input;
		int lxtype, bxtype;
		fscanf(file, "%f", &input);
		long double lx = (long double) input;
		fscanf(file, "%d", &lxtype);
		if(lxtype)
			lx *= 180.0/PI;
		fscanf(file, "%f", &input);
		long double bx = (long double) input;
		fscanf(file, "%d", &bxtype);
		if(bxtype)
			bx *= 180.0/PI;
		fscanf(file, "%d", &type);
		Bset(type);
		fscanf(file, "%d", &Bringmul);
		if(type>2)
			Bringmul = 0;
		fscanf(file, "%d", &multi);
		fscanf(file, "%d", &N);
		fscanf(file, "%f", &input);
		long double obsx = (long double) input;
		fscanf(file, "%f", &input);
		long double obsy = (long double) input;
		fscanf(file, "%f", &input);
		long double obsz = (long double) input;
		obs = {obsx,obsy,obsz};
		fscanf(file, "%d", &omegaif);
        omegaifset(omegaif);
		fscanf(file, "%d", &ABSIF);
		Deltaa=-c*0.5*pow(10,-15)*pow(ma,2.0)/omega;
		fscanf(file, "%d", &BCOMPARISON);
		fscanf(file, "%d", &N_gauleg);
		fscanf(file, "%f", &input);
		ga_f = (long double) input;
		ga_f *= pow(10,-11);
		fscanf(file, "%f", &input);
		ma_f = (long double) input;
		fscanf(file, "%f", &input);
		omega_f = (long double) input;
		fscanf(file, "%d", &energymeasure_f);
		omega_f*=pow(10,3*(energymeasure_f-2));
		fscanf(file, "%d", &sampling_type);
		fscanf(file, "%d", &N1);
		fscanf(file, "%d", &plot_if);
		fscanf(file, "%d", &adornot);
		fclose(file);
		if(ABSIF){
			if(!gauleg(-1.0,1.0,N_gauleg))
				return 0;
			dngamma_dE_evaluator();
			Gammaabs_value(omega, N_gauleg);
		}
		if(single_multi==1){
			// Calculation of the probability one a single line
			if(adornot)
				ProbSingleLine(lx, bx, d, g, omega, Deltaa, mai, ABSIF);
			else
				ProbSingleLine_nonadaptive(lx, bx, d, g, omega, Deltaa, mai, ABSIF);
		}
		else{
			if(single_multi==2){
				long double ma_i = ma;
				// Calculation of the probability with respect to the mass
				long double deltaa = -c*0.5*pow(10,-15)/omega;
				cout<<single_multi<<endl;
				for(int i = 0; i<=N1; i++){
					if(sampling_type==0)
						ma = ma_i*pow(10,float(i)/float(N1) * log10(ma_f/ma_i));
					else
						ma = ma_i + float(i)/float(N1) * (ma_f-ma_i);
					Deltaa=deltaa*pow(ma,2.0);
					if(adornot)
						cout<<ma<<"   "<<Prob(lx, bx, d, g, omega, Deltaa, ma + 1, ABSIF)<<endl;
					else
						cout<<ma<<"   "<<Prob_nonadaptive(lx, bx, d, g, omega, Deltaa, ma + 1, ABSIF)<<endl;
				}
			}
			else if(single_multi==3){
				long double omega_i = omega;
				// Calculation of the probability with respect to the energy
				long double deltaa = -c*0.5*pow(10,-15)*pow(ma,2.0);
				cout<<single_multi<<"   "<<energymeasure<<endl;
				for(int i = 0; i<=N1; i++){
					if(sampling_type==0)
						omega = omega_i*pow(10,float(i)/float(N1) * log10(omega_f/omega_i));
					else
						omega = omega_i + float(i)/float(N1) * (omega_f-omega_i);
					Deltaa = deltaa/omega;
					if(adornot)
						cout<<omega<<"   "<<Prob(lx, bx, d, g, omega, Deltaa, mai, ABSIF)<<endl;
					else
						cout<<omega<<"   "<<Prob_nonadaptive(lx, bx, d, g, omega, Deltaa, mai, ABSIF)<<endl;
				}
			}
			else if(single_multi==4){
				long double ga_i = g;
				// Calculation of the probability with respect to the g
				cout<<single_multi<<"   "<<energymeasure<<endl;
				for(int i = 0; i<=N1; i++){
					if(sampling_type==0)
						g = ga_i*pow(10,float(i)/float(N1) * log10(ga_f/ga_i));
					else
						g = ga_i + float(i)/float(N1) * (ga_f-ga_i);
					if(adornot)
						cout<<g<<"    "<<Prob(lx, bx, d, g, omega, Deltaa, mai, ABSIF)<<endl;
					else
						cout<<g<<"    "<<Prob_nonadaptive(lx, bx, d, g, omega, Deltaa, mai, ABSIF)<<endl;
				}
			}
			else if(single_multi==5){
				gauleg(-1.0,1.0,N_gauleg);
				// Calculation of the probability with non perturbative theory
				double pb = (double)Non_pert_Prob(lx,bx,d, g, omega, ma, N_gauleg);
				cout<<pb<<endl;
			}
            else{
                cout<<"Wrong input"<<endl;
            }
		}
	}
	return 0;
}
