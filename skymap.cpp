//program to generate a skymaps and axion-photon conversion probability
#include <iostream>
#include <vector>
#include <math.h>
#include <chrono>
#include <algorithm>
#include <omp.h>
#include <stdio.h>
#include "ALP_photon_conversion.h"
#include "B_field.h"
#include "matrix.h"
using namespace std;

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
    //! Prints a progress bar to the console
    //@param percentage: percentage of the progress
    //@return void
    //TO DO: change the color of the progress bar
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
int main(){
    //! Main function
    //@return 0
    //TO DO: add more parameters choices
	// Initialize the parameters
    long double basetime = 15.0, g, ma, mai, omega, d, Deltaa;
	int divmappa, energymeasure, single_multi, maxprobscale, minprobscale, scalesetting, type, Bringmul, omegaif;
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
		fscanf(file, "%s", filename);
		fclose(file);
		if(type == 11)
			customBfielevaluator(filename);
		if(ABSIF){
			gauleg(-1.0,1.0);
			dngamma_dE_evaluator();
			Gammaabs_value(omega);
		}
		if(multi){
			#pragma omp parallel shared(cpucount)
			{
				cpucount = omp_get_num_threads();
			}
			cout<<"Number of available threads: "<<cpucount<<endl;
			cout<<"How many of them would you like to use?"<<endl;
			cin >> multi;
			while(multi<=0 || multi>cpucount){
				cout<<"Expected integer value between 1 and "<<cpucount<<'\n'<<"Please insert a reasonable value: "<<endl;
				cin >> multi;
			}
		}
		cpucount = multi;
		cout<<"Skymap calculation in progress..." <<endl;
		long double li = -180.0;				//<------- initial latitude
		long double lf = 180.0;					//<------- final latitude
		long double bi = -90.0;					//<------- initial longitude
		long double bf = 90.0;					//<------- final longitude
		long double dl = (lf-li)/divmappa;		//<------- lattice step for latitude
		long double db = (bf-bi)/divmappa;		//<------- lattice step for longitude
		int divmappa1 = divmappa + 1;
		int N1 = pow(divmappa1,2), count = 0;
		vector <long double> R(N1,0);			//<------- initialize the map of the values of the skymap
		long double estimatime = (16.0/cpucount)*basetime*pow(divmappa/180.0,2)*N/200.0*(1.0+ma)*d/10.0;
		if(estimatime<=1)
			cout<<"Estimate time for the completion of the evaluation: "<<estimatime*60<<" s"<<endl;
		else
			cout<<"Estimate time for the completion of the evaluation: "<<estimatime<<" min"<<endl;
		cout<<"Please note that this time may vary depending on the machine on which this application is running."<<endl;
		// Skymap calculation
		auto ti=std::chrono::steady_clock::now();
		printProgress(0.0);
		#pragma omp parallel for num_threads(cpucount) schedule(dynamic)
		for(int i = 0; i<N1; i++){
			R[i] = Prob(li + (i/divmappa1)*dl, bi + (i%divmappa1)*db, d, g, omega, Deltaa, mai, ABSIF);
			count += 1;
			if((count+1)%(N1/100)==0)
				printProgress(float(i)/float(N1-1));
		}
		printProgress(1.0);
		// end of skymap calculation
		cout << "\n\n";
		auto tf=std::chrono::steady_clock::now();
		if (std::chrono::duration_cast <std::chrono::seconds>(tf-ti).count()>60.0)
			cout<<"Evaluation completed in "<<round(std::chrono::duration_cast <std::chrono::seconds>(tf-ti).count()/6.0)/10.0 <<" min" <<endl;
		else
			cout<<"Evaluation completed in "<<round(std::chrono::duration_cast <std::chrono::seconds>(tf-ti).count()*10.0)/10.0<<" s" <<endl;
		if((file=fopen("skymap.txt","w+"))==NULL)
			cout<<"Error opening the file"<<endl;
		int Digs = 16;
		fprintf(file, "%d\n", scalesetting);
		fprintf(file, "%f    %d\n", minprob, minprobscale);
		fprintf(file, "%f    %d\n", maxprob, maxprobscale);
		for(int i=0;i<N1;i++){
			fprintf(file,"%.*e    %.*e    %.*e\n", Digs ,  double(R[i]), Digs , double(li + (i/divmappa1)*dl), Digs , double(bi + (i%divmappa1)*db));
		}
		fclose(file);
		system("Mollweide_plot.exe");
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
		fclose(file);
		if(ABSIF){
			gauleg(-1.0,1.0);
			dngamma_dE_evaluator();
			Gammaabs_value(omega);
		}
		cout<<"Calculation in progress..." <<endl;
		if(single_multi==1){
			// Calculation of the probability one a single line
			ProbSingleLine(lx, bx, d, g, omega, Deltaa, mai, ABSIF);
			// End of the calculation of the probability one a single line
			system("1D_plot.exe");
		}
		else{
			int N1 = 1000;
			vector <long double> xR(N1,0);
			vector <long double> R(N1,0);
			if(single_multi==2){
				long double ma_i = 0;
				long double ma_f = 10;
                if(multi){
                    #pragma omp parallel shared(cpucount)
                    {
                        cpucount = omp_get_num_threads();
                    }
                    cout<<"Number of available threads: "<<cpucount<<endl;
                    cout<<"How many of them would you like to use?"<<endl;
                    cin >> multi;
                    while(multi<=0 || multi>cpucount){
                        cout<<"Expected integer value between 1 and "<<cpucount<<'\n'<<"Please insert a reasonable value: "<<endl;
                        cin >> multi;
                    }
                }
                cpucount = multi;
				// Calculation of the probability with respect to the mass
				printProgress(0.0);
				#pragma omp parallel for num_threads(cpucount) schedule(dynamic)
				for(int i = 0; i<N1; i++){
					ma = ma_i + i*(ma_f-ma_i)/N1;
					mai = ma + 1;
					xR[i] = ma;
					R[i] = Prob(lx, bx, d, g, omega, Deltaa, mai, ABSIF);
					if((i+1)%(N1/100)==0)
						printProgress(float(i)/float(N1-1));
				}
				printProgress(1.0);
				// End of the calculation of the probability with respect to the mass
				cout << "\n";
				if((file=fopen("1D_plot.txt","w+"))==NULL)
					cout<<"Error opening the file"<<endl;
				double x_component, probability;
				int Digs = 16;
				fprintf(file, "%d\n",  single_multi);
				for(int i=0;i<N1;i++){
					x_component=double(xR[i]);
					probability=double(R[i]);
					fprintf(file,"%.*e    %.*e\n", Digs ,  x_component, Digs ,  probability);
				}
				fclose(file);
				cout<<"Your plot is saved in: 1D_plot.txt and in 1D_plot.pdf in this folder."<<endl;
				system("1D_plot.exe");
			}
			else if(single_multi==3){
				long double omegai = 1*pow(10,3*(energymeasure-2));
				long double omegaf = 1000*pow(10,3*(energymeasure-2));
                if(multi){
                    #pragma omp parallel shared(cpucount)
                    {
                        cpucount = omp_get_num_threads();
                    }
                    cout<<"Number of available threads: "<<cpucount<<endl;
                    cout<<"How many of them would you like to use?"<<endl;
                    cin >> multi;
                    while(multi<=0 || multi>cpucount){
                        cout<<"Expected integer value between 1 and "<<cpucount<<'\n'<<"Please insert a reasonable value: "<<endl;
                        cin >> multi;
                    }
                }
                cpucount = multi;
				// Calculation of the probability with respect to the energy
				printProgress(0.0);
				#pragma omp parallel for num_threads(cpucount) schedule(dynamic)
				for(int i = 0; i<N1; i++){
					omega = omegai + i*(omegaf-omegai)/N1;
					xR[i] = omega;
					R[i] = Prob(lx, bx, d, g, omega, Deltaa, mai, ABSIF);
					if((i+1)%(N1/100)==0)
						printProgress(float(i)/float(N1-1));
				}
				printProgress(1.0);
				// End of the calculation of the probability with respect to the energy
				cout<<"\n";
				if((file=fopen("1D_plot.txt","w+"))==NULL)
					cout<<"Error opening the file"<<endl;
				double x_component, probability;
				int Digs = 16;
				fprintf(file, "%d     %d\n",  single_multi, energymeasure);
				for(int i=0;i<N1;i++){
					x_component=double(xR[i]);
					probability=double(R[i]);
					fprintf(file,"%.*e    %.*e\n", Digs ,  x_component, Digs ,  probability);
				}
				fclose(file);
				cout<<"Your plot is saved in: 1D_plot.txt and in 1D_plot.pdf in this folder."<<endl;
				system("1D_plot.exe");
			}
			else if(single_multi==4){
				long double gai = g/10.0;
				long double gaf = g*10.0;
                if(multi){
                    #pragma omp parallel shared(cpucount)
                    {
                        cpucount = omp_get_num_threads();
                    }
                    cout<<"Number of available threads: "<<cpucount<<endl;
                    cout<<"How many of them would you like to use?"<<endl;
                    cin >> multi;
                    while(multi<=0 || multi>cpucount){
                        cout<<"Expected integer value between 1 and "<<cpucount<<'\n'<<"Please insert a reasonable value: "<<endl;
                        cin >> multi;
                    }
                }
		        cpucount = multi;
				// Calculation of the probability with respect to the g
				printProgress(0.0);
                #pragma omp parallel for num_threads(cpucount) schedule(dynamic)
				for(int i = 0; i<N1; i++){
					g = gai + i*(gaf-gai)/N1;
					xR[i] = g;
					R[i] = Prob(lx, bx, d, g, omega, Deltaa, mai, ABSIF);
					if((i+1)%(N1/100)==0)
						printProgress(float(i)/float(N1-1));
				}
				printProgress(1.0);
				// End of the calculation of the probability with respect to the g
				cout << "\n";
				if((file=fopen("1D_plot.txt","w+"))==NULL)
					cout<<"Error opening the file"<<endl;
				double x_component, probability;
				int Digs = 16;
				fprintf(file, "%d     %d\n",  single_multi, energymeasure);
				for(int i=0;i<N1;i++){
					x_component=double(xR[i]);
					probability=double(R[i]);
					fprintf(file,"%.*e    %.*e\n", Digs ,  x_component, Digs ,  probability);
				}
				fclose(file);
				cout<<"Your plot is saved in: 1D_plot.txt and in 1D_plot.pdf in this folder."<<endl;
				system("1D_plot.exe");
			}
			else if(single_multi==5){
				gauleg(-1.0,1.0);
				// Calculation of the probability with non perturbative theory
				double pb = (double)Non_pert_Prob(lx,bx,d, g, omega, ma);
				cout<<pb<<endl;
				// End of the calculation of the probability with non perturbative theory
				int Digs = 16;
				if((file=fopen("1D_plot.txt","w+"))==NULL)
					cout<<"Error opening the file"<<endl;
				fprintf(file,"%.*e\n", Digs ,  pb);
				fclose(file);
				cout<<"Your plot is saved in: 1D_plot.txt and in 1D_plot.pdf in this folder."<<endl;
			}
            else{
                cout<<"Wrong input"<<endl;
            }
		}
	}
	return 0;
}
