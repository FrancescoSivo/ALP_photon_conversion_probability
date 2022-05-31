#ifndef ALP_PHOTON_CONVERSION_H
#define ALP_PHOTON_CONVERSION_H
#include <iostream>
#include <vector>
#include <math.h>
#include "B_field.h"
#include "matrix.h"

//function that initializes a complex matrix with all elements equal to 0.0 + 0.0i but the [2][2] element equal to 1.0 + 0.0i
mat_complex ALPs_initial_conditions(){
    //! initializes a complex matrix with all elements equal to 0.0 + 0.0i but the [2][2] element equal to 1.0 + 0.0i
    //@return: the initialized matrix
    mat_complex mat(3);
    for(int i = 0; i < 3; i++){
        mat[i].resize(3);
        for(int j = 0; j < 3; j++){
            if(i == 2 && j == 2){
                mat[i][j] = complex<long double>(1.0, 0.0);
            }
            else{
                mat[i][j] = complex<long double>(0.0, 0.0);
            }
        }
    }
    return mat;
}
//function that impose the required conditions on the density matrix
//function that returns a complex matrix with simmetric real elements an antisimmetric imaginary elements
mat_complex ALPs_density_matrix_condition(mat_complex mat){
	//! imposes the required conditions on the density matrix
	//! function that returns a complex matrix with simmetric real elements an antisimmetric imaginary elements
	//@param mat: the density matrix
	//@return: the density matrix with the required conditions
    mat[1][0] = complex<long double>(mat[0][1].real(), -mat[0][1].imag());
	mat[2][0] = complex<long double>(mat[0][2].real(), -mat[0][2].imag());
	mat[2][1] = complex<long double>(mat[1][2].real(), -mat[1][2].imag());
    return mat;
}

//Global variables
long double Gammaabs;										//< ----   Initial inputs (coupling, axion mass, skymap distance, energy) + other stuff
long double c=1.56*pow(10,17);								//< ----   Multiplicative constant of the model
vec_real obs = {-8.5,0.0,0.0};								//< ----   Observer coordinates ({-8.5,0.0,0.0} for the Earth)
long double conv=PI/180.0;									//< ----   Degree to radian conversion factor
int N;		    											//< ----   Average number of iterations per kpc for runge-kutta
long double BPRINT = 0.0;

//Runge Kutta coefficients
//coefficients of order 5
long double AA_RK8[8][7]={{0.0,0.0,0.0,0.0,0.0,0.0,0.0},{1.0/6.0,0.0,0.0,0.0,0.0,0.0,0.0},{4.0/75.0,16.0/75.0,0.0,0.0,0.0,0.0,0.0},{5.0/6.0,-8.0/3.0,5.0/2.0,0.0,0.0,0.0,0.0},{-8.0/5.0,144.0/25.0,-4.0,16.0/25.0,0.0,0.0,0.0},{361.0/320.0,-18.0/5.0,407.0/128.0,-11.0/80.0,55.0/128.0,0.0,0.0},{-11.0/640.0,0.0,11.0/256.0,-11.0/160.0,11.0/256.0,0.0,0.0},{93.0/640.0,-18.0/5.0,803.0/256.0,-11.0/160.0,99.0/256.0,0.0,1.0}};
vec_real BB_RK8={31.0/384.0,0.0,1125.0/2816.0,9.0/32.0,125.0/768.0,5.0/66.0,0.0,0.0};
vec_real BBSTAR_RK8={7.0/1408.0,0.0,1125.0/2816.0,9.0/32.0,125.0/768.0,0.0,5.0/66.0,5.0/66.0};
vec_real CC_RK8={0.0,1.0/6.0,4.0/15.0,2.0/3.0,4.0/5.0,1.0,0.0,1.0};
int order = 8;
int errororder = order - 4;
int N_GAULEG = 500;
vector<long double> XGAUSS(N_GAULEG), WGAUSS(N_GAULEG);
void gauleg(long double x1, long double x2){
    //! This function calculates the Gauss-Legendre quadrature points and weights
    //! for the integration. It does not return anything, but it sets the global
	//! variables XGAUSS and WGAUSS which are respectively the quadrature points
	//! and the weights.
	//! This function needs to be called only once, when during the process the user
	//! wants to use absorpion or non-perturbative calculations of the probability.
    //@param x1: lower integration limit (set to -1.0)
    //@param x2: upper integration limit (set to 1.0)
    //@return: none
    //TO DO: let the user choose the number of points
	long double epsilon_gauleg = 3.0*pow(10.0,-14.0);
	for(int i = 0; i < (N_GAULEG)/2; i++){
		XGAUSS[i] = 0.0;
		WGAUSS[i] = 0.0;
	}
	int mm = (N_GAULEG + 1)/2;
	long double xm = 0.5*(x2 + x1);
	long double xl = 0.5*(x2 - x1);
	long double z, p1, p2, p3, pp, z1;
	vector<long double> xx(N_GAULEG+1), ww(N_GAULEG+1);
	for(int i = 0; i < mm; i++){
		z = cos(PI*((i+1)*1.0-0.25)/(N_GAULEG*1.0+0.5));
		z1 = z + 1.0;
		while((abs(z-z1)>epsilon_gauleg)){
			p1 = 1.0;
			p2 = 0.0;
			for(int j = 1; j < N_GAULEG+1; j++){
				p3 = p2;
				p2 = p1;
				p1 = ((2.0*j - 1.0)*z*p2 - (j-1.0)*p3)/j;
			}
			pp = N_GAULEG*(z*p1 - p2)/(pow(z,2.0)-1.0);
			z1 = z;
			z = z1 - p1/pp;
		}
		xx[i] = xm - xl*z;
		xx[N_GAULEG - 1 - i] = xm + xl*z;
		ww[i] = 2.0*xl/((1.0 - pow(z,2.0))*pow(pp,2.0));
		ww[N_GAULEG - 1 - i] = ww[i];

	}
	for(int i = 0; i < (N_GAULEG)/2+1; i++){
		XGAUSS[N_GAULEG - 2*i] = xx[i - 1];
		XGAUSS[N_GAULEG - 2*(i - 1) - 1] = xx[N_GAULEG - i];
		WGAUSS[N_GAULEG - 2*i] = ww[i - 1];
		WGAUSS[N_GAULEG - 2*(i - 1) - 1] = ww[N_GAULEG - i];
	}
}


//Function definitions of physical objects
//Definitions for the absorption
const int N_dngamma = 1000000;
long double dngamma_dE_tab[N_dngamma][2],DeltaE_interpol, m_e = 0.511, sigma_0 = 1.3203125*pow(10.0,-74.0);
void dngamma_dE_evaluator(){
    //! initialize the table of dngamma/dE
    //@param: none
    //@return: none    
	FILE* file;
	float x,y;
	if((file=fopen("dngamma_dE.dat","r"))==NULL)
			cout<<"Error opening the file"<<endl;
	for(int i = 0; i < N_dngamma; i++){
		fscanf(file, "%f", &x);
		fscanf(file, "%f", &y);
		dngamma_dE_tab[i][0] = (long double)x;
		dngamma_dE_tab[i][1] = (long double)y;
	}
	fclose(file);
	DeltaE_interpol = (dngamma_dE_tab[N_dngamma-1][0] - dngamma_dE_tab[0][0])/N_dngamma;
}
long double dngamma_dE(long double e){
    //! calculate the value of dngamma/dE at energy e with a linear interpolation
    //@param: e (MeV), energy of the photon
    //@return: dngamma/dE
	if(e < dngamma_dE_tab[0][0])
		return 0;
	else{
		int i = floor((e - dngamma_dE_tab[0][0])/DeltaE_interpol);
		return dngamma_dE_tab[i][1] + (dngamma_dE_tab[i+1][1]-dngamma_dE_tab[i][1])/(DeltaE_interpol)*(e - dngamma_dE_tab[i][0]);
	}
}
void Gammaabs_value(long double omega){
    //! calculate the value of Gammaabs
    //@param: omega (MeV), energy of the ALP
    //@return: none
	long double diff1,somm1,diff2,somm2, xmax1, xmin1, xmax2, xmin2, res = 0.0, res_i = 0.0, eps_i, eps_j, beta_j, sigma_gg;
	long double beta_j2;
	xmin1 = pow(m_e,2.0)/omega;
	xmax1 = dngamma_dE_tab[N_dngamma - 1][0];
	diff1 = (xmax1 - xmin1)/2.0;
	somm1 = (xmax1 + xmin1)/2.0;
	for(int i = 0; i < N_GAULEG; i++){
		eps_i = diff1*XGAUSS[i] + somm1;
		res_i = 0.0;
		xmax2 = 1.0 - 2.0*pow(m_e,2.0)/(omega * eps_i);
		xmin2 = -1.0;
		diff2 = (xmax2 - xmin2)/2.0;
		somm2 = (xmax2 + xmin2)/2.0;
		for(int j = 0; j < N_GAULEG; j++){
			eps_j = diff2*XGAUSS[j] + somm2;
			beta_j = sqrt(1.0 - 2.0*pow(m_e,2.0)/(omega*eps_i*(1.0 - eps_j)));
			beta_j2 = pow(beta_j,2.0);
			sigma_gg = sigma_0 * (1.0 - beta_j2)*(2.0 * beta_j*(beta_j2 - 2.0) + (3.0 - pow(beta_j,4.0))*log((1+beta_j)/(1-beta_j)));
			res_i += sigma_gg*(1.0 - eps_j)/2.0*WGAUSS[j];
		}
		res_i *= diff2;
		res += res_i*dngamma_dE(eps_i)*WGAUSS[i];
	}
	Gammaabs = -res*diff1*0.5;
}

mat_complex Mk(long double z , vec_real kdir, long double g, long double omega, long double Deltaa, bool absif){ 	
    //! calculate the matrix Mk for the Von Neumann commutation equation
    //@param: z (kpc), the z-component of the ALP
    //@param: kdir (rad,rad,rad), the direction of the ALP
    //@param: g (GeV^-1), the coupling constant of the ALP with the photon
    //@param: omega (MeV), the energy of the photon
    //@param: Deltaa (MeV^-1), the Delta_a related to the mass of the ALP
	//@param: absif (bool), if true, the absorption is considered in the calculation
    //@return: Mk
    long double xi = obs[0] + z*kdir[0];
	long double yi = obs[1] + z*kdir[1];
	long double zi = obs[2] + z*kdir[2];
	long double mgz2 = omegapl2(xi,yi,zi);
	vec_real Bvecz = Bfield(xi,yi,zi);
	if(BCOMPARISON == 1)
		BPRINT = sqrt(pow(Bvecz[0],2.0)+pow(Bvecz[1],2.0)+pow(Bvecz[2],2.0));
	long double Bpar = Bvecz[0]*kdir[0] + Bvecz[1]*kdir[1] + Bvecz[2]*kdir[2];
	vec_real BTz = {Bvecz[0]-Bpar*kdir[0],Bvecz[1]-Bpar*kdir[1],Bvecz[2]-Bpar*kdir[2]};
	long double Bmodz = sqrt(pow(Bvecz[0],2.0)+pow(Bvecz[1],2.0)+pow(Bvecz[2],2.0));
	long double BTmodz = sqrt(pow(BTz[0],2.0)+pow(BTz[1],2.0)+pow(BTz[2],2.0));
	long double xiz = 2.64114 * pow(10,-44)*pow(Bmodz,2.0);
	vec_real ort = {BTz[0],0.0,BTz[2]};
	long double normort = sqrt(pow(ort[0],2.0)+pow(ort[2],2.0));
	long double sps = 1.0;
	long double cps = 1.0;
	if(BCOMPARISON == 2)
		BPRINT = BTmodz;
	if (normort!=0 && BTmodz!=0){
		if (BTz[0]>0.0)
			sps = (BTz[0]*ort[0] + BTz[1]*ort[1] + BTz[2]*ort[2])/(BTmodz*normort);
		else
			sps = -(BTz[0]*ort[0] + BTz[1]*ort[1] + BTz[2]*ort[2])/(BTmodz*normort);
	}
	if (BTmodz!=0.0)
		cps = BTz[1]/BTmodz;
	long double sps2 = pow(sps,2);
	long double cps2 = pow(cps,2);
	long double D1=0.0;
	if(Bmodz!=0.0)
		D1 = c*pow(10,15)*xiz*omega*(1.0-pow(Bpar/Bmodz,2.0));
	long double D2 = -c*0.5*pow(10,-15)*mgz2/omega + 8.0*pow(10,-11)*omega;
	long double Dpar = 3.5*D1 + D2;
	long double Dperp = 2.0*D1 + D2;
	long double Dag = c*g*BTmodz*0.975*pow(10,-8);
	long double Dxy = (Dpar - Dperp)*sps*cps;
	mat_complex result = init_mat_complex(3,3);
	result[0][0] = Dpar*sps2 + Dperp*cps2;
    if(absif){
        result[0][0] += Gammaabs;
	}
	result[0][1] = Dxy;
	result[0][2] = Dag*sps;
	result[1][0] = Dxy;
	result[1][1] = Dperp*sps2 + Dpar*cps2;
    if(absif){
        result[1][1] += Gammaabs;
	}
	result[1][2] = Dag*cps;
	result[2][0] = Dag*sps;
	result[2][1] = Dag*cps;
	result[2][2]= Deltaa;
	return result;
}

long double Prob(long double ldir, long double bdir, long double distz, long double g, long double omega, long double Deltaa, long double mai, bool absif){
	//! calculate the probability of the photon to be absorbed by the ALP
	//@param: ldir (deg), the l direction of the ALP from the observer
	//@param: bdir (deg), the b direction of the ALP from the observer
	//@param: distz (kpc), the initial distance between ALP and the observer
	//@param: g (GeV^-1), the coupling constant of the ALP with the photon
	//@param: omega (MeV), the energy of the ALP
	//@param: Deltaa (MeV^-1), the Delta_a related to the mass of the ALP
	//@param: mai (neV), the mass of the ALP incremented by 1 neV
	//@param: absif (bool), if true, the absorption is considered in the calculation
	//@return: Prob
	ldir*=conv;
	bdir*=conv;
	long double Delta0, zi, dz0, dz, Delta1, Delta18;
	//I apply the initial conditions to the density matrix
	mat_complex rho;
	mat_complex Erho;
	mat_complex rhoidzak;
	mat_complex Mkzi;
	mat_complex ak;
	mat_complex k_RK[8];
	long double cosbdir = cos(bdir);
	vec_real kdir = {cos(ldir)*cosbdir , sin(ldir)*cosbdir , sin(bdir)};		//Direction vector definition
	long double div = N*distz;													//Number of average divisions desired as a function of the distance to be covered
	long double p = 0.0;
	int i = 1;
	while ((p>0.1 && i<=3)||(p>0.5 && i<=6 && i>3 )|| p<=0 || isnan(p)){
		rho = ALPs_initial_conditions();
		Erho = ALPs_initial_conditions();
		zi = distz;
		dz0 = -distz/(div*i*mai);      //Average step for the Runge Kutta
		dz = dz0;
		Delta1 = pow(abs(dz0),errororder);    //Estimated desired error as a function of the order of the RK method and the base average step
		Delta18 = Delta1/8.0;
		while (zi>0){				//I solve the Von Neuman commutator equation with Runge-Kutta method
			//initialize all k_RK with init_mat_complex
            for (int j=0; j<8; j++){
                k_RK[j] = init_mat_complex(3,3);
            }
			for(int j=0; j<order;j++){
				Mkzi = Mk(zi + CC_RK8[j]*dz, kdir, g, omega, Deltaa, absif);
				ak = init_mat_complex(3,3);
				for(int k = 0; k < j; k++){
                        ak = sum_mat(ak , mult_mat_scalar(k_RK[k],AA_RK8[j][k]));
                }
                rhoidzak = sum_mat(rho,mult_mat_scalar(ak,dz));
                k_RK[j] = mult_mat_scalar(commutator(Mkzi,rhoidzak),(0.0,-dz));  //!TO DO: check the warning
				}
            for(int j=0; j<order; j++){
                rho = sum_mat(rho,mult_mat_scalar(k_RK[j],BB_RK8[j]));
                Erho = sum_mat(Erho,mult_mat_scalar(k_RK[j],BB_RK8[j]));
            }
			rho = ALPs_density_matrix_condition(rho);
			Erho = ALPs_density_matrix_condition(Erho);
			Delta0=fmax(max_value(matrix_abs(sub_mat(rho,Erho))),0.0);
			zi+=dz;
			if (Delta0<Delta18)
				dz = dz0*1.8;
			else
				dz = 0.9*dz0*fmin(fmax(0.3 , pow(Delta1/(2.0*Delta0) , 0.5 )),2.0);
		}
		p = 1.0 - rho[2][2].real();
		i+=1;
	}
	return p;
}

long double Non_pert_Prob(long double ldir, long double bdir, long double distz, long double g, long double omega, long double ma){
	//! calculate the probability ALP photon conversion with a non-perturbative calculation
	//@param: ldir (deg), the l direction of the ALP from the observer
	//@param: bdir (deg), the b direction of the ALP from the observer
	//@param: distz (kpc), the initial distance between ALP and the observer
	//@param: g (GeV^-1), the coupling constant of the ALP with the photon
	//@param: omega (MeV), the energy of the ALP
	//@param: ma (neV), the mass of the ALP
	//@return: Prob
	ldir*=conv;
	bdir*=conv;
	long double cosbdir = cos(bdir);
	vec_real kdir = {cos(ldir)*cosbdir , sin(ldir)*cosbdir , sin(bdir)};
	long double eta = 156.0*pow(ma,2.0)/(2.0*omega);
	long double are = 0.0;
	long double aim = 0.0;
	long double are2 = 0.0;
	long double aim2 = 0.0;
	long double diff = distz/2.0;
	long double somm = diff;
	long double rr,xi,yi,zi,Bdotk,BTmod,stheta,ctheta,sgn, ortmod,Bdoty,con,P1,P2;
	vec_real Bvecz= {0.0,0.0,0.0}, BT = {0.0,0.0,0.0} , ort= {0.0,0.0,0.0};
	for(int i = 0; i < N_GAULEG; i++){
		rr = diff*XGAUSS[i]+somm;
		xi = obs[0]+rr*kdir[0];
		yi = obs[1]+rr*kdir[1];
		zi = obs[2]+rr*kdir[2];
		Bvecz = Bfield(xi,yi,zi);
		Bdotk = Bvecz[0]*kdir[0] + Bvecz[1]*kdir[1] + Bvecz[2]*kdir[2];
		BT[0] = Bvecz[0] - Bdotk*kdir[0];
		BT[1] = Bvecz[1] - Bdotk*kdir[1];
		BT[2] = Bvecz[2] - Bdotk*kdir[2];
		BTmod = sqrt(pow(BT[0],2.0)+pow(BT[1],2.0)+pow(BT[2],2.0));
		if(BTmod < pow(10,-10))
			ctheta = 0.0;
		else
			ctheta = BT[1]/BTmod;
		Bdoty = BT[1];
		ort[0] = BT[0];
		ort[1] = BT[1]-Bdoty;
		ort[2] = BT[2];
		ortmod = sqrt(pow(ort[0],2.0)+pow(ort[1],2.0)+pow(ort[2],2.0));
		if(ortmod < pow(10,-10))
			stheta = 0.0;
		else{
			sgn = BT[0]/abs(BT[0]);
			stheta = sgn*(BT[0]*ort[0] + BT[1]*ort[1] + BT[2]*ort[2])/(BTmod*ortmod);
		}
		are += cos(eta*(rr-distz))*BTmod*ctheta*WGAUSS[i];
		aim -= sin(eta*(rr-distz))*BTmod*ctheta*WGAUSS[i];
		are2 += cos(eta*(rr-distz))*BTmod*stheta*WGAUSS[i];
		aim2 -= sin(eta*(rr-distz))*BTmod*stheta*WGAUSS[i];
	}
	con = 0.5 * g * pow(10.0,10.0) * diff * 0.304504;
	are *= con;
	aim *= con;
	are2 *= con;
	aim2 *= con;
	P1 = pow(are,2.0) + pow(aim,2.0);
	P2 = pow(are2,2.0) + pow(aim2,2.0);
	return P1 + P2;
}

long double ProbSingleLine(long double ldir, long double bdir, long double distz, long double g, long double omega, long double Deltaa, long double mai, bool absif){
	//! calculate the probability of the photon to be absorbed by the ALP
	//@param: ldir (deg), the l direction of the ALP from the observer
	//@param: bdir (deg), the b direction of the ALP from the observer
	//@param: distz (kpc), the initial distance between ALP and the observer
	//@param: g (GeV^-1), the coupling constant of the ALP with the photon
	//@param: omega (MeV), the energy of the ALP
	//@param: Deltaa (MeV^-1), the Delta_a related to the mass of the ALP
	//@param: mai (neV), the mass of the ALP incremented by 1 neV
	//@param: absif (bool), if true, the photon is absorbed by the ALP
	//@return: Prob
	ldir*=conv;
	bdir*=conv;
	FILE* file;
	long double Delta0, zi, dz0, dz, Delta1, Delta18;
	//I apply the initial conditions to the density matrix
	mat_complex rho;
	mat_complex Erho;
	mat_complex rhoidzak;
	mat_complex Mkzi;
	mat_complex ak;
	mat_complex k_RK[8];
	double scrz, scrp;
	int Digs = 16;
	long double cosbdir = cos(bdir);
	vec_real kdir = {cos(ldir)*cosbdir , sin(ldir)*cosbdir , sin(bdir)};		//Direction vector definition
	long double div = N*distz;													//Number of average divisions desired as a function of the distance to be covered
	long double p = 0.0;
	int i = 1;
	while ((p>0.1 && i<=3)||(p>0.5 && i<=6 && i>3 )|| p<=0 || isnan(p)){
		if((file=fopen("1D_plot.txt","w+"))==NULL)
			cout<<"Error opening the file"<<endl;
		if(BCOMPARISON==0)
			fprintf(file, "1\n");
		else
			fprintf(file, "1    %d\n", BCOMPARISON);
		rho = ALPs_initial_conditions();
		Erho = ALPs_initial_conditions();
		zi = distz;
		dz0 = -distz/(div*i*mai);      //Average step for the Runge Kutta
		dz = dz0;
		Delta1 = pow(abs(dz0),errororder);    //Estimated desired error as a function of the order of the RK method and the base average step
		Delta18 = Delta1/8.0;
		while (zi>0){				//I solve the Von Neuman commutator equation with Runge-Kutta method
			//initialize all k_RK with init_mat_complex
            for (int j=0; j<8; j++){
                k_RK[j] = init_mat_complex(3,3);
            }
			for(int j=0; j<order;j++){
				Mkzi = Mk(zi + CC_RK8[j]*dz, kdir, g, omega, Deltaa, absif);
				ak = init_mat_complex(3,3);
				for(int k = 0; k < j; k++){
                        ak = sum_mat(ak , mult_mat_scalar(k_RK[k],AA_RK8[j][k]));
                }
                rhoidzak = sum_mat(rho,mult_mat_scalar(ak,dz));
                k_RK[j] = mult_mat_scalar(commutator(Mkzi,rhoidzak),(0.0,-dz));  //!TO DO: check the warning
				}
            for(int j=0; j<order; j++){
                rho = sum_mat(rho,mult_mat_scalar(k_RK[j],BB_RK8[j]));
                Erho = sum_mat(Erho,mult_mat_scalar(k_RK[j],BB_RK8[j]));
            }
			rho = ALPs_density_matrix_condition(rho);
			Erho = ALPs_density_matrix_condition(Erho);
			Delta0 = fmax(max_value(matrix_abs(sub_mat(rho,Erho))),0.0);
			scrz = double(zi);
			p = 1.0 - rho[2][2].real();
			if(p>0)
				scrp=double(p);
			else
				scrp=0.0;
			if(BCOMPARISON==0)
				fprintf(file,"%.*e   %.*e\n", Digs ,  scrz, Digs ,  scrp);
			else
				fprintf(file,"%.*e   %.*e   %.*e\n", Digs ,  scrz, Digs ,  scrp, Digs , double(BPRINT));
			zi+=dz;
			if (Delta0<Delta18)
				dz = dz0*1.8;
			else
				dz = 0.9*dz0*fmin(fmax(0.3 , pow(Delta1/(2.0*Delta0) , 0.5 )),2.0);
		}
		p = 1.0 - rho[2][2].real();
		i+=1;
		fclose(file);
	}
	return p;
}
#endif //ALP_PHOTON_CONVERSION_H
//end of alp_photon_conversion.h