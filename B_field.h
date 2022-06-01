//header file for B_field.cpp
#ifndef B_FIELD_H
#define B_FIELD_H

//library for the calculation of different models of galactic magnetic field
//include the header file for the calculation of the magnetic field
#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>

using namespace std;

//declare global variables
const long double PI = 3.14159265358979324;
const long double IFAC = 11.5*PI/180.0;
const long double SIFAC = sin(IFAC);
const long double CIFAC = cos(IFAC);
const long double BB = tan((PI/2)-IFAC);
const long double BRING = 0.1;
const vector<long double> BCOMP = {0.1, 2.0, -0.9, 2.0, -3.0, -3.5, 0.0, 1.92}; //!DIFFERENZA CON GAMMAALPS [0.1, 3., -0.9, -0.8, -2.0, -4.2, 0., np.nan]
const vector<long double> RX = {5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5};
const int LENRX = RX.size();
const long double HDISK = 0.4;
const long double WDISK = 0.27;
const long double BN = 1.0;
const long double BS = -0.8;
const long double RN = 9.22;
const long double RS = 16.7;
const long double WH = 0.2;
const long double Z0 = 5.3;
const long double RXP = 2.9;
const long double RXC = 4.8;
const long double BBX = 3.0;
const long double THETAX0 = 49.0*PI/180.0;
const long double TANTHETAX0 = tan(THETAX0);
const long double DELTA = 3.0;
const long double R_SPIRAL1 = 5.0;
const long double R_SPIRAL2 = 20.0;
const long double BSTR = 3.0;
const long double RMAG = 5.0;
const long double ZMAG = 7.0;
const long double ZMIN = 0.1;
const long double ZS = 1.0;
const long double BETA0 = 2.0*(RXC* TANTHETAX0)/ZS;
const long double RSC = RXC+ZS/TANTHETAX0;
int BTYPE = 0;
int BFERMIBUBBLE = 0;
int BRINGMUL = 0;
int OMEGAIF = 0;
int ABSIF = 0;
int BCOMPARISON = 0;

long double selec(int i, long double r, long double phi){
    //! This function select the correct spiral arm for the calculation of the magnetic field
    //! for the JF12 model.
    //@param i: the spiral arm number
    //@param r: the radial coordinate
    //@param phi: the azimuthal coordinate
    //@return: the minimum distance between the point and the spiral arm center
	long double r0=abs(RX[i]*exp((phi-5*PI)/BB)-r);
	long double r1=abs(RX[i]*exp((phi-3*PI)/BB)-r);
	long double r2=abs(RX[i]*exp((phi-PI)/BB)-r);
	long double r3=abs(RX[i]*exp((phi+PI)/BB)-r);
	long double r4=abs(RX[i]*exp((phi+3*PI)/BB)-r);
	return min(min(min(min(r0,r1),r2),r3),r4);
}

long double thetaxfun(long double z, long double r, long double rp){
    //! This function calculates the angle between the magnetic field and the x-axis.
    //@param z: the vertical coordinate
    //@param r: the radial coordinate
    //@param rp: the radial coordinate of the spiral arm
    //@return: the angle between the magnetic field and the x-axis
	if (z!=0 || rp<0)
		return atan(abs(z)/abs(r-rp));
	else if(rp>0)
		return 0.0;
	else
		return PI/2.0;
}

long double BcompK(long double r ,long double phi){
    //! This function calculates a projection of the magnetic field spiral component.
    //@param r: the radial coordinate
    //@param phi: the azimuthal coordinate
    //@return: the projection of the magnetic field spiral component in \mu G
	long double test = 1000;
	int k = 0;
	long double s;
	for(int i=0; i<LENRX; i++){
		s = selec(i,r,phi);
		if (test > s){
			test = s;
			k = i;
		}
	}
	return BCOMP[k];
}

long double Hs(long double r, long double phi, long double absz){
    //! This function calculates the integral of the magnetic field over a circular path in
    //! the  space around the center of the galaxy where there is the spiral component of
    //! the magnetic field.
    //@param r: the radial coordinate
    //@param phi: the azimuthal coordinate
    //@param absz: the absolute value of the vertical coordinate
    //@return: integral over the specified path
	if(phi==0)
		return 0;
	else{
		long double integral = 0;
		if(phi<0)
			phi+=2*PI;
		vector <long double> a = {0,0,0,0,0,0,0,0};
		for(int i=0; i<LENRX-1; i++){
			a[i] = BB * (log((2*r)/(RX[i]+RX[i+1])))-PI;
			while(a[i]<0 || a[i]>2*PI){
				if(a[i]<0)
					a[i]+=2*PI;
				else
					a[i]-=2*PI;
			}
		}
		a[LENRX-1] = BB * (log((2*r)/(RX[LENRX-1]*exp(-PI/BB)+RX[0]*exp(PI/BB))));
		while(a[LENRX-1]<0 || a[LENRX-1]>2*PI){
			if(a[LENRX-1]<0)
				a[LENRX-1]+=2*PI;
			else
				a[LENRX-1]-=2*PI;
		}
		sort(a.begin(), a.end());
		if(phi<a[0])
			integral = BcompK(r,phi)*phi;
		else if(phi > a[LENRX-1]){
			integral = BcompK(r, phi)*(phi-a[LENRX-1]);
			for(int i = 0; i<LENRX-1; i++){
				integral += BcompK(r, (a[i]+a[i+1])/2)* (a[i+1]-a[i]);
			}
			integral += BcompK(r, 0)*a[0];
		}
		else{
			integral = BcompK(r, 0)*a[0];
			int i = 0;
			while(phi >= a[i+1]){
				integral += BcompK(r, (a[i]+a[i+1])/2) * (a[i+1]-a[i]);
				i++;
			}
			integral += BcompK(r, phi)*(phi-a[i+1]);
		}
		return integral*(1-1/(1 + exp(-2*(absz-HDISK)/WDISK)));
	}
}

vector<long double> JF12(long double x, long double y, long double z){
    //! This function calculates the magnetic field components of the JF12 model.
    //@param x (kpc): the x-coordinate
    //@param y (kpc): the y-coordinate
    //@param z (kpc): the z-coordinate
    //@return: the magnetic field components in \mu G
	long double r = sqrt(( pow (x,2.0) + pow (y,2.0)));
	long double rz2 = pow(r,2.0) + pow(z,2.0);
	long double absz = abs(z);
	long double phi = 0.0;
	vector<long double> bcentral = {0.0,0.0,0.0};
	vector<long double> bspiral = {0.0,0.0,0.0};
	vector<long double> btoroidal = {0.0,0.0,0.0};
	vector<long double> bX = {0.0,0.0,0.0};
	vector<long double> btot = {0.0,0.0,0.0};
	long double thetax = 0.0;
	long double Bpol = 0.0;
	if(x<0)
		phi = PI + atan (y/x);
	else if(y<0.0){
		if (x==0.0)
			phi = 3.0*PI/2.0;
		else
			phi = 2.0*PI + atan(y/x);
	}
	else {
		if (x==0.0)
			phi = PI/2.0;
		else
			phi = atan(y/x);
	}
	//Calculation of central correction from https : // arxiv.org/pdf/1506.05334.pdf
	if (rz2<=1){
		long double bcent = 4*exp(-absz);
		bcentral = {bcent,0.0,bcent};
	}
	//Calculation of the spiral component
	if(r>3.0 && r<=5.0 && BRINGMUL==1)
		bspiral = {0.0,BRING*(1.0-1.0/(1.0 + exp(-2.0*(absz-HDISK)/WDISK))),0.0};
	else if(r>5 && r<20){
		long double bspir = BcompK(r,phi)*(1.0-1.0/(1.0 + exp(-2.0*(absz-HDISK)/WDISK)))*5.0/r;
		bspiral = {SIFAC*bspir , CIFAC*bspir , 0.0};
	}
	//Calculation of the toroidal component
	if(z>0.0)
		btoroidal={0.0,BN*(1.0-1.0/(1.0+exp(-2.0*(r-RN)/WH)))*exp(-absz/Z0)/(1.0+exp(-2.0*(absz-HDISK)/WDISK)),0.0};
	else if(z<0.0)
		btoroidal={0.0,BS*(1.0-1.0/(1.0+exp(-2.0*(r-RS)/WH)))*exp(-absz/Z0)/(1.0+exp(-2.0*(absz-HDISK)/WDISK)),0.0};
	//Calculation of the X component
	if (rz2>1){
		long double rp = r*RXC/(RXC + absz/TANTHETAX0);
		if (rp>RXC){
			thetax = THETAX0;
			Bpol = BBX*exp((absz/TANTHETAX0-r)/RXP)*(r-absz/TANTHETAX0)/r;
		}
		else{
			thetax =thetaxfun(z,r,rp);
			Bpol = BBX*exp(-rp/RXP)*pow((rp/r),2);
		}
	}

	//Calculation of bX
	if(z>0)
		bX={Bpol*cos(thetax),0.0,Bpol*sin(thetax)};
	else
		bX={-Bpol*cos(thetax),0.0,Bpol*sin(thetax)};
	//Calculation of btot
	btot={bcentral[0]+bspiral[0]+bX[0],bspiral[1]+btoroidal[1],bcentral[2]+bX[2]};
	long double cphi=cos(phi);
	long double sphi=sin(phi);
	vector<long double> result={btot[0]*cphi-btot[1]*sphi,btot[0]*sphi+btot[1]*cphi,btot[2]};
	return result;
}

vector<long double> JF12CORRECTED(long double x, long double y, long double z){
    //! This function calculates the magnetic field components of the JF12 model with solenoidal corrections.
    //@param x (kpc): the x-coordinate
    //@param y (kpc): the y-coordinate
    //@param z (kpc): the z-coordinate
    //@return: the magnetic field components in \mu G
    long double r = sqrt(( pow (x,2.0) + pow (y,2.0)));
	long double rz2 = pow(r,2.0) + pow(z,2.0);
	long double absz = abs(z);
	long double phi = 0.0;
	vector<long double> bcentral = {0.0,0.0,0.0};
	vector<long double> bspiral = {0.0,0.0,0.0};
	vector<long double> btoroidal={0.0,0.0,0.0};
	vector<long double> bX={0.0,0.0,0.0};
	vector<long double> btot={0.0,0.0,0.0};
	long double thetax=THETAX0;
	long double Bpol=0.0;
	//Calculation of phi
	if(x<0)
		phi = PI + atan (y/x);
	else if(y<0.0){
		if (x==0.0)
			phi = 3.0*PI/2.0;
		else
			phi = 2.0*PI + atan(y/x);
	}
	else {
		if (x==0.0)
			phi = PI/2.0;
		else
			phi = atan(y/x);
	}
	//Calculation of central correction from https : // arxiv.org/pdf/1506.05334.pdf
	if (rz2<=1){
		long double bcent = 4*exp(-absz);
		bcentral = {bcent,0.0,bcent};
	}
	//Calculation of the plane component
	if(r>3.0 && r<=R_SPIRAL1 && BRINGMUL==1)
		bspiral = {0.0,BRING*(1.0-1.0/(1.0 + exp(-2.0*(absz-HDISK)/WDISK))),0.0};
	else if(r>R_SPIRAL1 && r<R_SPIRAL1+DELTA){
		long double ralpha = R_SPIRAL1;
		long double rbeta = R_SPIRAL1 + DELTA;
		long double pdr = R_SPIRAL1/rbeta*(2.0 - r/rbeta + (ralpha/rbeta - 2.0)* pow((r - rbeta)/(ralpha - rbeta),2.0));
		long double bspiral1 = BcompK(r,phi)*(1.0-1.0/(1.0 + exp(-2.0*(absz-HDISK)/WDISK)))*pdr;
		long double qdr = R_SPIRAL1/rbeta * (2.0 - 2.0*r/rbeta + (ralpha/rbeta-2.0)*( (3.0*pow(r,2.0) - 4.0*rbeta*r + pow(rbeta,2.0)) / pow((ralpha - rbeta),2.0) ));
		long double H = Hs(r,phi,absz);
		bspiral = {SIFAC*bspiral1 , CIFAC*bspiral1 - qdr * H * SIFAC, 0.0};
	}
	else if(r>R_SPIRAL2-DELTA && r<R_SPIRAL2){
		long double ralpha = R_SPIRAL2;
		long double rbeta = R_SPIRAL2 - DELTA;
		long double pdr = R_SPIRAL1/rbeta*(2.0 - r/rbeta + (ralpha/rbeta - 2.0)* pow((r - rbeta)/(ralpha - rbeta),2.0));
		long double bspiral1 = BcompK(r,phi)*(1.0-1.0/(1.0 + exp(-2.0*(absz-HDISK)/WDISK)))*pdr;
		long double qdr = R_SPIRAL1/rbeta * (2.0 - 2.0*r/rbeta + (ralpha/rbeta-2.0)*( (3.0*pow(r,2.0) - 4.0*rbeta*r + pow(rbeta,2.0)) / pow((ralpha - rbeta),2.0) ));
		long double H = Hs(r,phi,absz);
		bspiral = {SIFAC*bspiral1 , CIFAC*bspiral1 - qdr * H * SIFAC, 0.0};
	}
	else if(r>=R_SPIRAL1+DELTA && r<=R_SPIRAL2-DELTA){
		long double bspir = BcompK(r,phi)*(1.0-1.0/(1.0 + exp(-2.0*(absz-HDISK)/WDISK)))*5.0/r;
		bspiral = {SIFAC*bspir , CIFAC*bspir , 0.0};
	}
	//Calculation of the toroidal component
	if(z>0.0)
		btoroidal={0.0,BN*(1.0-1.0/(1.0+exp(-2.0*(r-RN)/WH)))*exp(-absz/Z0)/(1.0+exp(-2.0*(absz-HDISK)/WDISK)),0.0};
	else if(z<0.0)
		btoroidal={0.0,BS*(1.0-1.0/(1.0+exp(-2.0*(r-RS)/WH)))*exp(-absz/Z0)/(1.0+exp(-2.0*(absz-HDISK)/WDISK)),0.0};
	//Calculation of the X component
	if (fabs(z)<ZS){
		long double rs = r*pow(1.0 - (1.0/(2.0 + BETA0))*(1-pow(z/ZS,2.0)) , -1.0);
		long double rp = rs*RXC/(RXC + ZS/TANTHETAX0);
		if (rs>RSC){
			Bpol = BBX*exp((absz/TANTHETAX0-rs)/RXP)*(rs-absz/TANTHETAX0)/rs;
		}
		else{
			thetax = thetaxfun(ZS,rs,rp);
			if(rs==0.0)
				Bpol = 1.0;
			else
				Bpol = BBX*exp(-rp/RXP)*pow((rp/rs),2);
		}
		long double F = pow(1.0 - 1/(2.0 + BETA0)*(1-pow(z/ZS,2.0)) , -2.0);
		if(z>0)
			bX={(z/ZS)*Bpol*cos(thetax)*F,0.0,Bpol*sin(thetax)*F};
		else
			bX={-(z/ZS)*Bpol*cos(thetax)*F,0.0,Bpol*sin(thetax)*F};
	}
	else{
		if (rz2>1){
			long double rp = r*RXC/(RXC + absz/TANTHETAX0);
			if (rp>RXC){
				Bpol = BBX*exp((absz/TANTHETAX0-r)/RXP)*(r-absz/TANTHETAX0)/r;
			}
			else{
				thetax = thetaxfun(z,r,rp);
				Bpol = BBX*exp(-rp/RXP)*pow((rp/r),2);
			}
		}
		if(z>0)
			bX={Bpol*cos(thetax),0.0,Bpol*sin(thetax)};
		else
			bX={-Bpol*cos(thetax),0.0,Bpol*sin(thetax)};
	}
	//Calculation of btot
	btot={bcentral[0]+bspiral[0]+bX[0],bspiral[1]+btoroidal[1],bcentral[2]+bX[2]};
	long double cphi=cos(phi);
	long double sphi=sin(phi);
	vector<long double> result={btot[0]*cphi-btot[1]*sphi,btot[0]*sphi+btot[1]*cphi,btot[2]};
	return result;
}
//constants for the second group of galactic magnetic field models
const long double BPSH = 1 / tan(IFAC);
const long double DASS = -0.6;
const long double DBSS = -0.6;
const long double Z0ASS = 1.0;
const long double Z0BSS = 1.0;
const long double Z0SUN = 1.0;
const long double PFACASS = -5*PI/180.0;
const long double PFACBSS = -6*PI/180.0;
const long double SINPFACASS = sin(PFACASS);
const long double SINPFACBSS = sin(PFACBSS);
const long double COSPFACASS = cos(PFACASS);
const long double COSPFACBSS = cos(PFACBSS);
const long double B0ASS = 2.0;
const long double B0BSS = 2.0;
const long double B0SUN = 2.0;
const long double RODOT = 8.5;
const long double PHIASS = BPSH*log(1+DASS/RODOT)-PI/2;
const long double PHIBSS = BPSH*log(1+DBSS/RODOT)-PI/2;
const long double COSPHIASS = cos(PHIASS);
const long double COSPHIBSS = cos(PHIBSS);
const long double RCPSH = 5.0;
const long double RCSUN = 5.0;
const long double B0HNASS = 4.0;
const long double B0HSASS = 2.0;
const long double B0HNBSS = 4.0;
const long double B0HSBSS = 2.0;
const long double Z0H = 1.3;
const long double Z1H = 0.25;
const long double Z2H = 0.4;
const long double R0H = 8.0;
const long double R0sun = 10.0;


vector<long double> PshASS(long double x, long double y, long double z){
	//! This function returns the galactical magnetic field components for the ASS model of the Pshirkov model
	//@param x (kpc): the x-coordinate
    //@param y (kpc): the y-coordinate
    //@param z (kpc): the z-coordinate
    //@return: the magnetic field components in \mu G
	long double r = sqrt(( pow (x,2.0) + pow (y,2.0)));
	long double absz = fabs(z);
	long double phi;
	//Calculation of phi
	if(x<0)
		phi = PI + atan (y/x);
	else if(y<0.0){
		if (x==0.0)
			phi = 3.0*PI/2.0;
		else
			phi = 2.0*PI + atan(y/x);
	}
	else {
		if (x==0.0)
			phi = PI/2.0;
		else
			phi = atan(y/x);
	}
	//B spiral component
	long double Br;
	if(r>=RCPSH)
		Br = B0ASS*RODOT/(r*COSPHIASS);
	else
		Br = B0ASS*RODOT/(RCPSH*COSPHIASS);
	long double Bspir = Br*fabs(cos(phi - BPSH*log(r/RODOT)+PHIASS))*exp(-fabs(z)/Z0ASS);
	//B halo component
	long double BH;
	if(z>0){
		if(absz<Z0H)
			BH = B0HNASS*pow( 1 + pow( (absz-Z0H)/Z1H , 2.0) , -1.0)* r/R0H * exp(1 - r/R0H);
		else
			BH = B0HNASS*pow( 1 + pow( (absz-Z0H)/Z2H , 2.0) , -1.0)* r/R0H * exp(1 - r/R0H);
	}
	else{
		if(absz<Z0H)
			BH = B0HSASS*pow( 1 + pow( (absz-Z0H)/Z1H , 2.0) , -1.0)* r/R0H * exp(1 - r/R0H);
		else
			BH = B0HSASS*pow( 1 + pow( (absz-Z0H)/Z2H , 2.0) , -1.0)* r/R0H * exp(1 - r/R0H);
	}
	vector<long double> btot = {(Bspir + BH)* COSPFACASS, (Bspir + BH)* SINPFACASS, 0.0};
	long double cphi=cos(phi);
	long double sphi=sin(phi);
	vector<long double> result={btot[0]*cphi-btot[1]*sphi,btot[0]*sphi+btot[1]*cphi,0.0};
	return result;
}
vector<long double> PshBSS(long double x, long double y, long double z){
	//! This function returns the galactical magnetic field components for the BSS model of the Pshirkov model
	//@param x (kpc): the x-coordinate
    //@param y (kpc): the y-coordinate
    //@param z (kpc): the z-coordinate
    //@return: the magnetic field components in \mu G
	long double r = sqrt(( pow (x,2.0) + pow (y,2.0)));
	long double absz = fabs(z);
	long double phi;
	//Calculation of phi
	if(x<0)
		phi = PI + atan (y/x);
	else if(y<0.0){
		if (x==0.0)
			phi = 3.0*PI/2.0;
		else
			phi = 2.0*PI + atan(y/x);
	}
	else {
		if (x==0.0)
			phi = PI/2.0;
		else
			phi = atan(y/x);
	}
	//B spiral component
	long double Br;
	if(r>=RCPSH)
		Br = B0BSS*RODOT/(r*COSPHIBSS);
	else
		Br = B0BSS*RODOT/(RCPSH*COSPHIBSS);
	long double Bspir = Br*cos(phi - BPSH*log(r/RODOT)+PHIBSS)*exp(-fabs(z)/Z0BSS);
	//B halo component
	long double BH;
	if(z>0){
		if(absz<Z0H)
			BH = B0HNBSS*pow( 1 + pow( (absz-Z0H)/Z1H , 2.0) , -1.0)* r/R0H * exp(1 - r/R0H);
		else
			BH = B0HNBSS*pow( 1 + pow( (absz-Z0H)/Z2H , 2.0) , -1.0)* r/R0H * exp(1 - r/R0H);
	}
	else{
		if(absz<Z0H)
			BH = B0HSBSS*pow( 1 + pow( (absz-Z0H)/Z1H , 2.0) , -1.0)* r/R0H * exp(1 - r/R0H);
		else
			BH = B0HSBSS*pow( 1 + pow( (absz-Z0H)/Z2H , 2.0) , -1.0)* r/R0H * exp(1 - r/R0H);
	}
	vector<long double> btot = {(Bspir + BH)* COSPFACBSS, (Bspir + BH)* SINPFACBSS, 0.0};
	long double cphi=cos(phi);
	long double sphi=sin(phi);
	vector<long double> result={btot[0]*cphi-btot[1]*sphi,btot[0]*sphi+btot[1]*cphi,0.0};
	return result;
}
vector<long double> Sun(long double x, long double y, long double z){
	//! This function returns the galactical magnetic field components for the Sun model
	//@param x (kpc): the x-coordinate
    //@param y (kpc): the y-coordinate
    //@param z (kpc): the z-coordinate
    //@return: the magnetic field components in \mu G
	long double r = sqrt(( pow (x,2.0) + pow (y,2.0)));
	long double absz = fabs(z);
	long double phi;
	//Calculation of phi
	if(x<0)
		phi = PI + atan (y/x);
	else if(y<0.0){
		if (x==0.0)
			phi = 3.0*PI/2.0;
		else
			phi = 2.0*PI + atan(y/x);
	}
	else {
		if (x==0.0)
			phi = PI/2.0;
		else
			phi = atan(y/x);
	}
	//B spiral component
	long double D1, D2, Bspir;
	if(r>RCSUN)
		D1 = B0SUN*exp(-(r-RODOT)/R0sun - absz/Z0SUN);
	else
		D1 = B0SUN*exp(- absz/Z0SUN);
	if((r>7.5) || (r>5 && r<=6))
		D2 = 1;
	else
		D2 = -1;
	Bspir = D1*D2;
	//B halo component
	long double BH;
	if(z>0){
		if(absz<Z0H)
			BH = B0HNASS*pow( 1 + pow( (absz-Z0H)/Z1H , 2.0) , -1.0)* r/R0H * exp(1 - r/R0H);
		else
			BH = B0HNASS*pow( 1 + pow( (absz-Z0H)/Z2H , 2.0) , -1.0)* r/R0H * exp(1 - r/R0H);
	}
	else{
		if(absz<Z0H)
			BH = B0HSASS*pow( 1 + pow( (absz-Z0H)/Z1H , 2.0) , -1.0)* r/R0H * exp(1 - r/R0H);
		else
			BH = B0HSASS*pow( 1 + pow( (absz-Z0H)/Z2H , 2.0) , -1.0)* r/R0H * exp(1 - r/R0H);
	}
	vector<long double> btot = {(Bspir + BH)* COSPFACASS, (Bspir + BH)* SINPFACASS, 0.0};
	long double cphi=cos(phi);
	long double sphi=sin(phi);
	vector<long double> result={btot[0]*cphi-btot[1]*sphi,btot[0]*sphi+btot[1]*cphi,0.0};
	return result;
}

vector<long double> BBubble(long double x, long double y, long double z){
	//! This function returns the galactical magnetic field components of the Fermi Bubble that needs to be added 
	//! to the other models to see the global effect of the Fermi Bubble
	//@param x (kpc): the x-coordinate
    //@param y (kpc): the y-coordinate
    //@param z (kpc): the z-coordinate
	//@return: Fermi Bubble Magnetic field components in \mu G
	vector<long double> B = {0,0,0};
	long double phi;
	//Calculation of phi
	if(x<0)
		phi = PI + atan (y/x);
	else if(y<0.0){
		if (x==0.0)
			phi = 3.0*PI/2.0;
		else
			phi = 2.0*PI + atan(y/x);
	}
	else {
		if (x==0.0)
			phi = PI/2.0;
		else
			phi = atan(y/x);
	}
	long double absz = abs(z);
	vector<long double> bubble = {0.0 , BSTR*exp(-absz/ZMAG)*exp(-ZMIN/absz)*exp(-sqrt(pow(x,2.0) + pow(y,2.0) + pow(z,2.0))/RMAG), 0.0};
	B[0] = - bubble[1]*sin(phi);
	B[1] = bubble[1]*cos(phi);
	return B;
}
void customBfielevaluator(string a){
	//! This function is used initialize the magnetic field model from a file called a
	//@param a: the name of the file
	//@return: none
	//do stuff
	//TODO: implement this
	cout << "This function is not implemented yet" << endl;
	cout << a << endl;
}
vector<long double> customBfield(long double x, long double y, long double z){
	//! This function returns the galactical magnetic field components for the custom model
	//@param x (kpc): the x-coordinate
    //@param y (kpc): the y-coordinate
    //@param z (kpc): the z-coordinate
    //@return: the magnetic field components in \mu G
	//TODO: implement this
	vector<long double> B = {0,0,0};
	B = {x,y,z};
	return B;
}
void Bset(int TYPE){
	//! This function sets the magnetic field model
	//@param TYPE: the type of the model
	//@return: none
	switch(TYPE){
	case 1:
		BTYPE = 0;
		BFERMIBUBBLE = 0;
		break;
	case 2:
		BTYPE = 0;
		BFERMIBUBBLE = 1;
		break;
	case 3:
		BTYPE = 1;
		BFERMIBUBBLE = 0;
		break;
	case 4:
		BTYPE = 1;
		BFERMIBUBBLE = 1;
		break;
	case 5:
		BTYPE = 2;
		BFERMIBUBBLE = 0;
		break;
	case 6:
		BTYPE = 2;
		BFERMIBUBBLE = 1;
		break;
	case 7:
		BTYPE = 3;
		BFERMIBUBBLE = 0;
		break;
	case 8:
		BTYPE = 3;
		BFERMIBUBBLE = 1;
		break;
	case 9:
		BTYPE = 4;
		BFERMIBUBBLE = 0;
		break;
	case 10:
		BTYPE = 4;
		BFERMIBUBBLE = 1;
		break;
	case 11:
		BTYPE = 5;
		BFERMIBUBBLE = 0;
		break;
	default:
		BTYPE = 0;
		BFERMIBUBBLE = 0;
	}
}
vector<long double> Bfield(long double xi, long double yi, long double zi){
	//! This function returns the galactical magnetic field components
	//@param x (kpc): the x-coordinate
    //@param y (kpc): the y-coordinate
    //@param z (kpc): the z-coordinate
    //@return: the magnetic field components in \mu G
	vector<long double> Bvecz;
	switch(BTYPE){
	case 0:
		Bvecz = JF12(xi,yi,zi);
		break;
	case 1:
		Bvecz = JF12CORRECTED(xi,yi,zi);
		break;
	case 2:
		Bvecz = PshASS(xi,yi,zi);
		break;
	case 3:
		Bvecz = PshBSS(xi,yi,zi);
		break;
	case 4:
		Bvecz = Sun(xi,yi,zi);
		break;
	case 5:
		Bvecz = customBfield(xi,yi,zi);
		break;
	}
	if(BFERMIBUBBLE){
		vector<long double> bubble = BBubble(xi,yi,zi);
		Bvecz[0] += bubble[0];
		Bvecz[1] += bubble[1];
	}
	return Bvecz;
}
void omegaifset(bool omega_if){
	//! This function sets the value of omega_if
	//@param omega_if: the value of omega_if
	//@return: none
	OMEGAIF = omega_if;
}
long double omegapl2(long double x, long double y, long double z){
	//! This function returns the value of omega_plasma^2
	//@param x (kpc): the x-coordinate
    //@param y (kpc): the y-coordinate
    //@param z (kpc): the z-coordinate
    //@return: the magnetic field components in \mu 
	//@return: omega_plasma^2
	if (OMEGAIF){
		long double r = sqrt(( pow(x,2.0) + pow (y,2.0)));
		long double gthin = exp (- pow((r-3.7)/1.8,2.0));
		long double galnel = 0.0;
		if (r>17)
			galnel = 0.09*gthin/pow(cosh(z/0.14),2.0);
		else{
			long double gthick = cos(PI*r/34.0)/cos(PI*8.5/34);
			galnel = 0.09*gthin/pow(cosh(z/0.14),2.0)+0.035*gthick/pow(cosh(z/0.95),2.0);
		}
		if (galnel>pow(10,-7))
			return 13.6161*pow(10,-4)*galnel;
		else
			return 13.6161*pow(10,-11);
	}
	else
		return 0;
}
#endif //B_FIELD_H
//end of B_field.h