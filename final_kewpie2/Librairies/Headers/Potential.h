#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <iostream>
#include <cmath>
#include <complex>
#include <Constantes.h>
#include <Readinput.h>
#include <cstdlib>

using namespace std;

class Potential {

	private:
	
		double Rmin;	// Starting point for the IWBC, minimum of the potential pocket 
		double Rcoul;	// Coulomb barrier position independent of lb or jb
		double Bmin;	// Minimum of the potential pocket 
		double Bcoul;	// Coulomb barrier height independent of lb or jb

	public:

		complex<double> OpticalPot(Readinput*,double,int,int,int,int,double,double,double);
		double FusPot(double,int,int,int,int,double);
		
		void PotShape(Readinput*,int,int,int,int,double,double,double);
		void CoulombShape(Readinput*,int,int,int,int,double);
		
		double getRmin();
		double getRcoul();
		double getBmin();
		double getBcoul();
		
		double Fusion(int,int,int,int,double);
		
		double frra(double,double,double);
		double dfrradr(double,double,double);
		double ddfrradrdr(double,double,double);
		double Vcoul(Readinput*,double,int,int,int,int);		
		double dVcouldr(Readinput*,double,int,int,int,int);
		double Vcentri(double,double,double);		
		double dVcentridr(double,double,double);
		double dV2_avrdr(double,int,int,double);
		double dV2_vardr(double,int,int,int,double,double,double);
		double dV2_kondr(double,int,int,int,double,double,double);
		double Vtot(Readinput*,double,int,int,int,int,double,double,double);
		double dVtotdr(Readinput*,double,int,int,int,int,double,double,double);
		
		double V2_avr(double,int,int,double);
		double W2_avr(double,int,int,double);
		double V2_var(double,int,int,int,double,double,double);
		double W2_var(double,int,int,int,double,double,double);
		double V2_kon(double,int,int,int,double,double,double);
		double W2_kon(double,int,int,int,double,double,double);
		
		Potential();
		~Potential();

};

#endif
