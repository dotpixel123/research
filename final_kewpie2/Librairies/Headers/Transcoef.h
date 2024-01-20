#ifndef TRANSCOEF_H
#define TRANSCOEF_H

#include <Constantes.h>
#include <Potential.h>
#include <Readinput.h>
#include <gsl/gsl_sf_coulomb.h>
#include <iostream>
#include <complex>
#include <cstdlib>
#include <cmath>

using namespace std;

class Transcoef {

	private:
		double transcoef;
		Potential pot;
		complex<double> Psi,Psi0,Psi1;

	public:
		double getTransCoef();
		void TransCal(Readinput*,int,int,int,int,double,double,double);
		void Numerov(Readinput*,double,double,double,double,double,int,int,int,int,int,double);
	Transcoef();
	~Transcoef();
};

#endif
