#ifndef FUSION_H
#define FUSION_H

#include <iomanip>
#include <cstring>
#include <complex>
#include <iostream>
#include <Potential.h>
#include <Constantes.h>
#include <cstdlib>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>

using namespace std;

class Fusion {

	private:
	
		Potential* Pot;
		double* SigFus;
		double* weight;
		double* point;
		size_t Npoint;

	public:
	
		double getWeight(int i) { return weight[i]; };
		double getPoint(int i) { return point[i]; };
		int getNpoint() { return Npoint;};
		double getSigfus(int);
		void Parabo(double*,double*,double*,double&,double&);
		double TWKB(int,int,int,int,double,double,double,double,double&,double&);
		double EBD_Fusion(int,int,int,int,double,double,double,double);
		void WKB_Fusion_J(int,int,int,int,double,double);
		void EBD_Fusion_J(int,int,int,int,double,double,double,double);
		Fusion(size_t);
		~Fusion();
};

#endif
