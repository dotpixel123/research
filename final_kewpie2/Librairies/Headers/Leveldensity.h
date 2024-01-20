#ifndef LEVELDENSITY_H
#define LEVELDENSITY_H

#include <Constantes.h>
#include <Readinput.h>
#include <Element.h>
#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

class Leveldensity {

	private:
		double levdensj;
		double levdens;
		int* Neutron;
		int* Proton;
	public:
		double beta_zero(double,double,double);
		double getLevDens();
		double getLevelDensityParameter();
		double parameter1(int);
		double parameter2(double,int,double,double,double,int,int,char);
		double parameter3(double,int,double,double,double,int,int,char);
		double parameter4(double,int,double,double,double,int,int,char);
		double CollEnhFac(double,double,double,double,int,int);
		int Levdensj(Readinput*,Element*,double,double,char,char);
	Leveldensity();
	~Leveldensity();
};

#endif
