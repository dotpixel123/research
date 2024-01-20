#ifndef GAMMA_H
#define GAMMA_H

#include <Readinput.h>
#include <Element.h>
#include <Leveldensity.h>
#include <Constantes.h>
#include <cmath>
#include <cstdlib>

using namespace std;

class Gamma {

	private:
		double gamwidth;
		Leveldensity* lev;
		double* coefe;
		double* coefm;

	public:
		double getGamWidth();
		double getfe(int);
		double getfm(int);
		
		double RSFE1EGLO(int,int,double,double,double,double); // Radiation strength function for E1 resonance within the EGLO model
		double RSFE1SMLO(int,int,double,double,double,double); // Radiation strength function for E1 resonance within the SMLO model
		double RSFE2(int,int,double); // Radiation strength function for E2 resonance
		double RSFM1(int,int,double); // Radiation strength function for M1 resonance
		
		void RadStrFuncs(Readinput*,Element*,double,double);
		void GamWidth_J(Readinput*,Element*,double,double,double,double,double);
		void GamWidth(Readinput*,Element*,double,double,double);

		Gamma();
		~Gamma();
};

#endif
