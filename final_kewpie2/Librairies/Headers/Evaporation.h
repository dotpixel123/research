#ifndef EVAPORATION_H
#define EVAPORATION_H

#include <Readinput.h>
#include <Element.h>
#include <Particle.h>
#include <Leveldensity.h>
#include <Transcoef.h>
#include <cstdlib>
#include <Constantes.h>

using namespace std;

class Evaporation {
	
	private:
		Leveldensity* lev;
		double evawidth;

	public:
		double getEvaWidth();
		double SeCal(double,double,double);
		double BcCal(int,int,int,int);
		void Weisskopf_Ewing(Readinput*,Element*,Element*,Particle*,double,double,double,double);
		void Hauser_Feshbach(Readinput*,Element*,Element*,double,double,double,double,double,double);
		Evaporation();
		~Evaporation();
};

#endif
