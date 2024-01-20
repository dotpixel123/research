#ifndef FISSION_H
#define FISSION_H

#include <InitMassArray.h>
#include <Readinput.h>
#include <Element.h>
#include <Constantes.h>
#include <Leveldensity.h>
#include <iostream>
#include <cstdlib>

using namespace std;

class Fission {

	private:
		double fiswidth;
		Leveldensity* lev;

	public:
		double getFisWidth();
		void FisWidth(Readinput*,Element*,double,double,double,double);
		Fission();
		~Fission();
};

#endif
