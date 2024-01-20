#ifndef DECAY_H
#define DECAY_H

#include <Readinput.h>
#include <Fission.h>
#include <Gamma.h>
#include <Constantes.h>
#include <Evaporation.h>
#include <Element.h>
#include <Particle.h>
#include <iostream>
#include <cstdlib>

using namespace std;

class Decay {
	
	private:
		Fission* fiswidth;
		Evaporation* evawidth;
		Gamma* gamwidth;
		
		double* Se; // separation energy
		double* Bc; // Coulomb barrier
		double*** gammaeva;
		double*** Width;
		bool* tabbool;
		
		int number;
		int maxi;
		int maxj;
		double FisTime;

	public:
		bool getBool(int);
		double getFisSum();
		void WidthCall(Readinput*,Element*,Element*,Particle*,double*,double,bool,int);
		void InitGammaEva();
		void InitWidth();
		
	Decay(int,int,int);
	~Decay();
};

#endif
