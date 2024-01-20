#ifndef CASCADE_H
#define CASCADE_H

#include <Readinput.h>
#include <InitMassArray.h>
#include <Element.h>
#include <Decay.h>
#include <Synthesis.h>
#include <Particle.h>
#include <cstdlib>

using namespace std;

class Cascade {

	private:
		Element* residuetab;
		Element* tabelem;
		Decay* Dec;
		Particle* evap;
		int* Neu;
		int* Pro;
		double* gammatot;
		
		int n;
		int p;
		int gamma;
		int numpart;
		int Number;
		int MaxI;
		int MaxT;
		int MaxJ;
		int comp;
		double FissionTime;
		double TimeDist;

	public:
		int getComp() {	return(comp); };
		double getFissionTime() { return(FissionTime); };
		double getTimeDist() { return(TimeDist); };
		
		void CascadeDecay(Readinput*,Table*,Element*);
		void Synth(Synthesis*);
		
		Cascade(Readinput*,int,int,int);
		~Cascade();
};

#endif
