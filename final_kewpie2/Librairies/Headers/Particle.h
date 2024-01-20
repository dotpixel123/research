#ifndef PARTICLE_H
#define PARTICLE_H

#include <Constantes.h>
#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

class Particle {

	private:
		int A,Z;
		double Spin,Mass,Radius,CoulRadius;

	public:
		int getA();
		void setA(int);
		int getZ();
		void setZ(int);
		
		double getSpin();
		void setSpin(double);
		
		double getMass();
		void setMass(double);

		void InitParticle(int,int);
	
	Particle();
	~Particle();
};

#endif
