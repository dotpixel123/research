#include <Particle.h>

int Particle::getA() {
	return(A);
}

void Particle::setA(int A1) {
	A = A1;
}

int Particle::getZ() {
	return(Z);
}

void Particle::setZ(int Z1) {
	Z = Z1;
}

double Particle::getSpin() {
	return(Spin);
}

void Particle::setSpin(double Spin1) {
	Spin = Spin1;
}

double Particle::getMass() {
	return(Mass);
}

void Particle::setMass(double Mass1) {
	Mass = Mass1;
}

void Particle::InitParticle(int A1,int Z1) {

	A = A1;
	Z = Z1;

	double MassExcess = 0.;
	switch (A) {
		case 1 :
    		if (Z==1) MassExcess = getBp(); // Proton
    		else MassExcess = getBn(); // Neutron
    		Spin = 0.5;
    		break;
  		case 4 :
    		MassExcess = getBa(); // Alpha
    		Spin = 0.;
    		break;
  		default : 
    		cerr << "Please add the mass and spin of your emitted particle in the 'Particle.cpp' file." << endl;
    		cin.get();
    		exit(EXIT_FAILURE);
	}
  
	Mass = Z*getBp()+(A-Z)*getBn()-MassExcess+Z*getBe(1)-getBe(Z); // Nuclear binding energy
	if (A==1) Mass = 0.; // for Neutrons and Protons
}

Particle::Particle() {
}

Particle::~Particle() {
}
