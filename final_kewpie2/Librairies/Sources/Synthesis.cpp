#include <Synthesis.h>

int Synthesis::getZ() {
	return(Z);
}
void Synthesis::setZ(int Z1) {
	Z = Z1;
}

int Synthesis::getA() {
	return(A);
}

void Synthesis::setA(int A1) {
	A = A1;
}

double Synthesis::getTaux() {
	return(Taux);
}

void Synthesis::setTaux(double Taux1) {
	Taux = Taux1;
}

Synthesis::Synthesis() {
	A = 0;
	Z = 0;
	Taux = 0.;
}

Synthesis::~Synthesis() {
}
