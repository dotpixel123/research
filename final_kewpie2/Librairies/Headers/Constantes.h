#ifndef CONSTANTES_H
#define CONSTANTES_H

#include <cmath>

inline double getBe(int Z1) {
	return(14.4381e-6*pow(Z1*1.,2.39)+1.55468e-12*pow(Z1*1.,5.35));
}

inline double max_value(double a,double b) {
	return((a<b)?b:a);
}

inline double min_value(double a,double b) {
	return((a<b)?a:b);
}

inline double getSqrtTwo() {
	return(1.41421356237309504880); // sqrt(2) = 1.41421356237309504880
}

inline double getSqrtThree() {
	return(1.73205080756887729353); // sqrt(3) = 1.73205080756887729353
}

inline double getInvSqrtTwo() {
	return(0.707106781186547524401); // 1/sqrt(2) = 1/1.41421356237309504880
}

inline double getInvSqrtThree() {
	return(0.57735026918962576451); // 1/sqrt(3) = 0.57735026918962576451
}

inline double getE() {
	return(2.71828182845904523536); // e = 2.71828182845904523536
}

inline double getPi() {
	return(3.14159265358979323846); // pi = 3.14159265358979323846 
}

inline double getSqrtPi() {
	return(1.7724538509055160273); // sqrt(pi) = 1.7724538509055160273
}

inline double getInvSqrtPi() {
	return(0.56418958354775628695); // 2/sqrt(pi)*0.5 = 1/sqrt(pi) 
}

inline double getLevdens() {
	return(10.);
}

inline double getAmu() {
	return(931.494061);  // MeV/c^2
}

inline double getHbarc() {
	return(197.3269718);  // Mev*fm
}

inline double getHbar() {
	return(6.58211928e-22);  // Mev*s
}

inline double getR0() {
	return(1.17);
}

inline double getR0Coul() {
	return(1.44);
}

inline double getBp() { // Mass excess of proton
	return(7.288970591);
}

inline double getBn() { // Mass excess of neutron
	return(8.071317144);
}

inline double getCoul() {
	return(1.4399);  // e^2/4*Pi*esp0 Mev.fm
}

inline double getPion() {
	return(2.);  // Pion compton wave lenght
}

inline double getSurf() {
	return(14.);  // MeV = 4*pi*(r0)^2
}

inline double getBa() {
	return(2.4249156);  // Mass excess of alpha particle
}

#endif
