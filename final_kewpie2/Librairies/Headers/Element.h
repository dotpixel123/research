#ifndef ELEMENT_H
#define ELEMENT_H

#include <cmath>
#include <Constantes.h>
#include <Transcoef.h>
#include <Readinput.h>
#include <string>
#include <cstdlib>
#include <iostream>

using namespace std;

class Element {

	private:
		int A,Z;
		int maxt,maxi,maxj;
		int nbmax,jbmax,emax;
		string Name;
		double Radius,CoulRadius,Mass,MassExcess,x,Esh,Inertias,Inertiag,Pairing,Betag,Betas,Pop,Bf;
		
		int* lbmax; // Cut-off values of lb for neutron, proton and alpha 
		double*** Spectrum; // Energy spectrum [t][i][j]
		double**** TabCoef; // Table of transmission coefficients
		Transcoef* trans;

	public:
		int getA();
		void setA(int);
		int getZ();
		void setZ(int);
		
		string getName();
		
		double getRadius();
		double getCoulRadius();
		double getMass();
		double getMassExcess();
		double getEsh();		
		void setEsh(double);
		
		double getx();
		double Bfcal(int,int,int,double);
		double FisBar1(int,int);
		double FisBar2(int,int);
		double getFisBarrier();

		double getInertias();
		double getInertiag();
		
		double getPairing();
		
		double getBetag();
		double getBetas();
		
		double getPop();
		void setPop(double);
		
		void ElementName(int);
		double Fissility(int,int);
		double InertiaMoment(int,double);
		double PairingEnergy(int,int);
		
		void allocspectrum(int,int,int);
		void deallocspectrum();		
		
		int getMaxt();
		void setMaxt(int);
		int getMaxi();
		void setMaxi(int);
		int getMaxj();
		void setMaxj(int);
		
		double getSpectrum(int,int,int);
		void setSpectrum(int,int,int,double);
		
		void Nucleus(Readinput*,int,int,double,double,double,double);
		
		void InitElement();
		void InitSpectrum();
		
		int getlbmax(int);
		double getTransCoefArray(int,int,int,int);
		void InitTransCoef(Readinput*);
		void CutAngularMom(Readinput*);
		void dealloctranscoef();

	Element();
	Element& operator+=(Element&);  
	Element& operator=(Element&);
	Element& operator=(Element*&);
	~Element();
};

#endif
