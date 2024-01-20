#ifndef READINPUT_H
#define READINPUT_H

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>

using namespace std;

class Readinput {
 
	private:

		ifstream macro;
		size_t NumPoint;
		int Zpro,Apro,Ztar,Atar,Shell,LevDensPar,MollerTable,Pair;
		int Eva,OptPot;
		int Kramers,CollEnhFac;
		int FisBar,FisTrans,Number;
		int Jmin,Jmax,n,p,Gamma,Time,Fuse;
		int GammaModel;
		double KinStep,afan,Friction,Omegags,Omegasd,Emin,Emax,Estep,ShellFactor,CstShell,CstDefor;
		double Delay,TimeFact,CutTcoef,CutWidth,Eview,FacInerMom;
		double DeltaB0,DeltaBvar;
		string reader,Output,FusInput,MassTable,FisBarTable;
	
	public:

		int getZpro();
		int getApro();
		int getZtar();
		int getAtar();
		int getShell();
		int getGammaModel();
		double getShellFactor();
		double getKinStep();
		int getLevDensPar();
		double getCstShell();
		double getafan();
		int getMollerTable();
		int getPair();
		int getEva();
		int getOptPot();
		int getKramers();
		int getCollEnhFac();
		double getFacInerMom();
		double getFriction();
		double getOmegags();
		double getOmegasd();
		int getFisBar();
		int getFisTrans();
		double getEmin();
		double getEmax();
		double getEstep();
		int getJmin();
		int getJmax();
		int getn();
		int getp();
		int getNumber();
		int getGamma();
		int getTime();
		double getTimeFact();
		double getDelay();
		double getCutTcoef();
		double getCutWidth();
		double getEview();
		int getFusion();
		double getDeltaB0();
		double getDeltaBvar();
		size_t getNumberPoint();
		
		string getFusInput();
		string getOutput();
		string getMassTable();
		string getFisBarTable();
		void read();
		int check();

		Readinput();
		~Readinput();

};

#endif
