#ifndef INITMASSARRAY_H
#define INITMASSARRAY_H

#include <fstream>
#include <iostream>
#include <cmath>
#include <Constantes.h>
#include <Readinput.h>
#include <cstdlib>

using namespace std;

class Table {

	private:
		int A,Z,N,P;
		double** ShellArray;
		double** MassArray;
		double** BetaArray;
		double** FisBarArray;
		string reader;

	public:
		double getShellArray(int,int);
		double getMassArray(int,int);
		double getBetaArray(int,int);
		double getFisBarArray(int,int);
		void InitTable(Readinput*);

	Table(int,int,int,int);
	~Table();
};

#endif
