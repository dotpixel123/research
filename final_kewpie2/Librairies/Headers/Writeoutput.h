#ifndef WRITEOUTPUT_H
#define WRITEOUTPUT_H

#include <fstream>
#include <iostream>
#include <string>
#include <Readinput.h>
#include <cstdlib>

using namespace std;

class Writeoutput {
 
	private:
		ofstream macro;

	public:
		void write(double);
		void writeendl();
		void writestr(string);
		void writepara(Readinput*);

		Writeoutput();
		Writeoutput(string);
		~Writeoutput();
};

#endif
