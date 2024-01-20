#ifndef SYNTHESIS_H
#define SYNTHESIS_H

class Synthesis {

	private:
		int A,Z;
		double Taux;

	public:
		int getZ();
		void setZ(int);
		int getA();
		void setA(int);
		double getTaux();
		void setTaux(double);

	Synthesis();
	~Synthesis();
};

#endif
