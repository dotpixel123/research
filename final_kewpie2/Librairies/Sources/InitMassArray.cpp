#include <InitMassArray.h>

double Table::getShellArray(int A1,int Z1) { 
	return(ShellArray[(A-A1)-(Z-Z1)][Z-Z1]); 
}

double Table::getMassArray(int A1,int Z1) { 
	return(MassArray[(A-A1)-(Z-Z1)][Z-Z1]); 
}

double Table::getBetaArray(int A1,int Z1) { 
	return(BetaArray[(A-A1)-(Z-Z1)][Z-Z1]); 
}

double Table::getFisBarArray(int A1,int Z1) { 
	return(FisBarArray[(A-A1)-(Z-Z1)][Z-Z1]); 
}

void Table::InitTable(Readinput* path) {

	ifstream macro;
	
	if (path->getFisBar()==0) {
		reader = "../Input/FissBar/"+path->getFisBarTable();
		macro.open(reader.c_str());
		if (macro.fail()) {
			cerr << "Warning: There is a problem to open the fission barrier file." << endl;
			macro.close();
			cin.get();
			exit(EXIT_FAILURE);
		} else {
			do {macro>>reader;} while(reader!="Bf");
			do {
				double Bfaux;
				int a,z;
				macro >> z >> a >> Bfaux;
				if (z<=Z && z>=(Z-P) && (a-z)<=(A-Z) && (a-z)>=(A-Z-N))
					FisBarArray[(A-Z)-(a-z)][Z-z] = Bfaux;
			} while (!macro.eof());
		}
		macro.close();
	}

	switch (path->getMollerTable()) {
		case 0:
			macro.open("../Input/MassTab/MollerNix.dat");
			if (macro.fail()) {
				cerr << "Warning: There is a problem to open 'Moller.dat' file." << endl;
				macro.close();
				exit(EXIT_FAILURE);
			} else {
				do {macro>>reader;} while(reader!="beta2");			
				do {
					double shell,mass,beta2;
					int a,z;
					macro >> z >> a >> mass >> shell >> beta2;
					if (z<=Z && z>=(Z-P) && (a-z)<=(A-Z) && (a-z)>=(A-Z-N)) { 
						MassArray[(A-Z)-(a-z)][Z-z] = mass;
						ShellArray[(A-Z)-(a-z)][Z-z] = shell;
						BetaArray[(A-Z)-(a-z)][Z-z] = beta2;
					}    
				} while (!macro.eof());
			}
			macro.close();	
			break;
		case 1:
			reader = "../Input/MassTab/"+path->getMassTable();
			macro.open(reader.c_str());
			if (macro.fail()) {
				cerr << "Warning: There is a problem to open the nuclear data file." << endl;
				macro.close();
				cin.get();
				exit(EXIT_FAILURE);
			} else {
				do {macro>>reader;} while(reader!="beta2");							
				do {
					double shell,mass,beta2;
					int a,z;
					macro >> z >> a >> mass >> shell >> beta2;
					if (z<=Z && z>=(Z-P) && (a-z)<=(A-Z) && (a-z)>=(A-Z-N)) { 
						MassArray[(A-Z)-(a-z)][Z-z] = mass;
						ShellArray[(A-Z)-(a-z)][Z-z] = shell;
						BetaArray[(A-Z)-(a-z)][Z-z] = beta2;
					} 
				} while (!macro.eof());
			}
			macro.close();			
			break;		
		default:
			cerr << "Warning: There is a problem with the mass table option in the input file." << endl;
			macro.close();
			cin.get();
			exit(EXIT_FAILURE);
	}
}

Table::Table(int Afus,int Zfus,int n,int p) {

	A = Afus;
	Z = Zfus;

	N = n+2; // two more for calculating the alpha separation energy of the last element 
	P = p+2; // two more for calculating the alpha separation energy of the last element

	// one more for taking into account the compound nucleus [0][0]
	ShellArray = new double*[N+1];
	MassArray = new double*[N+1];
	BetaArray = new double*[N+1];
	FisBarArray = new double*[N+1];
	for (int i=0;i<N+1;i++) {
		ShellArray[i] = new double[P+1];
		MassArray[i] = new double[P+1];
		BetaArray[i] = new double[P+1];
		FisBarArray[i] = new double[P+1];
		for (int j=0;j<P+1;j++) {
			ShellArray[i][j] = 0.;
			MassArray[i][j] = 0.;
			BetaArray[i][j] = 50.;
			FisBarArray[i][j] = 0.;
		}
	}
}

Table::~Table() {

	for (int i=0;i<N+1;i++) {
		delete[] ShellArray[i];
		delete[] MassArray[i];
		delete[] BetaArray[i];
		delete[] FisBarArray[i];
	}

	delete[] ShellArray;
	delete[] MassArray;
	delete[] BetaArray;
	delete[] FisBarArray;
	
	ShellArray = NULL;
	MassArray = NULL;
	BetaArray = NULL;
	FisBarArray = NULL;
}
