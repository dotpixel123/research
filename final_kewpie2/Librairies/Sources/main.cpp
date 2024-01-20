#include <Constantes.h>
#include <Element.h>
#include <Readinput.h>
#include <InitMassArray.h>
#include <Writeoutput.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Cascade.h>
#include <Synthesis.h>
#include <Fusion.h>
#include <string>
#include <ctime>
#include <cmath>
#include <cstdlib>

using namespace std;

int main() {

	time_t heure_debut = time(0);
	
	// Reading parameters given in the input file
	Readinput* path = new Readinput();
	path->read();

	// check the input file
	int error = path->check();
	if (error==0) {
		cerr << "Warning: Please check the input file." << endl;
		cin.get();
		exit(EXIT_FAILURE); 
	}
  
	// Declaration of the output file
	Writeoutput* out = new Writeoutput(path->getOutput());
	out->writepara(path);
	string reader = "# Estar\tEcm\tElab\tFusion";
	char* new_name = new char[80];
	if (path->getTime()==1) reader+= "\tFissTime";
	
	for(int i=0;i<=path->getn();i++) {
		for(int j=0;j<=path->getp();j++){
			sprintf(new_name,"\t%dn%dp",i,j);
			reader += new_name;
		}
	}
	delete[] new_name;
	out->writestr(reader);
	out->writeendl();				
  
	// Declaration for the fission time calculation
	int MaxT = 0;	
	if (path->getTime()==1) MaxT = 400; // Time bins 0-400 leading to an error of 10 percent for single neutron chain
  
	// Declaration of the nuclei
	Element* fused = new Element();
  
	int Atar = path->getAtar();
	int Ztar = path->getZtar();
	int Apro = path->getApro();
	int Zpro = path->getZpro();
	int Afus = Atar+Apro;
	int Zfus = Ztar+Zpro;
  
	Table* NuclTab = new Table(Afus,Zfus,path->getn(),path->getp());
  
	// Initialisation of the nuclear data Table
	NuclTab->InitTable(path);

	// Initialization of the compound nucleus
	double massaux = NuclTab->getMassArray(Afus,Zfus);
	double shellaux = NuclTab->getShellArray(Afus,Zfus);
	double beta2aux = NuclTab->getBetaArray(Afus,Zfus);
	double fissbaraux = NuclTab->getFisBarArray(Afus,Zfus);
	fused->Nucleus(path,Afus,Zfus,massaux,shellaux,beta2aux,fissbaraux);
  
	// Visualisation of the differents parameters of the compound nucleus
	cout << endl; 
	cout << fixed << scientific << setprecision(5) << "Compound nucleus: " << endl;
	cout << " A=" << fused->getA() << " Z=" << fused->getZ() << " Name=" << fused->getName() << endl; 
	cout << " R=" << fused->getRadius() << " Rc=" << fused->getCoulRadius() << endl; 
	cout << " Mass=" << fused->getMass() << " Esh=" << fused->getEsh() << endl; 
	cout << " x=" << fused->getx()  << " Is=" << fused->getInertias() << " Ig=" << fused->getInertiag() << endl; 
	cout << " Beta2g=" << fused->getBetag() << " Beta2s=" << fused->getBetas() << endl;
	cout << endl;

	double beta2t=0.;
	double beta2p=0.;
	double Masstar = 0.;
	double Masspro = 0.;

	ifstream macro_mass;
	
	macro_mass.open("../Input/MassTab/MollerNix.dat");
	if (macro_mass.fail()) {
		cerr << "Warning: Problem to open mass table file for projectile and target nuclei." << endl;
		macro_mass.close();
		cin.get();
		exit(EXIT_FAILURE);
	} else {
		do {macro_mass>>reader;} while(reader!="beta2");	
		do {
			int z=0,a=0;
			double mass=0.,beta2=0.,shell=0.;
			macro_mass >> z >> a >> mass >> shell >> beta2; // mass is the mass excess	
			if (z==Ztar && a==Atar) {
				Masstar = Ztar*getBp()+(Atar-Ztar)*getBn()-mass+Ztar*getBe(1)-getBe(Ztar); // Masstar is the Binding Energy
				if (Atar==1) Masstar = 0.;
				beta2t = beta2;
			}
			if (z==Zpro && a==Apro) {
				Masspro = Zpro*getBp()+(Apro-Zpro)*getBn()-mass+Zpro*getBe(1)-getBe(Zpro); // Masspro is the Binding Energy
				if (Apro==1) Masspro = 0.;
				beta2p = beta2; 
			}
		} while (!macro_mass.eof());
	}
	
	// parameters within the barrier-distribution capture model
	double Rt = 1.15*cbrt(Atar*1.);
	double Rp = 1.15*cbrt(Apro*1.);
	double zfus = Zpro*Ztar/(cbrt(1.*Atar)+cbrt(1.*Apro));
	double delta0 = 0.531;
	double C = 0.0421;
	double deltat = Rt*Rt*beta2t*beta2t/(4.*getPi());
	double deltap = Rp*Rp*beta2p*beta2p/(4.*getPi());
	double B0 = 0.853315*zfus + 0.0011695*zfus*zfus - 0.000001544*zfus*zfus*zfus;
	double Bvar = B0*C*sqrt(deltap*deltap+deltat*deltat+delta0*delta0);
	double Rfus_WKB = 1.126*(cbrt(1.*Atar)+cbrt(1.*Apro)); // contact point for WKB
	double Rfus_EBD = 1.16*(cbrt(1.*Atar)+cbrt(1.*Apro)); // contact point for EBD 
	
	if (path->getFusion()==1) {
		cout << "Default barrier height [MeV]: " << B0 << endl;
		cout << "Default barrier width [MeV]: " << Bvar << endl;
		if (path->getDeltaB0()!=0.) {
			B0=path->getDeltaB0();
			cout << "Adjusted barrier height [MeV]: " << B0 << endl;
		}
		if (path->getDeltaBvar()!=0.) {
			Bvar=path->getDeltaBvar();
			cout << "Adjusted barrier width [MeV]: " << Bvar << endl;
		}
		cout << endl;
	}
	macro_mass.close();

	double Qvalue = fused->getMass()-Masstar-Masspro;

	cout << "Binding energies and Q-value calculation: " << endl;
	cout << " Bcomp=" << fused->getMass() << endl;
	cout << " Btar=" << Masstar << endl;
	cout << " Bpro=" << Masspro << endl;
	cout << " Qvalue=" << Qvalue << endl << endl;

	// Declaration and initialisation of the table for the residue cross-sections
	double** pri;
	pri = new double*[(path->getp()+1)];
	for (int i=0;i<(path->getp()+1);i++) {
		pri[i] = new double[(path->getn()+1)];
	}
  
	// Declaration of differents variables
	double* sigfus;
	int JmaxFus = 1000; // Maximal number of partial waves
	sigfus = new double[JmaxFus]; // Declaration of partial cross-sections
	for (int i=0;i<JmaxFus;i++) 
		sigfus[i] = 0.; // Initialization of partial cross-sections
  

	Fusion* fuse = new Fusion(path->getNumberPoint());
  
	// Opening the fusion data file
	ifstream macro;
	if (path->getFusion()==2) {
		reader = "../Input/Fusion/"+path->getFusInput();
		cout << "reader= " << reader << endl;
		macro.open(reader.c_str());
		if (macro.fail()) {
			cerr << "Warning: Cannot open the fusion input file" << endl;
			macro.close();
			cin.get();
			exit(EXIT_FAILURE);
		}
	}
	
	// Beginning of the loop on the excitation energy
	double estar = path->getEmin();
			
	int Jmin = 0;
	int Jmax = 0;
	if (path->getEva()==0) {
		Jmin = path->getJmin();
		Jmax = path->getJmax();
	}
	    
	double fact = double(Atar)/double(Atar+Apro); // factor for converting Elab to Ecm

	do {

		// Calculation of the center-of-mass (Ecm) and laboratory (Elab) energies
		double ecm = estar-Qvalue;
		double elab = ecm/fact;

		sigfus[0] = 0.;
		switch (path->getFusion()) {
			case 0:      
				fuse->WKB_Fusion_J(Apro,Zpro,Atar,Ztar,ecm,Rfus_WKB);
				for (int i=0;i<JmaxFus;i++) {
					if (path->getEva()==1) 
						sigfus[0] += fuse->getSigfus(i);
					else sigfus[i] = fuse->getSigfus(i);
				}        
            	break;
			case 1:
				if (path->getEva()==1) 
					sigfus[0] = fuse->EBD_Fusion(Apro,Zpro,Atar,Ztar,ecm,B0,Bvar,Rfus_EBD);
				else {
					fuse->EBD_Fusion_J(Apro,Zpro,Atar,Ztar,ecm,B0,Bvar,Rfus_EBD);	
					for (int i=0;i<JmaxFus;i++) 
						sigfus[i] = fuse->getSigfus(i);
				}        
            	break;
			case 2:
				// Reading and saving the fusion data file
				macro >> estar >> Jmin >> Jmax;
				for (int i=Jmin;i<=Jmax;i++) {	
					if (path->getEva()==1) {
						double aux = 0.;
						macro >> aux;
						sigfus[0] += aux;
					} else {
						if (i<=path->getJmax()) macro >> sigfus[i-1];
					}
				}
      			if (path->getEva()==1) {
					Jmin = 0;
					Jmax = 0;	
				} else {
					Jmin -= 1;
					Jmax -= 1;
				}
				break;      
			default :
				cerr << "Warning: Please check the fusion calculation option." << endl;
				cin.get();
				exit(EXIT_FAILURE);
		}

		// Conditions on the excitation energy
		if (estar<=path->getEmax() && estar>=path->getEmin() && estar>0.) {

			// Determination of the number of the bins for the energy spectrum of the compound nucleus
			int MaxI = int(estar/path->getEstep()+0.5);
			Cascade* cas = new Cascade(path,(MaxT+1),(MaxI+1),(Jmax+1));

			cout << "MaxI=" << MaxI << "  Estar=" << estar << endl;
		    
			for (int i=0;i<path->getp()+1;i++)
				for (int k=0;k<path->getn()+1;k++) 
					pri[i][k] = 0.;

			// initialization of the compound nucleus
			double sigsum = 0.;
			fused->allocspectrum((MaxT+1),(MaxI+1),(Jmax+1));
			for (int q=Jmin;q<=Jmax;q+=1) {
				fused->setSpectrum(0,MaxI,q,sigfus[q]);
				sigsum += sigfus[q];
			}
							
			// Beginning of the cascade  
			cas->CascadeDecay(path,NuclTab,fused);
			Synthesis* synth = new Synthesis[cas->getComp()];
			cas->Synth(synth);
		
			double popsum = 0.;
			for (int i=0;i<cas->getComp();i++) 
				popsum += synth[i].getTaux();
				
			double FissionTime = 0.;
			if (path->getTime()==1) 
				FissionTime += cas->getFissionTime()/(sigsum-popsum);
  
			// Treatment of the residue cross-sections
			for (int i=0;i<cas->getComp();i++) {
				if (synth[i].getTaux()>0.) {
					int p = fused->getZ()-synth[i].getZ(); // number of emitted protons
					int n = (fused->getA()-fused->getZ())-(synth[i].getA()-synth[i].getZ()); // number of emitted neutrons
					pri[p][n] += synth[i].getTaux();
				}
			}

			delete[] synth;
			fused->deallocspectrum();
	       
			// Writing the residue cross-sections
			cout << "Estar=" << estar << " Ecm=" << ecm << " Elab=" << elab << " Sigfus=" << sigsum << endl;
			out->write(estar);
			out->write(ecm);
			out->write(elab);
			out->write(sigsum);

			if (path->getTime()==1) {
				FissionTime *= getHbar();
				cout << " FissionTime is equal to " << FissionTime << endl;
				out->write(FissionTime);
			}

			for (int i=0;i<=path->getn();i++) {
				cout << " " << i << " " << pri[0][i] << endl;
				for (int j=0;j<=path->getp();j++) 
					out->write((pri[j][i]));
			}
       
			out->writeendl();
			cout << endl;
			delete cas;
		} // end on the excitation energy conditions
    	if (path->getFusion()==2) macro >> reader;
		else {
			if(estar<=path->getEmax()) reader = "continue";
			else reader = "finished";
		} // macro >> reader;
		estar += path->getEview();
	} while (reader=="continue" && estar<=path->getEmax());

	cout << "Computational time: " << difftime(time(0),heure_debut) << endl;
  
	// Destruction of the pointors
	for (int i=0;i<path->getp()+1;i++) 
		delete[] pri[i];
	delete[] pri;
	
	delete[] sigfus;
	delete fused;
	delete NuclTab;
	delete out;
	delete path;
	delete fuse;
  
	// Closing the fusion input data file
	macro.close();
  
	return(0); 
}
