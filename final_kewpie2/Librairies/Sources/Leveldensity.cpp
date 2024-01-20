#include <Leveldensity.h>

double Leveldensity::getLevDens() {
	return(levdens);
}

double Leveldensity::beta_zero(double Estarnew, double a, double acc) {

	double error = 0.;
	double ratio_new = 0.;
	double ratio_old = sqrt(a*Estarnew);
	
	do {
		ratio_new = sqrt(a*Estarnew*(1.-exp(-ratio_old)));
		error = abs(ratio_old-ratio_new)/ratio_old;
		ratio_old = ratio_new;
	} while (error>acc);

	double beta0=a/ratio_new;
	
	if (Estarnew==0.) beta0=1.e100; 
    return(beta0);
}

double Leveldensity::parameter1(int A) {

	// ********** Level Density Parameter ******************
	//            Simple Fermi Gas Model    
	// *****************************************************
  
	return(double(A)/getLevdens());  
}

double Leveldensity::parameter2(double CstShell,int Shell,double Esh,double Beta2,double Estarnew,int A,int Z,char state) {
  
	// ********** Level Density Parameter ******************
	//             J. Toke W.J. Swiatecki
	//           Nuc. Phys. A. (1981) 141-150
	// *****************************************************
 
	double beta2 = Beta2;
	double alpha2 = beta2*sqrt(5./(4.*getPi()));
	
	double F2 = 1.+2.*alpha2*alpha2/5.-4.*alpha2*alpha2*alpha2/105.-66.*alpha2*alpha2*alpha2*alpha2/175.;  // Bs
	double F3 = 1.+2.*alpha2*alpha2/5.+16.*alpha2*alpha2*alpha2/105.-82.*alpha2*alpha2*alpha2*alpha2/175.; // Bk

	double I = (A*1.-2.*Z)/(A*1.);
	double aux = (1./14.61)*(A+3.114*F2*cbrt(A*A*1.)+5.626*F3*cbrt(A*1.))*(1.-I*I/9.);

	if ((Shell==1)&&(state=='g')&&(Estarnew!=0.)) {
		double fact = 1.-exp(-Estarnew/CstShell);
		aux *= (1.+fact*Esh/Estarnew);
	}
	
	if ((Shell==1)&&(state=='g')&&(Estarnew==0.)) 
		aux *= (1.+Esh/CstShell);
  
	if (aux<0.) {
		aux = 0.;
		cerr << "Warning: The level density parameter cannot be negatif." << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
  
	return(aux);
}

double Leveldensity::parameter3(double CstShell,int Shell,double Esh,double Beta2,double Estarnew,int A,int Z,char state) {
  
	// ********** Level Density Parameter ******************
	//             		Reisdorf
	//     W. Reisdorf, Zeitschrift für Physik A Atoms and 
	//				Nuclei 300, 227 (1981).
	// *****************************************************
	
	double beta2 = Beta2;
	double alpha2 = beta2*sqrt(5./(4.*getPi()));

	double F2 = 1.+2.*alpha2*alpha2/5.-4.*alpha2*alpha2*alpha2/105.-66.*alpha2*alpha2*alpha2*alpha2/175.;
	double F3 = 1.+2.*alpha2*alpha2/5.+16.*alpha2*alpha2*alpha2/105.-82.*alpha2*alpha2*alpha2*alpha2/175.;

	double r0 = 1.153;
	double aux = 0.04543*r0*r0*r0*A+0.1315*r0*r0*F2*cbrt(A*A*1.)+0.1426*r0*F3*cbrt(A*1.);

	if ((Shell==1)&&(state=='g')&&(Estarnew!=0.)) {
		double fact = 1.-exp(-Estarnew/CstShell);
		aux *= (1.+fact*Esh/Estarnew);
	}
	
	if ((Shell==1)&&(state=='g')&&(Estarnew==0.)) 
		aux *= (1.+Esh/CstShell);
  
	if (aux<0.) {
		aux = 0.;
		cerr << "Warning: The level density parameter cannot be negatif." << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
  
	return(aux);
}

double Leveldensity::parameter4(double CstShell,int Shell,double Esh,double Beta2,double Estarnew,int A,int Z,char state) {
	
	// ********** Level Density Parameter ************************
	//     B. Nerlo-Pomorska, K. Pomorski and J. Bartel,
	// 			Phys. Rev. C 74, 034327 (2006).
	// ***********************************************************
	
	double beta2 = Beta2;
	double alpha2 = beta2*sqrt(5./(4.*getPi()));

	double F2 = 1.+2.*alpha2*alpha2/5.-4.*alpha2*alpha2*alpha2/105.-66.*alpha2*alpha2*alpha2*alpha2/175.;
	double F3 = 1.+2.*alpha2*alpha2/5.+16.*alpha2*alpha2*alpha2/105.-82.*alpha2*alpha2*alpha2*alpha2/175.;
	double F4 = 1.-1.*alpha2*alpha2/5.-4.*alpha2*alpha2*alpha2/105.+51.*alpha2*alpha2*alpha2*alpha2/245.;	// Bcoul
	
	double aux = 0.092*A+0.036*F2*cbrt(A*A*1.)+0.275*F3*cbrt(A*1.)-0.00146*Z*Z*F4/cbrt(A*1.);
	
	if ((Shell==1)&&(state=='g')&&(Estarnew!=0.)) {
		double fact = 1.-exp(-Estarnew/CstShell);
		aux *= (1.+fact*Esh/Estarnew);
	}
	
	if ((Shell==1)&&(state=='g')&&(Estarnew==0.)) 
		aux *= (1.+Esh/CstShell);
  
	if (aux<0.) {
		aux = 0.;
		cerr << "Warning: The level density parameter cannot be negatif." << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
  
	return(aux);
}

double Leveldensity::CollEnhFac(double Estarnew,double Inertia,double a,double Beta2,int A1,int Z1) {

	//***************************************************************************
	//  V.I. Zagrebaev et al., Phys. Rev., 2001, vol.C65, p.014607
	//***************************************************************************
	
	double beta2 = Beta2;
	double Ecr = 40.;
	double dcr = 10.;
	double b0 = 0.15;
	double deltab = 0.04;  
	double cut = 1./(1.+exp((Estarnew-Ecr)/dcr));
	double phi = 1./(1.+exp((b0-abs(beta2))/deltab));

	int deltaN = 1000;
	int deltaZ = 1000;
	for (int i=0;i<8;i++) { 
		int aux = 0;   
		int N1 = A1-Z1;
		aux = Neutron[i]-N1;
		if (abs(aux)<deltaN) 
			deltaN = abs(aux);
		aux = Proton[i]-Z1;
		if (abs(aux)<deltaZ) 
			deltaZ = abs(aux);
	}
 	
	double Fact = 1.;
	if (deltaN!=1000 && deltaZ!=1000) { 

		double delta2 = 0.022+0.003*deltaN+0.005*deltaZ;
		double Temp = sqrt(Estarnew/a);
		double K_rot0 = Temp*Inertia;            
		double K_vib0 = 25.*delta2*delta2*K_rot0;
//		K_vib0 = exp(0.0555*pow(A1*1.,2./3.)*pow(Temp,4./3.));

		double K_vib = 1.;
		double K_rot = 1.;
		if ((K_rot0<=1.)&&(K_rot0>=0.)) K_rot = 1.;
		if (K_rot0>1.) K_rot = (K_rot0-1.)*cut + 1.;
		if (K_rot0<0.) cerr << "Warning: The rotational enhancement factor cannot be negatif." << endl;
		if ((K_vib0<=1.)&&(K_vib0>=0.)) K_vib = 1.;
		if (K_vib0>1.) K_vib = (K_vib0-1.)*cut + 1.;
		if (K_vib0<0.) cerr << "Warning: The vibrational enhancement factor cannot be negatif." << endl;
		
		Fact = K_vib*(1.-phi) + K_rot*phi;
	} else {
		cerr << "Warning: There is problem with the collective enhancement factor. Please check it." << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}

	return(Fact);
}

int Leveldensity::Levdensj(Readinput* path,Element* elem,double Estarnew,double J,char state,char react) {

	// ********************************************************************
	//  Grossjean & Feldmeier's formula for two types of nucleons 
	//  Nucl. Phys. A 444 (1985) 113-132.
	// ********************************************************************

	int err = 1;
	
	double Inertia = 0.;
	if (state=='s') Inertia = elem->getInertias();
	if (state=='g') Inertia = elem->getInertiag();
	
	double Beta2 = 0.;
	if (state=='s') Beta2 = elem->getBetas();
	if (state=='g') Beta2 = elem->getBetag();	  

//	double Erot = J*(J+1.)/(2.*Inertia);
	double Estar = Estarnew-elem->getPairing();
	
	if (Estar>=0.) {
	
		double a = 0.;		
		switch (path->getLevDensPar()) {
			case 0:
				a = parameter1(elem->getA());
		    	break;
		    case 1:
				a = parameter2(path->getCstShell(),path->getShell(),elem->getEsh(),Beta2,Estar,elem->getA(),elem->getZ(),state);
				break;
			case 2:
				a = parameter3(path->getCstShell(),path->getShell(),elem->getEsh(),Beta2,Estar,elem->getA(),elem->getZ(),state);
				break;
			case 3:
				a = parameter4(path->getCstShell(),path->getShell(),elem->getEsh(),Beta2,Estar,elem->getA(),elem->getZ(),state);
				break;
			default:
				cerr << "Warning: There is a problem with the option for the level-density parameter." << endl;
				cin.get();
				exit(EXIT_FAILURE);
		}
	  
		if (state=='s' && react=='f') a*=path->getafan();
		
		if (Estar==0.) {
			levdens = getSqrtTwo()*getSqrtPi()*getE()*a*exp(a*Estar)/12.;
			if (path->getEva()==0) levdens = 0.;
		} else {
			double acc = 1.e-6;	// precision of the calculation of beta0
			double beta0 = beta_zero(Estar,a,acc);
			// double beta0 = sqrt(a/Estar); // high energy limit
	 		levdens = getSqrtPi()*exp(beta0*Estar+a/beta0)/(12.*sqrt(beta0*Estar*Estar*Estar))
	 				*(1.-exp(-a/beta0))/sqrt(1.-Estar*beta0*exp(-a/beta0)*0.5);
			if (path->getEva()==0) {
				double M = Inertia/beta0;
				levdens *= (2.*J+1.)*exp(-(J+0.5)*(J+0.5)/(2.*M))/(sqrt(2.*getPi()*M)*M);
			}
		} 	

		switch (path->getCollEnhFac()) {
		    case 0:
				levdens = levdens;
				break;
			case 1:
		        if (state=='s') levdens *= CollEnhFac(Estar,Inertia,a,elem->getBetas(),elem->getA(),elem->getZ());
		        if (state=='g') levdens *= CollEnhFac(Estar,Inertia,a,elem->getBetag(),elem->getA(),elem->getZ());
				break;
			default: 
				cerr << "Warning: Please check the option for the collective enhancement factor." << endl;
				cin.get();
				exit(EXIT_FAILURE);
		}
		err = 0;    
	}
  
	return(err);
}

Leveldensity::Leveldensity() {

	levdens = 0.;

	Neutron = new int[8];
	Neutron[0] = 2;
	Neutron[1] = 8;							 
	Neutron[2] = 20;
	Neutron[3] = 28;
	Neutron[4] = 50;
	Neutron[5] = 82;
	Neutron[6] = 126;
	Neutron[7] = 184;
  
	Proton = new int[8];
	Proton[0] = 2;
	Proton[1] = 8;							 
	Proton[2] = 20;
	Proton[3] = 28;
	Proton[4] = 50;
	Proton[5] = 82;
	Proton[6] = 114;
	Proton[7] = 184;
}

Leveldensity::~Leveldensity() {

	delete[] Neutron;
	delete[] Proton;
	
	Neutron = NULL;
	Proton = NULL;
}
