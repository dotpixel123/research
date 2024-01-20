#include <Transcoef.h>

double Transcoef::getTransCoef() {
	return(transcoef);
}

void Transcoef::Numerov(Readinput* path,double rmin,double dr,double j,double l,double elab,int Iter,
						int Ares,int Zres,int Aemi,int Zemi,double Fact1) {
	
// *************************************************************************** //  
// ********** Usual boundary conditions within OM **************************** //
// *************************************************************************** //  
/*   
	double ecm = elab/((1.*Ares+1.*Aemi)/(Ares*1.));
	complex<double> Uopt0,Uopt1,Uopt,aux;
	double Fact = Fact1/12.;
	
	if (l==1.) Psi0=-2./15.*pow(dr,(l+1.));
	else Psi0=0.;
    
    Psi1=pow(dr,(l+1.));
    
    double r = dr*2.;
	double r1 = dr;
	double r0 = 0.;   
    
    Uopt1=pot.OpticalPot(path,r1,Ares,Zres,Aemi,Zemi,j,l,elab);
	Uopt=pot.OpticalPot(path,r,Ares,Zres,Aemi,Zemi,j,l,elab);
		 
	aux = Fact*2.*getSqrtThree()*(Uopt1-ecm)+getSqrtThree();
	aux=(aux*aux-1.)*(1.-Fact/12.*(Uopt1-ecm))*Psi1;
	Psi=aux/(1.-Fact/12.*(Uopt-ecm));
    
	Psi0=Psi1;
	Psi1=Psi;
	
	for (int i=3;i<=(Iter+1);i++) {
		
		r = rmin+dr*i;
		r1 = rmin+dr*(i-1);
		r0 = rmin+dr*(i-2);

  		Uopt0=pot.OpticalPot(path,r0,Ares,Zres,Aemi,Zemi,j,l,elab);
		Uopt1=pot.OpticalPot(path,r1,Ares,Zres,Aemi,Zemi,j,l,elab);
		Uopt=pot.OpticalPot(path,r,Ares,Zres,Aemi,Zemi,j,l,elab);	
				 
		aux = Fact*2.*getSqrtThree()*(Uopt1-ecm)+getSqrtThree();
		aux=(aux*aux-1.)*(1.-Fact/12.*(Uopt1-ecm))*Psi1-(1.-Fact/12.*(Uopt0-ecm))*Psi0;
		Psi=aux/(1.-Fact/12.*(Uopt-ecm));
    
		if (i<(Iter+1)) {
			Psi0=Psi1;
			Psi1=Psi;
		}
	}
*/ 
// *************************************************************************** //
// ********** Incoming-wave boundary conditions with real part of OP ********* //
// *************************************************************************** //

	double ecm = elab/((1.*Ares+1.*Aemi)/(Ares*1.));	
	double Fact = Fact1/12.;
	
	complex<double> aux,Phi0;
	complex<double> ai (0.0,1.0);
	complex<double> ak1,ak2,ak3,ak4;
    complex<double> bk1,bk2,bk3,bk4;
	
	double rpp = rmin+dr;
	double rh = rmin+dr*0.5;    
	
	double Vrmin = ecm-pot.Vtot(path,rmin,Ares,Zres,Aemi,Zemi,j,l,elab);
	double Vrpp = ecm-pot.Vtot(path,rpp,Ares,Zres,Aemi,Zemi,j,l,elab);
	double Vrh = ecm-pot.Vtot(path,rh,Ares,Zres,Aemi,Zemi,j,l,elab);
    
	double mu = (Ares*Aemi*1.)/(Ares+Aemi)*getAmu();	
	
	double k = sqrt(2.*mu/(getHbarc()*getHbarc())*abs(Vrmin));
	Psi0 = exp(k*rmin);
	Phi0 = k*Psi0;
    
    if (Vrmin>=0.) {
		k = sqrt(2.*mu/(getHbarc()*getHbarc())*Vrmin);
		Psi0 = exp(-ai*k*rmin);
		Phi0 = -ai*k*Psi0;
    }
    
    double Fac=dr*(2.*mu/(getHbarc()*getHbarc())); 
    
	ak1 = ak1-Fac*Vrmin*Psi0;
	bk1 = dr*Phi0;
	ak2 = ak2-Fac*Vrh*(Psi0+0.5*bk1);
	bk2 = dr*(Phi0+0.5*ak1);
	ak3 = ak3-Fac*Vrh*(Psi0+0.5*bk2);
	bk3 = dr*(Phi0+0.5*ak2);
	ak4 = ak4-Fac*Vrpp*(Psi0+bk3);
	bk4 = dr*(Phi0+ak3);
	
	Psi1 = Psi0+(bk1+2.*bk2+2.*bk3+bk4)/6.;

	for (int i=2;i<=(Iter+1);i++) {			
		double r = rmin+dr*i;
		double r1 = rmin+dr*(i-1);
		double r0 = rmin+dr*(i-2);

  		double V0 = pot.Vtot(path,r0,Ares,Zres,Aemi,Zemi,j,l,elab);
		double V1 = pot.Vtot(path,r1,Ares,Zres,Aemi,Zemi,j,l,elab);
		double Vr = pot.Vtot(path,r,Ares,Zres,Aemi,Zemi,j,l,elab);	
				 
		aux = Fact*2.*getSqrtThree()*(V1-ecm)+getSqrtThree();
		aux = (aux*aux-1.)*(1.-Fact*(V1-ecm))*Psi1-(1.-Fact*(V0-ecm))*Psi0;
		Psi = aux/(1.-Fact*(Vr-ecm));
    
		if (i<Iter+1) {
			Psi0 = Psi1;
			Psi1 = Psi;
		}
	}
/**/
}

void Transcoef::TransCal(Readinput* path,int Ares,int Zres,int Aemi,int Zemi,double ecm,double j,double l) {
  
	transcoef = 0.;
	
	if (ecm>0.) {  
	
		double mu = (Ares*Aemi*1.)/(Ares*1.+Aemi*1.)*getAmu();
		double elab = ecm*(Ares+Aemi)/(Ares*1.);
		double k = sqrt(2.*mu*ecm)/getHbarc();
    
		pot.PotShape(path,Ares,Zres,Aemi,Zemi,j,l,elab); // Determination of rabs and Vrmin for IWBC
		
	//	double Rmin=0.; // for the OM calculation without IWBC
		double Rmin = pot.getRmin();
		double Vrmin = ecm-pot.getBmin();
		double eta = (Zres*Zemi*1./137.)*sqrt(mu/(2.*ecm));
		double dr = 0.05;
	    double Rmax = 50.;
	    double Fact = (dr*dr)*(2.*mu/(getHbarc()*getHbarc())); 
		int Iter = int((Rmax-Rmin)/dr);
	   	Rmax = Rmin+Iter*dr;
			
		gsl_sf_result fcw, gcw, fcwp, gcwp;
		double exp_F = 0.;
		double exp_G = 0.;
	   	int err = 1; // for Coulomb wave function calculations
		int gsl_sf_coulomb_wave_FG_e(const double, const double, const double, const int, 
									gsl_sf_result*, gsl_sf_result*, gsl_sf_result*, gsl_sf_result*, 
									double*, double*);

		double ap = dr*(Iter+1);
	    double kap = k*ap;
	    double am = dr*(Iter-1);
	    double kam = k*am;
	
	    Numerov(path,Rmin,dr,j,l,elab,Iter,Ares,Zres,Aemi,Zemi,Fact);
	    
	    err = gsl_sf_coulomb_wave_FG_e(eta,kam,l,0,&fcw,&fcwp,&gcw,&gcwp,&exp_F,&exp_G);
    	
	    if (err==1) {
			cerr << "Warning: Coulomb wave function are not well calculated." << endl;
			cin.get();
			exit(EXIT_FAILURE);
		}
	    
	    complex<double> zm(gcw.val,fcw.val);
	    complex<double> izm(gcw.val,-1.*fcw.val);
	    
	    err = gsl_sf_coulomb_wave_FG_e(eta,kap,l,0,&fcw,&fcwp,&gcw,&gcwp,&exp_F,&exp_G);
	
	    if (err==1) {
			cerr << "Warning: Coulomb wave function are not well calculated." << endl;
			cin.get();
			exit(EXIT_FAILURE);
		}
	    
		complex<double> zp(gcw.val,fcw.val);
		complex<double> izp(gcw.val,-1.*fcw.val);
		complex<double> tcplx = (izm*Psi-izp*Psi0);
		complex<double> tpcplx = (zm*Psi-zp*Psi0);
		
		complex<double> Smat = tcplx/tpcplx;
	
		transcoef = 1.-abs(Smat*Smat); // |a|^2 = |a*a|
					
		if ((transcoef<1.e-12)||(Vrmin<0.)) transcoef = 0.;
			
		if (transcoef>1. && transcoef<0.) {
			cerr << "Warning: Problem with optical transmission coefficient." << endl;
			cin.get();
			exit(EXIT_FAILURE);
		}
	}
}

Transcoef::Transcoef() {
}

Transcoef::~Transcoef() {
}
