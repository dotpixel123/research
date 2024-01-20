#include <Potential.h>

double Potential::getRmin() {
	return(Rmin);
}

double Potential::getRcoul() {
	return(Rcoul);
}

double Potential::getBmin() {
	return(Bmin);
}

double Potential::getBcoul() {
	return(Bcoul);
}

void Potential::CoulombShape(Readinput* path,int Ares,int Zres,int Aemi,int Zemi,double elab) {

	double lb0 = 0.;
	double jb0 = 0.;
	
	double rmin = 0.5; /* Minimum of rmin (the starting point for IWBC) */
	double rmax = 50.; 

// ***************** Determination of Rcoul and Bcoul for charged particles ******************************** //

	double u0 = 0.;
	double u1 = 0.;
	double u = 0.;
	double raux = rmax;

	if (Zemi==2) jb0 = 0.;
	else jb0 = 0.5;
	
	u0 = dVtotdr(path,rmax,Ares,Zres,Aemi,Zemi,jb0,lb0,elab);
	do {
		raux -= 1.;
	   	if (raux<0.) {
	   	    cerr << "Warning: Problem with the local search for the maximum of the potential." << endl;
			cin.get();
			exit(EXIT_FAILURE);        
	   	}
		u1 = dVtotdr(path,raux,Ares,Zres,Aemi,Zemi,jb0,lb0,elab);
	} while ((u0*u1>0.)&&(raux>rmin));

	double ra = raux+1.0;
	double rb = raux;
	double tolk = 1.e-6;
	int Iter = int(log10(abs(rb-ra)/tolk)/log10(2.)+0.5);
	
	for (int i=0;i<Iter;i++) {
		raux = (ra+rb)/2.;
		u = dVtotdr(path,raux,Ares,Zres,Aemi,Zemi,jb0,lb0,elab);
		if (u0*u<0.) rb = raux;
		else ra = raux;                
	}

	Rcoul = raux;
	Bcoul = Vtot(path,Rcoul,Ares,Zres,Aemi,Zemi,jb0,lb0,elab);
	
	double ddv00 = (dVtotdr(path,Rcoul+1.e-5,Ares,Zres,Aemi,Zemi,jb0,lb0,elab)
				-dVtotdr(path,Rcoul-1.e-5,Ares,Zres,Aemi,Zemi,jb0,lb0,elab))/2.e-5;
	if (ddv00>0.) {
       	cerr << "Warning: Something is strange when calculating the Coulomb barrier position." << endl;
		cin.get();
		exit(EXIT_FAILURE);        
	}		
}

void Potential::PotShape(Readinput* path,int Ares,int Zres,int Aemi,int Zemi,double jb,double lb,double elab) {

	double rmin = 0.5; /* Minimum of rmin (the starting point for IWBC) */
	double rmax = 50.; 

// ********************************** Determination of Rmin and Bmin *********************************************** //
	
	double u0 = 0.;
	double u1 = 0.;
	double u = 0.;
	double raux = rmax;
	
	u0 = dVtotdr(path,rmax,Ares,Zres,Aemi,Zemi,jb,lb,elab);	
	do {
	 	raux -= 0.5;    
    	if (raux<0.) {
	   	    cerr << "Warning: Problem with the local search for the maximum of the potential." << endl;
			cin.get();
			exit(EXIT_FAILURE);        
    	}
		u1 = dVtotdr(path,raux,Ares,Zres,Aemi,Zemi,jb,lb,elab);
	} while ((u0*u1>0.)&&(raux>rmin));
	
	double ra = raux+0.5;
	double rb = raux;
	double tolk = 1.e-6;
	int Iter = int(log10(abs(rb-ra)/tolk)/log10(2.)+0.5);
		
	for (int i=0;i<Iter;i++) {
		raux = (ra+rb)/2.;
		u = dVtotdr(path,raux,Ares,Zres,Aemi,Zemi,jb,lb,elab);
		if (u0*u<0.) rb = raux;
		else ra = raux;                
	}

	rb = raux;
	ra = rb-0.5;
	u0 = dVtotdr(path,ra,Ares,Zres,Aemi,Zemi,jb,lb,elab);

	rb = rmin;
    Iter = int(log10(abs(rb-ra)/tolk)/log10(2.)+0.5);
    	
    for (int i=0;i<Iter;i++) {
		raux = (ra+rb)/2.;
		double u = dVtotdr(path,raux,Ares,Zres,Aemi,Zemi,jb,lb,elab);
		if (u0*u<0.) rb = raux;
		else ra = raux;                
	}
		
	Rmin = raux;
	if (Rmin<rmin) Rmin = rmin;
	
	Bmin = Vtot(path,Rmin,Ares,Zres,Aemi,Zemi,jb,lb,elab);
}

// ******************************************************************************************************************************** //
// ****************************************** Determination of fusion potential *************************************************** //
// ******************************************************************************************************************************** //

double Potential::FusPot(double r,int Apro,int Zpro,int Atar,int Ztar,double jc) {

	double mu = (Atar*Apro*1.)/(1.*(Atar+Apro))*getAmu();
	double muhb2 = 2.*mu/(getHbarc()*getHbarc());
    double aux = Fusion(Atar,Ztar,Apro,Zpro,r)+getCoul()*Ztar*Zpro/r+Vcentri(r,jc,muhb2);

	return(aux);
}

// ******************************************************************************************************************************** //
// ***************************** Determination of real part of optical potential for IWBC ***************************************** //
// ******************************************************************************************************************************** //

double Potential::Vtot(Readinput* path,double r,int Ares,int Zres,int Aemi,int Zemi,double jb,double lb,double elab) {
	
	double mu = (Ares*Aemi*1.)/(1.*(Ares+Aemi))*getAmu();
	double muhb2 = 2.*mu/(getHbarc()*getHbarc());
	double aux = 0.;

	if (Zemi==2) aux = V2_avr(r,Ares,Aemi,elab);
	else {
		if (path->getOptPot()==0) aux = V2_var(r,Ares,Zres,Zemi,jb,lb,elab);
		if (path->getOptPot()==1) aux = V2_kon(r,Ares,Zres,Zemi,jb,lb,elab);		
	}
	    
	aux = aux+Vcoul(path,r,Zres,Zemi,Ares,Aemi)+Vcentri(r,lb,muhb2);
	
	return(aux);
}

double Potential::dVtotdr(Readinput* path,double r,int Ares,int Zres,int Aemi,int Zemi,double jb,double lb,double elab) {

	double mu = (Ares*Aemi*1.)/(1.*(Ares+Aemi))*getAmu();
	double muhb2 = 2.*mu/(getHbarc()*getHbarc());
	double aux = 0.;

	if (Zemi==2) aux = dV2_avrdr(r,Ares,Aemi,elab);
	else {
		if (path->getOptPot()==0) aux = dV2_vardr(r,Ares,Zres,Zemi,jb,lb,elab);
		if (path->getOptPot()==1) aux = dV2_kondr(r,Ares,Zres,Zemi,jb,lb,elab);
	}
	    
	aux = aux+dVcouldr(path,r,Zres,Zemi,Ares,Aemi)+dVcentridr(r,lb,muhb2);
	
	return(aux);
}

// ******************************************************************************************************************************** //
// ************************************* Determination of optical potential ******************************************************* //
// ******************************************************************************************************************************** //

complex<double> Potential::OpticalPot(Readinput* path,double r,int Ares,int Zres,int Aemi,int Zemi,double jb,double lb,double elab) {
  
	double mu = (Ares*Aemi*1.)/(1.*(Ares+Aemi))*getAmu();
	double muhb2 = 2.*mu/(getHbarc()*getHbarc());

	double img = 0.;
	double rea = 0.;

	if (Zemi==2) rea = V2_avr(r,Ares,Aemi,elab);
	else {
		if (path->getOptPot()==0) rea = V2_var(r,Ares,Zres,Zemi,jb,lb,elab);
		if (path->getOptPot()==1) rea = V2_kon(r,Ares,Zres,Zemi,jb,lb,elab);		
	}
	    
	rea = rea+Vcoul(path,r,Zres,Zemi,Ares,Aemi)+Vcentri(r,lb,muhb2);
    
	if (Zemi==2) img = W2_avr(r,Ares,Aemi,elab);
	else {
		if (path->getOptPot()==0) img = W2_var(r,Ares,Zres,Zemi,jb,lb,elab);
		if (path->getOptPot()==1) img = W2_kon(r,Ares,Zres,Zemi,jb,lb,elab);
	}
	
	complex<double> Uopt(rea,img);

	return(Uopt);
}

// ******************************************************************************************************************************** //
// ***************************** Determination of Fermi function as well as its derivatives *************************************** //
// ******************************************************************************************************************************** //

double Potential::frra(double r,double Rad,double a) {

	double aux=0.;
  
	if (a==0) {
		if ((r-Rad)>=0) aux = 0.;
		if ((r-Rad)<0) aux = 1.;
	} else aux = 1./(1.+exp((r-Rad)/a));

	return(aux);
}

double Potential::dfrradr(double r,double Rad,double a) {

	double aux = 0.;
	double Fact = 1./(1.+exp((r-Rad)/a));
	
	if (a==0) {
		if ((r-Rad)>=0) aux =0.;
		if ((r-Rad)<0) aux =1.;
	} else aux = -1.*(1./a)*exp((r-Rad)/a)*Fact*Fact;

	return(aux); 
}

double Potential::ddfrradrdr(double r,double Rad,double a) {

	double aux = 0.;
	double Fact = 1./(1.+exp((r-Rad)/a));
 
	if (a==0) {
		if ((r-Rad)>=0) aux = 0.;
		if ((r-Rad)<0) aux = 1.;
	} else aux = 2.*(1./a/a)*exp(2.*(r-Rad)/a)*Fact*Fact*Fact-1.*(1./a/a)*exp((r-Rad)/a)*Fact*Fact;

	return(aux); 
}

// ******************************************************************************************************************************** //
// ***************************** Determination of Avrigeanu potential for alpha particles ***************************************** //
// ******************************************************************************************************************************** //

double Potential::V2_avr(double r,int Ares,int Zres,double elab) {
  
	double a0 = 101.1;
	double a1 = 6.051;
	double a2 = -0.248;
	double c0 = 0.817;
	double c1 = -0.0085;
	double rv = 1.245;
	double av = c0+c1*cbrt(Ares*1.);
	double aux=-(a0+a1*Zres/cbrt(Ares*1.)+a2*elab)*frra(r,rv,av);

	return(aux);
}

double Potential::dV2_avrdr(double r,int Ares,int Zres,double elab) {

	double a0 = 101.1;
	double a1 = 6.051;
	double a2 = -0.248;
	double c0 = 0.817;
	double c1 = -0.0085;
	double rv = 1.245;
	double av = c0+c1*cbrt(Ares*1.);
	double aux = -(a0+a1*Zres/cbrt(Ares*1.)+a2*elab)*dfrradr(r,rv,av);

	return(aux);
}

double Potential::W2_avr(double r,int Ares,int Zres,double elab) {
  
	double b0 = 26.82;
	double b1 = -1.706;
	double b2 = 0.006;
	double d0 = 0.692;
	double d1 = 0.020;
	double rw = 1.570;
	double aw = d0+d1*cbrt(Ares*1.);

	if (elab<73.) {
		b0 = 12.64;
		b1 = -1.706;
		b2 = 0.20;
	}

	double aux = -(b0+b1*cbrt(Ares*1.)+b2*elab)*frra(r,rw,aw);

	return(aux);
}

// ******************************************************************************************************************************** //
// ***************************** Determination of Varner potential for neutrons and protons *************************************** //
// ******************************************************************************************************************************** //

double Potential::V2_var(double r,int Ares,int Zres,int Zemi,double j,double l,double elab) {
    
	double Ro  = 1.25*cbrt(Ares*1.)-0.225;
	double rso = 1.34*cbrt(Ares*1.)-1.2;
  
	double Ec = 0.;
	double Rc = 0.;
	double Fact = 1.;
	
	if (Zemi==1) {
		Rc = 1.238*cbrt(Ares*1.)+0.116;
		Ec = 6.*getCoul()*Zres/(5.*Rc);
		Fact = -1.;
	} else {
		Ec = 0.;
		Fact = 1.;
	}

	double vr = 52.9 - Fact*13.1*(Ares-2.*Zres)/Ares-(elab-Ec)*0.299;
	double vso = 5.9;
	double SpinOrbit = j*(j+1.)-l*(l+1.)-0.5*(0.5+1.);

	double aux = -vr*frra(r,Ro,0.69)+2.*vso*SpinOrbit*(1./r)*dfrradr(r,rso,0.63);

	return(aux);
}

double Potential::dV2_vardr(double r,int Ares,int Zres,int Zemi,double j,double l,double elab) {
  
	double Ro  = 1.25*cbrt(Ares*1.)-0.225;
	double rso = 1.34*cbrt(Ares*1.)-1.2;

	double Ec = 0.;
	double Rc = 0.;
	double Fact=1.;
	  
	if (Zemi==1) {
		Rc = 1.238*cbrt(Ares*1.)+0.116;
		Ec = 6.*getCoul()*Zres/(5.*Rc);
		Fact = -1.;
	} else {
		Ec = 0.;
		Fact = 1.;
	}

	double vr = 52.9 - Fact*13.1*(Ares-2.*Zres)/Ares-(elab-Ec)*0.299;
	double vso = 5.9;
	double SpinOrbit = j*(j+1.)-l*(l+1.)-0.5*(0.5+1.);

	double aux = -vr*dfrradr(r,Ro,0.69)-2.*vso*SpinOrbit*(1./(r*r))*dfrradr(r,rso,0.63)+2.*vso*SpinOrbit*(1./r)*ddfrradrdr(r,rso,0.63);

	return(aux);
}

double Potential::W2_var(double r,int Ares,int Zres,int Zemi,double j, double l,double elab) {
  
	double rw  = 1.33*cbrt(Ares*1.)-0.42;
	double aw  = 0.69;
	double rso = 1.34*cbrt(Ares*1.)-1.2;
   
	double Ec = 0.;
	double Rc = 0.; 
	double Fact=1.;
	
	if (Zemi==1) {
		Rc = 1.238*cbrt(Ares*1.)+0.116;
		Ec = 6.*getCoul()*Zres/(5.*Rc);
		Fact = -1.;
	} else {
		Ec = 0.;
		Fact = 1.;
	}

	double wv = 7.8/(1.+exp((35.-(elab-Ec))/16.));
	double wsf = 4.*aw*(10.-Fact*18.*(Ares-2.*Zres)/Ares)/(1.+exp(((elab-Ec)-36)/37.));
	double wso = 10.;
	double SpinOrbit = j*(j+1.)-l*(l+1.)-0.5*(0.5+1.);
	
	double aux = -wv*frra(r,rw,aw)+wsf*dfrradr(r,rw,aw)+wso*2.*SpinOrbit*(1./r)*dfrradr(r,rso,0.63);
	
	return(aux);
}

// ******************************************************************************************************************************** //
// ***************************** Determination of Koning potential for neutrons and protons *************************************** //
// ******************************************************************************************************************************** //

double Potential::V2_kon(double r,int Ares,int Zres,int Zemi,double j,double l,double elab) {
    
    double aux = 0.;
    
    if (Zemi==0) {	    
		double vn1 = 59.30 - 21.0*(Ares - 2.*Zres)/Ares - 0.024*Ares;
		double vn2 = 0.007228 - 1.48e-6*Ares; 
		double vn3 = 1.994e-5 - 2.0e-8*Ares;
		double vn4 = 7.e-9;
		double vnso1 = 5.922 + 0.0030*Ares;
		double vnso2 = 0.0040;
		double Enf = -11.2814 + 0.02646*Ares;
		double vv = vn1*(1.-vn2*(elab-Enf)+vn3*(elab-Enf)*(elab-Enf)-vn4*(elab-Enf)*(elab-Enf)*(elab-Enf));
		double vso = vnso1*exp(-vnso2*(elab-Enf));	
		double rv = 1.3039 - 0.4054/cbrt(Ares*1.);
		double av = 0.6778 - 1.487e-4*Ares;
		double rso = 1.1854 - 0.647/cbrt(Ares*1.);
		double aso = 0.59;
		double SpinOrbit = j*(j+1.)-l*(l+1.)-0.5*(0.5+1.);
	
		aux = -vv*frra(r,rv,av)+2.*vso*SpinOrbit*(1./r)*dfrradr(r,rso,aso);
    } else {
		double vp1 = 59.30 + 21.0*(Ares - 2.*Zres)/Ares - 0.024*Ares;
		double vp2 = 0.007067 + 4.23e-6*Ares;
		double vp3 = 1.729e-5 + 1.136e-8*Ares;
		double vp4 = 7.e-9;
		double vpso1 = 5.922 + 0.0030*Ares;
		double vpso2 = 0.0040;
		double Epf = -8.4075 + 0.01378*Ares;
		double rc = 1.198 + 0.697/cbrt(Ares*Ares*1.) 
				+ 12.994/cbrt(double(Ares)*double(Ares)*double(Ares)*double(Ares)*double(Ares));
		double vc = 1.73*Zres/(rc*cbrt(Ares*1.));
		double vv = vp1*(1.-vp2*(elab-Epf)+vp3*(elab-Epf)*(elab-Epf)-vp4*(elab-Epf)*(elab-Epf)*(elab-Epf))
					+vc*vp1*(vp2-2.*vp3*(elab-Epf)+3.*vp4*(elab-Epf)*(elab-Epf));
		double vso = vpso1*exp(-vpso2*(elab-Epf));	
		double rv = 1.3039 - 0.4054/cbrt(Ares*1.);
		double av = 0.6778 - 1.487e-4*Ares;
		double rso = 1.1854 - 0.647/cbrt(Ares*1.);
		double aso = 0.59;
		double SpinOrbit = j*(j+1.)-l*(l+1.)-0.5*(0.5+1.);
	
		aux = -vv*frra(r,rv,av)+2.*vso*SpinOrbit*(1./r)*dfrradr(r,rso,aso);
    }
    
	return(aux);		
}

double Potential::dV2_kondr(double r,int Ares,int Zres,int Zemi,double j,double l,double elab) {
  
	double aux = 0.;
    
    if (Zemi==0) {	    
		double vn1 = 59.30 - 21.0*(Ares - 2.*Zres)/Ares - 0.024*Ares;
		double vn2 = 0.007228 - 1.48e-6*Ares; 
		double vn3 = 1.994e-5 - 2.0e-8*Ares;
		double vn4 = 7.e-9;
		double vnso1 = 5.922 + 0.0030*Ares;
		double vnso2 = 0.0040;
		double Enf = -11.2814 + 0.02646*Ares;
		double vv = vn1*(1.-vn2*(elab-Enf)+vn3*(elab-Enf)*(elab-Enf)-vn4*(elab-Enf)*(elab-Enf)*(elab-Enf));
		double vso = vnso1*exp(-vnso2*(elab-Enf));	
		double rv = 1.3039 - 0.4054/cbrt(Ares*1.);
		double av = 0.6778 - 1.487e-4*Ares;
		double rso = 1.1854 - 0.647/cbrt(Ares*1.);
		double aso = 0.59;
		double SpinOrbit = j*(j+1.)-l*(l+1.)-0.5*(0.5+1.);
	
    	aux = -vv*dfrradr(r,rv,av)-2.*vso*SpinOrbit*(1./(r*r))*dfrradr(r,rso,aso)+2.*vso*SpinOrbit*(1./r)*ddfrradrdr(r,rso,aso);
    } else {
		double vp1 = 59.30 + 21.0*(Ares - 2.*Zres)/Ares - 0.024*Ares;
		double vp2 = 0.007067 + 4.23e-6*Ares;
		double vp3 = 1.729e-5 + 1.136e-8*Ares;
		double vp4 = 7.e-9;
		double vpso1 = 5.922 + 0.0030*Ares;
		double vpso2 = 0.0040;
		double Epf = -8.4075 + 0.01378*Ares;
		double rc = 1.198 + 0.697/cbrt(Ares*Ares*1.) 
				+ 12.994/cbrt(double(Ares)*double(Ares)*double(Ares)*double(Ares)*double(Ares));
		double vc = 1.73*Zres/(rc*cbrt(Ares*1.));
		double vv = vp1*(1.-vp2*(elab-Epf)+vp3*(elab-Epf)*(elab-Epf)-vp4*(elab-Epf)*(elab-Epf)*(elab-Epf))
					+vc*vp1*(vp2-2.*vp3*(elab-Epf)+3.*vp4*(elab-Epf)*(elab-Epf));
		double vso = vpso1*exp(-vpso2*(elab-Epf));	
		double rv = 1.3039 - 0.4054/cbrt(Ares*1.);
		double av = 0.6778 - 1.487e-4*Ares;
		double rso = 1.1854 - 0.647/cbrt(Ares*1.);
		double aso = 0.59;
		double SpinOrbit = j*(j+1.)-l*(l+1.)-0.5*(0.5+1.);
	
    	aux = -vv*dfrradr(r,rv,av)-2.*vso*SpinOrbit*(1./(r*r))*dfrradr(r,rso,aso)+2.*vso*SpinOrbit*(1./r)*ddfrradrdr(r,rso,aso);
    }

	return(aux);
}

double Potential::W2_kon(double r,int Ares,int Zres,int Zemi,double j, double l,double elab) {
  
	double aux = 0.;
    
    if (Zemi==0) {	    
		double wn1 = 12.195 + 0.0167*Ares;
		double wn2 = 73.55 + 0.0795*Ares;
		double dn1 = 16.0 - 16.0*(Ares - 2.*Zres)/Ares;
		double dn2 = 0.0180 + 0.003802/(1.0 + exp((Ares - 156.)/8.));
		double dn3 = 11.5;
		double wnso1 = -3.1;
		double wnso2 = 160.;
		double Enf = -11.2814 + 0.02646*Ares;
		double wv = wn1*(elab-Enf)*(elab-Enf)/((elab-Enf)*(elab-Enf)+wn2*wn2);
		double wd = dn1*(elab-Enf)*(elab-Enf)/((elab-Enf)*(elab-Enf)+dn3*dn3)*exp(-dn2*(elab-Enf));
		double wso = wnso1*(elab-Enf)*(elab-Enf)/((elab-Enf)*(elab-Enf)+wnso2*wnso2);
		double rv = 1.3039 - 0.4054/cbrt(Ares*1.);
		double av = 0.6778 - 1.487e-4*Ares;
		double rd = 1.3424 - 0.01585*cbrt(Ares*1.);
		double ad = 0.5446 - 1.656e-4*Ares;
		double wsf = 4.*ad*wd;
		double rso = 1.1854 - 0.647/cbrt(Ares*1.);
		double aso = 0.59;
		double SpinOrbit = j*(j+1.)-l*(l+1.)-0.5*(0.5+1.);
	
		aux = -wv*frra(r,rv,av)+wsf*dfrradr(r,rd,ad)+wso*2.*SpinOrbit*(1./r)*dfrradr(r,rso,aso);
    } else {
		double wp1 = 14.667 + 0.009629*Ares;
		double wp2 = 73.55 + 0.0795*Ares;
		double dp1 = 16.0 + 16.0*(Ares - 2.*Zres)/Ares;
		double dp2 = 0.0180 + 0.003802/(1.0 + exp((Ares - 156.)/8.));
		double dp3 = 11.5;
		double wpso1 = -3.1;
		double wpso2 = 160.;
		double Epf = -8.4075 + 0.01378*Ares;
		double wv = wp1*(elab-Epf)*(elab-Epf)/((elab-Epf)*(elab-Epf)+wp2*wp2);
		double wd = dp1*(elab-Epf)*(elab-Epf)/((elab-Epf)*(elab-Epf)+dp3*dp3)*exp(-dp2*(elab-Epf));
		double wso = wpso1*(elab-Epf)*(elab-Epf)/((elab-Epf)*(elab-Epf)+wpso2*wpso2);		
		double rv = 1.3039 - 0.4054/cbrt(Ares*1.);
		double av = 0.6778 - 1.487e-4*Ares;
		double rd = 1.3424 - 0.01585*cbrt(Ares*1.);
		double ad = 0.5187 + 5.205e-4*Ares;
		double wsf = 4.*ad*wd;
		double rso = 1.1854 - 0.647/cbrt(Ares*1.);
		double aso = 0.59;
		double SpinOrbit = j*(j+1.)-l*(l+1.)-0.5*(0.5+1.);
	
		aux = -wv*frra(r,rv,av)+wsf*dfrradr(r,rd,ad)+wso*2.*SpinOrbit*(1./r)*dfrradr(r,rso,aso);
    }
    
	return(aux);
}

// ******************************************************************************************************************************** //
// **************************************** Determination of Coulomb potential **************************************************** //
// ******************************************************************************************************************************** //

double Potential::Vcoul(Readinput* path,double r,int Zres,int Zemi,int Ares,int Aemi) {
  
	double aux = 0.;
	double rc = 0.;
	
	if (Zemi==0) return(0.);
	else if (Zemi==1) {
		if (path->getOptPot()==0) 
			rc = 1.238*cbrt(Ares*1.) + 0.116;;
		if (path->getOptPot()==1) 
			rc = 1.198*cbrt(Ares*1.) + 0.697/cbrt(Ares*1.) + 12.994/cbrt(double(Ares)*double(Ares)*double(Ares)*double(Ares));
	}
	else if (Zemi==2) rc = 1.2*(cbrt(Ares*1.) + cbrt(Aemi*1.));
	else {
		cerr << "Warning: Problem with the Coulomb radius." << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
	
	if (r<rc) aux = 1.5*getCoul()*Zres*Zemi*(1.-(r/rc)*(r/rc)/3.)/rc;
	if (r>=rc) aux = getCoul()*Zres*Zemi/r;
	       
	return(aux);
}

double Potential::dVcouldr(Readinput* path,double r,int Zres,int Zemi,int Ares,int Aemi) {
	
	double aux = 0.;
	double rc = 0.;
	
	if (Zemi==0) return(0.);
	else if (Zemi==1) {
		if (path->getOptPot()==0) 
			rc = 1.238*cbrt(Ares*1.) + 0.116;;
		if (path->getOptPot()==1) 
			rc = 1.198*cbrt(Ares*1.) + 0.697/cbrt(Ares*1.) + 12.994/cbrt(double(Ares)*double(Ares)*double(Ares)*double(Ares));
	}
	else if (Zemi==2) rc = 1.2*(cbrt(Ares*1.) + cbrt(Aemi*1.));
	else {
		cerr << "Warning: Problem with the Coulomb radius." << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
  
	if (r<rc) aux = 1.5*getCoul()*Zres*Zemi*(-2.*r/(3.*rc*rc*rc));
	if (r>=rc) aux = -getCoul()*Zres*Zemi/(r*r);
	       
	return(aux);
}

// ******************************************************************************************************************************** //
// **************************************** Determination of centrifuge potential ************************************************* //
// ******************************************************************************************************************************** //

double Potential::Vcentri(double r,double l,double muhb2) {

	return(l*(l+1.)/(muhb2*r*r));
}


double Potential::dVcentridr(double r,double l,double muhb2) {

	return(-2.*l*(l+1.)/(muhb2*r*r*r));
}

double Potential::Fusion(int Apro,int Zpro,int Atar,int Ztar,double r) {

	double At = 1.126*cbrt(Atar*1.);
	double Ap = 1.126*cbrt(Apro*1.);
	double S0 = r-At-Ap;
	double Zeta = S0/0.79;
	double Phi = 0.;
  
	if (Zeta<1.2511) {
		double Dum = Zeta-2.54;
		Phi = -0.5*Dum*Dum-0.0852*Dum*Dum*Dum;
	} else Phi = -3.437*exp(-Zeta/0.75);

	double aux = 59.*0.79*At*Ap/(Ap*1.+At*1.)*Phi;

/*	Myers and Swiatecki, PRC 62, 044610 (2000)

	double b=1.;
	double r0=1.14;
	double J=32.65;
	double c1=0.757895;
	double Q=35.4;
	double It=(Atar-2.*Ztar)/Atar, Ip=(Apro-2.*Zpro)/Apro;
	
	double dt=1.5*r0*(J*It-c1*Ztar*pow(1.*Atar,-1./3.)/12.)/(Q+2.25*J*pow(1.*Atar,-1./3.));
	double dp=1.5*r0*(J*Ip-c1*Zpro*pow(1.*Apro,-1./3.)/12.)/(Q+2.25*J*pow(1.*Apro,-1./3.));

	double Rt=1.24*pow(1.*Atar,1./3.)*(1.+1.646/Atar-0.191*It);
	Rt=Rt*(1.-7./(2*Rt*Rt)-49./(8.*Rt*Rt*Rt*Rt))+1.*(Atar-Ztar)*dt/Atar;
	double Rp=1.24*pow(1.*Apro,1./3.)*(1.+1.646/Apro-0.191*Ip);	
	Rp=Rp*(1.-7./(2*Rp*Rp)-49./(8.*Rp*Rp*Rp*Rp))+1.*(Apro-Zpro)*dp/Apro;
	
	double Req=Rt*Rp/(Rt+Rp);
	double gamma=(18.36-0.5*Q*(dt*dt+dp*dp)/(r0*r0))/(4.*getPi()*r0*r0);
	
	double Zeta = r-Rt-Rp;
	double Phi=0.;	
	
	double c[6];
	c[0]=-0.1886;
	c[1]=-0.2628;
	c[2]=-0.15216;
	c[3]=-0.04562;
	c[4]=0.069136;
	c[5]=-0.011454; 
	
	if (Zeta<2.5 && Zeta>0.) {
		double Dum = 2.5-Zeta;
		Phi=-0.1353;
		for (int i=0;i<6;i++) Phi+=c[i]*pow(Dum,i+1.)/(i+1.);
	} else Phi=-0.09551*exp((2.75-Zeta)/0.7176);
	
	double aux = 4.*getPi()*b*gamma*Req*Phi;
*/	
	return(aux);
}

Potential::Potential() {
	Rmin = 0.;	
	Rcoul = 0.;	
	Bmin = 0.;	
	Bcoul = 0.;
}

Potential::~Potential() {
}
