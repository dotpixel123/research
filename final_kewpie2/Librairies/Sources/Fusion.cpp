#include <Fusion.h>

double Fusion::getSigfus(int i) {
	return(SigFus[i]);
}

void Fusion::Parabo(double* X,double* Y,double* A,double& X0,double& Y0) { 

	// Local Parabolic Approximation: Finds A[] of V(R)=A[2]R^2+A[1]R+A[0]
	// for known three points (X,Y) or (R,V)

	X0 = 0.;
	Y0 = 0.;
	double aux = X[1]-X[0];
	if (aux==0.) return;
	double F01 = (Y[1]-Y[0])/aux;
	aux = X[2]-X[0];
	if (aux==0.) return;
	double F02 = (Y[2]-Y[0])/aux;
	aux = X[1]-X[2];
	if (aux==0.) return;
	double F012 = (F01-F02)/aux;
	A[1] = F01-F012*(X[0]+X[1]);
	A[2] = F012;
	A[0] = Y[0]-X[0]*(A[1]+X[0]*F012);
	if (A[2]==0) return;
	// (X0,Y0) is the maximum of the parabola
	X0 = -0.5*A[1]/A[2];	
	Y0 = -0.25*A[1]*A[1]/A[2]+A[0];
}

double Fusion::TWKB(int Ap,int Zp,int At,int Zt,double J,double Ecm,double Rfus,double Dr,double& Rin,double& Rout) {

	double Fac = Ap*At*1./(Ap+At);
	double F0 = sqrt(2.*Fac*getAmu()/(getHbarc()*getHbarc()));

	int Iup = int(4./Dr);
	int Iupp = int(8./Dr); 

	double R1 = 0.;
	double V1 = 0.;
	double R = 0.;
	double V = 0.;
	int i = 0;

// ******************************************************************************************************** //	
// ************************************** Determination of Rin ******************************************** //	
// ******************************************************************************************************** //	

	R = Rin;
	V = Pot->FusPot(R,Ap,Zp,At,Zt,J);
	i = -1;
	bool test1 = false;	
  	
	do {
		i++;
		R1 = R;
		V1 = V;
		R = R-Dr;
		V = Pot->FusPot(R,Ap,Zp,At,Zt,J);
		if (V<Ecm) test1 = true;
	} while (i<Iup && !test1);

	if (!test1) {
		if (J<3) 
			return(0.);
	} else {
		double R3 = 0.5*(R+R1);
		double V2 = Pot->FusPot(R3,Ap,Zp,At,Zt,J);
		if (V2<Ecm) {
			V = V2;
			R = R3;
		} else {
			V1 = V2;
			R1 = R3;
		}
		Rin = ((V1-Ecm)*R-(V-Ecm)*R1)/(V1-V); // position of the inner turning point   
	}
	
// ******************************************************************************************************** //	
// ************************************** Determination of Rout ******************************************** //	
// ******************************************************************************************************** //

	R = Rout;
	V = Pot->FusPot(Rout,Ap,Zp,At,Zt,J);
	i = -1;
	bool test2 = false;

	do {
		i++;
		R1 = R;
		V1 = V;
		R = R+Dr;
		V = Pot->FusPot(R,Ap,Zp,At,Zt,J);
		if (V<Ecm) test2 = true;
	} while (i<Iupp && !test2);

	if (!test2) {
		if (J<3)
			return(0.);
	} else {
		double R3 = 0.5*(R+R1);
		double V2 = Pot->FusPot(R3,Ap,Zp,At,Zt,J);
		if (V2<Ecm) {
			V = V2;
			R = R3;
		} else {
			V1 = V2;
			R1 = R3;
		}
		Rout = ((V1-Ecm)*R-(V-Ecm)*R1)/(V1-V); // position of the outer turning point
		if (Rfus<Rout) Rin = max(Rin,Rfus);
	}
	
// ******************************************************************************************************** //	
// *********************************** Calculation of integral ******************************************** //	
// ******************************************************************************************************** //

	double C = (Rout-Rin)/2.;
	double D = (Rout+Rin)/2.;
	double sum = 0.;
	for (size_t j=0;j<Npoint;j++) {
		double Xj = C*point[j] + D;
		double aux = Pot->FusPot(Xj,Ap,Zp,At,Zt,J)-Ecm;
		if (aux>0.) aux = sqrt(aux);
		else aux = 0.;
		sum += weight[j]*aux;
	}
	
	double T = exp(-2.*C*sum*F0);
	return(T);
}

void Fusion::WKB_Fusion_J(int Ap,int Zp,int At,int Zt,double Ecm,double R) {

	double Fac = Ap*At*1./(Ap+At);
	double HW0 = getHbarc()*getHbarc()/(getAmu()*Fac);
	double SIG0 = 15.708*HW0/Ecm;
	int Afus = Ap+At;

	double Dr = 0.1;
	double MaxExp = 100.;
	double Tmin = exp(-MaxExp);
	double V = Pot->FusPot(R,Ap,Zp,At,Zt,0.);
  
	R += Dr;
	double Jaux = 0.;	
	if (Afus%2==1) Jaux=0.5;

	double V1 = Pot->FusPot(R,Ap,Zp,At,Zt,Jaux);
	double V2 = 0.;
	int i = -1;
	bool test = false;
	do {
		i++;
		R += Dr;
		V2 = Pot->FusPot(R,Ap,Zp,At,Zt,Jaux);
		if (V2<V1 && V<V1) test = true;
		else {
			V = V1;
			V1 = V2;
		}
	} while (i<100 && !test);

	if (test) {

		double* X = new double[3];
		double* Y = new double[3];
		double* A = new double[3];
		
		for (int i=0;i<3;i++) {
			X[i] = 0.;
			Y[i] = 0.;
			A[i] = 0.;
		}
		
		double Rmax = 0.;
		double Vlb = 0.;
    
		X[2] = R;
		Y[2] = V2;
		X[1] = X[2]-Dr;
		Y[1] = V1;
		X[0] = X[1]-Dr;
		Y[0] = V;
    
		Parabo(X,Y,A,Rmax,Vlb);
		
		int J = -1;
		double T = 0.;
		double Tcut = 0.;
		double Rin = 0.;		
		double Rout = 0.;
		double Rinn = 0.;
		double Routt = 0.;
		double RMAX = 0.;
		double Tsave = 0.;

		do {
      
			J++;
			// For nuclei with odd A, J is increased by 1/2.
			if (Afus%2==1) J = J+0.5;
			
			if (J!=0) Rmax = RMAX;

			X[2] = Rmax+0.15;
			Y[2] = Pot->FusPot(X[2],Ap,Zp,At,Zt,J);
			X[1] = Rmax;
			Y[1] = Pot->FusPot(X[1],Ap,Zp,At,Zt,J);
			X[0] = Rmax-0.15;
			Y[0] = Pot->FusPot(X[0],Ap,Zp,At,Zt,J);
   
			Parabo(X,Y,A,Rmax,Vlb);

			bool test3 = false;

			if (J!=0) {
				if (Rmax>X[2]) {
					T = 0.;
					test3 = true;
				} else {
					if (Rmax<X[0]) {
						double Raux = Rmax-0.15;
						do {
							X[2] = X[1];
							Y[2] = Y[1];
							X[1] = X[0];
							Y[1] = Y[0];
							Raux = Rmax-0.15;  
							X[0] = Raux;
							Y[0] = Pot->FusPot(X[0],Ap,Zp,At,Zt,J);
							Parabo(X,Y,A,Rmax,Vlb);
						} while (Raux>0 && Rmax<X[2] && Rmax<X[0]);
						if (Rmax>X[2] || Raux==0.) {
							T = 0.;
							test3 = true;
						}
					}
				}
      		}
      		
			T = 0.;

			if (!test3) {
				double aux1 = 6.28319*(Vlb-Ecm)/sqrt(HW0*abs(2.*A[2]));
				if (A[2]<0. || J==0) {
					if (abs(aux1)<MaxExp) T = 1./(1.+exp(aux1));
					if (aux1>MaxExp) T = 0.;
					if (aux1<-MaxExp) T = 1.;
					if (J==0 || T>0.) Tsave = T;
				} else T = 0.;
			}
     
			RMAX = Rmax;

			if (T<0.49)
				if (T!=0.) {
					double Dr2 = 0.1;
					double aux = 0.;
					if (J==0) {
						Rin = Rmax;
						Rout = 1.44*Zp*Zt/Ecm-2.;
						if (Zp==0) Rout = sqrt(HW0/Ecm)-3.;
						Rout = max(Rmax,Rout);
	
						aux = TWKB(Ap,Zp,At,Zt,J,Ecm,0.,Dr2,Rin,Rout);
	    
						T = 0.;
						if (aux>0.) T = 1./(1.+(1./aux));
						Rinn = Rin;
						Routt = Rout;
					} else {
						if (Routt>0.) {
							if (Rinn>0.) {
								aux = (Routt-Rinn)*0.1;
								Dr = min(0.2,aux);
								Dr = max(0.02,Dr);
								Rin = Rinn+Dr;
								Rin = min(RMAX,Rin);
								Rout = Routt-Dr;
								Rout = max(Rout,RMAX);
							} else {
								Rin = Rmax;
								Rout = Rmax;
							}
						} else {
							Rin = Rmax;
							Rout = Rmax;
						}
					
						double Rfus = 0.;
						if (Zp==0) Rout = sqrt(HW0*J*(J+1)/(2.*Ecm))-0.5;
						else Rfus = -Rinn;
						
						aux = TWKB(Ap,Zp,At,Zt,J,Ecm,Rfus,Dr2,Rin,Rout);
	    
						if (aux>0.) T = 1./(1.+(1./aux));
						if (aux<=0. && Tsave>0.05) T = Tsave;
						Rinn = Rin;
						Routt = Rout;
					}
				}
			
			if (J==0) Tcut = 1.e-6*T;      
			if (T>Tcut) SigFus[J] = T*SIG0*(2.*J+1);
		} while (J<999 && T>Tmin && T>Tcut);

		delete[] X;
		delete[] Y;
		delete[] A;
	}
}

double Fusion::EBD_Fusion(int Ap,int Zp,int At,int Zt,double ecm,double B0,double Bvar,double R) {

	double Xfuse = (ecm-B0)/(Bvar*getSqrtTwo());
	double sigmax = 10.*getPi()*R*R*Bvar/(ecm*getSqrtTwo())*(Xfuse*(1.+erf(Xfuse))+exp(-Xfuse*Xfuse)*getInvSqrtPi()); // in mb
	
	return(sigmax);
}

void Fusion::EBD_Fusion_J(int Ap,int Zp,int At,int Zt,double ecm,double B0,double Bvar,double R) {

	double Fac = Ap*At*1./(Ap+At);
	double HW0 = getHbarc()*getHbarc()/(getAmu()*Fac);
	double muhb2= 2.*(getAmu()*Fac)/(getHbarc()*getHbarc());
	double SIG0 = 15.708*HW0/ecm; // in mb
	int Afus = Ap+At;
	
	int J = -1;
	do {
		J++; // For nuclei with odd A, J is increased by 1/2.
		double Jaux = 0.;
		if (Afus%2==1) 
			Jaux = J+0.5;
		double Beff = B0+Jaux*(Jaux+1.)/(muhb2*R*R); 
		double T = 0.5*(1.+erf((ecm-Beff)/(Bvar*getSqrtTwo())));
		SigFus[J] = SIG0*(2.*Jaux+1.)*T;
	} while (J<999); // number of partial waves equals 1000
}

Fusion::Fusion(size_t NumPoint) {

	weight = new double [NumPoint];
	point = new double [NumPoint];

// ************************************************************************************************************************* //
// ************************************** Initialization for Gauss-Legendre Quadrature ************************************* //
// ************************************************************************************************************************* //

	Npoint = NumPoint;
	
	gsl_integration_glfixed_table* tab;
	tab = gsl_integration_glfixed_table_alloc(NumPoint);
	size_t halfNPoint = tab->n/2;

	for (size_t j=0;j<halfNPoint;j++) {
		weight[j] = tab->w[j];
		point[j] = -tab->x[j];
		weight[j+halfNPoint] = weight[j];
		point[j+halfNPoint] = -point[j];
	}
	
	gsl_integration_glfixed_table_free(tab); // Free memory

	Pot = new Potential();
	SigFus = new double[1000];
	for (int i=0;i<1000;i++) SigFus[i] = 0.;
}

Fusion::~Fusion() {

	delete[] weight;
	delete[] point;
	delete[] SigFus;
	delete Pot;
	
	weight = NULL;
	point = NULL;
	SigFus = NULL;
	Pot = NULL;
}
