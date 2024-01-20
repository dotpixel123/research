#include <Element.h>

int Element::getA() {
	return(A);
}

void Element::setA(int A1) {
	A=A1;
}

int Element::getZ() {
	return(Z);
}

void Element::setZ(int Z1) {
	Z=Z1;
}

string Element::getName() { 
	return(Name);
}

double Element::getRadius() {
	return(Radius);
}

double Element::getCoulRadius() {
	return(CoulRadius);
}

double Element::getMass() {
	return(Mass);
}

double Element::getMassExcess() {
	return(MassExcess);
}

double Element::getx() {
	return(x);
}

double Element::getEsh() {
	return(Esh);
}

double Element::getFisBarrier() {
	return(Bf);
}

void Element::setEsh(double Esh1) {
	Esh=Esh1;
}

double Element::getInertias() {
	return(Inertias);
}

double Element::getInertiag() {
	return(Inertiag);
}

double Element::getPairing() {
	return(Pairing);
}

double Element::getBetag() {
	return(Betag);
}

double Element::getBetas() {
	return(Betas);
}

double Element::getPop() {
	return(Pop);
}

void Element::setPop(double Pop1) {
	Pop=Pop1;
}

int Element::getMaxt() {
	return(maxt);
}

void Element::setMaxt(int Maxt) {
	maxt = Maxt;
}

int Element::getMaxi() {
	return(maxi);
}

void Element::setMaxi(int Maxi) {
	maxi = Maxi;
}

int Element::getMaxj() {
	return(maxj);
}

void Element::setMaxj(int Maxj) {
	maxj = Maxj;
}

double Element::getSpectrum(int t,int i,int j) {
	return(Spectrum[t][i][j]);
}

void Element::setSpectrum(int t,int i,int j,double value) {
	Spectrum[t][i][j] = value;
}

void Element::ElementName(int Z1) {
  
	string mendel[140] = {"n  ","H  ","He ","Li ","Be ","B  ","C  ","N  ","O  ","F  ","Ne ","Na ","Mg ","Al ","Si ",
	"P  ","S  ","Cl ","Ar ","K  ","Ca ","Sc ","Ti ","V  ","Cr ","Mn ","Fe ","Co ","Ni ","Cu ","Zn ","Ga ","Ge ","As ",
	"Se ","Br ","Kr ","Rb ","Sr ","Y  ","Zr ","Nb ","Mo ","Tc ","Ru ","Rh ","Pd ","Ag ","Cd ","In ","Sn ","Sb ","Te ",
	"I  ","Xe ","Cs ","Ba ","La ","Ce ","Pr ","Nd ","Pm ","Sm ","Eu ","Gd ","Td ","Dy ","Ho ","Er ","Tm ","Yb ","Lu ",
	"Hf ","Ta ","W  ","Re ","Os ","Ir ","Pt ","Au ","Hg ","Tl ","Pb ","Bi ","Po ","At ","Rn ","Fr ","Ra ","Ac ","Th ",
	"Pa ","U  ","Np ","Pu ","Am ","Cm ","Bk ","Cf ","Es ","Fm ","Md ","No ","Lr ","Rf ","Db ","Sg ","Bh ","Hs ","Mt ",
	"Ds ","Rg ","Cn ","113 ","Fl ","115","Lv ","117 ","118 ","119 ","120 ","121 ","122 ","123 ","124 ","125 ","126 ",
	"127 ","128 ","129 ","130 ","131 ","132 ","133 ","134 ","135 ","136 ","137 ","138 ","139 "};
  
	Name = mendel[Z1];
}

double Element::Fissility(int A1,int Z1) {

	// *********** Calculation of the fissility ********************
	//        Reference K H Schmidt  W Moraveck (1991)  
	//          Rep. Prog. Phys. 54 (1991) 949-1003                
	// *************************************************************

	double I = (A1-2.*Z1)/(A1*1.);
	double aux = (Z1*Z1*1.)/(49.22*A1*(1.-0.3803*I*I-20.489*I*I*I*I));
	
	if (aux<0.) aux = 0.;
	if (aux>1.) aux = 1.;
	
	return(aux);
}

double Element::FisBar1(int A1,int Z1) {
 
  // ----- Calculation of the fission Barrier -----------------------------
  //           Myers et Swiatecki PRC 60 (1999) 014606
  // ----------------------------------------------------------------------

	double F = 0.;
	double X0 = 48.5428;
	double X1 = 34.15;
	double I = (A1*1.-2.*Z1)/A1;
	double k = 1.9+(Z1*1.-80.)/75.;
	double S = cbrt(A1*A1*1.)*(1.-k*I*I);
	double X = Z1*Z1*1./(A1*(1.-k*I*I));
	
	if (30.<=X && X<=X1) F = 0.595553-0.124136*(X-X1);
	if (X1<X && X<=X0) F = 0.000199749*(X0-X)*(X0-X)*(X0-X);
	
	double FisBar = F*S; // Fission barrier height in MeV
	
	if (FisBar<0.) FisBar = 0.;

	return(FisBar);
}

double Element::FisBar2(int A1,int Z1) {

  // ----- Calculation of the fission Barrier -----------------------------
  // M. Dahlinger, D. Vermeulen, and K.-H. Schmidt, Nucl. Phys. A 376, 94 (1982)                
  // ----------------------------------------------------------------------

/*  
	double u = 0.368-5.057*x+8.93*pow(x,2.)-8.71*pow(x,3.);
	double FisBar = 0.7322*pow(10.,u)*(pow((1.*Z),2.)/pow((1.*A),(1./3.)));
*/

  // ----- Calculation of the fission Barrier -----------------------------
  //  F. A. Ivanyuk and K. Pomorski, Phys. Rev. C 79, 054327 (2009).                
  // ----------------------------------------------------------------------

	double a0 = 723.255;
	double a1 = -19.60;
	double a2 = 17.70;
	double a3 = -5.328;
	double a4 = 0.2739;
	double a5 = -0.8633;
	double a6 = 0.09067;
	double a7 = 0.6911;
	double a8 = -0.661;

	if (Z<75) { 
		a0 = 23.92;
		a1 = 2.235;
		a2 = -3.816;
		a3 = 1.128;
		a4 = 0.2268;
		a5 = 4.340;
		a6 = 0.7851;
		a7 = -0.4141;
		a8 = -0.438;
	} 
	
	double Bmax = a0+a1*Z1+a2*Z1*Z1*1e-2+a3*double(Z1)*double(Z1)*double(Z1)*1e-4;
	double I = (A1*1.-2.*Z1)/A1;
	double I0 = a4+a5*Z1*1e-4;
	double DI = a6+a7*Z1*1e-2+a8*Z1*Z1*1e-4;
	double FisBar = Bmax*exp(-(I-I0)*(I-I0)/(DI*DI));
	
	if (FisBar<0.) FisBar = 0.;
	return(FisBar);  
}

double Element::Bfcal(int FisBar,int Asrc,int Zsrc,double DeltaShell) {

	double Bf = 0.;
	switch (FisBar) {
		case 1: 
			Bf = FisBar1(Asrc,Zsrc);
			break;
		case 2:
			Bf = FisBar2(Asrc,Zsrc);
			break;
		default: 
			cerr << "Warning: Please check the fission-barrier option in the input file." << endl;
    		cin.get();
			exit(EXIT_FAILURE);
	}
	return(Bf-DeltaShell);
}

double Element::InertiaMoment(int A1,double Beta2) {

	// *******************Calcultaion of I/hbar ************************
	//      reference  Bohr et Mottelson Book 2 equation 4-104
	//      V.I. Zagrebaev, et al., Phys. Rev., 2001, vol.C65, p.014607
	// *****************************************************************
 
	double inertia = 0.4*A1*getAmu()*getR0()*cbrt(A1*1.)*getR0()*cbrt(A1*1.)/(getHbarc()*getHbarc());
	
	inertia *= (1.+Beta2*sqrt(5./(16.*getPi()))+Beta2*Beta2*45./(28.*getPi())); // Taking the deformation into account
	
	return(inertia);
}

double Element::PairingEnergy(int A1,int Z1) {

	double Epair = 0.;

	if ((((A1-Z1)%2)==0)&&((Z1%2)==0)) Epair = 24./sqrt(A1*1.);
	if ((A1%2)==1) Epair = 12./sqrt(A1*1.);

	return(Epair); 
}

void Element::Nucleus(Readinput* path,int A1,int Z1,double MassArray,double ShellArray,double BetaArray,double FisBarArray) {
  
	A = A1;
	Z = Z1;
	
	ElementName(Z);
	Radius = getR0()*cbrt(A*1.); 
	CoulRadius = getR0Coul()*cbrt(A*1.);
	MassExcess = MassArray;

	Mass = Z*getBp()+(A-Z)*getBn()-MassArray+Z*getBe(1)-getBe(Z);
	if (MassArray==0.) {
		cerr << "Warning: No mass excess found for element Z=" << Z << " and A=" << A << ". Please Check it." << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
  
	x = 0.;  
	if (A>1) x = Fissility(A,Z);
  
	Pairing = 0.;  
	if (path->getPair()==1) Pairing = PairingEnergy(A,Z);

	Esh = 0.;
	if (path->getShell()==1) Esh = ShellArray*path->getShellFactor();
	
	if (path->getFisBar()==0) {
		Bf = FisBarArray;
		if (Bf<0.) Bf = 0.;
	} else Bf = Bfcal(path->getFisBar(),A1,Z1,Esh);	

	Betag = BetaArray;
	if (BetaArray==50.) {
		cerr << "Warning: No beta2 value found for element Z=" << Z << " A=" << A << ". Please Check it." << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}	
	double y = 1.-x;
	Betas = sqrt(4.*getPi()/5.)*(7.*y/3.-938.*y*y/765.+9.499768*y*y*y-8.050944*y*y*y*y);
	if (Betas<Betag) Betas = Betas + Betag;

	Pop = 0.;
   
	Inertias = InertiaMoment(A1,Betas)*path->getFacInerMom();
	Inertiag = InertiaMoment(A1,Betag)*path->getFacInerMom();
}

// ************************************************************************************************************* //
// ******************************** construction of the spectrum *********************************************** //
// ************************************************************************************************************* //

void Element::allocspectrum(int MaxT,int MaxI,int MaxJ) {

	maxt = MaxT; // number of bins for estimating fission time
	maxi = MaxI; // numer of bins for discretizing exitation energy
	maxj = MaxJ; // number of bins for discretizing angular momentum

	Spectrum = new double**[maxt];
	for (int t=0;t<maxt;t++) {
		Spectrum[t] = new double*[maxi];
		for (int i=0;i<maxi;i++) {
			Spectrum[t][i] = new double[maxj];
			for (int j=0;j<maxj;j++)
				Spectrum[t][i][j] = 0.;
		}
	}   
}

void Element::InitSpectrum() {
     
	for (int t=0;t<maxt;t++)
		for (int i=0;i<maxi;i++)
			for (int j=0;j<maxj;j++)
				Spectrum[t][i][j] = 0.;    
}

void Element::deallocspectrum() {

	if (maxt!=0 && maxi!=0 && maxj!=0) {
		for (int t=0;t<maxt;t++) {
			for (int i=0;i<maxi;i++) 
				delete[] Spectrum[t][i];
			delete[] Spectrum[t];
		}
		delete[] Spectrum;
	}

	maxt = 0;
	maxi = 0;
	maxj = 0;
	
	Spectrum = NULL;
}

void Element::InitElement() {

	A = 0;
	Z = 0;
	Mass = 0.;
	Name = " ";
	Radius = 0.;
	CoulRadius = 0.;
	x = 0.;
	Esh = 0.;
	Inertias = 0.;
	Inertiag = 0.;
	Pairing = 0.;
	Betag = 50.;
	Betas = 50.;
	Pop = 0.;
	Bf = 0.;
  
	deallocspectrum();
}

// ************************************************************************************************************* //
// ************************************************************************************************************* //
// ************************************************************************************************************* //

Element::Element() {

	maxi = 0;
	maxt = 0;
	maxj = 0;
	nbmax = 0;
	emax = 0;

	Spectrum = NULL;
	TabCoef = NULL;

	trans = new Transcoef();
}

Element &Element::operator+=(Element& source) {

	int Maxi = maxi;
	double aux[maxt][Maxi][maxj];
	
	for (int t=0;t<maxt;t++)
		for (int i=0;i<Maxi;i++) 
			for (int j=0;j<maxj;j++) 
				aux[t][i][j] = Spectrum[t][i][j];
	
	deallocspectrum();
	allocspectrum(source.getMaxt(),source.getMaxi(),source.getMaxj());
	
	for (int t=0;t<maxt;t++)
		for (int i=0;i<Maxi;i++)
			for (int j=0;j<maxj;j++) 
				Spectrum[t][i][j] = aux[t][i][j];
  
	if (maxi>=source.getMaxi()) Maxi = source.getMaxi();

	for (int t=0;t<maxt;t++)
		for (int i=0;i<Maxi;i++)
			for(int j=0;j<maxj;j++) 
				Spectrum[t][i][j] = Spectrum[t][i][j]+source.getSpectrum(t,i,j);

	return(*this);
}

Element &Element::operator=(Element& source) {

	A = source.getA();
	Z = source.getZ();
	Mass = source.getMass();
	Name = source.getName();
	Radius = source.getRadius();
	CoulRadius = source.getCoulRadius();
	x = source.getx();
	Esh = source.getEsh();
	Inertias = source.getInertias();
	Inertiag = source.getInertiag();
	Pairing = source.getPairing();
	Betag = source.getBetag();
	Betas = source.getBetas();
	Pop = source.getPop();
	Bf = source.getFisBarrier();
  	
	if (maxi!=source.getMaxi()||maxt!=source.getMaxt()||maxj!=source.getMaxj()) {
		deallocspectrum();
		maxi = source.getMaxi();
		maxt = source.getMaxt();
		maxj = source.getMaxj();
		allocspectrum(maxt,maxi,maxj);
	}
	
	for (int t=0;t<maxt;t++)
		for (int i=0;i<maxi;i++)
			for (int j=0;j<maxj;j++)
				Spectrum[t][i][j] = source.getSpectrum(t,i,j);

	return(*this);
}

Element &Element::operator=(Element*& source) {

	A = source->getA();
	Z = source->getZ();
	Mass = source->getMass();
	Name = source->getName();
	Radius = source->getRadius();
	CoulRadius = source->getCoulRadius();
	x = source->getx();
	Esh = source->getEsh();
	Inertias = source->getInertias();
	Inertiag = source->getInertiag();
	Pairing = source->getPairing();
	Betag = source->getBetag();
	Betas = source->getBetas();
	Pop = source->getPop();
	Bf = source->getFisBarrier();
	
	if (maxi!=source->getMaxi()||maxt!=source->getMaxt()||maxj!=source->getMaxj()) {
		deallocspectrum();
		maxi = source->getMaxi();
		maxt = source->getMaxt();
		maxj = source->getMaxj();
		allocspectrum(maxt,maxi,maxj);
	}
	
	for (int t=0;t<maxt;t++)
		for (int i=0;i<maxi;i++)
			for(int j=0;j<maxj;j++)
				Spectrum[t][i][j] = source->getSpectrum(t,i,j);

	return(*this);
}


Element::~Element() {

	dealloctranscoef();
	deallocspectrum();
	
	delete trans;
}

// ************************************************************************************************************* //
// ************************* Transmission coefficient calculation ********************************************** //
// ************************************************************************************************************* //

int Element::getlbmax(int nb){
	return(lbmax[nb]);
}

double Element::getTransCoefArray(int nb,int e,int lb,int jb) {
	return(TabCoef[nb][e][lb][jb]);
}

void Element::InitTransCoef(Readinput* path) {

	lbmax = new int[path->getNumber()];
	for (int nb=0;nb<path->getNumber();nb++) 
		lbmax[nb] = 0;
		
	int l = 0;
	int liml = 100;
	
	switch (path->getNumber()) {
		case 3:
			// Determination of lbmax for alpha
			l = 0;
			do {
				double j = (double(l)+0.5);
				trans->TransCal(path,A-4,Z-2,4,2,(maxi-1)*path->getEstep(),j,double(l));
				l++;
			} while (l<liml && trans->getTransCoef()>path->getCutTcoef());
			lbmax[2] = l;	
		case 2:
			// Determination of lbmax for proton
			l = 0;
			do {
				double j = (double(l)+0.5);
				trans->TransCal(path,A-1,Z-1,1,1,(maxi-1)*path->getEstep(),j,double(l));
				l++;
			} while (l<liml && trans->getTransCoef()>path->getCutTcoef());
			lbmax[1] = l;
		case 1:
			// Determination of lbmax for proton
			l = 0;
			do {
				double j = (double(l)+0.5);
				trans->TransCal(path,A-1,Z,1,0,(maxi-1)*path->getEstep(),j,double(l));
				l++;
			} while (l<liml && trans->getTransCoef()>path->getCutTcoef());
			lbmax[0] = l;
			break;
		default:
			cerr << "Warning: Please check the option for the type of evaporated particles." << endl;
			cin.get();
			exit(EXIT_FAILURE);
			break;
	}
	// lbmax could reach 15-20 in most cases. 

	jbmax = 3;
	nbmax = path->getNumber();
	emax = int((maxi-1)*path->getEstep()/path->getKinStep()+0.5);
//	emax = maxi;

	TabCoef = new double***[nbmax];
	for (int nb=0;nb<nbmax;nb++) {
		TabCoef[nb] = new double**[emax];
		for (int e=0;e<emax;e++) {
			TabCoef[nb][e] = new double*[lbmax[nb]];
			for (int lb=0;lb<lbmax[nb];lb++) {
				TabCoef[nb][e][lb] = new double[jbmax];
				for (int jb=0;jb<jbmax;jb++)
					TabCoef[nb][e][lb][jb] = 0.;	
			}
		}
	}

	switch (path->getNumber()) {
		case 3:
			for (int e=0;e<emax;e++) {
				double kin = (e+0.5)*path->getKinStep();
				int lb = 0.;
				do {
					double jb = double(lb); 
					trans->TransCal(path,A-4,Z-2,4,2,kin,jb,double(lb));
					TabCoef[2][e][lb][1] = trans->getTransCoef(); // for alpha one has (jb-lb+1)*2-1 = 1
					lb++;
				} while (lb<lbmax[2] && trans->getTransCoef()>path->getCutTcoef());
			}	
		case 2:
			for (int e=0;e<emax;e++) {
				double kin = (e+0.5)*path->getKinStep();
				int lb = 0;
				do {
					for (double jb=abs(double(lb)-0.5);jb<=(double(lb)+0.5);jb+=1.) {
						trans->TransCal(path,A-1,Z-1,1,1,kin,jb,double(lb));
						TabCoef[1][e][lb][int((jb-lb+1)*2-1)] = trans->getTransCoef(); // for nucleons one has (jb-lb+1)*2-1 = 0, 2
					}
					lb++;
				} while (lb<lbmax[1] && trans->getTransCoef()>path->getCutTcoef());
			}
		case 1:
			for (int e=0;e<emax;e++) {
				double kin = (e+0.5)*path->getKinStep();
				int lb = 0;
				do {	
					for (double jb=abs(double(lb)-0.5);jb<=(double(lb)+0.5);jb+=1.) {
						trans->TransCal(path,A-1,Z,1,0,kin,jb,double(lb)); 
						TabCoef[0][e][lb][int((jb-lb+1)*2-1)] = trans->getTransCoef(); // for nucleons one has (jb-lb+1)*2-1 = 0, 2
					}
					lb++;
				} while (lb<lbmax[0] && trans->getTransCoef()>path->getCutTcoef());
			}
			break;
		default:
			cerr << "Warning: Please check the option for the type of evaporated particles." << endl;
			cin.get();
			exit(EXIT_FAILURE);
			break;
	}
}

void Element::dealloctranscoef() {
/**/
	for (int nb=0;nb<nbmax;nb++) {
		for (int e=0;e<emax;e++) {
			for (int lb=0;lb<lbmax[nb];lb++) 
				delete[] TabCoef[nb][e][lb];
			delete[] TabCoef[nb][e];
		}
		delete[] TabCoef[nb];
	}
	delete[] TabCoef;
	TabCoef = NULL;

}
