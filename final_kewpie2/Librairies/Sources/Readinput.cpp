#include <Readinput.h>

int Readinput::getZpro() {
	return(Zpro);
}

int Readinput::getApro() {
	return(Apro);
}

int Readinput::getZtar() {
	return(Ztar);
}

int Readinput::getAtar() {
	return(Atar);
}

int Readinput::getShell() {
	return(Shell);
}

double Readinput::getShellFactor() {
	return(ShellFactor);
}

double Readinput::getKinStep() {
	return(KinStep);
}

int Readinput::getLevDensPar() {
	return(LevDensPar);
}

double Readinput::getCstShell() {
	return(CstShell);
}

double Readinput::getafan() {
	return(afan);
}

int Readinput::getMollerTable() {
	return(MollerTable);
}

int Readinput::getPair() {
	return(Pair);
}

int Readinput::getEva() {
	return(Eva);
}

int Readinput::getOptPot() {
	return(OptPot);
}

int Readinput::getKramers() {
	return(Kramers);
}

int Readinput::getCollEnhFac() {
	return(CollEnhFac);
}

double Readinput::getFriction() {
	return(Friction);
}

double Readinput::getOmegags() {
	return(Omegags);
}

double Readinput::getOmegasd() {
	return(Omegasd);
}

int Readinput::getFisBar() {
	return(FisBar);
}

int Readinput::getFisTrans() {
	return(FisTrans);
}

double Readinput::getEmin() {
	return(Emin);
}

double Readinput::getEmax() {
	return(Emax);
}

double Readinput::getEstep() {
	return(Estep);
}

int Readinput::getJmin() {
	return(Jmin);
}

int Readinput::getJmax() {
	return(Jmax);
}

int Readinput::getn() {
	return(n);
}

int Readinput::getp() {
	return(p);
}

int Readinput::getNumber() {
	return(Number);
}

int Readinput::getGamma() {
	return(Gamma);
}

int Readinput::getGammaModel() {
	return(GammaModel);
}

int Readinput::getTime() {
	return(Time);
}

double Readinput::getTimeFact() {
	return(TimeFact);
}

double Readinput::getDelay(){
  return(Delay);
}

double Readinput::getCutTcoef() {
	return(CutTcoef);
}

double Readinput::getCutWidth() {
	return(CutWidth);
}

double Readinput::getFacInerMom() {
	return(FacInerMom);
}

double Readinput::getEview() {
	return(Eview);
}

int Readinput::getFusion() {
	return(Fuse);
}

double Readinput::getDeltaBvar() {
	return(DeltaBvar);
}

double Readinput::getDeltaB0() {
	return(DeltaB0);
}

size_t Readinput::getNumberPoint() {
	return(NumPoint);
}

string Readinput::getFusInput() {
	return(FusInput);
}

string Readinput::getOutput() {
	return(Output); 
}

string Readinput::getMassTable() {
	return(MassTable); 
}

string Readinput::getFisBarTable() {
	return(FisBarTable); 
}

void Readinput::read() {
  
	macro>>Zpro;
	do {macro>>reader;} while(reader!="charge");
	macro>>Apro;
	do {macro>>reader;} while(reader!="mass");
	macro>>Ztar;
	do {macro>>reader;} while(reader!="charge");
	macro>>Atar;
	do {macro>>reader;} while(reader!="mass");
	macro>>Shell;
	do {macro>>reader;} while(reader!="used)");
	macro>>ShellFactor;
	do {macro>>reader;} while(reader!="1)");
	macro>>KinStep;  
	do {macro>>reader;} while(reader!="(MeV)");
	macro>>LevDensPar;
	do {macro>>reader;} while(reader!="README)");
	macro>>CstShell;
	do {macro>>reader;} while(reader!="(MeV)");
	macro>>afan;
	do {macro>>reader;} while(reader!="ratio");
	macro>>FacInerMom;
	do {macro>>reader;} while(reader!="inertia");
	macro>>MollerTable;
	do {macro>>reader;} while(reader!="file)");
	macro>>Pair;
	do {macro>>reader;} while(reader!="used)");
	macro>>Kramers;
	do {macro>>reader;} while(reader!="used)");
	macro>>CollEnhFac;
	do {macro>>reader;} while(reader!="used)");
	macro>>Friction;
	do {macro>>reader;} while(reader!="(s^-1)");
	macro>>Omegags;
	do {macro>>reader;} while(reader!="(MeV)");
	macro>>Omegasd;
	do {macro>>reader;} while(reader!="(MeV)");
	macro>>Eva;
	do {macro>>reader;} while(reader!="WE)");
	macro>>OptPot;
	do {macro>>reader;} while(reader!="Koning)");
	macro>>FisTrans;
	do {macro>>reader;} while(reader!="used)");
	macro>>FisBar;
	do {macro>>reader;} while(reader!="LSD)");
	macro>>Emin;
	do {macro>>reader;} while(reader!="(MeV)");
	macro>>Emax;
	do {macro>>reader;} while(reader!="(MeV)");
	macro>>Estep;
	do {macro>>reader;} while(reader!="(MeV)");
	macro>>Eview;
	do {macro>>reader;} while(reader!="(MeV)");
	macro>>Jmin;
	do {macro>>reader;} while(reader!="(hbar)");
	macro>>Jmax;
	do {macro>>reader;} while(reader!="(hbar)");
	macro>>n;
	do {macro>>reader;} while(reader!="neutrons");
	macro>>p;
	do {macro>>reader;} while(reader!="protons");
	macro>>Number;
	do {macro>>reader;} while(reader!="README)");
	macro>>Gamma;
	do {macro>>reader;} while(reader!="used)");
	macro>>GammaModel;
	do {macro>>reader;} while(reader!="SMLO)");
	macro>>Time;
	do {macro>>reader;} while(reader!="used)");
	macro>>TimeFact;
	do {macro>>reader;} while(reader!="calculation");
	macro>>Delay;
  	do {macro>>reader;} while(reader!="process");
	macro>>CutTcoef;
	do {macro>>reader;} while(reader!="coefficient");
	macro>>CutWidth;
	do {macro>>reader;} while(reader!="(MeV)");
	macro>>Fuse;
	do {macro>>reader;} while(reader!="file)");
	macro>>DeltaB0>>DeltaBvar;
	do {macro>>reader;} while(reader!="README)");
	macro>>NumPoint;
	do {macro>>reader;} while(reader!="README)");
	macro>>FusInput;
	do {macro>>reader;} while(reader!="used)");
	macro>>MassTable;
	do {macro>>reader;} while(reader!="used)");
	macro>>FisBarTable;
	do {macro>>reader;} while(reader!="README)");
	macro>>Output;
	do {macro>>reader;} while(reader!="file");
}

int Readinput::check() {

	int error = 1;
  
	if (Zpro<0 || Zpro>140 || Ztar<0 || Ztar>140) error = 0;
	if (Apro<0 || Apro>350 || Atar<0 || Atar>350) error = 0;
	if (Shell!=0 && Shell!=1) error = 0;
	if (LevDensPar!=0 && LevDensPar!=1 && LevDensPar!=2 && LevDensPar!=3) error = 0;
	if (MollerTable!=0 && MollerTable!=1 && MollerTable!=2) error = 0;
	if (FacInerMom<=0.) error = 0;
	if (Pair!=0 && Pair!=1) error = 0;
	if (Eva!=0 && Eva!=1) error = 0;
	if (OptPot!=0 && OptPot!=1) error = 0;
	if (GammaModel!=0 && GammaModel!=1) error = 0;
	if (Kramers!=0 && Kramers!=1) error = 0;
	if (CollEnhFac!=0 && CollEnhFac!=1) error = 0;
	if (FisTrans!=0 && FisTrans!=1) error = 0;
	if (FisBar!=0 && FisBar!=1 && FisBar!=2) error = 0;
	if (Number!=1 && Number!=2 && Number!=3) error = 0;
	if (Jmin<0 || Jmin>200 || Jmax<0 || Jmax>200) error = 0;
	if (Jmin>Jmax) error = 0;
	if (Gamma!=0 && Gamma!=1) error = 0;
	if (Time!=0 && Time!=1) error = 0;
	if (Emin<0.) error = 0;
	if (Emin>Emax) error = 0;
	if (ShellFactor<0.) error = 0;
	if (KinStep<=0.) error = 0;
	if (Estep<=0.) error = 0;
	if (KinStep>Estep) error = 0;
	if (Friction<0.) error = 0;
	if (CstShell<0.) error = 0;
	if (Delay<0.) error = 0;
	if (TimeFact<=0.) error = 0;
	if (Fuse!=0 && Fuse!=1 && Fuse!=2) error = 0;
	if ((NumPoint%2)==1) error = 0;
	if (n<0) error = 0;
	if (p<0) error = 0;

	return(error);
}

Readinput::Readinput() {

	macro.open("../Input/input.inp");
	if (macro.fail()) {
		cerr << "cannot open the input file" << endl;
		macro.close();
		cin.get();
		exit(EXIT_FAILURE);
	}

	Zpro = 0;
	Apro = 0;
	Ztar = 0;
	Atar = 0;
	Output = " ";
}

Readinput::~Readinput() {

	macro.close();
}
