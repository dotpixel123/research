#include <Fission.h>

double Fission::getFisWidth() {
	return(fiswidth);
}

void Fission::FisWidth(Readinput* path,Element* source,double Estarnew,double kin,double JC,double Bf) {
  
	fiswidth = 0.;  
  
	int errc = lev->Levdensj(path,source,Estarnew,JC,'g','f');
	if (errc==0) {
		double leveldeno = lev->getLevDens();
		double enr = Estarnew-Bf-kin;		
		int errf = lev->Levdensj(path,source,enr,JC,'s','f');
		if (errf==0 && leveldeno!=0.) {
			double levelnume = lev->getLevDens();
			fiswidth = levelnume/(leveldeno*2.*getPi());
		} else fiswidth = 0.;
	} else fiswidth = 0.;
	
	if (path->getFisTrans()==1) 
		fiswidth *= (1./(1.+exp(2.*getPi()*(-kin)/path->getOmegasd())));

	double Estar = Estarnew-source->getPairing();

	if (path->getKramers()==1) {		
		double a = 0.;
		switch (path->getLevDensPar()) {
			case 0:
				a = lev->parameter1(source->getA());
				break;
			case 1:
				a = lev->parameter2(path->getCstShell(),path->getShell(),source->getEsh(),source->getBetag(),
									Estar,source->getA(),source->getZ(),'g');
				break;
			case 2:
				a = lev->parameter3(path->getCstShell(),path->getShell(),source->getEsh(),source->getBetag(),
									Estar,source->getA(),source->getZ(),'g');
				break;
			case 3:
				a = lev->parameter4(path->getCstShell(),path->getShell(),source->getEsh(),source->getBetag(),
									Estar,source->getA(),source->getZ(),'g');
				break;
			default : 
				cerr << "Warning: There is a problem with the option for the level-density parameter." << endl;
				cin.get();
				exit(EXIT_FAILURE);
		}
		if (errc==0) {
			double temp = sqrt(Estar/a);
			double betaf = path->getFriction()*getHbar();
			double omegak = sqrt(path->getOmegasd()*path->getOmegasd()+betaf*betaf*0.25)-betaf*0.5; 
			if (omegak>0.) fiswidth *= (path->getOmegags()/path->getOmegasd())*(omegak/temp);
			if (Estar<=0.) fiswidth = 0.;
		} else fiswidth = 0.;
	}		
}

Fission::Fission() {

	lev = new Leveldensity();
}

Fission::~Fission() {
	
	delete lev;	
	lev = NULL;
}
