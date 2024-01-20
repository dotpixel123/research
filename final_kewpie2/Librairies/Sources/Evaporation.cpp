#include <Evaporation.h>

double Evaporation::getEvaWidth() {       
	return(evawidth);
}

double Evaporation::SeCal(double Mres,double Memi,double Msrc) {
	return(Msrc-Mres-Memi);
}

double Evaporation::BcCal(int Ares,int Zres,int Aemi,int Zemi) {
	
	if (Zemi==0 && Aemi==1) return(0.);
	else if (Zemi==1 && Aemi==1) {
		double Ri = 1.44;
		double Re = 1.81; 
		double Vc = Zemi*Zres*1.44/(Re*cbrt(Ares*1.)+Ri);
		return(Vc);
	}
	else if (Zemi==2 && Aemi==4) {
		double Ri = 2.53;
		double Re = 2.452-0.408*log10(1.*Zemi*Zres);
		double Vc = Zemi*Zres*1.44/(Re*cbrt(Ares*1.)+Ri);
		return(Vc);
	}
	else {
		cerr << "Warning: Problem with the type of emitted particle." << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
}

// ************************************************************************************************************* //
// **************************************** Hauser-Feshbach formalism ****************************************** //
// ************************************************************************************************************* //

void Evaporation::Hauser_Feshbach(Readinput* path,Element* source,Element* residue,double Estarnew,
									double kin,double JC,double IB,double Sn,double Trans) {

	evawidth = 0.;

	int errc = lev->Levdensj(path,source,Estarnew,JC,'g','e');
	if (errc==0) {  
		double leveldeno = lev->getLevDens();
		double enr = Estarnew-kin-Sn;
		int errf = lev->Levdensj(path,residue,enr,IB,'g','e'); 
		if (errf==0 && leveldeno!=0.) {
			double levelnume = lev->getLevDens();
			evawidth = levelnume*Trans/(leveldeno*2.*getPi());
		} else evawidth = 0.; // end condition on errf
	} else evawidth = 0.; // end condition on errc
}

// ************************************************************************************************************* //
// **************************************** Weisskopf-Ewing model ********************************************** //
// ************************************************************************************************************* //

void Evaporation::Weisskopf_Ewing(Readinput* path,Element* source,Element* residue,Particle* emit,double Estarnew,
						double kin, double JC,double Sn) {
	double sigma = 0.;
	double Ri = 0.;
	double Re = 0.;
  
	if (emit->getZ()==0) {
		sigma = (0.76+1.93/cbrt(residue->getA()*1.))*kin+(1.66/cbrt(residue->getA()*residue->getA()*1.)-0.05);
		sigma *= getPi()*1.7*cbrt(residue->getA()*1.)*1.7*cbrt(residue->getA()*1.);
		if (sigma<0.) sigma = 0.;
	} else {
		if (emit->getZ()==1) {
			Ri = 1.44;
			Re = 1.81; 
		}
		if (emit->getZ()==2) {
			Ri = 2.53;
			Re = 2.452-0.408*log10(1.*emit->getZ()*residue->getZ());
		}
		sigma = kin-emit->getZ()*residue->getZ()*1.44/(Re*cbrt(residue->getA()*1.)+Ri);
		sigma *= getPi()*(1.42*cbrt(residue->getA()*1.)+Ri)*(1.42*cbrt(residue->getA()*1.)+Ri);
		if (sigma<0.) sigma = 0.; // classical description of inverse cross-section
	} 

	evawidth = 0.;

	int errf = lev->Levdensj(path,source,Estarnew,JC,'g','e');	
	if (errf==0) {
		double rhoc = lev->getLevDens();
		double enr = Estarnew-kin-Sn;		
		int errc = lev->Levdensj(path,residue,enr,JC,'g','e');
		if (errc==0 && rhoc!=0.) {
			double rhof = lev->getLevDens();
			evawidth = rhof/rhoc;
		} else evawidth = 0.;
	} else evawidth = 0.;
	 
	double g = 2.*emit->getSpin()+1.;
	double mu = (residue->getA()*emit->getA()*1.)/(residue->getA()+emit->getA())*getAmu();
  	
	evawidth = g*mu*sigma/((getPi()*getHbarc())*(getPi()*getHbarc()))*evawidth;
}

Evaporation::Evaporation() {
	
	lev = new Leveldensity();
}

Evaporation::~Evaporation() {     
	
	delete lev;	
	lev = NULL;
}
