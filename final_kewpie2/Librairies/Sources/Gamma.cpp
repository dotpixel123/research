#include <Gamma.h>

double Gamma::getGamWidth() {
	return(gamwidth);
}

double Gamma::getfe(int i) {
	return(coefe[i]);
}

double Gamma::getfm(int i) {
	return(coefm[i]);
}

double Gamma::RSFE1SMLO(int Asrc,int Zsrc,double a,double beta2,double Estar,double e) {

	double fe1 = 0.;
		
	if (e<Estar) {
			
		double Tempi = sqrt(Estar/a);
		double Tempf = sqrt((Estar-e)/a);
		double lambdaf = 1./(1.-exp(-e/Tempf));
		double I = (Asrc-2.*Zsrc)/Asrc;
			
		double a1=28.69, a2=21.731, a3=0.33078, a4=1.2669, b1=0., b2=0.; // Parametrization of the resonance energy
		double er = a1*(1.+b1*I*I)/cbrt(1.*Asrc)+a2*(1.+b2*I*I)/cbrt(sqrt(1.*Asrc));
		double gamma = a3*er;
		double sigma = 120.*a4*(Asrc-1.*Zsrc)*Zsrc/(Asrc*getPi()*gamma);
		double alpha2 = beta2*sqrt(5./(4.*getPi()));
	
		// Calculation of RSFE1EGLO
		double a0 = (1.+alpha2)/cbrt(1.+0.6*alpha2*alpha2+2.*alpha2*alpha2*alpha2/35.);
		double b0 = (1.-0.5*alpha2)/cbrt(1.+0.6*alpha2*alpha2+2.*alpha2*alpha2*alpha2/35.);
		double er2 = er*(1.-1.51*0.01*(a0*a0-b0*b0))/b0;
		double er1 = er2/(0.911*a0/b0 + 0.089);
		double gamma1 = a3*er1;
		double gamma2 = a3*er2;
		double sigma1 = sigma/3.;
		double sigma2 = sigma*2./3.;
		double eta1 = gamma1/er1;
		double eta2 = gamma2/er2;
		double gammaT1 = eta1*Estar;
		double gammaT2 = eta2*Estar;

		if (e==0.) {
			fe1 = 8.674e-8*sigma1*gamma1*Tempi*gammaT1/(er1*er1*er1*er1);
			fe1 += 8.674e-8*sigma2*gamma2*Tempi*gammaT2/(er2*er2*er2*er2);
		} else {
			fe1 = 8.674e-8*sigma1*gamma1*lambdaf*e*gammaT1/((e*e-er1*er1)*(e*e-er1*er1)+e*e*gammaT1*gammaT1);
			fe1 += 8.674e-8*sigma2*gamma2*lambdaf*e*gammaT2/((e*e-er2*er2)*(e*e-er2*er2)+e*e*gammaT2*gammaT2);
		}	
	}
		
	return(fe1);
}

double Gamma::RSFE1EGLO(int Asrc,int Zsrc,double a,double beta2,double Estar,double e) {

	double fe1 = 0.;	
	
	if (e<Estar) {
		
		double Tempi = sqrt(Estar/a);
		double Tempf = sqrt((Estar-e)/a);
		double I = (Asrc-2.*Zsrc)/Asrc;
		
		double a1=27.469, a2=22.063, a3=0.02691, a4=1.2224, b1=0., b2=0.;  // Parametrization of the resonance energy
		double er = a1*(1.+b1*I*I)/cbrt(1.*Asrc)+a2*(1.+b2*I*I)/cbrt(sqrt(1.*Asrc));
		double gamma = a3*pow(er,1.91);
		double sigma = 120.*a4*(Asrc-1.*Zsrc)*Zsrc/(Asrc*getPi()*gamma);
	
		double kaux = 1.;
		if (Asrc>=148) kaux=1.+0.09*(Asrc-148.)*(Asrc-148.)*exp(-0.18*(Asrc-148.)); // according to BSFG 
		
		double alpha2 = beta2*sqrt(5./(4.*getPi()));
			
		// Calculation of RSFE1SMLO
		double a0 = (1.+alpha2)/cbrt(1.+0.6*alpha2*alpha2+2.*alpha2*alpha2*alpha2/35.);
		double b0 = (1.-0.5*alpha2)/cbrt(1.+0.6*alpha2*alpha2+2.*alpha2*alpha2*alpha2/35.);
		double er2 = er*(1.-1.51*0.01*(a0*a0-b0*b0))/b0;
		double er1 = er2/(0.911*a0/b0 + 0.089);
		double gamma1 = a3*pow(er1,1.91);
		double gamma2 = a3*pow(er2,1.91);
		double sigma1 = sigma/3.;
		double sigma2 = sigma*2./3.;
		double epsilon0 = 4.5;
		double ksi1 = kaux+(1.-kaux)*(e-epsilon0)/(er1-epsilon0);
		double ksi2 = kaux+(1.-kaux)*(e-epsilon0)/(er2-epsilon0);
			
		double gammaTf1 = gamma1*ksi1/er1/er1*(e*e+4.*getPi()*getPi()*Tempf*Tempf);
		double gammaTf2 = gamma2*ksi2/er2/er2*(e*e+4.*getPi()*getPi()*Tempf*Tempf);
		double gammaTi1 = gamma1*ksi1/er1/er1*(4.*getPi()*getPi()*Tempi*Tempi);
		double gammaTi2 = gamma2*ksi2/er2/er2*(4.*getPi()*getPi()*Tempi*Tempi);
		
		fe1 = 8.674e-8*sigma1*gamma1*(e*gammaTf1/((e*e-er1*er1)*(e*e-er1*er1)+e*e*gammaTf1*gammaTf1)+0.7*gammaTi1/er1/er1/er1);
		fe1 += 8.674e-8*sigma2*gamma2*(e*gammaTf2/((e*e-er2*er2)*(e*e-er2*er2)+e*e*gammaTf2*gammaTf2)+0.7*gammaTi2/er2/er2/er2);	
	}	

	return(fe1);
}

double Gamma::RSFE2(int Asrc,int Zsrc,double e) {
	
	double ere2 = 63./cbrt(1.*Asrc);
	double gammae2 = 6.11-0.012*Asrc;
	double sigmae2 = 1.5*0.0001*Zsrc*Zsrc*ere2*ere2/(cbrt(Asrc*1.)*gammae2);
	
	double fe2 = 5.2e-8*sigmae2*gammae2/e*gammae2/((e*e-ere2*ere2)*(e*e-ere2*ere2)+e*e*gammae2*gammae2);
	return(fe2);
}

double Gamma::RSFM1(int Asrc,int Zsrc,double e) {

	double erm1 = 41./cbrt(1.*Asrc);
	double gammam1 = 4.0;
	double sigmam1 = 1.0;

	double fm1 = 8.674e-8*sigmam1*gammam1*e*gammam1/((e*e-erm1*erm1)*(e*e-erm1*erm1)+e*e*gammam1*gammam1);
	return(fm1);
}

void Gamma::RadStrFuncs(Readinput* path,Element* source,double Estarnew,double e) {

	double Estar = Estarnew-source->getPairing();

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
			default:
				cerr << "Warning: There is a problem with the option of the level-density parameter in the input file." << endl;
				cin.get();
				exit(EXIT_FAILURE);
	}

	if (e<Estar) { // e cannot be greater than Estar

		double fe1 = 0.;
		switch (path->getGammaModel()) {
				case 0:
					fe1 = RSFE1EGLO(source->getA(),source->getZ(),a,source->getBetag(),Estar,e);
			    	break;
			    case 1:
					fe1 = RSFE1SMLO(source->getA(),source->getZ(),a,source->getBetag(),Estar,e);
					break;
				default:
					cerr << "Warning: There is a problem with the calculation of RSF." << endl;
					cin.get();
					exit(EXIT_FAILURE);
		}
	
		double fe2 = RSFE2(source->getA(),source->getZ(),e);
		double fm1 = RSFM1(source->getA(),source->getZ(),e);
		double fm2 = 0.307/cbrt(source->getA()*source->getA()*1.)*fe2;
		
		coefe[1] = fe1*2.*getPi()*e*e*e;
		coefe[2] = fe2*2.*getPi()*e*e*e*e*e;
		coefm[1] = fm1*2.*getPi()*e*e*e;
		coefm[2] = fm2*2.*getPi()*e*e*e*e*e;
		
//		coefe[1] = (1.8432E-10)*pow(source->getA(),2./3.);
//		coefe[2] = (8.9365E-12)*pow(source->getA(),4./3.);
	}
}

void Gamma::GamWidth_J(Readinput* path,Element* source,double Estarnew,double e,double JC,double IB,double Trans) {

	gamwidth = 0.;
	
	int errc = lev->Levdensj(path,source,Estarnew,JC,'g','e');
	if (errc==0) {
		double leveldeno = lev->getLevDens();
		double enr = Estarnew-e;
		int errf = lev->Levdensj(path,source,enr,IB,'g','e');
		if (errf==0 && leveldeno!=0.) {
			double levelnume = lev->getLevDens();
			gamwidth = levelnume*Trans/(leveldeno*2.*getPi());
		} else gamwidth = 0.;	
	} else gamwidth = 0.;
}

void Gamma::GamWidth(Readinput* path,Element* source,double Estarnew,double e,double JC) {

	double Estar = Estarnew-source->getPairing();

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
			default:
				cerr << "Warning: There is a problem with the option of the level-density parameter in the input file." << endl;
				cin.get();
				exit(EXIT_FAILURE);
	}
	
	gamwidth = 0.;

	if (e<Estar) { // e cannot be greater than Estar

		double fe1 = 0.;
		switch (path->getGammaModel()) {
				case 0:
					fe1 = RSFE1EGLO(source->getA(),source->getZ(),a,source->getBetag(),Estar,e);
			    	break;
			    case 1:
					fe1 = RSFE1SMLO(source->getA(),source->getZ(),a,source->getBetag(),Estar,e);
					break;
				default:
					cerr << "Warning: There is a problem with the calculation of RSF." << endl;
					cin.get();
					exit(EXIT_FAILURE);
		}
	
		double fe2 = RSFE2(source->getA(),source->getZ(),e);
		double fm1 = RSFM1(source->getA(),source->getZ(),e);
		double fm2 = 0.307/cbrt(source->getA()*source->getA()*1.)*fe2;
		
		coefe[1] = fe1*2.*getPi()*e*e*e;
		coefe[2] = fe2*2.*getPi()*e*e*e*e*e;
		coefm[1] = fm1*2.*getPi()*e*e*e;
		coefm[2] = fm2*2.*getPi()*e*e*e*e*e;

//		Weisskopf estimates:
//		coefe[1] = (1.8432E-10)*pow(source->getA(),2./3.); 
//		coefe[2] = (8.9365E-12)*pow(source->getA(),4./3.);
		
		int errc = lev->Levdensj(path,source,Estarnew,JC,'g','e');
		if (errc==0) {
			double leveldeno = lev->getLevDens();
			for (int l=1;l<=2;l++) {
				for (double IB=abs(JC-double(l));IB<=JC+double(l);IB++) {
					double enr = Estarnew-e;
					int errf = lev->Levdensj(path,source,enr,IB,'g','e');
					if (errf==0 && leveldeno!=0.) {
						double levelnume = lev->getLevDens();
						gamwidth += levelnume*(coefe[l]+coefm[l])/(leveldeno*2.*getPi());
					}	
				} 
			}
		} else gamwidth = 0.;
	} else gamwidth = 0.;
}

Gamma::Gamma() {

	coefe = new double[3];
	coefm = new double[3];
	for (int j=0;j<3;j++) {
		coefe[j] = 0.;
		coefm[j] = 0.;
	}
	
	lev = new Leveldensity();
}

Gamma::~Gamma() {

	delete[] coefe;
	delete[] coefm;
	delete lev;
	
	coefe = NULL;
	coefm = NULL;
}
