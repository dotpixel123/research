#include <Decay.h>

bool Decay::getBool(int i) {
	return(tabbool[i]);
}

double Decay::getFisSum() {
	return(FisTime);
}

void Decay::WidthCall(Readinput* path,Element* source,Element* residue,Particle* evap,double* gammatot,double DeltaT,bool Delay,int nbT) {
	
	if (path->getEva()==0)
		source->InitTransCoef(path);
//	cout << source->getZ() << " " << source->getA() << " " << source->getlbmax(0) << endl;
//	cin.get();

	int nbTres = nbT;
	if (path->getTime()==1) {
		nbTres = nbT+1;		
		FisTime = 0.;	
		InitWidth();
	}
	
	for (int q=0;q<number;q++) 
		tabbool[q] = false;
	
	double Bf = source->getFisBarrier();
	
	for (int q=0;q<path->getNumber();q++) {
		Se[q] = evawidth->SeCal(residue[q].getMass(),evap[q].getMass(),source->getMass());
		Bc[q] = evawidth->BcCal(residue[q].getA(),residue[q].getZ(),evap[q].getA(),evap[q].getZ());
	}

//	cout << source->getZ() << " " << source->getA() << " " << Bf << endl;
//	cin.get();
//	cout << source->getZ() << " " << source->getA() << " " << Se[0]<< endl;
//	cin.get();

//	if (source->getZ()==90 && source->getA()==224) Se[0] = 7.46349;
//	if (source->getZ()==90 && source->getA()==223) Se[0] = 5.88852;
//	if (source->getZ()==90 && source->getA()==222) Se[0] = 7.80635;
//	if (source->getZ()==90 && source->getA()==221) Se[0] = 5.80227;
//	if (source->getZ()==90 && source->getA()==220) Se[0] = 7.87508;	

	int jmin = 0;
	int jmax = 0;	
	if (path->getEva()==0) {
		jmin = path->getJmin();
		jmax = path->getJmax();
	}
	
	for (int i=0;i<source->getMaxi();i++) { // loop on excitation energy

		double width = 0.;			
		double evap_width = 0.;
		double fiss_width = 0.;
		double gamm_width = 0.;
		double excit = i*path->getEstep();
		
		for (int j=jmin;j<=jmax;j++) { // loop on angular momentum

			// Initialization of local variables
			double widthsum = 0.;
			double fissum = 0.;
			double JC = j*1.;
			if (source->getA()%2==1) 
				JC = j*1.+0.5;
				
			InitGammaEva();	
			
			if (source->getSpectrum(nbT,i,j)>0.) { 			

				int kmax = int(excit/path->getKinStep()+0.5); // number of bins for kinetic energy
				
				for (int k=0;k<kmax;k++) {
	      			
		      		double ener = (k+0.5)*path->getKinStep();		
	      			// Evaporation width calculation
	      			for (int q=0;q<path->getNumber();q++) {
						if (excit-ener-residue[q].getPairing()-Se[q]-Bc[q]-path->getKinStep()>=0.) {	
							int exi2 = int((excit-ener-Se[q]-Bc[q])/path->getEstep()+0.5); // kinetic energy above the Coulomb barrier
						//	cout << evap[q].getZ() << " " << evap[q].getA() << " " << Bc[q] << " " << Se[q] << endl;
						//	cout << maxi << " " << exi2 << endl;
						//	cin.get();
							if (path->getTime()==0) {
								label1: 
									if (path->getEva()==0) {
										width = 0.;
										for (int lb=0;lb<source->getlbmax(q);lb++) {
											double Tcoef = 0.;
											for (double jb=abs(double(lb)-evap[q].getSpin());jb<=double(lb)+evap[q].getSpin();jb+=1.) {
												Tcoef = source->getTransCoefArray(q,k,lb,int((jb-lb+1)*2-1));
												for (double IB=abs(JC-jb);IB<=JC+jb;IB+=1.) {
													evawidth->Hauser_Feshbach(path,source,&residue[q],excit,ener+Bc[q],JC,IB,Se[q],Tcoef);
													width += evawidth->getEvaWidth()*path->getKinStep();
													if (int(IB)<maxj) 
														gammaeva[q][exi2][int(IB)] += evawidth->getEvaWidth()*path->getKinStep();
												}
											}
										}
									} else {
										evawidth->Weisskopf_Ewing(path,source,&residue[q],&evap[q],excit,ener+Bc[q],JC,Se[q]);
										width = evawidth->getEvaWidth()*path->getKinStep();
										gammaeva[q][exi2][0] += width;
									}
									if (path->getTime()==1) 
										Width[q][i][k] = width;
							} else {
								if (Width[q][i][k]!=0.) 
									width = Width[q][i][k];
								else goto label1;
							} // end if(Time==0)
							if (width>0. && !tabbool[q]) tabbool[q] = true; // control for element loops in Cascade 
							widthsum += width;
							evap_width += width;
						} // end condition ener kinetic energy
					} // end evaporation loop q type of particles
					
	      			// Fission width calculation
					if (!Delay) {
						if ((excit-ener-Bf-source->getPairing()-path->getKinStep())>=0.) {
							if (path->getTime()==0) {
								label2:
									fiswidth->FisWidth(path,source,excit,ener,JC,Bf);
									width = fiswidth->getFisWidth()*path->getKinStep();
									if (path->getTime()==1) Width[path->getNumber()][i][k] = width;
							} else {
								if (Width[path->getNumber()][i][k]!=0.) 
									width = Width[path->getNumber()][i][k];
								else goto label2;
							}
						} else width = 0.; // end if (Estar-e-Ep>Bf)				
						widthsum += width;
						fiss_width += width;
						if (path->getTime()==1) 
							fissum += width;   
					} // end if (!Delay)
	      
	      			// Gamma width calculation
					if (path->getGamma()==1 && widthsum!=0.) {  						 
						if (excit-ener-source->getPairing()-path->getKinStep()>=0.) {
							int exi2 = int((excit-ener)/path->getEstep());  	  
							if (path->getTime()==0) {
								label3:
									if (path->getEva()==0) {
										width = 0.;
										gamwidth->RadStrFuncs(path,source,excit,ener);
										for (int le=1;le<=2;le++) {
											double Tcoef = gamwidth->getfe(le)+gamwidth->getfm(le);
											for (double IB=abs(JC-double(le));IB<=JC+double(le);IB++) {
												gamwidth->GamWidth_J(path,source,excit,ener,JC,IB,Tcoef);
												width += gamwidth->getGamWidth()*path->getKinStep();
												if (int(IB)<maxj) 
													gammaeva[number-1][exi2][int(IB)] += gamwidth->getGamWidth()*path->getKinStep();
											}
										}
									} else {
										gamwidth->GamWidth(path,source,excit,ener,JC);
										width = gamwidth->getGamWidth()*path->getKinStep();
										gammaeva[number-1][exi2][0] += width;
									}
									if (path->getTime()==1) 
										Width[number][i][k] = width;
							} else {
								if (Width[number][i][k]!=0.) 
									width = Width[number][i][k];
								else goto label3;
							}  // end if(Time==0)
							widthsum += width;
							gamm_width += width;
							if (width>0. && !tabbool[number-1]) 
								tabbool[number-1] = true; // control for element loops in Cascade 
						}
					}
				} // end loop on kinetic energy
	    
				gammatot[i] = widthsum;
				
				if (path->getTime()==1) FisTime += fissum*source->getSpectrum(nbT,i,j);
				
				if (widthsum!=0.) {
					for (int q=0;q<number;q++) {
						if (tabbool[q]) { // 
							for (int e=0;e<residue[q].getMaxi();e++) {
								for (int m=0;m<residue[q].getMaxj();m++) {
									// cout << residue[q].getMaxi() << "\t" << residue[q].getMaxj() << endl;
									double factor = 0.;
									if (path->getTime()==1) factor = exp(-gammatot[i]*DeltaT);
										width = gammaeva[q][e][m]*source->getSpectrum(nbT,i,j)*(1.-factor)/gammatot[i]
												+residue[q].getSpectrum(nbTres,e,m);
									if (width>path->getCutWidth()) 
										residue[q].setSpectrum(nbTres,e,m,width);
								// Extract spectrum of the daughter nucleus
								// if (residue[q].getA()==255) 
								//	cout << e*path->getEstep() << "\t" << m << "\t" << residue[q].getSpectrum(nbTres,e,m) << endl;
								// else if (residue[q].getA()==254) 
								//	cout << e*path->getEstep() << "\t" << m << "\t" << residue[q].getSpectrum(nbTres,e,m) << endl;
								// else 
								//	continue;
								}
							}
						}
					}
//	cout << source->getZ() << "\t" << source->getA() << "\t" << excit << "\t" << evap_width << "\t" << fiss_width << endl;			
//	cout << scientific << fiss_width << "\t" << sum_n << "\t" << sum_p << "\t" << sum_a << "\t" << gamm_width << endl;
//	exit(0);		
				} // end if (widthsum!=0)			
			} // end on the condition to start the calculation
		} // end loop on angular momentum
	} // end loop on excitation energy	
//	cin.get();
}

void Decay::InitGammaEva() {

	for (int q=0;q<number;q++)		
		for (int i=0;i<maxi;i++)
			for (int j=0;j<maxj;j++)
				gammaeva[q][i][j] = 0.;
}

void Decay::InitWidth() {

	for (int q=0;q<number+1;q++) // type of emitted particle
		for (int i=0;i<maxi;i++) // excitation energy
			for (int k=0;k<maxi;k++) // kinetic energy
				Width[q][i][k] = 0.;
}
 
Decay::Decay(int Number,int Maxi,int Maxj) {

	number = Number;
	maxi = Maxi;
	maxj = Maxj;
	
	FisTime = 0.;

	evawidth = new Evaporation();
	fiswidth = new Fission();
	gamwidth = new Gamma();
	Se = new double[number];
	Bc = new double[number];
	tabbool = new bool[number];

	gammaeva = new double**[number];
	for (int q=0;q<number;q++) {
		tabbool[q] = false;
		Se[q] = 0.;
		Bc[q] = 0.;
		gammaeva[q] = new double*[maxi];
		for (int i=0;i<maxi;i++) {
			gammaeva[q][i] = new double[maxj];
			for (int j=0;j<maxj;j++)	
				gammaeva[q][i][j] = 0.;
		}
	}

	Width = new double**[number+1];
	for (int q=0;q<number+1;q++) {
		Width[q] = new double*[maxi];
		for (int i=0;i<maxi;i++) {
		    Width[q][i] = new double[maxi];
			for (int k=0;k<maxi;k++) 
				Width[q][i][k] = 0.;	
		}
	}
}

Decay::~Decay() {

	delete[] Se;
	delete[] tabbool;
	delete evawidth;
	delete fiswidth;
	delete gamwidth;
	
	for (int q=0;q<number;q++) {
		for (int i=0;i<maxi;i++)
			delete[] gammaeva[q][i];
		delete[] gammaeva[q];
	}
	delete[] gammaeva;
	
	for (int q=0;q<number+1;q++) {
		for (int i=0;i<maxi;i++)
			delete[] Width[q][i];
		delete[] Width[q];
	}
	delete[] Width;
	
	Se = NULL;
	tabbool = NULL;
	evawidth = NULL;
	fiswidth = NULL;
	gamwidth = NULL;
	gammaeva = NULL;
	Width = NULL;
}
