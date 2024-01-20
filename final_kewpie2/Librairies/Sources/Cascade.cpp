#include <Cascade.h>

void Cascade::CascadeDecay(Readinput* path,Table* MassTable,Element* source) {

	int rang = 0; // rank of each step increased by one at each time
	
	// first element at the beginning of the decay chain
	tabelem[0].allocspectrum(MaxT,MaxI,MaxJ);
	tabelem[0] = source;
	for (int q=0;q<Number;q++) 
		residuetab[q].allocspectrum(MaxT,tabelem[0].getMaxi(),MaxJ);
	
	FissionTime = 0.;	
	comp = 1; // number of residual nuclei at each step
	
	int num;
	do {
	  
		TimeDist = 0.; 
		num = 0;
		
		for (int r=0;r<comp;r++) {
	    
			if (tabelem[r].getA()==source->getA()-rang) {	
	      
				for (int q=0;q<Number;q++) {
					Neu[q] = 0;
					Pro[q] = 0;
					residuetab[q].InitSpectrum();
					int Aaux = tabelem[r].getA()-evap[q].getA();
					int Zaux = tabelem[r].getZ()-evap[q].getZ();
					double massaux = MassTable->getMassArray((tabelem[r].getA()-evap[q].getA()),(tabelem[r].getZ()-evap[q].getZ()));
					double shellaux = MassTable->getShellArray((tabelem[r].getA()-evap[q].getA()),(tabelem[r].getZ()-evap[q].getZ()));
					double beta2aux = MassTable->getBetaArray((tabelem[r].getA()-evap[q].getA()),(tabelem[r].getZ()-evap[q].getZ()));
					double fissbaraux = MassTable->getFisBarArray((tabelem[r].getA()-evap[q].getA()),(tabelem[r].getZ()-evap[q].getZ()));
					residuetab[q].Nucleus(path,Aaux,Zaux,massaux,shellaux,beta2aux,fissbaraux);
				}

				bool Test = false;
				double TimeDist = 0.;
				double DeltaT = 0.01;
				double Time = DeltaT;
				
				int nbT = 0;
				int nbTnext = nbT;
				if (path->getTime()==1) 
					nbTnext = nbT+1;

				bool StillDecay;
				do {
					
					for (int i=0;i<MaxI;i++) 
						gammatot[i] = 0.;
			 
					bool Continue = false;
					StillDecay = false;
					
					bool Delay = false; // Delay for taking into account the transient effect  
					if (path->getTime()==1 && Time<(path->getDelay()/getHbar())) Delay = true;
			 
					for (int q=0;q<Number;q++) 
						residuetab[q].InitSpectrum();

					Dec->WidthCall(path,&tabelem[r],residuetab,evap,gammatot,DeltaT,Delay,nbT); // to feed the daughter nucleus

			 		if (path->getTime()==1) {
						FissionTime += Dec->getFisSum()*DeltaT*(Time-DeltaT*0.5);
						TimeDist += Dec->getFisSum()*DeltaT;
					} 
			 
					for (int i=0;i<tabelem[r].getMaxi();i++) {
						if (gammatot[i]!=0.) {
							double factor = 0.;
							if (path->getTime()==1) 
								factor = exp(-gammatot[i]*DeltaT);
							for (int j=0;j<tabelem[r].getMaxj();j++) {
								double width = tabelem[r].getSpectrum(nbT,i,j)*factor+tabelem[r].getSpectrum(nbTnext,i,j);
								if (width>path->getCutWidth()) 
									tabelem[r].setSpectrum(nbTnext,i,j,width);
								if (width>path->getCutWidth() && !Continue) 
									Continue = true;
								if (path->getTime()==1 && width>path->getCutWidth() && !Test) 
									Test = true;
								tabelem[r].setSpectrum(nbT,i,j,0.);
							}
						}
					} // Decay to next time step
			 
					if (path->getTime()==1 && !Test && !Continue) Continue = true;
				
					for (int q=0;q<Number;q++) {
						Neu[q] = (source->getA()-source->getZ())-(residuetab[q].getA()-residuetab[q].getZ()); // number of emitted neutrons 
						Pro[q] = source->getZ()-residuetab[q].getZ(); // number of emitted protons
						// if (q==Number-1) cout << Dec->getBool(q) << endl;
						if (Dec->getBool(q)) { // if the spectrum of the residual nucleus is not empty
							if (Neu[q]<=path->getn() && Pro[q]<=path->getp()) {
								bool Already = false;
								for (int nb=0;nb<comp+num;nb++) {
									if (residuetab[q].getA()==tabelem[nb].getA() && residuetab[q].getZ()==tabelem[nb].getZ()) {
									// if the mother and daughter nuclei are the same in the case of gamma-ray emission	
										tabelem[nb] += residuetab[q];
										Already = true;
									}
								}
								if (!Already) {
									tabelem[comp+num].allocspectrum(residuetab[q].getMaxt(),residuetab[q].getMaxi(),residuetab[q].getMaxj());
									tabelem[comp+num] = residuetab[q];
									num += 1;
								}
							}
							if (q==(Number-1) && path->getGamma()==1) StillDecay = true; // for gamma decay
						//	cout << "q=" << q << " " << Number-1 << " " << StillDecay << endl;
						}
					} // end loop on the residual nuclei
		 
			 		if (path->getTime()==1) {
						nbT++;
						if (path->getTime()==1) nbTnext = nbT+1;
						else nbTnext = nbT;
						DeltaT *= path->getTimeFact();
						Time += DeltaT;
					}	  

					if (path->getTime()==1 && Continue) StillDecay = true;
					
					if (path->getTime()==1 && nbTnext>MaxT-1) StillDecay = false;
					
				} while (StillDecay);
		      
				double sum = 0.;
				for (int t=0;t<tabelem[r].getMaxt();t++)
					for (int i=0;i<tabelem[r].getMaxi();i++)
						for (int j=0;j<tabelem[r].getMaxj();j++)
							sum += tabelem[r].getSpectrum(t,i,j);	
				
				// Normalization to the population of the rth residual nucleus
				tabelem[r].setPop(sum);
				tabelem[r].deallocspectrum();	
			
			} // if (Ares==Asrc-rang)
		} // end loop on rth element
	  
		comp += num;
		rang++;
		
	} while (num!=0);
	
	for (int q=0;q<Number;q++) 
		residuetab[q].deallocspectrum();
}

void Cascade::Synth(Synthesis* synth) {
	
	int nbsynth = 0;
	for (int i=0;i<comp;i++) {
		synth[nbsynth].setA(tabelem[i].getA());
		synth[nbsynth].setZ(tabelem[i].getZ());
		synth[nbsynth].setTaux(tabelem[i].getPop());
		nbsynth +=1;
	}
	for (int i=0;i<comp;i++) tabelem[i].InitElement();
}

Cascade::Cascade(Readinput* path,int maxt,int maxi,int maxj) {
  
	n = path->getn(); // number of emitted neutrons 
	p = path->getp(); // number of emitted protons
	numpart = path->getNumber();
	gamma = path->getGamma();  
	Number = numpart+gamma; // from 1 to 4 types of emitted particles (n,p,a,gamma)	
	
	MaxI = maxi; 
	MaxT = maxt; 
	MaxJ = maxj;           

	comp = 0;
	FissionTime = 0.;
	TimeDist = 0.;

	Dec = new Decay(Number,MaxI,MaxJ);
	residuetab = new Element[Number];
	tabelem = new Element[(Number+1)*(n+1)*(p+1)];
	evap = new Particle[Number];
	gammatot = new double[MaxI];
	Neu = new int[Number];
	Pro = new int[Number];

	switch (numpart) {
		case 3: 
			evap[2].InitParticle(4,2);
		case 2:
			evap[1].InitParticle(1,1);
		case 1:
			evap[0].InitParticle(1,0);
			break;
		default :
			cerr << "Warning: Please check the input file for the number of evaporated particles." << endl;
			cin.get();
			exit(EXIT_FAILURE);
	}
	
	if (gamma==1) {
		evap[Number-1].setA(0);
		evap[Number-1].setZ(0);
	}
}

Cascade::~Cascade() {

	delete[] residuetab;
	delete[] tabelem;
	delete[] evap;
	delete[] gammatot;
	delete[] Neu;
	delete[] Pro;
	delete Dec; 
	
	residuetab = NULL;
	tabelem = NULL;
	evap = NULL;
	gammatot = NULL;
	Neu = NULL;
	Pro = NULL;
	Dec = NULL;
}
