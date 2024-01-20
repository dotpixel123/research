#include <Writeoutput.h>

void Writeoutput::write(double value1) {
	macro << fixed << scientific << value1 << "   " ;
}

void Writeoutput::writeendl() {
	macro << endl;
}

void Writeoutput::writestr(string phrase) {
	macro << phrase << "   ";
}

void Writeoutput::writepara(Readinput* path) {

	macro << "# Parameters and options used in calculation" << endl;
	macro << "# " << path->getZpro() << "   Projectile charge" << endl;
	macro << "# " << path->getApro() << "   Projectile mass" << endl;
	macro << "# " << path->getZtar() << "   Target charge" << endl;
	macro << "# " << path->getAtar() << "   Target mass" << endl;
	macro << "# " << path->getShell() << "   Shell correction" << endl;
	macro << "# " << path->getShellFactor() << "   Shell-correction factor" << endl;
	macro << "# " << path->getKinStep() << "   Stepsize for the kinetic energy" << endl;
	macro << "# " << path->getLevDensPar() << "   Level-density parameter used" << endl;
	macro << "# " << path->getCstShell() << "   Damping-shell energy used" << endl;
	macro << "# " << path->getafan() << "   af/an ratio" << endl;
	macro << "# " << path->getFacInerMom() << "   Correction factor for the moment of inertia" << endl;
	macro << "# " << path->getMollerTable() << "   Mass table used" << endl;
	macro << "# " << path->getPair() << "   Pairing correction" << endl;
	macro << "# " << path->getKramers() << "   Kramers-Strutinsky factor" << endl;
	macro << "# " << path->getCollEnhFac() << "   Collective enhancement factor" << endl;
	macro << "# " << path->getFriction() << "   Reduced friction used" << endl;
	macro << "# " << path->getOmegags() << "   Potential curvature in the ground state (MeV)" << endl;
	macro << "# " << path->getOmegasd() << "   Potential curvature at the saddle point (MeV)" << endl;
	macro << "# " << path->getEva() << "   Evaporation model used" << endl;
	macro << "# " << path->getFisTrans() << "   Hill-Wheeler penetration factor used" << endl;
	macro << "# " << path->getFisBar() << "   Fission-barrier model used" << endl;
	macro << "# " << path->getEmin() << "   Minimum excitation energy" << endl;
	macro << "# " << path->getEmax() << "   Maximum excitation energy" << endl;
	macro << "# " << path->getEstep() << "   Stepsize for the spectral discretization" << endl;
	macro << "# " << path->getEview() << "   Stepsize for plotting the excitation function" << endl;
	macro << "# " << path->getJmin() << "   Minimum angular momentum" << endl;
	macro << "# " << path->getJmax() << "   Maximum angular momentum" << endl;
	macro << "# " << path->getn() << "   Number of emitted neutrons" << endl;
	macro << "# " << path->getp() << "   Number of emitted protons" << endl;
	macro << "# " << path->getNumber() << "   Types of emitted particles at each step" << endl;
	macro << "# " << path->getGamma() << "   Gamma-ray emission" << endl;
	macro << "# " << path->getGammaModel() << "   Gamma-ray emission model used" << endl;
	macro << "# " << path->getTime() << "   Fission-time calculation" << endl;
	macro << "# " << path->getTimeFact() << "   Time factor" << endl;
	macro << "# " << path->getCutTcoef() << "   Cut-off for the transmission coefficient" << endl;
	macro << "# " << path->getCutWidth() << "   Cut-off for the decay width" << endl;
	macro << "# " << path->getFusion() << "   Fusion model used" << endl;
	macro << "# " << path->getDeltaB0()<< ", " << path->getDeltaBvar() << "   Barrier-distribution parameters" << endl;
	macro << "# " << path->getNumberPoint() << "   Number of abscissas for the G-L integration" << endl;
	macro << "# " << path->getFusInput() << "   Name of the fusion data file" << endl;
	macro << "# " << path->getMassTable() << "   Name of the nuclear data file" << endl;
	macro << "# " << path->getFisBarTable() << "   Name of the fission-barrier data file" << endl;
	macro << "# " << path->getOutput() << "   Name of the output file" << endl;
		
}

Writeoutput::Writeoutput(string phrase) {

	phrase = "../Results/" + phrase;

	macro.open(phrase.c_str());
	if (macro.fail()) {
		cerr<<"Warning: Cannot open the output file! Please check it."<<endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
}

Writeoutput::~Writeoutput() {
	macro.close();
}
