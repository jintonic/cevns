static const double mm = 1;
static const double cm = 10*mm;
static const double cm2 = cm*cm;
static const double fermi = 1e-12*mm;

static const double  eV = 1e-6;
static const double keV = 1e-3;
static const double MeV = 1;
static const double GeV = 1e3;

static const double GF = 1.166364e-5/GeV/GeV;
static const double pi = 3.14159265358979323846;
static const double Qe = 1.602176487e-19; // electron charge in coulomb
static const double NA = 6.02214179e+23; // Avogadro constant
static const double hbarc= 197.327e-12*MeV*mm;
//______________________________________________________________________________
// F(q^2)
double F(double *x, double *parameter)
{
	double Er=x[0]*keV;
 	if (Er==0) return 1; // no momentum transfer
 	double M=parameter[0]; // GeV
 	double rn=parameter[1]; // fermi
	double q2=2*M*Er;
	double q = sqrt(q2);
	double r0= sqrt(rn*rn-5*fermi*fermi);
	double qr0 = q*r0/hbarc; // convert r0 from mm to natural unit MeV^-1
	double j1= (sin(qr0)-qr0*cos(qr0))/(qr0)/(qr0); // Levin & Smith, 1996
	double s2= fermi/hbarc * fermi/hbarc;
	return 3.*j1/qr0*exp(-q2*s2/2.);
}
//______________________________________________________________________________
// CEvNS differential cross section
double dXS(double Ev, double Er, unsigned short Z=53)
{
	unsigned short N = 74; // neutron number of I
	double M = 118.1846*GeV; // mass of I
	double R = 4.7500*fermi; // radius of I

	if (Z==11) { // Na
		N = 12;
		M = 21.4094*GeV;
		R = 2.9936*fermi;
	}

	// incident neutrinos must be energetic enough to cause a nuclear recoil
	double minEv = (Er+sqrt(2*M*Er))/2;
	if (Ev<=minEv) return 0;

	double sin2thetaw = 0.231;
	double qw = (N-Z) + 4*sin2thetaw * Z;

	double f;
	if (Er==0) {
		f=1; // no momentum transfer
	} else {
		double q2 = 2*M*Er;
		double q = sqrt(q2);
		double r0= sqrt(R*R-5*fermi*fermi);
		double qr0 = q*r0/hbarc; // convert r0 from mm to natural unit MeV^-1
		double j1= (sin(qr0)-qr0*cos(qr0))/(qr0)/(qr0); // Levin & Smith, 1996
		double s2= fermi/hbarc * fermi/hbarc;
		f = 3.*j1/qr0*exp(-q2*s2/2.);
	}

	double xs = GF*GF*M/8/pi * (1 + (1-Er/Ev)*(1-Er/Ev) - M*Er/Ev/Ev) *qw*qw*f*f;
	return xs*hbarc*hbarc; // 1/E^3 -> L^2/E
}
//______________________________________________________________________________
//
double dXSfunc(double *x, double *parameter)
{
	double Ev = x[0]*MeV; // neutrino energy
	double Er = x[1]*keV; // nuclear recoil energy
	// proton number of the target nucleus
	unsigned short Z = static_cast<unsigned short>(parameter[0]);
	return dXS(Ev, Er, Z)/(cm2/keV);
}
//______________________________________________________________________________
//
double dXSvsEr(double *x, double *parameter)
{
	double Er = x[0]*keV;
	double Ev = parameter[0]*MeV;
	unsigned short Z = static_cast<unsigned short>(parameter[1]);
	return dXS(Ev, Er, Z)/(1e-38*cm2/keV);
}

TF1 *fdXSvsEr = new TF1("fdXSvsEr", dXSvsEr, 0, 10, 2); // Er in keV
//______________________________________________________________________________
//
double dXSvsEv(double *x, double *parameter)
{
	double Ev = x[0]; fdXSvsEr->SetParameter(0, Ev);
	double Z = parameter[0]; fdXSvsEr->SetParameter(1, Z);
	if (Z==11) return fdXSvsEr->Integral(0, 2*Ev*Ev/21); 
	else return fdXSvsEr->Integral(0, 2*Ev*Ev/118); 
}
//______________________________________________________________________________
//
TF1 *fdNdE = new TF1("fdNdE", "exp(1.39575-0.549369*x+0.00136642*x*x"
			"-0.00665779*x*x*x)", 0,10);
//______________________________________________________________________________
//
double dXSxdNv(double *x, double *parameter)
{
	return dXSvsEv(x,parameter)*fdNdE->Eval(x[0])/MeV;
}
//______________________________________________________________________________
//
double dXSxdNvsEv(double *x, double *parameter)
{
	double Ev = x[0]*MeV;
	double Er = parameter[0]*keV;
	unsigned short Z = static_cast<unsigned short>(parameter[1]);
	return dXS(Ev, Er, Z)*fdNdE->Eval(x[0]); // L2/E * 1/E
}

TF1 *fdXSxdNvsEv = new TF1("fdXSxdNvsEv", dXSxdNvsEv, 0, 10/*MeV*/, 2);
//______________________________________________________________________________
//
TF1 *fEffI = new TF1("fEffI", "1-exp(-0.293746*x-0.1299*x*x)",0,10);
TF1 *fEffNa = new TF1("fEffNa", "1-exp(-0.596213*x-0.510011*x*x)",0,10);
//______________________________________________________________________________
//
double dXSxdNxEff(double *x, double *parameter)
{
	fdXSxdNvsEv->SetParameter(0, x[0]); // x[0] = Er in keV
	double Z = parameter[0]; fdXSxdNvsEv->SetParameter(1, Z);
	double eff = fEffI->Eval(x[0]);
	if (Z==11) eff = fEffNa->Eval(x[0]);
	double power=3400; // MW
	double Efission235=202; // MeV
	double Rfission235=power/Efission235/Qe; // # of fissions per second
	double mass=10000; // 10 kg
	double Nnuclei=mass/(23+127)*NA; // # of Na or I
	double distance=10000*mm;
	double duration=365*24*3600;
	double coefficient=Rfission235*duration*Nnuclei/4./pi/distance/distance;
	return fdXSxdNvsEv->Integral(0, 10) /*L2/E*/ * keV * eff * coefficient;
}
//______________________________________________________________________________
//
double dXSxdNdE(double *x, double *parameter)
{
	return dXSfunc(x,parameter) * fdNdE->Eval(x[0]);
}
//______________________________________________________________________________
//
double dXSxdNdExEff(double *x, double *parameter)
{
	if (parameter[0]==11) return dXSxdNdE(x,parameter) * fEffNa->Eval(x[1]);
	else return dXSxdNdE(x,parameter) * fEffI->Eval(x[1]);
}
//______________________________________________________________________________
//
TF2 *fdXSxdNdExEff = new TF2("fdXSxdNdE", dXSxdNdExEff,
			0, 10, // neutrino energy range in MeV
			0, 10, // nuclear recoil energy range in keV
			1); // 1 parameter
//______________________________________________________________________________
//
void drawReactorNeutrinoSpectra()
{
	// set drawing style
	gROOT->SetStyle("Plain"); // pick up an existing style to modify
	gStyle->SetMarkerStyle(kFullDotLarge);
	gStyle->SetMarkerSize(0.8);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFont(132);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetFitFormat(".2f"); // shorten fit result
	gStyle->SetStatX(0.99); // top right corner x
	gStyle->SetStatY(0.99); // top right corner y
	gStyle->SetStatFont(132);
	gStyle->SetLabelFont(132,"XYZ");
	gStyle->SetTitleFont(132,"H");
	gStyle->SetTitleFont(132,"XYZ");
	gStyle->SetLabelSize(0.05,"XYZ");
	gStyle->SetTitleSize(0.05,"XYZ");
	gStyle->SetTitleOffset(1.0,"Y");
	gStyle->SetTitleOffset(1.2,"Z");
	gStyle->SetPadRightMargin(0.01);
	gStyle->SetPadLeftMargin(0.1);
	gStyle->SetPadTopMargin(0.01);
	gStyle->SetPadBottomMargin(0.11);
	gStyle->SetTitleAlign(23);
	gStyle->SetTitleX(0.52);
	gStyle->SetTitleY(0.96);

	//____________________________________________________________________________
	// draw F(q^2)
	TF1 *fF = new TF1("fF", F, 0, 5 /* keV */, 2);
	fF->SetParameter(0,118.1846*GeV);
	fF->SetParameter(1,4.7500*fermi);
	TF1 *fI = fF->DrawCopy();
	fI->SetTitle(";E_{recoil};F(q^{2})");
	fF->SetParameter(0,21.4094*GeV);
	fF->SetParameter(1,2.9936*fermi);
	fF->SetLineColor(kRed);
	fF->Draw("same");
	
	TLegend *l1 = new TLegend(0.75, 0.6, 0.9, 0.78);
	l1->AddEntry(fF, "Na", "l");
	l1->AddEntry(fI, "I", "l");
	l1->Draw();

	gPad->SetGridx(); gPad->SetGridy();
	gPad->Print("F.png");

	//____________________________________________________________________________
	// draw reactor neutrino spectra normalized per fission
	// * data from Fig. 3, Hayes & Vogel, Ann. Rev. Nucl. Part. Sci. 66 (2016) 219
	// * formula from Sec. 4.1.1, Bernstein, et al., Sci. Glob. Sec. 18 (2010) 127

	// mapping between pixel and Ev in MeV
	const int np=4; double energy[np]={2, 4, 6, 8}, pixel[np]={133, 405, 675, 946};
	TGraph *gmp = new TGraph(np, pixel, energy);
	gmp->SetTitle(";pixel;E_{#nu} [MeV]");
	gmp->Draw("ap");

	TF1 *fmp = new TF1("fmp","pol1",100,1000);
	gmp->Fit("fmp");
	gPad->Print("EvVSpixel.png");

	// reproduce Fig. 3 in Hayes and Vogel
	const int n=11; double Ev5[n], Ev8[n], Ev9[n], Ev1[n];
	double dNdE[n]={0.003, 0.005, 0.01, 0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 1, 2};
	double pixel5[n]={895, 856, 813, 747, 647, 556, 404, 319, 257, 186, 33};
	double pixel8[n]={971, 930, 871, 815, 711, 627, 464, 377, 315, 232, 62};
	double pixel9[n]={848, 814, 749, 683, 574, 486, 358, 270, 219, 156, 22};
	double pixel1[n]={872, 841, 790, 726, 630, 544, 407, 325, 262, 192, 33};
	for (int i=0; i<n; i++) {
		Ev5[i] = fmp->Eval(pixel5[i]);
		Ev8[i] = fmp->Eval(pixel8[i]);
		Ev9[i] = fmp->Eval(pixel9[i]);
		Ev1[i] = fmp->Eval(pixel1[i]);
	}
	TGraph *gne5 = new TGraph(n, Ev5, dNdE);
	TGraph *gne8 = new TGraph(n, Ev8, dNdE);
	TGraph *gne9 = new TGraph(n, Ev9, dNdE);
	TGraph *gne1 = new TGraph(n, Ev1, dNdE);

	gne5->SetTitle(";E_{#nu} [MeV]; dN_{#nu}/dE_{#nu}");
	gne5->SetMarkerColor(kOrange); gne5->SetLineColor(kOrange);
	gne8->SetMarkerColor(kAzure+1); gne8->SetLineColor(kAzure+1);
	gne9->SetMarkerColor(8); gne9->SetLineColor(8);
	gne1->SetMarkerColor(kRed); gne1->SetLineColor(kRed);

	gne5->Draw("apc");
	gne8->Draw("pc");
	gne9->Draw("pc");
	gne1->Draw("pc");

	// empirical formula from Bernstein, et al.
	TF1 *f5 = new TF1("f5", "exp(-([0]+[1]*x+[2]*x*x+[3]*x*x*x))",1,9);
	TF1 *f8 = new TF1("f8", "exp(-([0]+[1]*x+[2]*x*x+[3]*x*x*x))",1,9);
	TF1 *f9 = new TF1("f9", "exp(-([0]+[1]*x+[2]*x*x+[3]*x*x*x))",1,9);
	TF1 *f1 = new TF1("f1", "exp(-([0]+[1]*x+[2]*x*x+[3]*x*x*x))",1,9);

	f5->SetLineColor(kOrange);
	f8->SetLineColor(kAzure+1);
	f9->SetLineColor(8);
	f1->SetLineColor(kRed);

	gne5->Fit("f5");
	gne8->Fit("f8");
	gne9->Fit("f9");
	gne1->Fit("f1");

	TLegend *leg = new TLegend(0.8, 0.6, 0.98, 0.98);
	leg->AddEntry(f5, "^{235}U", "l");
	leg->AddEntry(f8, "^{238}U", "l");
	leg->AddEntry(f9, "^{239}Pu", "l");
	leg->AddEntry(f1, "^{241}Pu", "l");
	leg->Draw();

	gPad->Print("dN_VS_dE.png");
	gPad->SetLogy();
	gPad->Print("logdN_VS_dE.png");

	//____________________________________________________________________________
	// draw CEvNS cross sections
	TF2 *fdXS = new TF2("fdXS", dXSfunc,
			0, 10, // neutrino energy range in MeV
			0, 10, // nuclear recoil energy range in keV
			1); // 1 parameter

	gPad->SetLogz();
	gPad->SetRightMargin(0.17);
	gPad->SetTopMargin(0.02);

	fdXS->SetParameter(0,11);
	fdXS->Draw("colz");
	fdXS->SetTitle(";E_{#nu} [MeV]; E_{r} [keV];"
			"d#sigma/dE_{r} [cm^{2}/keV]");
	gPad->Print("dXS.png");

	TF2 *fdXSxdNdE = new TF2("fdXSxdNdE", dXSxdNdE,
			0, 10, // neutrino energy range in MeV
			0, 10, // nuclear recoil energy range in keV
			1); // 1 parameter
	fdXSxdNdE->SetParameter(0,11);
	fdXSxdNdE->Draw("colz");
	fdXSxdNdE->SetTitle(";E_{#nu} [MeV];E_{r} [keV];"
			"d#sigma/dE_{r} [cm^{2}/keV] #times dN_{#nu}/dE_{#nu} [1/MeV]");
	gPad->Print("dXSxdNdE.png");

	fdXSxdNdExEff->SetParameter(0,11);
	fdXSxdNdExEff->Draw("colz");
	fdXSxdNdExEff->SetTitle(";E_{#nu} [MeV];E_{r} [keV];"
			"d#sigma/dE_{r} [cm^{2}/keV] #times dN_{#nu}/dE_{#nu} [1/MeV] #times #epsilon");
	gPad->Print("dXSxdNdExEff.png");

	gPad->Clear();
	gPad->SetLeftMargin(0.11);
	gPad->SetRightMargin(0.02);

	TF1 *fdXSvsEv = new TF1("fdXSvsEv", dXSvsEv, 0, 8 /* MeV */, 1);
	fdXSvsEv->SetParameter(0,53);
	TF1 *fc = fdXSvsEv->DrawCopy();
	fc->GetYaxis()->SetRangeUser(5e-6, 5e-1);
	fc->SetTitle(";E_{#nu} [MeV];#sigma(E_{#nu}) [10^{-38} cm^{2}]");
	fc->GetYaxis()->SetTitleOffset(1.1);

	fdXSvsEv->SetParameter(0,11);
	fdXSvsEv->SetLineColor(kRed);
	fdXSvsEv->Draw("same");

	l1->Draw();
	gPad->Print("dXSvsEv.png");

	//____________________________________________________________________________
	// draw (XS x Ev spectra)
	gPad->Clear();
	gPad->SetLogy(0);
	gPad->SetTopMargin(0.06);

	TF1 *fdXSxdNv = new TF1("fdXSxdNv", dXSxdNv, 0,10, 1);
	TF1 *fco = fdXSxdNv->DrawCopy();
	fco->SetTitle(";E_{#nu} [MeV];"
			"#sigma(E_{#nu}) #times dN_{#nu}/dE_{#nu} [arbitury unit]");
	fco->GetYaxis()->SetTitleOffset(1.1);

	fdXSxdNv->SetParameter(0,11);
	fdXSxdNv->SetLineColor(kRed);
	fdXSxdNv->Draw("same");

	l1->Draw();
	gPad->Print("dXSxdNv.png");

	//____________________________________________________________________________
	// draw (XS x Ev spectra x eff)
	TF1 *fdXSxdNxEff = new TF1("fdXSxdNxEff", dXSxdNxEff, 0, 5, 1);
	fdXSxdNxEff->SetParameter(0,11);

	gStyle->SetTitleY(1);
	gPad->SetLeftMargin(0.13);
	fdXSxdNxEff->Draw();
	fdXSxdNxEff->GetYaxis()->SetTitleOffset(1.3);
	//fdXSxdNxEff->GetYaxis()->SetRangeUser(0,21);
	fdXSxdNxEff->SetTitle("3.4 GW, 10 kg NaI, 10 m from core;E^{Na}_{recoil} [keV];"
			"Events / year / keV");
	//fdXSxdNxEff->GetHistogram()->Draw("hist");
	gPad->Print("dXSxdNxEffNa.png");
	cout<<fdXSxdNxEff->Integral(1.3,5)<<endl;

	fdXSxdNxEff->SetParameter(0,53);
	fdXSxdNxEff->Draw();
	fdXSxdNxEff->GetYaxis()->SetTitleOffset(1.3);
	fdXSxdNxEff->SetTitle("3.4 GW, 10 kg NaI, 10 m from core;E^{I}_{recoil} [keV];"
			"Events / year / keV");
	//fdXSxdNxEff->GetHistogram()->Draw("hist");
	gPad->Print("dXSxdNxEffI.png");
	cout<<fdXSxdNxEff->Integral(2.5,5)<<endl;

	double power=3400; // MW
	double Efission235=202; // MeV
	double Rfission235=power/Efission235/Qe; // # of fissions per second
	double mass=10000; // 10 kg
	double Nnuclei=mass/(23+127)*NA; // # of Na or I
	double distance=10000*mm/cm;
	double duration=365*24*3600;
	double coefficient=Rfission235*duration*Nnuclei/4./pi/distance/distance;

	fdXSxdNdExEff->SetParameter(0,11);
	cout<<"Rfission235: "<<Rfission235<<endl;
	cout<<"N_Na: "<<Nnuclei<<endl;
	cout<<"365*24*3600: "<<duration<<endl;
	cout<<"1/4/pi/1000/1000 cm^-2: "<<1/4./pi/distance/distance<<endl;
	cout<<"N of CEvNS on Na: "<<fdXSxdNdExEff->Integral(3.5,9, 1,6) * coefficient<<endl;
}
