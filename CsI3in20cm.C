double Lindhard(double *x, double *p)
{
   double z = p[0];
   double a = p[1];
   double keVee = x[0];

   double k = 0.133*pow(z,2./3)/sqrt(a);
   double epsilon = 11.5*keVee/pow(z,7./3);
   double g = 3*pow(epsilon,0.15) + 0.7*pow(epsilon,0.6) + epsilon;
   return k*g/(1+k*g);
}

const int nbin = 200; // max number of bins in the Poisson histogram
const int nexp = 10000; // max number of experiments to be performed

double Eff2PE(double keVee=0.1, double yield=30/*PE/keV*/, bool drawHist=false)
{
   double npe = yield*keVee;

   gStyle->SetOptStat(false);
   TH1F *hc, *htc;
   TH1F h("h",Form("LY=%.0fPE/keV;Number of PE;Entries",yield),nbin,0,nbin);
   for (int i=0; i<nexp; i++) h.Fill(gRandom->Poisson(npe));
   if (drawHist) {
      hc = (TH1F*) h.DrawClone();
      hc->GetXaxis()->SetRangeUser(0,18);
   }

   TH1F htrg("htrg","",nbin,0,nbin);
   for (int ipe=2; ipe<nbin; ipe++) {
      int ntrg = 0; // number of triggered events (both PMTs have hits)
      // n = number of experiments, in which ipe are observed
      int n = h.GetBinContent(ipe+1);
      for (int j=0; j<n; j++) {
         int n1=0, n2=0; // number of PE in PMT 1 & 2
         for (int idx=0; idx<ipe; idx++) { // idx of PE
            if (gRandom->Gaus(58.26,23.69)>30) { // above electronic noise
               if (gRandom->Rndm()<0.5) n1++;
               else n2++;
            }
         }
         if (n1!=0 && n2!=0) ntrg++;
      }
      htrg.SetBinContent(ipe+1,ntrg);
   }
   htrg.SetLineColor(kRed);
   if (drawHist) {
      htc = (TH1F*) htrg.DrawClone("same");
      TLegend *l = new TLegend(0.4,0.7,0.88,0.88);
      l->AddEntry(hc,"No coincidence","l");
      l->AddEntry(htc,"Requiring coincidence","l");
      l->Draw();
   }

   double rate = htrg.Integral()/h.GetEntries();
   cout<<"trigger rate @ "<<keVee<<" keVee: "<<rate*100<<"%"<<endl;
   return rate;
}

void CsI3in20cm()
{
   TF1 *qfNa = new TF1("qfNa",Lindhard,0,10,2);
   qfNa->SetParameters(11, 23);
   TF1 *qfI = new TF1("qfI", Lindhard, 0, 5, 2);
   qfI->SetParameters(53, 127);
   const int n = 120;
   double effee[n], effNa[n], effI[n], keV[n];
   for (int i=0; i<n; i++) {
      keV[i] = 0.01*(i+1);
      effee[i] = Eff2PE(keV[i]);
      effNa[i] = Eff2PE(keV[i]*qfNa->Eval(keV[i]));
      effI[i] = Eff2PE(keV[i]*qfI->Eval(keV[i]));
   }
   TGraph *g = new TGraph(n,keV,effee);
   TGraph *gNa = new TGraph(n,keV,effNa);
   TGraph *gI = new TGraph(n,keV,effI);
   g->SetTitle("");
   g->GetXaxis()->SetTitle("Recoil energy [keV]");
   g->GetYaxis()->SetTitle("Trigger efficiency (2 PE, 2 SiPM coincidence)");
   g->Draw("ap");
   gNa->SetMarkerColor(kRed);
   gNa->Draw("p");
   gI->SetMarkerColor(kBlue);
   gI->Draw("p");
   gPad->SetGridx(); gPad->SetGridy();

   TLegend *l = new TLegend(0.7,0.2,0.88,0.4);
   l->AddEntry(g,"Electron","p");
   l->AddEntry(gNa,"Na","p");
   l->AddEntry(gI,"I","p");
   l->Draw();

   gPad->Print("eff.png");
}
