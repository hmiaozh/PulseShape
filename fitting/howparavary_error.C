#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "RooGlobalFunc.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TRandom.h"
#include <string>
#include <vector>
#include <math.h>

//parameter vary with different charge, he hb separately


Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t yy;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         yy = x[0]-xx;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(yy,0.0,par[3]);

         xx = xupp - (i-.5) * step;
         yy = x[0]-xx;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(yy,0.0,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}



TF1 *langaufit(TProfile *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

//   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
//   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");

   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function

}

/*
Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

   // Seaches for the location (x value) at the maximum of the
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.

   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 10000;


   // Search for maximum

   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = langaufun(&x,params);

      if (l < lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-1);

   maxx = x;

   fy = l/2;


   // Search for right x location of fy

   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-2);

   fxr = x;


   // Search for left x location of fy

   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-3);


   fxl = x;

   FWHM = fxr - fxl;
   return (0);
}
*/

void setDataPulse(TGraphErrors gr, double *dataP) {
    for (Int_t i=0; i<10; i++) {
      double x,y;
      gr.GetPoint(i,x,y);
      dataP[i]=y;
    }
}

void setDataPulse2(TProfile pr, double *dataPu){
   for (Int_t i=1; i<11; i++) {
      double y;
      y=pr.GetBinContent(i);
      dataPu[i-1]=y;
    }
}

void howparavary_error() {
   // Fill Histogram
   

//   TFile* f1 = TFile::Open("/afs/cern.ch/work/h/hum/public/CMSSW_8_0_8/src/PulseShapes/thegraphs2.root");
  TFile* f1 = TFile::Open("/afs/cern.ch/work/h/hum/public/CMSSW_8_0_1/src/PulseShapes/Ped/datall.root");
  TFile* f2 = new TFile("/afs/cern.ch/work/h/hum/public/CMSSW_8_0_1/src/PulseShapes/fitting/LG_para_data_try.root","recreate");
  
  TGraphErrors* Chi2_c_hb = new TGraphErrors(58);
  Chi2_c_hb->SetName("Chi2_c_hb");
  TGraphErrors* LWidth_c_hb = new TGraphErrors(58);
  LWidth_c_hb->SetName("LWidth_c_hb");
  TGraphErrors* MP_c_hb = new TGraphErrors(58);
  MP_c_hb->SetName("MP_c_hb");
  TGraphErrors* Gsigma_c_hb = new TGraphErrors(58);
  Gsigma_c_hb->SetName("Gsigma_c_hb");
//  TGraphErrors* Area_c_hb = new TGraphErrors(58);
  
  TGraphErrors* Chi2_c_he = new TGraphErrors(58);
  Chi2_c_he->SetName("Chi2_c_he");
  TGraphErrors* LWidth_c_he = new TGraphErrors(58);
  LWidth_c_he->SetName("LWidth_c_he");
  TGraphErrors* MP_c_he = new TGraphErrors(58);
  MP_c_he->SetName("MP_c_he");
  TGraphErrors* Gsigma_c_he = new TGraphErrors(58);
  Gsigma_c_he->SetName("Gsigma_c_he");
//  TGraphErrors* Area_c_he = new TGraphErrors(58);

  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
  
  pllo[0]=0.05; pllo[1]=3.5; pllo[2]=1.0; pllo[3]=0.1;
  plhi[0]=0.5; plhi[1]=5.0; plhi[2]=1.0; plhi[3]=0.8;
  sv[0]=0.3; sv[1]=4.0; sv[2]=1.0; sv[3]=0.3;
  Double_t chisqr;
  Int_t    ndf;


  for(Int_t m=0; m<58; m++){ 
   char gname[50];
   char pname[50];
   sprintf(gname, "hb_prof_%i",m);
   sprintf(pname, "he_prof_%i",m);
  
   Int_t n = m*10+20;
 
   TProfile *ro = (TProfile*)f1->Get(gname); 
   TProfile *ra = (TProfile*)f1->Get(pname); 

//   TGraphErrors *g = (TGraphErrors*)f1->Get("q_100_110");

//   double dataPulse[10];
//   setDataPulse2(*g, &dataPulse[0]);

//   TH1F *ro = new TH1F("ro","pulseshape",10,-0.5,9.5);
//   for(Int_t i=0; i < 10; i++) {
//       cout << dataPulse[i] << endl;
//       ro->Fill(i,dataPulse[i]);
//       double j = (double)i;
//       ro->Fill(j,g->GetBinContent(j+1.0));
//       ro->SetBinContent(i,g->GetBinContent(i));
//   }
   
   // Fitting SNR histo
//   printf("Fitting...\n");

   // Setting fit range and start values
   fr[0]=0.3*ro->GetMean();
   fr[1]=3.0*ro->GetMean();

//   Double_t chisqr;
//   Int_t    ndf;
   TF1 *fitsnr = langaufit(ro,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  
   Chi2_c_hb->SetPoint(m,n,chisqr);

   LWidth_c_hb->SetPoint(m,n,fp[0]);
   LWidth_c_hb->SetPointError(m,0,fpe[0]);

   MP_c_hb->SetPoint(m,n,fp[1]);
   MP_c_hb->SetPointError(m,0,fpe[1]);

   Gsigma_c_hb->SetPoint(m,n,fp[3]); 
   Gsigma_c_hb->SetPointError(m,0,fpe[3]);

//   Area_c_hb->Fill(n,fp[2]);


   fr[0]=0.3*ra->GetMean();
   fr[1]=3.0*ra->GetMean();

   TF1 *fitsnr1 = langaufit(ra,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);

   Chi2_c_he->SetPoint(m,n,chisqr);
   LWidth_c_he->SetPoint(m,n,fp[0]);
   LWidth_c_he->SetPointError(m,0,fpe[0]);
   MP_c_he->SetPoint(m,n,fp[1]);
   MP_c_he->SetPointError(m,0,fpe[1]);
   Gsigma_c_he->SetPoint(m,n,fp[3]);
   Gsigma_c_he->SetPointError(m,0,fpe[3]);
//   Area_c_he->Fill(n,fp[2]);

 }

 Chi2_c_hb->Write();
 Chi2_c_he->Write();
 LWidth_c_hb->Write();
 LWidth_c_he->Write();
 MP_c_hb->Write();
 MP_c_he->Write();
 Gsigma_c_hb->Write();
 Gsigma_c_he->Write();



  f2->Write();
  f2->Close();
   
//   Double_t SNRPeak, SNRFWHM;
//   langaupro(fp,SNRPeak,SNRFWHM);

//   printf("Fitting done\nPlotting results...\n");

   // Global style settings
/*   gStyle->SetOptStat(1111);
   gStyle->SetOptFit(111);
   gStyle->SetLabelSize(0.03,"time slice");
   gStyle->SetLabelSize(0.03,"y");

   ro->GetXaxis()->SetRange(0,10);
   ro->Draw("");
   fitsnr->Draw("lsame");
*/
}
