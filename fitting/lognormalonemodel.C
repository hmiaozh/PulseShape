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

/*
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

      return (par[2] *  step * sum * invsq2pi / par[3]);
}

*/

TF1 *lognormalfit(TProfile *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

//   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
//   if (ffitold) delete ffitold;

//   TF1 *ffit = new TF1(FunName,"TMath::LogNormal(x,[0],[1],[2])",fitrange[0],fitrange[1],3);
   TF1 *ffit = new TF1(FunName,"TMath::LogNormal(x,[0],[1],[2])",fitrange[0],fitrange[1]);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Sigma","Theta","m");

   for (i=0; i<3; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<3; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function

}


void lognormalonemodel() {
   // Fill Histogram
   

//   TFile* f1 = TFile::Open("/afs/cern.ch/work/h/hum/public/CMSSW_8_0_8/src/PulseShapes/thegraphs2.root");
   TFile* f1 = TFile::Open("/afs/cern.ch/work/h/hum/public/CMSSW_8_0_1/src/PulseShapes/Ped/datall.root");
   TProfile *ro = (TProfile*)f1->Get("hb_prof_32");  



   // Fitting SNR histo
   printf("Fitting...\n");

   // Setting fit range and start values
   Double_t fr[2];
   Double_t sv[3], pllo[3], plhi[3], fp[3], fpe[3];
//   fr[0]=0.3*ro->GetMean();
//   fr[1]=3.0*ro->GetMean();

   fr[0]=0.0;
   fr[1]=9.0;

   pllo[0]=0.1; pllo[1]=0.0; pllo[2]=0.5;
   plhi[0]=1.0; plhi[1]=5.0; plhi[2]=2.0;
   sv[0]=0.5; sv[1]=4.0; sv[2]=1.0;

   Double_t chisqr;
   Int_t    ndf;
   TF1 *fitsnr = lognormalfit(ro,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);

   Double_t differ[7]={0.0};
   Double_t xv[7]={0.0};   


   for (Int_t i = 3; i<10; i++){
        xv[i-3] = (double)i;
        Double_t funcv = TMath::LogNormal(xv[i-3], fp[0],fp[1],fp[2]);
        Double_t pulv = ro->GetBinContent(i+1);
	differ[i-3] = (funcv - pulv)/pulv;
        cout << i << " pulse value " << pulv << ",function value " << funcv << ", difference is " << differ[i-3] << endl;   
   }

     


//   Double_t SNRPeak, SNRFWHM;
//   langaupro(fp,SNRPeak,SNRFWHM);

   printf("Fitting done\nPlotting results...\n");

/*   TF1 *gau = new TF1 ("gau","0.3989422804014 * TMath::Gaus(x,[0],[1])/[1]",0.,10.);
   gau->SetParameters(0.0,fp[3]);
   TF1 *lan = new TF1 ("lau","TMath::Landau(x,[0],[1])/[1]",0.,10.);
   lan->SetParameters(fp[1],fp[0]);
*/
   // Global style settings
   gStyle->SetOptStat(1111);
   gStyle->SetOptFit(111);
   gStyle->SetLabelSize(0.03,"time slice");
   gStyle->SetLabelSize(0.03,"y");

   ro->GetXaxis()->SetRange(0,10);
   ro->GetYaxis()->SetRangeUser(0,0.8);
   ro->Draw("C");
   fitsnr->Draw("lsame");
//   gau->SetLineColor(1); 
//   lan->SetLineColor(3);  
//   gau->Draw("same");
//   lan->Draw("same");    
}
