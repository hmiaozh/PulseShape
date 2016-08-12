#include "TROOT.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"
#include "RooGlobalFunc.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TRandom.h"



#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>

#include "HcalPulseShapes.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using namespace std;

//static const double tzero[3]= {23.960177, 13.307784, 9.109694};
//static const double slope[3] = {-3.178648,  -1.556668, -1.075824 };
//static const double tmax[3] = {16.00, 10.00, 6.25 };
//
//double delay(double fC, BiasSetting bias) {
//  double rawDelay=tzero[bias]+slope[bias]*log(fC);
//  return (rawDelay<0)?(0):((rawDelay>tmax[bias])?(tmax[bias]):(rawDelay));   
//}


void Analy() {
  
//  TFile* f1 = TFile::Open("/afs/cern.ch/work/h/hum/public/CMSSW_8_0_1/src/PulseShapes/Ped/datall.root");
//  TFile* f1 = TFile::Open("/afs/cern.ch/work/h/hum/public/CMSSW_8_0_1/src/PulseShapes/NewMC/newmc.root");
  TFile* f1 = TFile::Open("/afs/cern.ch/work/h/hum/public/CMSSW_8_0_1/src/PulseShapes/digishiftMCsample/digi_shift0_PulseShape.root");
//  TFile* f2 = new TFile("/afs/cern.ch/work/h/hum/public/CMSSW_8_0_1/src/PulseShapes/fitting/Results/fixed8para_data_shift0.root","recreate");
  TFile* f2 = new TFile("/afs/cern.ch/work/h/hum/public/CMSSW_8_0_1/src/PulseShapes/fitting/Results/8para_digi_shift0.root","recreate");

  vector<TProfile*> tsc_hb;
//  vector<TProfile*> tsc_he;
//  vector<TProfile*> tsc_15;
//  vector<TProfile*> tsc_10;

  char tname[50];
//  char sname[50];
//  char aname[50];
//  char bname[50];
  for (UInt_t i=3; i<10; i++) {
     sprintf(tname,"ts%i_c_hb",i);
//     sprintf(sname,"ts%i_c_he",i);
//     sprintf(aname,"ts%i_c_15",i);
//     sprintf(bname,"ts%i_c_10",i);
     tsc_hb.push_back(new TProfile(tname, "",28,15.0,305.0));
//     tsc_hb.push_back(new TProfile(tname, "",58,15.0,605.0));
//     tsc_he.push_back(new TProfile(sname, "",58,15.0,605.0));
//     tsc_15.push_back(new TProfile(aname, "",58,15.0,605.0));
//     tsc_10.push_back(new TProfile(bname, "",58,15.0,605.0));

  }


  float ts1, ts2, ts3, thpd, tpre, wd1, wd2, wd3;
  ts1=8.0; ts2=19.0; ts3=29.3; thpd=4.0; tpre=9.0; wd1=2.0; wd2=0.7; wd3=0.32;  //105
//  ts1=8.0; ts2=10.0; ts3=22.3; thpd=4.0; tpre=7.0; wd1=2.0; wd2=0.7; wd3=1.0;  //125
  HcalPulseShapes thePulses_;
  
 
  thePulses_.computeHPDShape(ts1,ts2,ts3,thpd,tpre,wd1,wd2,wd3, thePulses_.hpdShape_);


  for(Int_t m=0; m<28; m++){
//  for(Int_t m=0; m<58; m++){
   char gname[50];
   sprintf(gname, "hb_prof_%i",m);

   TProfile *ro = (TProfile*)f1->Get(gname);
  
   Int_t n = m*10+20;
 
   Double_t differ[7]={0.0};

   Double_t pul = ro->GetBinContent(5);

   Double_t rawtslew=13.307784-1.556668*log(pul*(n+5));
   Double_t tslew=(rawtslew<0)?(0):((rawtslew>10)?(10):(rawtslew));


   Double_t pulse[10];
   for (UInt_t i=0; i<10; i++) {pulse[i]=0;}

   Double_t sum = 0.0;
   for (UInt_t i=0; i<10; i++) {
     for (UInt_t j=0; j<25; j++) {
       pulse[i]+=thePulses_.hpdShape_(25*i+j-92.5-tslew);
     }
     sum += pulse[i];
   }

   for (UInt_t i=0; i<10; i++) {
       pulse[i]/=sum;
   }

   for (Int_t i = 3; i<10; i++){
        Double_t pulv = ro->GetBinContent(i+1);
       // differ[i-3] = (pulse[i] - pulv)/pulse[i];
        if (pulv == 0){differ[i-3] = 0;}
        else {differ[i-3] = (pulse[i] - pulv)/pulv;}
        tsc_hb[i-3]->Fill(n,differ[i-3]);
   }

 }
/*
  for(Int_t m=0; m<58; m++){
   char pname[50];
   sprintf(pname, "he_prof_%i",m);

   TProfile *ra = (TProfile*)f1->Get(pname);
  
   Int_t n = m*10+20;
 
   Double_t differ[7]={0.0};

   Double_t pul = ra->GetBinContent(5);


   Double_t rawtslew=13.307784-1.556668*log(pul*(n+5));
   Double_t tslew=(rawtslew<0)?(0):((rawtslew>10)?(10):(rawtslew));


   Double_t pulses[10];
   for (UInt_t i=0; i<10; i++) {pulses[i]=0;}

   Double_t sumhe = 0.0;
   for (UInt_t i=0; i<10; i++) {
     for (UInt_t j=0; j<25; j++) {
       pulses[i]+=thePulses_.hpdShape_(25*i+j-92.5-tslew);
     }
     sumhe += pulses[i];
   }

   for (UInt_t i=0; i<10; i++) {
       pulses[i]/=sumhe;
   }

   for (Int_t i = 3; i<10; i++){
        Double_t pulv = ra->GetBinContent(i+1);
       // differ[i-3] = (pulses[i] - pulv)/pulses[i];
        if (pulv == 0){differ[i-3] = 0;}
        else {differ[i-3] = (pulses[i] - pulv)/pulv;}
        tsc_he[i-3]->Fill(n,differ[i-3]);
   }
 }

 
  for(Int_t m=0; m<58; m++){
   char xname[50];
   sprintf(xname, "hb15_prof_%i",m);


   TProfile *rx = (TProfile*)f1->Get(xname);


   Int_t n = m*10+20;

   Double_t differ[7]={0.0};

   Double_t pul = rx->GetBinContent(5);


   Double_t rawtslew=13.307784-1.556668*log(pul*n);
   Double_t tslew=(rawtslew<0)?(0):((rawtslew>10)?(10):(rawtslew));



   Double_t pulsess[10];
   for (UInt_t i=0; i<10; i++) {pulsess[i]=0;}


   for (UInt_t i=0; i<10; i++) {
     for (UInt_t j=0; j<25; j++) {
       pulsess[i]+=thePulses_.hpdShape_(25*i+j-92.5-tslew);
     }
   }

   for (Int_t i = 3; i<10; i++){
        Double_t pulv = rx->GetBinContent(i+1);
        if (pulv == 0){differ[i-3] = 0;}
        else {differ[i-3] = (pulsess[i] - pulv)/pulv;}
        tsc_15[i-3]->Fill(n,differ[i-3]);
   }
  }
 
   
  for(Int_t m=0; m<58; m++){
   char yname[50];
   sprintf(yname, "hb10_prof_%i",m);


   TProfile *ry = (TProfile*)f1->Get(yname);


   Int_t n = m*10+20;

   Double_t differ[7]={0.0};

   Double_t pul = ry->GetBinContent(5);


   Double_t rawtslew=13.307784-1.556668*log(pul*n);
   Double_t tslew=(rawtslew<0)?(0):((rawtslew>10)?(10):(rawtslew));



   Double_t pulser[10];
   for (UInt_t i=0; i<10; i++) {pulser[i]=0;}


   for (UInt_t i=0; i<10; i++) {
     for (UInt_t j=0; j<25; j++) {
       pulser[i]+=thePulses_.hpdShape_(25*i+j-92.5-tslew);
     }
   }

   for (Int_t i = 3; i<10; i++){
        Double_t pulv = ry->GetBinContent(i+1);
        if (pulv == 0){differ[i-3] = 0;}
        else {differ[i-3] = (pulser[i] - pulv)/pulv;}
        tsc_10[i-3]->Fill(n,differ[i-3]);
   }
  }
 
*/


 f2->Write();
 f2->Close();
 
  
}
