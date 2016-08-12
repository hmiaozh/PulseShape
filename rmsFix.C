#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <iostream>
#include <iomanip>

#endif

double sigma(double ifC) {
  if(ifC < 75) return (0.577 + 0.0686*ifC)/3.;
  return (2.75  + 0.0373*ifC + 3e-6*ifC*ifC)/3.;
}

void rmsFix(TString infile="/afs/cern.ch/user/h/hum/work/public/CMSSW_8_0_7/src/HcalRecoTesting/Results/Output_shift12.root",
	    TString outfile="/afs/cern.ch/user/h/hum/work/public/CMSSW_8_0_1/src/PulseShapes/fitting/Results/digi_shift12.root") {

  const int nbin=58;
  int why=nbin;

  float binsize=10;
  float binstart=20;

  TChain *chain = new TChain("Events");
  chain->AddFile(infile);
  
  //cout << chain->GetEntries() << endl;

  Int_t iphi, ieta, depth;
  Double_t pulse[10], ped[10];

  chain->SetBranchAddress("iphi", &iphi);
  chain->SetBranchAddress("ieta", &ieta);
  chain->SetBranchAddress("depth", &depth);
  chain->SetBranchAddress("pulse", &pulse);
  chain->SetBranchAddress("ped", &ped);

  TFile *outf = new TFile(outfile, "recreate");
  TTree *outt = new TTree("Averages", "Averages");
 
 
  vector<TProfile*> prof;
//  vector<TProfile*> pieta;
  vector<TProfile*> hb_prof;
  vector<TProfile*> he_prof;
//  vector<TProfile*> transarea;
//  vector<TProfile*> endcap;
//  vector<TProfile*> ends;

  /*TProfile* pTS3 = new TProfile("pTS3", "pTS3",nbin,20,20+10*nbin,"s");
  TProfile* pTS4 = new TProfile("pTS4", "pTS4",nbin,20,20+10*nbin,"s");
  TProfile* pTS5 = new TProfile("pTS5", "pTS5",nbin,20,20+10*nbin,"s");
  TProfile* pTS6 = new TProfile("pTS6", "pTS6",nbin,20,20+10*nbin,"s");

  TProfile* pTS3_hb = new TProfile("pTS3_hb", "pTS3_hb",nbin,20,20+10*nbin,"s");
  TProfile* pTS4_hb = new TProfile("pTS4_hb", "pTS4_hb",nbin,20,20+10*nbin,"s");
  TProfile* pTS5_hb = new TProfile("pTS5_hb", "pTS5_hb",nbin,20,20+10*nbin,"s");
  TProfile* pTS6_hb = new TProfile("pTS6_hb", "pTS6_hb",nbin,20,20+10*nbin,"s");

  TProfile* pTS3_he = new TProfile("pTS3_he", "pTS3_he",nbin,20,20+10*nbin,"s");
  TProfile* pTS4_he = new TProfile("pTS4_he", "pTS4_he",nbin,20,20+10*nbin,"s");
  TProfile* pTS5_he = new TProfile("pTS5_he", "pTS5_he",nbin,20,20+10*nbin,"s");
  TProfile* pTS6_he = new TProfile("pTS6_he", "pTS6_he",nbin,20,20+10*nbin,"s");*/

//above originally commented.

/*  
  TProfile* pTS3 = new TProfile("pTS3", "pTS3",80,1,2.8,"s");
  TProfile* pTS4 = new TProfile("pTS4", "pTS4",80,1,2.8,"s");
  TProfile* pTS5 = new TProfile("pTS5", "pTS5",80,1,2.8,"s");
  TProfile* pTS6 = new TProfile("pTS6", "pTS6",80,1,2.8,"s");

  TProfile* pTS3_hb = new TProfile("pTS3_hb", "pTS3_hb",80,1,2.8,"s");
  TProfile* pTS4_hb = new TProfile("pTS4_hb", "pTS4_hb",80,1,2.8,"s");
  TProfile* pTS5_hb = new TProfile("pTS5_hb", "pTS5_hb",80,1,2.8,"s");
  TProfile* pTS6_hb = new TProfile("pTS6_hb", "pTS6_hb",80,1,2.8,"s");

  TProfile* pTS3_he = new TProfile("pTS3_he", "pTS3_he",80,1,2.8,"s");
  TProfile* pTS4_he = new TProfile("pTS4_he", "pTS4_he",80,1,2.8,"s");
  TProfile* pTS5_he = new TProfile("pTS5_he", "pTS5_he",80,1,2.8,"s");
  TProfile* pTS6_he = new TProfile("pTS6_he", "pTS6_he",80,1,2.8,"s");


  TProfile* he = new TProfile("he", "he",10,-0.5,9.5);
  TProfile* hb = new TProfile("hb", "hb",10,-0.5,9.5);
  TProfile* heo = new TProfile("heo", "heo",10,-0.5,9.5);
*/
  


  Int_t nPulses[nbin]={0};  
  Int_t ibin[10]={0,1,2,3,4,5,6,7,8,9};
  Double_t avgPulse[nbin][10], avgPulse2[nbin][10];//, adcUnc[nbin][10], totUnc[nbin][10]; 
  Double_t ucut[nbin], lcut[nbin];

  char pname[50];
  char gname[50];
//  char posname[50];
//  char negname[50];
  for (UInt_t i=0; i<nbin; i++) {
    sprintf(pname,"prof_%i",i);
                       //for (UInt_t j=0; j<10; j++) { 
                       //avgPulse[i][j]=0; avgPulse2[i][j]=0; } //adcUnc[i][j]=0; totUnc[i][j]=0;}
    lcut[i]=binstart+i*binsize;  ucut[i]=binstart+(i+1)*binsize;
    prof.push_back(new TProfile(pname, "",10,-0.5,9.5));
  }

//  for (Int_t i=-29; i<30; i++) {
//    sprintf(gname,"pieta_%i",i);
//    pieta.push_back(new TProfile(gname, "",10,-0.5,9.5));
//  }

  for (UInt_t i=0; i<nbin; i++) {
    sprintf(pname,"hb_prof_%i",i);
    hb_prof.push_back(new TProfile(pname, "",10,-0.5,9.5));
  }

  for (UInt_t i=0; i<nbin; i++) {
    sprintf(pname,"he_prof_%i",i);
    he_prof.push_back(new TProfile(pname, "",10,-0.5,9.5));
  }

/*  for (Int_t j=1; j<3; j++) {
    for (Int_t i=-26; i<-17; i++) {
      sprintf(negname,"endcap_ieta%i_d%i",i,j);
      endcap.push_back(new TProfile(negname, "",10,-0.5,9.5));
    }
    for (Int_t i=18; i<27; i++) {
      sprintf(posname,"endcap_ieta%i_d%i",i,j);
      endcap.push_back(new TProfile(posname, "",10,-0.5,9.5));
    }
  
  }

    for (Int_t j=1; j<4; j++) {
    for (Int_t i=-16; i<-14; i++) {
      sprintf(negname,"transarea_ieta%i_d%i",i,j);
      transarea.push_back(new TProfile(negname, "",10,-0.5,9.5));
    }
    for (Int_t i=15; i<17; i++) {
      sprintf(posname,"transarea_ieta%i_d%i",i,j);
      transarea.push_back(new TProfile(posname, "",10,-0.5,9.5));
    }
  }

  for (Int_t j=1; j<4; j++) {
    for (Int_t i=-29; i<-26; i++) {
      sprintf(negname,"ends_ieta%i_d%i",i,j);
      ends.push_back(new TProfile(negname, "",10,-0.5,9.5));
    }
    for (Int_t i=27; i<30; i++) {
      sprintf(posname,"ends_ieta%i_d%i",i,j);
      ends.push_back(new TProfile(posname, "",10,-0.5,9.5));
    }
  }
*/
  //cout << lcut[0] << ", " << ucut[nbin-1] << endl;

  TH1D *hQ = new TH1D("hq", "hq",nbin, binstart, binstart+nbin*binsize);


  //outt->Branch("why",       &why,       "why/I");
  outt->Branch("nPulses",   &nPulses,   "nPulses[58]/i");
  outt->Branch("ibin",      &ibin,      "ibin[10]/i");
  outt->Branch("ucut",      &ucut,      "ucut[58]/D");
  outt->Branch("lcut",      &lcut,      "lcut[58]/D");
  outt->Branch("avgPulse",  &avgPulse,  "avgPulse[58][10]/D");
  outt->Branch("ieta",      &ieta,      "ieta/i");
  //outt->Branch("avgPulse2", &avgPulse2, "avgPulse2[58][10]/D");
  //outt->Branch("adcUnc",    &adcUnc,    "adcUnc[58][10]/D");

  for (UInt_t i=0; i<chain->GetEntries(); i++) {
  //for (UInt_t i=0; i<1000; i++) {
    chain->GetEntry(i);

    Double_t sumQ=0, ped=0;
    for (UInt_t j=0; j<3; j++) { ped+=pulse[j];}
    for (UInt_t j=0; j<10; j++) { sumQ+=pulse[j];}

//    Double_t ratio[10];
//    for (UInt_t j=0; j<10; j++) {ratio[j]=sumQ/pulse[j]; cout << "the" << j << " ratio is " << ratio[j] << endl;}

    ped/=3;
    sumQ-=10*ped;

    hQ->Fill(sumQ);

    if (sumQ<lcut[0] || sumQ>ucut[nbin-1]) continue;
    //if (TMath::Log10(sumQ)<1 || sumQ>600) continue;

    Int_t ibin= floor( (sumQ-binstart)/binsize );
    
    //if (ibin>20) cout << i << ", " << sumQ << ", " << ibin << endl;
    nPulses[ibin]++;

    //cout << ieta << endl;
/*
    pTS3->Fill(TMath::Log10(sumQ),(pulse[3]-ped)/sumQ);
    pTS4->Fill(TMath::Log10(sumQ),(pulse[4]-ped)/sumQ);
    pTS5->Fill(TMath::Log10(sumQ),(pulse[5]-ped)/sumQ);
    pTS6->Fill(TMath::Log10(sumQ),(pulse[6]-ped)/sumQ);

    if (abs(ieta)<17) {
      pTS3_hb->Fill(TMath::Log10(sumQ),(pulse[3]-ped)/sumQ);
      pTS4_hb->Fill(TMath::Log10(sumQ),(pulse[4]-ped)/sumQ);
      pTS5_hb->Fill(TMath::Log10(sumQ),(pulse[5]-ped)/sumQ);
      pTS6_hb->Fill(TMath::Log10(sumQ),(pulse[6]-ped)/sumQ);
    }
    else {
      pTS3_he->Fill(TMath::Log10(sumQ),(pulse[3]-ped)/sumQ);
      pTS4_he->Fill(TMath::Log10(sumQ),(pulse[4]-ped)/sumQ);
      pTS5_he->Fill(TMath::Log10(sumQ),(pulse[5]-ped)/sumQ);
      pTS6_he->Fill(TMath::Log10(sumQ),(pulse[6]-ped)/sumQ);
    }
*/
      
    for (UInt_t j=0; j<10; j++) { 
      //cout << (pulse[j]-ped)*(pulse[j]-ped)/sumQ << ", ";
      //avgPulse[ibin][j]+=(pulse[j]-ped)/sumQ;
      //avgPulse2[ibin][j]+=(pulse[j]-ped)*(pulse[j]-ped)/sumQ;
      //adcUnc[ibin][j]+=(sigma(pulse[j])*sigma(pulse[j])+1)/(pulse[j]*pulse[j]);
      prof[ibin]->Fill(j,(pulse[j]-ped)/sumQ);
      
    }


    if (abs(ieta)<15){
      for (UInt_t j=0; j<10; j++) {
	hb_prof[ibin]->Fill(j,(pulse[j]-ped)/sumQ);
      }
    }

    if (abs(ieta)>17 && abs(ieta)<27){
      for (UInt_t j=0; j<10; j++) {
        he_prof[ibin]->Fill(j,(pulse[j]-ped)/sumQ);
      }
    }

//    if ( abs(ieta) < 30){
//       Int_t m = ieta+29;
//       for (UInt_t j=0; j<10; j++) {
//          pieta[m]->Fill(j,(pulse[j]-ped)/sumQ);
//       }
//    }

/*    if ( abs(ieta) > 14 && abs(ieta)<17) {
       Int_t d = (depth-1)*4;
       for (UInt_t j=0; j<10; j++) {
          if (ieta < 0){Int_t m = ieta + 16;  transarea[d+m]->Fill(j,(pulse[j]-ped)/sumQ);} 
          else if (ieta >0){ Int_t m = ieta - 13; transarea[d+m]->Fill(j,(pulse[j]-ped)/sumQ);}
       }
    }

    if ( abs(ieta) > 17 && abs(ieta)<27) {
       Int_t d = (depth-1)*18;
       for (UInt_t j=0; j<10; j++) {
          if (ieta < 0){ Int_t m = ieta + 26;  endcap[d+m]->Fill(j,(pulse[j]-ped)/sumQ); }
          else if (ieta > 0){ Int_t m = ieta - 9; endcap[d+m]->Fill(j,(pulse[j]-ped)/sumQ);}
       }
    }

    if ( abs(ieta) > 26 && abs(ieta)<30) {
       Int_t d = (depth-1)*6;
       for (UInt_t j=0; j<10; j++) {
          if (ieta < 0){ Int_t m = ieta + 29;  ends[d+m]->Fill(j,(pulse[j]-ped)/sumQ) ;}
          else if (ieta > 0){ Int_t m = ieta - 24; ends[d+m]->Fill(j,(pulse[j]-ped)/sumQ);}
       }
    }





    
    if (abs(ieta)<15) {
       for (UInt_t j=0; j<10; j++) {
          hb->Fill(j,(pulse[j]-ped)/sumQ);
       }
    }
    else if (abs(ieta)>18 && abs(ieta)<27) {
       for (UInt_t j=0; j<10; j++) {
          he->Fill(j,(pulse[j]-ped)/sumQ);
       }
    }
    else if (abs(ieta)>26 && abs(ieta)<30) {
       for (UInt_t j=0; j<10; j++) {
          heo->Fill(j,(pulse[j]-ped)/sumQ);
       }
    }
*/

    
  }

  //for (UInt_t i=0; i<2; i++) {
  //  cout << endl;
  //  cout << "i : " << i << ", " << nPulses[i] << endl;
  //
  //  for (UInt_t j=0; j<10; j++) {
  //    cout << j << ", " << avgPulse[i][j]/nPulses[i] << ", " << sqrt(avgPulse[i][j]*avgPulse[i][j] - avgPulse2[i][j])/nPulses[i] << endl;
  //    //cout << j << ", " << avgPulse[i][j] << ", " << avgPulse2[i][j] << endl;
  //  }
  //
  //}
  outt->Fill();
  outf->Write();
  cout << "mission completed!" << endl;
  outf->Close();

}
