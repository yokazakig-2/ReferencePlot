#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TText.h>
#include <TArc.h>
#include <TArrow.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <iostream>
#include "stdio.h"
#include <fstream>
#include <string>

#include <math.h>

#include "TVirtualFFT.h"

// Header file for the classes stored in the TTree if any.
#include "vector"

//Double_t uni( Double_t *x, Double_t *par ){
//  Double_t out = 0;
//  if( x[0] > par[0] && x[0] < par[1] ) out = par[2];
//  return out;
//}
Double_t uni( Double_t *x, Double_t *par ){// for time window
  Double_t out = 0;
  if( fabs(x[0]-par[0]) < par[1]/2 ) out = par[2];
  return out;
}
Double_t c2l( Double_t *x, Double_t *par ){// for TDC count -> Drift length
  Double_t out = (x[0] - par[0])/par[1];
  return out;
}
Double_t c2p( Double_t *x, Double_t *par ){// for TDC count -> drift-length
  // 40 - (c-l)/(u - l)*80 = (t_max - count )/(t_max - t_min)*80
  Double_t out = 40. - (par[0] - x[0])/(par[0] - par[1]) * 80.;
  return out;
}
Double_t sq( Double_t *x, Double_t *par ){// for beam waist
  Double_t out = par[0]*(x[0] - par[1])*(x[0] - par[1]) + par[2];
  return out;
}
Double_t line( Double_t *x, Double_t *par ){
  Double_t out = par[1]*x[0] + par[0];
  return out;
}
Double_t line2( Double_t *x, Double_t *par ){
  Double_t out;
  if( x[0] < par[3] )
    out = par[1]*x[0] + par[0];
  else
    out = par[2]*x[0] + par[1]*par[3] + par[0];
  return out;
}
Double_t gausc( Double_t *x, Double_t *par ){
  Double_t xx = x[0];

  Double_t h = par[0];
  Double_t tx = par[1];
  Double_t sigma = par[2];
  Double_t c = par[3];

  Double_t y = h*TMath::Gaus(xx, tx, sigma) + c;

  return y;
}
Double_t gausex( Double_t *x, Double_t *par ){
  Double_t xx = x[0];

  Double_t h = par[0];
  Double_t tx = par[1];
  Double_t sigma = par[2];
  Double_t hi = par[3];
  Double_t slp = par[4];

  Double_t y = h*TMath::Gaus(xx, tx, sigma) + hi *TMath::Exp( slp*xx );

  return y;
}
Double_t gauses( Double_t *x, Double_t *par ){
  Double_t xx = x[0];

  Double_t h = par[0];
  Double_t tx = par[1];
  Double_t sigmal = par[2];
  Double_t sigmar = par[3];

  Double_t y;
  if( xx < tx )
    y = h*TMath::Gaus(xx, tx, sigmal);
  else
    y = h*TMath::Gaus(xx, tx, sigmar);

  return y;
}
Double_t gausesc( Double_t *x, Double_t *par ){
  Double_t xx = x[0];

  Double_t h = par[0];
  Double_t tx = par[1];
  Double_t sigmal = par[2];
  Double_t sigmar = par[3];
  Double_t c = par[4];

  Double_t y;
  if( xx < tx )
    y = h*TMath::Gaus(xx, tx, sigmal);
  else
    y = h*TMath::Gaus(xx, tx, sigmar);

  y = y + c;
  
  return y;
}
Double_t gausescc( Double_t *x, Double_t *par ){
  Double_t xx = x[0];

  Double_t h = par[0];
  Double_t tx = par[1];
  Double_t sigmal = par[2];
  Double_t sigmar = par[3];
  Double_t c = par[4];
  Double_t ll = par[5];
  Double_t ul = par[6];

  Double_t y;
  if( xx > ll && xx < ul ){
    if( xx < tx )
      y = h*TMath::Gaus(xx, tx, sigmal);
    else
      y = h*TMath::Gaus(xx, tx, sigmar);
  }else y = 0;
  
  y = y + c;
  
  return y;
}
Double_t gausesex( Double_t *x, Double_t *par ){
  Double_t xx = x[0];

  Double_t h = par[0];
  Double_t tx = par[1];
  Double_t sigmal = par[2];
  Double_t sigmar = par[3];
  Double_t hi = par[4];
  Double_t slp = par[5];

  Double_t y;
  if( xx < tx )
    y = h*TMath::Gaus(xx, tx, sigmal);
  else
    y = h*TMath::Gaus(xx, tx, sigmar);

  y = y + hi *TMath::Exp( slp*xx );

  return y;
}
Double_t bgexp( Double_t *x, Double_t *par ){
  Double_t xx = x[0];

  Double_t hi = par[0];
  Double_t slp = par[1];

  Double_t y;
  y = hi *TMath::Exp( slp*xx );

  return y;
}
Double_t sggaus( Double_t *x, Double_t *par ){
  Double_t xx = x[0];

  Double_t h = par[0];
  Double_t tx = par[1];
  Double_t sigma = par[2];

  Double_t y = h*TMath::Gaus(xx, tx, sigma);

  return y;
}
Double_t gaus_2D( Double_t *x, Double_t *par ){
  Double_t xx = x[0];
  Double_t yy = x[1];

  Double_t h = fabs(par[0]);// height
  Double_t xc = par[1];// x center
  Double_t yc = par[2];// y center
  Double_t sx = par[3];// sigma x
  Double_t sy = par[4];// sigma y
  Double_t r = par[5];// correlation factor
  Double_t c = par[6];// baseline

  Double_t y = ( h/sqrt(1-r*r) )*exp( (xx-xc)*(xx-xc)/sx/sx + (yy-yc)*(yy-yc)/sy/sy +2*r*(xx-xc)*(yy-yc)/sx/sy/2/(1-r*r ) ) + c;

  return y;
}
Double_t ggs( Double_t *x, Double_t *par ){
  Double_t xx = x[0];

  Double_t h0 = par[0];
  Double_t tx0 = par[1];
  Double_t sigma0 = par[2];
  Double_t h1 = par[3];
  Double_t tx1 = par[4];
  Double_t sigma1l = par[5];
  Double_t sigma1r = par[6];

  Double_t y0 = h0*TMath::Gaus(xx, tx0, sigma0 );
  Double_t y1;
  if( xx < tx1 )
    y1 = h1*TMath::Gaus(xx, tx1, sigma1l );
  else
    y1 = h1*TMath::Gaus(xx, tx1, sigma1r );

  return y0+y1;
}
Double_t egg( Double_t *xx, Double_t *par ){
  Double_t x = xx[0];
  
  Double_t h0 = fabs(par[0]);
  Double_t c0 = fabs(par[1]);
  Double_t s0l = fabs(par[2]);
  Double_t s0r = fabs(par[3]);
  Double_t h1 = fabs(par[4]);
  Double_t c1 = fabs(par[5]);
  Double_t s1 = fabs(par[6]);
  Double_t n0 = fabs(par[7]);
  Double_t sl = -fabs(par[8]);
  
  Double_t y0;
  if( x < c0 ) y0 = h0*TMath::Gaus(x, c0, s0l );
  else  y0 = h0*TMath::Gaus(x, c0, s0r );
  Double_t y1 = h1*TMath::Gaus(x, c1, s1 );
  Double_t y2 = n0 *TMath::Exp( sl*x );
  return y0+y1+y2;
}
Double_t lgg( Double_t *xx, Double_t *par ){
  Double_t x = xx[0];
  
  Double_t h0 = fabs(par[0]);
  Double_t c0 = fabs(par[1]);
  Double_t s0l = fabs(par[2]);
  Double_t s0r = fabs(par[3]);
  Double_t h1 = fabs(par[4]);
  Double_t c1 = fabs(par[5]);
  Double_t s1 = fabs(par[6]);
  Double_t n0 = fabs(par[7]);
  Double_t sl = -fabs(par[8]);
  
  Double_t y0;
  if( x < c0 ) y0 = h0*TMath::Gaus(x, c0, s0l );
  else  y0 = h0*TMath::Gaus(x, c0, s0r );
  Double_t y1 = h1*TMath::Gaus(x, c1, s1 );
  Double_t y2 = TMath::Max( n0 + sl*x, 0. );
  return y0+y1+y2;
}

// Main routine
void tracking(TString runnumber="904", Int_t full = 1 ){

  Int_t runno0 = atoi(runnumber );
  Int_t runno1 = runno0+1;
  TString runnum0, runnum1;

  if( runno0 < 900 )
    if( (runno0/2)*2 != runno0 ){
      cout << "Use even number for run numer less than 900" << endl;
      return;
    }

  Int_t N_files = 1;
  if( runno0 == 991 ) N_files = 5;
  if( runno0 == 996 ) N_files = 3;
  
  Int_t LG_on = 1;
  Float_t LG_offset = 2.27292;
  Float_t LG_slope = 47.4055;
  Float_t LG_plot_offset = 2.28;
  Float_t LG_async_cut = 0.2;
  Float_t LG_1bin_energy = 0.12;
  if( runno0 < 900 ) LG_on = 0;

  Int_t ON_stg = 0;
  if( runno0 > 969 && runno0 < 1017 || runno0 > 3000 && runno0 < 3128 || runno0 > 3157 && runno0 < 3340 ) ON_stg = 1;
  if( runno0 > 3339 ) ON_stg = 2;
  if( runno0 > 3345 ) ON_stg = 3;
  
  Int_t No_tdc = 0;
  if( runno0 >929 && runno0 < 970 ) No_tdc = 1;

  Int_t Nmin_fit = 500;
  
  TString fname0;
  TString fname1;
  TString fname2;

  TString cname;
  cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/ddt2x3300.txt";
  if( runno0 > 3010 && runno0 < 3158 )
    cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/ddt2x3130.txt";
  if( runno0 > 4091 )
    cname = "/group/itdc/tbl/ARTBL_Reference_Data/2025Nov_AR5p0GeV/macro/cal/ddt2x4090.txt";
  /*cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x1017.txt";
  if( runno0 > 899 && runno0 < 970 )
    cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x900.txt";
  if( runno0 > 969 && runno0 < 1017 )
    cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x970.txt";
  if( runno0 < 800 )
    cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x622.txt";
  cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x3029.txt";
  if( runno0 > 3059 )
    cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x3060.txt";
  if( runno0 > 3066 )
    cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x3067.txt";
  if( runno0 > 3074 )
    cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x3075.txt";
  if( runno0 > 3101 )
    cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x3103.txt";
  if( runno0 > 3123 )
    cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x3124.txt";
  if( runno0 > 3130 )
    cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x3133.txt";
  if( runno0 > 3138 )
    cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x3139.txt";
  if( runno0 > 3159 )
  cname = "/group/itdc/tbl/2024ARTBL018/macro/cal/t2x3160.txt";*/

  
  runnum0 = to_string( runno0 );
  runnum1 = to_string( runno1 );

  TString dname;
  dname = "/group/itdc/tbl/2024ARTBL018/macro/cal/dtdc" + runnum0 + ".txt";

  TString zname;
  zname = "/group/itdc/tbl/2024ARTBL018/text/waist" + runnum0 + ".txt";

  if( runno0 < 900 ){
    fname0 = "/group/itdc/tbl/2024ARTBL018/rootfiles/outputs/run" + runnum0 + "_bin.root";  ;//down -> Jet0
    fname1 = "/group/itdc/tbl/2024ARTBL018/rootfiles/outputs/run" + runnum1 + "_bin.root";  ;//up   -> Jet1
    fname2 = "/group/itdc/tbl/2024ARTBL018/rootfiles/outputs/run" + runnum1 + "_bin.root";  ;//up   -> dummy
  }else{
    if( runno0 < 3300 ){
      //                                                                       ON_stg      0      1
      //                                                                                <970   >969
      //                                                                               >1016  <1017
      fname0 = "/group/itdc/tbl/2024ARTBL018/rootfiles/run" + runnumber + "_192.168.10.10_bin.root";  ;//down or up
      fname1 = "/group/itdc/tbl/2024ARTBL018/rootfiles/run" + runnumber + "_192.168.10.16_bin.root";  ;//  up or down
      fname2 = "/group/itdc/tbl/2024ARTBL018/rootfiles/run" + runnumber + "_192.168.10.4_bin.root";  ;//LG
    }else if( runno0 < 4000){
      fname0 = "/group/itdc/tbl/2025ARTBL005/rootfiles/run" + runnumber + "_192.168.10.10_bin.root";  ;//down or up
      fname1 = "/group/itdc/tbl/2025ARTBL005/rootfiles/run" + runnumber + "_192.168.10.16_bin.root";  ;//  up or down
      fname2 = "/group/itdc/tbl/2025ARTBL005/rootfiles/run" + runnumber + "_192.168.10.4_bin.root";  ;//LG
    }else{
      fname0 = " /group/itdc/tbl/ARTBL_Reference_Data/2025Nov_AR5p0GeV/rootfiles/run" + runnumber + "_192.168.10.10_bin.root";  ;//down or up
      fname1 = " /group/itdc/tbl/ARTBL_Reference_Data/2025Nov_AR5p0GeV/rootfiles/run" + runnumber + "_192.168.10.16_bin.root";  ;//  up or down
      fname2 = " /group/itdc/tbl/ARTBL_Reference_Data/2025Nov_AR5p0GeV/rootfiles/run" + runnumber + "_192.168.10.4_bin.root";  ;//LG
    }
   }

  TString pdf_waveform = "/group/itdc/tbl/2024ARTBL018/pdf/waveform/waveform_run" + runnumber + ".pdf(";
  TString pdf_orbit = "/group/itdc/tbl/2024ARTBL018/pdf/orbit/orbit_run" + runnumber + ".pdf(";
  TString pdf_tracking = "/group/itdc/tbl/2024ARTBL018/pdf/2D/tracking_run" + runnumber + ".pdf(";
  
  TFile *f0 = TFile::Open(fname0);//down
  TFile *f1 = TFile::Open(fname1);//up -> in front of LG(970~)
  TFile *f2 = TFile::Open(fname2);//LG

  TTree* tree0 = (TTree*)f0->Get("tree_ev");
  TTree* tree1 = (TTree*)f1->Get("tree_ev");
  TTree* tree2 = (TTree*)f2->Get("tree_ev");

  vector<int> *f0_tdc13;
  vector<int> *f0_tdc14;
  vector<int> *f0_tdc25;
  vector<int> *f0_tdc26;
  vector<int> *f0_tdc13_lt;
  vector<int> *f0_tdc14_lt;
  vector<int> *f0_tdc25_lt;
  vector<int> *f0_tdc26_lt;

  vector<int> *f1_tdc13;
  vector<int> *f1_tdc14;
  vector<int> *f1_tdc25;
  vector<int> *f1_tdc26;
  vector<int> *f1_tdc13_lt;
  vector<int> *f1_tdc14_lt;
  vector<int> *f1_tdc25_lt;
  vector<int> *f1_tdc26_lt;

  vector<int> *f2_adc30;
  vector<int> *f2_tdc30;
  vector<int> *f2_tdc30_lt;

  TBranch        *f0_b_tdc13;   //!
  TBranch        *f0_b_tdc14;   //!
  TBranch        *f0_b_tdc25;   //!
  TBranch        *f0_b_tdc26;   //!
  TBranch        *f0_b_tdc13_lt;   //!
  TBranch        *f0_b_tdc14_lt;   //!
  TBranch        *f0_b_tdc25_lt;   //!
  TBranch        *f0_b_tdc26_lt;   //!

  TBranch        *f1_b_tdc13;   //!
  TBranch        *f1_b_tdc14;   //!
  TBranch        *f1_b_tdc25;   //!
  TBranch        *f1_b_tdc26;   //!
  TBranch        *f1_b_tdc13_lt;   //!
  TBranch        *f1_b_tdc14_lt;   //!
  TBranch        *f1_b_tdc25_lt;   //!
  TBranch        *f1_b_tdc26_lt;   //!

  TBranch        *f2_b_adc30;   //!
  TBranch        *f2_b_tdc30;   //!
  TBranch        *f2_b_tdc30_lt;   //!

  f0_tdc13 = 0;
  f0_tdc14 = 0;
  f0_tdc25 = 0;
  f0_tdc26 = 0;

  f0_tdc13_lt = 0;
  f0_tdc14_lt = 0;
  f0_tdc25_lt = 0;
  f0_tdc26_lt = 0;

  f1_tdc13 = 0;
  f1_tdc14 = 0;
  f1_tdc25 = 0;
  f1_tdc26 = 0;

  f1_tdc13_lt = 0;
  f1_tdc14_lt = 0;
  f1_tdc25_lt = 0;
  f1_tdc26_lt = 0;

  f2_adc30 = 0;
  f2_tdc30 = 0;
  f2_tdc30_lt = 0;

  tree0->SetBranchAddress( "tdc13", &f0_tdc13, &f0_b_tdc13 );
  tree0->SetBranchAddress( "tdc14", &f0_tdc14, &f0_b_tdc14 );
  tree0->SetBranchAddress( "tdc25", &f0_tdc25, &f0_b_tdc25 );
  tree0->SetBranchAddress( "tdc26", &f0_tdc26, &f0_b_tdc26 );

  tree0->SetBranchAddress( "tdc13_lt", &f0_tdc13_lt, &f0_b_tdc13_lt );
  tree0->SetBranchAddress( "tdc14_lt", &f0_tdc14_lt, &f0_b_tdc14_lt );
  tree0->SetBranchAddress( "tdc25_lt", &f0_tdc25_lt, &f0_b_tdc25_lt );
  tree0->SetBranchAddress( "tdc26_lt", &f0_tdc26_lt, &f0_b_tdc26_lt );

  tree1->SetBranchAddress( "tdc13", &f1_tdc13, &f1_b_tdc13 );
  tree1->SetBranchAddress( "tdc14", &f1_tdc14, &f1_b_tdc14 );
  tree1->SetBranchAddress( "tdc25", &f1_tdc25, &f1_b_tdc25 );
  tree1->SetBranchAddress( "tdc26", &f1_tdc26, &f1_b_tdc26 );

  tree1->SetBranchAddress( "tdc13_lt", &f1_tdc13_lt, &f1_b_tdc13_lt );
  tree1->SetBranchAddress( "tdc14_lt", &f1_tdc14_lt, &f1_b_tdc14_lt );
  tree1->SetBranchAddress( "tdc25_lt", &f1_tdc25_lt, &f1_b_tdc25_lt );
  tree1->SetBranchAddress( "tdc26_lt", &f1_tdc26_lt, &f1_b_tdc26_lt );

  if( LG_on ){
    tree2->SetBranchAddress( "adc30", &f2_adc30, &f2_b_adc30 );
    tree2->SetBranchAddress( "tdc30", &f2_tdc30, &f2_b_tdc30 );
    tree2->SetBranchAddress( "tdc30_lt", &f2_tdc30_lt, &f2_b_tdc30_lt );
  }

  //const Int_t offtx0 = -15;
  //const Int_t offty0 = 82;
  //const Int_t offtx1 = -38;
  //const Int_t offty1 = 46;
  const Int_t offty0 = 19;
  const Int_t offtx0 = -39;
  const Int_t offty1 = 68;
  const Int_t offtx1 = -27;

  const Int_t n_tol = 1;
  const Int_t tolty0 = 57;//1314
  const Int_t toltx0 = 98;//2526
  const Int_t tolty1 = 82;//1314
  const Int_t toltx1 = 83;//2526
  
  const Int_t Nch = 8;
  const Int_t Nexp = 6;
  Int_t length;

  Int_t Nstop[Nch];
  Int_t NstopL[Nch];
  Int_t NstopT[Nch];
  Int_t Ntot[Nch];

  Int_t tdc_count[Nch][100];
  Int_t tdc_Lcount[Nch][100];
  Int_t tdc_Tcount[Nch][100];
  Int_t tot_count[Nch][100];
  //Int_t tdc_tot_count[Nch][100];

  Int_t tdc_LCorr[Nch][100];

  Int_t tdc_LCcount[Nch][100];
  Int_t NstopLC[Nch];
  Int_t tot_Ccount[Nch][100];
      
  Int_t xx[1024], yy[1024];
  for( Int_t i = 0; i < 1024; i++ ) xx[i] = i;

  // center and width of beam at LG surface QS on-off
  Float_t xc_on[5] = { -6.09, -7.621, -8.276, -12.8, 0 };// run 991
  Float_t xs_on[5] = { 23.58, 23.25, 22.83, 20.65, 0 };
  Float_t xc_off[5] = { -6.486, -6.656, -6.696, -6.109, 0 };// run 996
  Float_t xs_off[5] = { 23.56, 23.23, 23.07, 22.92, 0 };

  if( !ON_stg ){
    xc_on[0] = -0.7683;// run 1018
    xs_on[0] = 34.17;// run 1018
    xc_off[0] = 2.032;// run 1019
    xs_off[0] = 30.05;// run 1019
    xc_on[1] = -0.1083;// run 1020
    xs_on[1] = 30.67;// run 1020
    xc_off[1] = 3.227;// run 1021
    xs_off[1] = 31.06;// run 1021
    xc_on[2] = -0.689;// run 1022
    xs_on[2] = 28.62;// run 1022
    xc_off[2] = 3.803;// run 1023
    xs_off[2] = 31.12;// run 1023
    xc_on[3] = -2.025;// run 1024
    xs_on[3] = 22.83;// run 1024
    xc_off[3] = 4.713;// run 1025
    xs_off[3] = 30.82;// run 1025
    xc_on[4] = -1.675;// run 1026
    xs_on[4] = 19.44;// run 1026
    xc_off[4] = 8.018;// run 1027
    xs_off[4] = 30.15;// run 1027
  }
  
  //jet0, (25,26,13,14), jet1 (25,26,13,14)
  Double_t a[8];
  Double_t ea[8];
  Double_t b[8];
  Double_t eb[8];

  Double_t dum;
  /*ifstream cfin( cname );
  if( !cfin ){
    for( Int_t i = 0; i < 4; i++ ){
      //b[i] = 3800;
      //a[i] = 5850;
      b[i] = 3850;
      a[i] = 5650;
    }
    for( Int_t i = 4; i < 8; i++ ){
      //b[i] = 3800;
      //a[i] = 5850;
      b[i] = 3900;
      a[i] = 5770;
    }
  }else{
    for( Int_t i = 0; i < 8; i++ ){
      cfin >> dum;
      b[i] = dum;
      cfin >> dum;
      eb[i] = dum;
    }
    for( Int_t i = 0; i < 8; i++ ){
      cfin >> dum;
      a[i] = dum;
      cfin >> dum;
      ea[i] = dum;
    }
    }*/
  
  /*b[0] = 3840;
  b[1] = 3840;
  b[2] = 3845;
  b[3] = 3845;
  b[4] = 3775;
  if( runno0 > 3130 ) b[4] = 3965;
  if( runno0 > 3160 ) b[4] = 3850; 
  b[5] = 3825;
  b[6] = 3775;
  if( runno0 > 3130 ) b[6] = 3965;
  if( runno0 > 3160 ) b[6] = 3810; 
  b[7] = 3775;
  
  a[0] = 5630;
  a[1] = 5635;
  a[2] = 5655;
  a[3] = 5658;
  a[4] = 5610;
  if( runno0 > 3130 ) a[4] = 5800;
  if( runno0 > 3160 ) a[4] = 5610; 
  a[5] = 5610;
  if( runno0 > 3130 ) a[5] = 5750;
  if( runno0 > 3160 ) a[5] = 5600; 
  a[6] = 5600;
  if( runno0 > 3130 ) a[6] = 5600;
  if( runno0 > 3160 ) a[6] = 5600;
  a[7] = 5650;*/

  ifstream cfin( cname );
  if( cfin ){
    for( Int_t i = 0; i < 8; i++ ){
      cfin >> dum;
      b[i] = dum;
      cfin >> dum;
      eb[i] = dum;
    }
    for( Int_t i = 0; i < 8; i++ ){
      cfin >> dum;
      a[i] = dum;
      cfin >> dum;
      ea[i] = dum;
    }
  }


  Double_t dtmean[4];
  Double_t dtsigma[4];

  ifstream dfin( dname );
  if( !dfin ){
    cout << " No dtdc file" << endl;
    for( Int_t i = 0; i < 4; i++ ){
      dtmean[i] = 0;
      dtsigma[i] = 30;
    }
  }else{
    for( Int_t i = 0; i < 4; i++ ){
      dfin >> dum;
      dtmean[i] = dum;
      dfin >> dum;
      dtsigma[i] = dum;
    }
  }
  dfin.close();

  ofstream dfout( dname, ios::out );
    
  TF1 *twin = new TF1( "twin", "uni", 0, 65535, 3 );
  
  //TF1 *x25j0 = new TF1( "x25j0", "c2l", -100, 100, 2 );
  //TF1 *x26j0 = new TF1( "x26j0", "c2l", -100, 100, 2 );
  //TF1 *y13j0 = new TF1( "y13j0", "c2l", -100, 100, 2 );
  //TF1 *y14j0 = new TF1( "y14j0", "c2l", -100, 100, 2 );
  //TF1 *x25j1 = new TF1( "x25j1", "c2l", -100, 100, 2 );
  //TF1 *x26j1 = new TF1( "x26j1", "c2l", -100, 100, 2 );
  //TF1 *y13j1 = new TF1( "y13j1", "c2l", -100, 100, 2 );
  //TF1 *y14j1 = new TF1( "y14j1", "c2l", -100, 100, 2 );

  TF1 *x25j0 = new TF1( "x25j0", "c2p", -100, 100, 2 );
  TF1 *x26j0 = new TF1( "x26j0", "c2p", -100, 100, 2 );
  TF1 *y13j0 = new TF1( "y13j0", "c2p", -100, 100, 2 );
  TF1 *y14j0 = new TF1( "y14j0", "c2p", -100, 100, 2 );
  TF1 *x25j1 = new TF1( "x25j1", "c2p", -100, 100, 2 );
  TF1 *x26j1 = new TF1( "x26j1", "c2p", -100, 100, 2 );
  TF1 *y13j1 = new TF1( "y13j1", "c2p", -100, 100, 2 );
  TF1 *y14j1 = new TF1( "y14j1", "c2p", -100, 100, 2 );

  TF1 *Line = new TF1( "Line", "line", -10000, 10000, 2 );
  TF1 *Linelg = new TF1( "Linelg", "line", -100, 100, 2 );
  
  TF1 *fitpro = new TF1( "fitpro", "gausescc", -100, 100, 7 );
  //TF1 *fitpro = new TF1( "fitpro", "gausesc", -100, 100, 5 );
  //TF1 *fitLG = new TF1( "fitLG", "gauses", 0, 500, 4 );
  TF1 *fitLGT = new TF1( "fitLGT", "gausesex", 0, 500, 6 );
  //TF1 *fitLG = new TF1( "fitLG", "gausesex", 0, 500, 6 );
  //TF1 *bg = new TF1( "bg", "[0]*TMath::Exp( [1]*x )", 0, 500, 2);
  TF1 *bg[14];
  bg[0]= new TF1( "bg00", "bgexp", 0, 500, 2);
  bg[1]= new TF1( "bg01", "bgexp", 0, 500, 2);
  bg[2]= new TF1( "bg02", "bgexp", 0, 500, 2);
  bg[3]= new TF1( "bg03", "bgexp", 0, 500, 2);
  bg[4]= new TF1( "bg04", "bgexp", 0, 500, 2);
  bg[5]= new TF1( "bg05", "bgexp", 0, 500, 2);
  bg[6]= new TF1( "bg06", "bgexp", 0, 500, 2);
  bg[7]= new TF1( "bg07", "bgexp", 0, 500, 2);
  bg[8]= new TF1( "bg08", "bgexp", 0, 500, 2);
  bg[9]= new TF1( "bg09", "bgexp", 0, 500, 2);
  bg[10]= new TF1( "bg10", "bgexp", 0, 500, 2);
  bg[11]= new TF1( "bg11", "bgexp", 0, 500, 2);
  bg[12]= new TF1( "bg12", "bgexp", 0, 500, 2);
  bg[13]= new TF1( "bg13", "bgexp", 0, 500, 2);
  for( Int_t i = 0; i < 14; i++ ){
    bg[i]->SetLineColor(3);
    bg[i]->SetLineStyle(2);
  }
  TF1 *fitLG = new TF1( "fitLG", "lgg", 0, 1200, 9 );
  TF1 *fitintLG = new TF1( "fitintLG", "lgg", 0, 1200, 9 );
  
  TF2 *fit2D = new TF2( "fit2D", "gaus_2D", -100, 100, -100, 100, 7 );

  TF1 *para = new TF1( "para", "sq", -10000, 10000, 3 );

  TF1 *ctdcj0y = new TF1( "ctdcj0y", "ggs", -1000, 1000, 7 );
  TF1 *ctdcj0x = new TF1( "ctdcj0x", "ggs", -1000, 1000, 7 );
  TF1 *ctdcj1y = new TF1( "ctdcj1y", "ggs", -1000, 1000, 7 );
  TF1 *ctdcj1x = new TF1( "ctdcj1x", "ggs", -1000, 1000, 7 );
  
  y13j0->SetParameters( a[0], b[0] );
  y14j0->SetParameters( a[1], b[1] );
  x25j0->SetParameters( a[2], b[2] );
  x26j0->SetParameters( a[3], b[3] );
  y13j1->SetParameters( a[4], b[4] );
  y14j1->SetParameters( a[5], b[5] );
  x25j1->SetParameters( a[6], b[6] );
  x26j1->SetParameters( a[7], b[7] );

  //Float_t x[2], y[2], z[2] = {500, -4000};
  //Float_t x[2], y[2], z[2] = {180, -4000};// chamber position
  Float_t x[2], y[2], z[2] = {180, -3100};// chamber position
  if( ON_stg )
    z[1] = 4700;
  if( ON_stg == 2 ) z[1] = 1380;
  if( ON_stg == 3 ) z[1] = 750;

  // interpolate position
  Float_t px, py;
  //Float_t pz[7] = {-3475, -2950, -2425, -1900, -1375, -850, -325};
  Float_t pz[7] = {-2960, -2300, -1640, -980, -320, -220, 0};
  if( ON_stg ){
    pz[0] = 760; pz[1] = 1320; pz[2] = 1880; pz[3] = 2440;
    pz[4] = 3000; pz[5] = 3560; pz[6] = 4120;
  }
  // finer
  Float_t fpz[6] = { -2900, -2800, -2700, -2600, -2500, -2400 };// -3100 - -2300
  //{ -3868.75, -3737.5, -3606.25, -3343.75, -3212.5, -3081.25 };// -4000 - -2950
  if( ON_stg ){// *760        *1320     *1880
    ////fpz[0] = 1460; fpz[1] = 1600; fpz[2] = 1740;// += 140
    ////fpz[3] = 2020; fpz[4] = 2160; fpz[5] = 2300;
    //fpz[0] = 900; fpz[1] = 1040; fpz[2] = 1180;// += 140 pz[0] ~ pz[1]
    //fpz[3] = 1460; fpz[4] = 1600; fpz[5] = 1740;// pz[1] ~ pz[2]
    fpz[0] = 900; fpz[1] = 1040; fpz[2] = 1180;// += 140 pz[0] ~ pz[1]
    fpz[3] = 1460; fpz[4] = 1600; fpz[5] = 1740;// pz[1] ~ pz[2]
  }

  // extrapolation
  //LG z 4700 + 190 + 47
  Float_t x_lg, y_lg, z_lg = 4937;

  Float_t expz[6];
  //expz[0] = -400;// T4
  //expz[1] = -850;// QSD down
  //expz[2] = -3450;//QSF up
  //expz[3] = -5570;// T0
  //expz[0] = -235;// T4
  //expz[1] = -470;// QSD down
  expz[0] = -470;// QSD down
  expz[1] = -1490;// QSD up
  expz[2] = -2810;//QSF up
  //expz[3] = -4810;// shutter
  expz[3] = -5570;// T0
  expz[4] = expz[2] - 5361.;// Bend down = -8171 mm
  expz[5] = expz[4] - 1780.;// Bend up   = -9951 mm
  
  Float_t expz2[2];
  expz2[0] = -1640;// middle of QS magnets
  expz2[1] = -2980;// between up stream stand and QSF up

  Int_t p10 = 30;// p10 = p(GeV) * 10

  Float_t PHth[14] = { 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 320 };
  Float_t PHcth[3] = { 25, 88, 173 };// 3 GeV/c QS ON
  //Float_t PHcth[3] = { 320, 88, 173 };// 3 GeV/c QS ON temporary
  if( runno0 == 1023 || runno0 == 1004 || runno0 >954 && runno0 < 960 ){// 3 GeV/c QS OFF
    PHcth[1] = 78;
    PHcth[2] = 170;
  }
  if( runno0 == 1024 || runno0 >990 && runno0 < 996 || runno0 >939 && runno0 < 945 ){// 4 GeV/c QS ON
    PHcth[1] = 125;
    PHcth[2] = 225;
  }
  if( runno0 == 1025 || runno0 >995 && runno0 < 999 || runno0 >944 && runno0 < 950 ){// 4 GeV/c QS OFF
    PHcth[1] = 116;
    PHcth[2] = 222;
  }
  if( runno0 == 1020 || runno0 == 999 || runno0 > 959 && runno0 < 965 || runno0 > 969 && runno0 < 976 || runno0 > 908 && runno0 < 914 ){// 2 GeV/c QS ON
    PHcth[1] = 47;
    PHcth[2] = 118;
  }
  if( runno0 == 1021 || runno0 == 1000 || runno0 > 964 && runno0 < 970 || runno0 > 975 && runno0 < 981 ){// 2 GeV/c QS OFF
    PHcth[1] = 39;
    PHcth[2] = 115;
  }
  if( runno0 > 1016 && runno0 < 1019 || runno0 == 1001 || runno0 > 980 && runno0 < 986 || runno0 > 913 && runno0 < 935 ){// 1 GeV/c QS ON
    PHcth[1] = 11;
    PHcth[2] = 60;
  }
  if( runno0 == 1019 || runno0 == 1002 || runno0 > 985 && runno0 < 991 || runno0 > 934 && runno0 < 940 ){// 1 GeV/c QS OF
    PHcth[1] = 10;
    PHcth[2] = 59;
  }
  //  1          2          3          4        4.5
  // 11 (+36=)  47 (+41=)  88 (+37=) 125 (+19=) 144 on
  // 10 (+39=)  39 (+39=)  78 (+38=) 116 (+20=) 136 off
  // 60 (+58=) 118 (+55=) 173 (+52=) 225 (+28=) 253 on
  // 59 (+56=) 115 (+55=) 170 (+52=) 222 (+27=) 249 off
  if( runno0 == 1026 ){// 4.5 GeV/c QS ON
    PHcth[1] = 155;
    PHcth[2] = 249;
  }
  if( runno0 == 1027 ){// 4.5 GeV/c QS OFF
    PHcth[1] = 153;
    PHcth[2] = 249;
  }

  p10 = 30;//
  if( runno0 > 3000 ){ PHcth[1] = 136; PHcth[2] = 153; }
  //if( runno0 == 3139 ){ PHcth[1] = 136; PHcth[2] = 153; }// 3 GeV
  if( runno0 > 3047 && runno0 < 3051 || runno0 > 3062 && runno0 < 3065 || runno0 > 3070 && runno0 < 3073 || runno0 > 3094 && runno0 < 3097 || runno0 == 3119 || runno0 == 3133 || runno0 == 3140 || runno0 == 3149 || runno0 == 3155 || runno0 == 3165 ) p10 = 20;

  if( runno0 > 3050 && runno0 < 3054 || runno0 > 3096 && runno0 < 3099 || runno0 == 3120 || runno0 == 3134 || runno0 == 3141 || runno0 == 3150 || runno0 == 3156 || runno0 == 3166 ) p10 = 10;
  
  if( runno0 > 3055 && runno0 < 3058 || runno0 == 3100 || runno0 == 3122 || runno0 == 3136 || runno0 == 3137 || runno0 == 3142 || runno0 == 3146 || runno0 == 3151 || runno0 == 3160 || runno0 == 3161 ) p10 = 50;
  
  if( runno0 > 3053 && runno0 < 3056 || runno0 > 3064 && runno0 < 3067 || runno0 > 3072 && runno0 < 3075 || runno0 == 3099 || runno0 == 3121 || runno0 == 3135 || runno0 == 3143 || runno0 == 3147 || runno0 == 3152 || runno0 == 3163 ) p10 = 40;

  if( runno0 == 3162 ) p10 = 55;
  if( runno0 == 3167 ) p10 = 5;

  if( p10 == 20 ){ PHcth[1] =  89; PHcth[2] = 103; }// 2 GeV
  if( p10 == 10 ){ PHcth[1] =  42; PHcth[2] =  52; }// 1 GeV
  if( p10 == 50 ){ PHcth[1] = 227; PHcth[2] = 250; }// 5 GeV
  if( p10 == 40 ){ PHcth[1] = 182; PHcth[2] = 204; }// 4 GeV
  if( p10 == 55 ){ PHcth[1] = 251; PHcth[2] = 276; }// 5.5 GeV
  if( p10 == 5 ) { PHcth[1] = 19; PHcth[2] = 28; }// 0.5 GeV
  
  TH1F *hrndm = new TH1F( "random", "random", 1000, -0.5, 1.5 );

  TH1F *ntdc[Nch];
  TH1F *ntdcn[Nch];
  TH1F *ntdcL[Nch];
  TH1F *ntdcT[Nch];

  TH1F *nCorr[Nch];

  TH1F *ntdcLC[Nch];
  
  TH1F *tdc[Nch];
  TH1F *tdcL[Nch];
  TH1F *tdcT[Nch];

  TH1F *tdcLs[Nch];
  TH1F *tdcLm[Nch];

  TH1F *LGlength = new TH1F( "LGlength", "LGlength", 1024, 0, 1024 );
  TH1F *LGintlength = new TH1F( "LGintlength", "LGintlength", 20, 0, 20 );
  TH1F *LGintlength1 = new TH1F( "LGintlength1", "LGintlength1", 20, 0, 20 );
  TH1F *LGintlength2 = new TH1F( "LGintlength2", "LGintlength2", 20, 0, 20 );
  TH1F *LG = new TH1F( "LG", "LG", 80, 0, 320 );
  TH1F *LG_converted = new TH1F( "LG_converted", "LG_converted", 80, 0, 5. );
  TH1F *LG2 = new TH1F( "LG2", "LG2", 80, 0, 320 );
  TH1F *LGff = new TH1F( "LGff", "LGff", 80, 0, 320 );
  TH1F *LGtdc = new TH1F( "LGtdc", "LGtdc", 512, 0, 5120 );
  //TH1F *intLG = new TH1F( "intLG", "intLG", 320, 8000, 11200 );
  TH1F *intLG = new TH1F( "intLG", "intLG", 400, -300, 1300 );
  TH1F *intLGth = new TH1F( "intLGth", "intLGth", 300, 0, 1200 );
  TH1F *intLGpos = new TH1F( "intLGpos", "intLGpos", 240, 0, 960 );
  TH1F *intLGpos1 = new TH1F( "intLGpos1", "intLGpos1", 240, 0, 960 );
  TH1F *intLGpos2 = new TH1F( "intLGpos2", "intLGpos2", 240, 0, 960 );
  TH1F *LGpeak = new TH1F( "LGpeak", "LGpeak", 20, 0, 20 );
  TH1F *LGvalley = new TH1F( "LGvalley", "LGvalley", 20, 0, 20 );
  TH1F *LGpeak2 = new TH1F( "LGpeak2", "LGpeak2", 20, 0, 20 );
    
  TH1F *tot[Nch];
  TH1F *tot1st[Nch];
  TH1F *totsngl[Nch];
  TH1F *totsnglf[Nch];

  //TH1F *totj0p0[Nch];

  TH2F *tot_tdcL[Nch];
  TH2F *tot_tdcL_0[Nch];
  TH2F *tot_tdcL_n[Nch];

  TH2F *tot_ihit[Nch];
  TH2F *tot_ihit_c[Nch];
  TH2F *tot_ihit_nc[Nch];
  TH2F *tot_ihit_lc[Nch];

  TH1F *stat = new TH1F( "stat", "stat", 10, 0, 10 );
  
  TH2F *ncell[2];
  ncell[0] = new TH2F( "ncell0", "ncell0", 3, 0, 3, 3, 0, 3 );
  ncell[1] = new TH2F( "ncell1", "ncell1", 3, 0, 3, 3, 0, 3 );

  ntdc[0] = new TH1F( "ntdc13j0", "ntdc13j0", 20, 0, 20 );
  ntdc[1] = new TH1F( "ntdc14j0", "ntdc14j0", 20, 0, 20 );
  ntdc[2] = new TH1F( "ntdc25j0", "ntdc25j0", 20, 0, 20 );
  ntdc[3] = new TH1F( "ntdc26j0", "ntdc26j0", 20, 0, 20 );
  ntdc[4] = new TH1F( "ntdc13j1", "ntdc13j1", 20, 0, 20 );
  ntdc[5] = new TH1F( "ntdc14j1", "ntdc14j1", 20, 0, 20 );
  ntdc[6] = new TH1F( "ntdc25j1", "ntdc25j1", 20, 0, 20 );
  ntdc[7] = new TH1F( "ntdc26j1", "ntdc26j1", 20, 0, 20 );

  ntdcn[0] = new TH1F( "ntdcn13j0", "ntdcn13j0", 20, 0, 20 );
  ntdcn[1] = new TH1F( "ntdcn14j0", "ntdcn14j0", 20, 0, 20 );
  ntdcn[2] = new TH1F( "ntdcn25j0", "ntdcn25j0", 20, 0, 20 );
  ntdcn[3] = new TH1F( "ntdcn26j0", "ntdcn26j0", 20, 0, 20 );
  ntdcn[4] = new TH1F( "ntdcn13j1", "ntdcn13j1", 20, 0, 20 );
  ntdcn[5] = new TH1F( "ntdcn14j1", "ntdcn14j1", 20, 0, 20 );
  ntdcn[6] = new TH1F( "ntdcn25j1", "ntdcn25j1", 20, 0, 20 );
  ntdcn[7] = new TH1F( "ntdcn26j1", "ntdcn26j1", 20, 0, 20 );

  ntdcL[0] = new TH1F( "ntdcL13j0", "ntdcL13j0", 20, 0, 20 );
  ntdcL[1] = new TH1F( "ntdcL14j0", "ntdcL14j0", 20, 0, 20 );
  ntdcL[2] = new TH1F( "ntdcL25j0", "ntdcL25j0", 20, 0, 20 );
  ntdcL[3] = new TH1F( "ntdcL26j0", "ntdcL26j0", 20, 0, 20 );
  ntdcL[4] = new TH1F( "ntdcL13j1", "ntdcL13j1", 20, 0, 20 );
  ntdcL[5] = new TH1F( "ntdcL14j1", "ntdcL14j1", 20, 0, 20 );
  ntdcL[6] = new TH1F( "ntdcL25j1", "ntdcL25j1", 20, 0, 20 );
  ntdcL[7] = new TH1F( "ntdcL26j1", "ntdcL26j1", 20, 0, 20 );

  ntdcT[0] = new TH1F( "ntdcT13j0", "ntdcT13j0", 20, 0, 20 );
  ntdcT[1] = new TH1F( "ntdcT14j0", "ntdcT14j0", 20, 0, 20 );
  ntdcT[2] = new TH1F( "ntdcT25j0", "ntdcT25j0", 20, 0, 20 );
  ntdcT[3] = new TH1F( "ntdcT26j0", "ntdcT26j0", 20, 0, 20 );
  ntdcT[4] = new TH1F( "ntdcT13j1", "ntdcT13j1", 20, 0, 20 );
  ntdcT[5] = new TH1F( "ntdcT14j1", "ntdcT14j1", 20, 0, 20 );
  ntdcT[6] = new TH1F( "ntdcT25j1", "ntdcT25j1", 20, 0, 20 );
  ntdcT[7] = new TH1F( "ntdcT26j1", "ntdcT26j1", 20, 0, 20 );

  ntdcLC[0] = new TH1F( "ntdcLC13j0", "ntdcLC13j0", 20, 0, 20 );
  ntdcLC[1] = new TH1F( "ntdcLC14j0", "ntdcLC14j0", 20, 0, 20 );
  ntdcLC[2] = new TH1F( "ntdcLC25j0", "ntdcLC25j0", 20, 0, 20 );
  ntdcLC[3] = new TH1F( "ntdcLC26j0", "ntdcLC26j0", 20, 0, 20 );
  ntdcLC[4] = new TH1F( "ntdcLC13j1", "ntdcLC13j1", 20, 0, 20 );
  ntdcLC[5] = new TH1F( "ntdcLC14j1", "ntdcLC14j1", 20, 0, 20 );
  ntdcLC[6] = new TH1F( "ntdcLC25j1", "ntdcLC25j1", 20, 0, 20 );
  ntdcLC[7] = new TH1F( "ntdcLC26j1", "ntdcLC26j1", 20, 0, 20 );

  nCorr[0] = new TH1F( "nCorr13j0", "nCorr13j0", 20, 0, 20 );
  nCorr[1] = new TH1F( "nCorr14j0", "nCorr14j0", 20, 0, 20 );
  nCorr[2] = new TH1F( "nCorr25j0", "nCorr25j0", 20, 0, 20 );
  nCorr[3] = new TH1F( "nCorr26j0", "nCorr26j0", 20, 0, 20 );
  nCorr[4] = new TH1F( "nCorr13j1", "nCorr13j1", 20, 0, 20 );
  nCorr[5] = new TH1F( "nCorr14j1", "nCorr14j1", 20, 0, 20 );
  nCorr[6] = new TH1F( "nCorr25j1", "nCorr25j1", 20, 0, 20 );
  nCorr[7] = new TH1F( "nCorr26j1", "nCorr26j1", 20, 0, 20 );

  tdc[0] = new TH1F( "tdc13j0", "tdc13j0", 1500, 0, 15000 );
  tdc[1] = new TH1F( "tdc14j0", "tdc14j0", 1500, 0, 15000 );
  tdc[2] = new TH1F( "tdc25j0", "tdc25j0", 1500, 0, 15000 );
  tdc[3] = new TH1F( "tdc26j0", "tdc26j0", 1500, 0, 15000 );
  tdc[4] = new TH1F( "tdc13j1", "tdc13j1", 1500, 0, 15000 );
  tdc[5] = new TH1F( "tdc14j1", "tdc14j1", 1500, 0, 15000 );
  tdc[6] = new TH1F( "tdc25j1", "tdc25j1", 1500, 0, 15000 );
  tdc[7] = new TH1F( "tdc26j1", "tdc26j1", 1500, 0, 15000 );

  tdcL[0] = new TH1F( "tdcL13j0", "tdcL13j0", 1500, 0, 15000 );
  tdcL[1] = new TH1F( "tdcL14j0", "tdcL14j0", 1500, 0, 15000 );
  tdcL[2] = new TH1F( "tdcL25j0", "tdcL25j0", 1500, 0, 15000 );
  tdcL[3] = new TH1F( "tdcL26j0", "tdcL26j0", 1500, 0, 15000 );
  tdcL[4] = new TH1F( "tdcL13j1", "tdcL13j1", 1500, 0, 15000 );
  tdcL[5] = new TH1F( "tdcL14j1", "tdcL14j1", 1500, 0, 15000 );
  tdcL[6] = new TH1F( "tdcL25j1", "tdcL25j1", 1500, 0, 15000 );
  tdcL[7] = new TH1F( "tdcL26j1", "tdcL26j1", 1500, 0, 15000 );

  tdcT[0] = new TH1F( "tdcT13j0", "tdcT13j0", 1500, 0, 15000 );
  tdcT[1] = new TH1F( "tdcT14j0", "tdcT14j0", 1500, 0, 15000 );
  tdcT[2] = new TH1F( "tdcT25j0", "tdcT25j0", 1500, 0, 15000 );
  tdcT[3] = new TH1F( "tdcT26j0", "tdcT26j0", 1500, 0, 15000 );
  tdcT[4] = new TH1F( "tdcT13j1", "tdcT13j1", 1500, 0, 15000 );
  tdcT[5] = new TH1F( "tdcT14j1", "tdcT14j1", 1500, 0, 15000 );
  tdcT[6] = new TH1F( "tdcT25j1", "tdcT25j1", 1500, 0, 15000 );
  tdcT[7] = new TH1F( "tdcT26j1", "tdcT26j1", 1500, 0, 15000 );

  tdcLs[0] = new TH1F( "tdcLs13j0", "tdcLs13j0", 1500, 0, 15000 );
  tdcLs[1] = new TH1F( "tdcLs14j0", "tdcLs14j0", 1500, 0, 15000 );
  tdcLs[2] = new TH1F( "tdcLs25j0", "tdcLs25j0", 1500, 0, 15000 );
  tdcLs[3] = new TH1F( "tdcLs26j0", "tdcLs26j0", 1500, 0, 15000 );
  tdcLs[4] = new TH1F( "tdcLs13j1", "tdcLs13j1", 1500, 0, 15000 );
  tdcLs[5] = new TH1F( "tdcLs14j1", "tdcLs14j1", 1500, 0, 15000 );
  tdcLs[6] = new TH1F( "tdcLs25j1", "tdcLs25j1", 1500, 0, 15060 );
  tdcLs[7] = new TH1F( "tdcLs26j1", "tdcLs26j1", 1500, 0, 15000 );

  tdcLm[0] = new TH1F( "tdcLm13j0", "tdcLm13j0", 3600, 2900, 6500 );
  tdcLm[1] = new TH1F( "tdcLm14j0", "tdcLm14j0", 3600, 2900, 6500 );
  tdcLm[2] = new TH1F( "tdcLm25j0", "tdcLm25j0", 3600, 2900, 6500 );
  tdcLm[3] = new TH1F( "tdcLm26j0", "tdcLm26j0", 3600, 2900, 6500 );
  tdcLm[4] = new TH1F( "tdcLm13j1", "tdcLm13j1", 3600, 2900, 6500 );
  tdcLm[5] = new TH1F( "tdcLm14j1", "tdcLm14j1", 3600, 2900, 6500 );
  tdcLm[6] = new TH1F( "tdcLm25j1", "tdcLm25j1", 3600, 2900, 6500 );
  tdcLm[7] = new TH1F( "tdcLm26j1", "tdcLm26j1", 3600, 2900, 6500 );

  tot[0] = new TH1F( "tot13j0", "tot13j0", 100, -500, 500 );
  tot[1] = new TH1F( "tot14j0", "tot14j0", 100, -500, 500 );
  tot[2] = new TH1F( "tot25j0", "tot25j0", 100, -500, 500 );
  tot[3] = new TH1F( "tot26j0", "tot26j0", 100, -500, 500 );
  tot[4] = new TH1F( "tot13j1", "tot13j1", 100, -500, 500 );
  tot[5] = new TH1F( "tot14j1", "tot14j1", 100, -500, 500 );
  tot[6] = new TH1F( "tot25j1", "tot25j1", 100, -500, 500 );
  tot[7] = new TH1F( "tot26j1", "tot26j1", 100, -500, 500 );

  tot1st[0] = new TH1F( "tot1st13j0", "tot1st13j0", 100, -500, 500 );
  tot1st[1] = new TH1F( "tot1st14j0", "tot1st14j0", 100, -500, 500 );
  tot1st[2] = new TH1F( "tot1st25j0", "tot1st25j0", 100, -500, 500 );
  tot1st[3] = new TH1F( "tot1st26j0", "tot1st26j0", 100, -500, 500 );
  tot1st[4] = new TH1F( "tot1st13j1", "tot1st13j1", 100, -500, 500 );
  tot1st[5] = new TH1F( "tot1st14j1", "tot1st14j1", 100, -500, 500 );
  tot1st[6] = new TH1F( "tot1st25j1", "tot1st25j1", 100, -500, 500 );
  tot1st[7] = new TH1F( "tot1st26j1", "tot1st26j1", 100, -500, 500 );

  totsngl[0] = new TH1F( "totsngl13j0", "totsngl13j0", 160, -100, 1500 );
  totsngl[1] = new TH1F( "totsngl14j0", "totsngl14j0", 160, -100, 1500 );
  totsngl[2] = new TH1F( "totsngl25j0", "totsngl25j0", 160, -100, 1500 );
  totsngl[3] = new TH1F( "totsngl26j0", "totsngl26j0", 160, -100, 1500 );
  totsngl[4] = new TH1F( "totsngl13j1", "totsngl13j1", 160, -100, 1500 );
  totsngl[5] = new TH1F( "totsngl14j1", "totsngl14j1", 160, -100, 1500 );
  totsngl[6] = new TH1F( "totsngl25j1", "totsngl25j1", 160, -100, 1500 );
  totsngl[7] = new TH1F( "totsngl26j1", "totsngl26j1", 160, -100, 1500 );

  totsnglf[0] = new TH1F( "totsnglf13j0", "totsnglf13j0", 30, -10, 20 );
  totsnglf[1] = new TH1F( "totsnglf14j0", "totsnglf14j0", 30, -10, 20 );
  totsnglf[2] = new TH1F( "totsnglf25j0", "totsnglf25j0", 30, -10, 20 );
  totsnglf[3] = new TH1F( "totsnglf26j0", "totsnglf26j0", 30, -10, 20 );
  totsnglf[4] = new TH1F( "totsnglf13j1", "totsnglf13j1", 30, -10, 20 );
  totsnglf[5] = new TH1F( "totsnglf14j1", "totsnglf14j1", 30, -10, 20 );
  totsnglf[6] = new TH1F( "totsnglf25j1", "totsnglf25j1", 30, -10, 20 );
  totsnglf[7] = new TH1F( "totsnglf26j1", "totsnglf26j1", 30, -10, 20 );

  /*totj0p0[0] = new TH1F( "totj0p013j0", "totj0p013j0", 100, -500, 500 );
  totj0p0[1] = new TH1F( "totj0p014j0", "totj0p014j0", 100, -500, 500 );
  totj0p0[2] = new TH1F( "totj0p025j0", "totj0p025j0", 100, -500, 500 );
  totj0p0[3] = new TH1F( "totj0p026j0", "totj0p026j0", 100, -500, 500 );
  totj0p0[4] = new TH1F( "totj0p013j1", "totj0p013j1", 100, -500, 500 );
  totj0p0[5] = new TH1F( "totj0p014j1", "totj0p014j1", 100, -500, 500 );
  totj0p0[6] = new TH1F( "totj0p025j1", "totj0p025j1", 100, -500, 500 );
  totj0p0[7] = new TH1F( "totj0p026j1", "totj0p026j1", 100, -500, 500 );*/

  tot_tdcL[0] = new TH2F( "tot_tdc13j0","tot_tdc13j0", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL[1] = new TH2F( "tot_tdc14j0","tot_tdc14j0", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL[2] = new TH2F( "tot_tdc25j0","tot_tdc25j0", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL[3] = new TH2F( "tot_tdc26j0","tot_tdc26j0", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL[4] = new TH2F( "tot_tdc13j1","tot_tdc13j1", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL[5] = new TH2F( "tot_tdc14j1","tot_tdc14j1", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL[6] = new TH2F( "tot_tdc25j1","tot_tdc25j1", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL[7] = new TH2F( "tot_tdc26j1","tot_tdc26j1", 384, 2580, 6420, 50, 0, 500 );

  tot_tdcL_0[0] = new TH2F( "tot_tdc13j0_0","tot_tdc13j0_0", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_0[1] = new TH2F( "tot_tdc14j0_0","tot_tdc14j0_0", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_0[2] = new TH2F( "tot_tdc25j0_0","tot_tdc25j0_0", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_0[3] = new TH2F( "tot_tdc26j0_0","tot_tdc26j0_0", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_0[4] = new TH2F( "tot_tdc13j1_0","tot_tdc13j1_0", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_0[5] = new TH2F( "tot_tdc14j1_0","tot_tdc14j1_0", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_0[6] = new TH2F( "tot_tdc25j1_0","tot_tdc25j1_0", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_0[7] = new TH2F( "tot_tdc26j1_0","tot_tdc26j1_0", 384, 2580, 6420, 50, 0, 500 );

  tot_tdcL_n[0] = new TH2F( "tot_tdc13j0_n","tot_tdc13j0_n", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_n[1] = new TH2F( "tot_tdc14j0_n","tot_tdc14j0_n", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_n[2] = new TH2F( "tot_tdc25j0_n","tot_tdc25j0_n", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_n[3] = new TH2F( "tot_tdc26j0_n","tot_tdc26j0_n", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_n[4] = new TH2F( "tot_tdc13j1_n","tot_tdc13j1_n", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_n[5] = new TH2F( "tot_tdc14j1_n","tot_tdc14j1_n", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_n[6] = new TH2F( "tot_tdc25j1_n","tot_tdc25j1_n", 384, 2580, 6420, 50, 0, 500 );
  tot_tdcL_n[7] = new TH2F( "tot_tdc26j1_n","tot_tdc26j1_n", 384, 2580, 6420, 50, 0, 500 );

  tot_ihit[0] = new TH2F( "tot_i13j0","tot_i13j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit[1] = new TH2F( "tot_i14j0","tot_i14j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit[2] = new TH2F( "tot_i25j0","tot_i25j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit[3] = new TH2F( "tot_i26j0","tot_i26j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit[4] = new TH2F( "tot_i13j1","tot_i13j1", 10, 0, 10, 50, 0, 500 );
  tot_ihit[5] = new TH2F( "tot_i14j1","tot_i14j1", 10, 0, 10, 50, 0, 500 );
  tot_ihit[6] = new TH2F( "tot_i25j1","tot_i25j1", 10, 0, 10, 50, 0, 500 );
  tot_ihit[7] = new TH2F( "tot_i26j1","tot_i26j1", 10, 0, 10, 50, 0, 500 );

  tot_ihit_c[0] = new TH2F( "tot_i_c13j0","tot_i_c13j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit_c[1] = new TH2F( "tot_i_c14j0","tot_i_c14j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit_c[2] = new TH2F( "tot_i_c25j0","tot_i_c25j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit_c[3] = new TH2F( "tot_i_c26j0","tot_i_c26j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit_c[4] = new TH2F( "tot_i_c13j1","tot_i_c13j1", 10, 0, 10, 50, 0, 500 );
  tot_ihit_c[5] = new TH2F( "tot_i_c14j1","tot_i_c14j1", 10, 0, 10, 50, 0, 500 );
  tot_ihit_c[6] = new TH2F( "tot_i_c25j1","tot_i_c25j1", 10, 0, 10, 50, 0, 500 );
  tot_ihit_c[7] = new TH2F( "tot_i_c26j1","tot_i_c26j1", 10, 0, 10, 50, 0, 500 );

  tot_ihit_nc[0] = new TH2F( "tot_i_nc13j0","tot_i_nc13j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit_nc[1] = new TH2F( "tot_i_nc14j0","tot_i_nc14j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit_nc[2] = new TH2F( "tot_i_nc25j0","tot_i_nc25j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit_nc[3] = new TH2F( "tot_i_nc26j0","tot_i_nc26j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit_nc[4] = new TH2F( "tot_i_nc13j1","tot_i_nc13j1", 10, 0, 10, 50, 0, 500 );
  tot_ihit_nc[5] = new TH2F( "tot_i_nc14j1","tot_i_nc14j1", 10, 0, 10, 50, 0, 500 );
  tot_ihit_nc[6] = new TH2F( "tot_i_nc25j1","tot_i_nc25j1", 10, 0, 10, 50, 0, 500 );
  tot_ihit_nc[7] = new TH2F( "tot_i_nc26j1","tot_i_nc26j1", 10, 0, 10, 50, 0, 500 );

  tot_ihit_lc[0] = new TH2F( "tot_i_lc13j0","tot_i_lc13j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit_lc[1] = new TH2F( "tot_i_lc14j0","tot_i_lc14j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit_lc[2] = new TH2F( "tot_i_lc25j0","tot_i_lc25j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit_lc[3] = new TH2F( "tot_i_lc26j0","tot_i_lc26j0", 10, 0, 10, 50, 0, 500 );
  tot_ihit_lc[4] = new TH2F( "tot_i_lc13j1","tot_i_lc13j1", 10, 0, 10, 50, 0, 500 );
  tot_ihit_lc[5] = new TH2F( "tot_i_lc14j1","tot_i_lc14j1", 10, 0, 10, 50, 0, 500 );
  tot_ihit_lc[6] = new TH2F( "tot_i_lc25j1","tot_i_lc25j1", 10, 0, 10, 50, 0, 500 );
  tot_ihit_lc[7] = new TH2F( "tot_i_lc26j1","tot_i_lc26j1", 10, 0, 10, 50, 0, 500 );

  TH2F *totvsPH[4];
  totvsPH[0] = new TH2F( "tot13j0vsPH", "tot13j0vsPH", 80, 0, 320, 100, 0, 500 );
  totvsPH[1] = new TH2F( "tot25j0vsPH", "tot25j0vsPH", 80, 0, 320, 100, 0, 500 );
  totvsPH[2] = new TH2F( "tot13j1vsPH", "tot13j1vsPH", 80, 0, 320, 100, 0, 500 );
  totvsPH[3] = new TH2F( "tot25j1vsPH", "tot25j1vsPH", 80, 0, 320, 100, 0, 500 );

  TH1F *dtdc1313 = new TH1F( "dtdc1313","dtdc1313", 100, -500, 500 );
  TH1F *dtdc2525 = new TH1F( "dtdc2525","dtdc2525", 100, -500, 500 );
  TH1F *dtdc1313j0 = new TH1F( "dtdc1313j0","dtdc1313j0", 100, -500, 500 );
  TH1F *dtdc2525j0 = new TH1F( "dtdc2525j0","dtdc2525j0", 100, -500, 500 );
  TH1F *dtdc1313j1 = new TH1F( "dtdc1313j1","dtdc1313j1", 100, -500, 500 );
  TH1F *dtdc2525j1 = new TH1F( "dtdc2525j1","dtdc2525j1", 100, -500, 500 );
  TH1F *dtdc1314j0 = new TH1F( "dtdc1314j0","dtdc1314j0", 100, -200, 200 );
  TH1F *dtdc2526j0 = new TH1F( "dtdc2526j0","dtdc2526j0", 100, -200, 200 );
  TH1F *dtdc1314j1 = new TH1F( "dtdc1314j1","dtdc1314j1", 100, -200, 200 );
  TH1F *dtdc2526j1 = new TH1F( "dtdc2526j1","dtdc2526j1", 100, -200, 200 );

  TH2F *dtdc1314_tdc13j0 = new TH2F( "dtdc1314_tdc13j0","dtdc1314_tdc13j0", 384, 2580, 6420, 100, -200, 200 );
  TH2F *dtdc2526_tdc25j0 = new TH2F( "dtdc2526_tdc25j0","dtdc2526_tdc25j0", 384, 2580, 6420, 100, -200, 200 );
  TH2F *dtdc1314_tdc13j1 = new TH2F( "dtdc1314_tdc13j1","dtdc1314_tdc13j1", 384, 2580, 6420, 100, -200, 200 );
  TH2F *dtdc2526_tdc25j1 = new TH2F( "dtdc2526_tdc25j1","dtdc2526_tdc25j1", 384, 2580, 6420, 100, -200, 200 );

  
  TH1F *dtdc13m14mj0 = new TH1F( "dtdc13m14mj0","dtdc13m14mj0", 100, -200, 200 );
  TH1F *dtdc25m26mj0 = new TH1F( "dtdc25m26mj0","dtdc25m26mj0", 100, -200, 200 );
  TH1F *dtdc13m14mj1 = new TH1F( "dtdc13m14mj1","dtdc13m14mj1", 100, -200, 200 );
  TH1F *dtdc25m26mj1 = new TH1F( "dtdc25m26mj1","dtdc25m26mj1", 100, -200, 200 );

  TH2F *dtdc13m14m_tdc13mj0 = new TH2F( "dtdc13m14m_tdc13mj0","dtdc13m14m_tdc13mj0", 384, 2580, 6420, 100, -200, 200 );
  TH2F *dtdc25m26m_tdc25mj0 = new TH2F( "dtdc25m26m_tdc25mj0","dtdc25m26m_tdc25mj0", 384, 2580, 6420, 100, -200, 200 );
  TH2F *dtdc13m14m_tdc13mj1 = new TH2F( "dtdc13m14m_tdc13mj1","dtdc13m14m_tdc13mj1", 384, 2580, 6420, 100, -200, 200 );
  TH2F *dtdc25m26m_tdc25mj1 = new TH2F( "dtdc25m26m_tdc25mj1","dtdc25m26m_tdc25mj1", 384, 2580, 6420, 100, -200, 200 );



  TH2F *tdc25_13j0 = new TH2F( "tdc25_13j0","tdc25_13j0", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25_13j1 = new TH2F( "tdc25_13j1","tdc25_13j1", 384, 2580, 6420, 384, 2580, 6420 );

  TH2F *tdc25c_13cj0 = new TH2F( "tdc25c_13cj0","tdc25c_13cj0", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25c_13cj1 = new TH2F( "tdc25c_13cj1","tdc25c_13cj1", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25nc_13ncj0 = new TH2F( "tdc25nc_13ncj0","tdc25nc_13ncj0", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25nc_13ncj1 = new TH2F( "tdc25nc_13ncj1","tdc25nc_13ncj1", 384, 2580, 6420, 384, 2580, 6420 );

  TH2F *tdc25m_13mj0 = new TH2F( "tdc25m_13mj0","tdc25m_13mj0", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj1 = new TH2F( "tdc25m_13mj1","tdc25m_13mj1", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj0cx = new TH2F( "tdc25m_13mj0cx","tdc25m_13mj0cx", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj1cx = new TH2F( "tdc25m_13mj1cx","tdc25m_13mj1cx", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj0ncx = new TH2F( "tdc25m_13mj0ncx","tdc25m_13mj0ncx", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj1ncx = new TH2F( "tdc25m_13mj1ncx","tdc25m_13mj1ncx", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj0cy = new TH2F( "tdc25m_13mj0cy","tdc25m_13mj0cy", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj1cy = new TH2F( "tdc25m_13mj1cy","tdc25m_13mj1cy", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj0ncy = new TH2F( "tdc25m_13mj0ncy","tdc25m_13mj0ncy", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj1ncy = new TH2F( "tdc25m_13mj1ncy","tdc25m_13mj1ncy", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj0n1 = new TH2F( "tdc25m_13mj0n1","tdc25m_13mj0n1", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj0w1 = new TH2F( "tdc25m_13mj0w1","tdc25m_13mj0w1", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj1n0 = new TH2F( "tdc25m_13mj1n0","tdc25m_13mj1n0", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25m_13mj1w0 = new TH2F( "tdc25m_13mj1w0","tdc25m_13mj1w0", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25mc_13mcj0 = new TH2F( "tdc25mc_13mcj0","tdc25mc_13mcj0", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc25mc_13mcj1 = new TH2F( "tdc25mc_13mcj1","tdc25mc_13mcj1", 384, 2580, 6420, 384, 2580, 6420 );

  TH2F *xy25m_13mj0 = new TH2F( "xy25m_13mj0","xy25m_13mj0", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj1 = new TH2F( "xy25m_13mj1","xy25m_13mj1", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj0 = new TH2F( "xy26m_14mj0","xy26m_14mj0", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj1 = new TH2F( "xy26m_14mj1","xy26m_14mj1", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13mj0w1 = new TH2F( "xy25m_13mj0w1","xy25m_13mj0w1", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj1w0 = new TH2F( "xy25m_13mj1w0","xy25m_13mj1w0", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj0w1 = new TH2F( "xy26m_14mj0w1","xy26m_14mj0w1", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj1w0 = new TH2F( "xy26m_14mj1w0","xy26m_14mj1w0", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13mj0n1 = new TH2F( "xy25m_13mj0n1","xy25m_13mj0n1", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj1n0 = new TH2F( "xy25m_13mj1n0","xy25m_13mj1n0", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj0n1 = new TH2F( "xy26m_14mj0n1","xy26m_14mj0n1", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj1n0 = new TH2F( "xy26m_14mj1n0","xy26m_14mj1n0", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13mj0w1_25_13 = new TH2F( "xy25m_13mj0w1_25_13","xy25m_13mj0w1_25_13", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj1w0_25_13 = new TH2F( "xy25m_13mj1w0_25_13","xy25m_13mj1w0_25_13", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj0w1_26_14 = new TH2F( "xy26m_14mj0w1_26_14","xy26m_14mj0w1_26_14", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj1w0_26_14 = new TH2F( "xy26m_14mj1w0_26_14","xy26m_14mj1w0_26_14", 100, -50, 50, 100, -50, 50 );

  // 2D map outside of time window
  TH2F *xy25o_13oj0 = new TH2F( "xy25o_13oj0","xy25o_13oj0", 400, -100, 300, 400, -100, 300 );
  TH2F *xy25o_13oj1 = new TH2F( "xy25o_13oj1","xy25o_13oj1", 400, -100, 300, 400, -100, 300 );
  TH2F *xy26o_14oj0 = new TH2F( "xy26o_14oj0","xy26o_14oj0", 400, -100, 300, 400, -100, 300 );
  TH2F *xy26o_14oj1 = new TH2F( "xy26o_14oj1","xy26o_14oj1", 400, -100, 300, 400, -100, 300 );

  TH2F *xy25o_13oj0w1 = new TH2F( "xy25o_13oj0w1","xy25o_13oj0w1", 400, -100, 300, 400, -100, 300 );
  TH2F *xy25o_13oj1w0 = new TH2F( "xy25o_13oj1w0","xy25o_13oj1w0", 400, -100, 300, 400, -100, 300 );
  TH2F *xy26o_14oj0w1 = new TH2F( "xy26o_14oj0w1","xy26o_14oj0w1", 400, -100, 300, 400, -100, 300 );
  TH2F *xy26o_14oj1w0 = new TH2F( "xy26o_14oj1w0","xy26o_14oj1w0", 400, -100, 300, 400, -100, 300 );

  TH2F *xy25o_13oj0n1 = new TH2F( "xy25o_13oj0n1","xy25o_13oj0n1", 400, -100, 300, 400, -100, 300 );
  TH2F *xy25o_13oj1n0 = new TH2F( "xy25o_13oj1n0","xy25o_13oj1n0", 400, -100, 300, 400, -100, 300 );
  TH2F *xy26o_14oj0n1 = new TH2F( "xy26o_14oj0n1","xy26o_14oj0n1", 400, -100, 300, 400, -100, 300 );
  TH2F *xy26o_14oj1n0 = new TH2F( "xy26o_14oj1n0","xy26o_14oj1n0", 400, -100, 300, 400, -100, 300 );

  // 2D map poles effect
  TH2F *xy25m_13mj0p1 = new TH2F( "xy25m_13mj0p1","xy25m_13mj0p1", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj1p0 = new TH2F( "xy25m_13mj1p0","xy25m_13mj1p0", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj0p1 = new TH2F( "xy26m_14mj0p1","xy26m_14mj0p1", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj1p0 = new TH2F( "xy26m_14mj1p0","xy26m_14mj1p0", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13mj0c1 = new TH2F( "xy25m_13mj0c1","xy25m_13mj0c1", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj1c0 = new TH2F( "xy25m_13mj1c0","xy25m_13mj1c0", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj0c1 = new TH2F( "xy26m_14mj0c1","xy26m_14mj0c1", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj1c0 = new TH2F( "xy26m_14mj1c0","xy26m_14mj1c0", 100, -50, 50, 100, -50, 50 );

  // 2D map each pole effect
  TH2F *xy25m_13mj0p10 = new TH2F( "xy25m_13mj0p10","xy25m_13mj0p10", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj0p11 = new TH2F( "xy25m_13mj0p11","xy25m_13mj0p11", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj0p12 = new TH2F( "xy25m_13mj0p12","xy25m_13mj0p12", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj0p13 = new TH2F( "xy25m_13mj0p13","xy25m_13mj0p13", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj1p00 = new TH2F( "xy25m_13mj1p00","xy25m_13mj1p00", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj1p01 = new TH2F( "xy25m_13mj1p01","xy25m_13mj1p01", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj1p02 = new TH2F( "xy25m_13mj1p02","xy25m_13mj1p02", 100, -50, 50, 100, -50, 50 );
  TH2F *xy25m_13mj1p03 = new TH2F( "xy25m_13mj1p03","xy25m_13mj1p03", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj0p10 = new TH2F( "xy26m_14mj0p10","xy26m_14mj0p10", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj0p11 = new TH2F( "xy26m_14mj0p11","xy26m_14mj0p11", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj0p12 = new TH2F( "xy26m_14mj0p12","xy26m_14mj0p12", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj0p13 = new TH2F( "xy26m_14mj0p13","xy26m_14mj0p13", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj1p00 = new TH2F( "xy26m_14mj1p00","xy26m_14mj1p00", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj1p01 = new TH2F( "xy26m_14mj1p01","xy26m_14mj1p01", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj1p02 = new TH2F( "xy26m_14mj1p02","xy26m_14mj1p02", 100, -50, 50, 100, -50, 50 );
  TH2F *xy26m_14mj1p03 = new TH2F( "xy26m_14mj1p03","xy26m_14mj1p03", 100, -50, 50, 100, -50, 50 );

  // wire difference in same chambers
  TH2F *xy_dxdy_mj0 = new TH2F( "xy_dxdy_mj0","xy_dxdy_mj0", 200, -10, 10, 200, -10, 10 );
  TH2F *xy_dxdy_mj1 = new TH2F( "xy_dxdy_mj1","xy_dxdy_mj1", 200, -10, 10, 200, -10, 10 );

  TH2F *dxdy25m_13mj1_j0 = new TH2F( "dxdy25m_13mj1_j0","dxdy25m_13mj1_j0", 100, -50, 50, 100, -50, 50 );
  TH2F *xdx25m_13mj1_j0 = new TH2F( "xdx25m_13mj1_j0","xdx25m_13mj1_j0", 100, -50, 50, 100, -50, 50 );
  TH2F *ydy25m_13mj1_j0 = new TH2F( "ydy25m_13mj1_j0","ydy25m_13mj1_j0", 100, -50, 50, 100, -50, 50 );
  TH2F *xdy25m_13mj1_j0 = new TH2F( "xdy25m_13mj1_j0","xdy25m_13mj1_j0", 100, -50, 50, 100, -50, 50 );
  TH2F *ydx25m_13mj1_j0 = new TH2F( "ydx25m_13mj1_j0","ydx25m_13mj1_j0", 100, -50, 50, 100, -50, 50 );

  TH2F *LGvsEv_j0 = new TH2F( "LGvsEv_j0", "LGvsEv_j0", 30, 0, 300000, 16, 0, 320 );
  TH2F *LGvsEv_j1 = new TH2F( "LGvsEv_j1", "LGvsEv_j1", 30, 0, 300000, 16, 0, 320 ); 
    
  // interpolated 2D map
  TH2F *xy25m_13m_z[7];
  xy25m_13m_z[0]= new TH2F( "xy25m_13m_z0","xy25m_13m_z0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_z[1]= new TH2F( "xy25m_13m_z1","xy25m_13m_z1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_z[2]= new TH2F( "xy25m_13m_z2","xy25m_13m_z2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_z[3]= new TH2F( "xy25m_13m_z3","xy25m_13m_z3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_z[4]= new TH2F( "xy25m_13m_z4","xy25m_13m_z4", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_z[5]= new TH2F( "xy25m_13m_z5","xy25m_13m_z5", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_z[6]= new TH2F( "xy25m_13m_z6","xy25m_13m_z6", 100, -50, 50, 100, -50, 50 );

  // finer intepolated 2D map
  TH2F *xy25m_13m_fz[6];
  xy25m_13m_fz[0]= new TH2F( "xy25m_13m_fz0","xy25m_13m_fz0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_fz[1]= new TH2F( "xy25m_13m_fz1","xy25m_13m_fz1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_fz[2]= new TH2F( "xy25m_13m_fz2","xy25m_13m_fz2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_fz[3]= new TH2F( "xy25m_13m_fz3","xy25m_13m_fz3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_fz[4]= new TH2F( "xy25m_13m_fz4","xy25m_13m_fz4", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_fz[5]= new TH2F( "xy25m_13m_fz5","xy25m_13m_fz5", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_undetz[2];
  xy25m_13m_undetz[0]= new TH2F( "xy25m_13m_undetz0","xy25m_13m_undetz0", 400, -100, 300, 400, -100, 300 );
  xy25m_13m_undetz[1]= new TH2F( "xy25m_13m_undetz1","xy25m_13m_undetz1", 400, -100, 300, 400, -100, 300 );

  // extrapolated 2D map
  TH2F *xy25m_13m_detz[2];
  xy25m_13m_detz[0]= new TH2F( "xy25m_13m_detz0","xy25m_13m_detz0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_detz[1]= new TH2F( "xy25m_13m_detz1","xy25m_13m_detz1", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_detz_E0[2];
  xy25m_13m_detz_E0[0]= new TH2F( "xy25m_13m_detz0_e0","xy25m_13m_detz0_e0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_detz_E0[1]= new TH2F( "xy25m_13m_detz1_e0","xy25m_13m_detz1_e0", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_detz_E1[2];
  xy25m_13m_detz_E1[0]= new TH2F( "xy25m_13m_detz0_e1","xy25m_13m_detz0_e1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_detz_E1[1]= new TH2F( "xy25m_13m_detz1_e1","xy25m_13m_detz1_e1", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_detz_E2[2];
  xy25m_13m_detz_E2[0]= new TH2F( "xy25m_13m_detz0_e2","xy25m_13m_detz0_e2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_detz_E2[1]= new TH2F( "xy25m_13m_detz1_e2","xy25m_13m_detz1_e2", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_detz_E3[2];
  xy25m_13m_detz_E3[0]= new TH2F( "xy25m_13m_detz0_e3","xy25m_13m_detz0_e3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_detz_E3[1]= new TH2F( "xy25m_13m_detz1_e3","xy25m_13m_detz1_e3", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_exz_E0[6];
  xy25m_13m_exz_E0[0]= new TH2F( "xy25m_13m_exz0_e0","xy25m_13m_exz0_e0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_E0[1]= new TH2F( "xy25m_13m_exz1_e0","xy25m_13m_exz1_e0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_E0[2]= new TH2F( "xy25m_13m_exz2_e0","xy25m_13m_exz2_e0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_E0[3]= new TH2F( "xy25m_13m_exz3_e0","xy25m_13m_exz3_e0", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_E0[4]= new TH2F( "xy25m_13m_exz3_e0","xy25m_13m_exz4_e0", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_E0[5]= new TH2F( "xy25m_13m_exz3_e0","xy25m_13m_exz5_e0", 100, -100, 100, 100, -100, 100 );

  TH2F *xy25m_13m_exz_E1[6];
  xy25m_13m_exz_E1[0]= new TH2F( "xy25m_13m_exz0_e1","xy25m_13m_exz0_e1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_E1[1]= new TH2F( "xy25m_13m_exz1_e1","xy25m_13m_exz1_e1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_E1[2]= new TH2F( "xy25m_13m_exz2_e1","xy25m_13m_exz2_e1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_E1[3]= new TH2F( "xy25m_13m_exz3_e1","xy25m_13m_exz3_e1", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_E1[4]= new TH2F( "xy25m_13m_exz3_e0","xy25m_13m_exz4_e1", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_E1[5]= new TH2F( "xy25m_13m_exz3_e0","xy25m_13m_exz5_e1", 100, -100, 100, 100, -100, 100 );

  TH2F *xy25m_13m_exz_E2[6];
  xy25m_13m_exz_E2[0]= new TH2F( "xy25m_13m_exz0_e2","xy25m_13m_exz0_e2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_E2[1]= new TH2F( "xy25m_13m_exz1_e2","xy25m_13m_exz1_e2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_E2[2]= new TH2F( "xy25m_13m_exz2_e2","xy25m_13m_exz2_e2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_E2[3]= new TH2F( "xy25m_13m_exz3_e2","xy25m_13m_exz3_e2", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_E2[4]= new TH2F( "xy25m_13m_exz3_e0","xy25m_13m_exz4_e2", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_E2[5]= new TH2F( "xy25m_13m_exz3_e0","xy25m_13m_exz5_e2", 100, -100, 100, 100, -100, 100 );

  TH2F *xy25m_13m_exz_E3[6];
  xy25m_13m_exz_E3[0]= new TH2F( "xy25m_13m_exz0_e3","xy25m_13m_exz0_e3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_E3[1]= new TH2F( "xy25m_13m_exz1_e3","xy25m_13m_exz1_e3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_E3[2]= new TH2F( "xy25m_13m_exz2_e3","xy25m_13m_exz2_e3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_E3[3]= new TH2F( "xy25m_13m_exz3_e3","xy25m_13m_exz3_e3", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_E3[4]= new TH2F( "xy25m_13m_exz3_e0","xy25m_13m_exz4_e3", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_E3[5]= new TH2F( "xy25m_13m_exz3_e0","xy25m_13m_exz5_e3", 100, -100, 100, 100, -100, 100 );

  TH2F *xy25m_13m_detz_LGxneg_E0[2];
  xy25m_13m_detz_LGxneg_E0[0]= new TH2F( "xy25m_13m_detz_LGxneg0_e0","xy25m_13m_detz_LGxneg0_e0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_detz_LGxneg_E0[1]= new TH2F( "xy25m_13m_detz_LGxneg1_e0","xy25m_13m_detz_LGxneg1_e0", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_detz_LGxneg_E1[2];
  xy25m_13m_detz_LGxneg_E1[0]= new TH2F( "xy25m_13m_detz_LGxneg0_e1","xy25m_13m_detz_LGxneg0_e1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_detz_LGxneg_E1[1]= new TH2F( "xy25m_13m_detz_LGxneg1_e1","xy25m_13m_detz_LGxneg1_e1", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_detz_LGxneg_E2[2];
  xy25m_13m_detz_LGxneg_E2[0]= new TH2F( "xy25m_13m_detz_LGxneg0_e2","xy25m_13m_detz_LGxneg0_e2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_detz_LGxneg_E2[1]= new TH2F( "xy25m_13m_detz_LGxneg1_e2","xy25m_13m_detz_LGxneg1_e2", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_detz_LGxneg_E3[2];
  xy25m_13m_detz_LGxneg_E3[0]= new TH2F( "xy25m_13m_detz_LGxneg0_e3","xy25m_13m_detz_LGxneg0_e3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_detz_LGxneg_E3[1]= new TH2F( "xy25m_13m_detz_LGxneg1_e3","xy25m_13m_detz_LGxneg1_e3", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_exz_LGxneg_E0[6];
  xy25m_13m_exz_LGxneg_E0[0]= new TH2F( "xy25m_13m_exz_LGxneg0_e0","xy25m_13m_exz_LGxneg0_e0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_LGxneg_E0[1]= new TH2F( "xy25m_13m_exz_LGxneg1_e0","xy25m_13m_exz_LGxneg1_e0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_LGxneg_E0[2]= new TH2F( "xy25m_13m_exz_LGxneg2_e0","xy25m_13m_exz_LGxneg2_e0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_LGxneg_E0[3]= new TH2F( "xy25m_13m_exz_LGxneg3_e0","xy25m_13m_exz_LGxneg3_e0", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_LGxneg_E0[4]= new TH2F( "xy25m_13m_exz_LGxneg3_e0","xy25m_13m_exz_LGxneg4_e0", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_LGxneg_E0[5]= new TH2F( "xy25m_13m_exz_LGxneg3_e0","xy25m_13m_exz_LGxneg5_e0", 100, -100, 100, 100, -100, 100 );

  TH2F *xy25m_13m_exz_LGxneg_E1[6];
  xy25m_13m_exz_LGxneg_E1[0]= new TH2F( "xy25m_13m_exz_LGxneg0_e1","xy25m_13m_exz_LGxneg0_e1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_LGxneg_E1[1]= new TH2F( "xy25m_13m_exz_LGxneg1_e1","xy25m_13m_exz_LGxneg1_e1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_LGxneg_E1[2]= new TH2F( "xy25m_13m_exz_LGxneg2_e1","xy25m_13m_exz_LGxneg2_e1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_LGxneg_E1[3]= new TH2F( "xy25m_13m_exz_LGxneg3_e1","xy25m_13m_exz_LGxneg3_e1", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_LGxneg_E1[4]= new TH2F( "xy25m_13m_exz_LGxneg3_e0","xy25m_13m_exz_LGxneg4_e1", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_LGxneg_E1[5]= new TH2F( "xy25m_13m_exz_LGxneg3_e0","xy25m_13m_exz_LGxneg5_e1", 100, -100, 100, 100, -100, 100 );

  TH2F *xy25m_13m_exz_LGxneg_E2[6];
  xy25m_13m_exz_LGxneg_E2[0]= new TH2F( "xy25m_13m_exz_LGxneg0_e2","xy25m_13m_exz_LGxneg0_e2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_LGxneg_E2[1]= new TH2F( "xy25m_13m_exz_LGxneg1_e2","xy25m_13m_exz_LGxneg1_e2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_LGxneg_E2[2]= new TH2F( "xy25m_13m_exz_LGxneg2_e2","xy25m_13m_exz_LGxneg2_e2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_LGxneg_E2[3]= new TH2F( "xy25m_13m_exz_LGxneg3_e2","xy25m_13m_exz_LGxneg3_e2", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_LGxneg_E2[4]= new TH2F( "xy25m_13m_exz_LGxneg3_e0","xy25m_13m_exz_LGxneg4_e2", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_LGxneg_E2[5]= new TH2F( "xy25m_13m_exz_LGxneg3_e0","xy25m_13m_exz_LGxneg5_e2", 100, -100, 100, 100, -100, 100 );

  TH2F *xy25m_13m_exz_LGxneg_E3[6];
  xy25m_13m_exz_LGxneg_E3[0]= new TH2F( "xy25m_13m_exz_LGxneg0_e3","xy25m_13m_exz_LGxneg0_e3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_LGxneg_E3[1]= new TH2F( "xy25m_13m_exz_LGxneg1_e3","xy25m_13m_exz_LGxneg1_e3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_LGxneg_E3[2]= new TH2F( "xy25m_13m_exz_LGxneg2_e3","xy25m_13m_exz_LGxneg2_e3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz_LGxneg_E3[3]= new TH2F( "xy25m_13m_exz_LGxneg3_e3","xy25m_13m_exz_LGxneg3_e3", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_LGxneg_E3[4]= new TH2F( "xy25m_13m_exz_LGxneg3_e0","xy25m_13m_exz_LGxneg4_e3", 100, -100, 100, 100, -100, 100 );
  xy25m_13m_exz_LGxneg_E3[5]= new TH2F( "xy25m_13m_exz_LGxneg3_e0","xy25m_13m_exz_LGxneg5_e3", 100, -100, 100, 100, -100, 100 );
 
  TH2F *xy25m_13m_exz2_E0[4];
  xy25m_13m_exz2_E0[0]= new TH2F( "xy25m_13m_exz20_e0","xy25m_13m_exz20_e0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz2_E0[1]= new TH2F( "xy25m_13m_exz21_e0","xy25m_13m_exz21_e0", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_exz2_E1[4];
  xy25m_13m_exz2_E1[0]= new TH2F( "xy25m_13m_exz20_e1","xy25m_13m_exz20_e1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz2_E1[1]= new TH2F( "xy25m_13m_exz21_e1","xy25m_13m_exz21_e1", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_exz2_E2[4];
  xy25m_13m_exz2_E2[0]= new TH2F( "xy25m_13m_exz20_e2","xy25m_13m_exz20_e2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz2_E2[1]= new TH2F( "xy25m_13m_exz21_e2","xy25m_13m_exz21_e2", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_exz2_E3[4];
  xy25m_13m_exz2_E3[0]= new TH2F( "xy25m_13m_exz20_e3","xy25m_13m_exz20_e3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_exz2_E3[1]= new TH2F( "xy25m_13m_exz21_e3","xy25m_13m_exz21_e3", 100, -50, 50, 100, -50, 50 );

  
  /*  TH2F *xy25m_13m_j0p0_exz_E0[4];
  xy25m_13m_j0p0_exz_E0[0]= new TH2F( "xy25m_13m_j0p0_exz0_e0","xy25m_13m_j0p0_exz0_e0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_j0p0_exz_E0[1]= new TH2F( "xy25m_13m_j0p0_exz1_e0","xy25m_13m_j0p0_exz1_e0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_j0p0_exz_E0[2]= new TH2F( "xy25m_13m_j0p0_exz2_e0","xy25m_13m_j0p0_exz2_e0", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_j0p0_exz_E0[3]= new TH2F( "xy25m_13m_j0p0_exz3_e0","xy25m_13m_j0p0_exz3_e0", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_j0p0_exz_E1[4];
  xy25m_13m_j0p0_exz_E1[0]= new TH2F( "xy25m_13m_j0p0_exz0_e1","xy25m_13m_j0p0_exz0_e1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_j0p0_exz_E1[1]= new TH2F( "xy25m_13m_j0p0_exz1_e1","xy25m_13m_j0p0_exz1_e1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_j0p0_exz_E1[2]= new TH2F( "xy25m_13m_j0p0_exz2_e1","xy25m_13m_j0p0_exz2_e1", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_j0p0_exz_E1[3]= new TH2F( "xy25m_13m_j0p0_exz3_e1","xy25m_13m_j0p0_exz3_e1", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_j0p0_exz_E2[4];
  xy25m_13m_j0p0_exz_E2[0]= new TH2F( "xy25m_13m_j0p0_exz0_e2","xy25m_13m_j0p0_exz0_e2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_j0p0_exz_E2[1]= new TH2F( "xy25m_13m_j0p0_exz1_e2","xy25m_13m_j0p0_exz1_e2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_j0p0_exz_E2[2]= new TH2F( "xy25m_13m_j0p0_exz2_e2","xy25m_13m_j0p0_exz2_e2", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_j0p0_exz_E2[3]= new TH2F( "xy25m_13m_j0p0_exz3_e2","xy25m_13m_j0p0_exz3_e2", 100, -50, 50, 100, -50, 50 );

  TH2F *xy25m_13m_j0p0_exz_E3[4];
  xy25m_13m_j0p0_exz_E3[0]= new TH2F( "xy25m_13m_j0p0_exz0_e3","xy25m_13m_j0p0_exz0_e3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_j0p0_exz_E3[1]= new TH2F( "xy25m_13m_j0p0_exz1_e3","xy25m_13m_j0p0_exz1_e3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_j0p0_exz_E3[2]= new TH2F( "xy25m_13m_j0p0_exz2_e3","xy25m_13m_j0p0_exz2_e3", 100, -50, 50, 100, -50, 50 );
  xy25m_13m_j0p0_exz_E3[3]= new TH2F( "xy25m_13m_j0p0_exz3_e3","xy25m_13m_j0p0_exz3_e3", 100, -50, 50, 100, -50, 50 );*/

  TH2F *ax_xj1_25_13 = new TH2F( "ax_xj1_25_13","ax_xj1_25_13", 100, -50, 50, 100, -0.05, 0.05 );
  TH2F *ay_yj1_25_13 = new TH2F( "ay_yj1_25_13","ay_yj1_25_13", 100, -50, 50, 100, -0.05, 0.05 );
  TH2F *ax_xj1_26_14 = new TH2F( "ax_xj1_26_14","ax_xj1_26_14", 100, -50, 50, 100, -0.05, 0.05 );
  TH2F *ay_yj1_26_14 = new TH2F( "ay_yj1_26_14","ay_yj1_26_14", 100, -50, 50, 100, -0.05, 0.05 );

  TH2F *tdc_j0j1_x = new TH2F( "tdc_j0j1_x","tdc_j0j1_x", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc_j0j1_y = new TH2F( "tdc_j0j1_y","tdc_j0j1_y", 384, 2580, 6420, 384, 2580, 6420 );

  TH2F *tdc_j0mj1m_x = new TH2F( "tdc_j0mj1m_x","tdc_j0mj1m_x", 384, 2580, 6420, 384, 2580, 6420 );
  TH2F *tdc_j0mj1m_y = new TH2F( "tdc_j0mj1m_y","tdc_j0mj1m_y", 384, 2580, 6420, 384, 2580, 6420 );

  TH2F *LG_tdc_j0x = new TH2F( "LG_tdc_j0x","LG_tdc_j0x", 384, 2580, 6420, 320, 0, 320 );

  TH2F *xy_LG= new TH2F( "xy_LG","xy_LG", 100, -100, 100, 100, -100, 100 );
  TH2F *x_LG= new TH2F( "x_LG","x_LG", 100, -100, 100, 80, 0, 320 );
  TH2F *y_LG= new TH2F( "y_LG","y_LG", 100, -100, 100, 80, 0, 320 );
  TH1F *xmode_LG= new TH1F( "xmode_LG","xmode_LG", 80, 0, 320 );
  TH2F *x_intLGpos = new TH2F( "x_intLGpos", "x_intLGpos", 100, -100, 100, 240, 0, 960 );
  TH2F *xw_intLGpos= new TH2F( "xw_LG","xw_LG", 100, -200, 200, 240, 0, 960 );
  
  TH2F *xy_LG_E[4];
  xy_LG_E[0]= new TH2F( "xy_LG0","xy_LG0", 100, -100, 100, 100, -100, 100 );
  xy_LG_E[1]= new TH2F( "xy_LG1","xy_LG1", 100, -100, 100, 100, -100, 100 );
  xy_LG_E[2]= new TH2F( "xy_LG2","xy_LG2", 100, -100, 100, 100, -100, 100 );
  xy_LG_E[3]= new TH2F( "xy_LG3","xy_LG3", 100, -100, 100, 100, -100, 100 );

  TH2F *xy_LG_j1xneg_E[4];
  xy_LG_j1xneg_E[0]= new TH2F( "xy_LG_j1xneg_0","xy_LG_j1xneg_0", 100, -100, 100, 100, -100, 100 );
  xy_LG_j1xneg_E[1]= new TH2F( "xy_LG_j1xneg_1","xy_LG_j1xneg_1", 100, -100, 100, 100, -100, 100 );
  xy_LG_j1xneg_E[2]= new TH2F( "xy_LG_j1xneg_2","xy_LG_j1xneg_2", 100, -100, 100, 100, -100, 100 );
  xy_LG_j1xneg_E[3]= new TH2F( "xy_LG_j1xneg_3","xy_LG_j1xneg_3", 100, -100, 100, 100, -100, 100 );

  /*TH2F *xy_LG_j0p0_E[4];
  xy_LG_j0p0_E[0]= new TH2F( "xy_LG_j0p0_E0","xy_LG_j0p0_E0", 100, -100, 100, 100, -100, 100 );
  xy_LG_j0p0_E[1]= new TH2F( "xy_LG_j0p0_E1","xy_LG_j0p0_E1", 100, -100, 100, 100, -100, 100 );
  xy_LG_j0p0_E[2]= new TH2F( "xy_LG_j0p0_E2","xy_LG_j0p0_E2", 100, -100, 100, 100, -100, 100 );
  xy_LG_j0p0_E[3]= new TH2F( "xy_LG_j0p0_E3","xy_LG_j0p0_E3", 100, -100, 100, 100, -100, 100 );*/

  TH2F *x25m_j0w1_LG= new TH2F( "x25m_j0w1_LG","x25m_j0w1_LG", 100, -50, 50, 80, 0, 320 );
  TH2F *y13m_j0w1_LG= new TH2F( "y13m_j0w1_LG","y13m_j0w1_LG", 100, -50, 50, 80, 0, 320 );
  TH2F *x25m_j1w0_LG= new TH2F( "x25m_j1w0_LG","x25m_j1w0_LG", 100, -50, 50, 80, 0, 320 );
  TH2F *y13m_j1w0_LG= new TH2F( "y13m_j1w0_LG","y13m_j1w0_LG", 100, -50, 50, 80, 0, 320 );
  TH2F *x25m_j0n1_LG= new TH2F( "x25m_j0n1_LG","x25m_j0n1_LG", 100, -50, 50, 80, 0, 320 );
  TH2F *y13m_j0n1_LG= new TH2F( "y13m_j0n1_LG","y13m_j0n1_LG", 100, -50, 50, 80, 0, 320 );
  TH2F *x25m_j1n0_LG= new TH2F( "x25m_j1n0_LG","x25m_j1n0_LG", 100, -50, 50, 80, 0, 320 );
  TH2F *y13m_j1n0_LG= new TH2F( "y13m_j1n0_LG","y13m_j1n0_LG", 100, -50, 50, 80, 0, 320 );

  TH2F *xy25m_13mj0_E[13];
  TH2F *xy25m_13mj0n1_E[13];
  TH2F *xy25m_13mj0w1_E[13];
  TH2F *xy25m_13mj1_E[13];
  TH2F *xy25m_13mj1n0_E[13];
  TH2F *xy25m_13mj1w0_E[13];

  xy25m_13mj0_E[0] = new TH2F( "xy25m_13mj0_E00","xy25m_13mj0_E00", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0_E[1] = new TH2F( "xy25m_13mj0_E01","xy25m_13mj0_E01", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0_E[2] = new TH2F( "xy25m_13mj0_E02","xy25m_13mj0_E02", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0_E[3] = new TH2F( "xy25m_13mj0_E03","xy25m_13mj0_E03", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0_E[4] = new TH2F( "xy25m_13mj0_E04","xy25m_13mj0_E04", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0_E[5] = new TH2F( "xy25m_13mj0_E05","xy25m_13mj0_E05", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0_E[6] = new TH2F( "xy25m_13mj0_E06","xy25m_13mj0_E06", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0_E[7] = new TH2F( "xy25m_13mj0_E07","xy25m_13mj0_E07", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0_E[8] = new TH2F( "xy25m_13mj0_E08","xy25m_13mj0_E08", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0_E[9] = new TH2F( "xy25m_13mj0_E09","xy25m_13mj0_E09", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0_E[10] = new TH2F( "xy25m_13mj0_E10","xy25m_13mj0_E10", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0_E[11] = new TH2F( "xy25m_13mj0_E11","xy25m_13mj0_E11", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0_E[12] = new TH2F( "xy25m_13mj0_E12","xy25m_13mj0_E12", 100, -50, 50, 100, -50, 50 );

  xy25m_13mj0n1_E[0] = new TH2F( "xy25m_13mj0n1_E00","xy25m_13mj0n1_E00", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0n1_E[1] = new TH2F( "xy25m_13mj0n1_E01","xy25m_13mj0n1_E01", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0n1_E[2] = new TH2F( "xy25m_13mj0n1_E02","xy25m_13mj0n1_E02", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0n1_E[3] = new TH2F( "xy25m_13mj0n1_E03","xy25m_13mj0n1_E03", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0n1_E[4] = new TH2F( "xy25m_13mj0n1_E04","xy25m_13mj0n1_E04", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0n1_E[5] = new TH2F( "xy25m_13mj0n1_E05","xy25m_13mj0n1_E05", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0n1_E[6] = new TH2F( "xy25m_13mj0n1_E06","xy25m_13mj0n1_E06", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0n1_E[7] = new TH2F( "xy25m_13mj0n1_E07","xy25m_13mj0n1_E07", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0n1_E[8] = new TH2F( "xy25m_13mj0n1_E08","xy25m_13mj0n1_E08", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0n1_E[9] = new TH2F( "xy25m_13mj0n1_E09","xy25m_13mj0n1_E09", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0n1_E[10] = new TH2F( "xy25m_13mj0n1_E10","xy25m_13mj0n1_E10", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0n1_E[11] = new TH2F( "xy25m_13mj0n1_E11","xy25m_13mj0n1_E11", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0n1_E[12] = new TH2F( "xy25m_13mj0n1_E12","xy25m_13mj0n1_E12", 100, -50, 50, 100, -50, 50 );
  
  xy25m_13mj0w1_E[0] = new TH2F( "xy25m_13mj0w1_E00","xy25m_13mj0w1_E00", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_E[1] = new TH2F( "xy25m_13mj0w1_E01","xy25m_13mj0w1_E01", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_E[2] = new TH2F( "xy25m_13mj0w1_E02","xy25m_13mj0w1_E02", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_E[3] = new TH2F( "xy25m_13mj0w1_E03","xy25m_13mj0w1_E03", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_E[4] = new TH2F( "xy25m_13mj0w1_E04","xy25m_13mj0w1_E04", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_E[5] = new TH2F( "xy25m_13mj0w1_E05","xy25m_13mj0w1_E05", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_E[6] = new TH2F( "xy25m_13mj0w1_E06","xy25m_13mj0w1_E06", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_E[7] = new TH2F( "xy25m_13mj0w1_E07","xy25m_13mj0w1_E07", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_E[8] = new TH2F( "xy25m_13mj0w1_E08","xy25m_13mj0w1_E08", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_E[9] = new TH2F( "xy25m_13mj0w1_E09","xy25m_13mj0w1_E09", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_E[10] = new TH2F( "xy25m_13mj0w1_E10","xy25m_13mj0w1_E10", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_E[11] = new TH2F( "xy25m_13mj0w1_E11","xy25m_13mj0w1_E11", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_E[12] = new TH2F( "xy25m_13mj0w1_E12","xy25m_13mj0w1_E12", 100, -50, 50, 100, -50, 50 );
  
  xy25m_13mj1_E[0] = new TH2F( "xy25m_13mj1_E00","xy25m_13mj1_E00", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1_E[1] = new TH2F( "xy25m_13mj1_E01","xy25m_13mj1_E01", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1_E[2] = new TH2F( "xy25m_13mj1_E02","xy25m_13mj1_E02", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1_E[3] = new TH2F( "xy25m_13mj1_E03","xy25m_13mj1_E03", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1_E[4] = new TH2F( "xy25m_13mj1_E04","xy25m_13mj1_E04", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1_E[5] = new TH2F( "xy25m_13mj1_E05","xy25m_13mj1_E05", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1_E[6] = new TH2F( "xy25m_13mj1_E06","xy25m_13mj1_E06", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1_E[7] = new TH2F( "xy25m_13mj1_E07","xy25m_13mj1_E07", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1_E[8] = new TH2F( "xy25m_13mj1_E08","xy25m_13mj1_E08", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1_E[9] = new TH2F( "xy25m_13mj1_E09","xy25m_13mj1_E09", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1_E[10] = new TH2F( "xy25m_13mj1_E10","xy25m_13mj1_E10", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1_E[11] = new TH2F( "xy25m_13mj1_E11","xy25m_13mj1_E11", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1_E[12] = new TH2F( "xy25m_13mj1_E12","xy25m_13mj1_E12", 100, -50, 50, 100, -50, 50 );

  xy25m_13mj1n0_E[0] = new TH2F( "xy25m_13mj1n0_E00","xy25m_13mj1n0_E00", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1n0_E[1] = new TH2F( "xy25m_13mj1n0_E01","xy25m_13mj1n0_E01", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1n0_E[2] = new TH2F( "xy25m_13mj1n0_E02","xy25m_13mj1n0_E02", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1n0_E[3] = new TH2F( "xy25m_13mj1n0_E03","xy25m_13mj1n0_E03", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1n0_E[4] = new TH2F( "xy25m_13mj1n0_E04","xy25m_13mj1n0_E04", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1n0_E[5] = new TH2F( "xy25m_13mj1n0_E05","xy25m_13mj1n0_E05", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1n0_E[6] = new TH2F( "xy25m_13mj1n0_E06","xy25m_13mj1n0_E06", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1n0_E[7] = new TH2F( "xy25m_13mj1n0_E07","xy25m_13mj1n0_E07", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1n0_E[8] = new TH2F( "xy25m_13mj1n0_E08","xy25m_13mj1n0_E08", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1n0_E[9] = new TH2F( "xy25m_13mj1n0_E09","xy25m_13mj1n0_E09", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1n0_E[10] = new TH2F( "xy25m_13mj1n0_E10","xy25m_13mj1n0_E10", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1n0_E[11] = new TH2F( "xy25m_13mj1n0_E11","xy25m_13mj1n0_E11", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1n0_E[12] = new TH2F( "xy25m_13mj1n0_E12","xy25m_13mj1n0_E12", 100, -50, 50, 100, -50, 50 );
  
  xy25m_13mj1w0_E[0] = new TH2F( "xy25m_13mj1w0_E00","xy25m_13mj1w0_E00", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_E[1] = new TH2F( "xy25m_13mj1w0_E01","xy25m_13mj1w0_E01", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_E[2] = new TH2F( "xy25m_13mj1w0_E02","xy25m_13mj1w0_E02", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_E[3] = new TH2F( "xy25m_13mj1w0_E03","xy25m_13mj1w0_E03", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_E[4] = new TH2F( "xy25m_13mj1w0_E04","xy25m_13mj1w0_E04", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_E[5] = new TH2F( "xy25m_13mj1w0_E05","xy25m_13mj1w0_E05", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_E[6] = new TH2F( "xy25m_13mj1w0_E06","xy25m_13mj1w0_E06", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_E[7] = new TH2F( "xy25m_13mj1w0_E07","xy25m_13mj1w0_E07", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_E[8] = new TH2F( "xy25m_13mj1w0_E08","xy25m_13mj1w0_E08", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_E[9] = new TH2F( "xy25m_13mj1w0_E09","xy25m_13mj1w0_E09", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_E[10] = new TH2F( "xy25m_13mj1w0_E10","xy25m_13mj1w0_E10", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_E[11] = new TH2F( "xy25m_13mj1w0_E11","xy25m_13mj1w0_E11", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_E[12] = new TH2F( "xy25m_13mj1w0_E12","xy25m_13mj1w0_E12", 100, -50, 50, 100, -50, 50 );
  
  TH2F *xy25m_13mj0w1_cE[4];
  TH2F *xy25m_13mj1w0_cE[4];

  xy25m_13mj0w1_cE[0] = new TH2F( "xy25m_13mj0w1_cE0","xy25m_13mj0w1_cE0", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_cE[1] = new TH2F( "xy25m_13mj0w1_cE1","xy25m_13mj0w1_cE1", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_cE[2] = new TH2F( "xy25m_13mj0w1_cE2","xy25m_13mj0w1_cE2", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1_cE[3] = new TH2F( "xy25m_13mj0w1_cE3","xy25m_13mj0w1_cE3", 100, -50, 50, 100, -50, 50 );
  
  xy25m_13mj1w0_cE[0] = new TH2F( "xy25m_13mj1w0_cE0","xy25m_13mj1w0_cE0", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_cE[1] = new TH2F( "xy25m_13mj1w0_cE1","xy25m_13mj1w0_cE1", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_cE[2] = new TH2F( "xy25m_13mj1w0_cE2","xy25m_13mj1w0_cE2", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0_cE[3] = new TH2F( "xy25m_13mj1w0_cE3","xy25m_13mj1w0_cE3", 100, -50, 50, 100, -50, 50 );

  /*TH2F *xy25m_13mj0w1j0p0_cE[4];
  TH2F *xy25m_13mj1w0j0p0_cE[4];

  xy25m_13mj0w1j0p0_cE[0] = new TH2F( "xy25m_13mj0w1j0p0_cE0","xy25m_13mj0w1j0p0_cE0", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1j0p0_cE[1] = new TH2F( "xy25m_13mj0w1j0p0_cE1","xy25m_13mj0w1j0p0_cE1", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1j0p0_cE[2] = new TH2F( "xy25m_13mj0w1j0p0_cE2","xy25m_13mj0w1j0p0_cE2", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj0w1j0p0_cE[3] = new TH2F( "xy25m_13mj0w1j0p0_cE3","xy25m_13mj0w1j0p0_cE3", 100, -50, 50, 100, -50, 50 );
  
  xy25m_13mj1w0j0p0_cE[0] = new TH2F( "xy25m_13mj1w0j0p0_cE0","xy25m_13mj1w0j0p0_cE0", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0j0p0_cE[1] = new TH2F( "xy25m_13mj1w0j0p0_cE1","xy25m_13mj1w0j0p0_cE1", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0j0p0_cE[2] = new TH2F( "xy25m_13mj1w0j0p0_cE2","xy25m_13mj1w0j0p0_cE2", 100, -50, 50, 100, -50, 50 );
  xy25m_13mj1w0j0p0_cE[3] = new TH2F( "xy25m_13mj1w0j0p0_cE3","xy25m_13mj1w0j0p0_cE3", 100, -50, 50, 100, -50, 50 );*/

  TH1F *tdcL_E0[Nch];
  tdcL_E0[0] = new TH1F( "tdcL_E013j0", "tdcL_E013j0", 1500, 0, 15000 );
  tdcL_E0[1] = new TH1F( "tdcL_E014j0", "tdcL_E014j0", 1500, 0, 15000 );
  tdcL_E0[2] = new TH1F( "tdcL_E025j0", "tdcL_E025j0", 1500, 0, 15000 );
  tdcL_E0[3] = new TH1F( "tdcL_E026j0", "tdcL_E026j0", 1500, 0, 15000 );
  tdcL_E0[4] = new TH1F( "tdcL_E013j1", "tdcL_E013j1", 1500, 0, 15000 );
  tdcL_E0[5] = new TH1F( "tdcL_E014j1", "tdcL_E014j1", 1500, 0, 15000 );
  tdcL_E0[6] = new TH1F( "tdcL_E025j1", "tdcL_E025j1", 1500, 0, 15000 );
  tdcL_E0[7] = new TH1F( "tdcL_E026j1", "tdcL_E026j1", 1500, 0, 15000 );
  TH1F *tdcL_En0[Nch];
  tdcL_En0[0] = new TH1F( "tdcL_En013j0", "tdcL_En013j0", 1500, 0, 15000 );
  tdcL_En0[1] = new TH1F( "tdcL_En014j0", "tdcL_En014j0", 1500, 0, 15000 );
  tdcL_En0[2] = new TH1F( "tdcL_En025j0", "tdcL_En025j0", 1500, 0, 15000 );
  tdcL_En0[3] = new TH1F( "tdcL_En026j0", "tdcL_En026j0", 1500, 0, 15000 );
  tdcL_En0[4] = new TH1F( "tdcL_En013j1", "tdcL_En013j1", 1500, 0, 15000 );
  tdcL_En0[5] = new TH1F( "tdcL_En014j1", "tdcL_En014j1", 1500, 0, 15000 );
  tdcL_En0[6] = new TH1F( "tdcL_En025j1", "tdcL_En025j1", 1500, 0, 15000 );
  tdcL_En0[7] = new TH1F( "tdcL_En026j1", "tdcL_En026j1", 1500, 0, 15000 );

  TGraph *linex;
  TGraph *liney;

  const Int_t Nq = 4;
  //Float_t xcp = 43., xcn = -45., ycp = 40., ycn = -47., rp = 30.;
  Float_t xcp = 39., xcn = -45., ycp = 40., ycn = -44., rp = 30.;
  TArc *Q[Nq];
  Q[0] = new TArc( xcp, ycp, rp, 180, 270 );
  Q[1] = new TArc( xcp, ycn, rp, 90, 180 );
  Q[2] = new TArc( xcn, ycn, rp, 0, 90 );
  Q[3] = new TArc( xcn, ycp, rp, 270, 360 );
  for( Int_t i = 0; i < Nq; i++ ){
    //Q[i]->SetLineWidth(2);
    //Q[i]->SetLineColor(2);
    //Q[i]->SetFillStyle(0);
    Q[i]->SetFillColor(16);
  }
  Float_t r0, r1, r2, r3;

  TArrow *path;

  Int_t ncellx0 = 0;
  Int_t ncelly0 = 0;
  Int_t ncellx1 = 0;
  Int_t ncelly1 = 0;

  Float_t fp0x, fp1x, fp0y, fp1y;

  Float_t dx, dy;

  Float_t x0, y0, x1, y1;

  gRandom->SetSeed(0);

  Int_t Ndraw = 0;
  Int_t Nproc = 0;

  TCanvas *c1 = new TCanvas( "c1", "Canvas", 0, 0, 800, 800 );
  TCanvas *c2 = new TCanvas( "c2", "Canvas2", 0, 0, 800, 800 );
  TCanvas *c3 = new TCanvas( "c3", "Canvas3", 0, 0, 800, 800 );
  TCanvas *c4 = new TCanvas( "c4", "Canvas4", 0, 0, 800, 800 );
  
  c1->Divide(2,2);
  c2->Divide(3,3);
  c3->Divide(3,3);
  c4->Divide(1,1);

  gStyle->SetOptFit();
  gStyle->SetOptStat(111111);
  
  c4->cd(1);
  TH1F *frame = gPad->DrawFrame( -50, -50, 50, 50 );
  frame->GetXaxis()->SetTitle( "x (mm)" );
  frame->GetYaxis()->SetTitle( "y (mm)" );
  for( Int_t i = 0; i < Nq; i++ ){
    Q[i]->Draw("same");
  }

  TGraph *wf;

  Int_t Nwf = 0;
  Int_t peak;
  Int_t peak2;
  Int_t valley;
  Int_t ppeak, ppeak2;
  Int_t pvall;

  Int_t data_size = tree0->GetEntries();
  if( full == -1 ) data_size = 50000;
  if( full == 0 ) data_size = 200000;
  Int_t data_size_1 = tree1->GetEntries();
  Int_t data_size_2 = tree2->GetEntries();

  //stat->Fill( 0., data_size );
  //stat->Fill( 1., data_size_1 );
  //stat->Fill( 2., data_size_2 );
  stat->SetBinContent(1, data_size-data_size );
  stat->SetBinContent(2, data_size_1-data_size );
  stat->SetBinContent(3, data_size_2-data_size );
  
  const Int_t interval = data_size/200;
  
  Int_t j1tdc26;

  Int_t goodj0y13 = 0;
  Int_t goodj0y14 = 0;
  Int_t goodj0x25 = 0;
  Int_t goodj0x26 = 0;
  Int_t goodj1y13 = 0;
  Int_t goodj1y14 = 0;
  Int_t goodj1x25 = 0;
  Int_t goodj1x26 = 0;
  Int_t goodj0yy = 0;
  Int_t goodj0xx = 0;
  Int_t goodj1yy = 0;
  Int_t goodj1xx = 0;

  Int_t vgoodj0y13 = 0;
  Int_t vgoodj0x25 = 0;
  Int_t vgoodj1y13 = 0;
  Int_t vgoodj1x25 = 0;

  Int_t dtdcth = 30;

  Int_t integral;
  Int_t intpos = 0;
  Int_t intpos1;
  Int_t intpos2;
  Int_t intlen;
  Int_t intlen1;
  Int_t intlen2;
  
  //////////////////////////////////////////////////
  for( Int_t ii = 0; ii < N_files; ii++ ){
  for( Long64_t h = 0; h < data_size; h++ ){
    //for( Long64_t i = 0; i < 500000; i++ ){
    //cout << "Load 0 " << 
    tree0->LoadTree(h);
    //cout << "Load 1 " << 
    tree1->LoadTree(h);
    if( LG_on )
      tree2->LoadTree(h);

    tree0->GetEntry(h);
    tree1->GetEntry(h);
    if( LG_on )
      tree2->GetEntry(h);

    peak = 1024;
    peak2 = 1024;
    valley = 512;
    ppeak = -1;
    ppeak2 = -1;
    pvall = -1;
    if( LG_on ){
      length = f2_adc30->size();
      LGlength->Fill( length );
      intlen = 0;
      intlen1 = 0;
      intlen2 = 0;
      if( length > 0 ){
	peak = 1024;
	peak2 = 1024;
	valley = 512;
	// peak search
	integral = 0;
	intpos = 0;
	for( Int_t i = 0; i < length; i++ ){
	  yy[i] = f2_adc30->at(i);
	  if( yy[i] < peak ){ peak = yy[i]; ppeak = i; }	  
	  integral = integral + (512 - yy[i] );//( 512 - yy[i] );
	  if( yy[i] < 512 ){
	    intpos = intpos + ( 512 - yy[i] );
	    intlen++;
	  }
	}
	// valley search
	if( ppeak > -1 && ppeak < length-1){
	  for( Int_t i = ppeak+1; i < length; i++ ){
	    //if( 512 - yy[i] < 0 ){ pvall = i; break;}
	    if( yy[i] > valley ){ valley = yy[i]; pvall = i;}
	  }
	}
	// second peak search
	if( pvall > -1 && pvall < length-1 ){
	  for( Int_t i = pvall; i < length; i++ ){
	    if( yy[i] < peak2 ){ peak2 = yy[i]; ppeak2 = i; }	  
	  }
	}
	
	intpos1 = 0;
	intpos2 = 0;
	if( pvall > -1 && peak2 > -1 ){
	  for( Int_t i = 0; i < pvall; i++ ){
	    if( yy[i] < 512 ){
	      intpos1 = intpos1 + ( 512 - yy[i] );
	      intlen1++;
	    }
	  }
	  for( Int_t i = pvall; i < length; i++ ){
	    if( yy[i] < 512 ){
	      intpos2 = intpos2 + ( 512 - yy[i] );
	      intlen2++;
	    }
	  }
	}
	
	LGintlength->Fill( intlen );
	if( intlen1 != 0 ) LGintlength1->Fill( intlen1 );
	if( intlen2 != 0 ) LGintlength2->Fill( intlen2 );
	LG->Fill( 512 - peak );
    	LG_converted->Fill((float(512-peak)-2.27292)/47.4055);
	intLG->Fill( integral );
	intLGpos->Fill( intpos );
	intLGpos1->Fill( intpos1 );
	intLGpos2->Fill( intpos2 );
	LGpeak->Fill( ppeak );
	if( pvall > -1 ){
	  LGvalley->Fill( pvall );
	}
	if( ppeak2 > -1 ){
	  LG2->Fill( 512 - peak2 );
	  LGpeak2->Fill( ppeak2 );
	}
	if( integral > 50 ) intLGth->Fill( integral );
	// draw wave form
	if( Nwf < 100 ){
	  wf = new TGraph( length, xx, yy );
	  c2->cd(Nwf%9+1);
	  wf->Draw( "AL" );
	  if( Nwf%9 == 8 ){
	    // if( Nwf == 8 ) c2->Print("waveform.pdf(");
	    // else c2->Print("waveform.pdf");
	    if( Nwf == 8 ) pdf_waveform = "/group/itdc/tbl/2024ARTBL018/pdf/waveform/waveform_run" + runnumber + ".pdf(";
	    else pdf_waveform = "/group/itdc/tbl/2024ARTBL018/pdf/waveform/waveform_run" + runnumber + ".pdf";
	    c2->Print(pdf_waveform);
	  }
	}
	Nwf++;
      }
      //if( !(Nwf%1000) ) cout << "#wave form check " << Nwf << "   h = " << h << endl;
      if( !(Nwf%interval) ) cout << "#wave form check " << Nwf << "   h = " << h << endl;
      length = f2_tdc30->size();
      if( length )
	LGtdc->Fill( f2_tdc30->at(0) );
    }// LG on

    //if( tree0->LoadTree(i) > 0 ){
    Nstop[0] = f0_tdc13->size(); 
    Nstop[1] = f0_tdc14->size(); 
    Nstop[2] = f0_tdc25->size(); 
    Nstop[3] = f0_tdc26->size(); 
    Nstop[4] = f1_tdc13->size(); 
    Nstop[5] = f1_tdc14->size(); 
    Nstop[6] = f1_tdc25->size(); 
    Nstop[7] = f1_tdc26->size(); 

    // check #hit cells
    ncellx0 = 0;
    ncelly0 = 0;
    ncellx1 = 0;
    ncelly1 = 0;
    if( Nstop[0] > 0 ) ncelly0++;
    if( Nstop[1] > 0 ) ncelly0++;
    if( Nstop[2] > 0 ) ncellx0++;
    if( Nstop[3] > 0 ) ncellx0++;
    if( Nstop[4] > 0 ) ncelly1++;
    if( Nstop[5] > 0 ) ncelly1++;
    if( Nstop[6] > 0 ) ncellx1++;
    if( Nstop[7] > 0 ) ncellx1++;

    ncell[0]->Fill( ncellx0, ncelly0 );
    ncell[1]->Fill( ncellx1, ncelly1 );

    for( Int_t i = 0; i < Nch; i++ ){
      ntdc[i]->Fill( Nstop[i] );
      if( Nstop[i] > 0 ) ntdcn[i]->Fill( Nstop[i] );
    }

    // store tdc count
    for( Int_t i = 0; i < Nch; i++ ){
      for( Int_t j = 0; j < 100; j++ ){
	tdc_Lcount[i][j] = 0;
	tdc_Tcount[i][j] = 0;
      }
      NstopL[i] = 0;
      NstopT[i] = 0;
      for( Int_t j = 0; j < Nstop[i]; j++ ){
	if( i == 0 ){
	  if( f0_tdc13_lt->at(j) == 1 ){
	    tdc_Lcount[i][NstopL[i]] = f0_tdc13->at(j);
	    NstopL[i]++;
	    tdcL[i]->Fill( f0_tdc13->at(j) );
	  }
	  if( f0_tdc13_lt->at(j) == 0 ){
	    tdc_Tcount[i][NstopT[i]] = f0_tdc13->at(j);
	    NstopT[i]++;
	    tdcT[i]->Fill( f0_tdc13->at(j) );
	  }
	  tdc[i]->Fill( f0_tdc13->at(j) );
	}// i = 0
	if( i == 1 ){
	  if( f0_tdc14_lt->at(j) == 1 ){
	    tdc_Lcount[i][NstopL[i]] = f0_tdc14->at(j);
	    NstopL[i]++;
	    tdcL[i]->Fill( f0_tdc14->at(j) );
	  }
	  if( f0_tdc14_lt->at(j) == 0 ){
	    tdc_Tcount[i][NstopT[i]] = f0_tdc14->at(j);
	    NstopT[i]++;
	    tdcT[i]->Fill( f0_tdc14->at(j) );
	  }
	  tdc[i]->Fill( f0_tdc14->at(j) );
	}// i = 1
	if( i == 2 ){
	  if( f0_tdc25_lt->at(j) == 1 ){
	    tdc_Lcount[i][NstopL[i]] = f0_tdc25->at(j);
	    NstopL[i]++;
	    tdcL[i]->Fill( f0_tdc25->at(j) );
	  }
	  if( f0_tdc25_lt->at(j) == 0 ){
	    tdc_Tcount[i][NstopT[i]] = f0_tdc25->at(j);
	    NstopT[i]++;
	    tdcT[i]->Fill( f0_tdc25->at(j) );
	  }
	  tdc[i]->Fill( f0_tdc25->at(j) );
	}// i = 2
	if( i == 3 ){
	  if( f0_tdc26_lt->at(j) == 1 ){
	    tdc_Lcount[i][NstopL[i]] = f0_tdc26->at(j);
	    NstopL[i]++;
	    tdcL[i]->Fill( f0_tdc26->at(j) );
	  }
	  if( f0_tdc26_lt->at(j) == 0 ){
	    tdc_Tcount[i][NstopT[i]] = f0_tdc26->at(j);
	    NstopT[i]++;
	    tdcT[i]->Fill( f0_tdc26->at(j) );
	  }
	  tdc[i]->Fill( f0_tdc26->at(j) );
	}// i = 3
	if( i == 4 ){
	  if( f1_tdc13_lt->at(j) == 1 ){
	    tdc_Lcount[i][NstopL[i]] = f1_tdc13->at(j);
	    NstopL[i]++;
	    tdcL[i]->Fill( f1_tdc13->at(j) );
	  }
	  if( f1_tdc13_lt->at(j) == 0 ){
	    tdc_Tcount[i][NstopT[i]] = f1_tdc13->at(j);
	    NstopT[i]++;
	    tdcT[i]->Fill( f1_tdc13->at(j) );
	  }
	  tdc[i]->Fill( f1_tdc13->at(j) );
	}// i = 4
	if( i == 5 ){
	  if( f1_tdc14_lt->at(j) == 1 ){
	    tdc_Lcount[i][NstopL[i]] = f1_tdc14->at(j);
	    NstopL[i]++;
	    tdcL[i]->Fill( f1_tdc14->at(j) );
	  }
	  if( f1_tdc14_lt->at(j) == 0 ){
	    tdc_Tcount[i][NstopT[i]] = f1_tdc14->at(j);
	    NstopT[i]++;
	    tdcT[i]->Fill( f1_tdc14->at(j) );
	  }
	  tdc[i]->Fill( f1_tdc14->at(j) );
	}// i = 5
	if( i == 6 ){
	  if( f1_tdc25_lt->at(j) == 1 ){
	    tdc_Lcount[i][NstopL[i]] = f1_tdc25->at(j);
	    NstopL[i]++;
	    tdcL[i]->Fill( f1_tdc25->at(j) );
	  }
	  if( f1_tdc25_lt->at(j) == 0 ){
	    tdc_Tcount[i][NstopT[i]] = f1_tdc25->at(j);
	    NstopT[i]++;
	    tdcT[i]->Fill( f1_tdc25->at(j) );
	  }
	  tdc[i]->Fill( f1_tdc25->at(j) );
	}// i = 6
	if( i == 7 ){
	  if( f1_tdc26_lt->at(j) == 1 ){
	    tdc_Lcount[i][NstopL[i]] = f1_tdc26->at(j);
	    NstopL[i]++;
	    tdcL[i]->Fill( f1_tdc26->at(j) );
	  }
	  if( f1_tdc26_lt->at(j) == 0 ){
	    tdc_Tcount[i][NstopT[i]] = f1_tdc26->at(j);
	    NstopT[i]++;
	    tdcT[i]->Fill( f1_tdc26->at(j) );
	  }
	  tdc[i]->Fill( f1_tdc26->at(j) );
	}// i = 7
	/*if( i == 7 ){
	  j1tdc26 = f1_tdc26->at(j);
	  // 0000 0000 0000 0000
	  // 31
	  // 2684 21
	  // 7310 0052 1
	  // 6899 4215 2631
	  // 8426 8426 8426 8421
	  
	  //if( j1tdc26 & 0x0200 ) j1tdc26 = j1tdc26 & 0xfdff;
	  //if( !(j1tdc26 & 0x0200) ) j1tdc26 = j1tdc26 | 0x0200;
	  //if( j1tdc26 & 0x0100 ) j1tdc26 = j1tdc26 & 0xfeff;
	  //if( !(j1tdc26 & 0x0100) ) j1tdc26 = j1tdc26 | 0x0100;
	  if( f1_tdc26_lt->at(j) == 1 ){
	    //tdc_Lcount[i][NstopL[i]] = f1_tdc26->at(j);
	    tdc_Lcount[i][NstopL[i]] = j1tdc26;
	    NstopL[i]++;
	    //tdcL[i]->Fill( f1_tdc26->at(j) );
	    tdcL[i]->Fill( j1tdc26 );
	  }
	  if( f1_tdc26_lt->at(j) == 0 ){
	    //tdc_Tcount[i][NstopT[i]] = f1_tdc26->at(j);
	    tdc_Tcount[i][NstopT[i]] = j1tdc26;
	    NstopT[i]++;
	    //tdcT[i]->Fill( f1_tdc26->at(j) );
	    tdcT[i]->Fill( j1tdc26 );
	  }
	  //tdc[i]->Fill( f1_tdc26->at(j) );
	  tdc[i]->Fill( j1tdc26 );
	}*/// i = 7
      }//j
    }// i
    // end of store tdc count


    // ToT
    for( Int_t i = 0; i < Nch; i++ ){
      ntdcL[i]->Fill( NstopL[i] );
      ntdcT[i]->Fill( NstopT[i] );
      if( NstopL[i] == 1 ) tdcLs[i]->Fill( tdc_Lcount[i][0] );
      if( NstopL[i] == 1 && NstopT[i] == 1 ){
	totsngl[i]->Fill( tdc_Lcount[i][0] - tdc_Tcount[i][0] );
	totsnglf[i]->Fill( tdc_Lcount[i][0] - tdc_Tcount[i][0] );
      }

      if( NstopL[i] < NstopT[i] ){
	for( Int_t j = 0; j < NstopL[i]; j++ ){
	  tot_count[i][j] = tdc_Lcount[i][j] - tdc_Tcount[i][j];
	  tot[i]->Fill( tot_count[i][j] );
	  tot_tdcL[i]->Fill( tdc_Lcount[i][j], tot_count[i][j] );
	  tot_ihit[i]->Fill( j, tot_count[i][j] );
	  if( j == 0 ) tot_tdcL_0[i]->Fill( tdc_Lcount[i][j], tot_count[i][j] );
	  else tot_tdcL_n[i]->Fill( tdc_Lcount[i][j], tot_count[i][j] );
	}
	tot1st[i]->Fill( tot_count[i][NstopL[i]-1] );
      }
      
      tdcLm[i]->Fill( tdc_Lcount[i][NstopL[i]-1] );
    }// ToT

    if( LG_on ){
      totvsPH[0]->Fill( 512-peak, tot_count[0][0] );
      totvsPH[1]->Fill( 512-peak, tot_count[2][0] );
      totvsPH[2]->Fill( 512-peak, tot_count[4][0] );
      totvsPH[3]->Fill( 512-peak, tot_count[6][0] );
    }
    
    // TDC count 2D map (all hits)
    for( Int_t i = 0; i < NstopL[2]; i++ ){
      for( Int_t j = 0; j < NstopL[0]; j++ ){
	tdc25_13j0->Fill( tdc_Lcount[2][i], tdc_Lcount[0][j] );
      }
    }
    for( Int_t i = 0; i < NstopL[6]; i++ ){
      for( Int_t j = 0; j < NstopL[4]; j++ ){
	tdc25_13j1->Fill( tdc_Lcount[6][i], tdc_Lcount[4][j] );
      }
    }// 2D map

    goodj0y13 = 0;
    goodj0y14 = 0;
    goodj0x25 = 0;
    goodj0x26 = 0;
    goodj1y13 = 0;
    goodj1y14 = 0;
    goodj1x25 = 0;
    goodj1x26 = 0;
    goodj0yy = 0;
    goodj0xx = 0;
    goodj1yy = 0;
    goodj1xx = 0;

    vgoodj0y13 = 0;
    vgoodj0x25 = 0;
    vgoodj1y13 = 0;
    vgoodj1x25 = 0;

    /*if( NstopL[0] == 1 && NstopT[0] == 1 && tdc_Lcount[0][0] - tdc_Tcount[0][0] > 5 && tdc_Lcount[0][0] - tdc_Tcount[0][0] < 1000  ) goodj0y13 = 1;
    if( NstopL[1] == 1 && NstopT[1] == 1 && tdc_Lcount[1][0] - tdc_Tcount[1][0] > 5 && tdc_Lcount[1][0] - tdc_Tcount[1][0] < 1000  ) goodj0y14 = 1;
    if( NstopL[2] == 1 && NstopT[2] == 1 && tdc_Lcount[2][0] - tdc_Tcount[2][0] > 5 && tdc_Lcount[2][0] - tdc_Tcount[2][0] < 1000  ) goodj0x25 = 1;
    if( NstopL[3] == 1 && NstopT[3] == 1 && tdc_Lcount[3][0] - tdc_Tcount[3][0] > 5 && tdc_Lcount[3][0] - tdc_Tcount[3][0] < 1000  ) goodj0x26 = 1;
    if( NstopL[4] == 1 && NstopT[4] == 1 && tdc_Lcount[4][0] - tdc_Tcount[4][0] > 5 && tdc_Lcount[4][0] - tdc_Tcount[4][0] < 1000  ) goodj1y13 = 1;
    if( NstopL[5] == 1 && NstopT[5] == 1 && tdc_Lcount[5][0] - tdc_Tcount[5][0] > 5 && tdc_Lcount[5][0] - tdc_Tcount[5][0] < 1000  ) goodj1y14 = 1;
    if( NstopL[6] == 1 && NstopT[6] == 1 && tdc_Lcount[6][0] - tdc_Tcount[6][0] > 5 && tdc_Lcount[6][0] - tdc_Tcount[6][0] < 1000  ) goodj1x25 = 1;
    if( NstopL[7] == 1 && NstopT[7] == 1 && tdc_Lcount[7][0] - tdc_Tcount[7][0] > 5 && tdc_Lcount[7][0] - tdc_Tcount[7][0] < 1000  ) goodj1x26 = 1;*/

    Double_t py13j0 = y13j0->Eval(tdc_Lcount[0][0]);
    Double_t py14j0 = y14j0->Eval(tdc_Lcount[1][0]);
    Double_t px25j0 = x25j0->Eval(tdc_Lcount[2][0]);
    Double_t px26j0 = x26j0->Eval(tdc_Lcount[3][0]);
    Double_t py13j1 = y13j1->Eval(tdc_Lcount[4][0]);
    Double_t py14j1 = y14j1->Eval(tdc_Lcount[5][0]);
    Double_t px25j1 = x25j1->Eval(tdc_Lcount[6][0]);
    Double_t px26j1 = x26j1->Eval(tdc_Lcount[7][0]);

    Double_t fid_x = 44;
    Double_t fid_y = 44;
    if( NstopL[0] > 0 && fabs(py13j0) < fid_y ) goodj0y13 = 1;
    if( NstopL[1] > 0 && fabs(py14j0) < fid_y ) goodj0y14 = 1;
    if( NstopL[2] > 0 && fabs(px25j0) < fid_x ) goodj0x25 = 1;
    if( NstopL[3] > 0 && fabs(px26j0) < fid_x ) goodj0x26 = 1;
    if( NstopL[4] > 0 && fabs(py13j1) < fid_y ) goodj1y13 = 1;
    if( NstopL[5] > 0 && fabs(py14j1) < fid_y ) goodj1y14 = 1;
    if( NstopL[6] > 0 && fabs(px25j1) < fid_x ) goodj1x25 = 1;
    if( NstopL[7] > 0 && fabs(px26j1) < fid_x ) goodj1x26 = 1;

    if( NstopL[0] > 0 && fabs(py13j0) < fid_y/2. ) vgoodj0y13 = 1;
    if( NstopL[2] > 0 && fabs(px25j0) < fid_x/2. ) vgoodj0x25 = 1;
    if( NstopL[4] > 0 && fabs(py13j1) < fid_y/2. ) vgoodj1y13 = 1;
    if( NstopL[6] > 0 && fabs(px25j1) < fid_x/2. ) vgoodj1x25 = 1;

    if( goodj0y13 ){//&& goodj0y14 ){
      //if( fabs( py14j0 - py13j0 + 0.57 ) < 0.7 ) goodj0yy = 1;
      //if( abs( tdc_Lcount[1][0] - tdc_Lcount[0][0] - dtmean[0] ) < 3*dtsigma[0] )
	goodj0yy = 1;
    }
    if( goodj0x25 ){//&& goodj0x26 ){
      //if( fabs( px26j0 - px25j0 - 0.29 ) < 0.7 )
      //if( abs( tdc_Lcount[3][0] - tdc_Lcount[2][0] - dtmean[1] ) < 3*dtsigma[1] )
      goodj0xx = 1;
    }
    if( goodj1y13 ){//&& goodj1y14 ){
      //if( fabs( py14j1 - py13j1 + 0.35 ) < 0.7 ) goodj1yy = 1;
      //if( abs( tdc_Lcount[5][0] - tdc_Lcount[4][0] - dtmean[2] ) < 3*dtsigma[2] )
      goodj1yy = 1;
    }
    if( goodj1x25 ){//&& goodj1x26 ){
      //if( fabs( px26j1 - px25j1 ) < 0.7 ) goodj1xx = 1;
      //if( abs( tdc_Lcount[7][0] - tdc_Lcount[6][0] - dtmean[3] ) < 3*dtsigma[3] )
      goodj1xx = 1;
    }
    //if( goodj1x25 ){
    //goodj1xx = 1;
    //}
    
    if( vgoodj0y13 && vgoodj0x25 ) LGvsEv_j0->Fill( h, 512-peak );
    if( vgoodj1y13 && vgoodj1x25 ) LGvsEv_j1->Fill( h, 512-peak );

    if( goodj0y13 == 1 && goodj0x25 == 1 ){//&& ( tdc_Lcount[2][0] > a[2] || tdc_Lcount[0][0] > a[0] ) ){
      xy25o_13oj0->Fill( px25j0, py13j0 );
      if( ncellx1 > 0 && ncelly1 > 0 )
	xy25o_13oj0w1->Fill( px25j0, py13j0 );
      else
	xy25o_13oj0n1->Fill( px25j0, py13j0 );
    }
    if( goodj0y14 == 1 && goodj0x26 == 1 ){//&& ( tdc_Lcount[3][0] > a[3] || tdc_Lcount[1][0] > a[1] ) ){
      xy26o_14oj0->Fill( px26j0, py14j0 );
      if( ncellx1 > 0 && ncelly1 > 0 )
	xy26o_14oj0w1->Fill( px26j0, py14j0 );
      else
	xy26o_14oj0n1->Fill( px26j0, py14j0 );
    }
    if( goodj1y13 == 1 && goodj1x25 == 1 ){//&& ( tdc_Lcount[6][0] > a[6] || tdc_Lcount[4][0] > a[4] ) ){
      xy25o_13oj1->Fill( px25j1, py13j1 );
      if( ncellx0 > 0 && ncelly0 > 0 )
	xy25o_13oj1w0->Fill( px25j1, py13j1 );
      else
	xy25o_13oj1n0->Fill( px25j1, py13j1 );
    }
    if( goodj1y14 == 1 && goodj1x26 == 1 ){//&& ( tdc_Lcount[7][0] > a[7] || tdc_Lcount[5][0] > a[5] ) ){
      xy26o_14oj1->Fill( px26j1, py14j1 );
      if( ncellx0 > 0 && ncelly0  > 0)
	xy26o_14oj1w0->Fill( px26j1, py14j1 );
      else
      	xy26o_14oj1n0->Fill( px26j1, py14j1 );
    }
    if( LG_on ){
      for( Int_t jj = 0; jj < Nch; jj++ ){
	if( NstopL[jj] > 0 ){
	  if( 512-peak < PHcth[0] )
	    tdcL_E0[jj]->Fill( tdc_Lcount[jj][0] );
	  else
	    tdcL_En0[jj]->Fill( tdc_Lcount[jj][0] );
	}
      }
    }

    // TDC count 2D map & xy 2D map (single good hit)
    // Jet 0
    if( goodj0yy == 1 && goodj0xx == 1 ){
      tdc25m_13mj0->Fill( tdc_Lcount[2][0], tdc_Lcount[0][0] );
      xy25m_13mj0->Fill( px25j0, py13j0 );

      if( LG_on){
	for( Int_t i = 0; i < 13; i++ ){
	  if( (512 - peak ) >= PHth[i] && (512 - peak ) < PHth[i+1] ){
	    xy25m_13mj0_E[i]->Fill( px25j0, py13j0 );
	    break;
	  }
	}
      }
      // no point in Jet 1
      if( ncellx1 == 0 || ncelly1 == 0){
	tdc25m_13mj0n1->Fill( tdc_Lcount[2][0], tdc_Lcount[0][0] );
	xy25m_13mj0n1->Fill( px25j0, py13j0 );

	if( LG_on ){
	  x25m_j0n1_LG->Fill( px25j0, 512- peak );
	  y13m_j0n1_LG->Fill( py13j0, 512- peak );

	  for( Int_t i = 0; i < 13; i++ ){
	    if( (512 - peak ) >= PHth[i] && (512 - peak ) < PHth[i+1] ){
	      xy25m_13mj0n1_E[i]->Fill( px25j0, py13j0 );
	      break;
	    }
	  }
	}
      }else{
	tdc25m_13mj0w1->Fill( tdc_Lcount[2][0], tdc_Lcount[0][0] );
	xy25m_13mj0w1->Fill( px25j0, py13j0 );

      	if( LG_on ){
	  x25m_j0w1_LG->Fill( px25j0, 512- peak );
	  y13m_j0w1_LG->Fill( py13j0, 512- peak );

	  /*for( Int_t i = 0; i < 13; i++ ){
	    if( (512 - peak ) >= PHth[i] && (512 - peak ) < PHth[i+1] ){
	      xy25m_13mj0w1_E[i]->Fill( px25j0, py13j0 );
	      break;
	    }
	  }// for end
	  */
	int LGconvert = int(((float(512-peak)-LG_offset)/LG_slope-LG_plot_offset)/LG_1bin_energy);
	if((float(512-peak)-LG_offset)/LG_slope > 0.2){
	if(LGconvert < 0) LGconvert = 0;
	else if(LGconvert > 9) LGconvert = 9;
	xy25m_13mj0w1_E[LGconvert]->Fill( px25j0, py13j0 );
	}
	  if( (512 - peak ) < PHcth[0] )
	    xy25m_13mj0w1_cE[0]->Fill( px25j0, py13j0 );
	  if( (512 - peak ) >= PHcth[0] && (512 - peak ) < PHcth[1] )
	    xy25m_13mj0w1_cE[1]->Fill( px25j0, py13j0 );
	  if( (512 - peak ) >= PHcth[1] && (512 - peak ) < PHcth[2] )
	    xy25m_13mj0w1_cE[2]->Fill( px25j0, py13j0 );
	  if( (512 - peak ) >= PHcth[2] )
	    xy25m_13mj0w1_cE[3]->Fill( px25j0, py13j0 );
	}// LG on
      }
      tdc_j0mj1m_x->Fill( tdc_Lcount[2][0], tdc_Lcount[6][0] );
    }
    // Jet 1
    if( goodj1yy == 1 && goodj1xx == 1 ){
      tdc25m_13mj1->Fill( tdc_Lcount[6][0], tdc_Lcount[4][0] );
      xy25m_13mj1->Fill( px25j1, py13j1 );

      if( LG_on){
	for( Int_t i = 0; i < 13; i++ ){
	  if( (512 - peak ) >= PHth[i] && (512 - peak ) < PHth[i+1] ){
	    xy25m_13mj1_E[i]->Fill( px25j1, py13j1 );
	    break;
	  }
	}
      }
      // no point in Jet 0
      if( ncellx0 == 0 || ncelly0 == 0){
	tdc25m_13mj1n0->Fill( tdc_Lcount[6][0], tdc_Lcount[4][0] );
	xy25m_13mj1n0->Fill( px25j1, py13j1 );

	if( LG_on ){
	  x25m_j1n0_LG->Fill( px25j1, 512- peak );
	  y13m_j1n0_LG->Fill( py13j1, 512- peak );

	  for( Int_t i = 0; i < 13; i++ ){
	    if( (512 - peak ) >= PHth[i] && (512 - peak ) < PHth[i+1] ){
	      xy25m_13mj1n0_E[i]->Fill( px25j1, py13j1 );
	      break;
	    }
	  }
	}
      }else{
	tdc25m_13mj1w0->Fill( tdc_Lcount[6][0], tdc_Lcount[4][0] );
	xy25m_13mj1w0->Fill( px25j1, py13j1 );

	if( LG_on ){
	  x25m_j1w0_LG->Fill( px25j1, 512- peak );
	  y13m_j1w0_LG->Fill( py13j1, 512- peak );

	  /*for( Int_t i = 0; i < 13; i++ ){
	    if( (512 - peak ) >= PHth[i] && (512 - peak ) < PHth[i+1] ){
	      xy25m_13mj1w0_E[i]->Fill( px25j1, py13j1 );
	      break;
	    }
	  }// for end
	 */
	int LGconvert1 = int(((float(512-peak)-LG_offset)/LG_slope-LG_plot_offset)/LG_1bin_energy);
	if((float(512-peak)-LG_offset)/LG_slope > 0.2){
	if(LGconvert1 < 0) LGconvert1 = 0;
	else if(LGconvert1 > 9) LGconvert1 = 9;
	xy25m_13mj1w0_E[LGconvert1]->Fill( px25j1, py13j1 );
	}
	  if( (512 - peak ) < PHcth[0] )
	    xy25m_13mj1w0_cE[0]->Fill( px25j1, py13j1 );
	  if( (512 - peak ) >= PHcth[0] && (512 - peak ) < PHcth[1] )
	    xy25m_13mj1w0_cE[1]->Fill( px25j1, py13j1 );
	  if( (512 - peak ) >= PHcth[1] && (512 - peak ) < PHcth[2] )
	    xy25m_13mj1w0_cE[2]->Fill( px25j1, py13j1 );
	  if( (512 - peak ) >= PHcth[2] )
	    xy25m_13mj1w0_cE[3]->Fill( px25j1, py13j1 );

	}// LG
      }
      tdc_j0mj1m_y->Fill( tdc_Lcount[0][0], tdc_Lcount[4][0] );
    }

    // hits in ch26 and ch14
    if( goodj0y14 == 1 && goodj0x26 == 1){
      if( ncellx1 > 0 && ncelly1 > 0)
	xy26m_14mj0w1->Fill( px26j0, py14j0 );
    }
    if( goodj1y14 == 1 && goodj1x26 == 1){
      if( ncellx0 > 0 && ncelly0 > 0)
	xy26m_14mj1w0->Fill( px26j1, py14j1 );
    }
    // 2D map

    // Q-pole effects
    if( goodj0yy == 1 && goodj0xx == 1 && goodj1yy == 1 && goodj1xx == 1 ){
      xy25m_13mj0w1_25_13->Fill( px25j0, py13j0 );
      xy25m_13mj1w0_25_13->Fill( px25j1, py13j1 );
    }
    if( goodj0y14 == 1 && goodj0x26 && goodj1y14 == 1 && goodj1x26 == 1 ){
      xy26m_14mj0w1_26_14->Fill( px26j0, py14j0 );
      xy26m_14mj1w0_26_14->Fill( px26j1, py14j1 );
    }
      
    // hits in pole shade
    // Jet 0
    if( goodj0yy == 1 && goodj0xx == 1 && goodj1yy == 1 && goodj1xx == 1 ){
      x1 = px25j1;
      y1 = py13j1;
      r0 = sqrt( (x1 - xcp)*(x1 - xcp) + (y1 - ycp)*(y1 - ycp) );// right top
      r1 = sqrt( (x1 - xcp)*(x1 - xcp) + (y1 - ycn)*(y1 - ycn) );// right bottom
      r2 = sqrt( (x1 - xcn)*(x1 - xcn) + (y1 - ycn)*(y1 - ycn) );// left bottom
      r3 = sqrt( (x1 - xcn)*(x1 - xcn) + (y1 - ycp)*(y1 - ycp) );// left top
      if( r0 < rp || r1 < rp || r2 < rp || r3 < rp ){
	xy25m_13mj0p1->Fill( px25j0, py13j0 );
	if( r0 < rp )
	  xy25m_13mj0p10->Fill( px25j0, py13j0 );
	if( r1 < rp )
	  xy25m_13mj0p11->Fill( px25j0, py13j0 );
	if( r2 < rp )
	  xy25m_13mj0p12->Fill( px25j0, py13j0 );
	if( r3 < rp )
	  xy25m_13mj0p13->Fill( px25j0, py13j0 );
      }else
	xy25m_13mj0c1->Fill( px25j0, py13j0 );
      
      x0 = px25j0;
      y0 = py13j0;
      r0 = sqrt( (x0 - xcp)*(x0 - xcp) + (y0 - ycp)*(y0 - ycp) );
      r1 = sqrt( (x0 - xcp)*(x0 - xcp) + (y0 - ycn)*(y0 - ycn) );
      r2 = sqrt( (x0 - xcn)*(x0 - xcn) + (y0 - ycn)*(y0 - ycn) );
      r3 = sqrt( (x0 - xcn)*(x0 - xcn) + (y0 - ycp)*(y0 - ycp) );
      if( r0 < rp || r1 < rp || r2 < rp || r3 < rp ){
	xy25m_13mj1p0->Fill( px25j1, py13j1 );
	if( r0 < rp )
	  xy25m_13mj1p00->Fill( px25j1, py13j1 );
	if( r1 < rp )
	  xy25m_13mj1p01->Fill( px25j1, py13j1 );
	if( r2 < rp )
	  xy25m_13mj1p02->Fill( px25j1, py13j1 );
	if( r3 < rp )
	  xy25m_13mj1p03->Fill( px25j1, py13j1 );
      }else
	xy25m_13mj1c0->Fill( px25j1, py13j1 );
    }

    Int_t j0p0 = 0;
    Int_t j1p0 = 0;
    if( goodj0yy == 1 && goodj0xx == 1 && goodj1yy == 1 && goodj1xx == 1 ){
      x0 = px25j0;
      y0 = py13j0;
      r0 = sqrt( (x0 - xcp)*(x0 - xcp) + (y0 - ycp)*(y0 - ycp) );
      if( r0 < rp ) j0p0 = 1;
      
      x1 = px25j1;
      y1 = py13j1;
      r0 = sqrt( (x1 - xcp)*(x1 - xcp) + (y1 - ycp)*(y1 - ycp) );// right top
      if( r0 < rp ) j1p0 = 1;
    }

    //if( ON_stg && j1p0 || !ON_stg && j0p0 ){
    /*if( j0p0 ){
      for( Int_t i = 0; i < Nch; i++ ){
	totj0p0[i]->Fill( tot_count[i][NstopL[i]-1] );
      }
    }*/
    // Jet 0
    if( LG_on ){
      if( goodj0yy == 1 && goodj0xx == 1 && goodj1yy == 1 && goodj1xx == 1 ){// j0 hit with j1 hit
	/*if( j0p0 ){
	  if( (512 - peak ) < PHcth[0] ){
	    xy25m_13mj0w1j0p0_cE[0]->Fill( px25j0, py13j0 );
	    xy25m_13mj1w0j0p0_cE[0]->Fill( px25j1, py13j1 );
	  }
	  if( (512 - peak ) >= PHcth[0] && (512 - peak ) < PHcth[1] ){
	    xy25m_13mj0w1j0p0_cE[1]->Fill( px25j0, py13j0 );
	    xy25m_13mj1w0j0p0_cE[1]->Fill( px25j1, py13j1 );
	  }
	  if( (512 - peak ) >= PHcth[1] && (512 - peak ) < PHcth[2] ){
	    xy25m_13mj0w1j0p0_cE[2]->Fill( px25j0, py13j0 );
	    xy25m_13mj1w0j0p0_cE[2]->Fill( px25j1, py13j1 );
	  }
	  if( (512 - peak ) >= PHcth[2] ){
	    xy25m_13mj0w1j0p0_cE[3]->Fill( px25j0, py13j0 );
	    xy25m_13mj1w0j0p0_cE[3]->Fill( px25j1, py13j1 );
	  }
	}*/
      }
    }

    // dtdc neighbering cells for 2D map
    if( goodj0y13 == 1 && goodj0x25 == 1 && goodj0y14 == 1 && goodj0x26 == 1 ){
      dx = px26j0 - px25j0;
      dy = py14j0 - py13j0;
      xy_dxdy_mj0->Fill( dx, dy );
    }
    if( goodj1y13 == 1 && goodj1x25 == 1 && goodj1y14 == 1 && goodj1x26 == 1 ){
      dx = px26j1 - px25j1;
      dy = py14j1 - py13j1;
      xy_dxdy_mj1->Fill( dx, dy );
    }// dtdc

    
    // tracking 25&13 to 25&13
    if( goodj0y13 == 1 && goodj0x25 == 1 && goodj1y13 == 1 && goodj1x25 == 1 ){
      x[0] = px25j0;
      x[1] = px25j1;
      y[0] = py13j0;
      y[1] = py13j1;
      
      if( ON_stg ){
	dxdy25m_13mj1_j0->Fill( x[1]-x[0], y[1]-y[0] );
	xdx25m_13mj1_j0->Fill( x[0], x[1]-x[0] );
	ydy25m_13mj1_j0->Fill( y[0], y[1]-y[0] );
	xdy25m_13mj1_j0->Fill( x[0], y[1]-y[0] );
	ydx25m_13mj1_j0->Fill( y[0], x[1]-x[0] );
      }else{
	dxdy25m_13mj1_j0->Fill( x[0]-x[1], y[0]-y[1] );
	xdx25m_13mj1_j0->Fill( x[1], x[0]-x[1] );
	ydy25m_13mj1_j0->Fill( y[1], y[0]-y[1] );
	xdy25m_13mj1_j0->Fill( x[1], y[0]-y[1] );
	ydx25m_13mj1_j0->Fill( y[1], x[0]-x[1] );
      }
      
      xy25m_13m_detz[0]->Fill( x[0], y[0] );
      xy25m_13m_detz[1]->Fill( x[1], y[1] );
      if( 512-peak < PHcth[0] ){
	xy25m_13m_detz_E0[0]->Fill( x[0], y[0] );
	xy25m_13m_detz_E0[1]->Fill( x[1], y[1] );
      }
      if( 512-peak < PHcth[1] && 512-peak >= PHcth[0] ){
	xy25m_13m_detz_E1[0]->Fill( x[0], y[0] );
	xy25m_13m_detz_E1[1]->Fill( x[1], y[1] );
      }
      if( 512-peak < PHcth[2] && 512-peak >= PHcth[1] ){
	xy25m_13m_detz_E2[0]->Fill( x[0], y[0] );
	xy25m_13m_detz_E2[1]->Fill( x[1], y[1] );
      }
      if( 512-peak >= PHcth[2] ){
	xy25m_13m_detz_E3[0]->Fill( x[0], y[0] );
	xy25m_13m_detz_E3[1]->Fill( x[1], y[1] );
      }
	
      linex = new TGraph( 2, z, x );
      liney = new TGraph( 2, z, y );

      linex->Fit("Line", "Q");
      fp0x = Line->GetParameter(0);
      fp1x = Line->GetParameter(1);
      
      liney->Fit("Line", "Q");
      fp0y = Line->GetParameter(0);
      fp1y = Line->GetParameter(1);
      
      for( Int_t i = 0; i < 7; i++ ){// interpolate
	px = pz[i]*fp1x + fp0x;
	py = pz[i]*fp1y + fp0y;
	xy25m_13m_z[i]->Fill( px, py );
      }
      for( Int_t i = 0; i < 6; i++ ){// finer interpolate
	px = fpz[i]*fp1x + fp0x;
	py = fpz[i]*fp1y + fp0y;
	xy25m_13m_fz[i]->Fill( px, py );
      }
      for( Int_t i = 0; i < Nexp; i++ ){// extrapolate
	px = expz[i]*fp1x + fp0x;
	py = expz[i]*fp1y + fp0y;
	if( 512-peak < PHcth[0] )
	  xy25m_13m_exz_E0[i]->Fill( px, py );
	if( 512-peak < PHcth[1] && 512-peak >= PHcth[0] )
	  xy25m_13m_exz_E1[i]->Fill( px, py );
	if( 512-peak < PHcth[2] && 512-peak >= PHcth[1] )
	  xy25m_13m_exz_E2[i]->Fill( px, py );
	if( 512-peak >= PHcth[2] )
	  xy25m_13m_exz_E3[i]->Fill( px, py );
	
	/*if( j0p0 ){
	  if( 512-peak < PHcth[0] )
	    xy25m_13m_j0p0_exz_E0[i]->Fill( px, py );
	  if( 512-peak < PHcth[1] && 512-peak >= PHcth[0] )
	    xy25m_13m_j0p0_exz_E1[i]->Fill( px, py );
	  if( 512-peak < PHcth[2] && 512-peak >= PHcth[1] )
	    xy25m_13m_j0p0_exz_E2[i]->Fill( px, py );
	  if( 512-peak >= PHcth[2] )
	    xy25m_13m_j0p0_exz_E3[i]->Fill( px, py );
	    }*/	  
      }// extrapolate
      for( Int_t i = 0; i < 2; i++ ){// extrapolate2
	px = expz2[i]*fp1x + fp0x;
	py = expz2[i]*fp1y + fp0y;
	if( 512-peak < PHcth[0] )
	  xy25m_13m_exz2_E0[i]->Fill( px, py );
	if( 512-peak < PHcth[1] && 512-peak >= PHcth[0] )
	  xy25m_13m_exz2_E1[i]->Fill( px, py );
	if( 512-peak < PHcth[2] && 512-peak >= PHcth[1] )
	  xy25m_13m_exz2_E2[i]->Fill( px, py );
	if( 512-peak >= PHcth[2] )
	  xy25m_13m_exz2_E3[i]->Fill( px, py );
      }
	
      if( LG_on ){
	x_lg = z_lg*fp1x + fp0x;
	y_lg = z_lg*fp1y + fp0y;
	xy_LG->Fill( x_lg, y_lg );
	x_LG->Fill( x_lg, 512-peak );
	y_LG->Fill( y_lg, 512-peak );
	x_intLGpos->Fill( x_lg, intpos );
	xw_intLGpos->Fill( x_lg, intpos );
	if( fabs(px25j0) < 40 && fabs(py13j0)<40 &&  fabs(px25j1) < 40 && fabs(py13j1)<40 )
	  LGff->Fill( 512-peak );
	
	if( 512-peak < PHcth[0] ){
	  xy_LG_E[0]->Fill( x_lg, y_lg );
	  if( fabs( px25j1 + 30 ) < 10 ) xy_LG_j1xneg_E[0]->Fill( x_lg, y_lg );
	}
	if( 512-peak >= PHcth[0] && 512-peak < PHcth[1] ){
	  xy_LG_E[1]->Fill( x_lg, y_lg );
	  if( fabs( px25j1 + 30 ) < 10 ) xy_LG_j1xneg_E[1]->Fill( x_lg, y_lg );
	}
	if( 512-peak >= PHcth[1] && 512-peak < PHcth[2] ){
	  xy_LG_E[2]->Fill( x_lg, y_lg );
	  if( fabs( px25j1 + 30 ) < 10 ) xy_LG_j1xneg_E[2]->Fill( x_lg, y_lg );
	}
	if( 512-peak >= PHcth[2] ){
	  xy_LG_E[3]->Fill( x_lg, y_lg );
	  if( fabs( px25j1 + 30 ) < 10 ) xy_LG_j1xneg_E[3]->Fill( x_lg, y_lg );
	}

	if( fabs( x_lg + 52.5 ) < 2.5 ){
	  if( 512-peak < PHcth[0] ){
	    xy25m_13m_detz_LGxneg_E0[0]->Fill( x[0], y[0] );
	    xy25m_13m_detz_LGxneg_E0[1]->Fill( x[1], y[1] );
	  }
	  if( 512-peak < PHcth[1] && 512-peak >= PHcth[0] ){
	    xy25m_13m_detz_LGxneg_E1[0]->Fill( x[0], y[0] );
	    xy25m_13m_detz_LGxneg_E1[1]->Fill( x[1], y[1] );
	  }
	  if( 512-peak < PHcth[2] && 512-peak >= PHcth[1] ){
	    xy25m_13m_detz_LGxneg_E2[0]->Fill( x[0], y[0] );
	    xy25m_13m_detz_LGxneg_E2[1]->Fill( x[1], y[1] );
	  }
	  if( 512-peak >= PHcth[2] ){
	    xy25m_13m_detz_LGxneg_E3[0]->Fill( x[0], y[0] );
	    xy25m_13m_detz_LGxneg_E3[1]->Fill( x[1], y[1] );
	  }
	  for( Int_t i = 0; i < Nexp; i++ ){// extrapolate
	    px = expz[i]*fp1x + fp0x;
	    py = expz[i]*fp1y + fp0y;
	    if( 512-peak < PHcth[0] )
	      xy25m_13m_exz_LGxneg_E0[i]->Fill( px, py );
	    if( 512-peak < PHcth[1] && 512-peak >= PHcth[0] )
	      xy25m_13m_exz_LGxneg_E1[i]->Fill( px, py );
	    if( 512-peak < PHcth[2] && 512-peak >= PHcth[1] )
	      xy25m_13m_exz_LGxneg_E2[i]->Fill( px, py );
	    if( 512-peak >= PHcth[2] )
	      xy25m_13m_exz_LGxneg_E3[i]->Fill( px, py );
	  }
	}
	/*if( j0p0 ){
	  if( 512-peak < PHcth[0] )
	    xy_LG_j0p0_E[0]->Fill( x_lg, y_lg );
	  if( 512-peak >= PHcth[0] && 512-peak < PHcth[1] )
	    xy_LG_j0p0_E[1]->Fill( x_lg, y_lg );
	  if( 512-peak >= PHcth[1] && 512-peak < PHcth[2] )
	    xy_LG_j0p0_E[2]->Fill( x_lg, y_lg );
	  if( 512-peak >= PHcth[2] )
	    xy_LG_j0p0_E[3]->Fill( x_lg, y_lg );
	}*/
	
	// QS on
	if( (runno0 == 1001 || runno0 == 1018 || runno0 == 1017 ) && fabs(x_lg - xc_on[0]) < xs_on[0] ) xmode_LG->Fill( 512 - peak );// 1 GeV/c
	if( (runno0 ==  999 || runno0 == 1020 )&& fabs(x_lg - xc_on[1]) < xs_on[1] ) xmode_LG->Fill( 512 - peak );// 2 GeV/c
	if( (runno0 == 1003 || runno0 == 1022 || runno0 > 1004 && runno0 < 1017 )&& fabs(x_lg - xc_on[2]) < xs_on[2] ) xmode_LG->Fill( 512 - peak );// 3 GeV/c
	if( (runno0 ==  991 || runno0 == 1024 ) && fabs(x_lg - xc_on[3]) < xs_on[3] ) xmode_LG->Fill( 512 - peak );// 4 GeV/c
	if( (runno0 == 1026) && fabs(x_lg - xc_on[4]) < xs_on[4] ) xmode_LG->Fill( 512 - peak );// 4.5 GeV/c
	
	// QS off
	if( (runno0 == 1002 || runno0 == 1019) && fabs(x_lg - xc_off[0]) < xs_off[0] ) xmode_LG->Fill( 512 - peak );// 1 GeV/c
	if( (runno0 == 1000 || runno0 == 1021) && fabs(x_lg - xc_off[1]) < xs_off[1] ) xmode_LG->Fill( 512 - peak );// 2 GeV/c
	if( (runno0 == 1004 || runno0 == 1023) && fabs(x_lg - xc_off[2]) < xs_off[2] ) xmode_LG->Fill( 512 - peak );// 3 GeV/c
	if( (runno0 ==  996 || runno0 == 1025) && fabs(x_lg - xc_off[3]) < xs_off[3] ) xmode_LG->Fill( 512 - peak );// 4 GeV/c
	if( (runno0 == 1027) && fabs(x_lg - xc_off[3]) < xs_off[4] ) xmode_LG->Fill( 512 - peak );// 4.5 GeV/c	  
      }
      
      ax_xj1_25_13->Fill( x[1], fp1x );
      ay_yj1_25_13->Fill( y[1], fp1y );
      
      Float_t rndm = gRandom->Uniform(0,1.);
      hrndm->Fill( rndm );
      if(  rndm < 0.001 ){
	c3->cd(Ndraw%9+1);
	TH1F *frame = gPad->DrawFrame( -50, -50, 50, 50 );
	path = new TArrow( x[1], y[1], x[0], y[0], 0.03, ">" );
	path->Draw();
	c4->cd(1);
	path->Draw("same");
	Ndraw++;
	if( (Ndraw%9) == 0 ){
	  //if( Ndraw == 9 ) c3->Print("orbit.pdf(" );
	  //else c3->Print("orbit.pdf" );
	  if( Ndraw == 9 ) pdf_orbit = "/group/itdc/tbl/2024ARTBL018/pdf/orbit/orbit_run" + runnumber + ".pdf(";
	  else pdf_orbit = "/group/itdc/tbl/2024ARTBL018/pdf/orbit/orbit_run" + runnumber + ".pdf";
	  c3->Print( pdf_orbit );
	  c3->Clear();
	  c3->Divide(3,3);
	  cout << "#Draw path " << Ndraw << " #processed " << Nproc << ", rate = "  << float(Ndraw)/float(Nproc) << endl;
	}
      }
    }else{// tracking
      xy25m_13m_undetz[0]->Fill( px25j0, py13j0 );
      xy25m_13m_undetz[1]->Fill( px25j1, py13j1 );
    }
    
    // tracking 26&14 to 26&14
    if( goodj0y14 == 1 && goodj0x26 && goodj1y14 == 1 && goodj1x26 == 1 ){
      x[0] = px26j0;
      x[1] = px26j1;
      y[0] = py14j0;
      y[1] = py14j1;
	
      linex = new TGraph( 2, z, x );
      liney = new TGraph( 2, z, y );
	
      linex->Fit("Line", "Q");
      fp0x = Line->GetParameter(0);
      fp1x = Line->GetParameter(1);
      
      liney->Fit("Line", "Q");
      fp0y = Line->GetParameter(0);
      fp1y = Line->GetParameter(1);
      
      ax_xj1_26_14->Fill( x[1], fp1x );
      ay_yj1_26_14->Fill( y[1], fp1y );
    }// tracking
    
    for( Int_t i = 0; i < NstopL[2]; i++ ){
      for( Int_t j = 0; j < NstopL[6]; j++ ){
	tdc_j0j1_x->Fill( tdc_Lcount[2][i], tdc_Lcount[6][j] );
      }
    }
    for( Int_t i = 0; i < NstopL[0]; i++ ){
      for( Int_t j = 0; j < NstopL[4]; j++ ){
	tdc_j0j1_y->Fill( tdc_Lcount[0][i], tdc_Lcount[4][j] );
      }
    }

    dtdc1313->Fill( tdc_Lcount[0][0]-tdc_Lcount[4][0] );
    dtdc2525->Fill( tdc_Lcount[2][0]-tdc_Lcount[6][0] );
    if(  abs( tdc_Lcount[0][0]-tdc_Lcount[4][0] + 100 ) < 200 ){
      tdc25m_13mj0cx->Fill( tdc_Lcount[2][0], tdc_Lcount[0][0] );
      tdc25m_13mj1cx->Fill( tdc_Lcount[6][0], tdc_Lcount[4][0] );
    }else{
      tdc25m_13mj0ncx->Fill( tdc_Lcount[2][0], tdc_Lcount[0][0] );
      tdc25m_13mj1ncx->Fill( tdc_Lcount[6][0], tdc_Lcount[4][0] );
    }
    if(  abs( tdc_Lcount[2][0]-tdc_Lcount[6][0] + 42 ) < 200 ){
      tdc25m_13mj0cy->Fill( tdc_Lcount[2][0], tdc_Lcount[0][0] );
      tdc25m_13mj1cy->Fill( tdc_Lcount[6][0], tdc_Lcount[4][0] );
    }else{
      tdc25m_13mj0ncy->Fill( tdc_Lcount[2][0], tdc_Lcount[0][0] );
      tdc25m_13mj1ncy->Fill( tdc_Lcount[6][0], tdc_Lcount[4][0] );
    }

    for( Int_t i = 0; i < Nch; i++ )
      for( Int_t j = 0; j < NstopL[i]; j++ )
      tdc_LCorr[i][j] = 0;

    for( Int_t i = 0; i < NstopL[0]-1; i++ ){
      for( Int_t j = i+1; j < NstopL[0]; j++ ){
	dtdc1313j0->Fill( tdc_Lcount[0][j]-tdc_Lcount[0][i] );
      }
    }
    for( Int_t i = 0; i < NstopL[2]-1; i++ ){
      for( Int_t j = i+1; j < NstopL[2]; j++ ){
	dtdc2525j0->Fill( tdc_Lcount[2][j]-tdc_Lcount[2][i] );
      }
    }
    for( Int_t i = 0; i < NstopL[4]-1; i++ ){
      for( Int_t j = i+1; j < NstopL[4]; j++ ){
	dtdc1313j1->Fill( tdc_Lcount[4][j]-tdc_Lcount[4][i] );
      }
    }
    for( Int_t i = 0; i < NstopL[6]-1; i++ ){
      for( Int_t j = i+1; j < NstopL[6]; j++ ){
	dtdc2525j1->Fill( tdc_Lcount[6][j]-tdc_Lcount[6][i] );
      }
    }
    
    for( Int_t i = 0; i < NstopL[0]; i++ ){
      for( Int_t j = 0; j < NstopL[1]; j++ ){
	dtdc1314j0->Fill( tdc_Lcount[0][i]-tdc_Lcount[1][j] );
	dtdc1314_tdc13j0->Fill( tdc_Lcount[0][i], tdc_Lcount[0][i]-tdc_Lcount[1][j] );
	if( abs (tdc_Lcount[0][i]-tdc_Lcount[1][j] - offty0 ) < n_tol*tolty0 ){
	  tdc_LCorr[0][i]++;
	  tdc_LCorr[1][j]++;
	}
      }
    }
    for( Int_t i = 0; i < NstopL[2]; i++ ){
      for( Int_t j = 0; j < NstopL[3]; j++ ){
	dtdc2526j0->Fill( tdc_Lcount[2][i]-tdc_Lcount[3][j] );
	dtdc2526_tdc25j0->Fill( tdc_Lcount[2][i], tdc_Lcount[2][i]-tdc_Lcount[3][j] );
	if( abs (tdc_Lcount[2][i]-tdc_Lcount[2][j] - offtx0 ) < n_tol*toltx0 ){
	  tdc_LCorr[2][i]++;
	  tdc_LCorr[3][j]++;
	}
      }
    }
    for( Int_t i = 0; i < NstopL[4]; i++ ){
      for( Int_t j = 0; j < NstopL[5]; j++ ){
	dtdc1314j1->Fill( tdc_Lcount[4][i]-tdc_Lcount[5][j] );
	dtdc1314_tdc13j1->Fill( tdc_Lcount[4][i], tdc_Lcount[4][i]-tdc_Lcount[5][j] );
	if( abs (tdc_Lcount[4][i]-tdc_Lcount[5][j] - offty1 ) < n_tol*tolty1 ){
	  tdc_LCorr[4][i]++;
	  tdc_LCorr[5][j]++;
	}
      }
    }
    for( Int_t i = 0; i < NstopL[6]; i++ ){
      for( Int_t j = 0; j < NstopL[7]; j++ ){
	dtdc2526j1->Fill( tdc_Lcount[6][i]-tdc_Lcount[7][j] );
	dtdc2526_tdc25j1->Fill( tdc_Lcount[6][i], tdc_Lcount[6][i]-tdc_Lcount[7][j] );
	if( abs (tdc_Lcount[6][i]-tdc_Lcount[7][j] - offtx1 ) < n_tol*toltx1 ){
	  tdc_LCorr[6][i]++;
	  tdc_LCorr[7][j]++;
	}
      }
    }

    for( Int_t i = 0; i < Nch; i++ ){
      for( Int_t j = 0; j < NstopL[i]; j++ ){
	nCorr[i]->Fill( tdc_LCorr[i][j] );
      }
    }
    
    if( NstopL[0] == 1 && NstopT[0] == 1 && NstopL[1] == 1 && NstopT[1] == 1 && tdc_Lcount[0][0] - tdc_Tcount[0][0] > 5 && tdc_Lcount[1][0] - tdc_Tcount[1][0] > 5  ){
      dtdc13m14mj0->Fill( tdc_Lcount[0][0]-tdc_Lcount[1][0] );
      dtdc13m14m_tdc13mj0->Fill( tdc_Lcount[0][0], tdc_Lcount[0][0]-tdc_Lcount[1][0] );
    }
    if( NstopL[2] == 1 && NstopT[2] == 1 && NstopL[3] == 1 && NstopT[3] == 1 && tdc_Lcount[2][0] - tdc_Tcount[2][0] > 5 && tdc_Lcount[3][0] - tdc_Tcount[3][0] > 5  ){
      dtdc25m26mj0->Fill( tdc_Lcount[2][0]-tdc_Lcount[3][0] );
      dtdc25m26m_tdc25mj0->Fill( tdc_Lcount[2][0], tdc_Lcount[2][0]-tdc_Lcount[3][0] );
    }
    if( NstopL[4] == 1 && NstopT[4] == 1 && NstopL[5] == 1 && NstopT[5] == 1 && tdc_Lcount[4][0] - tdc_Tcount[4][0] > 5 && tdc_Lcount[5][0] - tdc_Tcount[5][0] > 5  ){
      dtdc13m14mj1->Fill( tdc_Lcount[4][0]-tdc_Lcount[5][0] );
      dtdc13m14m_tdc13mj1->Fill( tdc_Lcount[4][0], tdc_Lcount[4][0]-tdc_Lcount[5][0] );
    }
    if( NstopL[6] == 1 && NstopT[6] == 1 && NstopL[7] == 1 && NstopT[7] == 1 && tdc_Lcount[6][0] - tdc_Tcount[6][0] > 5 && tdc_Lcount[7][0] - tdc_Tcount[7][0] > 5  ){
      dtdc25m26mj1->Fill( tdc_Lcount[6][0]-tdc_Lcount[7][0] );
      dtdc25m26m_tdc25mj1->Fill( tdc_Lcount[6][0], tdc_Lcount[6][0]-tdc_Lcount[7][0] );
    }

    for( Int_t i = 0; i < Nch; i++ ){
      NstopLC[i] = 0;
	for( Int_t j = 0; j < NstopL[i]; j++ ){
	  if( tdc_LCorr[i] ){
	    tdc_LCcount[i][NstopLC[i]] = tdc_Lcount[i][j];
	    if( NstopL[i] < NstopT[i] ){
	      tot_Ccount[i][NstopLC[i]] = tot_count[i][j];
	      tot_ihit_lc[i]->Fill( j, tot_Ccount[i][j] );
	    }
	    NstopLC[i]++;
	  }
	}
    }
    for( Int_t i = 0; i < Nch; i++ ){
      ntdcLC[i]->Fill( NstopLC[i] );
    }

    if( NstopLC[0] > 0 && NstopLC[2] > 0 ){
      tdc25mc_13mcj0->Fill( tdc_LCcount[2][0], tdc_LCcount[0][0] );
    }
    if( NstopLC[4] > 0 && NstopLC[6] > 0 ){
      tdc25mc_13mcj1->Fill( tdc_LCcount[6][0], tdc_LCcount[4][0] );
    }
    

    // TDC count 2D map (all hits)
    for( Int_t i = 0; i < NstopL[2]; i++ ){
      for( Int_t j = 0; j < NstopL[0]; j++ ){
	if( tdc_LCorr[2][i] && tdc_LCorr[0][j] )
	  tdc25c_13cj0->Fill( tdc_Lcount[2][i], tdc_Lcount[0][j] );
	else
	  tdc25nc_13ncj0->Fill( tdc_Lcount[2][i], tdc_Lcount[0][j] );
      }
    }
    for( Int_t i = 0; i < NstopL[6]; i++ ){
      for( Int_t j = 0; j < NstopL[4]; j++ ){
	if( tdc_LCorr[6][i] && tdc_LCorr[4][j] )
	  tdc25c_13cj1->Fill( tdc_Lcount[6][i], tdc_Lcount[4][j] );
	else
	  tdc25nc_13ncj1->Fill( tdc_Lcount[6][i], tdc_Lcount[4][j] );
      }
    }// 2D map

    for( Int_t i =0 ; i < Nch; i++ ){
      if( NstopL[i] < NstopT[i] ){
	for( Int_t j = 0; j < NstopL[i]; j++ ){
	  if( tdc_LCorr[i][j] )
	    tot_ihit_c[i]->Fill( j, tot_count[i][j] );
	  else
	    tot_ihit_nc[i]->Fill( j, tot_count[i][j] );
	}
      }
    }
    if( LG_on ) LG_tdc_j0x->Fill( tdc_Lcount[2][0], 512-peak );
    
    Nproc++;
  }// h   entry loop  (= end of data read)
  }// ii
  
  TText *t0, *t1;
  if( ON_stg ){
    t0 = new TText( 0, 0.5, "Up stream" );
    t1 = new TText( 0, 0.5, "Down stream" );
  }else{
    t0 = new TText( 0, 0.5, "Down stream" );
    t1 = new TText( 0, 0.5, "Up stream" );
  }
  
  //TCanvas *c1 = new TCanvas( "c1", "Canvas", 0, 0, 800, 800 );
  //c2->Print("waveform.pdf)" );
  pdf_waveform = "/group/itdc/tbl/2024ARTBL018/pdf/waveform/waveform_run" + runnumber + ".pdf)";
  c2->Print(pdf_waveform);

  //if( Ndraw%4 != 0 ) c3->Print("orbit.pdf)" );
  //pdf_orbit = "/group/itdc/tbl/2024ARTBL018/pdf/orbit/orbit_run" + runnumber + ".pdf)";
  if( Ndraw%4 != 0 ) c3->Print( pdf_orbit );

  pdf_orbit = "/group/itdc/tbl/2024ARTBL018/pdf/orbit/orbit_run" + runnumber + ".pdf)";
  c4->Print( pdf_orbit );

  c1->Clear();

  // #hit cells & random 
  c1->Divide( 1, 2 );
  c1->cd(1);
  LGvsEv_j0->Draw();
  c1->cd(2);
  LGvsEv_j1->Draw();

  pdf_tracking = "/group/itdc/tbl/2024ARTBL018/pdf/2D/tracking_run" + runnumber + ".pdf(";
  c1->Print(pdf_tracking);
  c1->Clear();

  // #hit cells & random 
  c1->Divide( 2, 2 );
  c1->cd(1);
  ncell[0]->Draw("box");
  c1->cd(2);
  ncell[1]->Draw("box");
  c1->cd(3);
  hrndm->Draw("box");
  c1->cd(4);
  stat->Draw();

  //pdf_tracking = "/group/itdc/tbl/2024ARTBL018/pdf/2D/tracking_run" + runnumber + ".pdf(";
  pdf_tracking = "/group/itdc/tbl/2024ARTBL018/pdf/2D/tracking_run" + runnumber + ".pdf";
  c1->Print(pdf_tracking);
  c1->Clear();

  // #hits per cell
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    gPad->SetLogy(1);
    ntdc[i]->Draw();
  }

  pdf_tracking = "/group/itdc/tbl/2024ARTBL018/pdf/2D/tracking_run" + runnumber + ".pdf";
  c1->Print(pdf_tracking);
  c1->Clear();

  // #hits (>0) per cell
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    gPad->SetLogy(1);
    ntdcn[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // #leading hits per cell
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    gPad->SetLogy(1);
    ntdcL[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // #trailing hits per cell
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    gPad->SetLogy(1);
    ntdcT[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // #leading hits correlate to other cell per cell
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    gPad->SetLogy(1);
    ntdcLC[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // #correlated hits cell per cell
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    gPad->SetLogy(1);
    nCorr[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tdc[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC of leading
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tdcL[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  //  TDC of trailing
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tdcT[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC of leading (one hit cell)
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tdcLs[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC of leading fit
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tdcL[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC of leading (1st hit)
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tdcLm[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // ToT
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tot[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // ToT 1st
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tot1st[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // ToT single hit
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    totsngl[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // ToT single hit
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    gPad->SetLogy(1);
    totsngl[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // ToT single hit
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    gPad->SetLogy(1);
    totsnglf[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // ToT j0p0
  /*c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    totj0p0[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();*/

  // ToT vs TDC
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tot_tdcL[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // ToT vs hit number
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tot_ihit[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // ToT vs TDC of 1st hit
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tot_tdcL_0[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // ToT vs TDC of n-th hits
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tot_tdcL_n[i]->Draw();
  }
  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC with no-Energy in LG
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tdcL_E0[i]->Draw();
  }
  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC with no-Energy in LG
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    gPad->SetLogy(1);
    tdcL_E0[i]->Draw();
  }
  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC with non0-Energy in LG
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tdcL_En0[i]->Draw();
  }
  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC with non0-Energy in LG
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    gPad->SetLogy(1);
    tdcL_En0[i]->Draw();
  }
  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC-2D, Jet0 TDC vs Jet1 TDC (x, y)
  c1->Divide( 2, 2 );
  c1->cd(1);
  tdc25_13j0->Draw();
  c1->cd(2);
  tdc25_13j1->Draw();
  c1->cd(3);
  tdc_j0j1_x->Draw();
  c1->cd(4);
  tdc_j0j1_y->Draw();
  
  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC-2D of 1st hits, Jet0 TDC 1st vs Jet1 TDC 1st (x, y)
  c1->Divide( 2, 2 );
  c1->cd(1);
  tdc25m_13mj0->Draw();
  c1->cd(2);
  tdc25m_13mj1->Draw();
  c1->cd(3);
  tdc_j0mj1m_x->Draw();
  c1->cd(4);
  tdc_j0mj1m_y->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC diff. Jet1-Jet0 (y,x), dy(ch14-ch13) vs dx(ch26-ch25) (Jet0, Jet1)
  c1->Divide( 2, 2 );
  c1->cd(1);
  //gPad->SetLogy(1);
  dtdc1313j0->Draw();
  c1->cd(2);
  dtdc2525j0->Draw();
  c1->cd(3);
  dtdc1313j1->Draw();
  c1->cd(4);
  dtdc2525j1->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC diff. Jet1-Jet0 (y,x), dy(ch14-ch13) vs dx(ch26-ch25) (Jet0, Jet1)
  c1->Divide( 2, 2 );
  c1->cd(1);
  //gPad->SetLogy(1);
  dtdc1313->Draw();
  c1->cd(2);
  dtdc2525->Draw();
  c1->cd(3);
  xy_dxdy_mj0->Draw();
  c1->cd(4);
  xy_dxdy_mj1->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC diff. ch14-ch13 Jet0, ch26-ch26 Jet0, ch14-ch13 Jet1, ch26-ch26 Jet1
  c1->Divide( 2, 2 );
  c1->cd(1);
  dtdc1314j0->Draw();
  c1->cd(2);
  dtdc2526j0->Draw();
  c1->cd(3);
  dtdc1314j1->Draw();
  c1->cd(4);
  dtdc2526j1->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC diff. ch14-ch13 Jet0, ch26-ch26 Jet0, ch14-ch13 Jet1, ch26-ch26 Jet1
  c1->Divide( 2, 2 );

  ctdcj0y->SetParameters( 20000, 0, 100, 1000, 0, 300, 300 );
  ctdcj0x->SetParameters( 20000, 0, 100, 1000, 0, 300, 300 );
  ctdcj1y->SetParameters( 20000, 0, 100, 1000, 0, 300, 300 );
  ctdcj1x->SetParameters( 20000, 0, 100, 1000, 0, 300, 300 );

  c1->cd(1);
  dtdc1314j0->Fit( "ctdcj0y", "ML" );
  c1->cd(2);
  dtdc2526j0->Fit( "ctdcj0x", "ML" );
  c1->cd(3);
  dtdc1314j1->Fit( "ctdcj1y", "ML" );
  c1->cd(4);
  dtdc2526j1->Fit( "ctdcj1x", "ML" );

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC diff. ch14-ch13 Jet0, ch26-ch26 Jet0, ch14-ch13 Jet1, ch26-ch26 Jet1
  c1->Divide( 2, 2 );
  c1->cd(1);
  dtdc1314_tdc13j0->Draw();
  c1->cd(2);
  dtdc2526_tdc25j0->Draw();
  c1->cd(3);
  dtdc1314_tdc13j1->Draw();
  c1->cd(4);
  dtdc2526_tdc25j1->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  Double_t dmean, stddev;
  // TDC diff. ch14-ch13 Jet0, ch26-ch26 Jet0, ch14-ch13 Jet1, ch26-ch26 Jet1
  c1->Divide( 2, 2 );
  c1->cd(1);
  dtdc13m14mj0->Draw();
  dmean = dtdc13m14mj0->GetMean();
  stddev = dtdc13m14mj0->GetStdDev();
  dfout << dmean << endl;
  dfout << stddev << endl;
  c1->cd(2);
  dtdc25m26mj0->Draw();
  dmean = dtdc25m26mj0->GetMean();
  stddev = dtdc25m26mj0->GetStdDev();
  dfout << dmean << endl;
  dfout << stddev << endl;
  c1->cd(3);
  dtdc13m14mj1->Draw();
  dmean = dtdc13m14mj1->GetMean();
  stddev = dtdc13m14mj1->GetStdDev();
  dfout << dmean << endl;
  dfout << stddev << endl;
  c1->cd(4);
  dtdc25m26mj1->Draw();
  dmean = dtdc25m26mj1->GetMean();
  stddev = dtdc25m26mj1->GetStdDev();
  dfout << dmean << endl;
  dfout << stddev << endl;

  c1->Print(pdf_tracking);
  c1->Clear();
  dfout.close();

  // TDC diff. ch14-ch13 Jet0, ch26-ch26 Jet0, ch14-ch13 Jet1, ch26-ch26 Jet1
  c1->Divide( 2, 2 );
  c1->cd(1);
  dtdc13m14m_tdc13mj0->Draw();
  c1->cd(2);
  dtdc25m26m_tdc25mj0->Draw();
  c1->cd(3);
  dtdc13m14m_tdc13mj1->Draw();
  c1->cd(4);
  dtdc25m26m_tdc25mj1->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC-2D, 13&14, 25&26Jet0 TDC vs Jet1 TDC (x, y)
  c1->Divide( 2, 2 );
  c1->cd(1);
  tdc25c_13cj0->Draw();
  c1->cd(2);
  tdc25c_13cj1->Draw();
  c1->cd(3);
  tdc25nc_13ncj0->Draw();
  c1->cd(4);
  tdc25nc_13ncj1->Draw();
  
  c1->Print(pdf_tracking);
  c1->Clear();

  // ToT vs hit number
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tot_ihit_c[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // ToT vs hit number
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tot_ihit_nc[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // ToT vs hit number
  c1->Divide( 2, 4 );
  for( Int_t i = 0; i < Nch; i++ ){
    c1->cd(i+1);
    tot_ihit_lc[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC-2D of 1st hits, Jet0 TDC 1st vs Jet1 TDC 1st (x, y)
  c1->Divide( 2, 2 );
  c1->cd(1);
  tdc25mc_13mcj0->Draw();
  c1->cd(2);
  tdc25mc_13mcj1->Draw();
  
  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC-2D dTDC(x) of Jet0 and Jet1 is small or large
  c1->Divide( 2, 2 );
  c1->cd(1);
  tdc25m_13mj0cx->Draw();
  c1->cd(2);
  tdc25m_13mj0ncx->Draw();
  c1->cd(3);
  tdc25m_13mj1cx->Draw();
  c1->cd(4);
  tdc25m_13mj1ncx->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC-2D dTDC(y) of Jet0 and Jet1 is small or large
  c1->Divide( 2, 2 );
  c1->cd(1);
  tdc25m_13mj0cy->Draw();
  c1->cd(2);
  tdc25m_13mj0ncy->Draw();
  c1->cd(3);
  tdc25m_13mj1cy->Draw();
  c1->cd(4);
  tdc25m_13mj1ncy->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // TDC-2D Jet0(no Jet1 hit/with Jet1 hit) vice versa
  c1->Divide( 2, 2 );
  c1->cd(1);
  tdc25m_13mj0n1->Draw();
  c1->cd(2);
  tdc25m_13mj0w1->Draw();
  c1->cd(3);
  tdc25m_13mj1n0->Draw();
  c1->cd(4);
  tdc25m_13mj1w0->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D-map of Jet0, all, with Jet1 hit or no Jet1 hit
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy25m_13mj0->Draw();
  c1->cd(3);
  xy25m_13mj0n1->Draw();
  c1->cd(4);
  xy25m_13mj0w1->Draw();
  c1->cd(2);
  t0->Draw();
  
  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D-map of Jet1, all, with Jet0 hit or no Jet0 hit
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy25m_13mj1->Draw();
  c1->cd(3);
  xy25m_13mj1n0->Draw();
  c1->cd(4);
  xy25m_13mj1w0->Draw();
  c1->cd(2);
  t1->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map and projection Jet0 (first hit, no ather cut)
  c1->Divide( 2, 2 );
  c1->cd(2);
  xy25m_13mj0->Draw();
  c1->cd(4);
  TH1D *xy25m_13mj0_prox = xy25m_13mj0->ProjectionX( "xy25m_13mj0_prox", 1, 100 );
  xy25m_13mj0_prox->Draw();
  c1->cd(1);
  TH1D *xy25m_13mj0_proy = xy25m_13mj0->ProjectionY( "xy25m_13mj0_proy", 1, 100 );
  xy25m_13mj0_proy->Draw();
  c1->cd(3);
  t0->Draw();
    
  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map and projection Jet1 (first hit, no ather cut)
  c1->Divide( 2, 2 );
  c1->cd(2);
  xy25m_13mj1->Draw();
  c1->cd(4);
  TH1D *xy25m_13mj1_prox = xy25m_13mj1->ProjectionX( "xy25m_13mj1_prox", 1, 100 );
  xy25m_13mj1_prox->Draw();
  c1->cd(1);
  TH1D *xy25m_13mj1_proy = xy25m_13mj1->ProjectionY( "xy25m_13mj1_proy", 1, 100 );
  xy25m_13mj1_proy->Draw();
  c1->cd(3);
  t1->Draw();
    
  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map and projection Jet0 (first hit, no hit in Jet1 )
  c1->Divide( 2, 2 );
  c1->cd(2);
  xy25m_13mj0n1->Draw();
  c1->cd(4);
  TH1D *xy25m_13mj0n1_prox = xy25m_13mj0n1->ProjectionX( "xy25m_13mj0n1_prox", 1, 100 );
  xy25m_13mj0n1_prox->Draw();
  c1->cd(1);
  TH1D *xy25m_13mj0n1_proy = xy25m_13mj0n1->ProjectionY( "xy25m_13mj0n1_proy", 1, 100 );
  xy25m_13mj0n1_proy->Draw();
  c1->cd(3);
  t0->Draw();
    
  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map and projection Jet0 (first hit, wit hit in Jet1 )
  c1->Divide( 2, 2 );
  c1->cd(2);
  xy25m_13mj0w1->Draw();
  c1->cd(4);
  TH1D *xy25m_13mj0w1_prox = xy25m_13mj0w1->ProjectionX( "xy25m_13mj0w1_prox", 1, 100 );
  xy25m_13mj0w1_prox->Draw();
  c1->cd(1);
  TH1D *xy25m_13mj0w1_proy = xy25m_13mj0w1->ProjectionY( "xy25m_13mj0w1_proy", 1, 100 );
  xy25m_13mj0w1_proy->Draw();
  c1->cd(3);
  t0->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map and projection Jet1 (first hit, no hit in Jet0 )
  c1->Divide( 2, 2 );
  c1->cd(2);
  xy25m_13mj1n0->Draw();
  c1->cd(4);
  TH1D *xy25m_13mj1n0_prox = xy25m_13mj1n0->ProjectionX( "xy25m_13mj1n0_prox", 1, 100 );
  xy25m_13mj1n0_prox->Draw();
  c1->cd(1);
  TH1D *xy25m_13mj1n0_proy = xy25m_13mj1n0->ProjectionY( "xy25m_13mj1n0_proy", 1, 100 );
  xy25m_13mj1n0_proy->Draw();
  c1->cd(3);
  t1->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map and projection Jet1 (first hit, wit hit in Jet0 )
  c1->Divide( 2, 2 );
  c1->cd(2);
  xy25m_13mj1w0->Draw();
  c1->cd(4);
  TH1D *xy25m_13mj1w0_prox = xy25m_13mj1w0->ProjectionX( "xy25m_13mj1w0_prox", 1, 100 );
  xy25m_13mj1w0_prox->Draw();
  c1->cd(1);
  TH1D *xy25m_13mj1w0_proy = xy25m_13mj1w0->ProjectionY( "xy25m_13mj1w0_proy", 1, 100 );
  xy25m_13mj1w0_proy->Draw();
  c1->cd(3);
  t1->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D difference map
  c1->Divide( 2, 2 );
  c1->cd(2);
  dxdy25m_13mj1_j0->Draw();
  c1->cd(4);
  xdx25m_13mj1_j0->Draw();
  c1->cd(1);
  ydy25m_13mj1_j0->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D difference map
  c1->Divide( 2, 2 );
  c1->cd(1);
  TH1D *projdx_xdx = xdx25m_13mj1_j0->ProjectionY( "xdx_projdx", 1, 100 );
  projdx_xdx->Draw();

  c1->cd(2);
  TH1D *projdx_ydx = ydx25m_13mj1_j0->ProjectionY( "ydx_projdx", 1, 100 );
  projdx_ydx->Draw();

  c1->cd(3);
  TH1D *slidx_ydx = ydx25m_13mj1_j0->ProjectionY( "ydx_slidx", 41, 60 );
  slidx_ydx->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D difference map
  c1->Divide( 2, 2 );
  c1->cd(1);
  xdx25m_13mj1_j0->Draw();
  c1->cd(2);
  ydy25m_13mj1_j0->Draw();
  c1->cd(3);
  xdy25m_13mj1_j0->Draw();
  c1->cd(4);
  ydx25m_13mj1_j0->Draw();
  
  c1->Print(pdf_tracking);
  c1->Clear();
  
  // 2D-map outside of time window
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy25o_13oj0->Draw();
  c1->cd(2);
  xy26o_14oj0->Draw();
  c1->cd(3);
  xy25o_13oj1->Draw();
  c1->cd(4);
  xy26o_14oj1->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();
  
  // 2D-map outside of time window
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy25o_13oj0w1->Draw();
  c1->cd(2);
  xy26o_14oj0w1->Draw();
  c1->cd(3);
  xy25o_13oj1w0->Draw();
  c1->cd(4);
  xy26o_14oj1w0->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();
  
  // 2D-map outside of time window
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy25o_13oj0n1->Draw();
  c1->cd(2);
  xy26o_14oj0n1->Draw();
  c1->cd(3);
  xy25o_13oj1n0->Draw();
  c1->cd(4);
  xy26o_14oj1n0->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();
  
  // 2D map Jet0 with Jet1 hit, 2D map Jet1 with Jet0 hit
  // slope of 2hits connecting line 25,13
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy25m_13mj0w1->Draw();
  c1->cd(2);
  xy25m_13mj1w0->Draw();
  c1->cd(3);
  ax_xj1_25_13->Draw();
  c1->cd(4);
  ay_yj1_25_13->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map Jet0 with Jet1 hit, 2D map Jet1 with Jet0 hit
  // slope of 2hits connecting line 26,14
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy26m_14mj0w1->Draw();
  c1->cd(2);
  xy26m_14mj1w0->Draw();
  c1->cd(3);
  ax_xj1_26_14->Draw();
  c1->cd(4);
  ay_yj1_26_14->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map Jet0 with Jet1 hit, 2D map Jet1 with Jet0 hit
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy25m_13mj0w1->Draw();
  for( Int_t i = 0; i < Nq; i ++ )
    Q[i]->Draw("only");
  c1->cd(2);
  xy25m_13mj1w0->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map Jet0 with Jet1, Jet1 hit is hide of the pole or center
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy25m_13mj0w1_25_13->Draw();
  c1->cd(3);
  xy25m_13mj0p1->Draw();
  c1->cd(4);
  xy25m_13mj0c1->Draw();
  c1->cd(2);
  t0->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map Jet0 with Jet1, Jet1 hit is hide of the each pole
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy25m_13mj0p13->Draw();
  c1->cd(2);
  xy25m_13mj0p10->Draw();
  c1->cd(3);
  xy25m_13mj0p12->Draw();
  c1->cd(4);
  xy25m_13mj0p11->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map Jet1 with Jet1, Jet0 hit is hide of the pole or center
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy25m_13mj1w0_25_13->Draw();
  c1->cd(3);
  xy25m_13mj1p0->Draw();
  c1->cd(4);
  xy25m_13mj1c0->Draw();
  c1->cd(2);
  t1->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map Jet1 with Jet0, Jet1 hit is hide of the each pole
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy25m_13mj1p03->Draw();
  c1->cd(2);
  xy25m_13mj1p00->Draw();
  c1->cd(3);
  xy25m_13mj1p02->Draw();
  c1->cd(4);
  xy25m_13mj1p01->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D map 
  c1->Divide( 2, 2 );
  c1->cd(1);
  xy25m_13m_detz[0]->Draw();
  c1->cd(2);
  xy25m_13m_undetz[0]->Draw();
  c1->cd(3);
  xy25m_13m_detz[1]->Draw();
  c1->cd(4);
  xy25m_13m_undetz[1]->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // 2D-map of Jet0 -interpolation- Jet1
  c1->Divide( 3, 3 );
  c1->cd(1);
  if( ON_stg )
    //xy25m_13mj0w1->Draw();
    xy25m_13m_detz[0]->Draw();
  else
    //xy25m_13mj1w0->Draw();
    xy25m_13m_detz[1]->Draw();
  for( Int_t i = 0; i < 7; i++ ){
    c1->cd(i+2);
    xy25m_13m_z[i]->Draw();
  }
  c1->cd(9);
  if( ON_stg )
    //xy25m_13mj1w0->Draw();
    xy25m_13m_detz[1]->Draw();
  else
    //xy25m_13mj0w1->Draw();
    xy25m_13m_detz[0]->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  // projection of 2D-map on x-axis Jet0 -interpolation- Jet1
  c1->Divide( 3, 3 );
  TH1D *twoD_prox[9];
  TH1D *twoD_proy[9];
  TString xname;
  TString yname;
  Float_t mean;
  Float_t sgm;

  Float_t xcen[9];
  Float_t ycen[9];
  Float_t excen[9];
  Float_t eycen[9];
  Float_t xwid_r[9];
  Float_t xwid_l[9];
  Float_t ywid_u[9];
  Float_t ywid_l[9];
  Float_t exwid_r[9];
  Float_t exwid_l[9];
  Float_t eywid_u[9];
  Float_t eywid_l[9];
  Float_t xc_max;
  Float_t yc_max;
  Float_t xc_min;
  Float_t yc_min;
  Float_t xw_max;
  Float_t yw_max;
  Float_t xw_min;
  Float_t yw_min;

  Float_t xfrom, xto;
  Float_t yfrom, yto;
  Int_t binfrom, binto;
  Float_t total;
  Float_t xwcen[9];
  Float_t ywcen[9];
  Float_t xwwid[9];
  Float_t ywwid[9];

  Float_t ave;
  Float_t avv;
  
  c1->cd(1);
  fit2D->SetParameters( 1000, 0, 0, 10, 10, 0, 1 );
  xname = "25m_13m_z0_prox";
  yname = "25m_13m_z0_proy";
  if( ON_stg ){
    //twoD_prox[0] = xy25m_13mj0w1->ProjectionX( xname, 1, 100 );
    //twoD_proy[0] = xy25m_13mj0w1->ProjectionY( yname, 1, 100 );
    twoD_prox[0] = xy25m_13m_detz[0]->ProjectionX( xname, 1, 100 );
    twoD_proy[0] = xy25m_13m_detz[0]->ProjectionY( yname, 1, 100 );
  }else{
    //twoD_prox[0] = xy25m_13mj1w0->ProjectionX( xname, 1, 100 );
    //twoD_proy[0] = xy25m_13mj1w0->ProjectionY( yname, 1, 100 );
    twoD_prox[0] = xy25m_13m_detz[1]->ProjectionX( xname, 1, 100 );
    twoD_proy[0] = xy25m_13m_detz[1]->ProjectionY( yname, 1, 100 );
  }
  //mean = twoD_prox[0]->GetMean();
  //sgm = twoD_prox[0]->GetStdDev();
  //fitpro->SetParameters( 1000, mean, sgm, sgm, 1, -40, 40 );
  //fitpro->SetParameters( 1000, 0, 40, 40, 1, -40, 40 );
  //twoD_prox[0]->Fit( "fitpro" );
  //xcen[0] = fitpro->GetParameter(1);
  //xwid_l[0] = fitpro->GetParameter(2);
  //xwid_r[0] = fitpro->GetParameter(3);
  //excen[0] = fitpro->GetParError(1);
  //exwid_l[0] = fitpro->GetParError(2);
  //exwid_r[0] = fitpro->GetParError(3);
  twoD_prox[0]->Draw();
  xcen[0] = twoD_prox[0]->GetMean();
  xwid_l[0] = twoD_prox[0]->GetStdDev();
  excen[0] = twoD_prox[0]->GetMeanError();
  exwid_l[0] = twoD_prox[0]->GetStdDevError();
  xc_max = xcen[0];
  xc_min = xcen[0];
  xw_max = xwid_l[0];
  xw_min = xwid_l[0];
  // calc width at rate center
  xfrom = xcen[0] - 20;// mm
  xto = xcen[0] + 20;// mm
  binfrom = xfrom + 51;// -50 mm - 50 mm/100 bins = 1 mm/bin. 1st bin = -50 + 51
  binto = xto + 51;
  ave = 0;
  avv = 0;
  total = 0;
  for( Int_t j = binfrom; j < binto; j++ ){
    ave = ave + (twoD_prox[0]->GetBinContent( j ))*( j-51 );
    avv = avv + (twoD_prox[0]->GetBinContent( j ))*( j-51 )*(j - 51);
    total = total + (twoD_prox[0]->GetBinContent( j ));
  }
  ave = ave/total;
  avv = avv/total;
  avv = sqrt(avv-ave*ave);
  xwcen[0] = ave;
  xwwid[0] = avv;
  //
  for( Int_t i = 0; i < 7; i++ ){
    c1->cd(i+2);
    xname = "25m_13m_z" + to_string( i+1 ) + "_prox";
    yname = "25m_13m_z" + to_string( i+1 ) + "_proy";
    twoD_prox[i+1] = xy25m_13m_z[i]->ProjectionX( xname, 1, 100 );
    twoD_proy[i+1] = xy25m_13m_z[i]->ProjectionY( yname, 1, 100 );
    //mean = twoD_prox[i+1]->GetMean();
    //sgm = twoD_prox[i+1]->GetStdDev();
    //fitpro->SetParameters( 1000, mean, sgm, sgm, 1, -40, 40 );
    //fitpro->SetParameters( 1000, 0, 40, 40, 1, -40, 40 );
    //fitpro->FixParameter( 5, -40 );
    //fitpro->FixParameter( 6, 40 );
    //twoD_prox[i+1]->Fit( "fitpro", "ME" );
    //xcen[i+1] = fitpro->GetParameter(1);
    //xwid_l[i+1] = fitpro->GetParameter(2);
    //xwid_r[i+1] = fitpro->GetParameter(3);
    //excen[i+1] = fitpro->GetParError(1);
    //exwid_l[i+1] = fitpro->GetParError(2);
    //exwid_r[i+1] = fitpro->GetParError(3);
    twoD_prox[i+1]->Draw();
    xcen[i+1] = twoD_prox[i+1]->GetMean();
    xwid_l[i+1] = twoD_prox[i+1]->GetStdDev();
    excen[i+1] = twoD_prox[i+1]->GetMeanError();
    exwid_l[i+1] = twoD_prox[i+1]->GetStdDevError();
    if( xc_max < xcen[i+1] ) xc_max = xcen[i+1];
    if( xc_min > xcen[i+1] ) xc_min = xcen[i+1];
    if( xw_max < xwid_l[i+1] ) xw_max = xwid_l[i+1];
    if( xw_min > xwid_l[i+1] ) xw_min = xwid_l[i+1];
    // calc width at rate center
    xfrom = xcen[i+1] - 20;// mm
    xto = xcen[i+1] + 20;// mm
    binfrom = xfrom + 51;// -50 mm - 50 mm/100 bins = 1 mm/bin. 1st bin = -50 + 51
    binto = xto + 51;
    ave = 0;
    avv = 0;
    total = 0;
    for( Int_t j = binfrom; j < binto; j++ ){
      ave = ave + (twoD_prox[i+1]->GetBinContent( j ))*( j-51 );
      avv = avv + (twoD_prox[i+1]->GetBinContent( j ))*( j-51 )*(j - 51);
      total = total + (twoD_prox[i+1]->GetBinContent( j ));
    }
    ave = ave/total;
    avv = avv/total;
    avv = sqrt(avv-ave*ave);
    xwcen[i+1] = ave;
    xwwid[i+1] = avv;
    //
  }
  c1->cd(9);
  //fit2D->SetParameters( 1000, 0, 0, 10, 10, 0, 1 );
  xname = "25m_13m_z8_prox";
  yname = "25m_13m_z8_proy";
  if( ON_stg ){
    //twoD_prox[8] = xy25m_13mj1w0->ProjectionX( xname, 1, 100 );
    //twoD_proy[8] = xy25m_13mj1w0->ProjectionY( yname, 1, 100 );
    twoD_prox[8] = xy25m_13m_detz[1]->ProjectionX( xname, 1, 100 );
    twoD_proy[8] = xy25m_13m_detz[1]->ProjectionY( yname, 1, 100 );
  }else{
    //twoD_prox[8] = xy25m_13mj0w1->ProjectionX( xname, 1, 100 );
    //twoD_proy[8] = xy25m_13mj0w1->ProjectionY( yname, 1, 100 );
    twoD_prox[8] = xy25m_13m_detz[0]->ProjectionX( xname, 1, 100 );
    twoD_proy[8] = xy25m_13m_detz[0]->ProjectionY( yname, 1, 100 );
  }
  twoD_prox[8]->Draw();
  xcen[8] = twoD_prox[8]->GetMean();
  xwid_l[8] = twoD_prox[8]->GetStdDev();
  excen[8] = twoD_prox[8]->GetMeanError();
  exwid_l[8] = twoD_prox[8]->GetStdDevError();
  if( xc_max < xcen[8] ) xc_max = xcen[8];
  if( xc_min > xcen[8] ) xc_min = xcen[8];
  if( xw_max < xwid_l[8] ) xw_max = xwid_l[8];
  if( xw_min > xwid_l[8] ) xw_min = xwid_l[8];
  // calc width at rate center
  xfrom = xcen[8] - 20;// mm
  xto = xcen[8] + 20;// mm
  binfrom = xfrom + 51;// -50 mm - 50 mm/100 bins = 1 mm/bin. 1st bin = -50 + 51
  binto = xto + 51;
  ave = 0;
  avv = 0;
  total = 0;
  for( Int_t j = binfrom; j < binto; j++ ){
    ave = ave + (twoD_prox[8]->GetBinContent( j ))*( j-51 );
    avv = avv + (twoD_prox[8]->GetBinContent( j ))*( j-51 )*(j - 51);
    total = total + (twoD_prox[8]->GetBinContent( j ));
  }
  ave = ave/total;
  avv = avv/total;
  avv = sqrt(avv-ave*ave);
  xwcen[8] = ave;
  xwwid[8] = avv;
  //

  c1->Print(pdf_tracking);
  c1->Clear();

  // projection of 2D-map on y-axis Jet0 -interpolation- Jet1
  c1->Divide( 3, 3 );
  yc_max = ycen[0];
  yc_min = ycen[0];
  yw_max = ywid_l[0];
  yw_min = ywid_l[0];
  for( Int_t i = 0; i < 9; i++ ){
    c1->cd(i+1);
    twoD_proy[i]->Draw();
    ycen[i] = twoD_proy[i]->GetMean();
    ywid_l[i] = twoD_proy[i]->GetStdDev();
    eycen[i] = twoD_proy[i]->GetMeanError();
    eywid_l[i] = twoD_proy[i]->GetStdDevError();
    //cout << "ywid " << ywid_l[i] << endl;
    if( i == 0 ){
      yc_max = ycen[i];
      yc_min = ycen[i];
      yw_max = ywid_l[i];
      yw_min = ywid_l[i];
    }else{
      if( yc_max < ycen[i] ) yc_max = ycen[i];
      if( yc_min > ycen[i] ) yc_min = ycen[i];
      if( yw_max < ywid_l[i] ) yw_max = ywid_l[i];
      if( yw_min > ywid_l[i] ) yw_min = ywid_l[i];
    }
    // calc width at rate center
    yfrom = ycen[i] - 20;// mm
    yto = ycen[i] + 20;// mm
    binfrom = yfrom + 51;// -50 mm - 50 mm/100 bins = 1 mm/bin. 1st bin = -50 + 51
    binto = yto + 51;
    ave = 0;
    avv = 0;
    total = 0;
    for( Int_t j = binfrom; j < binto; j++ ){
      ave = ave + (twoD_proy[i]->GetBinContent( j ))*( j-51 );
      avv = avv + (twoD_proy[i]->GetBinContent( j ))*( j-51 )*(j - 51);
      total = total + (twoD_proy[i]->GetBinContent( j ));
    }
    ave = ave/total;
    avv = avv/total;
    avv = sqrt(avv-ave*ave);
    ywcen[i] = ave;
    ywwid[i] = avv;
    //
}

  c1->Print(pdf_tracking);
  c1->Clear();

  Float_t wz[9];
  Float_t ewz[9] = { 20, 20, 20, 20, 20, 20, 20, 20, 20 };
  if( ON_stg )
    wz[0] = z[0];
  else
    wz[0] = z[1];
  for( Int_t i = 0; i < 7; i++ ) wz[i+1] = pz[i];
  if( ON_stg )
    wz[8] = z[1];
  else
    wz[8] = z[0];

  // width of projection on x-axis
  c1->Divide( 1, 2 );
  c1->cd(1);
  TGraphErrors *xwaist = new TGraphErrors( 9, wz, xwid_l, ewz, exwid_l );
  Float_t xw_wid = xw_max - xw_min;
  Float_t z_wid = fabs( z[0] - z[1] );
  TH1F *frameBxw = gPad->DrawFrame( wz[0]-z_wid*0.1, xw_min-xw_wid*0.5, wz[8]+z_wid*0.1, xw_max+xw_wid*0.5 );
  frameBxw->GetXaxis()->SetTitle( "z (mm)" );
  frameBxw->GetYaxis()->SetTitle( "StdDev(x) (mm)" );
  xwaist->Draw("same");
  
  c1->cd(2);
  TGraphErrors *ywaist = new TGraphErrors( 9, wz, ywid_l, ewz, eywid_l );
  Float_t yw_wid = yw_max - yw_min;
  TH1F *frameByw = gPad->DrawFrame( wz[0]-z_wid*0.1, yw_min-yw_wid*0.5, wz[8]+z_wid*0.1, yw_max+yw_wid*0.5 );
  frameByw->GetXaxis()->SetTitle( "z (mm)" );
  frameByw->GetYaxis()->SetTitle( "StdDev(y) (mm)" );
  ywaist->Draw("same");
  
  c1->Print(pdf_tracking);
  c1->Clear();

  // width of projection on x-axis fixed frame
  c1->Divide( 1, 2 );
  c1->cd(1);
  TH1F *frameBxwf = gPad->DrawFrame( wz[0]-z_wid*0.1, 0, wz[8]+z_wid*0.1, 40 );
  frameBxwf->GetXaxis()->SetTitle( "z (mm)" );
  frameBxwf->GetYaxis()->SetTitle( "StdDev(x) (mm)" );
  xwaist->Draw("same");
  
  c1->cd(2);
  TH1F *frameBywf = gPad->DrawFrame( wz[0]-z_wid*0.1, 0, wz[8]+z_wid*0.1, 25 );
  frameBywf->GetXaxis()->SetTitle( "z (mm)" );
  frameBywf->GetYaxis()->SetTitle( "StdDev(y) (mm)" );
  ywaist->Draw("same");
  
  c1->Print(pdf_tracking);
  c1->Clear();

  ofstream zout( zname, ios::out );
  // width of projection on x-axis
  c1->Divide( 1, 2 );
  c1->cd(1);
  para->SetParameter(2, xw_min );
  if( ON_stg ){
    para->SetParameter(0, (xw_max-xw_min)/3000/3000 );
    para->SetParameter(1, 1000 );
  }else{
    para->SetParameter(0, (xw_max-xw_min)/4000/4000 );
    para->SetParameter(1, -1000 );
  }
  xwaist->Draw("AP");
  xwaist->Fit("para");
  zout << para->GetParameter( 1 ) << endl;
  zout << para->GetParError( 1 ) << endl;
    
  c1->cd(2);
  para->SetParameter(2, yw_min );
  if( ON_stg ){
    para->SetParameter(0, (yw_max-yw_min)/3000/3000 );
    para->SetParameter(1, 1000 );
  }else{
    para->SetParameter(0, (xw_max-xw_min)/4000/4000 );
    para->SetParameter(1, -1000 );
  }
  ywaist->Draw("AP");
  ywaist->Fit("para");
  zout << para->GetParameter( 1 ) << endl;
  zout << para->GetParError( 1 ) << endl;
  zout.close();

  c1->Print(pdf_tracking);
  c1->Clear();

  // width of projection on x-axis fixed frame
  c1->Divide( 1, 2 );
  c1->cd(1);
  TH1F *frameBxwf0 = gPad->DrawFrame( wz[0]-z_wid*0.1, 0, wz[8]+z_wid*0.1, 40 );
  frameBxwf0->GetXaxis()->SetTitle( "z (mm)" );
  frameBxwf0->GetYaxis()->SetTitle( "StdDev(x) (mm)" );
  xwaist->Draw("Psame");
  
  c1->cd(2);
  TH1F *frameBywf0 = gPad->DrawFrame( wz[0]-z_wid*0.1, 0, wz[8]+z_wid*0.1, 25 );
  frameBywf0->GetXaxis()->SetTitle( "z (mm)" );
  frameBywf0->GetYaxis()->SetTitle( "StdDev(y) (mm)" );
  ywaist->Draw("Psame");
  
  c1->Print(pdf_tracking);
  c1->Clear();

  // mean of projection on x-axis
  c1->Divide( 1, 2 );
  c1->cd(1);
  TGraphErrors *xcent = new TGraphErrors( 9, wz, xcen, ewz, excen );
  Float_t xc_wid = xc_max - xc_min;
  TH1F *frameBxc = gPad->DrawFrame( wz[0]-z_wid*0.1, xc_min-xc_wid*0.5, wz[8]+z_wid*0.1, xc_max+xc_wid*0.5 );
  frameBxc->GetXaxis()->SetTitle( "z (mm)" );
  frameBxc->GetYaxis()->SetTitle( "Mean x (mm)" );
  xcent->Draw("same");
  
  c1->cd(2);
  TGraphErrors *ycent = new TGraphErrors( 9, wz, ycen, ewz, eycen );
  Float_t yc_wid = yc_max - yc_min;
  TH1F *frameByc = gPad->DrawFrame( wz[0]-z_wid*0.1, yc_min-yc_wid*0.5, wz[8]+z_wid*0.1, yc_max+yc_wid*0.5 );
  frameByc->GetXaxis()->SetTitle( "z (mm)" );
  frameByc->GetYaxis()->SetTitle( "Mean y (mm)" );
  ycent->Draw("same");
  
  c1->Print(pdf_tracking);
  c1->Clear();

  // mean of projection on x-axis fixed frame
  c1->Divide( 1, 2 );
  c1->cd(1);
  TH1F *frameBxcf = gPad->DrawFrame( wz[0]-z_wid*0.1, -12, wz[8]+z_wid*0.1, 6 );
  frameBxcf->GetXaxis()->SetTitle( "z (mm)" );
  frameBxcf->GetYaxis()->SetTitle( "Mean x (mm)" );
  xcent->Draw("same");
  
  c1->cd(2);
  TH1F *frameBycf = gPad->DrawFrame( wz[0]-z_wid*0.1, -8, wz[8]+z_wid*0.1, 6 );
  frameByc->GetXaxis()->SetTitle( "z (mm)" );
  frameByc->GetYaxis()->SetTitle( "Mean y (mm)" );
  ycent->Draw("same");
  
  c1->Print(pdf_tracking);
  c1->Clear();

  // calculated width of projection on x-axis fixed frame
  c1->Divide( 1, 2 );
  c1->cd(1);
  TGraph *xwwaist = new TGraph( 9, wz, xwwid );
  xwwaist->SetMarkerStyle(20);
  //TH1F *frameBxwfw = gPad->DrawFrame( wz[0]-z_wid*0.1, 0, wz[8]+z_wid*0.1, 40 );
  //frameBxwfw->GetXaxis()->SetTitle( "z (mm)" );
  //frameBxwfw->GetYaxis()->SetTitle( "StdDev(x) (mm)" );
  //xwwaist->Draw("Psame");
  xwwaist->Draw("APL");
  
  c1->cd(2);
  TGraph *ywwaist = new TGraph( 9, wz, ywwid );
  ywwaist->SetMarkerStyle(8);
  //TH1F *frameBywfw = gPad->DrawFrame( wz[0]-z_wid*0.1, 0, wz[8]+z_wid*0.1, 25 );
  //frameBywfw->GetXaxis()->SetTitle( "z (mm)" );
  //frameBywfw->GetYaxis()->SetTitle( "StdDev(y) (mm)" );
  //ywwaist->Draw("Psame");
  ywwaist->Draw("APL");
    
  c1->Print(pdf_tracking);
  c1->Clear();

  // finer interpolation
  // 2D-map of Jet0 -finer interpolation- Jet1
  c1->Divide( 3, 3 );
  c1->cd(1);
  if( ON_stg ){
    xy25m_13m_z[0]->Draw();

    xcen[0] = xy25m_13m_z[0]->GetMean(1);
    ycen[0] = xy25m_13m_z[0]->GetMean(2);
    excen[0] = xy25m_13m_z[0]->GetMeanError(1);
    eycen[0] = xy25m_13m_z[0]->GetMeanError(2);
    xwid_l[0] = xy25m_13m_z[0]->GetStdDev(1);
    ywid_l[0] = xy25m_13m_z[0]->GetStdDev(2);
    exwid_l[0] = xy25m_13m_z[0]->GetStdDevError(1);
    eywid_l[0] = xy25m_13m_z[0]->GetStdDevError(2);

    wz[0] = pz[0];
  }else{
    xy25m_13mj1w0->Draw();

    xcen[0] = xy25m_13mj1w0->GetMean(1);
    ycen[0] = xy25m_13mj1w0->GetMean(2);
    excen[0] = xy25m_13mj1w0->GetMeanError(1);
    eycen[0] = xy25m_13mj1w0->GetMeanError(2);
    xwid_l[0] = xy25m_13mj1w0->GetStdDev(1);
    ywid_l[0] = xy25m_13mj1w0->GetStdDev(2);
    exwid_l[0] = xy25m_13mj1w0->GetStdDevError(1);
    eywid_l[0] = xy25m_13mj1w0->GetStdDevError(2);

    wz[0] = z[1];
  }
  xw_max = xwid_l[0];
  yw_max = ywid_l[0];
  xw_min = xwid_l[0];
  yw_min = ywid_l[0];
  //cout << wz[0] << endl;
  //cout << xwid_l[0] << endl;
  //cout << ywid_l[0] << endl;
  
  for( Int_t i = 0; i < 3; i++ ){
    c1->cd(i+2);
    xy25m_13m_fz[i]->Draw();

    xcen[i+1] = xy25m_13m_fz[i]->GetMean(1);
    ycen[i+1] = xy25m_13m_fz[i]->GetMean(2);
    excen[i+1] = xy25m_13m_fz[i]->GetMeanError(1);
    eycen[i+1] = xy25m_13m_fz[i]->GetMeanError(2);
    xwid_l[i+1] = xy25m_13m_fz[i]->GetStdDev(1);
    ywid_l[i+1] = xy25m_13m_fz[i]->GetStdDev(2);
    exwid_l[i+1] = xy25m_13m_fz[i]->GetStdDevError(1);
    eywid_l[i+1] = xy25m_13m_fz[i]->GetStdDevError(2);

    wz[i+1] = fpz[i];

    if( xw_max < xwid_l[i+1] )
      xw_max = xwid_l[i+1];
    if( yw_max < ywid_l[i+1] )
      yw_max = ywid_l[i+1];
    if( xw_min > xwid_l[i+1] )
      xw_min = xwid_l[i+1];
    if( yw_min > ywid_l[i+1] )
      yw_min = ywid_l[i+1];
    //cout << wz[i+1] << endl;
    //cout << xwid_l[i+1] << endl;
    //cout << ywid_l[i+1] << endl;
  }

  c1->cd(5);
  if( ON_stg ){
    xy25m_13m_z[1]->Draw();

    xcen[4] = xy25m_13m_z[1]->GetMean(1);
    ycen[4] = xy25m_13m_z[1]->GetMean(2);
    excen[4] = xy25m_13m_z[1]->GetMeanError(1);
    eycen[4] = xy25m_13m_z[1]->GetMeanError(2);
    xwid_l[4] = xy25m_13m_z[1]->GetStdDev(1);
    ywid_l[4] = xy25m_13m_z[1]->GetStdDev(2);
    exwid_l[4] = xy25m_13m_z[1]->GetStdDevError(1);
    eywid_l[4] = xy25m_13m_z[1]->GetStdDevError(2);

    wz[4] = pz[1];
  }else{
    xy25m_13m_z[0]->Draw();

    xcen[4] = xy25m_13m_z[0]->GetMean(1);
    ycen[4] = xy25m_13m_z[0]->GetMean(2);
    excen[4] = xy25m_13m_z[0]->GetMeanError(1);
    eycen[4] = xy25m_13m_z[0]->GetMeanError(2);
    xwid_l[4] = xy25m_13m_z[0]->GetStdDev(1);
    ywid_l[4] = xy25m_13m_z[0]->GetStdDev(2);
    exwid_l[4] = xy25m_13m_z[0]->GetStdDevError(1);
    eywid_l[4] = xy25m_13m_z[0]->GetStdDevError(2);

    wz[4] = pz[0];
  }
  if( xw_max < xwid_l[4] )
    xw_max = xwid_l[4];
  if( yw_max < ywid_l[4] )
    yw_max = ywid_l[4];
  if( xw_min > xwid_l[4] )
    xw_min = xwid_l[4];
  if( yw_min > ywid_l[4] )
    yw_min = ywid_l[4];
  //cout << wz[4] << endl;
  //cout << xwid_l[4] << endl;
  //cout << ywid_l[4] << endl;
  
  for( Int_t i = 3; i < 6; i++ ){
    c1->cd(i+3);
    xy25m_13m_fz[i]->Draw();

    xcen[i+2] = xy25m_13m_fz[i]->GetMean(1);
    ycen[i+2] = xy25m_13m_fz[i]->GetMean(2);
    excen[i+2] = xy25m_13m_fz[i]->GetMeanError(1);
    eycen[i+2] = xy25m_13m_fz[i]->GetMeanError(2);
    xwid_l[i+2] = xy25m_13m_fz[i]->GetStdDev(1);
    ywid_l[i+2] = xy25m_13m_fz[i]->GetStdDev(2);
    exwid_l[i+2] = xy25m_13m_fz[i]->GetStdDevError(1);
    eywid_l[i+2] = xy25m_13m_fz[i]->GetStdDevError(2);

    wz[i+2] = fpz[i];

    if( xw_max < xwid_l[i+2] )
      xw_max = xwid_l[i+2];
    if( yw_max < ywid_l[i+2] )
      yw_max = ywid_l[i+2];
    if( xw_min > xwid_l[i+2] )
      xw_min = xwid_l[i+2];
    if( yw_min > ywid_l[i+2] )
      yw_min = ywid_l[i+2];
    //cout << wz[i+2] << endl;
    //cout << xwid_l[i+2] << endl;
    //cout << ywid_l[i+2] << endl;
  }

  c1->cd(9);
  if( ON_stg ){
    xy25m_13m_z[2]->Draw();

    xcen[8] = xy25m_13m_z[2]->GetMean(1);
    ycen[8] = xy25m_13m_z[2]->GetMean(2);
    excen[8] = xy25m_13m_z[2]->GetMeanError(1);
    eycen[8] = xy25m_13m_z[2]->GetMeanError(2);
    xwid_l[8] = xy25m_13m_z[2]->GetStdDev(1);
    ywid_l[8] = xy25m_13m_z[2]->GetStdDev(2);
    exwid_l[8] = xy25m_13m_z[2]->GetStdDevError(1);
    eywid_l[8] = xy25m_13m_z[2]->GetStdDevError(2);

    wz[8] = pz[2];
  }else{
    xy25m_13m_z[1]->Draw();

    xcen[8] = xy25m_13m_z[1]->GetMean(1);
    ycen[8] = xy25m_13m_z[1]->GetMean(2);
    excen[8] = xy25m_13m_z[1]->GetMeanError(1);
    eycen[8] = xy25m_13m_z[1]->GetMeanError(2);
    xwid_l[8] = xy25m_13m_z[1]->GetStdDev(1);
    ywid_l[8] = xy25m_13m_z[1]->GetStdDev(2);
    exwid_l[8] = xy25m_13m_z[1]->GetStdDevError(1);
    eywid_l[8] = xy25m_13m_z[1]->GetStdDevError(2);

    wz[8] = pz[1];
  }
  if( xw_max < xwid_l[8] )
    xw_max = xwid_l[8];
  if( yw_max < ywid_l[8] )
    yw_max = ywid_l[8];
  if( xw_min > xwid_l[8] )
    xw_min = xwid_l[8];
  if( yw_min > ywid_l[8] )
    yw_min = ywid_l[8];
  //cout << wz[8] << endl;
  //cout << xwid_l[8] << endl;
  //cout << ywid_l[8] << endl;

  c1->Print(pdf_tracking);
  c1->Clear();

  // width of projection on x-axis
  c1->Divide( 1, 2 );
  c1->cd(1);
  TGraphErrors *xwaistf = new TGraphErrors( 9, wz, xwid_l, ewz, exwid_l );
  //Float_t xw_wid = xw_max - xw_min;
  //Float_t z_wid = fabs( z[0] - z[1] );
  //TH1F *frameBxw = gPad->DrawFrame( wz[0]-z_wid*0.1, xw_min-xw_wid*0.5, wz[8]+z_wid*0.1, xw_max+xw_wid*0.5 );
  //frameBxw->GetXaxis()->SetTitle( "z (mm)" );
  //frameBxw->GetYaxis()->SetTitle( "StdDev(x) (mm)" );
  xwaistf->Draw("AL");
  
  c1->cd(2);
  TGraphErrors *ywaistf = new TGraphErrors( 9, wz, ywid_l, ewz, eywid_l );
  //Float_t yw_wid = yw_max - yw_min;
  //TH1F *frameByw = gPad->DrawFrame( wz[0]-z_wid*0.1, yw_min-yw_wid*0.5, wz[8]+z_wid*0.1, yw_max+yw_wid*0.5 );
  //frameByw->GetXaxis()->SetTitle( "z (mm)" );
  //frameByw->GetYaxis()->SetTitle( "StdDev(y) (mm)" );
  ywaistf->Draw("AL");
  
  c1->Print(pdf_tracking);
  c1->Clear();
  
  // width of projection on x-axis
  c1->Divide( 1, 2 );
  c1->cd(1);
  para->SetParameter(2, xw_min );
  if( ON_stg ){
    para->SetParameter(0, (xw_max-xw_min)/1000/1000 );
    para->SetParameter(1, 1000 );
  }else{
    para->SetParameter(0, (xw_max-xw_min)/1000/1000 );
    para->SetParameter(1, -3500 );
  }
  xwaistf->Draw("AP");
  xwaistf->Fit("para");
  
  c1->cd(2);
  para->SetParameter(2, yw_min );
  if( ON_stg ){
    para->SetParameter(0, (yw_max-yw_min)/1000/1000 );
    para->SetParameter(1, 1000 );
  }else{
    para->SetParameter(0, (yw_max-yw_min)/1000/1000 );
    para->SetParameter(1, -3500 );
  }
  ywaistf->Draw("AP");
  ywaistf->Fit("para");
  

  for( Int_t i = 0; i < Nexp; i++ ){
    c1->Print(pdf_tracking);
    c1->Clear();

    // extrapolate to upper stream
    c1->Divide( 2, 2 );
    c1->cd(1);
    xy25m_13m_exz_E0[i]->Draw();
    c1->cd(2);
    xy25m_13m_exz_E1[i]->Draw();
    c1->cd(3);
    xy25m_13m_exz_E2[i]->Draw();
    c1->cd(4);
    xy25m_13m_exz_E3[i]->Draw();
  }
  
  // pole shade of dounw stream p0
  /*for( Int_t i = 0; i < 4; i++ ){
    c1->Print(pdf_tracking);
    c1->Clear();

    // extrapolate to upper stream
    c1->Divide( 2, 2 );
    c1->cd(1);
    xy25m_13m_j0p0_exz_E0[i]->Draw();
    c1->cd(2);
    xy25m_13m_j0p0_exz_E1[i]->Draw();
    c1->cd(3);
    xy25m_13m_j0p0_exz_E2[i]->Draw();
    c1->cd(4);
    xy25m_13m_j0p0_exz_E3[i]->Draw();
    }*/
  
  c1->Print(pdf_tracking);
  c1->Clear();

  Int_t ID1st;
  // ON_stg == 1   Bu  Bl  T0  QSF   QSD  T4   Jet0 Jet1  LG
  // ON_stg == 0   Bu  Bl  T0  Jet1  QSF  QSD  T4   Jet0  LG
  //               1   2   3   4     5    6    7    8     9
  c1->Divide( 3, 3 );// pedestal
  c1->cd(1);
  xy25m_13m_exz_E0[5]->Draw();
  c1->cd(2);
  xy25m_13m_exz_E0[4]->Draw();
  c1->cd(3);
  xy25m_13m_exz_E0[3]->Draw();
  if( ON_stg ) ID1st = 6;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(4);
    //xy25m_13mj1w0_cE[0]->Draw();
    xy25m_13m_detz_E0[1]->Draw();
    ID1st = 7;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_exz_E0[i]->Draw();
  }
  if( ON_stg ) c1->cd(7);
  else c1->cd(8);
  //xy25m_13mj0w1_cE[0]->Draw();
  xy25m_13m_detz_E0[0]->Draw();
  if( ON_stg ){
    c1->cd(8);
    //xy25m_13mj1w0_cE[0]->Draw();
    xy25m_13m_detz_E0[1]->Draw();
  }
  c1->cd(9);
  //TH1F *frame0 = gPad->DrawFrame( -50, -50, 50, 50 );
  //xy_LG_E[0]->Draw("same");
  c1->cd(9);
  xy_LG_E[0]->Draw();
  
  c1->Print(pdf_tracking);
  c1->Clear();
  
  c1->Divide( 3, 3 );// low energy
  c1->cd(1);
  xy25m_13m_exz_E1[5]->Draw();
  c1->cd(2);
  xy25m_13m_exz_E1[4]->Draw();
  c1->cd(3);
  xy25m_13m_exz_E1[3]->Draw();
  if( ON_stg ) ID1st = 6;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(4);
    //xy25m_13mj1w0_cE[1]->Draw();
    xy25m_13m_detz_E1[1]->Draw();
    ID1st = 7;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_exz_E1[i]->Draw();
  }
  if( ON_stg ) c1->cd(7);
  else c1->cd(8);
  //xy25m_13mj0w1_cE[1]->Draw();
  xy25m_13m_detz_E1[0]->Draw();
  if( ON_stg ){
    c1->cd(8);
    //xy25m_13mj1w0_cE[1]->Draw();
    xy25m_13m_detz_E1[1]->Draw();
  }
  c1->cd(9);
  //TH1F *frame1 = gPad->DrawFrame( -50, -50, 50, 50 );
  //xy_LG_E[1]->Draw("same");
  c1->cd(9);
  xy_LG_E[1]->Draw();
  
  c1->Print(pdf_tracking);
  c1->Clear();
  
  c1->Divide( 3, 3 );// target energy
  c1->cd(1);
  xy25m_13m_exz_E2[5]->Draw();
  c1->cd(2);
  xy25m_13m_exz_E2[4]->Draw();
  c1->cd(3);
  xy25m_13m_exz_E2[3]->Draw();
  if( ON_stg ) ID1st = 6;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(4);
    //xy25m_13mj1w0_cE[2]->Draw();
    xy25m_13m_detz_E2[1]->Draw();
    ID1st = 7;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_exz_E2[i]->Draw();
  }
  if( ON_stg ) c1->cd(7);
  else c1->cd(8);
  //xy25m_13mj0w1_cE[2]->Draw();
  xy25m_13m_detz_E2[0]->Draw();
  if( ON_stg ){
    c1->cd(8);
    //xy25m_13mj1w0_cE[2]->Draw();
    xy25m_13m_detz_E2[1]->Draw();
  }
  c1->cd(9);
  //TH1F *frame2 = gPad->DrawFrame( -50, -50, 50, 50 );
  //xy_LG_E[2]->Draw("same");
  c1->cd(9);
  xy_LG_E[2]->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();
  
  c1->Divide( 3, 3 );// high energy
  c1->cd(1);
  xy25m_13m_exz_E3[5]->Draw();
  c1->cd(2);
  xy25m_13m_exz_E3[4]->Draw();
  c1->cd(3);
  xy25m_13m_exz_E3[3]->Draw();
  if( ON_stg ) ID1st = 6;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(4);
    //xy25m_13mj1w0_cE[3]->Draw();
    xy25m_13m_detz_E3[1]->Draw();
    ID1st = 7;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_exz_E3[i]->Draw();
  }
  if( ON_stg ) c1->cd(7);
  else c1->cd(8);
  //xy25m_13mj0w1_cE[3]->Draw();
  xy25m_13m_detz_E3[0]->Draw();
  if( ON_stg ){
    c1->cd(8);
    //xy25m_13mj1w0_cE[3]->Draw();
    xy25m_13m_detz_E3[1]->Draw();
  }
  c1->cd(9);
  //TH1F *frame3 = gPad->DrawFrame( -50, -50, 50, 50 );
  //xy_LG_E[3]->Draw("same");
  c1->cd(9);
  xy_LG_E[3]->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();

  //////////
  c1->Divide( 3, 3 );// pedestal
  c1->cd(1);
  xy25m_13m_exz_LGxneg_E0[5]->Draw();
  c1->cd(2);
  xy25m_13m_exz_LGxneg_E0[4]->Draw();
  c1->cd(3);
  xy25m_13m_exz_LGxneg_E0[3]->Draw();
  if( ON_stg ) ID1st = 6;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(4);
    //xy25m_13mj1w0_cE[0]->Draw();
    xy25m_13m_detz_LGxneg_E0[1]->Draw();
    ID1st = 7;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_exz_LGxneg_E0[i]->Draw();
  }
  if( ON_stg ) c1->cd(7);
  else c1->cd(8);
  //xy25m_13mj0w1_cE[0]->Draw();
  xy25m_13m_detz_LGxneg_E0[0]->Draw();
  if( ON_stg ){
    c1->cd(8);
    //xy25m_13mj1w0_cE[0]->Draw();
    xy25m_13m_detz_LGxneg_E0[1]->Draw();
  }
  c1->cd(9);
  //TH1F *frame0 = gPad->DrawFrame( -50, -50, 50, 50 );
  //xy_LG_E[0]->Draw("same");
  c1->cd(9);
  xy_LG_E[0]->Draw();
  
  c1->Print(pdf_tracking);
  c1->Clear();
  
  c1->Divide( 3, 3 );// low energy
  c1->cd(1);
  xy25m_13m_exz_LGxneg_E1[5]->Draw();
  c1->cd(2);
  xy25m_13m_exz_LGxneg_E1[4]->Draw();
  c1->cd(3);
  xy25m_13m_exz_LGxneg_E1[3]->Draw();
  if( ON_stg ) ID1st = 6;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(4);
    //xy25m_13mj1w0_cE[1]->Draw();
    xy25m_13m_detz_LGxneg_E1[1]->Draw();
    ID1st = 7;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_exz_LGxneg_E1[i]->Draw();
  }
  if( ON_stg ) c1->cd(7);
  else c1->cd(8);
  //xy25m_13mj0w1_cE[1]->Draw();
  xy25m_13m_detz_LGxneg_E1[0]->Draw();
  if( ON_stg ){
    c1->cd(8);
    //xy25m_13mj1w0_cE[1]->Draw();
    xy25m_13m_detz_LGxneg_E1[1]->Draw();
  }
  c1->cd(9);
  //TH1F *frame1 = gPad->DrawFrame( -50, -50, 50, 50 );
  //xy_LG_E[1]->Draw("same");
  c1->cd(9);
  xy_LG_E[1]->Draw();
  
  c1->Print(pdf_tracking);
  c1->Clear();
  
  c1->Divide( 3, 3 );// target energy
  c1->cd(1);
  xy25m_13m_exz_LGxneg_E2[5]->Draw();
  c1->cd(2);
  xy25m_13m_exz_LGxneg_E2[4]->Draw();
  c1->cd(3);
  xy25m_13m_exz_LGxneg_E2[3]->Draw();
  if( ON_stg ) ID1st = 6;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(4);
    //xy25m_13mj1w0_cE[2]->Draw();
    xy25m_13m_detz_LGxneg_E2[1]->Draw();
    ID1st = 7;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_exz_LGxneg_E2[i]->Draw();
  }
  if( ON_stg ) c1->cd(7);
  else c1->cd(8);
  //xy25m_13mj0w1_cE[2]->Draw();
  xy25m_13m_detz_LGxneg_E2[0]->Draw();
  if( ON_stg ){
    c1->cd(8);
    //xy25m_13mj1w0_cE[2]->Draw();
    xy25m_13m_detz_LGxneg_E2[1]->Draw();
  }
  c1->cd(9);
  //TH1F *frame2 = gPad->DrawFrame( -50, -50, 50, 50 );
  //xy_LG_E[2]->Draw("same");
  c1->cd(9);
  xy_LG_E[2]->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();
  
  c1->Divide( 3, 3 );// high energy
  c1->cd(1);
  xy25m_13m_exz_LGxneg_E3[5]->Draw();
  c1->cd(2);
  xy25m_13m_exz_LGxneg_E3[4]->Draw();
  c1->cd(3);
  xy25m_13m_exz_LGxneg_E3[3]->Draw();
  if( ON_stg ) ID1st = 6;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(4);
    //xy25m_13mj1w0_cE[3]->Draw();
    xy25m_13m_detz_LGxneg_E3[1]->Draw();
    ID1st = 7;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_exz_LGxneg_E3[i]->Draw();
  }
  if( ON_stg ) c1->cd(7);
  else c1->cd(8);
  //xy25m_13mj0w1_cE[3]->Draw();
  xy25m_13m_detz_LGxneg_E3[0]->Draw();
  if( ON_stg ){
    c1->cd(8);
    //xy25m_13mj1w0_cE[3]->Draw();
    xy25m_13m_detz_LGxneg_E3[1]->Draw();
  }
  c1->cd(9);
  //TH1F *frame3 = gPad->DrawFrame( -50, -50, 50, 50 );
  //xy_LG_E[3]->Draw("same");
  c1->cd(9);
  xy_LG_E[3]->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();
  /////////

  //////////////]
  // ON_stg == 1   T0  QSF  QSD   T4 Jet0 Jet1 LG
  // ON_stg == 0   T0 Jet1  QSF  QSD   T4 Jet0 LG
  //               1     2    3    4    5    6  7
  /*c1->Divide( 3, 3 );// pedestal
  c1->cd(1);
  xy25m_13m_exz_E0[3]->Draw();
  if( ON_stg ) ID1st = 4;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(2);
    xy25m_13m_detz_E0[1]->Draw();
    ID1st = 5;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_exz_E0[i]->Draw();
  }
  if( ON_stg ) c1->cd(5);
  else c1->cd(6);
  xy25m_13m_detz_E0[0]->Draw();
  if( ON_stg ){
    c1->cd(6);
    xy25m_13m_detz_E0[1]->Draw();
  }
  for( Int_t i = 0; i < 2; i++ ){
    c1->cd(8+i);
    xy25m_13m_exz2_E0[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();
  
  c1->Divide( 3, 3 );// low energy
  c1->cd(1);
  xy25m_13m_exz_E1[3]->Draw();
  if( ON_stg ) ID1st = 4;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(2);
    xy25m_13m_detz_E1[1]->Draw();
    ID1st = 5;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_exz_E1[i]->Draw();
  }
  if( ON_stg ) c1->cd(5);
  else c1->cd(6);
  xy25m_13m_detz_E1[0]->Draw();
  if( ON_stg ){
    c1->cd(6);
    xy25m_13m_detz_E1[1]->Draw();
  }
  for( Int_t i = 0; i < 2; i++ ){
    c1->cd(8+i);
    xy25m_13m_exz2_E1[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();
  
  c1->Divide( 3, 3 );// target energy
  c1->cd(1);
  xy25m_13m_exz_E2[3]->Draw();
  if( ON_stg ) ID1st = 4;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(2);
    xy25m_13m_detz_E2[1]->Draw();
    ID1st = 5;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_exz_E2[i]->Draw();
  }
  if( ON_stg ) c1->cd(5);
  else c1->cd(6);
  xy25m_13m_detz_E2[0]->Draw();
  if( ON_stg ){
    c1->cd(6);
    xy25m_13m_detz_E2[1]->Draw();
  }
  for( Int_t i = 0; i < 2; i++ ){
    c1->cd(8+i);
    xy25m_13m_exz2_E2[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();
  
  c1->Divide( 3, 3 );// high energy
  c1->cd(1);
  xy25m_13m_exz_E3[3]->Draw();
  if( ON_stg ) ID1st = 4;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(2);
    xy25m_13m_detz_E3[1]->Draw();
    ID1st = 5;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_exz_E3[i]->Draw();
  }
  if( ON_stg ) c1->cd(5);
  else c1->cd(6);
  xy25m_13m_detz_E3[0]->Draw();
  if( ON_stg ){
    c1->cd(6);
    xy25m_13m_detz_E3[1]->Draw();
  }
  for( Int_t i = 0; i < 2; i++ ){
    c1->cd(8+i);
    xy25m_13m_exz2_E3[i]->Draw();
  }

  c1->Print(pdf_tracking);
  c1->Clear();*/
  
  //////////////
  // pole shade of down stream p0
  // ON_stg == 1   T0  QSF  QSD   T4 Jet0 Jet1 LG
  // ON_stg == 0   T0 Jet1  QSF  QSD   T4 Jet0 LG
  //               1     2    3    4    5    6  7
  /*c1->Divide( 3, 3 );
  c1->cd(1);
  xy25m_13m_j0p0_exz_E0[3]->Draw();
  if( ON_stg ) ID1st = 4;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(2);
    xy25m_13mj1w0j0p0_cE[0]->Draw();
    ID1st = 5;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_j0p0_exz_E0[i]->Draw();
  }
  if( ON_stg ) c1->cd(5);
  else c1->cd(6);
  xy25m_13mj0w1j0p0_cE[0]->Draw();
  if( ON_stg ){
    c1->cd(6);
    xy25m_13mj1w0j0p0_cE[0]->Draw();
  }
  c1->cd(7);
  TH1F *frame0dp = gPad->DrawFrame( -50, -50, 50, 50 );
  xy_LG_j0p0_E[0]->Draw("same");
  c1->cd(9);
  xy_LG_j0p0_E[0]->Draw();
  
  c1->Print(pdf_tracking);
  c1->Clear();
  
  c1->Divide( 3, 3 );// low energy
  c1->cd(1);
  xy25m_13m_j0p0_exz_E1[3]->Draw();
  if( ON_stg ) ID1st = 4;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(2);
    xy25m_13mj1w0j0p0_cE[1]->Draw();
    ID1st = 5;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_j0p0_exz_E1[i]->Draw();
  }
  if( ON_stg ) c1->cd(5);
  else c1->cd(6);
  xy25m_13mj0w1j0p0_cE[1]->Draw();
  if( ON_stg ){
    c1->cd(6);
    xy25m_13mj1w0j0p0_cE[1]->Draw();
  }
  c1->cd(7);
  TH1F *frame1dp = gPad->DrawFrame( -50, -50, 50, 50 );
  xy_LG_j0p0_E[1]->Draw("same");
  c1->cd(9);
  xy_LG_j0p0_E[1]->Draw();
  
  c1->Print(pdf_tracking);
  c1->Clear();
  
  c1->Divide( 3, 3 );// target energy
  c1->cd(1);
  xy25m_13m_j0p0_exz_E2[3]->Draw();
  if( ON_stg ) ID1st = 4;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(2);
    xy25m_13mj1w0j0p0_cE[2]->Draw();
    ID1st = 5;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_j0p0_exz_E2[i]->Draw();
  }
  if( ON_stg ) c1->cd(5);
  else c1->cd(6);
  xy25m_13mj0w1j0p0_cE[2]->Draw();
  if( ON_stg ){
    c1->cd(6);
    xy25m_13mj1w0j0p0_cE[2]->Draw();
  }
  c1->cd(7);
  TH1F *frame2dp = gPad->DrawFrame( -50, -50, 50, 50 );
  xy_LG_j0p0_E[2]->Draw("same");
  c1->cd(9);
  xy_LG_j0p0_E[2]->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();
  
  c1->Divide( 3, 3 );// high energy
  c1->cd(1);
  xy25m_13m_j0p0_exz_E3[3]->Draw();
  if( ON_stg ) ID1st = 4;// i=2, 2; i=1, 3; i=0, 4
  else{
    c1->cd(2);
    xy25m_13mj1w0j0p0_cE[3]->Draw();
    ID1st = 5;// i=2, 3; i=1, 2; i=0, 5
  }
  for( Int_t i = 2; i >= 0; i-- ){
    c1->cd(ID1st - i);
    xy25m_13m_j0p0_exz_E3[i]->Draw();
  }
  if( ON_stg ) c1->cd(5);
  else c1->cd(6);
  xy25m_13mj0w1j0p0_cE[3]->Draw();
  if( ON_stg ){
    c1->cd(6);
    xy25m_13mj1w0j0p0_cE[3]->Draw();
  }
  c1->cd(7);
  TH1F *frame3dp = gPad->DrawFrame( -50, -50, 50, 50 );
  xy_LG_j0p0_E[3]->Draw("same");
  c1->cd(9);
  xy_LG_j0p0_E[3]->Draw();

  c1->Print(pdf_tracking);
  c1->Clear();*/
  

  /////// LG //////
  if( LG_on ){
    //c1->Print(pdf_tracking);
    //c1->Clear();

    // 2D map on LG surface
    c1->Divide( 2, 2 );
    c1->cd(1);
    LG->Draw();
    Int_t lg_entry;
    Float_t lg_mean = LG->GetMean();
    Float_t lg_rms = LG->GetStdDev();
    Float_t ll = lg_mean - lg_rms;
    Float_t ul = lg_mean + 3*lg_rms;

    //fitLGT->SetParameter( 1, lg_mean*2 );
    //fitLGT->SetParameter( 2, lg_rms/2 );
    //fitLGT->SetParameter( 3, lg_rms/2 );
    //fitLGT->SetParameter( 4, 20000 );
    //fitLGT->SetParameter( 5, -0.01 );

    //LG->Fit( "gaus", "", "", ll, ul );
    //LG->Fit( "fitLG", "", "", ll, ul );
    //LG->Fit( "fitLG", "LLE", "", 20, ul );
    //LG->Fit( "fitLGT", "L", "", 7, ul );
    
    //c1->cd(2);
    //LG_tdc_j0x->Draw();

    c1->cd(4);
    xy_LG->Draw();
    TLine *lgb[4];
    lgb[0] = new TLine( -60, -100, -60, 100 );
    lgb[1] = new TLine( 60, -100, 60, 100 );
    lgb[2] = new TLine( -100, -60, 100, -60 );
    lgb[3] = new TLine( -100, 60, 100, 60 );
    for( Int_t i = 0; i < 4; i++ ){
      lgb[i]->SetLineWidth(2);
      lgb[i]->Draw("same");
    }
    
    c1->cd(2);
    x_LG->Draw();
    //c1->cd(3);
    //LGtdc->Draw();
    c1->cd(3);
    y_LG->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    //LG xy 
    c1->Divide( 2, 2 );
    for( Int_t i = 0; i < 4; i++ ){
      c1->cd(i+1);
      xy_LG_E[i]->Draw();
    }
    
    c1->Print(pdf_tracking);
    c1->Clear();

    //LG xy 
    c1->Divide( 2, 2 );
    for( Int_t i = 0; i < 4; i++ ){
      c1->cd(i+1);
      xy_LG_j1xneg_E[i]->Draw();
    }
    
    c1->Print(pdf_tracking);
    c1->Clear();

    // fit x-sliced LG PH 
    c1->Divide( 4, 4 );
    c1->cd(1);
    x_LG->Draw();

    TH1D *x_LG_sl[14];
    x_LG_sl[0] = x_LG->ProjectionY( "25m_13m_-70_-60", 16, 20 );
    x_LG_sl[1] = x_LG->ProjectionY( "25m_13m_-60_-50", 21, 25 );
    x_LG_sl[2] = x_LG->ProjectionY( "25m_13m_-50_-40", 26, 30 );
    x_LG_sl[3] = x_LG->ProjectionY( "25m_13m_-40_-30", 31, 35 );
    x_LG_sl[4] = x_LG->ProjectionY( "25m_13m_-30_-20",   36, 40 );
    x_LG_sl[5] = x_LG->ProjectionY( "25m_13m_-20_-10",  41, 45) ;
    x_LG_sl[6] = x_LG->ProjectionY( "25m_13m_-10_0", 46, 50 );
    x_LG_sl[7] = x_LG->ProjectionY( "25m_13m_0_10", 51, 55 );
    x_LG_sl[8] = x_LG->ProjectionY( "25m_13m_10_20", 56, 60 );
    x_LG_sl[9] = x_LG->ProjectionY( "25m_13m_20_30", 61, 65 );
    x_LG_sl[10] = x_LG->ProjectionY( "25m_13m_30_40", 66, 70 );
    x_LG_sl[11] = x_LG->ProjectionY( "25m_13m_40_50", 71, 75 );
    x_LG_sl[12] = x_LG->ProjectionY( "25m_13m_50_60", 76, 80 );
    x_LG_sl[13] = x_LG->ProjectionY( "25m_13m_60_70", 81, 85 );

    TH1D *x_LG_slc = x_LG->ProjectionY( "25m_13m_-10_10", 46, 55 );
    
    Float_t x_pos[14] = { -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65 };
    Float_t ex_pos[14] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    Float_t ph_lg[14];
    Float_t eph_lg[14];
    Float_t x_pos_fit[14];
    Float_t ex_pos_fit[14];
    Float_t ph_lg_fit[14];
    Float_t eph_lg_fit[14];
    Int_t n_fit = 0;
    
    Int_t lg_n = 0;
    Float_t lg_max = 0;
    Float_t max_lmax = 0;
    Float_t min_lmax = 0;
    Int_t fit_r;
    //Int_t ent[80];
    // (0-320)/80 1st:0-3, 2nd:3-7, 3rd:8-11, 4th:12-15, 5th:16-19 
    Int_t ent3, ent6, entp, ent, pp;
    // ent1 = e0, ent5 = e0*exp(d*20) ent5/ent1 = exp(d*20) ln(ent5/ent1) = d*20
    Float_t dump;
    Float_t slope;
    for( Int_t i = 0; i < 14; i++ ){
      c1->cd( i+2 );
      ph_lg[i] = 0;
      eph_lg[i] = 0;
      lg_entry = x_LG_sl[i]->GetEntries();
      if( lg_entry > Nmin_fit ){
	lg_mean = x_LG_sl[i]->GetMean();
	lg_rms = x_LG_sl[i]->GetStdDev();
	ul = lg_mean + 3*lg_rms;
	entp = 0;
	pp = 0;
	for( Int_t j = 3; j < 80; j++ ){
	  if( entp < x_LG_sl[i]->GetBinContent(j) ){ entp = x_LG_sl[i]->GetBinContent(j); pp = j * 4; }
	}
	ent3 = x_LG_sl[i]->GetBinContent(3);
	ent6 = x_LG_sl[i]->GetBinContent(6);
	if( ent6 == 0 ) dump = -10;
	else dump = log(((float)ent6)/((float)ent3))/12;

	fitLGT->SetParameter( 0, entp );
	fitLGT->SetParameter( 1, pp );
	fitLGT->SetParameter( 2, lg_rms/4 );
	fitLGT->SetParameter( 3, lg_rms/5 );
	fitLGT->SetParameter( 4, ent3*1.2 );
	fitLGT->SetParameter( 5, dump );
      
	fit_r = x_LG_sl[i]->Fit( "fitLGT", "", "", 10, ul );

	bg[i]->SetParameter( 0, fitLGT->GetParameter( 4 ) );
	bg[i]->SetParameter( 1, fitLGT->GetParameter( 5 ) );
	bg[i]->Draw("same");
	
	ph_lg[i] = fitLGT->GetParameter(1);
	eph_lg[i] = fitLGT->GetParError(1);
	//if( ph_lg[i]/eph_lg[i] < 0.1 ){
	lg_max+= ph_lg[i];
	lg_n++;
	//}
	//if( ph_lg[i] > 0 ){
	if( fit_r == 0 ){
	  if( max_lmax < ph_lg[i] ) max_lmax = ph_lg[i];
	  if( min_lmax == 0 ) min_lmax = ph_lg[i];
	  else if( min_lmax > ph_lg[i] ) min_lmax = ph_lg[i];
	}
	if( fit_r == 0 ){
	  x_pos_fit[n_fit] = x_pos[i];
	  ex_pos_fit[n_fit] = ex_pos[i];
	  ph_lg_fit[n_fit] = ph_lg[i];
	  eph_lg_fit[n_fit] = eph_lg[i];
	  n_fit++;
	}
      }else{
	x_LG_sl[i]->Draw();
      }
    }
    if( lg_n == 0 )
      lg_max = 1;
    else
      lg_max = lg_max/lg_n;

    c1->cd(16);
    TGraphErrors *lg_xdep = new TGraphErrors( 14, x_pos, ph_lg, ex_pos, eph_lg );
    //lg_xdep->Draw("AL" );
    TH1F *frameLGN = gPad->DrawFrame( -100, 0, 100, 320 );
    frameLGN->GetXaxis()->SetTitle( "x (mm)" );
    frameLGN->GetYaxis()->SetTitle( "LG Pulse Height (count)" );
    lg_xdep->Draw("same" );

    c1->Print(pdf_tracking);
    c1->Clear();

    // LG
    c1->Divide( 2, 2 );
    c1->cd(1);
    LG->Draw();
    c1->cd(2);
    LGff->Draw();

    c1->cd(3);
    gPad->SetLogy(1);
    LG->Draw();
    c1->cd(4);
    gPad->SetLogy(1);
    LGff->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    // LG
    c1->Divide( 2, 2 );
    c1->cd(1);
    LG->Draw();
    c1->cd(2);
    LG2->Draw();

    c1->cd(3);
    gPad->SetLogy(1);
    LG->Draw();
    c1->cd(4);
    gPad->SetLogy(1);
    LG2->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    // LG
    c1->Divide( 2, 2 );
    c1->cd(1);
    LGpeak->Draw();

    c1->cd(3);
    LGvalley->Draw();

    c1->cd(4);
    LGpeak2->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    //gStyle->SetOptFit(1111);
    // LG PH fit
    c1->Divide( 1, 1 );
    c1->cd(1);
    LG->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1, 1 );
    c1->cd(1);
    gPad->SetLogy(1);
    LG->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1, 1 );
    c1->cd(1);
    LG->Draw();

    lg_mean = LG->GetMean();
    lg_rms = LG->GetStdDev();
    ul = lg_mean + 3*lg_rms;
    entp = 0;
    pp = 0;
    for( Int_t j = 3; j < 80; j++ ){
      if( entp < LG->GetBinContent(j) ){ entp = LG->GetBinContent(j); pp = j * 4;}
    }
    ent3 = LG->GetBinContent(3);
    ent6 = LG->GetBinContent(6);
    if( ent6 == 0 ) dump = -10;
    else dump = log(((float)ent6)/((float)ent3))/12;

    /*fitLGT->SetParameter( 0, entp );
    fitLGT->SetParameter( 1, pp );
    fitLGT->SetParameter( 2, lg_rms/4 );
    fitLGT->SetParameter( 3, lg_rms/5 );
    //fitLGT->SetParameter( 4, 200 );
    //fitLGT->SetParameter( 5, -0.01 );
    fitLGT->SetParameter( 4, ent3*1.2 );
    fitLGT->SetParameter( 5, dump );*/
    fitLG->SetParameter( 0, entp );
    fitLG->SetParameter( 1, pp );
    fitLG->SetParameter( 2, lg_rms/4 );
    fitLG->SetParameter( 3, lg_rms/5 );
    fitLG->SetParameter( 4, entp/1000 );
    fitLG->SetParameter( 5, pp*2 );
    fitLG->SetParameter( 6, lg_rms/4 );
    fitLG->SetParameter( 7, ent3*1.2 );
    //fitLG->SetParameter( 8, dump );
    slope = ent3*1.2/pp;
    fitLG->SetParameter( 8, slope );

    //LG->Fit( "gaus", "", "", ll, ul );
    //LG->Fit( "fitLG", "", "", ll, ul );
    //LG->Fit( "fitLG", "LLE", "", 20, ul );
    //LG->Fit( "fitLGT", "L", "", 20, ul );
    //LG->Fit( "fitLGT", "", "", 10, ul );
    LG->Fit( "fitLG", "", "", 10, 320 );

    //bg[0]->SetParameter( 0, fitLGT->GetParameter( 4 ) );
    //bg[0]->SetParameter( 1, fitLGT->GetParameter( 5 ) );
    bg[0]->SetParameter( 0, fitLG->GetParameter( 7 ) );
    bg[0]->SetParameter( 1, fitLG->GetParameter( 8 ) );
    bg[0]->Draw("same");

    Int_t maxc = 0;
    for( Int_t i = 0; i < 100; ++i ){
      if( maxc < LG->GetBinContent(i) ) maxc = LG->GetBinContent(i);
    }
    
    /*Float_t pa1 = fitLGT->GetParameter(1);
    Float_t pa2 = fitLGT->GetParameter(2);
    Float_t pa3 = fitLGT->GetParameter(3);
    Float_t epa1 = fitLGT->GetParError(1);
    Float_t epa2 = fitLGT->GetParError(2);
    Float_t epa3 = fitLGT->GetParError(3);
    TString spa1 = to_string( pa1 );
    TString spa2 = to_string( pa2 );
    TString spa3 = to_string( pa3 );
    TString espa1 = to_string( epa1 );
    TString espa2 = to_string( epa2 );
    TString espa3 = to_string( epa3 );
    TString peak = "Peak = " + spa1 + " +- " + espa1;
    TString sigl = "SigL = " + spa2 + " +- " + espa2;
    TString sigr = "SigR = " + spa3 + " +- " + espa3;
    TText *tpeak = new TText( 50, ((float)maxc)*0.9, peak );
    TText *tsigl = new TText( 50, ((float)maxc)*0.8, sigl );
    TText *tsigr = new TText( 50, ((float)maxc)*0.7, sigr );
    tpeak->Draw("same");
    tsigl->Draw("same");
    tsigr->Draw("same");*/

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1, 1 );
    c1->cd(1);
    gPad->SetLogy(1);
    LG->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1, 1 );
    c1->cd(1);
    x_LG_slc->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1, 1 );
    c1->cd(1);
    gPad->SetLogy(1);
    x_LG_slc->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    lg_mean = x_LG_slc->GetMean();
    lg_rms = x_LG_slc->GetStdDev();
    ul = lg_mean + 3*lg_rms;
    entp = 0;
    pp =0;
    for( Int_t j = 3; j < 80; j++ ){
      if( entp < x_LG_slc->GetBinContent(j) ){ entp = x_LG_slc->GetBinContent(j); pp = j * 4; }
    }
    ent3 = x_LG_slc->GetBinContent(3);
    ent6 = x_LG_slc->GetBinContent(6);
    if( ent6 == 0 ) dump = -10;
    else dump = log(((float)ent6)/((float)ent3))/12;

    /*fitLGT->SetParameter( 0, entp );
    fitLGT->SetParameter( 1, pp );
    fitLGT->SetParameter( 2, lg_rms/4 );
    fitLGT->SetParameter( 3, lg_rms/5 );
    fitLGT->SetParameter( 4, ent3*1.2 );
    fitLGT->SetParameter( 5, dump );*/
    fitLG->SetParameter( 0, entp );
    fitLG->SetParameter( 1, pp );
    fitLG->SetParameter( 2, lg_rms/4 );
    fitLG->SetParameter( 3, lg_rms/5 );
    fitLG->SetParameter( 4, entp/1000 );
    fitLG->SetParameter( 5, pp*2 );
    fitLG->SetParameter( 6, lg_rms/4 );
    fitLG->SetParameter( 7, ent3*1.2 );
    fitLG->SetParameter( 8, dump );
    slope = ent3*1.2/pp;
    fitLG->SetParameter( 8, slope );

    //x_LG_slc->Fit( "fitLGT", "", "", 10, ul );
    x_LG_slc->Fit( "fitLG", "", "", 10, 320 );

    //bg[0]->SetParameter( 0, fitLGT->GetParameter( 4 ) );
    //bg[0]->SetParameter( 1, fitLGT->GetParameter( 5 ) );
    bg[0]->SetParameter( 0, fitLGT->GetParameter( 7 ) );
    bg[0]->SetParameter( 1, fitLGT->GetParameter( 8 ) );
    bg[0]->Draw("same");

    maxc = 0;
    for( Int_t i = 0; i < 100; ++i ){
      if( maxc < x_LG_slc->GetBinContent(i) ) maxc = x_LG_slc->GetBinContent(i);
    }
    
    /*Float_t pa1c = fitLGT->GetParameter(1);
    Float_t pa2c = fitLGT->GetParameter(2);
    Float_t pa3c = fitLGT->GetParameter(3);
    Float_t epa1c = fitLGT->GetParError(1);
    Float_t epa2c = fitLGT->GetParError(2);
    Float_t epa3c = fitLGT->GetParError(3);
    TString spa1c = to_string( pa1c );
    TString spa2c = to_string( pa2c );
    TString spa3c = to_string( pa3c );
    TString espa1c = to_string( epa1c );
    TString espa2c = to_string( epa2c );
    TString espa3c = to_string( epa3c );
    TString peakc = "Peak = " + spa1c + " +- " + espa1c;
    TString siglc = "SigL = " + spa2c + " +- " + espa2c;
    TString sigrc = "SigR = " + spa3c + " +- " + espa3c;
    TText *tpeakc = new TText( 50, ((float)maxc)*0.9, peakc );
    TText *tsiglc = new TText( 50, ((float)maxc)*0.8, siglc );
    TText *tsigrc = new TText( 50, ((float)maxc)*0.7, sigrc );
    tpeakc->Draw("same");
    tsiglc->Draw("same");
    tsigrc->Draw("same");*/

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1, 1 );
    c1->cd(1);
    gPad->SetLogy(1);
    x_LG_slc->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1, 1 );
    c1->cd(1);
    xmode_LG->Draw();

    lg_mean = xmode_LG->GetMean();
    lg_rms = xmode_LG->GetStdDev();
    ul = lg_mean + 3*lg_rms;

    fitLGT->SetParameter( 1, lg_mean );
    fitLGT->SetParameter( 2, lg_rms/4 );
    fitLGT->SetParameter( 3, lg_rms/4 );
    fitLGT->SetParameter( 4, 200 );
    fitLGT->SetParameter( 5, -0.01 );

    xmode_LG->Fit( "fitLGT", "", "", 25, ul );

    bg[0]->SetParameter( 0, fitLGT->GetParameter( 4 ) );
    bg[0]->SetParameter( 1, fitLGT->GetParameter( 5 ) );
    bg[0]->Draw("same");

    maxc = 0;
    for( Int_t i = 0; i < 100; ++i ){
      if( maxc < xmode_LG->GetBinContent(i) ) maxc = xmode_LG->GetBinContent(i);
    }
    
    /*Float_t pa1m = fitLGT->GetParameter(1);
    Float_t pa2m = fitLGT->GetParameter(2);
    Float_t pa3m = fitLGT->GetParameter(3);
    Float_t epa1m = fitLGT->GetParError(1);
    Float_t epa2m = fitLGT->GetParError(2);
    Float_t epa3m = fitLGT->GetParError(3);
    TString spa1m = to_string( pa1m );
    TString spa2m = to_string( pa2m );
    TString spa3m = to_string( pa3m );
    TString espa1m = to_string( epa1m );
    TString espa2m = to_string( epa2m );
    TString espa3m = to_string( epa3m );
    TString peakm = "Peak = " + spa1m + " +- " + espa1m;
    TString siglm = "SigL = " + spa2m + " +- " + espa2m;
    TString sigrm = "SigR = " + spa3m + " +- " + espa3m;
    TText *tpeakm = new TText( 50, ((float)maxc)*0.9, peakm );
    TText *tsiglm = new TText( 50, ((float)maxc)*0.8, siglm );
    TText *tsigrm = new TText( 50, ((float)maxc)*0.7, sigrm );
    tpeakm->Draw("same");
    tsiglm->Draw("same");
    tsigrm->Draw("same");*/

    c1->Print(pdf_tracking);
    c1->Clear();

    //gStyle->SetOptFit(0);
    // x-dependence of LG PH
    c1->Divide( 2, 2 );
    c1->cd(1);
    x_LG->Draw();
    
    c1->cd(2);
    x_LG->Draw();
    lg_xdep->SetLineWidth(2);
    lg_xdep->Draw("same" );
    
    c1->cd(3);
    lg_xdep->SetLineColor(1);
    TH1F *frameLG = gPad->DrawFrame( -100, 0, 100, 320 );
    frameLG->GetXaxis()->SetTitle( "x (mm)" );
    frameLG->GetYaxis()->SetTitle( "LG Pulse Height (count)" );
    lg_xdep->Draw("same P" );

    c1->cd(4);
    Int_t y_wid = max_lmax - min_lmax;
    if( y_wid == 0 ) y_wid = 1;
    TH1F *frameLGZ = gPad->DrawFrame( -100, min_lmax-y_wid*0.5, 100, max_lmax+y_wid*0.5 );
    frameLGZ->GetXaxis()->SetTitle( "x (mm)" );
    frameLGZ->GetYaxis()->SetTitle( "LG Pulse Height (count)" );
    lg_xdep->Draw("same P" );

    c1->Print(pdf_tracking);
    c1->Clear();

    // fit x-dep of LGPH with linear function
    TGraphErrors *lg_xfit = new TGraphErrors( n_fit, x_pos_fit, ph_lg_fit, ex_pos_fit, eph_lg_fit );
    c1->Divide( 2, 2 );
    c1->cd(1);
    x_LG->Draw();

    c1->cd(3);
    TH1F *frameLG2 = gPad->DrawFrame( -100, 0, 100, 320 );
    frameLG2->GetXaxis()->SetTitle( "x (mm)" );
    frameLG2->GetYaxis()->SetTitle( "LG Pulse Height (count)" );
    lg_xfit->Draw("same P" );
    lg_xfit->Fit("Linelg" );
    Float_t p0 = Linelg->GetParameter(0);
    Float_t p1 = Linelg->GetParameter(1);
    Float_t ep0 = Linelg->GetParError(0);
    Float_t ep1 = Linelg->GetParError(1);
    TString sp0 = to_string( p0 );
    TString sp1 = to_string( p1 );
    TString esp0 = to_string( ep0 );
    TString esp1 = to_string( ep1 );
    TString form = "PH = " + sp1+" x + " + sp0;
    TString eform = "dp1 = " + esp1+", dp0 = " + esp0;
    TText *tform = new TText( -80, 280, form );
    TText *terror = new TText( -80, 250, eform );
    tform->Draw("same");
    terror->Draw("same");

    c1->cd(4);
    y_wid = max_lmax - min_lmax;
    if( y_wid == 0 ) y_wid = 1;
    TH1F *frameLGZ2 = gPad->DrawFrame( -100, min_lmax-y_wid*0.5, 100, max_lmax+y_wid*0.5 );
    frameLGZ2->GetXaxis()->SetTitle( "x (mm)" );
    frameLGZ2->GetYaxis()->SetTitle( "LG Pulse Height (count)" );
    lg_xfit->Draw("same P" );

    c1->Print(pdf_tracking);
    c1->Clear();

    // LG ph not both hits in chamber
    c1->Divide( 2, 2 );
    if( ON_stg ){
      c1->cd(1);
      x25m_j0n1_LG->Draw();
      c1->cd(2);
      y13m_j0n1_LG->Draw();
      c1->cd(3);
      x25m_j1n0_LG->Draw();
      c1->cd(4);
      y13m_j1n0_LG->Draw();
    }else{
      c1->cd(3);
      x25m_j0n1_LG->Draw();
      c1->cd(4);
      y13m_j0n1_LG->Draw();
      c1->cd(1);
      x25m_j1n0_LG->Draw();
      c1->cd(2);
      y13m_j1n0_LG->Draw();
    }

    c1->Print(pdf_tracking);
    c1->Clear();

    // LG ph both hits in chamber
    c1->Divide( 2, 2 );
    if( ON_stg ){
      c1->cd(1);
      x25m_j0w1_LG->Draw();
      c1->cd(2);
      y13m_j0w1_LG->Draw();
      c1->cd(3);
      x25m_j1w0_LG->Draw();
      c1->cd(4);
      y13m_j1w0_LG->Draw();
    }else{
      c1->cd(3);
      x25m_j0w1_LG->Draw();
      c1->cd(4);
      y13m_j0w1_LG->Draw();
      c1->cd(1);
      x25m_j1w0_LG->Draw();
      c1->cd(2);
      y13m_j1w0_LG->Draw();
    }

    c1->Print(pdf_tracking);
    c1->Clear();
    
    // map Energy divided
    c1->Divide( 4,4 );

    for( Int_t i = 0; i < 13; i++ ){
      c1->cd(i+1);
      if( ON_stg )
	xy25m_13mj0_E[i]->Draw();
      else
	xy25m_13mj1_E[i]->Draw();
    }
    c1->cd(14
	   );
      if( ON_stg )
	t0->Draw();
      else
	t1->Draw();
    
    c1->Print(pdf_tracking);
    c1->Clear();
    
    // map Energy divided
    c1->Divide( 4,4 );

    for( Int_t i = 0; i < 13; i++ ){
      c1->cd(i+1);
      if( ON_stg )
	xy25m_13mj0n1_E[i]->Draw();
      else
	xy25m_13mj1n0_E[i]->Draw();
    }
    c1->cd(14);
      if( ON_stg )
	t0->Draw();
      else
	t1->Draw();
    
    c1->Print(pdf_tracking);
    c1->Clear();
    
    // map Energy divided
    c1->Divide( 4,4 );

    for( Int_t i = 0; i < 13; i++ ){
      c1->cd(i+1);
      if( ON_stg )
	xy25m_13mj0w1_E[i]->Draw();
      else
	xy25m_13mj1w0_E[i]->Draw();
    }
    c1->cd(14);
      if( ON_stg )
	t0->Draw();
      else
	t1->Draw();
    
    c1->Print(pdf_tracking);
    c1->Clear();
    
    // map Energy divided
    c1->Divide( 4,4 );

    for( Int_t i = 0; i < 13; i++ ){
      c1->cd(i+1);
      if( ON_stg )
	xy25m_13mj1_E[i]->Draw();
      else
	xy25m_13mj0_E[i]->Draw();
    }
    c1->cd(14);
      if( ON_stg )
	t1->Draw();
      else
	t0->Draw();
    
    c1->Print(pdf_tracking);
    c1->Clear();
    
    // map Energy divided
    c1->Divide( 4,4 );

    for( Int_t i = 0; i < 13; i++ ){
      c1->cd(i+1);
      if( ON_stg )
	xy25m_13mj1n0_E[i]->Draw();
      else
	xy25m_13mj0n1_E[i]->Draw();
    }
    c1->cd(14);
      if( ON_stg )
	t1->Draw();
      else
	t0->Draw();
    
    c1->Print(pdf_tracking);
    c1->Clear();
    
    // map Energy divided
    c1->Divide( 4,4 );

    for( Int_t i = 0; i < 13; i++ ){
      c1->cd(i+1);
      if( ON_stg )
	xy25m_13mj1w0_E[i]->Draw();
      else
	xy25m_13mj0w1_E[i]->Draw();
    }
    c1->cd(14);
      if( ON_stg )
	t1->Draw();
      else
	t0->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();
    
    // map Energy course divided
    c1->Divide( 2, 2 );

    for( Int_t i = 0; i < 4; i++ ){
      c1->cd(i+1);
      if( ON_stg )
	xy25m_13mj0w1_cE[i]->Draw();
      else
	xy25m_13mj1w0_cE[i]->Draw();
    }
    c1->cd(14);
      if( ON_stg )
	t0->Draw();
      else
	t1->Draw();
    
    c1->Print(pdf_tracking);
    c1->Clear();

    // map Energy course divided
    c1->Divide( 2,2 );

    for( Int_t i = 0; i < 4; i++ ){
      c1->cd(i+1);
      if( ON_stg )
	xy25m_13mj1w0_cE[i]->Draw();
      else
	xy25m_13mj0w1_cE[i]->Draw();
    }
    c1->cd(4);
    if( ON_stg )
      t1->Draw();
    else
      t0->Draw();
    

    c1->Print(pdf_tracking);
    c1->Clear();

    // ToT vs LG
    c1->Divide( 2,2 );

    for( Int_t i = 0; i < 4; i++ ){
      c1->cd(i+1);
      totvsPH[i]->Draw();
    }

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );

    c1->cd(1);
    LGlength->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );

    c1->cd(1);
    LGintlength->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 2,2 );

    c1->cd(1);
    LGintlength->Draw();

    c1->cd(3);
    LGintlength1->Draw();

    c1->cd(4);
    LGintlength2->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );

    c1->cd(1);
    intLG->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );

    c1->cd(1);
    gPad->SetLogy(1);
    intLG->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );

    c1->cd(1);
    intLGth->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );

    c1->cd(1);
    gPad->SetLogy(1);
    intLGth->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );

    c1->cd(1);
    intLGpos->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );

    c1->cd(1);
    gPad->SetLogy(1);
    intLGpos->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 2,2 );

    c1->cd(1);
    gPad->SetLogy(1);
    intLGpos->Draw();

    c1->cd(3);
    gPad->SetLogy(1);
    intLGpos1->Draw();

    c1->cd(4);
    gPad->SetLogy(1);
    intLGpos2->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );

    c1->cd(1);
    gPad->SetLogy(1);
    fitintLG->SetParameter( 7, 2000 );// n0
    fitintLG->SetParameter( 8, -0.007 );// slope
    fitintLG->SetParameter( 0, 40000 );// h0
    fitintLG->SetParameter( 4, 120 );// h1
    fitintLG->SetParameter( 2, 20 );// s0l
    fitintLG->SetParameter( 3, 20 );// s0r
    fitintLG->SetParameter( 6, 15 );// s1
    if( p10 == 30 ){
      fitintLG->SetParameter( 1, 350 );// c0
      fitintLG->SetParameter( 5, 650 );// c1
    }
    if( p10 == 20 ){
      fitintLG->SetParameter( 1, 250 );
      fitintLG->SetParameter( 5, 450 );
    }
    if( p10 <= 10 ){
      fitintLG->SetParameter( 1, 150 );
      fitintLG->SetParameter( 5, 200 );
    }
    if( p10 >= 50 ){
      fitintLG->SetParameter( 1, 550 );
      fitintLG->SetParameter( 5, 900 );
    }
    if( p10 == 40 ){
      fitintLG->SetParameter( 1, 450 );
      fitintLG->SetParameter( 5, 800 );
    }

    intLGpos->Fit( "fitintLG", "", "", 50, 1200 );

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );

    c1->cd(1);
    x_intLGpos->Draw();

    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );

    TH1D *x_intLGpos_slc = x_intLGpos->ProjectionY( "25m_13m_int_-10_10", 46, 55 );
    
    x_intLGpos_slc->Draw();
    
    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );
    gPad->SetLogy(1);
    x_intLGpos_slc->Draw();
    
    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );
    gPad->SetLogy(1);
    fitintLG->SetParameter( 0, 10000 );
    fitintLG->SetParameter( 4, 20 );
    fitintLG->SetParameter( 2, 20 );
    fitintLG->SetParameter( 3, 20 );
    fitintLG->SetParameter( 6, 10 );
    if( p10 == 30 ){
      fitintLG->SetParameter( 1, 350 );
      fitintLG->SetParameter( 5, 650 );
      fitintLG->SetParameter( 7, 100 );
      slope = 100/350;
      fitintLG->SetParameter( 8, slope );
    }
    if( p10 == 20 ){
      fitintLG->SetParameter( 1, 250 );
      fitintLG->SetParameter( 5, 450 );
      fitintLG->SetParameter( 7, 200 );
      fitintLG->SetParameter( 8, -1. );
      slope = 200/250;
      fitintLG->SetParameter( 8, slope );
    }
    if( p10 <= 10 ){
      fitintLG->SetParameter( 1, 110 );
      fitintLG->SetParameter( 5, 220 );
      fitintLG->SetParameter( 7, 400 );
      fitintLG->SetParameter( 8, -2. );
      slope = 400/110;
      fitintLG->SetParameter( 8, slope );
    }
    if( p10 >= 50 ){
      fitintLG->SetParameter( 1, 550 );
      fitintLG->SetParameter( 5, 850 );
      fitintLG->SetParameter( 7, 40 );
      fitintLG->SetParameter( 8, -0.08 );
      slope = 40/550;
      fitintLG->SetParameter( 8, slope );
    }
    if( p10 == 40 ){
      fitintLG->SetParameter( 1, 450 );
      fitintLG->SetParameter( 5, 800 );
      fitintLG->SetParameter( 7, 100 );
      fitintLG->SetParameter( 8, -0.1 );
      slope = 100/450;
      fitintLG->SetParameter( 8, slope );
    }

    x_intLGpos_slc->Fit( "fitintLG", "", "", 50, 1200 );
    
    c1->Print(pdf_tracking);
    c1->Clear();

    c1->Divide( 1,1 );
    gPad->SetLogy(0);
    xw_intLGpos->Draw();
    
}// LG_on

  pdf_tracking = "/group/itdc/tbl/2024ARTBL018/pdf/2D/tracking_run" + runnumber + ".pdf)";
  c1->Print(pdf_tracking);
  TString outroot_tracking = "output.root";
  Float_t SF = 1.;
  if(runnumber == "4161"){
	  SF = 4878.42/200000;
	  outroot_tracking = "Nov_QAFOFF.root";
	  pdf_tracking = "Nov_QAFOFF.pdf";
  }else if(runnumber == "4180"){
	  SF = 7412.87/200000;
	  outroot_tracking = "Nov_QAFON.root";
	  pdf_tracking = "Nov_QAFON.pdf";
  }else if(runnumber == "3305"){
	  SF = 4878.42/300000;
	  outroot_tracking = "3p0GeV_GapOpen.root";
	  pdf_tracking = "3p0GeV_GapOpen.pdf";
	  //outroot_tracking = "3p0GeV_GapOpen.root";
	  //pdf_tracking = "3p0GeV_GapOpen.pdf";
  }else if(runnumber == "3306"){
	  SF = 12521.9/300000;
	  outroot_tracking = "2p5GeV_GapOpen.root";
	  pdf_tracking = "2p5GeV_GapOpen.pdf";
  }else if(runnumber == "3307"){
	  SF = 13308.8/300000;
	  outroot_tracking = "2p0GeV_GapOpen.root";
	  pdf_tracking = "2p0GeV_GapOpen.pdf";
  }else if(runnumber == "3308"){
	  SF = 13010.8/300000;
	  outroot_tracking = "1p5GeV_GapOpen.root";
	  pdf_tracking = "1p5GeV_GapOpen.pdf";
  }else if(runnumber == "3309"){
	  SF = 10652.6/300000;
	  outroot_tracking = "1p0GeV_GapOpen.root";
	  pdf_tracking = "1p0GeV_GapOpen.pdf";
  }else if(runnumber == "3310"){
	  SF = 5229.9/300000;
	  outroot_tracking = "0p5GeV_GapOpen.root";
	  pdf_tracking = "0p5GeV_GapOpen.pdf";
  }else if(runnumber == "3311"){
	  SF = 3867.0/300000;
	  outroot_tracking = "0p4GeV_GapOpen.root";
	  pdf_tracking = "0p4GeV_GapOpen.pdf";
  }else if(runnumber == "3312"){
	  SF = 2524.7/300000;
	  outroot_tracking = "0p3GeV_GapOpen.root";
	  pdf_tracking = "0p3GeV_GapOpen.pdf";
  }else if(runnumber == "3313"){
	  SF = 1284.0/300000;
	  outroot_tracking = "0p2GeV_GapOpen.root";
	  pdf_tracking = "0p2GeV_GapOpen.pdf";
  }else if(runnumber == "3314"){
	  SF = 359.3/100000;
	  outroot_tracking = "0p1GeV_GapOpen.root";
	  pdf_tracking = "0p1GeV_GapOpen.pdf";
  }else if(runnumber == "3316"){
	  SF = 1286.1/300000;
	  outroot_tracking = "5p8GeV_GapOpen.root";
	  pdf_tracking = "5p8GeV_GapOpen.pdf";
  }else if(runnumber == "3318"){
	  SF = 2008.0/300000;
	  outroot_tracking = "5p5GeV_GapOpen.root";
	  pdf_tracking = "5p5GeV_GapOpen.pdf";
  }else if(runnumber == "3319"){
	  SF = 3422.0/300000;
	  outroot_tracking = "5p0GeV_GapOpen.root";
	  pdf_tracking = "5p0GeV_GapOpen.pdf";
  }else if(runnumber == "3320"){
	  SF = 5339.5/300000;
	  outroot_tracking = "4p5GeV_GapOpen.root";
	  pdf_tracking = "4p5GeV_GapOpen.pdf";
  }else if(runnumber == "3321"){
	  SF = 7269.0/300000;
	  outroot_tracking = "4p0GeV_GapOpen.root";
	  pdf_tracking = "4p0GeV_GapOpen.pdf";
  }else if(runnumber == "3322"){
	  SF = 9043.4/300000;
	  outroot_tracking = "3p5GeV_GapOpen.root";
	  pdf_tracking = "3p5GeV_GapOpen.pdf";
  }
  TFile *fout = new TFile(outroot_tracking,"RECREATE");
  xy25m_13mj0->Scale(SF);
  xy25m_13mj1->Scale(SF);
  LG_converted->Scale(SF);
  xdx25m_13mj1_j0->Scale(SF);
  ydx25m_13mj1_j0->Scale(SF);
  xdy25m_13mj1_j0->Scale(SF);
  ydy25m_13mj1_j0->Scale(SF);
  xy25m_13mj0->GetXaxis()->SetTitle("x [mm]");
  xy25m_13mj0->GetYaxis()->SetTitle("y [mm]");
  xy25m_13mj0->GetZaxis()->SetTitle("Hz");
  xy25m_13mj1->GetXaxis()->SetTitle("x [mm]");
  xy25m_13mj1->GetYaxis()->SetTitle("y [mm]");
  xy25m_13mj1->GetZaxis()->SetTitle("Hz");
  LG_converted->GetXaxis()->SetTitle("GeV");
  LG_converted->GetYaxis()->SetTitle("Hz");
  xdx25m_13mj1_j0->GetXaxis()->SetTitle("x [mm] at 180 mm downstream from the magnet");
  xdx25m_13mj1_j0->GetYaxis()->SetTitle("Difference in x [mm] (4.7 m along the beam axis)");
  xdx25m_13mj1_j0->GetZaxis()->SetTitle("Hz");
  ydx25m_13mj1_j0->GetXaxis()->SetTitle("y [mm] at 180 mm downstream from the magnet");
  ydx25m_13mj1_j0->GetYaxis()->SetTitle("Difference in x [mm] (4.7 m along the beam axis)");
  ydx25m_13mj1_j0->GetZaxis()->SetTitle("Hz");
  xdy25m_13mj1_j0->GetXaxis()->SetTitle("x [mm] at 180 mm downstream from the magnet");
  xdy25m_13mj1_j0->GetYaxis()->SetTitle("Difference in y [mm] (4.7 m along the beam axis)");
  xdy25m_13mj1_j0->GetZaxis()->SetTitle("Hz");
  ydy25m_13mj1_j0->GetXaxis()->SetTitle("y [mm] at 180 mm downstream from the magnet");
  ydy25m_13mj1_j0->GetYaxis()->SetTitle("Difference in y [mm] (4.7 m along the beam axis)");
  ydy25m_13mj1_j0->GetZaxis()->SetTitle("Hz");
  xy25m_13mj0->SetTitle("x:y at 180 mm downstream from the magnet");
  xy25m_13mj1->SetTitle("x:y at 4880 mm downstream from the magnet");
  LG_converted->SetTitle("Momentum");
  xdx25m_13mj1_j0->SetTitle("x:dx");
  xdy25m_13mj1_j0->SetTitle("x:dy");
  ydx25m_13mj1_j0->SetTitle("y:dx");
  ydy25m_13mj1_j0->SetTitle("y:dy");
  xy25m_13mj0->Write();
  xy25m_13mj1->Write();
  LG_converted->Write();
  xdx25m_13mj1_j0->Write();
  ydx25m_13mj1_j0->Write();
  xdy25m_13mj1_j0->Write();
  ydy25m_13mj1_j0->Write();
  TCanvas *cout = new TCanvas();
  cout->Divide( 2,2 );
  cout->cd(1)->SetRightMargin(0.2);
  xy25m_13mj0->SetStats(0);
  xy25m_13mj0->Draw("COLZ");
  cout->cd(2)->SetRightMargin(0.2);
  xy25m_13mj1->SetStats(0);
  xy25m_13mj1->Draw("COLZ");
  cout->cd(3);
  LG_converted->SetStats(0);
  LG_converted->Draw("hist");
  cout->cd(4)->SetRightMargin(0.2);
  xdx25m_13mj1_j0->SetStats(0);
  xdx25m_13mj1_j0->Draw("COLZ");
  xy25m_13mj0w1->Write();
  xy25m_13mj1w0->Write();
  for(int ii = 0; ii < 13; ii++){
	  xy25m_13mj0w1_E[ii]->Write();
	  xy25m_13mj1w0_E[ii]->Write();
  }
//cout->Print(pdf_tracking);

  fout->Close();

}
