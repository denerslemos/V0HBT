#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include <TRandom1.h>
#include <vector>
#include <TLorentzVector.h>
#include "THnSparse.h"
#include "TRandom3.h"
#include <cstring>
#include <ctime>
#include <ctime>
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <vector>
#include <map>
#include "TFrame.h"
#include "TH1F.h"
#include "TBenchmark.h"
#include "TSystem.h"
#include "THnBase.h"


#define PI ( 4.*atan(1.) )


//ctes
#define k0s_mass 0.497614
#define la_mass 1.1156
#define pimass 0.1396
#define pro_mass 0.93827
#define electronMass 0.000511
#define kaon_mass 0.493677
#define xi_mass 1.32171
#define om_mass 1.67245

using namespace std;


//constants
 
//kinematics

double ptminK0s=0.3, ptmaxK0s=8.5, massK0s=0.4976, sigmaK0s=0.006;
double ptminLAL=0.5, ptmaxLAL=8.5, massLAL=1.1156, sigmaLAL=0.0025;
double ptminXAX=0.5, ptmaxXAX=8.5, massXAX=1.322, sigmaXAX=0.005;
double ptminOAO=0.5, ptmaxOAO=8.5, massOAO=1.672, sigmaOAO=0.005;
//double ptminJet=0.0, ptmaxJet=10.0; // for low pT jet
//double ptminJet=20.0, ptmaxJet=200.0; // for high pT jet
//double ptminJet=10.0, ptmaxJet=200.0; // for all jets
double ptminJet=0.0, ptmaxJet=200.0; // for all jets

double ctaucutmin = 0.0, ctaucutmax = 100000000.;
//double ctaucutmin = 0.0, ctaucutmax = 10.;
//double ctaucutmin = 0.0, ctaucutmax = 6.;
//double ctaucutmin = 0.0, ctaucutmax = 15.;

//rapidity range
double rapvar = 50.0; //probably i'll not cut on rapidity for pPb

 //for MC matching
double drdpt = 0.03; //do not change anymore
 
//declare variables

//cut variables
//general
double vtxzcutmax, vtxzcutmin, vtxrhocut;
double nsigmapeak, nsigmaside;
double Dchi2;
double pixelhits;
double nhitss;

//misidentified
double misK0s, misLam, misee; 

//V0s
double dxyz, cospoint, trkdca, decaylength; 


//Xi
double Xi_VTrkP3DIpSigValue_val = 3.; //3D impact parameter significance of proton track from Λ decay with respect to the primary vertex (>) sys (4)
double Xi_VTrkPi3DIpSigValue_val = 4.; //3D impact parameter significance of π track from Λ decay with respect to the primary vertex(>) sys (5)
double Xi_casPi3DIpSigValue_val = 5.; //3D impact parameter significance of π (K) track from Ξ− (Ω−) decay with respect to the primary vertex(>) sys (3)
double Xi_cas3DIpSigValue_val = 2.5; //3D impact parameter significance of the Ξ− (Ω−) candidate with respect to the primary vertex (<) sys (1)
double Xi_distanceSigValue_val = 12.; //3D separation significance between Λ vertex and primary vertex (>) sys (6)
double Xi_casFlightSigValue_val = 3.; //3D separation significance between Ξ− (Ω−) vertex and primary vertex (>) sys (2)

//Om
double Om_VTrkP3DIpSigValue_val = 2.; //3D impact parameter significance of proton track from Λ decay with respect to the primary vertex (>) sys (4)
double Om_VTrkPi3DIpSigValue_val = 3.; //3D impact parameter significance of π track from Λ decay with respect to the primary vertex(>) sys (5)
double Om_casPi3DIpSigValue_val = 4.; //3D impact parameter significance of π (K) track from Ξ− (Ω−) decay with respect to the primary vertex(>) sys (3)
double Om_cas3DIpSigValue_val = 3.0; //3D impact parameter significance of the Ξ− (Ω−) candidate with respect to the primary vertex (<) sys (1)
double Om_distanceSigValue_val = 10.; //3D separation significance between Λ vertex and primary vertex (>) sys (6)
double Om_casFlightSigValue_val = 2.; //3D separation significance between Ξ− (Ω−) vertex and primary vertex (>) sys (2)

double misXi = 0.015, misOm = 0.015;



//systematic for pPb - double check

TH2D* effhisto_K0s = nullptr;
TH2D* effhisto_Lam = nullptr;
TH2D* effhisto_ALam = nullptr;
TH2D* effhisto_LAL = nullptr;
TH1D* Vz_w_hist  = nullptr;
TH1D* Ntk_w_hist  = nullptr;
TH2D* NtkVz_w_hist  = nullptr;

string get_sys(int p){

string syst;

//related to corr

if(p==0){ //standard
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0std";
} else if (p==1){ //tight
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.2;
cospoint = 0.9995;
trkdca = 0.5;
decaylength = 6.5;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0tight";
} else if (p==2){ //loose
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.997;
trkdca = 0.5;
decaylength = 3.5;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0loose";
} else if (p==7){ //dca plus
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 1.0;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0dcaplus";
} else if (p==8){ //dca minus
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.3;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0dcaminus";
} else if (p==9){ //misidentified plus
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.012;
misLam = 0.022;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0misplus";
} else if (p==10){ //misidentified minus
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.007;
misLam = 0.015;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0misminus";
} else if (p==11){ //misidentified ee plus
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.020;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0miseeplus";
} else if (p==12){ //misidentified ee minus
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.010;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0miseeminus";
} else if (p==13){ //rho vertex plus
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.20; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0rhoplus";
} else if (p==14){ //z vertex 3 to 15
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0zvxt3to15";
} else if (p==15){ //z vertex 3
//general
vtxzcutmax = 3.0;
vtxzcutmin = -3.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0zvxt3";
} else if (p==16){ //sigma peak minus (1.5)
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 1.8;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0sigmapeakminus";
} else if (p==17){ //sigma peak plus (2.5)
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.2;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0sigmapeakplus";
} else if (p==18){ //sigma side plus (3.5)
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.2;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0sigmasideplus";
} else if (p==19){ //sigma side minus (2.5)
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 2.8;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0sigmasideminus";
} else if (p==20){//EPOS eff corr
TString effile;
effile = "/eos/cms/store/group/phys_heavyions/ddesouza/DATA/EPOS_Eff_V0.root";
TFile *fileeff = TFile::Open(Form("%s",effile.Data()));
if(!fileeff->IsOpen()){cout << "cannot find the file of efficiency" << endl; }else{cout << "eff file founded" << endl;}
fileeff->GetObject("K0s_eff", effhisto_K0s);
if(!effhisto_K0s){cout << "cannot find K0s efficiency histogram" << endl; }else{cout << "K0s eff hist founded" << endl;}
fileeff->GetObject("Lam_eff", effhisto_Lam);
if(!effhisto_Lam){cout << "cannot find Lam efficiency histogram" << endl; }else{cout << "Lam eff hist founded" << endl;}
fileeff->GetObject("ALam_eff", effhisto_ALam);
if(!effhisto_ALam){cout << "cannot find ALam efficiency histogram" << endl; }else{cout << "ALam eff hist founded" << endl;}
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0epos";
} else if (p==21){//HIJING eff corr
TString effile;
effile = "/afs/cern.ch/work/d/ddesouza/Eff/Eff003/HIJING_Eff_V0.root";
TFile *fileeff = TFile::Open(Form("%s",effile.Data()));
if(!fileeff->IsOpen()){cout << "cannot find the file of efficiency" << endl; }else{cout << "eff file founded" << endl;}
fileeff->GetObject("K0s_eff", effhisto_K0s);
if(!effhisto_K0s){cout << "cannot find K0s efficiency histogram" << endl; }else{cout << "K0s eff hist founded" << endl;}
fileeff->GetObject("Lam_eff", effhisto_Lam);
if(!effhisto_Lam){cout << "cannot find Lam efficiency histogram" << endl; }else{cout << "Lam eff hist founded" << endl;}
fileeff->GetObject("ALam_eff", effhisto_ALam);
if(!effhisto_ALam){cout << "cannot find ALam efficiency histogram" << endl; }else{cout << "ALam eff hist founded" << endl;}
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0hijing";
} else if(p==23){ //PU pileUpFilter_pPb8TeV_Gplus
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0Gplus";
} else if(p==24){ //PU pileUpFilter_pPb8TeV_vtx1
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0vtx1";
} else if(p==25){ //no PU filter
//general
vtxzcutmax = 15.0;
vtxzcutmin = -15.0;
vtxrhocut = 0.15; //0.20 (sys)
nsigmapeak = 2.0;
nsigmaside = 3.0;
// misidenfied
misK0s = 0.009;
misLam = 0.018;
misee = 0.015;
// V0 selection
dxyz = 1.;
cospoint = 0.999;
trkdca = 0.5;
decaylength = 5.0;
//chi2
Dchi2 = 0.000001;
pixelhits = 0;
nhitss = 4;
syst = "V0noPU";
}

return syst;

} //DONE
//----------------------------------------------------------------------------------------------------------------------


//identical particles

//for MC weight
std::vector<double> weight_K0sK0s;
std::vector<double> weight_LL;
std::vector<double> weight_ALAL;
std::vector<double> weight_LAL;
std::vector<double> weight_K0sL;
std::vector<double> weight_K0sAL;

//vectors for mixing

//K0s - 2 per event

std::vector<double> ev_z_vtxK0s;
std::vector<int> ev_ntrkoff_vecK0s;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0s;
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vecK0s;
std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vecK0s;

//Lam - 2 per event

std::vector<double> ev_z_vtxLam;
std::vector<int> ev_ntrkoff_vecLam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLam;
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vecLam;
std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vecLam;

//ALam - 2 per event

std::vector<double> ev_z_vtxALam;
std::vector<int> ev_ntrkoff_vecALam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALam;
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vecALam;
std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vecALam;


//non-identical femtoscopy

//Lam-ALam
std::vector<double> ev_z_vtxLamALam;
std::vector<int> ev_ntrkoff_vecLamALam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamALam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamLam;
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vecLamALam;
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vecALamLam;
std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vecLamALam;

//K0s-Lam
std::vector<double> ev_z_vtxK0sLam;
std::vector<int> ev_ntrkoff_vecK0sLam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0sLam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamK0s;
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vecK0sLam;
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vecLamK0s;
std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vecK0sLam;

//K0s-ALam
std::vector<double> ev_z_vtxK0sALam;
std::vector<int> ev_ntrkoff_vecK0sALam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0sALam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamK0s;
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vecK0sALam;
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vecALamK0s;
std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vecK0sALam;

//for feed down

std::vector<double> weight_LLFD;
std::vector<double> ev_z_vtx_Lam_LamFD;
std::vector<int> ev_ntrkoff_vec_Lam_LamFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_Lam_LamFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_LamFD_Lam;

std::vector<double> weight_LFDLFD;
std::vector<double> ev_z_vtx_LamFD_LamFD;
std::vector<int> ev_ntrkoff_vec_LamFD_LamFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_LamFD_LamFD;

std::vector<double> weight_ALALFD;
std::vector<double> ev_z_vtx_ALam_ALamFD;
std::vector<int> ev_ntrkoff_vec_ALam_ALamFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_ALam_ALamFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_ALamFD_ALam;

std::vector<double> weight_ALFDALFD;
std::vector<double> ev_z_vtx_ALamFD_ALamFD;
std::vector<int> ev_ntrkoff_vec_ALamFD_ALamFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_ALamFD_ALamFD;

std::vector<double> weight_LFDALFD;
std::vector<double> ev_z_vtx_LFDALFD;
std::vector<int> ev_ntrkoff_vec_LFDALFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_LFDALFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_ALFDLFD;

std::vector<double> weight_LALFD;
std::vector<double> ev_z_vtx_LALFD;
std::vector<int> ev_ntrkoff_vec_LALFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_LALFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_LFDAL;

std::vector<double> weight_ALLFD;
std::vector<double> ev_z_vtx_ALLFD;
std::vector<int> ev_ntrkoff_vec_ALLFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_ALLFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_ALFDL;

std::vector<double> weight_KLFD;
std::vector<double> ev_z_vtx_KLFD;
std::vector<int> ev_ntrkoff_vec_KLFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_KLFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_LFDK;

std::vector<double> weight_KALFD;
std::vector<double> ev_z_vtx_KALFD;
std::vector<int> ev_ntrkoff_vec_KALFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_KALFD;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_ALFDK;

//for shape studies

//K0sK0s

std::vector<double> ev_z_vtxK0s_TM;
std::vector<int> ev_ntrkoff_vecK0s_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0s_TM;

std::vector<double> ev_z_vtxK0s_TU;
std::vector<int> ev_ntrkoff_vecK0s_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0s_TU;

std::vector<double> ev_z_vtxK0s_TM_TU;
std::vector<int> ev_ntrkoff_vecK0s_TM_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0s_TM_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0s_TU_TM;

//LamLam

std::vector<double> ev_z_vtxLam_TM;
std::vector<int> ev_ntrkoff_vecLam_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLam_TM;

std::vector<double> ev_z_vtxLam_TU;
std::vector<int> ev_ntrkoff_vecLam_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLam_TU;

std::vector<double> ev_z_vtxLam_TM_TU;
std::vector<int> ev_ntrkoff_vecLam_TM_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLam_TM_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLam_TU_TM;

//ALamALam

std::vector<double> ev_z_vtxALam_TM;
std::vector<int> ev_ntrkoff_vecALam_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALam_TM;

std::vector<double> ev_z_vtxALam_TU;
std::vector<int> ev_ntrkoff_vecALam_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALam_TU;

std::vector<double> ev_z_vtxALam_TM_TU;
std::vector<int> ev_ntrkoff_vecALam_TM_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALam_TM_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALam_TU_TM;

//LAL

std::vector<double> ev_z_vtxLAL_TM;
std::vector<int> ev_ntrkoff_vecLAL_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLAL_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALL_TM;

std::vector<double> ev_z_vtxLAL_TU;
std::vector<int> ev_ntrkoff_vecLAL_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLAL_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALL_TU;

std::vector<double> ev_z_vtxLAL_TM_TU;
std::vector<int> ev_ntrkoff_vecLAL_TM_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLAL_TM_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALL_TM_TU;

std::vector<double> ev_z_vtxLAL_TU_TM;
std::vector<int> ev_ntrkoff_vecLAL_TU_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLAL_TU_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALL_TU_TM;

//KL

std::vector<double> ev_z_vtxKL_TM;
std::vector<int> ev_ntrkoff_vecKL_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecKL_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLK_TM;

std::vector<double> ev_z_vtxKL_TU;
std::vector<int> ev_ntrkoff_vecKL_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecKL_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLK_TU;

std::vector<double> ev_z_vtxKL_TM_TU;
std::vector<int> ev_ntrkoff_vecKL_TM_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecKL_TM_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLK_TM_TU;

std::vector<double> ev_z_vtxKL_TU_TM;
std::vector<int> ev_ntrkoff_vecKL_TU_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecKL_TU_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLK_TU_TM;

//KAL

std::vector<double> ev_z_vtxKAL_TM;
std::vector<int> ev_ntrkoff_vecKAL_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecKAL_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALK_TM;

std::vector<double> ev_z_vtxKAL_TU;
std::vector<int> ev_ntrkoff_vecKAL_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecKAL_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALK_TU;

std::vector<double> ev_z_vtxKAL_TM_TU;
std::vector<int> ev_ntrkoff_vecKAL_TM_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecKAL_TM_TU;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALK_TM_TU;

std::vector<double> ev_z_vtxKAL_TU_TM;
std::vector<int> ev_ntrkoff_vecKAL_TU_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecKAL_TU_TM;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALK_TU_TM;


//working with jets

//pure jet

std::vector<double> ev_z_vtxJet;
std::vector<int> ev_ntrkoff_vecJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecJet;

//new V0Jet

std::vector<double> ev_z_vtxK0sJet;
std::vector<int> ev_ntrkoff_vecK0sJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0sJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecJetK0s;

std::vector<double> ev_z_vtxLamJet;
std::vector<int> ev_ntrkoff_vecLamJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecJetLam;

std::vector<double> ev_z_vtxALamJet;
std::vector<int> ev_ntrkoff_vecALamJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecJetALam;


std::vector<double> ev_z_vtxK0snojetJet;
std::vector<int> ev_ntrkoff_vecK0snojetJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0snojetJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecJetK0snojet;

std::vector<double> ev_z_vtxLamnojetJet;
std::vector<int> ev_ntrkoff_vecLamnojetJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamnojetJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecJetLamnojet;

std::vector<double> ev_z_vtxALamnojetJet;
std::vector<int> ev_ntrkoff_vecALamnojetJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamnojetJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecJetALamnojet;


std::vector<double> ev_z_vtxK0sfromjetJet;
std::vector<int> ev_ntrkoff_vecK0sfromjetJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0sfromjetJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecJetK0sfromjet;

std::vector<double> ev_z_vtxLamfromjetJet;
std::vector<int> ev_ntrkoff_vecLamfromjetJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamfromjetJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecJetLamfromjet;

std::vector<double> ev_z_vtxALamfromjetJet;
std::vector<int> ev_ntrkoff_vecALamfromjetJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamfromjetJet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecJetALamfromjet;



//no jet

std::vector<double> ev_z_vtxK0s_nojet;
std::vector<int> ev_ntrkoff_vecK0s_nojet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0s_nojet;

std::vector<double> ev_z_vtxLam_nojet;
std::vector<int> ev_ntrkoff_vecLam_nojet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLam_nojet;

std::vector<double> ev_z_vtxALam_nojet;
std::vector<int> ev_ntrkoff_vecALam_nojet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALam_nojet;

std::vector<double> ev_z_vtxLamALam_nojet;
std::vector<int> ev_ntrkoff_vecLamALam_nojet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamALam_nojet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamLam_nojet;

std::vector<double> ev_z_vtxK0sLam_nojet;
std::vector<int> ev_ntrkoff_vecK0sLam_nojet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0sLam_nojet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamK0s_nojet;

std::vector<double> ev_z_vtxK0sALam_nojet;
std::vector<int> ev_ntrkoff_vecK0sALam_nojet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0sALam_nojet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamK0s_nojet;

//1 from jet and 1 not

std::vector<double> ev_z_vtxK0s_only1jet;
std::vector<int> ev_ntrkoff_vecK0s_only1jet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0s_only1jet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0s_only1jetX;

std::vector<double> ev_z_vtxLam_only1jet;
std::vector<int> ev_ntrkoff_vecLam_only1jet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLam_only1jet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLam_only1jetX;

std::vector<double> ev_z_vtxALam_only1jet;
std::vector<int> ev_ntrkoff_vecALam_only1jet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALam_only1jet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALam_only1jetX;

std::vector<double> ev_z_vtxLamALam_only1jet;
std::vector<int> ev_ntrkoff_vecLamALam_only1jet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamALam_only1jet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamLam_only1jet;

std::vector<double> ev_z_vtxLamALam_only1jetX;
std::vector<int> ev_ntrkoff_vecLamALam_only1jetX;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamALam_only1jetX;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamLam_only1jetX;


std::vector<double> ev_z_vtxK0sLam_only1jet;
std::vector<int> ev_ntrkoff_vecK0sLam_only1jet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0sLam_only1jet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamK0s_only1jet;

std::vector<double> ev_z_vtxK0sLam_only1jetX;
std::vector<int> ev_ntrkoff_vecK0sLam_only1jetX;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0sLam_only1jetX;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamK0s_only1jetX;


std::vector<double> ev_z_vtxK0sALam_only1jet;
std::vector<int> ev_ntrkoff_vecK0sALam_only1jet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0sALam_only1jet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamK0s_only1jet;

std::vector<double> ev_z_vtxK0sALam_only1jetX;
std::vector<int> ev_ntrkoff_vecK0sALam_only1jetX;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0sALam_only1jetX;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamK0s_only1jetX;

//for same or different jetstudies

std::vector<double> ev_z_vtxK0s_samediffjet;
std::vector<int> ev_ntrkoff_vecK0s_samediffjet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0s_samediffjet;
std::vector<std::vector<Int_t>> ev_njet_vecK0s_samediffjet;


std::vector<double> ev_z_vtxLam_samediffjet;
std::vector<int> ev_ntrkoff_vecLam_samediffjet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLam_samediffjet;
std::vector<std::vector<Int_t>> ev_njet_vecLam_samediffjet;

std::vector<double> ev_z_vtxALam_samediffjet;
std::vector<int> ev_ntrkoff_vecALam_samediffjet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALam_samediffjet;
std::vector<std::vector<Int_t>> ev_njet_vecALam_samediffjet;

std::vector<double> ev_z_vtxLamALam_samediffjet;
std::vector<int> ev_ntrkoff_vecLamALam_samediffjet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamALam_samediffjet;
std::vector<std::vector<Int_t>> ev_njet_vecLamALam_samediffjet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamLam_samediffjet;
std::vector<std::vector<Int_t>> ev_njet_vecALamLam_samediffjet;

std::vector<double> ev_z_vtxK0sLam_samediffjet;
std::vector<int> ev_ntrkoff_vecK0sLam_samediffjet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0sLam_samediffjet;
std::vector<std::vector<Int_t>> ev_njet_vecK0sLam_samediffjet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecLamK0s_samediffjet;
std::vector<std::vector<Int_t>> ev_njet_vecLamK0s_samediffjet;

std::vector<double> ev_z_vtxK0sALam_samediffjet;
std::vector<int> ev_ntrkoff_vecK0sALam_samediffjet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecK0sALam_samediffjet;
std::vector<std::vector<Int_t>> ev_njet_vecK0sALam_samediffjet;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecALamK0s_samediffjet;
std::vector<std::vector<Int_t>> ev_njet_vecALamK0s_samediffjet;

std::vector<double> ev_z_vtxK0sGen;
std::vector<int> ev_ntrkoff_vecGenK0s;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenK0s;

std::vector<double> ev_z_vtxLamGen;
std::vector<int> ev_ntrkoff_vecGenLam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenLam;

std::vector<double> ev_z_vtxALamGen;
std::vector<int> ev_ntrkoff_vecGenALam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenALam;

std::vector<double> ev_z_vtxLamALamGen;
std::vector<int> ev_ntrkoff_vecGenLamALam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenLamALam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenALamLam;

std::vector<double> ev_z_vtxK0sLamGen;
std::vector<int> ev_ntrkoff_vecGenK0sLam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenK0sLam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenLamK0s;

std::vector<double> ev_z_vtxK0sALamGen;
std::vector<int> ev_ntrkoff_vecGenK0sALam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenK0sALam;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenALamK0s;


//true

std::vector<double> ev_z_vtxK0sGenT;
std::vector<int> ev_ntrkoff_vecGenK0sT;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenK0sT;

std::vector<double> ev_z_vtxLamGenT;
std::vector<int> ev_ntrkoff_vecGenLamT;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenLamT;

std::vector<double> ev_z_vtxALamGenT;
std::vector<int> ev_ntrkoff_vecGenALamT;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenALamT;

std::vector<double> ev_z_vtxLamALamGenT;
std::vector<int> ev_ntrkoff_vecGenLamALamT;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenLamALamT;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenALamLamT;

std::vector<double> ev_z_vtxK0sLamGenT;
std::vector<int> ev_ntrkoff_vecGenK0sLamT;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenK0sLamT;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenLamK0sT;

std::vector<double> ev_z_vtxK0sALamGenT;
std::vector<int> ev_ntrkoff_vecGenK0sALamT;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenK0sALamT;
std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vecGenALamK0sT;


//--------------------------------------------------------------------------------------------------------------


//tree variables

Int_t N_tkoff;
Int_t N_tkoffGen;
Double_t N_tkoffGenwc;
Int_t pileUpFilter_pPb8TeV_vtx1;
Int_t pileUpFilter_pPb8TeV_Gplus;
Int_t olvFilter_pPb8TeV_dz1p0;


//Event selection filters
Int_t phffilter_1tw;
Int_t pprimaryvertexfilter;
Int_t pnoscrampingfilter;

//Event selection filters
Int_t validPV;
Int_t tracksizePV;
Int_t minnTowersTh3HF;
Int_t trkColl_noscr;
Double_t fraction_noscr;

Double_t vtx_z;
Double_t vtx_rho;
Double_t etHFtowerSum;
Double_t Eta_weight;
Double_t Eta_weightGen;
Double_t Eta_weightGenwc;


//---------------------------  K0s  -------------------------
//daughter 1 (pi+)
std::vector<double> *K0s_dxy1 = 0; 
std::vector<double> *K0s_dz1 = 0; 
std::vector<double> *K0s_chi21 = 0; 
std::vector<double> *K0s_d1px = 0; 
std::vector<double> *K0s_d1py = 0; 
std::vector<double> *K0s_d1pz = 0; 
std::vector<double> *K0s_d1M = 0; 
std::vector<double> *K0s_d1Nhit = 0;
std::vector<double> *K0s_d1Pix = 0;
std::vector<double> *K0s_d1pterr = 0; 

//------------------------------------------
//daughter 2 (pi-)
std::vector<double> *K0s_dxy2 = 0; 
std::vector<double> *K0s_dz2 = 0; 
std::vector<double> *K0s_chi22 = 0; 
std::vector<double> *K0s_d2px = 0; 
std::vector<double> *K0s_d2py = 0; 
std::vector<double> *K0s_d2pz = 0; 
std::vector<double> *K0s_d2M = 0; 
std::vector<double> *K0s_d2Nhit = 0;
std::vector<double> *K0s_d2Pix = 0;
std::vector<double> *K0s_d2pterr = 0; 

//------------------------------------------
//Mother (K0s)
std::vector<double> *K0s_3Dagl = 0; 
std::vector<double> *K0s_3Ddl = 0; 
std::vector<double> *K0s_3Ddca = 0; 
std::vector<double> *K0s_pt = 0; 
std::vector<double> *K0s_eta = 0; 
std::vector<double> *K0s_phi = 0; 
std::vector<double> *K0s_mass = 0; 
std::vector<double> *K0s_ctau = 0; 
std::vector<double> *K0s_dca = 0; 
std::vector<double> *K0s_vtx = 0; 
//------------------------------------------
//---------------------------  K0s  -------------------------

//---------------------------  Lambda  -------------------------
//daughter 1 (p)
std::vector<double> *Lam_dxy1 = 0; 
std::vector<double> *Lam_dz1 = 0; 
std::vector<double> *Lam_chi21 = 0; 
std::vector<double> *Lam_d1px = 0; 
std::vector<double> *Lam_d1py = 0; 
std::vector<double> *Lam_d1pz = 0; 
std::vector<double> *Lam_d1M = 0;
std::vector<double> *Lam_d1Nhit = 0;
std::vector<double> *Lam_d1Pix = 0;
std::vector<double> *Lam_d1pterr = 0;  
//------------------------------------------
//daughter 2 (pi-)
std::vector<double> *Lam_dxy2 = 0; 
std::vector<double> *Lam_dz2 = 0; 
std::vector<double> *Lam_chi22 = 0; 
std::vector<double> *Lam_d2px = 0; 
std::vector<double> *Lam_d2py = 0; 
std::vector<double> *Lam_d2pz = 0; 
std::vector<double> *Lam_d2M = 0;
std::vector<double> *Lam_d2Nhit = 0;
std::vector<double> *Lam_d2Pix = 0;
std::vector<double> *Lam_d2pterr = 0; 
 
//------------------------------------------
//Mother (Lambda)
std::vector<double> *Lam_3Dagl = 0; 
std::vector<double> *Lam_3Ddl = 0; 
std::vector<double> *Lam_3Ddca = 0; 
std::vector<double> *Lam_pt = 0; 
std::vector<double> *Lam_eta = 0; 
std::vector<double> *Lam_phi = 0; 
std::vector<double> *Lam_mass = 0; 
std::vector<double> *Lam_ctau = 0; 
std::vector<double> *Lam_dca = 0; 
std::vector<double> *Lam_vtx = 0; 
std::vector<double> *Lam_id = 0; 
//------------------------------------------
//---------------------------  Lambda  -------------------------


//gen
//Int_t gk0size;

std::vector<double> *gK0s_pt = 0; 
std::vector<double> *gK0s_eta = 0; 
std::vector<double> *gK0s_phi = 0; 
std::vector<double> *gK0s_mass = 0; 
std::vector<double> *gK0s_mom1 = 0; 
std::vector<double> *gK0s_mom2 = 0;
std::vector<double> *gK0s_stat = 0;
std::vector<double> *gK0s_statmom1 = 0; 
std::vector<double> *gK0s_statmom2 = 0;


//Int_t gLamsize;
std::vector<double> *gLam_pt = 0; 
std::vector<double> *gLam_eta = 0; 
std::vector<double> *gLam_phi = 0; 
std::vector<double> *gLam_mass = 0; 
std::vector<double> *gLam_id = 0; 
std::vector<double> *gLam_mom1 = 0; 
std::vector<double> *gLam_mom2 = 0;
std::vector<double> *gLam_stat = 0;
std::vector<double> *gLam_statmom1 = 0; 
std::vector<double> *gLam_statmom2 = 0;




//for feed down case


//---------------------------  Xi  -------------------------
//daughter 1 (Lambda)
std::vector<double> *Xi_d1pt = 0; 
std::vector<double> *Xi_d1eta = 0; 
std::vector<double> *Xi_d1phi = 0; 
std::vector<double> *Xi_d1mass = 0; 

//Lambda daughters
//proton
std::vector<double> *Xi_chi21_1 = 0; 
std::vector<double> *Xi_d1pt_1 = 0; 
std::vector<double> *Xi_d1eta_1 = 0;
std::vector<double> *Xi_d1phi_1 = 0;
std::vector<double> *Xi_d1mass_1 = 0;
std::vector<double> *Xi_d1Nhit_1 = 0;  
std::vector<double> *Xi_d1Pix_1 = 0;  

//pion
std::vector<double> *Xi_chi21_2 = 0; 
std::vector<double> *Xi_d1pt_2 = 0; 
std::vector<double> *Xi_d1eta_2 = 0;
std::vector<double> *Xi_d1phi_2 = 0;
std::vector<double> *Xi_d1mass_2 = 0;
std::vector<double> *Xi_d1Nhit_2 = 0;  
std::vector<double> *Xi_d1Pix_2 = 0;  

//daughter 2 (pion)
std::vector<double> *Xi_chi22 = 0; 
std::vector<double> *Xi_d2pt = 0; 
std::vector<double> *Xi_d2eta = 0;
std::vector<double> *Xi_d2phi = 0;
std::vector<double> *Xi_d2mass = 0;
std::vector<double> *Xi_d2Nhit = 0;  
std::vector<double> *Xi_d2Pix = 0;  
 
//------------------------------------------
//Mother (Xi)
std::vector<double> *Xi_cas3DIpSigValue = 0; 
std::vector<double> *Xi_casPi3DIpSigValue = 0; 
std::vector<double> *Xi_VTrkPi3DIpSigValue = 0; 
std::vector<double> *Xi_VTrkP3DIpSigValue = 0; 
std::vector<double> *Xi_casFlightSigValue = 0; 
std::vector<double> *Xi_distanceSigValue = 0; 
std::vector<double> *Xi_pt = 0; 
std::vector<double> *Xi_eta = 0; 
std::vector<double> *Xi_phi = 0; 
std::vector<double> *Xi_mass = 0; 
std::vector<double> *Xi_id = 0; 
//---------------------------  Xi  -------------------------


//---------------------------  Omega  -------------------------
//daughter 1 (Lambda)
std::vector<double> *Om_d1pt = 0; 
std::vector<double> *Om_d1eta = 0; 
std::vector<double> *Om_d1phi = 0; 
std::vector<double> *Om_d1mass = 0; 

//Lambda daughters
//proton
std::vector<double> *Om_chi21_1 = 0; 
std::vector<double> *Om_d1pt_1 = 0; 
std::vector<double> *Om_d1eta_1 = 0;
std::vector<double> *Om_d1phi_1 = 0;
std::vector<double> *Om_d1mass_1 = 0;
std::vector<double> *Om_d1Nhit_1 = 0;  
std::vector<double> *Om_d1Pix_1 = 0;  

//pion
std::vector<double> *Om_chi21_2 = 0; 
std::vector<double> *Om_d1pt_2 = 0; 
std::vector<double> *Om_d1eta_2 = 0;
std::vector<double> *Om_d1phi_2 = 0;
std::vector<double> *Om_d1mass_2 = 0;
std::vector<double> *Om_d1Nhit_2 = 0;  
std::vector<double> *Om_d1Pix_2 = 0;  

//daughter 2 (pion)
std::vector<double> *Om_chi22 = 0; 
std::vector<double> *Om_d2pt = 0; 
std::vector<double> *Om_d2eta = 0;
std::vector<double> *Om_d2phi = 0;
std::vector<double> *Om_d2mass = 0;
std::vector<double> *Om_d2Nhit = 0;  
std::vector<double> *Om_d2Pix = 0;  
 
//------------------------------------------
//Mother (Om)
std::vector<double> *Om_cas3DIpSigValue = 0; 
std::vector<double> *Om_casPi3DIpSigValue = 0; 
std::vector<double> *Om_VTrkPi3DIpSigValue = 0; 
std::vector<double> *Om_VTrkP3DIpSigValue = 0; 
std::vector<double> *Om_casFlightSigValue = 0; 
std::vector<double> *Om_distanceSigValue = 0; 
std::vector<double> *Om_pt = 0; 
std::vector<double> *Om_eta = 0; 
std::vector<double> *Om_phi = 0; 
std::vector<double> *Om_mass = 0; 
std::vector<double> *Om_id = 0; 
//---------------------------  Omega  -------------------------


std::vector<double> *Jet_pt = 0; 
std::vector<double> *Jet_eta = 0; 
std::vector<double> *Jet_phi = 0; 
std::vector<double> *Jet_mass = 0; 


void read_tree(TTree *tree, bool isMC){

tree->AddFriend("K0SAnalysis/my_treeK0s");
if(isMC)tree->AddFriend("K0SAnalysis/my_treegK0s");
tree->AddFriend("K0SAnalysis/my_treeLam");
if(isMC)tree->AddFriend("K0SAnalysis/my_treegLam");
tree->AddFriend("K0SAnalysis/my_treeXi");
//if(isMC)tree->AddFriend("K0SAnalysis/my_treegXi");
tree->AddFriend("K0SAnalysis/my_treeOm");
//if(isMC)tree->AddFriend("K0SAnalysis/my_treegOm");
//jets
if(isMC)tree->AddFriend("K0SAnalysis/my_treeJet");
//if(isMC)tree->AddFriend("K0SAnalysis/my_treegJet");

//commom variables

tree->SetBranchStatus("*",0);     // disable all branches - this is important for big files

//enable branches

tree->SetBranchStatus("vtx_z",1); //z vertex of the detector
tree->SetBranchStatus("vtx_rho",1);//rho vertex of the detector
tree->SetBranchStatus("N_tkoff",1);//number of tracks generated per event

tree->SetBranchAddress("vtx_z",&vtx_z);
tree->SetBranchAddress("vtx_rho",&vtx_rho);
tree->SetBranchAddress("N_tkoff",&N_tkoff);

tree->SetBranchStatus("Eta_weight",1);
tree->SetBranchAddress("Eta_weight",&Eta_weight);

tree->SetBranchStatus("validPV",1); //z vertex of the detector
tree->SetBranchStatus("tracksizePV",1); //z vertex of the detector
tree->SetBranchStatus("minnTowersTh3HF",1); //z vertex of the detector
tree->SetBranchStatus("trkColl_noscr",1); //z vertex of the detector
tree->SetBranchStatus("fraction_noscr",1); //z vertex of the detector
tree->SetBranchAddress("validPV",&validPV); //z vertex of the detector
tree->SetBranchAddress("tracksizePV",&tracksizePV); //z vertex of the detector
tree->SetBranchAddress("minnTowersTh3HF",&minnTowersTh3HF); //z vertex of the detector
tree->SetBranchAddress("trkColl_noscr",&trkColl_noscr); //z vertex of the detector
tree->SetBranchAddress("fraction_noscr",&fraction_noscr); //z vertex of the detector

tree->SetBranchStatus("phffilter_1tw",1); //z vertex of the detector
tree->SetBranchStatus("pprimaryvertexfilter",1); //z vertex of the detector
tree->SetBranchStatus("pnoscrampingfilter",1); //z vertex of the detector
tree->SetBranchAddress("phffilter_1tw",&phffilter_1tw); //z vertex of the detector
tree->SetBranchAddress("pprimaryvertexfilter",&pprimaryvertexfilter); //z vertex of the detector
tree->SetBranchAddress("pnoscrampingfilter",&pnoscrampingfilter); //z vertex of the detector


tree->SetBranchStatus("pileUpFilter_pPb8TeV_vtx1",1);
tree->SetBranchStatus("pileUpFilter_pPb8TeV_Gplus",1);
tree->SetBranchStatus("olvFilter_pPb8TeV_dz1p0",1);

tree->SetBranchAddress("pileUpFilter_pPb8TeV_vtx1",&pileUpFilter_pPb8TeV_vtx1);
tree->SetBranchAddress("pileUpFilter_pPb8TeV_Gplus",&pileUpFilter_pPb8TeV_Gplus);
tree->SetBranchAddress("olvFilter_pPb8TeV_dz1p0",&olvFilter_pPb8TeV_dz1p0);

//K0s

//daughter1
tree->SetBranchStatus("K0s_dxy1",1);//xy impact parameter of the daughter 1 (pi+)
tree->SetBranchStatus("K0s_dz1",1);//z impact parameter of the daughter 1 (pi+)
tree->SetBranchStatus("K0s_chi21",1);//chi2/ndof parameter of the daughter 1 (pi+)
tree->SetBranchStatus("K0s_d1px",1);//daughter 1 (pi+) Px
tree->SetBranchStatus("K0s_d1py",1);//daughter 1 (pi+) Py
tree->SetBranchStatus("K0s_d1pz",1);//daughter 1 (pi+) Pz
tree->SetBranchStatus("K0s_d1M",1);//daughter 1 (pi+) Mass
tree->SetBranchStatus("K0s_d1Nhit",1);//daughter 1 (pi+) Nhits
tree->SetBranchStatus("K0s_d1Pix",1);//daughter 1 (pi+) Npixelhits
tree->SetBranchStatus("K0s_d1pterr",1);//daughter 1 (pi+) pterr
//daughter2
tree->SetBranchStatus("K0s_dxy2",1);//xy impact parameter of the daughter 2 (pi-)
tree->SetBranchStatus("K0s_dz2",1);//z impact parameter of the daughter 2 (pi-)
tree->SetBranchStatus("K0s_chi22",1);//chi2/ndof parameter of the daughter 2 (pi-)
tree->SetBranchStatus("K0s_d2px",1);//daughter 2 (pi-) Px
tree->SetBranchStatus("K0s_d2py",1);//daughter 2 (pi-) Py
tree->SetBranchStatus("K0s_d2pz",1);//daughter 2 (pi-) Pz
tree->SetBranchStatus("K0s_d2M",1);//daughter 2 (pi-) Mass
tree->SetBranchStatus("K0s_d2Nhit",1);//daughter 2 (pi+) Nhits
tree->SetBranchStatus("K0s_d2Pix",1);//daughter 2 (pi+) Npixelhits
tree->SetBranchStatus("K0s_d2pterr",1);//daughter 2 (pi+) pterr
// V0s
tree->SetBranchStatus("K0s_3Dagl",1);//cosx theta point (theta = angle between the momenta of K0s and primary vertex)
tree->SetBranchStatus("K0s_3Ddl",1);//decay length over the error of primary vertex and K0s
tree->SetBranchStatus("K0s_3Ddca",1);//distance of closest approuch between two daughters
tree->SetBranchStatus("K0s_pt",1);//K0s pT
tree->SetBranchStatus("K0s_eta",1);//K0s eta
tree->SetBranchStatus("K0s_phi",1);//K0s phi
tree->SetBranchStatus("K0s_mass",1);//K0s mass
tree->SetBranchStatus("K0s_ctau",1);//K0s ctau
tree->SetBranchStatus("K0s_dca",1);//K0s dca
tree->SetBranchStatus("K0s_vtx",1);//K0s vtx

//set branchs to variables

//daughter1
tree->SetBranchAddress("K0s_dxy1",&K0s_dxy1);
tree->SetBranchAddress("K0s_dz1",&K0s_dz1);
tree->SetBranchAddress("K0s_chi21",&K0s_chi21);
tree->SetBranchAddress("K0s_d1px",&K0s_d1px);
tree->SetBranchAddress("K0s_d1py",&K0s_d1py);
tree->SetBranchAddress("K0s_d1pz",&K0s_d1pz);
tree->SetBranchAddress("K0s_d1M",&K0s_d1M);
tree->SetBranchAddress("K0s_d1Nhit",&K0s_d1Nhit);
tree->SetBranchAddress("K0s_d1Pix",&K0s_d1Pix);
tree->SetBranchAddress("K0s_d1pterr",&K0s_d1pterr);


//daughter2
tree->SetBranchAddress("K0s_dxy2",&K0s_dxy2);
tree->SetBranchAddress("K0s_dz2",&K0s_dz2);
tree->SetBranchAddress("K0s_chi22",&K0s_chi22);
tree->SetBranchAddress("K0s_d2px",&K0s_d2px);
tree->SetBranchAddress("K0s_d2py",&K0s_d2py);
tree->SetBranchAddress("K0s_d2pz",&K0s_d2pz);
tree->SetBranchAddress("K0s_d2M",&K0s_d2M);
tree->SetBranchAddress("K0s_d2Nhit",&K0s_d2Nhit);
tree->SetBranchAddress("K0s_d2Pix",&K0s_d2Pix);
tree->SetBranchAddress("K0s_d2pterr",&K0s_d2pterr);

// K0s
tree->SetBranchAddress("K0s_3Dagl",&K0s_3Dagl);
tree->SetBranchAddress("K0s_3Ddl",&K0s_3Ddl);
tree->SetBranchAddress("K0s_3Ddca",&K0s_3Ddca);
tree->SetBranchAddress("K0s_pt",&K0s_pt);
tree->SetBranchAddress("K0s_eta",&K0s_eta);
tree->SetBranchAddress("K0s_phi",&K0s_phi);
tree->SetBranchAddress("K0s_mass",&K0s_mass);
tree->SetBranchAddress("K0s_ctau",&K0s_ctau);
tree->SetBranchAddress("K0s_dca",&K0s_dca);
tree->SetBranchAddress("K0s_vtx",&K0s_vtx);

if(isMC){
//Gen===============================================

tree->SetBranchStatus("gK0s_pt",1);
tree->SetBranchStatus("gK0s_eta",1);
tree->SetBranchStatus("gK0s_phi",1);
tree->SetBranchStatus("gK0s_mass",1);
tree->SetBranchStatus("gK0s_mom1",1);
tree->SetBranchStatus("gK0s_mom2",1);
tree->SetBranchStatus("gK0s_stat",1);
tree->SetBranchStatus("gK0s_statmom1",1);
tree->SetBranchStatus("gK0s_statmom2",1);

tree->SetBranchAddress("gK0s_pt",&gK0s_pt);
tree->SetBranchAddress("gK0s_eta",&gK0s_eta);
tree->SetBranchAddress("gK0s_phi",&gK0s_phi);
tree->SetBranchAddress("gK0s_mass",&gK0s_mass);
tree->SetBranchAddress("gK0s_mom1",&gK0s_mom1);
tree->SetBranchAddress("gK0s_mom2",&gK0s_mom2);
tree->SetBranchAddress("gK0s_stat",&gK0s_stat);
tree->SetBranchAddress("gK0s_statmom1",&gK0s_statmom1);
tree->SetBranchAddress("gK0s_statmom2",&gK0s_statmom2);

//==================================================
}

//daughter1
tree->SetBranchStatus("Lam_dxy1",1);//xy impact parameter of the daughter 1 (pi+)
tree->SetBranchStatus("Lam_dz1",1);//z impact parameter of the daughter 1 (pi+)
tree->SetBranchStatus("Lam_chi21",1);//chi2/ndof parameter of the daughter 1 (pi+)
tree->SetBranchStatus("Lam_d1px",1);//daughter 1 (pi+) Px
tree->SetBranchStatus("Lam_d1py",1);//daughter 1 (pi+) Py
tree->SetBranchStatus("Lam_d1pz",1);//daughter 1 (pi+) Pz
tree->SetBranchStatus("Lam_d1M",1);//daughter 1 (pi+) Mass
tree->SetBranchStatus("Lam_d1Nhit",1);//daughter 1 (pi+) Nhits
tree->SetBranchStatus("Lam_d1Pix",1);//daughter 1 (pi+) Npixelhits
tree->SetBranchStatus("Lam_d1pterr",1);//daughter 1 (pi+) pterr
//daughter2
tree->SetBranchStatus("Lam_dxy2",1);//xy impact parameter of the daughter 2 (pi-)
tree->SetBranchStatus("Lam_dz2",1);//z impact parameter of the daughter 2 (pi-)
tree->SetBranchStatus("Lam_chi22",1);//chi2/ndof parameter of the daughter 2 (pi-)
tree->SetBranchStatus("Lam_d2px",1);//daughter 2 (pi-) Px
tree->SetBranchStatus("Lam_d2py",1);//daughter 2 (pi-) Py
tree->SetBranchStatus("Lam_d2pz",1);//daughter 2 (pi-) Pz
tree->SetBranchStatus("Lam_d2M",1);//daughter 2 (pi-) Mass
tree->SetBranchStatus("Lam_d2Nhit",1);//daughter 2 (pi+) Nhits
tree->SetBranchStatus("Lam_d2Pix",1);//daughter 2 (pi+) Npixelhits
tree->SetBranchStatus("Lam_d2pterr",1);//daughter 2 (pi+) pterr
// V0s
tree->SetBranchStatus("Lam_3Dagl",1);//cosx theta point (theta = angle between the momenta of Lam and primary vertex)
tree->SetBranchStatus("Lam_3Ddl",1);//decay length over the error of primary vertex and Lam
tree->SetBranchStatus("Lam_3Ddca",1);//distance of closest approuch between two daughters
tree->SetBranchStatus("Lam_pt",1);//Lam pT
tree->SetBranchStatus("Lam_eta",1);//Lam eta
tree->SetBranchStatus("Lam_phi",1);//Lam phi
tree->SetBranchStatus("Lam_mass",1);//Lam mass
tree->SetBranchStatus("Lam_ctau",1);//Lam ctau
tree->SetBranchStatus("Lam_dca",1);//Lam dca
tree->SetBranchStatus("Lam_vtx",1);//Lam vtx
tree->SetBranchStatus("Lam_id",1);//Lam id


//set branchs to variables

//tree->SetBranchAddress("Lamsize",&Lamsize);

//daughter1
tree->SetBranchAddress("Lam_dxy1",&Lam_dxy1);
tree->SetBranchAddress("Lam_dz1",&Lam_dz1);
tree->SetBranchAddress("Lam_chi21",&Lam_chi21);
tree->SetBranchAddress("Lam_d1px",&Lam_d1px);
tree->SetBranchAddress("Lam_d1py",&Lam_d1py);
tree->SetBranchAddress("Lam_d1pz",&Lam_d1pz);
tree->SetBranchAddress("Lam_d1M",&Lam_d1M);
tree->SetBranchAddress("Lam_d1Nhit",&Lam_d1Nhit);
tree->SetBranchAddress("Lam_d1Pix",&Lam_d1Pix);
tree->SetBranchAddress("Lam_d1pterr",&Lam_d1pterr);

//daughter2
tree->SetBranchAddress("Lam_dxy2",&Lam_dxy2);
tree->SetBranchAddress("Lam_dz2",&Lam_dz2);
tree->SetBranchAddress("Lam_chi22",&Lam_chi22);
tree->SetBranchAddress("Lam_d2px",&Lam_d2px);
tree->SetBranchAddress("Lam_d2py",&Lam_d2py);
tree->SetBranchAddress("Lam_d2pz",&Lam_d2pz);
tree->SetBranchAddress("Lam_d2M",&Lam_d2M);
tree->SetBranchAddress("Lam_d2Nhit",&Lam_d2Nhit);
tree->SetBranchAddress("Lam_d2Pix",&Lam_d2Pix);
tree->SetBranchAddress("Lam_d2pterr",&Lam_d2pterr);

// Lam
tree->SetBranchAddress("Lam_3Dagl",&Lam_3Dagl);
tree->SetBranchAddress("Lam_3Ddl",&Lam_3Ddl);
tree->SetBranchAddress("Lam_3Ddca",&Lam_3Ddca);
tree->SetBranchAddress("Lam_pt",&Lam_pt);
tree->SetBranchAddress("Lam_eta",&Lam_eta);
tree->SetBranchAddress("Lam_phi",&Lam_phi);
tree->SetBranchAddress("Lam_mass",&Lam_mass);
tree->SetBranchAddress("Lam_ctau",&Lam_ctau);
tree->SetBranchAddress("Lam_dca",&Lam_dca);
tree->SetBranchAddress("Lam_vtx",&Lam_vtx);
tree->SetBranchAddress("Lam_id",&Lam_id);//Lam id

if(isMC){
//Gen===============================================
tree->SetBranchStatus("gLam_pt",1);
tree->SetBranchStatus("gLam_eta",1);
tree->SetBranchStatus("gLam_phi",1);
tree->SetBranchStatus("gLam_mass",1);
tree->SetBranchStatus("gLam_id",1);
tree->SetBranchStatus("gLam_mom1",1);
tree->SetBranchStatus("gLam_mom2",1);
tree->SetBranchStatus("gLam_stat",1);
tree->SetBranchStatus("gLam_statmom1",1);
tree->SetBranchStatus("gLam_statmom2",1);

tree->SetBranchAddress("gLam_pt",&gLam_pt);
tree->SetBranchAddress("gLam_eta",&gLam_eta);
tree->SetBranchAddress("gLam_phi",&gLam_phi);
tree->SetBranchAddress("gLam_mass",&gLam_mass);
tree->SetBranchAddress("gLam_id",&gLam_id);
tree->SetBranchAddress("gLam_mom1",&gLam_mom1);
tree->SetBranchAddress("gLam_mom2",&gLam_mom2);
tree->SetBranchAddress("gLam_stat",&gLam_stat);
tree->SetBranchAddress("gLam_statmom1",&gLam_statmom1);
tree->SetBranchAddress("gLam_statmom2",&gLam_statmom2);
//==================================================

}


//daughter1 - Lambda
tree->SetBranchStatus("Xi_d1pt",1);//daughter pt
tree->SetBranchStatus("Xi_d1eta",1);//daughter eta
tree->SetBranchStatus("Xi_d1phi",1);//daughter phi
tree->SetBranchStatus("Xi_d1mass",1);//daughter Mass

//lambda daughters
//proton
tree->SetBranchStatus("Xi_chi21_1",1);//daughter chi2/ndf
tree->SetBranchStatus("Xi_d1pt_1",1);//daughter pt
tree->SetBranchStatus("Xi_d1eta_1",1);//daughter eta
tree->SetBranchStatus("Xi_d1phi_1",1);//daughter phi
tree->SetBranchStatus("Xi_d1mass_1",1);//daughter  Mass
tree->SetBranchStatus("Xi_d1Nhit_1",1);//daughter Nhits
tree->SetBranchStatus("Xi_d1Pix_1",1);//daughter Npixelhits

//pion
tree->SetBranchStatus("Xi_chi21_2",1);//daughter chi2/ndf
tree->SetBranchStatus("Xi_d1pt_2",1);//daughter pt
tree->SetBranchStatus("Xi_d1eta_2",1);//daughter eta
tree->SetBranchStatus("Xi_d1phi_2",1);//daughter phi
tree->SetBranchStatus("Xi_d1mass_2",1);//daughter  Mass
tree->SetBranchStatus("Xi_d1Nhit_2",1);//daughter Nhits
tree->SetBranchStatus("Xi_d1Pix_2",1);//daughter Npixelhits


//daughter2 - pion
tree->SetBranchStatus("Xi_chi22",1);//chi2/ndof parameter of the daughter 2
tree->SetBranchStatus("Xi_d2pt",1);//daughter pt
tree->SetBranchStatus("Xi_d2eta",1);//daughter eta
tree->SetBranchStatus("Xi_d2phi",1);//daughter phi
tree->SetBranchStatus("Xi_d2mass",1);//daughter Mass
tree->SetBranchStatus("Xi_d2Nhit",1);//daughter Nhits
tree->SetBranchStatus("Xi_d2Pix",1);//daughter Npixelhits


// Xi
tree->SetBranchStatus("Xi_cas3DIpSigValue",1); //3D impact parameter significance of the Ξ− (Ω−) candidate with respect to the primary vertex
tree->SetBranchStatus("Xi_casPi3DIpSigValue",1); //3D impact parameter significance of π (K) track from Ξ− (Ω−) decay with respect to the primary vertex
tree->SetBranchStatus("Xi_VTrkPi3DIpSigValue",1); //3D impact parameter significance of π track from Λ decay with respect to the primary vertex
tree->SetBranchStatus("Xi_VTrkP3DIpSigValue",1);//3D impact parameter significance of proton track from Λ decay with respect to the primary vertex 
tree->SetBranchStatus("Xi_casFlightSigValue",1);//3D separation significance between Ξ− (Ω−) vertex and primary vertex
tree->SetBranchStatus("Xi_distanceSigValue",1); //3D separation significance between Λ vertex and primary vertex
tree->SetBranchStatus("Xi_pt",1);
tree->SetBranchStatus("Xi_eta",1);
tree->SetBranchStatus("Xi_phi",1);
tree->SetBranchStatus("Xi_mass",1);
tree->SetBranchStatus("Xi_id",1);

//tree->SetBranchAddress("Xisize",&Xisize);

//daughter1 - Lambda
tree->SetBranchAddress("Xi_d1pt",&Xi_d1pt);//daughter pt
tree->SetBranchAddress("Xi_d1eta",&Xi_d1eta);//daughter eta
tree->SetBranchAddress("Xi_d1phi",&Xi_d1phi);//daughter phi
tree->SetBranchAddress("Xi_d1mass",&Xi_d1mass);//daughter Mass

//lambda daughters
//proton
tree->SetBranchAddress("Xi_chi21_1",&Xi_chi21_1);//daughter chi2/ndf
tree->SetBranchAddress("Xi_d1pt_1",&Xi_d1pt_1);//daughter pt
tree->SetBranchAddress("Xi_d1eta_1",&Xi_d1eta_1);//daughter eta
tree->SetBranchAddress("Xi_d1phi_1",&Xi_d1phi_1);//daughter phi
tree->SetBranchAddress("Xi_d1mass_1",&Xi_d1mass_1);//daughter  Mass
tree->SetBranchAddress("Xi_d1Nhit_1",&Xi_d1Nhit_1);//daughter Nhits
tree->SetBranchAddress("Xi_d1Pix_1",&Xi_d1Pix_1);//daughter Npixelhits

//pion
tree->SetBranchAddress("Xi_chi21_2",&Xi_chi21_2);//daughter chi2/ndf
tree->SetBranchAddress("Xi_d1pt_2",&Xi_d1pt_2);//daughter pt
tree->SetBranchAddress("Xi_d1eta_2",&Xi_d1eta_2);//daughter eta
tree->SetBranchAddress("Xi_d1phi_2",&Xi_d1phi_2);//daughter phi
tree->SetBranchAddress("Xi_d1mass_2",&Xi_d1mass_2);//daughter  Mass
tree->SetBranchAddress("Xi_d1Nhit_2",&Xi_d1Nhit_2);//daughter Nhits
tree->SetBranchAddress("Xi_d1Pix_2",&Xi_d1Pix_2);//daughter Npixelhits


//daughter2 - pion
tree->SetBranchAddress("Xi_chi22",&Xi_chi22);//chi2/ndof parameter of the daughter 2
tree->SetBranchAddress("Xi_d2pt",&Xi_d2pt);//daughter pt
tree->SetBranchAddress("Xi_d2eta",&Xi_d2eta);//daughter eta
tree->SetBranchAddress("Xi_d2phi",&Xi_d2phi);//daughter phi
tree->SetBranchAddress("Xi_d2mass",&Xi_d2mass);//daughter Mass
tree->SetBranchAddress("Xi_d2Nhit",&Xi_d2Nhit);//daughter Nhits
tree->SetBranchAddress("Xi_d2Pix",&Xi_d2Pix);//daughter Npixelhits


// Xi
tree->SetBranchAddress("Xi_cas3DIpSigValue",&Xi_cas3DIpSigValue);
tree->SetBranchAddress("Xi_casPi3DIpSigValue",&Xi_casPi3DIpSigValue);
tree->SetBranchAddress("Xi_VTrkPi3DIpSigValue",&Xi_VTrkPi3DIpSigValue);
tree->SetBranchAddress("Xi_VTrkP3DIpSigValue",&Xi_VTrkP3DIpSigValue);
tree->SetBranchAddress("Xi_casFlightSigValue",&Xi_casFlightSigValue);
tree->SetBranchAddress("Xi_distanceSigValue",&Xi_distanceSigValue);
tree->SetBranchAddress("Xi_pt",&Xi_pt);//Lambda pT
tree->SetBranchAddress("Xi_eta",&Xi_eta);//Lambda eta
tree->SetBranchAddress("Xi_phi",&Xi_phi);
tree->SetBranchAddress("Xi_mass",&Xi_mass);
tree->SetBranchAddress("Xi_id",&Xi_id);

//daughter1 - Lambda
tree->SetBranchStatus("Om_d1pt",1);//daughter pt
tree->SetBranchStatus("Om_d1eta",1);//daughter eta
tree->SetBranchStatus("Om_d1phi",1);//daughter phi
tree->SetBranchStatus("Om_d1mass",1);//daughter Mass

//lambda daughters
//proton
tree->SetBranchStatus("Om_chi21_1",1);//daughter chi2/ndf
tree->SetBranchStatus("Om_d1pt_1",1);//daughter pt
tree->SetBranchStatus("Om_d1eta_1",1);//daughter eta
tree->SetBranchStatus("Om_d1phi_1",1);//daughter phi
tree->SetBranchStatus("Om_d1mass_1",1);//daughter  Mass
tree->SetBranchStatus("Om_d1Nhit_1",1);//daughter Nhits
tree->SetBranchStatus("Om_d1Pix_1",1);//daughter Npixelhits

//pion
tree->SetBranchStatus("Om_chi21_2",1);//daughter chi2/ndf
tree->SetBranchStatus("Om_d1pt_2",1);//daughter pt
tree->SetBranchStatus("Om_d1eta_2",1);//daughter eta
tree->SetBranchStatus("Om_d1phi_2",1);//daughter phi
tree->SetBranchStatus("Om_d1mass_2",1);//daughter  Mass
tree->SetBranchStatus("Om_d1Nhit_2",1);//daughter Nhits
tree->SetBranchStatus("Om_d1Pix_2",1);//daughter Npixelhits


//daughter2 - pion
tree->SetBranchStatus("Om_chi22",1);//chi2/ndof parameter of the daughter 2
tree->SetBranchStatus("Om_d2pt",1);//daughter pt
tree->SetBranchStatus("Om_d2eta",1);//daughter eta
tree->SetBranchStatus("Om_d2phi",1);//daughter phi
tree->SetBranchStatus("Om_d2mass",1);//daughter Mass
tree->SetBranchStatus("Om_d2Nhit",1);//daughter Nhits
tree->SetBranchStatus("Om_d2Pix",1);//daughter Npixelhits


// Om
tree->SetBranchStatus("Om_cas3DIpSigValue",1);  //3D impact parameter significance of the Ξ− (Ω−) candidate with respect to the primary vertex
tree->SetBranchStatus("Om_casPi3DIpSigValue",1);//3D impact parameter significance of π (K) track from Ξ− (Ω−) decay with respect to the primary vertex
tree->SetBranchStatus("Om_VTrkPi3DIpSigValue",1);//3D impact parameter significance of π track from Λ decay with respect to the primary vertex
tree->SetBranchStatus("Om_VTrkP3DIpSigValue",1);//3D impact parameter significance of proton track from Λ decay with respect to the primary vertex 
tree->SetBranchStatus("Om_casFlightSigValue",1);//3D separation significance between Ξ− (Ω−) vertex and primary vertex
tree->SetBranchStatus("Om_distanceSigValue",1); //3D separation significance between Λ vertex and primary vertex
tree->SetBranchStatus("Om_pt",1);
tree->SetBranchStatus("Om_eta",1);
tree->SetBranchStatus("Om_phi",1);
tree->SetBranchStatus("Om_mass",1);
tree->SetBranchStatus("Om_id",1);



//daughter1 - Lambda
tree->SetBranchAddress("Om_d1pt",&Om_d1pt);//daughter pt
tree->SetBranchAddress("Om_d1eta",&Om_d1eta);//daughter eta
tree->SetBranchAddress("Om_d1phi",&Om_d1phi);//daughter phi
tree->SetBranchAddress("Om_d1mass",&Om_d1mass);//daughter Mass

//lambda daughters
//proton
tree->SetBranchAddress("Om_chi21_1",&Om_chi21_1);//daughter chi2/ndf
tree->SetBranchAddress("Om_d1pt_1",&Om_d1pt_1);//daughter pt
tree->SetBranchAddress("Om_d1eta_1",&Om_d1eta_1);//daughter eta
tree->SetBranchAddress("Om_d1phi_1",&Om_d1phi_1);//daughter phi
tree->SetBranchAddress("Om_d1mass_1",&Om_d1mass_1);//daughter  Mass
tree->SetBranchAddress("Om_d1Nhit_1",&Om_d1Nhit_1);//daughter Nhits
tree->SetBranchAddress("Om_d1Pix_1",&Om_d1Pix_1);//daughter Npixelhits

//pion
tree->SetBranchAddress("Om_chi21_2",&Om_chi21_2);//daughter chi2/ndf
tree->SetBranchAddress("Om_d1pt_2",&Om_d1pt_2);//daughter pt
tree->SetBranchAddress("Om_d1eta_2",&Om_d1eta_2);//daughter eta
tree->SetBranchAddress("Om_d1phi_2",&Om_d1phi_2);//daughter phi
tree->SetBranchAddress("Om_d1mass_2",&Om_d1mass_2);//daughter  Mass
tree->SetBranchAddress("Om_d1Nhit_2",&Om_d1Nhit_2);//daughter Nhits
tree->SetBranchAddress("Om_d1Pix_2",&Om_d1Pix_2);//daughter Npixelhits


//daughter2 - Kaon
tree->SetBranchAddress("Om_chi22",&Om_chi22);//chi2/ndof parameter of the daughter 2
tree->SetBranchAddress("Om_d2pt",&Om_d2pt);//daughter pt
tree->SetBranchAddress("Om_d2eta",&Om_d2eta);//daughter eta
tree->SetBranchAddress("Om_d2phi",&Om_d2phi);//daughter phi
tree->SetBranchAddress("Om_d2mass",&Om_d2mass);//daughter Mass
tree->SetBranchAddress("Om_d2Nhit",&Om_d2Nhit);//daughter Nhits
tree->SetBranchAddress("Om_d2Pix",&Om_d2Pix);//daughter Npixelhits


// Om
tree->SetBranchAddress("Om_cas3DIpSigValue",&Om_cas3DIpSigValue);
tree->SetBranchAddress("Om_casPi3DIpSigValue",&Om_casPi3DIpSigValue);
tree->SetBranchAddress("Om_VTrkPi3DIpSigValue",&Om_VTrkPi3DIpSigValue);
tree->SetBranchAddress("Om_VTrkP3DIpSigValue",&Om_VTrkP3DIpSigValue);
tree->SetBranchAddress("Om_casFlightSigValue",&Om_casFlightSigValue);
tree->SetBranchAddress("Om_distanceSigValue",&Om_distanceSigValue);
tree->SetBranchAddress("Om_pt",&Om_pt);//Lambda pT
tree->SetBranchAddress("Om_eta",&Om_eta);//Lambda eta
tree->SetBranchAddress("Om_phi",&Om_phi);
tree->SetBranchAddress("Om_mass",&Om_mass);
tree->SetBranchAddress("Om_id",&Om_id);

if(isMC){

tree->SetBranchStatus("Jet_pt",1);
tree->SetBranchStatus("Jet_eta",1);
tree->SetBranchStatus("Jet_phi",1);
tree->SetBranchStatus("Jet_mass",1);
tree->SetBranchAddress("Jet_pt",&Jet_pt);
tree->SetBranchAddress("Jet_eta",&Jet_eta);
tree->SetBranchAddress("Jet_phi",&Jet_phi);
tree->SetBranchAddress("Jet_mass",&Jet_mass);

}


} //DONE

//----------------------------------------------------------------------------------------------------------------------

//Important

//----------------------------------------------------------------------------------------------------------------------

double getTrkCorrWeight(TH2 *effhisto, double pT, double eta){
double efficiency = effhisto->GetBinContent(effhisto->GetXaxis()->FindBin(pT),effhisto->GetYaxis()->FindBin(eta) );
if(efficiency > 0.01){return 1./efficiency;}else{return 0;}
} //DONE

const TLorentzVector InvertPVector( TLorentzVector &vec){
   TLorentzVector ovec = vec;
   ovec.SetPx(-vec.Px());
   ovec.SetPy(-vec.Py());
   ovec.SetPz(-vec.Pz());
   return ovec;
} //DONE

const TLorentzVector InvertXYVector( TLorentzVector &vec){
  TLorentzVector ovec = vec;
  ovec.SetX(-vec.X());
  ovec.SetY(-vec.Y());
  return ovec;
} //DONE

bool splitcomb(TLorentzVector &vec1,TLorentzVector &vec2){
   bool issplit=false;
   Double_t cosa = TMath::Abs(vec1.Px()*vec2.Px() + vec1.Py()*vec2.Py() + vec1.Pz()*vec2.Pz())/(vec1.P()*vec2.P());
   Double_t deltapt = TMath::Abs(vec1.Pt() - vec2.Pt());if ( (cosa > 0.99996) && (deltapt < 0.04)) { issplit = true;} //same used by Sunil in pPb
   return issplit;
} //DONE

Double_t GetQ(const TLorentzVector &p1, const TLorentzVector &p2){
  TLorentzVector pSum = p1 + p2;
  TLorentzVector pDiff = p1 - p2;
  TLorentzVector qinvVec = (pDiff - (pDiff * pSum/pSum.Mag2())*pSum);
  return -1.*qinvVec.Mag();
} //DONE

bool etaMixSort(std::pair<Double_t,std::vector<TLorentzVector>> &a, std::pair<Double_t,std::vector<TLorentzVector>> &b){

   return a.first < b.first ? true : false ;

} //DONE

bool etaMixSort_ntrkoff(std::pair<Double_t,int> &a, std::pair<Double_t,int> &b){

   return a.first < b.first ? true : false ;

}  //DONE


//----------------------------------------------------------------------------------------------------------------------

//mass dep

//----------------------------------------------------------------------------------------------------------------------

void Mass_dep(std::vector<TLorentzVector> GoodTrackFourVector, Double_t aux_N_tk_offline, THnSparse* h1, double weight){
if(GoodTrackFourVector.size()>=2){ 
//just for quantities mass kT dep.                        
for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){ 
if(itk1!=0){break;}
Double_t mmmm = GoodTrackFourVector[itk1].M();
for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
    Double_t mmmm2 = GoodTrackFourVector[itk2].M();
    Double_t q = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
    TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
    Double_t ktX=(psum2X.Pt())/2.;
    Double_t x3DX2[4]={mmmm2,ktX,(Double_t)aux_N_tk_offline,q};
    h1->Fill(x3DX2,weight);
    if(itk2==1){Double_t x3DX[4]={mmmm,ktX,(Double_t)aux_N_tk_offline,q};
    h1->Fill(x3DX,weight);}
}}}}


void Mass_dep_sep(std::vector<TLorentzVector> GoodTrackFourVector, Double_t aux_N_tk_offline, THnSparse* h1, THnSparse* h2, double weight){
if(GoodTrackFourVector.size()>=2){ 
//just for quantities mass kT qinv and ntk dep.                        
for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){ 
if(itk1!=0){break;}
Double_t mmmm = GoodTrackFourVector[itk1].M();
for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
    Double_t mmmm2 = GoodTrackFourVector[itk2].M();
    Double_t q = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
    TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
    Double_t ktX=(psum2X.Pt())/2.;
    Double_t x3DX2[4]={mmmm2,ktX,(Double_t)aux_N_tk_offline,q};
    h2->Fill(x3DX2,weight);
    if(itk2==1){Double_t x3DX[4]={mmmm,ktX,(Double_t)aux_N_tk_offline,q};
    h1->Fill(x3DX,weight);}
}}}}


void Mass_dep_OS(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVector2, Double_t aux_N_tk_offline, THnSparse* h1, THnSparse* h2,  THnSparse* h3, double weight){
if(GoodTrackFourVector.size()>=1 && GoodTrackFourVector2.size()>=1){ 
//just for quantities mass kT dep.                        
for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){ 
Double_t mmmm = GoodTrackFourVector[itk1].M();
for(unsigned int itk2=0; itk2<GoodTrackFourVector2.size();itk2++){
    Double_t mmmm2 = GoodTrackFourVector2[itk2].M();
    Double_t q = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector2[itk2]);
    TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector2[itk2];
    Double_t ktX=(psum2X.Pt())/2.;
    Double_t x3DX2[4]={mmmm2,ktX,(Double_t)aux_N_tk_offline,q};
    Double_t x3DX[4]={mmmm,ktX,(Double_t)aux_N_tk_offline,q};
    if(itk1==0){h2->Fill(x3DX2,weight);h3->Fill(x3DX2,weight);}
    if(itk2==0){h1->Fill(x3DX,weight);h3->Fill(x3DX,weight);}
}}}}


//----------------------------------------------------------------------------------------------------------------------

//chi2 checks and plots - do not change

//----------------------------------------------------------------------------------------------------------------------

void get_chi2plot(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<double> GoodTrackchi21, std::vector<double> GoodTrackchi22, TH1* h_11, TH1* h_22, TH1* h_1221, double weight){
for(Int_t ik=0; ik<GoodTrackFourVector.size(); ik++){
double chi2d1 = GoodTrackchi21[ik];
double chi2d2 = GoodTrackchi22[ik];
for(Int_t iik=ik+1; iik<GoodTrackFourVector.size(); iik++){
double chi2d1a = GoodTrackchi21[iik];
double chi2d2a = GoodTrackchi22[iik]; 
h_11->Fill(fabs(chi2d1 - chi2d1a),weight);
h_22->Fill(fabs(chi2d2 - chi2d2a),weight);
h_1221->Fill(fabs(chi2d1 - chi2d2a),weight);
h_1221->Fill(fabs(chi2d1a - chi2d2),weight);
}}}

void get_chi2plot_TF(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<double> GoodTrackchi21, std::vector<double> GoodTrackchi22, std::vector<TLorentzVector> GoodTrackFourVectorX, std::vector<double> GoodTrackchi21X, std::vector<double> GoodTrackchi22X, TH1* h_11, TH1* h_22, TH1* h_1221, double weight){
for(Int_t ik=0; ik<GoodTrackFourVector.size(); ik++){
double chi2d1 = GoodTrackchi21[ik];
double chi2d2 = GoodTrackchi22[ik];
for(Int_t iik=0; iik<GoodTrackFourVectorX.size(); iik++){
double chi2d1a = GoodTrackchi21X[iik];
double chi2d2a = GoodTrackchi22X[iik]; 
h_11->Fill(fabs(chi2d1 - chi2d1a),weight);
h_22->Fill(fabs(chi2d2 - chi2d2a),weight);
h_1221->Fill(fabs(chi2d1 - chi2d2a),weight);
h_1221->Fill(fabs(chi2d1a - chi2d2),weight);
}}}

void get_chi2plot_TF_new(std::vector<double> GoodTrackchi21, std::vector<double> GoodTrackchi22, std::vector<double> GoodTrackchi21X, std::vector<double> GoodTrackchi22X, TH1* h_11, TH1* h_1221, double weight){
for(Int_t ik=0; ik<GoodTrackchi21.size(); ik++){
double chi2d1 = GoodTrackchi21[ik];
double chi2d2 = GoodTrackchi22[ik];
for(Int_t iik=0; iik<GoodTrackchi21X.size(); iik++){
double chi2d1a = GoodTrackchi21X[iik];
double chi2d2a = GoodTrackchi22X[iik]; 
h_11->Fill(fabs(chi2d1 - chi2d1a),weight);
h_11->Fill(fabs(chi2d2 - chi2d2a),weight);
h_1221->Fill(fabs(chi2d1 - chi2d2a),weight);
h_1221->Fill(fabs(chi2d1a - chi2d2),weight);
}}}



//until here do not need correction and it is done


//-----------------------------------------------------------------------------
//HBT
//-----------------------------------------------------------------------------


void hbt(std::vector<TLorentzVector> GoodTrackFourVector, Double_t aux_N_tk_offline, THnSparse* h1){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
       Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        h1->Fill(x3DX);
}}}


void hbtOS(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVector2, Double_t aux_N_tk_offline, THnSparse* h1){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    for(unsigned int itk2=0; itk2<GoodTrackFourVector2.size();itk2++){
       Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector2[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector2[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        h1->Fill(x3DX);
}}}


void hbt_reco(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVectord1, std::vector<TLorentzVector> GoodTrackFourVectord2, Double_t aux_N_tk_offline, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h1x, THnSparse* h11, THnSparse* h12, THnSparse* h12x, THnSparse* h1rot, THnSparse* h1inv, THnSparse* h2, THnSparse* h2rot, THnSparse* h2inv, THnSparse* h3, THnSparse* h3rot, THnSparse* h3inv, THnSparse* hbbL, THnSparse* hbbR, THnSparse* hsbL, THnSparse* hsbR, int p, TH2* effhist, TH1* nevss, TH1* nevbb, TH1* nevsb, double weight){
      
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    double eff1 = 1.0;
    if(p==20 || p==21 || p==22){eff1 = getTrkCorrWeight(effhist,GoodTrackFourVector[itk1].Pt(),GoodTrackFourVector[itk1].Rapidity());}
      for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
        double eff2 = 1.0;
        if(p==20 || p==21 || p==22){eff2 = getTrkCorrWeight(effhist,GoodTrackFourVector[itk2].Pt(),GoodTrackFourVector[itk2].Rapidity());}
        double toteff = eff1*eff2;
        //V0 qinv
        Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        //daughter1 qinv
        Double_t qX1 = GetQ(GoodTrackFourVectord1[itk1], GoodTrackFourVectord1[itk2]);
    	TLorentzVector psum2X1 = GoodTrackFourVectord1[itk1] + GoodTrackFourVectord1[itk2];
        Double_t ktX1=(psum2X1.Pt())/2.;
        Double_t x3DX1[3]={qX1,ktX1,aux_N_tk_offline};
        //daughter2 qinv
        Double_t qX2 = GetQ(GoodTrackFourVectord2[itk1], GoodTrackFourVectord2[itk2]);
    	TLorentzVector psum2X2 = GoodTrackFourVectord2[itk1] + GoodTrackFourVectord2[itk2];
        Double_t ktX2=(psum2X2.Pt())/2.;
        Double_t x3DX2[3]={qX2,ktX2,aux_N_tk_offline};
        //daughter12 qinv
        Double_t qX3 = GetQ(GoodTrackFourVectord1[itk1], GoodTrackFourVectord2[itk2]);
    	TLorentzVector psum2X3 = GoodTrackFourVectord1[itk1] + GoodTrackFourVectord2[itk2];
        Double_t ktX3=(psum2X3.Pt())/2.;
        Double_t x3DX3[3]={qX3,ktX3,aux_N_tk_offline};
         //daughter21 qinv       
        Double_t qX4 = GetQ(GoodTrackFourVectord2[itk1], GoodTrackFourVectord1[itk2]);
    	TLorentzVector psum2X4 = GoodTrackFourVectord2[itk1] + GoodTrackFourVectord1[itk2];
        Double_t ktX4=(psum2X4.Pt())/2.;
        Double_t x3DX4[3]={qX4,ktX4,aux_N_tk_offline};
        
        //Rotated XY
        Double_t q_rot = GetQ(GoodTrackFourVector[itk1], InvertXYVector(GoodTrackFourVector[itk2]));        
        TLorentzVector psum2_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector[itk2]);
        Double_t kt_rot = (psum2_rot.Pt())/2.;
        Double_t x3D_rot[3]={q_rot,kt_rot,(Double_t)aux_N_tk_offline};
       
        //Inverted PxPyPz
        Double_t q_inv = GetQ(GoodTrackFourVector[itk1], InvertPVector(GoodTrackFourVector[itk2]));
        TLorentzVector psum2_inv = GoodTrackFourVector[itk1] + InvertPVector(GoodTrackFourVector[itk2]);
        Double_t kt_inv = (psum2_inv.Pt())/2.;
        Double_t x3D_inv[3]={q_inv,kt_inv,(Double_t)aux_N_tk_offline};             
        
        if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma)){
        h1->Fill(x3DX,toteff*weight);
        h11->Fill(x3DX1,toteff*weight);    
        h12->Fill(x3DX2,toteff*weight);   
        h12x->Fill(x3DX3,toteff*weight);
        h12x->Fill(x3DX4,toteff*weight);
        h1rot->Fill(x3D_rot,toteff*weight);
        h1inv->Fill(x3D_inv,toteff*weight);

        if(itk1==0 && itk2==1)nevss->Fill(1);
        if(!splitcomb(GoodTrackFourVectord1[itk1],  GoodTrackFourVectord1[itk2]) && !splitcomb(GoodTrackFourVectord2[itk1],  GoodTrackFourVectord2[itk2]) && !splitcomb(GoodTrackFourVectord1[itk1],  GoodTrackFourVectord2[itk2]) && !splitcomb(GoodTrackFourVectord2[itk1],  GoodTrackFourVectord1[itk2])){h1x->Fill(x3DX,toteff*weight);}
        } //ss
        
        if((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk2].M() > mass + nsigmaside*sigma)){
        h2->Fill(x3DX,toteff*weight);
        h2rot->Fill(x3D_rot,toteff*weight);
        h2inv->Fill(x3D_inv,toteff*weight);
        if(itk1==0 && itk2==1)nevbb->Fill(1);
        } //bb
        
        if(((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk2].M() > mass + nsigmaside*sigma)) || ((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma))){
        h3->Fill(x3DX,toteff*weight);
        h3rot->Fill(x3D_rot,toteff*weight);
        h3inv->Fill(x3D_inv,toteff*weight);    
        if(itk1==0 && itk2==1)nevsb->Fill(1);
        } //sb

        if((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() < mass - nsigmaside*sigma)){hbbL->Fill(x3DX,toteff*weight);}
        
        if((GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() > mass + nsigmaside*sigma)){hbbR->Fill(x3DX,toteff*weight);}

        if(((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() < mass - nsigmaside*sigma)) || ((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma))){hsbL->Fill(x3DX,toteff*weight);}
        
        if(((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() > mass + nsigmaside*sigma)) || ((GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma))){hsbR->Fill(x3DX,toteff*weight);}
       
    }}}

//removed and bias test

void hbt_reco_fake(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVectord1, std::vector<TLorentzVector> GoodTrackFourVectord2, Double_t aux_N_tk_offline, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h11, THnSparse* h12, THnSparse* h12x, THnSparse* h2, THnSparse* h3, THnSparse* hX, THnSparse* hX1, THnSparse* hX2, THnSparse* hX12, double weight){
    
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
       //V0 qinv
        Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        //daughter1 qinv
        Double_t qX1 = GetQ(GoodTrackFourVectord1[itk1], GoodTrackFourVectord1[itk2]);
    	TLorentzVector psum2X1 = GoodTrackFourVectord1[itk1] + GoodTrackFourVectord1[itk2];
        Double_t ktX1=(psum2X1.Pt())/2.;
        Double_t x3DX1[3]={qX1,ktX1,aux_N_tk_offline};
        //daughter2 qinv
        Double_t qX2 = GetQ(GoodTrackFourVectord2[itk1], GoodTrackFourVectord2[itk2]);
    	TLorentzVector psum2X2 = GoodTrackFourVectord2[itk1] + GoodTrackFourVectord2[itk2];
        Double_t ktX2=(psum2X2.Pt())/2.;
        Double_t x3DX2[3]={qX2,ktX2,aux_N_tk_offline};
        //daughter12 qinv
        Double_t qX3 = GetQ(GoodTrackFourVectord1[itk1], GoodTrackFourVectord2[itk2]);
    	TLorentzVector psum2X3 = GoodTrackFourVectord1[itk1] + GoodTrackFourVectord2[itk2];
        Double_t ktX3=(psum2X3.Pt())/2.;
        Double_t x3DX3[3]={qX3,ktX3,aux_N_tk_offline};
         //daughter21 qinv       
        Double_t qX4 = GetQ(GoodTrackFourVectord2[itk1], GoodTrackFourVectord1[itk2]);
    	TLorentzVector psum2X4 = GoodTrackFourVectord2[itk1] + GoodTrackFourVectord1[itk2];
        Double_t ktX4=(psum2X4.Pt())/2.;
        Double_t x3DX4[3]={qX4,ktX4,aux_N_tk_offline};    
        
        if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma)){
        h1->Fill(x3DX,weight);
        h11->Fill(x3DX1,weight);    
        h12->Fill(x3DX2,weight);   
        h12x->Fill(x3DX3,weight);
        h12x->Fill(x3DX4,weight);
        hX->Fill(x3DX,weight);
        hX1->Fill(x3DX1,weight);
        hX2->Fill(x3DX2,weight);
        hX12->Fill(x3DX3,weight);
        hX12->Fill(x3DX4,weight);
        } 
        
        if((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk2].M() > mass + nsigmaside*sigma)){h2->Fill(x3DX,weight);}        
        if(((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk2].M() > mass + nsigmaside*sigma)) || ((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma))){h3->Fill(x3DX,weight);}

    }}
    
}


void hbt_reco_truefake(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVectord1, std::vector<TLorentzVector> GoodTrackFourVectord2, std::vector<TLorentzVector> GoodTrackFourVectorX, std::vector<TLorentzVector> GoodTrackFourVectord1X, std::vector<TLorentzVector> GoodTrackFourVectord2X, Double_t aux_N_tk_offline, double mass, double sigma, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h11, THnSparse* h12, THnSparse* h12x, THnSparse* h2, THnSparse* h3, THnSparse* hX, THnSparse* hX1, THnSparse* hX2, THnSparse* hX12, double weight){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    for(unsigned int itk2=0; itk2<GoodTrackFourVectorX.size();itk2++){
       //V0 qinv
        Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVectorX[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVectorX[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        //daughter1 qinv
        Double_t qX1 = GetQ(GoodTrackFourVectord1[itk1], GoodTrackFourVectord1X[itk2]);
    	TLorentzVector psum2X1 = GoodTrackFourVectord1[itk1] + GoodTrackFourVectord1X[itk2];
        Double_t ktX1=(psum2X1.Pt())/2.;
        Double_t x3DX1[3]={qX1,ktX1,aux_N_tk_offline};
        //daughter2 qinv
        Double_t qX2 = GetQ(GoodTrackFourVectord2[itk1], GoodTrackFourVectord2X[itk2]);
    	TLorentzVector psum2X2 = GoodTrackFourVectord2[itk1] + GoodTrackFourVectord2X[itk2];
        Double_t ktX2=(psum2X2.Pt())/2.;
        Double_t x3DX2[3]={qX2,ktX2,aux_N_tk_offline};
        //daughter12 qinv
        Double_t qX3 = GetQ(GoodTrackFourVectord1[itk1], GoodTrackFourVectord2X[itk2]);
    	TLorentzVector psum2X3 = GoodTrackFourVectord1[itk1] + GoodTrackFourVectord2X[itk2];
        Double_t ktX3=(psum2X3.Pt())/2.;
        Double_t x3DX3[3]={qX3,ktX3,aux_N_tk_offline};
         //daughter21 qinv       
        Double_t qX4 = GetQ(GoodTrackFourVectord2[itk1], GoodTrackFourVectord1X[itk2]);
    	TLorentzVector psum2X4 = GoodTrackFourVectord2[itk1] + GoodTrackFourVectord1X[itk2];
        Double_t ktX4=(psum2X4.Pt())/2.;
        Double_t x3DX4[3]={qX4,ktX4,aux_N_tk_offline};    
        
        if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVectorX[itk2].M() > mass2 - nsigmapeak*sigma2 && GoodTrackFourVectorX[itk2].M() < mass2 + nsigmapeak*sigma2)){
        h1->Fill(x3DX,weight);
        h11->Fill(x3DX1,weight);    
        h12->Fill(x3DX2,weight);   
        h12x->Fill(x3DX3,weight);
        h12x->Fill(x3DX4,weight);
        hX->Fill(x3DX,weight);
        hX1->Fill(x3DX1,weight);
        hX2->Fill(x3DX2,weight);
        hX12->Fill(x3DX3,weight);
        hX12->Fill(x3DX4,weight);
        } 
        
        if((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVectorX[itk2].M() < mass2 - nsigmaside*sigma2 || GoodTrackFourVectorX[itk2].M() > mass2 + nsigmaside*sigma2)){h2->Fill(x3DX,weight);}        
        if(((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVectorX[itk2].M() < mass2 - nsigmaside*sigma2 || GoodTrackFourVectorX[itk2].M() > mass2 + nsigmaside*sigma2)) || ((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVectorX[itk2].M() > mass2 - nsigmapeak*sigma2 && GoodTrackFourVectorX[itk2].M() < mass2 + nsigmapeak*sigma2))){h3->Fill(x3DX,weight);}

    }}
    
}

//match && unmatch

void hbt_MM_UU(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVectord1, std::vector<TLorentzVector> GoodTrackFourVectord2, Double_t aux_N_tk_offline, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h11, THnSparse* h12, THnSparse* h12x, THnSparse* hX, THnSparse* hX1, THnSparse* hX2, THnSparse* hX12, double weight){
    
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
       //V0 qinv
        Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        //daughter1 qinv
        Double_t qX1 = GetQ(GoodTrackFourVectord1[itk1], GoodTrackFourVectord1[itk2]);
    	TLorentzVector psum2X1 = GoodTrackFourVectord1[itk1] + GoodTrackFourVectord1[itk2];
        Double_t ktX1=(psum2X1.Pt())/2.;
        Double_t x3DX1[3]={qX1,ktX1,aux_N_tk_offline};
        //daughter2 qinv
        Double_t qX2 = GetQ(GoodTrackFourVectord2[itk1], GoodTrackFourVectord2[itk2]);
    	TLorentzVector psum2X2 = GoodTrackFourVectord2[itk1] + GoodTrackFourVectord2[itk2];
        Double_t ktX2=(psum2X2.Pt())/2.;
        Double_t x3DX2[3]={qX2,ktX2,aux_N_tk_offline};
        //daughter12 qinv
        Double_t qX3 = GetQ(GoodTrackFourVectord1[itk1], GoodTrackFourVectord2[itk2]);
    	TLorentzVector psum2X3 = GoodTrackFourVectord1[itk1] + GoodTrackFourVectord2[itk2];
        Double_t ktX3=(psum2X3.Pt())/2.;
        Double_t x3DX3[3]={qX3,ktX3,aux_N_tk_offline};
         //daughter21 qinv       
        Double_t qX4 = GetQ(GoodTrackFourVectord2[itk1], GoodTrackFourVectord1[itk2]);
    	TLorentzVector psum2X4 = GoodTrackFourVectord2[itk1] + GoodTrackFourVectord1[itk2];
        Double_t ktX4=(psum2X4.Pt())/2.;
        Double_t x3DX4[3]={qX4,ktX4,aux_N_tk_offline};    
        
        if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma)){
        h1->Fill(x3DX,weight);
        h11->Fill(x3DX1,weight);    
        h12->Fill(x3DX2,weight);   
        h12x->Fill(x3DX3,weight);
        h12x->Fill(x3DX4,weight);
        hX->Fill(x3DX,weight);
        hX1->Fill(x3DX1,weight);
        hX2->Fill(x3DX2,weight);
        hX12->Fill(x3DX3,weight);
        hX12->Fill(x3DX4,weight);
        } 
        
    }}
    
}

//cross

void hbt_reco_cross(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVectord1, std::vector<TLorentzVector> GoodTrackFourVectord2,std::vector<TLorentzVector> GoodTrackFourVectorX, std::vector<TLorentzVector> GoodTrackFourVectord1X, std::vector<TLorentzVector> GoodTrackFourVectord2X, Double_t aux_N_tk_offline, double mass, double sigma, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h1x, THnSparse* h11, THnSparse* h12, THnSparse* h12x, THnSparse* h1rot, THnSparse* h1inv, THnSparse* h2, THnSparse* h2rot, THnSparse* h2inv, THnSparse* h3, THnSparse* h3rot, THnSparse* h3inv, THnSparse* hbbL, THnSparse* hbbR, THnSparse* hsbL, THnSparse* hsbR, int p, TH2* effhist, TH2* effhist2, TH1* nevss, TH1* nevbb, TH1* nevsb, double weight){
    
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    double eff1 = 1.0;
    if(p==20 || p==21 || p==22)eff1 = getTrkCorrWeight(effhist,GoodTrackFourVector[itk1].Pt(),GoodTrackFourVector[itk1].Rapidity());
      for(unsigned int itk2=0; itk2<GoodTrackFourVectorX.size();itk2++){
        double eff2 = 1.0;
        if(p==20 || p==21 || p==22)eff2 = getTrkCorrWeight(effhist2,GoodTrackFourVectorX[itk2].Pt(),GoodTrackFourVectorX[itk2].Rapidity());
        double toteff = eff1*eff2;
        //V0 qinv
        Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVectorX[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVectorX[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        //daughter1 qinv
        Double_t qX1 = GetQ(GoodTrackFourVectord1[itk1], GoodTrackFourVectord1X[itk2]);
    	TLorentzVector psum2X1 = GoodTrackFourVectord1[itk1] + GoodTrackFourVectord1X[itk2];
        Double_t ktX1=(psum2X1.Pt())/2.;
        Double_t x3DX1[3]={qX1,ktX1,aux_N_tk_offline};
        //daughter2 qinv
        Double_t qX2 = GetQ(GoodTrackFourVectord2[itk1], GoodTrackFourVectord2X[itk2]);
    	TLorentzVector psum2X2 = GoodTrackFourVectord2[itk1] + GoodTrackFourVectord2X[itk2];
        Double_t ktX2=(psum2X2.Pt())/2.;
        Double_t x3DX2[3]={qX2,ktX2,aux_N_tk_offline};
        //daughter12 qinv
        Double_t qX3 = GetQ(GoodTrackFourVectord1[itk1], GoodTrackFourVectord2X[itk2]);
    	TLorentzVector psum2X3 = GoodTrackFourVectord1[itk1] + GoodTrackFourVectord2X[itk2];
        Double_t ktX3=(psum2X3.Pt())/2.;
        Double_t x3DX3[3]={qX3,ktX3,aux_N_tk_offline};
         //daughter21 qinv       
        Double_t qX4 = GetQ(GoodTrackFourVectord2[itk1], GoodTrackFourVectord1X[itk2]);
    	TLorentzVector psum2X4 = GoodTrackFourVectord2[itk1] + GoodTrackFourVectord1X[itk2];
        Double_t ktX4=(psum2X4.Pt())/2.;
        Double_t x3DX4[3]={qX4,ktX4,aux_N_tk_offline};
        
        //Rotated XY
        Double_t q_rot = GetQ(GoodTrackFourVector[itk1], InvertXYVector(GoodTrackFourVectorX[itk2]));        
        TLorentzVector psum2_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVectorX[itk2]);
        Double_t kt_rot = (psum2_rot.Pt())/2.;
        Double_t x3D_rot[3]={q_rot,kt_rot,(Double_t)aux_N_tk_offline};
       
        //Inverted PxPyPz
        Double_t q_inv = GetQ(GoodTrackFourVector[itk1], InvertPVector(GoodTrackFourVectorX[itk2]));
        TLorentzVector psum2_inv = GoodTrackFourVector[itk1] + InvertPVector(GoodTrackFourVectorX[itk2]);
        Double_t kt_inv = (psum2_inv.Pt())/2.;
        Double_t x3D_inv[3]={q_inv,kt_inv,(Double_t)aux_N_tk_offline};             
      
        if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVectorX[itk2].M() > mass2 - nsigmapeak*sigma2 && GoodTrackFourVectorX[itk2].M() < mass2 + nsigmapeak*sigma2)){
        h1->Fill(x3DX,toteff*weight);
        h11->Fill(x3DX1,toteff*weight);    
        h12->Fill(x3DX2,toteff*weight);   
        h12x->Fill(x3DX3,toteff*weight);
        h12x->Fill(x3DX4,toteff*weight);
        h1rot->Fill(x3D_rot,toteff*weight);
        h1inv->Fill(x3D_inv,toteff*weight);
        if(itk1==0 && itk2==1)nevss->Fill(1);
        if(splitcomb(GoodTrackFourVectord1[itk1],  GoodTrackFourVectord1X[itk2]))continue;
        if(splitcomb(GoodTrackFourVectord2[itk1],  GoodTrackFourVectord2X[itk2]))continue;
        if(splitcomb(GoodTrackFourVectord1[itk1],  GoodTrackFourVectord2X[itk2]))continue;
        if(splitcomb(GoodTrackFourVectord2[itk1],  GoodTrackFourVectord1X[itk2]))continue;
        h1x->Fill(x3DX,toteff*weight);
        }
        
        if((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVectorX[itk2].M() < mass2 - nsigmaside*sigma2 || GoodTrackFourVectorX[itk2].M() > mass2 + nsigmaside*sigma2)){
        h2->Fill(x3DX,toteff*weight);
        h2rot->Fill(x3D_rot,toteff*weight);
        h2inv->Fill(x3D_inv,toteff*weight);
        if(itk1==0 && itk2==1)nevbb->Fill(1);
        }
        
        if(((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVectorX[itk2].M() < mass2 - nsigmaside*sigma2 || GoodTrackFourVectorX[itk2].M() > mass2 + nsigmaside*sigma2)) || ((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVectorX[itk2].M() > mass2 - nsigmapeak*sigma2 && GoodTrackFourVectorX[itk2].M() < mass2 + nsigmapeak*sigma2))){
        h3->Fill(x3DX,toteff*weight);
        h3rot->Fill(x3D_rot,toteff*weight);
        h3inv->Fill(x3D_inv,toteff*weight);    
        if(itk1==0 && itk2==1)nevsb->Fill(1);
        }

        if((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma) && (GoodTrackFourVectorX[itk2].M() < mass2 - nsigmaside*sigma2)){hbbL->Fill(x3DX,toteff*weight);}
        
        if((GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVectorX[itk2].M() > mass2 + nsigmaside*sigma2)){hbbR->Fill(x3DX,toteff*weight);}

        if(((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVectorX[itk2].M() < mass2 - nsigmaside*sigma2)) || ((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma) && (GoodTrackFourVectorX[itk2].M() > mass2 - nsigmapeak*sigma2 && GoodTrackFourVectorX[itk2].M() < mass2 + nsigmapeak*sigma2))){hsbL->Fill(x3DX,toteff*weight);}
        
        if(((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVectorX[itk2].M() > mass2 + nsigmaside*sigma2)) || ((GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVectorX[itk2].M() > mass2 - nsigmapeak*sigma2 && GoodTrackFourVectorX[itk2].M() < mass2 + nsigmapeak*sigma2))){hsbR->Fill(x3DX,toteff*weight);}

        
    }}}


void hbt_MM_UU_cross(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVectord1, std::vector<TLorentzVector> GoodTrackFourVectord2, std::vector<TLorentzVector> GoodTrackFourVectorX, std::vector<TLorentzVector> GoodTrackFourVectord1X, std::vector<TLorentzVector> GoodTrackFourVectord2X, Double_t aux_N_tk_offline, double mass, double sigma, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h11, THnSparse* h12, THnSparse* h12x, THnSparse* hX, THnSparse* hX1, THnSparse* hX2, THnSparse* hX12, double weight){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    for(unsigned int itk2=0; itk2<GoodTrackFourVectorX.size();itk2++){
       //V0 qinv
        Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVectorX[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVectorX[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        //daughter1 qinv
        Double_t qX1 = GetQ(GoodTrackFourVectord1[itk1], GoodTrackFourVectord1X[itk2]);
    	TLorentzVector psum2X1 = GoodTrackFourVectord1[itk1] + GoodTrackFourVectord1X[itk2];
        Double_t ktX1=(psum2X1.Pt())/2.;
        Double_t x3DX1[3]={qX1,ktX1,aux_N_tk_offline};
        //daughter2 qinv
        Double_t qX2 = GetQ(GoodTrackFourVectord2[itk1], GoodTrackFourVectord2X[itk2]);
    	TLorentzVector psum2X2 = GoodTrackFourVectord2[itk1] + GoodTrackFourVectord2X[itk2];
        Double_t ktX2=(psum2X2.Pt())/2.;
        Double_t x3DX2[3]={qX2,ktX2,aux_N_tk_offline};
        //daughter12 qinv
        Double_t qX3 = GetQ(GoodTrackFourVectord1[itk1], GoodTrackFourVectord2X[itk2]);
    	TLorentzVector psum2X3 = GoodTrackFourVectord1[itk1] + GoodTrackFourVectord2X[itk2];
        Double_t ktX3=(psum2X3.Pt())/2.;
        Double_t x3DX3[3]={qX3,ktX3,aux_N_tk_offline};
         //daughter21 qinv       
        Double_t qX4 = GetQ(GoodTrackFourVectord2[itk1], GoodTrackFourVectord1X[itk2]);
    	TLorentzVector psum2X4 = GoodTrackFourVectord2[itk1] + GoodTrackFourVectord1X[itk2];
        Double_t ktX4=(psum2X4.Pt())/2.;
        Double_t x3DX4[3]={qX4,ktX4,aux_N_tk_offline};    
        
        if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVectorX[itk2].M() > mass2 - nsigmapeak*sigma2 && GoodTrackFourVectorX[itk2].M() < mass2 + nsigmapeak*sigma2)){
        h1->Fill(x3DX,weight);
        h11->Fill(x3DX1,weight);    
        h12->Fill(x3DX2,weight);   
        h12x->Fill(x3DX3,weight);
        h12x->Fill(x3DX4,weight);
        hX->Fill(x3DX,weight);
        hX1->Fill(x3DX1,weight);
        hX2->Fill(x3DX2,weight);
        hX12->Fill(x3DX3,weight);
        hX12->Fill(x3DX4,weight);
        } 
        
    }}
    
}


//simples

void hbt_simples(std::vector<TLorentzVector> GoodTrackFourVector, Double_t aux_N_tk_offline, double mass, double sigma, double nsigmapeak, THnSparse* h1, double weight){
for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
Double_t ktX=(psum2X.Pt())/2.;
Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma))
{h1->Fill(x3DX,weight);}}}}

void hbt_simples_cross(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVector2, Double_t aux_N_tk_offline, double mass, double sigma, double mass2, double sigma2, double nsigmapeak, THnSparse* h1, double weight){
for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
for(unsigned int itk2=0; itk2<GoodTrackFourVector2.size();itk2++){
Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector2[itk2]);
TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector2[itk2];
Double_t ktX=(psum2X.Pt())/2.;
Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector2[itk2].M() > mass2 - nsigmapeak*sigma2 && GoodTrackFourVector2[itk2].M() < mass2 + nsigmapeak*sigma2))
{h1->Fill(x3DX,weight);}}}}

//----------------------------------------------------------------------------------------------------------------------

//mixing

//----------------------------------------------------------------------------------------------------------------------

void MixEvents(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, THnSparse* hbbL, THnSparse* hbbR, THnSparse* hsbL, THnSparse* hsbR, int p, TH2* effhist, TH1* nevss, TH1* nevbb, TH1* nevsb, std::vector<double> weight){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
   int takeAssociated = 0;
   double weightev1 = weight[nevt1];
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
      if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      takeAssociated++;
      double weightev2 = weight[nevt_assoc];
      double weightprod = weightev1;
      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
      for(int imix=0; imix<nMixmult_nevt1; imix++){
	     double eff1 = 1.0;
 	     if(p==20 || p==21 || p==22)eff1 = getTrkCorrWeight(effhist,ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Pt(),ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Rapidity());
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
		    double eff2 = 1.0;
	        if(p==20 || p==21 || p==22)eff2 = getTrkCorrWeight(effhist,ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Pt(),ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Rapidity());
            double toteff = eff1*eff2;
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){h1->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevss->Fill(1);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)){h2->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevbb->Fill(1);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){h3->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevsb->Fill(1);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma)){hbbL->Fill(xmix3D,toteff*weightprod);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)){hbbR->Fill(xmix3D,toteff*weightprod);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){hsbL->Fill(xmix3D,toteff*weightprod);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){hsbR->Fill(xmix3D,toteff*weightprod);}

}}

}}} //DONE

void MixEvents_eta(int ntrkoff_min, int ntrkoff_max,   std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vec, std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vec, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, THnSparse* hbbL, THnSparse* hbbR, THnSparse* hsbL, THnSparse* hsbR, int p, TH2* effhist, TH1* nevss, TH1* nevbb, TH1* nevsb, std::vector<double> weight){
///sort vectors of maps for each ntrkoffline range
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > aux_ev_GoodTrackFourVector_etaMixWeight_vec;
std::vector<std::pair<Double_t,int> > aux_ev_ntrkoff_etaMixWeight_vec;
std::vector<double> w;
int aux_n_evts = ev_GoodTrackFourVector_etaMixWeight_vec.size();
for(int ievt=0; ievt<aux_n_evts; ievt++){
   if((ev_ntrkoff_etaMixWeight_vec[ievt]).second < ntrkoff_min || ntrkoff_max < (ev_ntrkoff_etaMixWeight_vec[ievt]).second)continue; //only events in a given range
   aux_ev_GoodTrackFourVector_etaMixWeight_vec.push_back(ev_GoodTrackFourVector_etaMixWeight_vec[ievt]);
   aux_ev_ntrkoff_etaMixWeight_vec.push_back(ev_ntrkoff_etaMixWeight_vec[ievt]);
   w.push_back(weight[ievt]);
}
std::sort( aux_ev_GoodTrackFourVector_etaMixWeight_vec.begin(), aux_ev_GoodTrackFourVector_etaMixWeight_vec.end(), etaMixSort );
std::sort( aux_ev_ntrkoff_etaMixWeight_vec.begin(), aux_ev_ntrkoff_etaMixWeight_vec.end(), etaMixSort_ntrkoff );
int aux_n_evts_inNtrkoffRange = aux_ev_GoodTrackFourVector_etaMixWeight_vec.size();
//auxiliar vectors for mixing
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec;
int nMixmult_nevt1;
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec;
int nMixmult_nevt_assoc;
for(int nevt=0; nevt+1<aux_n_evts_inNtrkoffRange; nevt+=2) {
   ev_GoodTrackFourVectorTemp_nevt1_vec= (aux_ev_GoodTrackFourVector_etaMixWeight_vec[nevt]).second;
   nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   ev_GoodTrackFourVectorTemp_nevt_assoc_vec= (aux_ev_GoodTrackFourVector_etaMixWeight_vec[nevt+1]).second;
   nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
   double weightev1 = w[nevt];
   double weightev2 = w[nevt+1];
   double weightprod = weightev1;   
   for(int imix=0; imix<nMixmult_nevt1; imix++){
	     double eff1 = 1.0;
 	     if(p==20 || p==21 || p==22)eff1 = getTrkCorrWeight(effhist,ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Pt(),ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Rapidity());
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
		    double eff2 = 1.0;
	        if(p==20 || p==21 || p==22)eff2 = getTrkCorrWeight(effhist,ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Pt(),ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Rapidity());
            double toteff = eff1*eff2;
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            TLorentzVector psum2 = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
            Double_t kt=(psum2.Pt())/2.;
            Double_t xmix3D[3]={qmix,kt,(Double_t)(aux_ev_ntrkoff_etaMixWeight_vec[nevt]).second};      
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){h1->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevss->Fill(1);}        
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)){h2->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevbb->Fill(1);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){h3->Fill(xmix3D,toteff*weightprod);  if(imix==0 && iimix==0)nevsb->Fill(1);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma)){hbbL->Fill(xmix3D,toteff*weightprod);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)){hbbR->Fill(xmix3D,toteff*weightprod);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){hsbL->Fill(xmix3D,toteff*weightprod);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){hsbR->Fill(xmix3D,toteff*weightprod);}

    }
   }
}//end first for loop
//clear vectors for next ntrkoff range
aux_ev_GoodTrackFourVector_etaMixWeight_vec.clear();
aux_ev_ntrkoff_etaMixWeight_vec.clear();
w.clear();
}//DONE

void MixEvents_cross(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, double mass, double sigma, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector2, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, THnSparse* hbbL, THnSparse* hbbR, THnSparse* hsbL, THnSparse* hsbR, int p, TH2* effhist, TH2* effhist2, TH1* nevss, TH1* nevbb, TH1* nevsb, std::vector<double> weight){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec2= ev_GoodTrackFourVector2[nevt1];	   
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   int nMixmult_nevt12=ev_GoodTrackFourVectorTemp_nevt1_vec2.size();
   if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
   int takeAssociated = 0;
   double weightev1 = weight[nevt1];
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
      if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec2= ev_GoodTrackFourVector2[nevt_assoc];
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      int nMixmult_nevt_assoc2=ev_GoodTrackFourVectorTemp_nevt_assoc_vec2.size();
      double weightev2 = weight[nevt_assoc];
      double weightprod = weightev1;
      takeAssociated++;
      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
      for(int imix=0; imix<nMixmult_nevt1; imix++){
	     double eff1 = 1.0;
 	     if(p==20 || p==21 || p==22)eff1 = getTrkCorrWeight(effhist,ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Pt(),ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Rapidity());
         for(int iimix=0; iimix<nMixmult_nevt_assoc2; iimix++){
		    double eff2 = 1.0;
	        if(p==20 || p==21 || p==22)eff2 = getTrkCorrWeight(effhist2,ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].Pt(),ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].Rapidity());
            double toteff = eff1*eff2;
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2)){h1->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevss->Fill(1);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 - nsigmaside*sigma2 || ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 + nsigmaside*sigma2)){h2->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevbb->Fill(1);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 - nsigmaside*sigma2 || ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 + nsigmaside*sigma2)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2))){h3->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevsb->Fill(1);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 - nsigmaside*sigma2)){hbbL->Fill(xmix3D,toteff*weightprod);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 + nsigmaside*sigma2)){hbbR->Fill(xmix3D,toteff*weightprod);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 - nsigmaside*sigma2)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2))){hsbL->Fill(xmix3D,toteff*weightprod);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 + nsigmaside*sigma2)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2))){hsbR->Fill(xmix3D,toteff*weightprod);}
		}}
      for(int imix=0; imix<nMixmult_nevt12; imix++){
	     double eff1 = 1.0;
 	     if(p==20 || p==21 || p==22)eff1 = getTrkCorrWeight(effhist2,ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].Pt(),ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].Rapidity());
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
		    double eff2 = 1.0;
	        if(p==20 || p==21 || p==22)eff2 = getTrkCorrWeight(effhist,ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Pt(),ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Rapidity());
            double toteff = eff1*eff2;
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec2[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec2[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
            if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){h1->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevss->Fill(1);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 - nsigmaside*sigma2 || ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 + nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)){h2->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevbb->Fill(1);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 - nsigmaside*sigma2 || ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 + nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){h3->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevsb->Fill(1);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 - nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma)){hbbL->Fill(xmix3D,toteff*weightprod);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 + nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)){hbbR->Fill(xmix3D,toteff*weightprod);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 - nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){hsbL->Fill(xmix3D,toteff*weightprod);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 + nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){hsbR->Fill(xmix3D,toteff*weightprod);}
		}}

}}}

void MixEvents_eta_cross(int ntrkoff_min, int ntrkoff_max,   std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vec, double mass, double sigma,  std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vec2, std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vec, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, THnSparse* hbbL, THnSparse* hbbR, THnSparse* hsbL, THnSparse* hsbR, int p, TH2* effhist,TH2* effhist2, TH1* nevss, TH1* nevbb, TH1* nevsb, std::vector<double> weight){
///sort vectors of maps for each ntrkoffline range
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > aux_ev_GoodTrackFourVector_etaMixWeight_vec;
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > aux_ev_GoodTrackFourVector_etaMixWeight_vec2;
std::vector<std::pair<Double_t,int> > aux_ev_ntrkoff_etaMixWeight_vec;
std::vector<double> w;
int aux_n_evts = ev_GoodTrackFourVector_etaMixWeight_vec.size();
for(int ievt=0; ievt<aux_n_evts; ievt++){
   if((ev_ntrkoff_etaMixWeight_vec[ievt]).second < ntrkoff_min || ntrkoff_max < (ev_ntrkoff_etaMixWeight_vec[ievt]).second)continue; //only events in a given range
   aux_ev_GoodTrackFourVector_etaMixWeight_vec.push_back(ev_GoodTrackFourVector_etaMixWeight_vec[ievt]);
   aux_ev_GoodTrackFourVector_etaMixWeight_vec2.push_back(ev_GoodTrackFourVector_etaMixWeight_vec2[ievt]);
   aux_ev_ntrkoff_etaMixWeight_vec.push_back(ev_ntrkoff_etaMixWeight_vec[ievt]);
   w.push_back(weight[ievt]);
}
std::sort( aux_ev_GoodTrackFourVector_etaMixWeight_vec.begin(), aux_ev_GoodTrackFourVector_etaMixWeight_vec.end(), etaMixSort );
std::sort( aux_ev_GoodTrackFourVector_etaMixWeight_vec2.begin(), aux_ev_GoodTrackFourVector_etaMixWeight_vec2.end(), etaMixSort );
std::sort( aux_ev_ntrkoff_etaMixWeight_vec.begin(), aux_ev_ntrkoff_etaMixWeight_vec.end(), etaMixSort_ntrkoff );
int aux_n_evts_inNtrkoffRange = aux_ev_GoodTrackFourVector_etaMixWeight_vec.size();
//auxiliar vectors for mixing
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec;
int nMixmult_nevt1;
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec;
int nMixmult_nevt_assoc;
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec2;
int nMixmult_nevt12;
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec2;
int nMixmult_nevt_assoc2;
for(int nevt=0; nevt+1<aux_n_evts_inNtrkoffRange; nevt+=2) {
   ev_GoodTrackFourVectorTemp_nevt1_vec= (aux_ev_GoodTrackFourVector_etaMixWeight_vec[nevt]).second;
   nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   ev_GoodTrackFourVectorTemp_nevt_assoc_vec= (aux_ev_GoodTrackFourVector_etaMixWeight_vec[nevt+1]).second;
   nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
   ev_GoodTrackFourVectorTemp_nevt1_vec2= (aux_ev_GoodTrackFourVector_etaMixWeight_vec2[nevt]).second;
   nMixmult_nevt12=ev_GoodTrackFourVectorTemp_nevt1_vec2.size();
   ev_GoodTrackFourVectorTemp_nevt_assoc_vec2= (aux_ev_GoodTrackFourVector_etaMixWeight_vec2[nevt+1]).second;
   nMixmult_nevt_assoc2=ev_GoodTrackFourVectorTemp_nevt_assoc_vec2.size();
   double weightev1 = w[nevt];
   double weightev2 = w[nevt+1];
   double weightprod = weightev1;   
      for(int imix=0; imix<nMixmult_nevt1; imix++){
	     double eff1 = 1.0;
 	     if(p==20 || p==21 || p==22)eff1 = getTrkCorrWeight(effhist,ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Pt(),ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Rapidity());
         for(int iimix=0; iimix<nMixmult_nevt_assoc2; iimix++){
		    double eff2 = 1.0;
	        if(p==20 || p==21 || p==22)eff2 = getTrkCorrWeight(effhist2,ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].Pt(),ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].Rapidity());
            double toteff = eff1*eff2;
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)(aux_ev_ntrkoff_etaMixWeight_vec[nevt]).second};
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2)){h1->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevss->Fill(1);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 - nsigmaside*sigma2 || ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 + nsigmaside*sigma2)){h2->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevbb->Fill(1);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 - nsigmaside*sigma2 || ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 + nsigmaside*sigma2)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2))){h3->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevsb->Fill(1);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 - nsigmaside*sigma2)){hbbL->Fill(xmix3D,toteff*weightprod);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 + nsigmaside*sigma2)){hbbR->Fill(xmix3D,toteff*weightprod);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 - nsigmaside*sigma2)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2))){hsbL->Fill(xmix3D,toteff*weightprod);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 + nsigmaside*sigma2)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2))){hsbR->Fill(xmix3D,toteff*weightprod);}
		}}

      for(int imix=0; imix<nMixmult_nevt12; imix++){
	     double eff1 = 1.0;
 	     if(p==20 || p==21 || p==22)eff1 = getTrkCorrWeight(effhist2,ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].Pt(),ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].Rapidity());
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
		    double eff2 = 1.0;
	        if(p==20 || p==21 || p==22)eff2 = getTrkCorrWeight(effhist,ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Pt(),ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Rapidity());
            double toteff = eff1*eff2;
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec2[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec2[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)(aux_ev_ntrkoff_etaMixWeight_vec[nevt]).second};
            if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){h1->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevss->Fill(1);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 - nsigmaside*sigma2 || ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 + nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)){h2->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevbb->Fill(1);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 - nsigmaside*sigma2 || ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 + nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){h3->Fill(xmix3D,toteff*weightprod); if(imix==0 && iimix==0)nevsb->Fill(1);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 - nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma)){hbbL->Fill(xmix3D,toteff*weightprod);}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 + nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)){hbbR->Fill(xmix3D,toteff*weightprod);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 - nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){hsbL->Fill(xmix3D,toteff*weightprod);}
            if(((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 + nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){hsbR->Fill(xmix3D,toteff*weightprod);}
		}}


}//end first for loop

//clear vectors for next ntrkoff range
aux_ev_GoodTrackFourVector_etaMixWeight_vec.clear();
aux_ev_GoodTrackFourVector_etaMixWeight_vec2.clear();
aux_ev_ntrkoff_etaMixWeight_vec.clear();
w.clear();


}//end MixEvents_eta function


//for shape studies
void MixEvents_test(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1, std::vector<double> weight){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
int takeAssociated = 0;
double weightev1 = weight[nevt1];
for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec = ev_GoodTrackFourVector[nevt_assoc];
int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
takeAssociated++;
double weightev2 = weight[nevt_assoc];
double weightprod = weightev1;
if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
for(int imix=0; imix<nMixmult_nevt1; imix++){
for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
Double_t ktmix=(psum2mix.Pt())/2.;
Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){h1->Fill(xmix3D,weightprod);}
}}}}} //DONE

void MixEvents_cross_test(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, double mass, double sigma, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector2, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1, std::vector<double> weight){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec2= ev_GoodTrackFourVector2[nevt1];	   
int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
int nMixmult_nevt12=ev_GoodTrackFourVectorTemp_nevt1_vec2.size();
if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
int takeAssociated = 0;
double weightev1 = weight[nevt1];
for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec2= ev_GoodTrackFourVector2[nevt_assoc];
int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
int nMixmult_nevt_assoc2=ev_GoodTrackFourVectorTemp_nevt_assoc_vec2.size();
takeAssociated++;
double weightev2 = weight[nevt_assoc];
double weightprod = weightev1;
if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.

for(int imix=0; imix<nMixmult_nevt1; imix++){
for(int iimix=0; iimix<nMixmult_nevt_assoc2; iimix++){
Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix]);
TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix];
Double_t ktmix=(psum2mix.Pt())/2.;
Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2)){h1->Fill(xmix3D,weightprod);}
}}

for(int imix=0; imix<nMixmult_nevt12; imix++){
for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec2[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec2[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
Double_t ktmix=(psum2mix.Pt())/2.;
Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){h1->Fill(xmix3D,weightprod);}
}}

}}}


//-----------------------------------------------------------------------------
//call mix
//-----------------------------------------------------------------------------

void call_mix(int m, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, THnSparse* hbbL, THnSparse* hbbR, THnSparse* hsbL, THnSparse* hsbR, int p, TH2* effhist, TH1* nevss, TH1* nevbb, TH1* nevsb, std::vector<double> weight){

cout << "Random Mixing" << endl;

cout << "Nmix = " << nEvt_to_mix << endl;  
    
if(m == 0){ //[0-120]
    
MixEvents(0,4,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(5,9,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(10,14,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(15,19,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(20,24,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(25,29,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(30,34,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(35,39,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(40,44,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(45,49,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(50,54,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(55,59,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(60,64,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(65,69,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(70,74,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(75,79,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(80,84,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(85,89,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(90,94,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(95,99,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(100,104,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(105,109,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(110,114,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(115,119,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
/*
MixEvents(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);

MixEvents(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
*/

}else if(m == 1){
    
MixEvents(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
    
}else if(m == 2){

MixEvents(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);

}else if(m == 3){


MixEvents(185,189,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(190,194,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(195,199,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(200,204,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(205,209,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(210,214,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(215,219,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(220,224,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(225,229,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(230,234,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(235,239,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(240,244,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents(245,249,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    


}else if(m == 4){

MixEvents(250,254,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(255,259,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(260,264,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(265,269,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(270,274,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(275,279,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(280,284,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(285,289,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(290,294,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(295,299,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(300,304,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(305,309,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(310,314,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(315,319,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(320,324,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(325,329,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(330,334,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(335,339,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(340,344,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(345,349,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(350,354,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(355,359,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(360,364,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(365,369,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(370,374,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(375,379,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(380,384,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(385,389,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(390,394,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents(395,399,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents(400,404,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents(405,409,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents(410,414,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents(415,419,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents(420,424,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents(425,429,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents(430,434,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents(435,439,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents(440,444,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents(445,449,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    


}//pPb

} //DONE

void call_mix_eta(int m, std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vec, std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vec, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, THnSparse* hbbL, THnSparse* hbbR, THnSparse* hsbL, THnSparse* hsbR, int p, TH2* effhist, TH1* nevss, TH1* nevbb, TH1* nevsb, std::vector<double> weight){

cout << "Eta Mixing" << endl;
    
if(m == 0){ //[0-120]
    
MixEvents_eta(0,4,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(5,9,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(10,14,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(15,19,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(20,24,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(25,29,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(30,34,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(35,39,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(40,44,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(45,49,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(50,54,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(55,59,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(60,64,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(65,69,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(70,74,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(75,79,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(80,84,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(85,89,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(90,94,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(95,99,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(100,104,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(105,109,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(110,114,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(115,119,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
/*
MixEvents_eta(120,124,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(125,129,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(130,134,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(135,139,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(140,144,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(145,149,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);

MixEvents_eta(150,154,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(155,159,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(160,164,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(165,169,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(170,174,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(175,179,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(180,184,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
*/
}else if(m == 1){
    
MixEvents_eta(120,124,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(125,129,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(130,134,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(135,139,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(140,144,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(145,149,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
    
}else if(m == 2){

MixEvents_eta(150,154,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(155,159,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(160,164,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(165,169,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(170,174,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(175,179,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(180,184,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);

}else if(m == 3){

MixEvents_eta(185,189,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(190,194,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(195,199,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(200,204,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(205,209,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(210,214,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(215,219,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(220,224,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(225,229,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(230,234,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(235,239,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(240,244,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents_eta(245,249,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    

}else if(m == 4){

MixEvents_eta(250,254,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(255,259,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(260,264,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(265,269,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(270,274,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(275,279,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(280,284,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(285,289,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(290,294,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(295,299,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(300,304,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(305,309,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(310,314,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(315,319,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(320,324,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(325,329,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(330,334,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(335,339,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(340,344,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(345,349,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(350,354,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(355,359,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(360,364,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(365,369,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(370,374,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(375,379,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(380,384,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(385,389,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(390,394,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(395,399,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);    
MixEvents_eta(400,404,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(405,409,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(410,414,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(415,419,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(420,424,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(425,429,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(430,434,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(435,439,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(440,444,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);
MixEvents_eta(445,449,ev_GoodTrackFourVector_etaMixWeight_vec,ev_ntrkoff_etaMixWeight_vec,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,nevss,nevbb,nevsb,weight);

}//pPb

}//DONE

void call_mix_cross(int m, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector2, double mass, double sigma, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, THnSparse* hbbL, THnSparse* hbbR, THnSparse* hsbL, THnSparse* hsbR, int p, TH2* effhist, TH2* effhist2, TH1* nevss, TH1* nevbb, TH1* nevsb, std::vector<double> weight){

cout << "Random Mixing" << endl;

cout << "Nmix = " << nEvt_to_mix << endl;
    
if(m == 0){ //[0-120]
    
MixEvents_cross(0,4,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(5,9,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(10,14,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(15,19,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(20,24,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(25,29,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(30,34,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(35,39,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(40,44,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(45,49,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(50,54,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(55,59,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(60,64,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(65,69,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(70,74,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(75,79,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(80,84,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(85,89,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(90,94,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(95,99,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(100,104,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(105,109,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(110,114,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(115,119,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
/*
MixEvents_cross(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);

MixEvents_cross(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
*/
}else if(m == 1){
    
MixEvents_cross(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
    
}else if(m == 2){

MixEvents_cross(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);

}else if(m == 3){

MixEvents_cross(185,189,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(190,194,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(195,199,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(200,204,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(205,209,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(210,214,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(215,219,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(220,224,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(225,229,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(230,234,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(235,239,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(240,244,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);    
MixEvents_cross(245,249,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);    

}else if(m == 4){

MixEvents_cross(250,254,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(255,259,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(260,264,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(265,269,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(270,274,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(275,279,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(280,284,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(285,289,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(290,294,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(295,299,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(300,304,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(305,309,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(310,314,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(315,319,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(320,324,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(325,329,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(330,334,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(335,339,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(340,344,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(345,349,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(350,354,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(355,359,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(360,364,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(365,369,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(370,374,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(375,379,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(380,384,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(385,389,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(390,394,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(395,399,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(400,404,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(405,409,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(410,414,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(415,419,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(420,424,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(425,429,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(430,434,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(435,439,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(440,444,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_cross(445,449,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
    
}//pPb

} //DONE

void call_mix_eta_cross(int m, std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vec,  std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vec2, std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vec, double mass, double sigma, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, THnSparse* hbbL, THnSparse* hbbR, THnSparse* hsbL, THnSparse* hsbR, int p, TH2* effhist, TH2* effhist2, TH1* nevss, TH1* nevbb, TH1* nevsb, std::vector<double> weight){

cout << "Eta Mixing" << endl;
    
if(m == 0){ //[0-120]
    
MixEvents_eta_cross(0,4,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(5,9,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(10,14,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(15,19,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(20,24,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(25,29,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(30,34,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(35,39,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(40,44,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(45,49,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(50,54,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(55,59,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(60,64,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(65,69,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(70,74,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(75,79,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(80,84,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(85,89,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(90,94,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(95,99,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(100,104,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(105,109,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(110,114,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(115,119,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
/*
MixEvents_eta_cross(120,124,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(125,129,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(130,134,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(135,139,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(140,144,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(145,149,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);

MixEvents_eta_cross(150,154,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(155,159,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(160,164,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(165,169,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(170,174,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(175,179,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(180,184,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
*/
}else if(m == 1){
    
MixEvents_eta_cross(120,124,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(125,129,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(130,134,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(135,139,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(140,144,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(145,149,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
    
}else if(m == 2){

MixEvents_eta_cross(150,154,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(155,159,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(160,164,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(165,169,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(170,174,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(175,179,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(180,184,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);

}else if(m == 3){

MixEvents_eta_cross(185,189,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(190,194,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(195,199,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(200,204,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(205,209,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(210,214,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(215,219,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(220,224,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(225,229,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(230,234,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(235,239,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(240,244,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);    
MixEvents_eta_cross(245,249,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);    

}else if(m == 4){

MixEvents_eta_cross(250,254,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(255,259,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(260,264,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(265,269,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(270,274,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(275,279,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(280,284,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(285,289,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(290,294,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(295,299,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(300,304,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(305,309,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(310,314,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(315,319,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(320,324,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(325,329,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(330,334,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(335,339,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(340,344,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(345,349,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(350,354,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(355,359,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(360,364,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(365,369,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(370,374,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(375,379,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(380,384,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(385,389,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(390,394,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(395,399,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(400,404,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(405,409,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(410,414,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(415,419,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(420,424,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(425,429,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(430,434,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(435,439,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(440,444,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);
MixEvents_eta_cross(445,449,ev_GoodTrackFourVector_etaMixWeight_vec,mass,sigma,ev_GoodTrackFourVector_etaMixWeight_vec2,ev_ntrkoff_etaMixWeight_vec,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,hbbL,hbbR,hsbL,hsbR,p,effhist,effhist2,nevss,nevbb,nevsb,weight);

}//pPb

} //DONE


//for feed-down studies

void call_mix_test(int m, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1,std::vector<double> weight){

cout << "Random Mixing" << endl;

cout << "Nmix = " << nEvt_to_mix << endl;  
    
if(m == 0){ //[0-120]
    
MixEvents_test(0,4,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(5,9,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(10,14,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(15,19,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(20,24,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(25,29,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(30,34,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(35,39,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(40,44,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(45,49,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(50,54,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(55,59,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(60,64,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(65,69,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(70,74,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(75,79,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(80,84,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(85,89,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(90,94,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(95,99,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(100,104,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(105,109,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(110,114,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(115,119,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
/*
MixEvents_test(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);

MixEvents_test(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
*/
}else if(m == 1){
    
MixEvents_test(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
    
}else if(m == 2){

MixEvents_test(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);

}else if(m == 3){

MixEvents_test(185,189,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(190,194,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(195,199,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(200,204,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(205,209,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(210,214,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(215,219,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(220,224,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(225,229,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(230,234,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(235,239,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(240,244,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);    
MixEvents_test(245,249,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);    

}else if(m == 4){

MixEvents_test(250,254,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(255,259,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(260,264,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(265,269,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(270,274,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(275,279,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(280,284,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(285,289,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(290,294,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(295,299,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(300,304,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(305,309,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(310,314,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(315,319,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(320,324,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(325,329,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(330,334,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(335,339,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(340,344,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(345,349,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(350,354,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(355,359,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(360,364,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(365,369,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(370,374,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(375,379,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(380,384,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(385,389,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(390,394,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(395,399,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);    
MixEvents_test(400,404,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(405,409,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(410,414,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(415,419,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(420,424,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(425,429,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(430,434,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(435,439,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(440,444,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);
MixEvents_test(445,449,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,weight);

}//pPb

} //DONE

void call_mix_cross_test(int m, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector2, double mass, double sigma, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1,std::vector<double> weight){

cout << "Random Mixing" << endl;

cout << "Nmix = " << nEvt_to_mix << endl;
    
if(m == 0){ //[0-120]
    
MixEvents_cross_test(0,4,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(5,9,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(10,14,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(15,19,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(20,24,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(25,29,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(30,34,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(35,39,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(40,44,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(45,49,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(50,54,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(55,59,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(60,64,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(65,69,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(70,74,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(75,79,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(80,84,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(85,89,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(90,94,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(95,99,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(100,104,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(105,109,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(110,114,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(115,119,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
/*
MixEvents_cross_test(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);

MixEvents_cross_test(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
*/
}else if(m == 1){
    
MixEvents_cross_test(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
    
}else if(m == 2){

MixEvents_cross_test(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);

}else if(m == 3){

MixEvents_cross_test(185,189,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(190,194,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(195,199,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(200,204,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(205,209,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(210,214,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(215,219,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(220,224,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(225,229,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(230,234,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(235,239,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(240,244,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);    
MixEvents_cross_test(245,249,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);    

}else if(m == 4){

MixEvents_cross_test(250,254,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(255,259,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(260,264,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(265,269,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(270,274,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(275,279,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(280,284,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(285,289,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(290,294,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(295,299,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(300,304,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(305,309,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(310,314,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(315,319,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(320,324,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(325,329,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(330,334,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(335,339,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(340,344,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(345,349,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(350,354,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(355,359,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(360,364,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(365,369,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(370,374,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(375,379,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(380,384,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(385,389,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(390,394,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(395,399,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(400,404,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(405,409,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(410,414,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(415,419,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(420,424,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(425,429,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(430,434,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(435,439,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(440,444,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
MixEvents_cross_test(445,449,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,weight);
    
}//pPb

} //DONE

//for shape studies

void MixEvents_testxxx(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
int takeAssociated = 0;
for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec = ev_GoodTrackFourVector[nevt_assoc];
int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
takeAssociated++;
if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
for(int imix=0; imix<nMixmult_nevt1; imix++){
for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
Double_t ktmix=(psum2mix.Pt())/2.;
Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){h1->Fill(xmix3D);}
}}}}} //DONE

void MixEvents_cross_testxxx(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, double mass, double sigma, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector2, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec2= ev_GoodTrackFourVector2[nevt1];	   
int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
int nMixmult_nevt12=ev_GoodTrackFourVectorTemp_nevt1_vec2.size();
if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
int takeAssociated = 0;
for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec2= ev_GoodTrackFourVector2[nevt_assoc];
int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
int nMixmult_nevt_assoc2=ev_GoodTrackFourVectorTemp_nevt_assoc_vec2.size();
takeAssociated++;
if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.

for(int imix=0; imix<nMixmult_nevt1; imix++){
for(int iimix=0; iimix<nMixmult_nevt_assoc2; iimix++){
Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix]);
TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix];
Double_t ktmix=(psum2mix.Pt())/2.;
Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2)){h1->Fill(xmix3D);}
}}

for(int imix=0; imix<nMixmult_nevt12; imix++){
for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec2[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec2[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
Double_t ktmix=(psum2mix.Pt())/2.;
Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){h1->Fill(xmix3D);}
}}

}}}


void call_mix_testxxx(int m, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1){

cout << "Random Mixing" << endl;

cout << "Nmix = " << nEvt_to_mix << endl;  
    
if(m == 0){ //[0-120]
    
MixEvents_testxxx(0,4,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(5,9,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(10,14,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(15,19,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(20,24,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(25,29,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(30,34,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(35,39,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(40,44,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(45,49,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(50,54,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(55,59,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(60,64,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(65,69,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(70,74,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(75,79,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(80,84,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(85,89,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(90,94,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(95,99,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(100,104,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(105,109,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(110,114,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(115,119,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
/*
MixEvents_testxxx(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);

MixEvents_testxxx(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
*/
}else if(m == 1){
    
MixEvents_testxxx(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
    
}else if(m == 2){

MixEvents_testxxx(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);

}else if(m == 3){

MixEvents_testxxx(185,189,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(190,194,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(195,199,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(200,204,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(205,209,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(210,214,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(215,219,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(220,224,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(225,229,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(230,234,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(235,239,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(240,244,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);    
MixEvents_testxxx(245,249,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);    

}else if(m == 4){

MixEvents_testxxx(250,254,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(255,259,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(260,264,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(265,269,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(270,274,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(275,279,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(280,284,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(285,289,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(290,294,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(295,299,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(300,304,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(305,309,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(310,314,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(315,319,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(320,324,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(325,329,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(330,334,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(335,339,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(340,344,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(345,349,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(350,354,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(355,359,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(360,364,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(365,369,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(370,374,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(375,379,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(380,384,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(385,389,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(390,394,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(395,399,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);    
MixEvents_testxxx(400,404,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(405,409,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(410,414,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(415,419,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(420,424,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(425,429,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(430,434,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(435,439,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(440,444,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);
MixEvents_testxxx(445,449,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1);

}//pPb

} //DONE

void call_mix_cross_testxxx(int m, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector2, double mass, double sigma, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1){

cout << "Random Mixing" << endl;

cout << "Nmix = " << nEvt_to_mix << endl;
    
if(m == 0){ //[0-120]
    
MixEvents_cross_testxxx(0,4,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(5,9,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(10,14,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(15,19,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(20,24,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(25,29,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(30,34,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(35,39,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(40,44,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(45,49,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(50,54,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(55,59,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(60,64,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(65,69,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(70,74,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(75,79,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(80,84,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(85,89,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(90,94,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(95,99,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(100,104,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(105,109,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(110,114,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(115,119,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
/*
MixEvents_cross_testxxx(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);

MixEvents_cross_testxxx(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
*/

}else if(m == 1){
    
MixEvents_cross_testxxx(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
    
}else if(m == 2){

MixEvents_cross_testxxx(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);

}else if(m == 3){

MixEvents_cross_testxxx(185,189,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(190,194,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(195,199,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(200,204,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(205,209,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(210,214,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(215,219,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(220,224,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(225,229,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(230,234,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(235,239,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(240,244,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);    
MixEvents_cross_testxxx(245,249,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);    

}else if(m == 4){

MixEvents_cross_testxxx(250,254,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(255,259,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(260,264,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(265,269,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(270,274,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(275,279,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(280,284,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(285,289,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(290,294,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(295,299,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(300,304,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(305,309,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(310,314,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(315,319,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(320,324,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(325,329,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(330,334,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(335,339,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(340,344,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(345,349,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(350,354,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(355,359,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(360,364,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(365,369,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(370,374,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(375,379,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(380,384,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(385,389,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(390,394,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(395,399,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(400,404,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(405,409,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(410,414,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(415,419,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(420,424,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(425,429,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(430,434,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(435,439,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(440,444,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
MixEvents_cross_testxxx(445,449,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1);
    
}//pPb

} //DONE


//until here is done


//H-Dibaryon


void hdibaryon(std::vector<TLorentzVector> GoodTrackFourVector, Double_t aux_N_tk_offline, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h1_rot, THnSparse* h1_inv, THnSparse* h2, THnSparse* h2_rot, THnSparse* h2_inv, THnSparse* h3, THnSparse* h3_rot, THnSparse* h3_inv, double weight){
	for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
		for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){

			//Default
			TLorentzVector psum = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
			Double_t mass_hdibaryon=psum.M();
			Double_t pt_hdibaryon=psum.Pt();
            
			Double_t y_hdibaryon=psum.Rapidity();  
            if(fabs(y_hdibaryon) > 1.0)continue;

            Double_t x3D[3]={mass_hdibaryon,pt_hdibaryon,aux_N_tk_offline};

	        //Rotated XY
	        TLorentzVector psum_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector[itk2]);
			Double_t mass_hdibaryon_rot=psum_rot.M();
			Double_t pt_hdibaryon_rot=psum_rot.Pt();
	        Double_t x3D_rot[3]={mass_hdibaryon_rot,pt_hdibaryon_rot,(Double_t)aux_N_tk_offline};
       
	        //Inverted PxPyPz
	        TLorentzVector psum_inv = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector[itk2]);
			Double_t mass_hdibaryon_inv=psum_inv.M();
			Double_t pt_hdibaryon_inv=psum_inv.Pt();
	        Double_t x3D_inv[3]={mass_hdibaryon_inv,pt_hdibaryon_inv,(Double_t)aux_N_tk_offline};

			if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma)){
				h1->Fill(x3D,weight);
				h1_rot->Fill(x3D_rot,weight);
				h1_inv->Fill(x3D_inv,weight);
			}	
        	if((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk2].M() > mass + nsigmaside*sigma)){
				h2->Fill(x3D,weight);
				h2_rot->Fill(x3D_rot,weight);
				h2_inv->Fill(x3D_inv,weight);
			}	
	        if(((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk2].M() > mass + nsigmaside*sigma)) || ((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma))){
 		       h3->Fill(x3D,weight);
	           h3_rot->Fill(x3D_rot,weight);
			   h3_inv->Fill(x3D_inv,weight);    
	        } //sb
		}
	}
}

void hdibaryon_OS(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVector2, Double_t aux_N_tk_offline, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, double weight){
	for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
		for(unsigned int itk2=0; itk2<GoodTrackFourVector2.size();itk2++){
			TLorentzVector psum = GoodTrackFourVector[itk1] + GoodTrackFourVector2[itk2];
			Double_t mass_hdibaryon=psum.M();
			Double_t pt_hdibaryon=psum.Pt();
            
			Double_t y_hdibaryon=psum.Rapidity();  
            if(fabs(y_hdibaryon) > 1.0)continue;


			Double_t x3D[3]={mass_hdibaryon,pt_hdibaryon,aux_N_tk_offline};
			if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector2[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector2[itk2].M() < mass + nsigmapeak*sigma)){h1->Fill(x3D,weight);}
        	if((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk2].M() > mass + nsigmaside*sigma)){h2->Fill(x3D,weight);}	
	        if(((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk2].M() > mass + nsigmaside*sigma)) || ((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma))){h3->Fill(x3D,weight);}
		}
	}
}

void hdibaryon_OS2(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVector2, Double_t aux_N_tk_offline, double mass, double sigma, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h1_rot, THnSparse* h1_inv, THnSparse* h2, THnSparse* h2_rot, THnSparse* h2_inv, THnSparse* h3, THnSparse* h3_rot, THnSparse* h3_inv, double weight){
	for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
		for(unsigned int itk2=0; itk2<GoodTrackFourVector2.size();itk2++){

			TLorentzVector psum = GoodTrackFourVector[itk1] + GoodTrackFourVector2[itk2];
			Double_t mass_hdibaryon=psum.M();
			Double_t pt_hdibaryon=psum.Pt();
            
			Double_t y_hdibaryon=psum.Rapidity();  
            if(fabs(y_hdibaryon) > 1.0)continue;            

			Double_t x3D[3]={mass_hdibaryon,pt_hdibaryon,aux_N_tk_offline};

	        //Rotated XY
	        TLorentzVector psum_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector[itk2]);
			Double_t mass_hdibaryon_rot=psum_rot.M();
			Double_t pt_hdibaryon_rot=psum_rot.Pt();
	        Double_t x3D_rot[3]={mass_hdibaryon_rot,pt_hdibaryon_rot,(Double_t)aux_N_tk_offline};
       
	        //Inverted PxPyPz
	        TLorentzVector psum_inv = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector[itk2]);
			Double_t mass_hdibaryon_inv=psum_inv.M();
			Double_t pt_hdibaryon_inv=psum_inv.Pt();
	        Double_t x3D_inv[3]={mass_hdibaryon_inv,pt_hdibaryon_inv,(Double_t)aux_N_tk_offline};


			if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector2[itk2].M() > mass2 - nsigmapeak*sigma2 && GoodTrackFourVector2[itk2].M() < mass2 + nsigmapeak*sigma2)){h1->Fill(x3D,weight); h1_rot->Fill(x3D_rot,weight); h1_inv->Fill(x3D_inv,weight);}
        	if((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() < mass2 - nsigmaside*sigma2 || GoodTrackFourVector[itk2].M() > mass2 + nsigmaside*sigma2)){h2->Fill(x3D,weight); h2_rot->Fill(x3D_rot,weight); h2_inv->Fill(x3D_inv,weight);}	
	        if(((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() < mass2 - nsigmaside*sigma2 || GoodTrackFourVector[itk2].M() > mass2 + nsigmaside*sigma2)) || ((GoodTrackFourVector[itk1].M() < mass - nsigmaside*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmaside*sigma) && (GoodTrackFourVector[itk2].M() > mass2 - nsigmapeak*sigma2 && GoodTrackFourVector[itk2].M() < mass2 + nsigmapeak*sigma2))){h3->Fill(x3D,weight); h3_rot->Fill(x3D_rot,weight); h3_inv->Fill(x3D_inv,weight);}
		}
	}
}

void MixEvents_hdibaryon(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, std::vector<double> weight){
	int aux_n_evts = (int)ev_ntrkoff.size();
	for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
	   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
	   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
	   if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
	   int takeAssociated = 0;
	   double weightev1 = weight[nevt1];
	   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
	      if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
	      if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
	      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
	      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
	      takeAssociated++;
          double weightev2 = weight[nevt_assoc];
          double weightprod = weightev1;
	      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
	      for(int imix=0; imix<nMixmult_nevt1; imix++){
	         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
	            TLorentzVector psummix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
				Double_t mass_hdibaryon=psummix.M();
				Double_t pt_hdibaryon=psummix.Pt();
            
                Double_t y_hdibaryon=psummix.Rapidity();  
                if(fabs(y_hdibaryon) > 1.0)continue;

	            Double_t xmix3D[3]={mass_hdibaryon,pt_hdibaryon,(Double_t)ev_ntrkoff[nevt1]};
	            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){h1->Fill(xmix3D,weightprod);}
	        	if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)){h2->Fill(xmix3D,weightprod);}	
		        if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){h3->Fill(xmix3D,weightprod);}
	         }
	      }
	   }
	}
} //DONE

void call_mix_hdibaryon(int m, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, double mass, double sigma, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, std::vector<double> weight){

cout << "Random Mixing for H-Dibaryon" << endl;

cout << "Nmix = " << nEvt_to_mix << endl;  
    
if(m == 0){ //[0-120]
    
MixEvents_hdibaryon(0,4,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(5,9,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(10,14,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(15,19,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(20,24,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(25,29,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(30,34,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(35,39,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(40,44,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(45,49,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(50,54,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(55,59,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(60,64,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(65,69,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(70,74,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(75,79,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(80,84,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(85,89,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(90,94,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(95,99,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(100,104,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(105,109,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(110,114,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(115,119,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
/*
MixEvents_hdibaryon(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);

MixEvents_hdibaryon(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
*/
}else if(m == 1){
    
MixEvents_hdibaryon(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
    
}else if(m == 2){

MixEvents_hdibaryon(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);

}else if(m == 3){

MixEvents_hdibaryon(185,189,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(190,194,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(195,199,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(200,204,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(205,209,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(210,214,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(215,219,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(220,224,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(225,229,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(230,234,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(235,239,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(240,244,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);    
MixEvents_hdibaryon(245,249,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);    

}else if(m == 4){

MixEvents_hdibaryon(250,254,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(255,259,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(260,264,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(265,269,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(270,274,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(275,279,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(280,284,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(285,289,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(290,294,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(295,299,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(300,304,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(305,309,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(310,314,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(315,319,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(320,324,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(325,329,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(330,334,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(335,339,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(340,344,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(345,349,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(350,354,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(355,359,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(360,364,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(365,369,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(370,374,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(375,379,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(380,384,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(385,389,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(390,394,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(395,399,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);    
MixEvents_hdibaryon(400,404,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(405,409,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(410,414,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(415,419,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(420,424,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(425,429,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(430,434,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(435,439,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(440,444,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon(445,449,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,nsigmapeak,nsigmaside,h1,h2,h3,weight);

}//pPb

} //DONE

void MixEvents_hdibaryon_OS(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, double mass, double sigma, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector2, double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, std::vector<double> weight){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec2= ev_GoodTrackFourVector2[nevt1];	   
int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
int nMixmult_nevt12=ev_GoodTrackFourVectorTemp_nevt1_vec2.size();
if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
int takeAssociated = 0;
double weightev1 = weight[nevt1];
for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec2= ev_GoodTrackFourVector2[nevt_assoc];
int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
int nMixmult_nevt_assoc2=ev_GoodTrackFourVectorTemp_nevt_assoc_vec2.size();
takeAssociated++;
double weightev2 = weight[nevt_assoc];
double weightprod = weightev1;
if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.

for(int imix=0; imix<nMixmult_nevt1; imix++){
for(int iimix=0; iimix<nMixmult_nevt_assoc2; iimix++){
TLorentzVector psummix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix];
Double_t mass_hdibaryon=psummix.M();
Double_t pt_hdibaryon=psummix.Pt();
            
Double_t y_hdibaryon=psummix.Rapidity();  
if(fabs(y_hdibaryon) > 1.0)continue;

Double_t xmix3D[3]={mass_hdibaryon,pt_hdibaryon,(Double_t)ev_ntrkoff[nevt1]};
if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2)){h1->Fill(xmix3D,weightprod);}
if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 - nsigmaside*sigma2 || ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 + nsigmaside*sigma2)){h2->Fill(xmix3D,weightprod);}
if(((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 - nsigmaside*sigma2 || ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 + nsigmaside*sigma2)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmaside*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2))){h3->Fill(xmix3D,weightprod);}
}}

for(int imix=0; imix<nMixmult_nevt12; imix++){
for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec2[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec2[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
Double_t ktmix=(psum2mix.Pt())/2.;
Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){h1->Fill(xmix3D,weightprod);}
if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 - nsigmaside*sigma2 || ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 + nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)){h2->Fill(xmix3D,weightprod);}
if(((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmaside*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmaside*sigma)) || ((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 - nsigmaside*sigma2 || ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 + nsigmaside*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma))){h3->Fill(xmix3D,weightprod);}
}}

}}}

void call_mix_hdibaryon_OS(int m, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector2, double mass, double sigma,  double mass2, double sigma2, double nsigmapeak, double nsigmaside, THnSparse* h1, THnSparse* h2, THnSparse* h3, std::vector<double> weight){

cout << "Random Mixing for H-Dibaryon" << endl;

cout << "Nmix = " << nEvt_to_mix << endl;  
    
if(m == 0){ //[0-120]
    
MixEvents_hdibaryon_OS(0,4,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(5,9,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(10,14,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(15,19,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(20,24,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(25,29,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(30,34,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(35,39,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(40,44,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(45,49,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(50,54,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(55,59,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(60,64,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(65,69,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(70,74,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(75,79,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(80,84,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(85,89,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(90,94,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(95,99,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(100,104,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(105,109,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(110,114,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(115,119,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
/*
MixEvents_hdibaryon_OS(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);

MixEvents_hdibaryon_OS(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
*/
}else if(m == 1){
    
MixEvents_hdibaryon_OS(120,124,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(125,129,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(130,134,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(135,139,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(140,144,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(145,149,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
    
}else if(m == 2){

MixEvents_hdibaryon_OS(150,154,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(155,159,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(160,164,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(165,169,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(170,174,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(175,179,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(180,184,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);

}else if(m == 3){

MixEvents_hdibaryon_OS(185,189,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(190,194,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(195,199,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(200,204,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(205,209,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(210,214,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(215,219,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(220,224,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(225,229,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(230,234,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(235,239,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(240,244,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);    
MixEvents_hdibaryon_OS(245,249,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);    

}else if(m == 4){

MixEvents_hdibaryon_OS(250,254,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(255,259,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(260,264,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(265,269,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(270,274,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(275,279,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(280,284,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(285,289,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(290,294,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(295,299,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(300,304,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(305,309,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(310,314,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(315,319,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(320,324,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(325,329,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(330,334,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(335,339,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(340,344,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(345,349,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(350,354,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(355,359,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(360,364,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(365,369,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(370,374,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(375,379,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(380,384,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(385,389,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(390,394,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(395,399,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);    
MixEvents_hdibaryon_OS(400,404,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(405,409,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(410,414,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(415,419,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(420,424,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(425,429,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(430,434,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(435,439,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(440,444,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);
MixEvents_hdibaryon_OS(445,449,nEvt_to_mix,ev_ntrkoff,pv_vtx_z,vzcut,ev_GoodTrackFourVector,mass,sigma,ev_GoodTrackFourVector2,mass2,sigma2,nsigmapeak,nsigmaside,h1,h2,h3,weight);

}//pPb

} //DONE


//after here is done

//jets
//for jet+jet 

void hbt_gen_jet(std::vector<TLorentzVector> GoodTrackFourVector, Double_t aux_N_tk_offline, THnSparse* h1, THnSparse* h1inv, THnSparse* h1rot){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
       Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        
        //Rotated XY
        Double_t q_rot = GetQ(GoodTrackFourVector[itk1], InvertXYVector(GoodTrackFourVector[itk2]));        
        TLorentzVector psum2_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector[itk2]);
        Double_t kt_rot = (psum2_rot.Pt())/2.;
        Double_t x3D_rot[3]={q_rot,kt_rot,(Double_t)aux_N_tk_offline};
       
        //Inverted PxPyPz
        Double_t q_inv = GetQ(GoodTrackFourVector[itk1], InvertPVector(GoodTrackFourVector[itk2]));
        TLorentzVector psum2_inv = GoodTrackFourVector[itk1] + InvertPVector(GoodTrackFourVector[itk2]);
        Double_t kt_inv = (psum2_inv.Pt())/2.;
        Double_t x3D_inv[3]={q_inv,kt_inv,(Double_t)aux_N_tk_offline};          

        h1->Fill(x3DX);
        h1rot->Fill(x3D_rot);
        h1inv->Fill(x3D_inv);
        
}}}

void MixEventsJet(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, THnSparse* h1){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
   int takeAssociated = 0;
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
      if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      takeAssociated++;
      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
      for(int imix=0; imix<nMixmult_nevt1; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
            h1->Fill(xmix3D);
}}}}}

void call_mix_Jet(std::vector<int> ev_ntrkoff_vec, std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec, int Evt_mix, THnSparse *h1, std::vector<double> pv_vtx_z, double vzcut){

cout << "Random Mixing " << endl;

cout << "RECO Level" << endl;

cout << "Default: Nmix = " << Evt_mix << endl;

MixEventsJet(0,4,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(5,9,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(10,14,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(15,19,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(20,24,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(25,29,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(30,34,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(35,39,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(40,44,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(45,49,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(50,54,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(55,59,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(60,64,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(65,69,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(70,74,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(75,79,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(80,84,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(85,89,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(90,94,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(95,99,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(100,104,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(105,109,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(110,114,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsJet(115,119,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);

}


//for V0+Jet

void hbt_gen_jet_OS(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVector2, Double_t aux_N_tk_offline, THnSparse* h1, THnSparse* h1inv, THnSparse* h1rot, double mass, double sigma, double nsigmapeak){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    for(unsigned int itk2=0; itk2<GoodTrackFourVector2.size();itk2++){
       Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector2[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector2[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        
        //Rotated XY
        Double_t q_rot = GetQ(GoodTrackFourVector[itk1], InvertXYVector(GoodTrackFourVector2[itk2]));        
        TLorentzVector psum2_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector2[itk2]);
        Double_t kt_rot = (psum2_rot.Pt())/2.;
        Double_t x3D_rot[3]={q_rot,kt_rot,(Double_t)aux_N_tk_offline};
       
        //Inverted PxPyPz
        Double_t q_inv = GetQ(GoodTrackFourVector[itk1], InvertPVector(GoodTrackFourVector2[itk2]));
        TLorentzVector psum2_inv = GoodTrackFourVector[itk1] + InvertPVector(GoodTrackFourVector2[itk2]);
        Double_t kt_inv = (psum2_inv.Pt())/2.;
        Double_t x3D_inv[3]={q_inv,kt_inv,(Double_t)aux_N_tk_offline};          

		if(GoodTrackFourVector[itk1].M() < mass - nsigmapeak*sigma || GoodTrackFourVector[itk1].M() > mass + nsigmapeak*sigma) continue;

        h1->Fill(x3DX);
        h1rot->Fill(x3D_rot);
        h1inv->Fill(x3D_inv);
        
}}}

void MixEventsgen_jet_OS(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector2, THnSparse* h1, double mass, double sigma, double nsigmapeak){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec2= ev_GoodTrackFourVector2[nevt1];	   
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   int nMixmult_nevt12=ev_GoodTrackFourVectorTemp_nevt1_vec2.size();
   if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
   int takeAssociated = 0;
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
      if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec2= ev_GoodTrackFourVector2[nevt_assoc];
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      int nMixmult_nevt_assoc2=ev_GoodTrackFourVectorTemp_nevt_assoc_vec2.size();
      takeAssociated++;
      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
      for(int imix=0; imix<nMixmult_nevt1; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc2; iimix++){
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass - nsigmapeak*sigma || ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass + nsigmapeak*sigma) ) continue;
            h1->Fill(xmix3D);
	  }}
      for(int imix=0; imix<nMixmult_nevt12; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec2[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec2[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
            if((ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass - nsigmapeak*sigma || ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass + nsigmapeak*sigma) ) continue;
            h1->Fill(xmix3D);
	  }}

}}}

void call_mix_gen_OS_basics(std::vector<int> ev_ntrkoff_vec, std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec, std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec2, int Evt_mix, THnSparse *h1, std::vector<double> pv_vtx_z, double vzcut, double mass, double sigma, double nsigmapeak){

cout << "Random Mixing " << endl;

cout << "RECO Level" << endl;

cout << "Default: Nmix = " << Evt_mix << endl;

MixEventsgen_jet_OS(0,4,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(5,9,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(10,14,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(15,19,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(20,24,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(25,29,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(30,34,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(35,39,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(40,44,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(45,49,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(50,54,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(55,59,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(60,64,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(65,69,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(70,74,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(75,79,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(80,84,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(85,89,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(90,94,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(95,99,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(100,104,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(105,109,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(110,114,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);
MixEventsgen_jet_OS(115,119,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,nsigmapeak);

}


//for V0NotaJet+V0NotaJet

void hbt_gen_V0jet(std::vector<TLorentzVector> GoodTrackFourVector, Double_t aux_N_tk_offline, THnSparse* h1, THnSparse* h1inv, THnSparse* h1rot, double mass, double sigma, double nsigmapeak){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
       Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        
        //Rotated XY
        Double_t q_rot = GetQ(GoodTrackFourVector[itk1], InvertXYVector(GoodTrackFourVector[itk2]));        
        TLorentzVector psum2_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector[itk2]);
        Double_t kt_rot = (psum2_rot.Pt())/2.;
        Double_t x3D_rot[3]={q_rot,kt_rot,(Double_t)aux_N_tk_offline};
       
        //Inverted PxPyPz
        Double_t q_inv = GetQ(GoodTrackFourVector[itk1], InvertPVector(GoodTrackFourVector[itk2]));
        TLorentzVector psum2_inv = GoodTrackFourVector[itk1] + InvertPVector(GoodTrackFourVector[itk2]);
        Double_t kt_inv = (psum2_inv.Pt())/2.;
        Double_t x3D_inv[3]={q_inv,kt_inv,(Double_t)aux_N_tk_offline};          

		if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma)){
        h1->Fill(x3DX);
        h1rot->Fill(x3D_rot);
        h1inv->Fill(x3D_inv);
		}        
}}}

void MixEventsgen_V0jet(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, THnSparse* h1, double mass, double sigma, double nsigmapeak){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
   int takeAssociated = 0;
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
      if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      takeAssociated++;
      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
      for(int imix=0; imix<nMixmult_nevt1; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){
            h1->Fill(xmix3D);
            }
}}}}}

void call_mix_gen_V0jet(std::vector<int> ev_ntrkoff_vec, std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec, int Evt_mix, THnSparse *h1, std::vector<double> pv_vtx_z, double vzcut, double mass, double sigma, double nsigmapeak){

cout << "Random Mixing " << endl;

cout << "RECO Level" << endl;

cout << "Default: Nmix = " << Evt_mix << endl;

MixEventsgen_V0jet(0,4,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(5,9,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(10,14,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(15,19,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(20,24,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(25,29,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(30,34,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(35,39,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(40,44,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(45,49,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(50,54,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(55,59,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(60,64,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(65,69,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(70,74,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(75,79,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(80,84,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(85,89,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(90,94,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(95,99,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(100,104,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(105,109,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(110,114,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);
MixEventsgen_V0jet(115,119,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1,mass,sigma,nsigmapeak);

}



//For V0NotaJet+V0Jet and V0NotaJet+V0NotaJet

void hbt_gen_V0jet_OS(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<TLorentzVector> GoodTrackFourVector2, Double_t aux_N_tk_offline, THnSparse* h1, THnSparse* h1inv, THnSparse* h1rot, double mass, double sigma, double mass2, double sigma2, double nsigmapeak){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    for(unsigned int itk2=0; itk2<GoodTrackFourVector2.size();itk2++){
       Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector2[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector2[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        
        //Rotated XY
        Double_t q_rot = GetQ(GoodTrackFourVector[itk1], InvertXYVector(GoodTrackFourVector2[itk2]));        
        TLorentzVector psum2_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector2[itk2]);
        Double_t kt_rot = (psum2_rot.Pt())/2.;
        Double_t x3D_rot[3]={q_rot,kt_rot,(Double_t)aux_N_tk_offline};
       
        //Inverted PxPyPz
        Double_t q_inv = GetQ(GoodTrackFourVector[itk1], InvertPVector(GoodTrackFourVector2[itk2]));
        TLorentzVector psum2_inv = GoodTrackFourVector[itk1] + InvertPVector(GoodTrackFourVector2[itk2]);
        Double_t kt_inv = (psum2_inv.Pt())/2.;
        Double_t x3D_inv[3]={q_inv,kt_inv,(Double_t)aux_N_tk_offline};          

        if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector2[itk2].M() > mass2 - nsigmapeak*sigma2 && GoodTrackFourVector2[itk2].M() < mass2 + nsigmapeak*sigma2)){

        h1->Fill(x3DX);
        h1rot->Fill(x3D_rot);
        h1inv->Fill(x3D_inv);
        
		}
}}}

void MixEventsgen_V0jet_OS(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector2, THnSparse* h1, double mass, double sigma, double mass2, double sigma2, double nsigmapeak){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec2= ev_GoodTrackFourVector2[nevt1];	   
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   int nMixmult_nevt12=ev_GoodTrackFourVectorTemp_nevt1_vec2.size();
   if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
   int takeAssociated = 0;
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
      if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec2= ev_GoodTrackFourVector2[nevt_assoc];
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      int nMixmult_nevt_assoc2=ev_GoodTrackFourVectorTemp_nevt_assoc_vec2.size();
      takeAssociated++;
      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
      for(int imix=0; imix<nMixmult_nevt1; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc2; iimix++){
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
	        if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2)){
            h1->Fill(xmix3D);
            }
	  }}
      for(int imix=0; imix<nMixmult_nevt12; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec2[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec2[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
	        if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){
            h1->Fill(xmix3D);
			}
	  }}

}}}

void call_mix_gen_V0jet_OS(std::vector<int> ev_ntrkoff_vec, std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec, std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec2, int Evt_mix, THnSparse *h1, std::vector<double> pv_vtx_z, double vzcut, double mass, double sigma, double mass2, double sigma2, double nsigmapeak){

cout << "Random Mixing " << endl;

cout << "RECO Level" << endl;

cout << "Default: Nmix = " << Evt_mix << endl;

MixEventsgen_V0jet_OS(0,4,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(5,9,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(10,14,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(15,19,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(20,24,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(25,29,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(30,34,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(35,39,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(40,44,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(45,49,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(50,54,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(55,59,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(60,64,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(65,69,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(70,74,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(75,79,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(80,84,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(85,89,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(90,94,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(95,99,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(100,104,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(105,109,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(110,114,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_V0jet_OS(115,119,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1,mass,sigma,mass2,sigma2,nsigmapeak);

}


//For VJet+V0Jet

void hbt_gen_jet_V0(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<Int_t> NJet, Double_t aux_N_tk_offline, THnSparse* h1, THnSparse* h1inv, THnSparse* h1rot, bool samejet, double mass, double sigma, double nsigmapeak){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
        Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        
        //Rotated XY
        Double_t q_rot = GetQ(GoodTrackFourVector[itk1], InvertXYVector(GoodTrackFourVector[itk2]));        
        TLorentzVector psum2_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector[itk2]);
        Double_t kt_rot = (psum2_rot.Pt())/2.;
        Double_t x3D_rot[3]={q_rot,kt_rot,(Double_t)aux_N_tk_offline};
       
        //Inverted PxPyPz
        Double_t q_inv = GetQ(GoodTrackFourVector[itk1], InvertPVector(GoodTrackFourVector[itk2]));
        TLorentzVector psum2_inv = GoodTrackFourVector[itk1] + InvertPVector(GoodTrackFourVector[itk2]);
        Double_t kt_inv = (psum2_inv.Pt())/2.;
        Double_t x3D_inv[3]={q_inv,kt_inv,(Double_t)aux_N_tk_offline};          

		if(NJet[itk1] < 0 || NJet[itk2] < 0)continue;

		if(samejet){if(NJet[itk1]!=NJet[itk2])continue;}else{if(NJet[itk1]==NJet[itk2])continue;}

		if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVector[itk2].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk2].M() < mass + nsigmapeak*sigma)){

        h1->Fill(x3DX);
        h1rot->Fill(x3D_rot);
        h1inv->Fill(x3D_inv);

		}
        
}}}

void MixEventsgen_jet_V0(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, std::vector<std::vector<Int_t>> NJet, THnSparse* h1, bool samejet, double mass, double sigma, double nsigmapeak){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
   std::vector<Int_t>  ev_NJet_nevt1_vec = NJet[nevt1];
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
   int takeAssociated = 0;
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
      if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
      std::vector<Int_t> ev_NJet_nevt_assoc_vec = NJet[nevt_assoc];
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      takeAssociated++;
      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
      for(int imix=0; imix<nMixmult_nevt1; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};

	  		if(ev_NJet_nevt1_vec[imix] < 0 || ev_NJet_nevt_assoc_vec[iimix] < 0)continue;
			if(samejet){if(ev_NJet_nevt1_vec[imix]!=ev_NJet_nevt_assoc_vec[iimix])continue;}else{if(ev_NJet_nevt1_vec[imix]==ev_NJet_nevt_assoc_vec[iimix])continue;}
            if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){
            h1->Fill(xmix3D);
            }
}}}}}

void call_mix_gen_V0_jet(std::vector<int> ev_ntrkoff_vec, std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec, std::vector<std::vector<Int_t>> NJet, int Evt_mix, THnSparse *h1, std::vector<double> pv_vtx_z, double vzcut, bool samejet, double mass, double sigma, double nsigmapeak){

cout << "Random Mixing" << endl;

cout << "RECO Level" << endl;

cout << "Default: Nmix = " << Evt_mix << endl;

MixEventsgen_jet_V0(0,4,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(5,9,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(10,14,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(15,19,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(20,24,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(25,29,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(30,34,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(35,39,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(40,44,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(45,49,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(50,54,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(55,59,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(60,64,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(65,69,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(70,74,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(75,79,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(80,84,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(85,89,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(90,94,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(95,99,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(100,104,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(105,109,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(110,114,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);
MixEventsgen_jet_V0(115,119,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,h1,samejet,mass,sigma,nsigmapeak);

}

//OS

void hbt_gen_jet_V0_OS(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<Int_t> NJet, std::vector<TLorentzVector> GoodTrackFourVectorX, std::vector<Int_t> NJetX, Double_t aux_N_tk_offline, THnSparse* h1, THnSparse* h1inv, THnSparse* h1rot, bool samejet, double mass, double sigma, double mass2, double sigma2, double nsigmapeak){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    for(unsigned int itk2=0; itk2<GoodTrackFourVectorX.size();itk2++){
        Double_t qX = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVectorX[itk2]);
    	TLorentzVector psum2X = GoodTrackFourVector[itk1] + GoodTrackFourVectorX[itk2];
        Double_t ktX=(psum2X.Pt())/2.;
        Double_t x3DX[3]={qX,ktX,aux_N_tk_offline};
        
        //Rotated XY
        Double_t q_rot = GetQ(GoodTrackFourVector[itk1], InvertXYVector(GoodTrackFourVectorX[itk2]));        
        TLorentzVector psum2_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVectorX[itk2]);
        Double_t kt_rot = (psum2_rot.Pt())/2.;
        Double_t x3D_rot[3]={q_rot,kt_rot,(Double_t)aux_N_tk_offline};
       
        //Inverted PxPyPz
        Double_t q_inv = GetQ(GoodTrackFourVector[itk1], InvertPVector(GoodTrackFourVectorX[itk2]));
        TLorentzVector psum2_inv = GoodTrackFourVector[itk1] + InvertPVector(GoodTrackFourVectorX[itk2]);
        Double_t kt_inv = (psum2_inv.Pt())/2.;
        Double_t x3D_inv[3]={q_inv,kt_inv,(Double_t)aux_N_tk_offline};          

		if(NJet[itk1] < 0 || NJetX[itk2] < 0)continue;

		if(samejet){if(NJet[itk1]!=NJetX[itk2])continue;}else{if(NJet[itk1]==NJetX[itk2])continue;}

        if((GoodTrackFourVector[itk1].M() > mass - nsigmapeak*sigma && GoodTrackFourVector[itk1].M() < mass + nsigmapeak*sigma) && (GoodTrackFourVectorX[itk2].M() > mass2 - nsigmapeak*sigma2 && GoodTrackFourVectorX[itk2].M() < mass2 + nsigmapeak*sigma2)){

        h1->Fill(x3DX);
        h1rot->Fill(x3D_rot);
        h1inv->Fill(x3D_inv);

		}
        
}}}

void MixEventsgen_jet_V0_OS(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, std::vector<std::vector<Int_t>> NJet, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector2, std::vector<std::vector<Int_t>> NJet2, THnSparse* h1, bool samejet, double mass, double sigma, double mass2, double sigma2, double nsigmapeak){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec2= ev_GoodTrackFourVector2[nevt1];	   
   std::vector<Int_t> ev_NJet_nevt1_vec = NJet[nevt1];
   std::vector<Int_t> ev_NJet_nevt1_vec2 = NJet2[nevt1];
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   int nMixmult_nevt12=ev_GoodTrackFourVectorTemp_nevt1_vec2.size();
   if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
   int takeAssociated = 0;
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
      if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec2= ev_GoodTrackFourVector2[nevt_assoc];
      std::vector<Int_t> ev_NJet_nevt_assoc_vec = NJet[nevt_assoc];
      std::vector<Int_t> ev_NJet_nevt_assoc_vec2 = NJet2[nevt_assoc];
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      int nMixmult_nevt_assoc2=ev_GoodTrackFourVectorTemp_nevt_assoc_vec2.size();
      takeAssociated++;
      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
      for(int imix=0; imix<nMixmult_nevt1; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc2; iimix++){
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
	  		if(ev_NJet_nevt1_vec[imix] < 0 || ev_NJet_nevt_assoc_vec2[iimix] < 0)continue;
			if(samejet){if(ev_NJet_nevt1_vec[imix]!=ev_NJet_nevt_assoc_vec2[iimix])continue;}else{if(ev_NJet_nevt1_vec[imix]==ev_NJet_nevt_assoc_vec2[iimix])continue;}
	        if((ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt1_vec[imix].M() < mass + nsigmapeak*sigma) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix].M() < mass2 + nsigmapeak*sigma2)){
            h1->Fill(xmix3D);
            }
	  }}
      for(int imix=0; imix<nMixmult_nevt12; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec2[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec2[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
	  		if(ev_NJet_nevt1_vec2[imix] < 0 || ev_NJet_nevt_assoc_vec[iimix] < 0)continue;
			if(samejet){if(ev_NJet_nevt1_vec2[imix]!=ev_NJet_nevt_assoc_vec[iimix])continue;}else{if(ev_NJet_nevt1_vec2[imix]==ev_NJet_nevt_assoc_vec[iimix])continue;}
	        if((ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() > mass2 - nsigmapeak*sigma2 && ev_GoodTrackFourVectorTemp_nevt1_vec2[imix].M() < mass2 + nsigmapeak*sigma2) && (ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() > mass - nsigmapeak*sigma && ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].M() < mass + nsigmapeak*sigma)){
            h1->Fill(xmix3D);
            }
	  }}

}}}

void call_mix_gen_V0_jet_OS(std::vector<int> ev_ntrkoff_vec, std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec, std::vector<std::vector<Int_t>> NJet, std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec2, std::vector<std::vector<Int_t>> NJet2, int Evt_mix, THnSparse *h1, std::vector<double> pv_vtx_z, double vzcut, bool samejet, double mass, double sigma, double mass2, double sigma2, double nsigmapeak){

cout << "Random Mixing" << endl;

cout << "RECO Level" << endl;

cout << "Default: Nmix = " << Evt_mix << endl;

MixEventsgen_jet_V0_OS(0,4,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(5,9,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(10,14,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(15,19,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(20,24,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(25,29,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(30,34,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(35,39,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(40,44,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(45,49,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(50,54,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(55,59,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(60,64,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(65,69,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(70,74,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(75,79,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(80,84,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(85,89,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(90,94,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(95,99,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(100,104,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(105,109,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(110,114,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);
MixEventsgen_jet_V0_OS(115,119,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,NJet,ev_GoodTrackFourVector_vec2,NJet2,h1,samejet,mass,sigma,mass2,sigma2,nsigmapeak);

}


//for gen tests


void MixEventsgen(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, THnSparse* h1){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
   int takeAssociated = 0;
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
      if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      takeAssociated++;
      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
      for(int imix=0; imix<nMixmult_nevt1; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
            h1->Fill(xmix3D);
}}}}}

void call_mix_gen_SS(std::vector<int> ev_ntrkoff_vec, std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec, int Evt_mix, THnSparse *h1, std::vector<double> pv_vtx_z, double vzcut){

cout << "Random Mixing" << endl;

cout << "RECO Level" << endl;

cout << "Default: Nmix = " << Evt_mix << endl;

MixEventsgen(0,4,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(5,9,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(10,14,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(15,19,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(20,24,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(25,29,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(30,34,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(35,39,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(40,44,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(45,49,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(50,54,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(55,59,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(60,64,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(65,69,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(70,74,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(75,79,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(80,84,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(85,89,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(90,94,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(95,99,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(100,104,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(105,109,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(110,114,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);
MixEventsgen(115,119,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,h1);

}

void MixEventsgen_OS(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector2, THnSparse* h1){
int aux_n_evts = (int)ev_ntrkoff.size();
for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector[nevt1];
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec2= ev_GoodTrackFourVector2[nevt1];	   
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   int nMixmult_nevt12=ev_GoodTrackFourVectorTemp_nevt1_vec2.size();
   if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
   int takeAssociated = 0;
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
      if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector[nevt_assoc];
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec2= ev_GoodTrackFourVector2[nevt_assoc];
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      int nMixmult_nevt_assoc2=ev_GoodTrackFourVectorTemp_nevt_assoc_vec2.size();
      takeAssociated++;
      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
      for(int imix=0; imix<nMixmult_nevt1; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc2; iimix++){
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec2[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
            h1->Fill(xmix3D);
	  }}
      for(int imix=0; imix<nMixmult_nevt12; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec2[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec2[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
    	    Double_t ktmix=(psum2mix.Pt())/2.;
            Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
            h1->Fill(xmix3D);
	  }}

}}}

void call_mix_gen_OS(std::vector<int> ev_ntrkoff_vec, std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec, std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec2, int Evt_mix, THnSparse *h1, std::vector<double> pv_vtx_z, double vzcut){

cout << "Random Mixing" << endl;

cout << "RECO Level" << endl;

cout << "Default: Nmix = " << Evt_mix << endl;

MixEventsgen_OS(0,4,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(5,9,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(10,14,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(15,19,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(20,24,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(25,29,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(30,34,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(35,39,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(40,44,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(45,49,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(50,54,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(55,59,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(60,64,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(65,69,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(70,74,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(75,79,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(80,84,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(85,89,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(90,94,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(95,99,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(100,104,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(105,109,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(110,114,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);
MixEventsgen_OS(115,119,Evt_mix,ev_ntrkoff_vec,pv_vtx_z,vzcut,ev_GoodTrackFourVector_vec,ev_GoodTrackFourVector_vec2,h1);

}


void cplots(std::vector<TLorentzVector> GoodTrackFourVector, TH1D* hpt, TH1D* heta, TH1D* hphi, TH1D* hmass, double weight){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
        hpt->Fill(GoodTrackFourVector[itk1].Pt(),weight);
        heta->Fill(GoodTrackFourVector[itk1].Eta(),weight);
        hphi->Fill(GoodTrackFourVector[itk1].Phi(),weight);
        hmass->Fill(GoodTrackFourVector[itk1].M(),weight);
}}

void cplots_peak(std::vector<TLorentzVector> GoodTrackFourVector, TH1D* hpt, TH1D* heta, TH1D* hphi, TH1D* hmass, double mass, double sigma, double weight){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
    if(GoodTrackFourVector[itk1].M() > mass + 2.*sigma)continue;
    if(GoodTrackFourVector[itk1].M() < mass - 2.*sigma)continue;
        hpt->Fill(GoodTrackFourVector[itk1].Pt(),weight);
        heta->Fill(GoodTrackFourVector[itk1].Eta(),weight);
        hphi->Fill(GoodTrackFourVector[itk1].Phi(),weight);
        hmass->Fill(GoodTrackFourVector[itk1].M(),weight);
}}


//-----------------------------------------------------------------------------
//histograms
//-----------------------------------------------------------------------------

//events
  
  TH1D *nocutev=new TH1D("nocutev","",2,0,2);  

  TH1D *nocutPV=new TH1D("nocutPV","",2,0,2);  
  TH1D *nocutPU=new TH1D("nocutPU","",2,0,2);  
  TH1D *nocutHF=new TH1D("nocutHF","",2,0,2);  
  TH1D *nocutSC=new TH1D("nocutSC","",2,0,2);  
  TH1D *nocutALLX=new TH1D("nocutALLX","",2,0,2);  
  TH1D *nocutALL=new TH1D("nocutALL","",2,0,2);  


  TH1D *nev_K0s_ini=new TH1D("nev_K0s_ini","",2,0,2);  
  TH1D *nev_Lam_ini=new TH1D("nev_Lam_ini","",2,0,2);  
  TH1D *nev_ALam_ini=new TH1D("nev_ALam_ini","",2,0,2);  

  TH1D *nev_K0s_AS=new TH1D("nev_K0s_AS","",2,0,2);  
  TH1D *nev_Lam_AS=new TH1D("nev_Lam_AS","",2,0,2);  
  TH1D *nev_ALam_AS=new TH1D("nev_ALam_AS","",2,0,2);  

  
//K0s

  TH1D *nev_K0sT=new TH1D("nev_K0sT","",2,0,2);  
  TH1D *nev_K0sF=new TH1D("nev_K0sF","",2,0,2);  
  TH1D *nev_K0sTF=new TH1D("nev_K0sTF","",2,0,2);  
  TH1D *nev_K0s_ssT=new TH1D("nev_K0s_ssT","",2,0,2);  
  TH1D *nev_K0s_bbT=new TH1D("nev_K0s_bbT","",2,0,2);  
  TH1D *nev_K0s_sbT=new TH1D("nev_K0s_sbT","",2,0,2);  
  TH1D *nev_K0s_ssT_mix=new TH1D("nev_K0s_ssT_mix","",2,0,2);  
  TH1D *nev_K0s_bbT_mix=new TH1D("nev_K0s_bbT_mix","",2,0,2);  
  TH1D *nev_K0s_sbT_mix=new TH1D("nev_K0s_sbT_mix","",2,0,2); 
  TH1D *nev_K0s_ssT_etamix=new TH1D("nev_K0s_ssT_etamix","",2,0,2);  
  TH1D *nev_K0s_bbT_etamix=new TH1D("nev_K0s_bbT_etamix","",2,0,2);  
  TH1D *nev_K0s_sbT_etamix=new TH1D("nev_K0s_sbT_etamix","",2,0,2); 
  

//Lam

  TH1D *nev_LamT=new TH1D("nev_LamT","",2,0,2);  
  TH1D *nev_LamF=new TH1D("nev_LamF","",2,0,2);  
  TH1D *nev_LamTF=new TH1D("nev_LamTF","",2,0,2);  
  TH1D *nev_Lam_ssT=new TH1D("nev_Lam_ssT","",2,0,2);  
  TH1D *nev_Lam_bbT=new TH1D("nev_Lam_bbT","",2,0,2);  
  TH1D *nev_Lam_sbT=new TH1D("nev_Lam_sbT","",2,0,2);
  TH1D *nev_Lam_ssT_mix=new TH1D("nev_Lam_ssT_mix","",2,0,2);  
  TH1D *nev_Lam_bbT_mix=new TH1D("nev_Lam_bbT_mix","",2,0,2);  
  TH1D *nev_Lam_sbT_mix=new TH1D("nev_Lam_sbT_mix","",2,0,2);   
  TH1D *nev_Lam_ssT_etamix=new TH1D("nev_Lam_ssT_etamix","",2,0,2);  
  TH1D *nev_Lam_bbT_etamix=new TH1D("nev_Lam_bbT_etamix","",2,0,2);  
  TH1D *nev_Lam_sbT_etamix=new TH1D("nev_Lam_sbT_etamix","",2,0,2);   
  
  
//ALam

  TH1D *nev_ALamT=new TH1D("nev_ALamT","",2,0,2);  
  TH1D *nev_ALamF=new TH1D("nev_ALamF","",2,0,2);  
  TH1D *nev_ALamTF=new TH1D("nev_ALamTF","",2,0,2);  
  TH1D *nev_ALam_ssT=new TH1D("nev_ALam_ssT","",2,0,2);  
  TH1D *nev_ALam_bbT=new TH1D("nev_ALam_bbT","",2,0,2);  
  TH1D *nev_ALam_sbT=new TH1D("nev_ALam_sbT","",2,0,2);  
  TH1D *nev_ALam_ssT_mix=new TH1D("nev_ALam_ssT_mix","",2,0,2);  
  TH1D *nev_ALam_bbT_mix=new TH1D("nev_ALam_bbT_mix","",2,0,2);  
  TH1D *nev_ALam_sbT_mix=new TH1D("nev_ALam_sbT_mix","",2,0,2);  
  TH1D *nev_ALam_ssT_etamix=new TH1D("nev_ALam_ssT_etamix","",2,0,2);  
  TH1D *nev_ALam_bbT_etamix=new TH1D("nev_ALam_bbT_etamix","",2,0,2);  
  TH1D *nev_ALam_sbT_etamix=new TH1D("nev_ALam_sbT_etamix","",2,0,2);  
  

//LAL

  TH1D *nev_LALT=new TH1D("nev_LALT","",2,0,2);  
  TH1D *nev_LALF=new TH1D("nev_LALF","",2,0,2);  
  TH1D *nev_LALTF=new TH1D("nev_LALTF","",2,0,2);  
  TH1D *nev_LAL_ssT=new TH1D("nev_LAL_ssT","",2,0,2);  
  TH1D *nev_LAL_bbT=new TH1D("nev_LAL_bbT","",2,0,2);  
  TH1D *nev_LAL_sbT=new TH1D("nev_LAL_sbT","",2,0,2);  
  TH1D *nev_LAL_ssT_mix=new TH1D("nev_LAL_ssT_mix","",2,0,2);  
  TH1D *nev_LAL_bbT_mix=new TH1D("nev_LAL_bbT_mix","",2,0,2);  
  TH1D *nev_LAL_sbT_mix=new TH1D("nev_LAL_sbT_mix","",2,0,2);  
  TH1D *nev_LAL_ssT_etamix=new TH1D("nev_LAL_ssT_etamix","",2,0,2);  
  TH1D *nev_LAL_bbT_etamix=new TH1D("nev_LAL_bbT_etamix","",2,0,2);  
  TH1D *nev_LAL_sbT_etamix=new TH1D("nev_LAL_sbT_etamix","",2,0,2);  
  


//KL

  TH1D *nev_KLT=new TH1D("nev_KLT","",2,0,2); 
  TH1D *nev_KLF=new TH1D("nev_KLF","",2,0,2);  
  TH1D *nev_KLTF=new TH1D("nev_KLTF","",2,0,2);  
  TH1D *nev_KL_ssT=new TH1D("nev_KL_ssT","",2,0,2);  
  TH1D *nev_KL_bbT=new TH1D("nev_KL_bbT","",2,0,2);  
  TH1D *nev_KL_sbT=new TH1D("nev_KL_sbT","",2,0,2);  
  TH1D *nev_KL_ssT_mix=new TH1D("nev_KL_ssT_mix","",2,0,2);  
  TH1D *nev_KL_bbT_mix=new TH1D("nev_KL_bbT_mix","",2,0,2);  
  TH1D *nev_KL_sbT_mix=new TH1D("nev_KL_sbT_mix","",2,0,2);  
  TH1D *nev_KL_ssT_etamix=new TH1D("nev_KL_ssT_etamix","",2,0,2);  
  TH1D *nev_KL_bbT_etamix=new TH1D("nev_KL_bbT_etamix","",2,0,2);  
  TH1D *nev_KL_sbT_etamix=new TH1D("nev_KL_sbT_etamix","",2,0,2);  

//KAL

  TH1D *nev_KALT=new TH1D("nev_KALT","",2,0,2);  
  TH1D *nev_KALF=new TH1D("nev_KALF","",2,0,2);  
  TH1D *nev_KALTF=new TH1D("nev_KALTF","",2,0,2);  
  TH1D *nev_KAL_ssT=new TH1D("nev_KAL_ssT","",2,0,2);  
  TH1D *nev_KAL_bbT=new TH1D("nev_KAL_bbT","",2,0,2);  
  TH1D *nev_KAL_sbT=new TH1D("nev_KAL_sbT","",2,0,2);  
  TH1D *nev_KAL_ssT_mix=new TH1D("nev_KAL_ssT_mix","",2,0,2);  
  TH1D *nev_KAL_bbT_mix=new TH1D("nev_KAL_bbT_mix","",2,0,2);  
  TH1D *nev_KAL_sbT_mix=new TH1D("nev_KAL_sbT_mix","",2,0,2);  
  TH1D *nev_KAL_ssT_etamix=new TH1D("nev_KAL_ssT_etamix","",2,0,2);  
  TH1D *nev_KAL_bbT_etamix=new TH1D("nev_KAL_bbT_etamix","",2,0,2);  
  TH1D *nev_KAL_sbT_etamix=new TH1D("nev_KAL_sbT_etamix","",2,0,2);  
  

  TH1D *cone_K0s=new TH1D("cone_K0s","",1000,0.,7.);
  TH1D *cone_K0sD1=new TH1D("cone_K0sD1","",1000,0.,7.);
  TH1D *cone_K0sD2=new TH1D("cone_K0sD2","",1000,0.,7.);
  TH1D *cone_K0sT=new TH1D("cone_K0sT","",1000,0.,7.);
  TH1D *cone_K0sTD1=new TH1D("cone_K0sTD1","",1000,0.,7.);
  TH1D *cone_K0sTD2=new TH1D("cone_K0sTD2","",1000,0.,7.);

  TH1D *cone_Lam=new TH1D("cone_Lam","",1000,0.,7.);
  TH1D *cone_LamT=new TH1D("cone_LamT","",1000,0.,7.);
  TH1D *cone_LamD1=new TH1D("cone_LamD1","",1000,0.,7.);
  TH1D *cone_LamD2=new TH1D("cone_LamD2","",1000,0.,7.);
  TH1D *cone_LamTD1=new TH1D("cone_LamTD1","",1000,0.,7.);
  TH1D *cone_LamTD2=new TH1D("cone_LamTD2","",1000,0.,7.);

  TH1D *cone_Jet=new TH1D("cone_Jet","",1000,0.,7.);


  //vertex / ntrk / centrality  
  TH1D *h_vtx_z=new TH1D("h_vtx_z","",400,-20.,20.);  
  TH1D *h_vtx_rho=new TH1D("h_vtx_rho","",200,0.,2.); 
  TH1D *h_ntrk_cent=new TH1D("h_ntrk_cent","",500,0.,500.); 

  TH1D *h_ntrk_cent_1Ks=new TH1D("h_ntrk_cent_1Ks","",500,0.,500.); 
  TH1D *h_ntrk_cent_1Lam=new TH1D("h_ntrk_cent_1Lam","",500,0.,500.); 
  TH1D *h_ntrk_cent_1ALam=new TH1D("h_ntrk_cent_1ALam","",500,0.,500.); 

  
  //mass  

  TH1D *K0s_Mass=new TH1D("K0s_Mass","",400,0.4,0.6); 
  TH1D *Lam_Mass=new TH1D("Lam_Mass","",400,1.0,1.2);   
  TH1D *ALam_Mass=new TH1D("ALam_Mass","",400,1.0,1.2);  
  TH1D *LAL_Mass=new TH1D("LAL_Mass","",400,1.0,1.2);    
  
  
 TH1D *K0sctau=new TH1D("K0sctau","",500,0.,50.);
 TH1D *K0sdca3D=new TH1D("K0sdca3D","",200,0.,2.);
 TH1D *K0sDL=new TH1D("K0sDL","",10000,0.,1000.);
 TH1D *K0strkdca=new TH1D("K0strkdca","",200,0.,2.);
 TH1D *K0scostheta=new TH1D("K0scostheta","",2000,0.99,1.01);
 TH1D *K0s_Vtx=new TH1D("K0s_vtx","", 100,0.,10.); 


 TH1D *Lamctau=new TH1D("Lamctau","",500,0.,50.);
 TH1D *Lamdca3D=new TH1D("Lamdca3D","",200,0.,2.);
 TH1D *LamDL=new TH1D("LamDL","",10000,0.,1000.);
 TH1D *Lamtrkdca=new TH1D("Lamtrkdca","",200,0.,2.);
 TH1D *Lamcostheta=new TH1D("Lamcostheta","",2000,0.99,1.01);
 TH1D *Lam_Vtx=new TH1D("Lam_vtx","",100,0.,10.); 


 TH1D *ALamctau=new TH1D("ALamctau","",500,0.,50.);
 TH1D *ALamdca3D=new TH1D("ALamdca3D","",200,0.,2.);
 TH1D *ALamDL=new TH1D("ALamDL","",10000,0.,1000.);
 TH1D *ALamtrkdca=new TH1D("ALamtrkdca","",200,0.,2.);
 TH1D *ALamcostheta=new TH1D("ALamcostheta","",2000,0.99,1.01);
 TH1D *ALam_Vtx=new TH1D("ALam_vtx","",100,0.,10.); 

 TH1D *K0s_d1_Nhits=new TH1D("K0s_d1_Nhits","",100,0.,100.);
 TH1D *K0s_d2_Nhits=new TH1D("K0s_d2_Nhits","",100,0.,100.);
 TH1D *Lam_d1_Nhits=new TH1D("Lam_d1_Nhits","",100,0.,100.);
 TH1D *Lam_d2_Nhits=new TH1D("Lam_d2_Nhits","",100,0.,100.);
 TH1D *ALam_d1_Nhits=new TH1D("ALam_d1_Nhits","",100,0.,100.);
 TH1D *ALam_d2_Nhits=new TH1D("ALam_d2_Nhits","",100,0.,100.);
 
 TH1D *K0s_d1_dxy=new TH1D("K0s_d1_dxy","",10000,0.,100.);
 TH1D *K0s_d2_dxy=new TH1D("K0s_d2_dxy","",10000,0.,100.);
 TH1D *K0s_d1_dz=new TH1D("K0s_d1_dz","",10000,0.,100.);
 TH1D *K0s_d2_dz=new TH1D("K0s_d2_dz","",10000,0.,100.);
 TH1D *Lam_d1_dxy=new TH1D("Lam_d1_dxy","",10000,0.,100.);
 TH1D *Lam_d2_dxy=new TH1D("Lam_d2_dxy","",10000,0.,100.);
 TH1D *Lam_d1_dz=new TH1D("Lam_d1_dz","",10000,0.,100.);
 TH1D *Lam_d2_dz=new TH1D("Lam_d2_dz","",10000,0.,100.);
 TH1D *ALam_d1_dxy=new TH1D("ALam_d1_dxy","",10000,0.,100.);
 TH1D *ALam_d2_dxy=new TH1D("ALam_d2_dxy","",10000,0.,100.);
 TH1D *ALam_d1_dz=new TH1D("ALam_d1_dz","",10000,0.,100.);
 TH1D *ALam_d2_dz=new TH1D("ALam_d2_dz","",10000,0.,100.);

  TH1D *DR_K0s=new TH1D("DR_K0s","",1000,0.,7.);
  TH1D *DR_Lam=new TH1D("DR_Lam","",1000,0.,7.);
  TH1D *DR_ALam=new TH1D("DR_ALam","",1000,0.,7.);

  TH1D *DPT_K0s=new TH1D("DPT_K0s","",1000,0.,1.);
  TH1D *DPT_Lam=new TH1D("DPT_Lam","",1000,0.,1.);
  TH1D *DPT_ALam=new TH1D("DPT_ALam","",1000,0.,1.);


double pla[21] = {0.0,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,5.0,5.6,6.6,8.5}; 

double mla[401] = {1.    ,  1.0005,  1.001 ,  1.0015,  1.002 ,  1.0025,  1.003 , 1.0035,  1.004 ,  1.0045,  1.005 ,  1.0055,  1.006 ,  1.0065,  1.007 ,  1.0075,  1.008 ,  1.0085,  1.009 ,  1.0095,  1.01  ,  1.0105,  1.011 ,  1.0115,  1.012 ,  1.0125,  1.013 ,  1.0135, 1.014 ,  1.0145,  1.015 ,  1.0155,  1.016 ,  1.0165,  1.017 , 1.0175,  1.018 ,  1.0185,  1.019 ,  1.0195,  1.02  ,  1.0205, 1.021 ,  1.0215,  1.022 ,  1.0225,  1.023 ,  1.0235,  1.024 , 1.0245,  1.025 ,  1.0255,  1.026 ,  1.0265,  1.027 ,  1.0275, 1.028 ,  1.0285,  1.029 ,  1.0295,  1.03  ,  1.0305,  1.031 ,  1.0315,  1.032 ,  1.0325,  1.033 ,  1.0335,  1.034 ,  1.0345,  1.035 ,  1.0355,  1.036 ,  1.0365,  1.037 ,  1.0375,  1.038 , 1.0385,  1.039 ,  1.0395,  1.04  ,  1.0405,  1.041 ,  1.0415, 1.042 ,  1.0425,  1.043 ,  1.0435,  1.044 ,  1.0445,  1.045 , 1.0455,  1.046 ,  1.0465,  1.047 ,  1.0475,  1.048 ,  1.0485, 1.049 ,  1.0495,  1.05  ,  1.0505,  1.051 ,  1.0515,  1.052 , 1.0525,  1.053 ,  1.0535,  1.054 ,  1.0545,  1.055 ,  1.0555, 1.056 ,  1.0565,  1.057 ,  1.0575,  1.058 ,  1.0585,  1.059 , 1.0595,  1.06  ,  1.0605,  1.061 ,  1.0615,  1.062 ,  1.0625, 1.063 ,  1.0635,  1.064 ,  1.0645,  1.065 ,  1.0655,  1.066 , 1.0665,  1.067 ,  1.0675,  1.068 ,  1.0685,  1.069 ,  1.0695,  1.07  ,  1.0705,  1.071 ,  1.0715,  1.072 ,  1.0725,  1.073 , 1.0735,  1.074 ,  1.0745,  1.075 ,  1.0755,  1.076 ,  1.0765,  1.077 ,  1.0775,  1.078 ,  1.0785,  1.079 ,  1.0795,  1.08  , 1.0805,  1.081 ,  1.0815,  1.082 ,  1.0825,  1.083 ,  1.0835, 1.084 ,  1.0845,  1.085 ,  1.0855,  1.086 ,  1.0865,  1.087 , 1.0875,  1.088 ,  1.0885,  1.089 ,  1.0895,  1.09  ,  1.0905,  1.091 ,  1.0915,  1.092 ,  1.0925,  1.093 ,  1.0935,  1.094 , 1.0945,  1.095 ,  1.0955,  1.096 ,  1.0965,  1.097 ,  1.0975, 1.098 ,  1.0985,  1.099 ,  1.0995,  1.1   ,  1.1005,  1.101 ,  1.1015,  1.102 ,  1.1025,  1.103 ,  1.1035,  1.104 ,  1.1045, 1.105 ,  1.1055,  1.106 ,  1.1065,  1.107 ,  1.1075,  1.108 , 1.1085,  1.109 ,  1.1095,  1.11  ,  1.1105,  1.111 ,  1.1115, 1.112 ,  1.1125,  1.113 ,  1.1135,  1.114 ,  1.1145,  1.115 ,  1.1155,  1.116 ,  1.1165,  1.117 ,  1.1175,  1.118 ,  1.1185, 1.119 ,  1.1195,  1.12  ,  1.1205,  1.121 ,  1.1215,  1.122 , 1.1225,  1.123 ,  1.1235,  1.124 ,  1.1245,  1.125 ,  1.1255, 1.126 ,  1.1265,  1.127 ,  1.1275,  1.128 ,  1.1285,  1.129 , 1.1295,  1.13  ,  1.1305,  1.131 ,  1.1315,  1.132 ,  1.1325, 1.133 ,  1.1335,  1.134 ,  1.1345,  1.135 ,  1.1355,  1.136 ,  1.1365,  1.137 ,  1.1375,  1.138 ,  1.1385,  1.139 ,  1.1395, 1.14  ,  1.1405,  1.141 ,  1.1415,  1.142 ,  1.1425,  1.143 ,  1.1435,  1.144 ,  1.1445,  1.145 ,  1.1455,  1.146 ,  1.1465,  1.147 ,  1.1475,  1.148 ,  1.1485,  1.149 ,  1.1495,  1.15  , 1.1505,  1.151 ,  1.1515,  1.152 ,  1.1525,  1.153 ,  1.1535, 1.154 ,  1.1545,  1.155 ,  1.1555,  1.156 ,  1.1565,  1.157 , 1.1575,  1.158 ,  1.1585,  1.159 ,  1.1595,  1.16  ,  1.1605, 1.161 ,  1.1615,  1.162 ,  1.1625,  1.163 ,  1.1635,  1.164 , 1.1645,  1.165 ,  1.1655,  1.166 ,  1.1665,  1.167 ,  1.1675, 1.168 ,  1.1685,  1.169 ,  1.1695,  1.17  ,  1.1705,  1.171 ,  1.1715,  1.172 ,  1.1725,  1.173 ,  1.1735,  1.174 ,  1.1745,  1.175 ,  1.1755,  1.176 ,  1.1765,  1.177 ,  1.1775,  1.178 , 1.1785,  1.179 ,  1.1795,  1.18  ,  1.1805,  1.181 ,  1.1815,  1.182 ,  1.1825,  1.183 ,  1.1835,  1.184 ,  1.1845,  1.185 , 1.1855,  1.186 ,  1.1865,  1.187 ,  1.1875,  1.188 ,  1.1885,  1.189 ,  1.1895,  1.19  ,  1.1905,  1.191 ,  1.1915,  1.192 ,  1.1925,  1.193 ,  1.1935,  1.194 ,  1.1945,  1.195 ,  1.1955,  1.196 ,  1.1965,  1.197 ,  1.1975,  1.198 ,  1.1985,  1.199 , 1.1995,  1.2 };

double mxi[201] = {1.2  , 1.201, 1.202, 1.203, 1.204, 1.205, 1.206, 1.207, 1.208, 1.209, 1.21 , 1.211, 1.212, 1.213, 1.214, 1.215, 1.216, 1.217,
1.218, 1.219, 1.22 , 1.221, 1.222, 1.223, 1.224, 1.225, 1.226, 1.227, 1.228, 1.229, 1.23 , 1.231, 1.232, 1.233, 1.234, 1.235, 1.236, 1.237, 1.238, 1.239, 1.24 , 1.241, 1.242, 1.243, 1.244,1.245, 1.246, 1.247, 1.248, 1.249, 1.25 , 1.251, 1.252, 1.253, 1.254, 1.255, 1.256, 1.257, 1.258, 1.259, 1.26 , 1.261, 1.262, 1.263, 1.264, 1.265, 1.266, 1.267, 1.268, 1.269, 1.27 , 1.271, 1.272, 1.273, 1.274, 1.275, 1.276, 1.277, 1.278, 1.279, 1.28 ,
1.281, 1.282, 1.283, 1.284, 1.285, 1.286, 1.287, 1.288, 1.289, 1.29 , 1.291, 1.292, 1.293, 1.294, 1.295, 1.296, 1.297, 1.298, 1.299, 1.3  , 1.301, 1.302, 1.303, 1.304, 1.305, 1.306, 1.307, 1.308, 1.309, 1.31 , 1.311, 1.312, 1.313, 1.314, 1.315, 1.316, 1.317, 1.318, 1.319, 1.32 , 1.321, 1.322, 1.323, 1.324, 1.325, 1.326, 1.327, 1.328, 1.329, 1.33 , 1.331, 1.332, 1.333, 1.334, 1.335, 1.336, 1.337, 1.338, 1.339, 1.34 , 1.341, 1.342, 1.343, 1.344, 1.345, 1.346, 1.347, 1.348, 1.349, 1.35 , 1.351, 1.352, 1.353, 1.354, 1.355, 1.356, 1.357, 1.358, 1.359, 1.36 , 1.361, 1.362, 1.363, 1.364, 1.365, 1.366, 1.367, 1.368, 1.369, 1.37 , 1.371, 1.372, 1.373, 1.374, 1.375, 1.376, 1.377, 1.378, 1.379, 1.38 , 1.381, 1.382, 1.383, 1.384, 1.385, 1.386, 1.387, 1.388, 1.389, 1.39 , 1.391, 1.392, 1.393, 1.394, 1.395, 1.396, 1.397, 1.398, 1.399,1.4};

double mOm[201] = {1.6, 1.601, 1.602, 1.603, 1.604, 1.605, 1.606, 1.607, 1.608, 1.609, 1.61 , 1.611, 1.612, 1.613, 1.614, 1.615, 1.616, 1.617, 1.618, 1.619, 1.62 , 1.621, 1.622, 1.623, 1.624, 1.625, 1.626, 1.627, 1.628, 1.629, 1.63 , 1.631, 1.632, 1.633, 1.634, 1.635, 1.636, 1.637, 1.638, 1.639, 1.64 , 1.641, 1.642, 1.643, 1.644, 1.645, 1.646, 1.647, 1.648, 1.649, 1.65 , 1.651, 1.652, 1.653, 1.654, 1.655, 1.656, 1.657, 1.658, 1.659, 1.66 , 1.661, 1.662, 1.663, 1.664, 1.665, 1.666, 1.667, 1.668, 1.669, 1.67 , 1.671, 1.672, 1.673, 1.674, 1.675, 1.676, 1.677, 1.678, 1.679, 1.68, 1.681, 1.682, 1.683, 1.684, 1.685, 1.686, 1.687, 1.688, 1.689, 1.69 , 1.691, 1.692, 1.693, 1.694, 1.695, 1.696, 1.697, 1.698, 1.699, 1.7  , 1.701, 1.702, 1.703, 1.704, 1.705, 1.706, 1.707, 1.708, 1.709, 1.71 , 1.711, 1.712, 1.713, 1.714, 1.715, 1.716, 1.717, 1.718, 1.719, 1.72 , 1.721, 1.722, 1.723, 1.724, 1.725, 1.726, 1.727, 1.728, 1.729, 1.73 , 1.731, 1.732, 1.733, 1.734, 1.735, 1.736, 1.737, 1.738, 1.739, 1.74 , 1.741, 1.742, 1.743,
1.744, 1.745, 1.746, 1.747, 1.748, 1.749, 1.75 , 1.751, 1.752, 1.753, 1.754, 1.755, 1.756, 1.757, 1.758, 1.759, 1.76 , 1.761, 1.762, 1.763, 1.764, 1.765, 1.766, 1.767, 1.768, 1.769, 1.77, 1.771, 1.772, 1.773, 1.774, 1.775, 1.776, 1.777, 1.778, 1.779, 1.78 , 1.781, 1.782, 1.783, 1.784, 1.785, 1.786, 1.787, 1.788, 1.789, 1.79 , 1.791, 1.792, 1.793, 1.794, 1.795, 1.796, 1.797, 1.798, 1.799, 1.8};

//////MC

  TH1D* pT_Lamb_from_Xi_T=new TH1D("pT_Lamb_from_Xi_T","",20,pla); //
  TH1D* pT_Lamb_from_Om_T=new TH1D("pT_Lamb_from_Om_T","",20,pla); //
  TH1D* pT_Lamb_from_prompt_T=new TH1D("pT_Lamb_from_prompt_T","",20,pla); //done

  TH1D* pT_Lamb_from_Xi_TF=new TH1D("pT_Lamb_from_Xi_TF","",20,pla); //done
  TH1D* pT_Lamb_from_Om_TF=new TH1D("pT_Lamb_from_Om_TF","",20,pla);//done
  TH1D* pT_Lamb_from_prompt_TF=new TH1D("pT_Lamb_from_prompt_TF","",20,pla);//done

  TH1D* pT_Lamb_daughter_Xi=new TH1D("pT_Lamb_daughter_Xi","",20,pla);//done
  TH1D* pT_Lamb_daughter_Om=new TH1D("pT_Lamb_daughter_Om","",20,pla);//done
  TH1D* pT_Lamb_daughter_prompt=new TH1D("pT_Lamb_daughter_prompt","",20,pla);//done
  
  TH2D* pT_Xi=new TH2D("pT_Xi","",20,pla,200,mxi); //
  TH2D* pT_Om=new TH2D("pT_Om","",20,pla,200,mOm); //
  TH2D* pT_Lamb=new TH2D("pT_Lamb","",20,pla,400,mla); //
  
//  THnSparseD *h_Xi_Xi  = new THnSparseD("h_Xi_Xi","h_Xi_Xi",3,bins3D,xmin3D,xmax3D);
//  THnSparseD *h_Om_Om  = new THnSparseD("h_Om_Om","h_Om_Om",3,bins3D,xmin3D,xmax3D);


  //-----------------------------------------
  
  //HBT histograms
  
  //------------------------------------------------------------------------

//for mass  
  Int_t bins3DM_K0s[4]            =   {400 , 15  ,  100    , 100};
  Double_t xmin3DM_K0s[4]    =   {0.4  , 0.0 ,  0.0    , 0.0};
  Double_t xmax3DM_K0s[4]   =   {0.6  , 3.0 ,  500.0 , 4.0};

  Int_t bins3DM_LAL[4]            =   {400 , 15  ,  100     , 100};
  Double_t xmin3DM_LAL[4]    =   {1.0  , 0.0 ,  0.0    , 0.0};
  Double_t xmax3DM_LAL[4]   =   {1.2  , 3.0 ,  500.0, 4.0};

  
  //DONE
  Int_t bins3DM_KL[4]            =   {4  , 6   , 5     ,  5};
  Double_t xmin3DM_KL[4]    =   {0.4 , 0.0, 0.0    ,  0.0};
  Double_t xmax3DM_KL[4]   =   {1.2 , 3.0, 500.0,  5.0};
  
// for qinv -> DONE

  Int_t bins3D[3]            =   {500 , 30 ,   100};
  Double_t xmin3D[3]    =   {0.0  , 0.0 ,   0.0};
  Double_t xmax3D[3]   =   {10.0, 3.0,   500.};

//for bias  test -> DONE
  
  Int_t bins3Db[3]            =  {200, 6   , 25 };
  Double_t xmin3Db[3]    =  {0.0 , 0.0, 0.0 };
  Double_t xmax3Db[3]   =  {1.0 , 3.0, 500.};
  
  
THnSparseD *hGen_K0s_K0s_Matched  = new THnSparseD("hGen_K0s_K0s_Matched","hGen_K0s_K0s_Matched",3,bins3D,xmin3D,xmax3D);
THnSparseD *hGen_K0s_K0s_Matched_mix  = new THnSparseD("hGen_K0s_K0s_Matched_mix","hGen_K0s_K0s_Matched_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hGen_Lam_Lam_Matched  = new THnSparseD("hGen_Lam_Lam_Matched","hGen_Lam_Lam_Matched",3,bins3D,xmin3D,xmax3D);
THnSparseD *hGen_Lam_Lam_Matched_mix  = new THnSparseD("hGen_Lam_Lam_Matched_mix","hGen_Lam_Lam_Matched_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hGen_ALam_ALam_Matched  = new THnSparseD("hGen_ALam_ALam_Matched","hGen_ALam_ALam_Matched",3,bins3D,xmin3D,xmax3D);
THnSparseD *hGen_ALam_ALam_Matched_mix  = new THnSparseD("hGen_ALam_ALam_Matched_mix","hGen_ALam_ALam_Matched_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hGen_Lam_ALam_Matched  = new THnSparseD("hGen_Lam_ALam_Matched","hGen_Lam_ALam_Matched",3,bins3D,xmin3D,xmax3D);
THnSparseD *hGen_Lam_ALam_Matched_mix  = new THnSparseD("hGen_Lam_ALam_Matched_mix","hGen_Lam_ALam_Matched_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hGen_K0s_Lam_Matched  = new THnSparseD("hGen_K0s_Lam_Matched","hGen_K0s_Lam_Matched",3,bins3D,xmin3D,xmax3D);
THnSparseD *hGen_K0s_Lam_Matched_mix  = new THnSparseD("hGen_K0s_Lam_Matched_mix","hGen_K0s_Lam_Matched_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hGen_K0s_ALam_Matched  = new THnSparseD("hGen_K0s_ALam_Matched","hGen_K0s_ALam_Matched",3,bins3D,xmin3D,xmax3D);
THnSparseD *hGen_K0s_ALam_Matched_mix  = new THnSparseD("hGen_K0s_ALam_Matched_mix","hGen_K0s_ALam_Matched_mix",3,bins3D,xmin3D,xmax3D);


THnSparseD *hGen_K0s_K0s_MatchedT  = new THnSparseD("hGen_K0s_K0s_MatchedT","hGen_K0s_K0s_MatchedT",3,bins3D,xmin3D,xmax3D);
THnSparseD *hGen_K0s_K0s_MatchedT_mix  = new THnSparseD("hGen_K0s_K0s_MatchedT_mix","hGen_K0s_K0s_MatchedT_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hGen_Lam_Lam_MatchedT  = new THnSparseD("hGen_Lam_Lam_MatchedT","hGen_Lam_Lam_MatchedT",3,bins3D,xmin3D,xmax3D);
THnSparseD *hGen_Lam_Lam_MatchedT_mix  = new THnSparseD("hGen_Lam_Lam_MatchedT_mix","hGen_Lam_Lam_MatchedT_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hGen_ALam_ALam_MatchedT  = new THnSparseD("hGen_ALam_ALam_MatchedT","hGen_ALam_ALam_MatchedT",3,bins3D,xmin3D,xmax3D);
THnSparseD *hGen_ALam_ALam_MatchedT_mix  = new THnSparseD("hGen_ALam_ALam_MatchedT_mix","hGen_ALam_ALam_MatchedT_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hGen_Lam_ALam_MatchedT  = new THnSparseD("hGen_Lam_ALam_MatchedT","hGen_Lam_ALam_MatchedT",3,bins3D,xmin3D,xmax3D);
THnSparseD *hGen_Lam_ALam_MatchedT_mix  = new THnSparseD("hGen_Lam_ALam_MatchedT_mix","hGen_Lam_ALam_MatchedT_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hGen_K0s_Lam_MatchedT  = new THnSparseD("hGen_K0s_Lam_MatchedT","hGen_K0s_Lam_MatchedT",3,bins3D,xmin3D,xmax3D);
THnSparseD *hGen_K0s_Lam_MatchedT_mix  = new THnSparseD("hGen_K0s_Lam_MatchedT_mix","hGen_K0s_Lam_MatchedT_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hGen_K0s_ALam_MatchedT  = new THnSparseD("hGen_K0s_ALam_MatchedT","hGen_K0s_ALam_MatchedT",3,bins3D,xmin3D,xmax3D);
THnSparseD *hGen_K0s_ALam_MatchedT_mix  = new THnSparseD("hGen_K0s_ALam_MatchedT_mix","hGen_K0s_ALam_MatchedT_mix",3,bins3D,xmin3D,xmax3D);
  
  
//K0s
//True

THnSparseD *hMass_K0s_K0s_T  = new THnSparseD("hMass_K0s_K0s_T","hMass_K0s_K0s_T",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_K0s_1_T  = new THnSparseD("hMass_K0s_1_T","hMass_K0s_1_T",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_K0s_2_T  = new THnSparseD("hMass_K0s_2_T","hMass_K0s_2_T",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   

  
THnSparseD *hpeak_K0s_K0s_T  = new THnSparseD("hpeak_K0s_K0s_T","hpeak_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakcos_K0s_K0s_T  = new THnSparseD("hpeakcos_K0s_K0s_T","hpeakcos_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_rot_K0s_K0s_T  = new THnSparseD("hpeak_rot_K0s_K0s_T","hpeak_rot_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_inv_K0s_K0s_T  = new THnSparseD("hpeak_inv_K0s_K0s_T","hpeak_inv_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_K0s_K0s_T  = new THnSparseD("hside_K0s_K0s_T","hside_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideL_K0s_K0s_T  = new THnSparseD("hsideL_K0s_K0s_T","hsideL_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideR_K0s_K0s_T  = new THnSparseD("hsideR_K0s_K0s_T","hsideR_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_rot_K0s_K0s_T  = new THnSparseD("hside_rot_K0s_K0s_T","hside_rot_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_inv_K0s_K0s_T  = new THnSparseD("hside_inv_K0s_K0s_T","hside_inv_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_K0s_K0s_T  = new THnSparseD("hpeakside_K0s_K0s_T","hpeakside_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideL_K0s_K0s_T  = new THnSparseD("hpeaksideL_K0s_K0s_T","hpeaksideL_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideR_K0s_K0s_T  = new THnSparseD("hpeaksideR_K0s_K0s_T","hpeaksideR_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_rot_K0s_K0s_T  = new THnSparseD("hpeakside_rot_K0s_K0s_T","hpeakside_rot_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_inv_K0s_K0s_T  = new THnSparseD("hpeakside_inv_K0s_K0s_T","hpeakside_inv_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_T  = new THnSparseD("hpeakd1_K0s_K0s_T","hpeakd1_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_T  = new THnSparseD("hpeakd2_K0s_K0s_T","hpeakd2_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_T  = new THnSparseD("hpeakd12_K0s_K0s_T","hpeakd12_K0s_K0s_T",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_K0s_K0s_T_mix  = new THnSparseD("hpeak_K0s_K0s_T_mix","hpeak_K0s_K0s_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_K0s_K0s_T_mix  = new THnSparseD("hside_K0s_K0s_T_mix","hside_K0s_K0s_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideL_K0s_K0s_T_mix  = new THnSparseD("hsideL_K0s_K0s_T_mix","hsideL_K0s_K0s_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideR_K0s_K0s_T_mix  = new THnSparseD("hsideR_K0s_K0s_T_mix","hsideR_K0s_K0s_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_K0s_K0s_T_mix  = new THnSparseD("hpeakside_K0s_K0s_T_mix","hpeakside_K0s_K0s_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideL_K0s_K0s_T_mix  = new THnSparseD("hpeaksideL_K0s_K0s_T_mix","hpeaksideL_K0s_K0s_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideR_K0s_K0s_T_mix  = new THnSparseD("hpeaksideR_K0s_K0s_T_mix","hpeaksideR_K0s_K0s_T_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_K0s_K0s_T_etamix  = new THnSparseD("hpeak_K0s_K0s_T_etamix","hpeak_K0s_K0s_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_K0s_K0s_T_etamix  = new THnSparseD("hside_K0s_K0s_T_etamix","hside_K0s_K0s_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideL_K0s_K0s_T_etamix  = new THnSparseD("hsideL_K0s_K0s_T_etamix","hsideL_K0s_K0s_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideR_K0s_K0s_T_etamix  = new THnSparseD("hsideR_K0s_K0s_T_etamix","hsideR_K0s_K0s_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_K0s_K0s_T_etamix  = new THnSparseD("hpeakside_K0s_K0s_T_etamix","hpeakside_K0s_K0s_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideL_K0s_K0s_T_etamix  = new THnSparseD("hpeaksideL_K0s_K0s_T_etamix","hpeaksideL_K0s_K0s_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideR_K0s_K0s_T_etamix  = new THnSparseD("hpeaksideR_K0s_K0s_T_etamix","hpeaksideR_K0s_K0s_T_etamix",3,bins3D,xmin3D,xmax3D);

TH1D *V0chi2d1_diff_T_K0s=new TH1D("V0chi2d1_diff_T_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_T_K0s=new TH1D("V0chi2d2_diff_T_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_T_K0s=new TH1D("V0chi2d12_diff_T_K0s","",10000,-0.01,0.11);

//False

THnSparseD *hMass_K0s_K0s_F  = new THnSparseD("hMass_K0s_K0s_F","hMass_K0s_K0s_F",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hpeak_K0s_K0s_F  = new THnSparseD("hpeak_K0s_K0s_F","hpeak_K0s_K0s_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_K0s_K0s_F  = new THnSparseD("hside_K0s_K0s_F","hside_K0s_K0s_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_K0s_K0s_F  = new THnSparseD("hpeakside_K0s_K0s_F","hpeakside_K0s_K0s_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_F  = new THnSparseD("hpeakd1_K0s_K0s_F","hpeakd1_K0s_K0s_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_F  = new THnSparseD("hpeakd2_K0s_K0s_F","hpeakd2_K0s_K0s_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_F  = new THnSparseD("hpeakd12_K0s_K0s_F","hpeakd12_K0s_K0s_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_F_bias  = new THnSparseD("hpeak_K0s_K0s_F_bias","hpeak_K0s_K0s_F_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_F_biasd1  = new THnSparseD("hpeak_K0s_K0s_F_biasd1","hpeak_K0s_K0s_F_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_F_biasd2  = new THnSparseD("hpeak_K0s_K0s_F_biasd2","hpeak_K0s_K0s_F_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_F_biasd12  = new THnSparseD("hpeak_K0s_K0s_F_biasd12","hpeak_K0s_K0s_F_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_F_K0s=new TH1D("V0chi2d1_diff_F_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_F_K0s=new TH1D("V0chi2d2_diff_F_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_F_K0s=new TH1D("V0chi2d12_diff_F_K0s","",10000,-0.01,0.11);

  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
THnSparseD *hpeak_K0s_K0s_TF  = new THnSparseD("hpeak_K0s_K0s_TF","hpeak_K0s_K0s_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_TF  = new THnSparseD("hpeakd1_K0s_K0s_TF","hpeakd1_K0s_K0s_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_TF  = new THnSparseD("hpeakd2_K0s_K0s_TF","hpeakd2_K0s_K0s_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_TF  = new THnSparseD("hpeakd12_K0s_K0s_TF","hpeakd12_K0s_K0s_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_K0s_K0s_TF  = new THnSparseD("hside_K0s_K0s_TF","hside_K0s_K0s_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_K0s_K0s_TF  = new THnSparseD("hpeakside_K0s_K0s_TF","hpeakside_K0s_K0s_TF",3,bins3D,xmin3D,xmax3D);  
THnSparseD *hpeak_K0s_K0s_TF_bias  = new THnSparseD("hpeak_K0s_K0s_TF_bias","hpeak_K0s_K0s_TF_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TF_biasd1  = new THnSparseD("hpeak_K0s_K0s_TF_biasd1","hpeak_K0s_K0s_TF_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TF_biasd2  = new THnSparseD("hpeak_K0s_K0s_TF_biasd2","hpeak_K0s_K0s_TF_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TF_biasd12  = new THnSparseD("hpeak_K0s_K0s_TF_biasd12","hpeak_K0s_K0s_TF_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TF_K0s=new TH1D("V0chi2d1_diff_TF_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TF_K0s=new TH1D("V0chi2d2_diff_TF_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TF_K0s=new TH1D("V0chi2d12_diff_TF_K0s","",10000,-0.01,0.11);
  
/////// matching hist ////// 
  
THnSparseD *hMass_K0s_K0s_TM  = new THnSparseD("hMass_K0s_K0s_TM","hMass_K0s_K0s_TM",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_K0s_K0s_TU  = new THnSparseD("hMass_K0s_K0s_TU","hMass_K0s_K0s_TU",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_K0s_K0s_FM  = new THnSparseD("hMass_K0s_K0s_FM","hMass_K0s_K0s_FM",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_K0s_K0s_FU  = new THnSparseD("hMass_K0s_K0s_FU","hMass_K0s_K0s_FU",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
 
THnSparseD *hpeak_K0s_K0s_TM  = new THnSparseD("hpeak_K0s_K0s_TM","hpeak_K0s_K0s_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_TM  = new THnSparseD("hpeakd1_K0s_K0s_TM","hpeakd1_K0s_K0s_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_TM  = new THnSparseD("hpeakd2_K0s_K0s_TM","hpeakd2_K0s_K0s_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_TM  = new THnSparseD("hpeakd12_K0s_K0s_TM","hpeakd12_K0s_K0s_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_TM_bias  = new THnSparseD("hpeak_K0s_K0s_TM_bias","hpeak_K0s_K0s_TM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TM_biasd1  = new THnSparseD("hpeak_K0s_K0s_TM_biasd1","hpeak_K0s_K0s_TM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TM_biasd2  = new THnSparseD("hpeak_K0s_K0s_TM_biasd2","hpeak_K0s_K0s_TM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TM_biasd12  = new THnSparseD("hpeak_K0s_K0s_TM_biasd12","hpeak_K0s_K0s_TM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_K0s=new TH1D("V0chi2d1_diff_TM_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_K0s=new TH1D("V0chi2d2_diff_TM_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_K0s=new TH1D("V0chi2d12_diff_TM_K0s","",10000,-0.01,0.11);

THnSparseD *hpeak_K0s_K0s_TU  = new THnSparseD("hpeak_K0s_K0s_TU","hpeak_K0s_K0s_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_TU  = new THnSparseD("hpeakd1_K0s_K0s_TU","hpeakd1_K0s_K0s_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_TU  = new THnSparseD("hpeakd2_K0s_K0s_TU","hpeakd2_K0s_K0s_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_TU  = new THnSparseD("hpeakd12_K0s_K0s_TU","hpeakd12_K0s_K0s_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_TU_bias  = new THnSparseD("hpeak_K0s_K0s_TU_bias","hpeak_K0s_K0s_TU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TU_biasd1  = new THnSparseD("hpeak_K0s_K0s_TU_biasd1","hpeak_K0s_K0s_TU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TU_biasd2  = new THnSparseD("hpeak_K0s_K0s_TU_biasd2","hpeak_K0s_K0s_TU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TU_biasd12  = new THnSparseD("hpeak_K0s_K0s_TU_biasd12","hpeak_K0s_K0s_TU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_K0s=new TH1D("V0chi2d1_diff_TU_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_K0s=new TH1D("V0chi2d2_diff_TU_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_K0s=new TH1D("V0chi2d12_diff_TU_K0s","",10000,-0.01,0.11);

THnSparseD *hpeak_K0s_K0s_FM  = new THnSparseD("hpeak_K0s_K0s_FM","hpeak_K0s_K0s_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_FM  = new THnSparseD("hpeakd1_K0s_K0s_FM","hpeakd1_K0s_K0s_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_FM  = new THnSparseD("hpeakd2_K0s_K0s_FM","hpeakd2_K0s_K0s_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_FM  = new THnSparseD("hpeakd12_K0s_K0s_FM","hpeakd12_K0s_K0s_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_FM_bias  = new THnSparseD("hpeak_K0s_K0s_FM_bias","hpeak_K0s_K0s_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_FM_biasd1  = new THnSparseD("hpeak_K0s_K0s_FM_biasd1","hpeak_K0s_K0s_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_FM_biasd2  = new THnSparseD("hpeak_K0s_K0s_FM_biasd2","hpeak_K0s_K0s_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_FM_biasd12  = new THnSparseD("hpeak_K0s_K0s_FM_biasd12","hpeak_K0s_K0s_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FM_K0s=new TH1D("V0chi2d1_diff_FM_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FM_K0s=new TH1D("V0chi2d2_diff_FM_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FM_K0s=new TH1D("V0chi2d12_diff_FM_K0s","",10000,-0.01,0.11);

THnSparseD *hpeak_K0s_K0s_FU  = new THnSparseD("hpeak_K0s_K0s_FU","hpeak_K0s_K0s_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_FU  = new THnSparseD("hpeakd1_K0s_K0s_FU","hpeakd1_K0s_K0s_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_FU  = new THnSparseD("hpeakd2_K0s_K0s_FU","hpeakd2_K0s_K0s_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_FU  = new THnSparseD("hpeakd12_K0s_K0s_FU","hpeakd12_K0s_K0s_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_FU_bias  = new THnSparseD("hpeak_K0s_K0s_FU_bias","hpeak_K0s_K0s_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_FU_biasd1  = new THnSparseD("hpeak_K0s_K0s_FU_biasd1","hpeak_K0s_K0s_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_FU_biasd2  = new THnSparseD("hpeak_K0s_K0s_FU_biasd2","hpeak_K0s_K0s_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_FU_biasd12  = new THnSparseD("hpeak_K0s_K0s_FU_biasd12","hpeak_K0s_K0s_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FU_K0s=new TH1D("V0chi2d1_diff_FU_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FU_K0s=new TH1D("V0chi2d2_diff_FU_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FU_K0s=new TH1D("V0chi2d12_diff_FU_K0s","",10000,-0.01,0.11);

//true matched + true unmatched

THnSparseD *hpeak_K0s_K0s_TM_TU  = new THnSparseD("hpeak_K0s_K0s_TM_TU","hpeak_K0s_K0s_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_TM_TU  = new THnSparseD("hpeakd1_K0s_K0s_TM_TU","hpeakd1_K0s_K0s_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_TM_TU  = new THnSparseD("hpeakd2_K0s_K0s_TM_TU","hpeakd2_K0s_K0s_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_TM_TU  = new THnSparseD("hpeakd12_K0s_K0s_TM_TU","hpeakd12_K0s_K0s_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_TM_TU_bias  = new THnSparseD("hpeak_K0s_K0s_TM_TU_bias","hpeak_K0s_K0s_TM_TU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TM_TU_biasd1  = new THnSparseD("hpeak_K0s_K0s_TM_TU_biasd1","hpeak_K0s_K0s_TM_TU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TM_TU_biasd2  = new THnSparseD("hpeak_K0s_K0s_TM_TU_biasd2","hpeak_K0s_K0s_TM_TU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TM_TU_biasd12  = new THnSparseD("hpeak_K0s_K0s_TM_TU_biasd12","hpeak_K0s_K0s_TM_TU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_TU_K0s=new TH1D("V0chi2d1_diff_TM_TU_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_TU_K0s=new TH1D("V0chi2d2_diff_TM_TU_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_TU_K0s=new TH1D("V0chi2d12_diff_TM_TU_K0s","",10000,-0.01,0.11);

//fake matched + fake unmatched

THnSparseD *hpeak_K0s_K0s_FM_FU  = new THnSparseD("hpeak_K0s_K0s_FM_FU","hpeak_K0s_K0s_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_FM_FU  = new THnSparseD("hpeakd1_K0s_K0s_FM_FU","hpeakd1_K0s_K0s_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_FM_FU  = new THnSparseD("hpeakd2_K0s_K0s_FM_FU","hpeakd2_K0s_K0s_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_FM_FU  = new THnSparseD("hpeakd12_K0s_K0s_FM_FU","hpeakd12_K0s_K0s_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_FM_FU_bias  = new THnSparseD("hpeak_K0s_K0s_FM_FU_bias","hpeak_K0s_K0s_FM_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_FM_FU_biasd1  = new THnSparseD("hpeak_K0s_K0s_FM_FU_biasd1","hpeak_K0s_K0s_FM_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_FM_FU_biasd2  = new THnSparseD("hpeak_K0s_K0s_FM_FU_biasd2","hpeak_K0s_K0s_FM_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_FM_FU_biasd12  = new THnSparseD("hpeak_K0s_K0s_FM_FU_biasd12","hpeak_K0s_K0s_FM_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FM_FU_K0s=new TH1D("V0chi2d1_diff_FM_FU_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FM_FU_K0s=new TH1D("V0chi2d2_diff_FM_FU_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FM_FU_K0s=new TH1D("V0chi2d12_diff_FM_FU_K0s","",10000,-0.01,0.11);


//true matched + fake matched

THnSparseD *hpeak_K0s_K0s_TM_FM  = new THnSparseD("hpeak_K0s_K0s_TM_FM","hpeak_K0s_K0s_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_TM_FM  = new THnSparseD("hpeakd1_K0s_K0s_TM_FM","hpeakd1_K0s_K0s_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_TM_FM  = new THnSparseD("hpeakd2_K0s_K0s_TM_FM","hpeakd2_K0s_K0s_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_TM_FM  = new THnSparseD("hpeakd12_K0s_K0s_TM_FM","hpeakd12_K0s_K0s_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_TM_FM_bias  = new THnSparseD("hpeak_K0s_K0s_TM_FM_bias","hpeak_K0s_K0s_TM_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TM_FM_biasd1  = new THnSparseD("hpeak_K0s_K0s_TM_FM_biasd1","hpeak_K0s_K0s_TM_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TM_FM_biasd2  = new THnSparseD("hpeak_K0s_K0s_TM_FM_biasd2","hpeak_K0s_K0s_TM_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TM_FM_biasd12  = new THnSparseD("hpeak_K0s_K0s_TM_FM_biasd12","hpeak_K0s_K0s_TM_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_FM_K0s=new TH1D("V0chi2d1_diff_TM_FM_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_FM_K0s=new TH1D("V0chi2d2_diff_TM_FM_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_FM_K0s=new TH1D("V0chi2d12_diff_TM_FM_K0s","",10000,-0.01,0.11);

//true matched + fake unmatched

THnSparseD *hpeak_K0s_K0s_TM_FU  = new THnSparseD("hpeak_K0s_K0s_TM_FU","hpeak_K0s_K0s_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_TM_FU  = new THnSparseD("hpeakd1_K0s_K0s_TM_FU","hpeakd1_K0s_K0s_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_TM_FU  = new THnSparseD("hpeakd2_K0s_K0s_TM_FU","hpeakd2_K0s_K0s_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_TM_FU  = new THnSparseD("hpeakd12_K0s_K0s_TM_FU","hpeakd12_K0s_K0s_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_TM_FU_bias  = new THnSparseD("hpeak_K0s_K0s_TM_FU_bias","hpeak_K0s_K0s_TM_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TM_FU_biasd1  = new THnSparseD("hpeak_K0s_K0s_TM_FU_biasd1","hpeak_K0s_K0s_TM_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TM_FU_biasd2  = new THnSparseD("hpeak_K0s_K0s_TM_FU_biasd2","hpeak_K0s_K0s_TM_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TM_FU_biasd12  = new THnSparseD("hpeak_K0s_K0s_TM_FU_biasd12","hpeak_K0s_K0s_TM_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_FU_K0s=new TH1D("V0chi2d1_diff_TM_FU_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_FU_K0s=new TH1D("V0chi2d2_diff_TM_FU_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_FU_K0s=new TH1D("V0chi2d12_diff_TM_FU_K0s","",10000,-0.01,0.11);

//true ummatched + fake unmatched

THnSparseD *hpeak_K0s_K0s_TU_FU  = new THnSparseD("hpeak_K0s_K0s_TU_FU","hpeak_K0s_K0s_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_TU_FU  = new THnSparseD("hpeakd1_K0s_K0s_TU_FU","hpeakd1_K0s_K0s_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_TU_FU  = new THnSparseD("hpeakd2_K0s_K0s_TU_FU","hpeakd2_K0s_K0s_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_TU_FU  = new THnSparseD("hpeakd12_K0s_K0s_TU_FU","hpeakd12_K0s_K0s_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_TU_FU_bias  = new THnSparseD("hpeak_K0s_K0s_TU_FU_bias","hpeak_K0s_K0s_TU_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TU_FU_biasd1  = new THnSparseD("hpeak_K0s_K0s_TU_FU_biasd1","hpeak_K0s_K0s_TU_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TU_FU_biasd2  = new THnSparseD("hpeak_K0s_K0s_TU_FU_biasd2","hpeak_K0s_K0s_TU_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TU_FU_biasd12  = new THnSparseD("hpeak_K0s_K0s_TU_FU_biasd12","hpeak_K0s_K0s_TU_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_FU_K0s=new TH1D("V0chi2d1_diff_TU_FU_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_FU_K0s=new TH1D("V0chi2d2_diff_TU_FU_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_FU_K0s=new TH1D("V0chi2d12_diff_TU_FU_K0s","",10000,-0.01,0.11);

//true unmatched + fake matched

THnSparseD *hpeak_K0s_K0s_TU_FM  = new THnSparseD("hpeak_K0s_K0s_TU_FM","hpeak_K0s_K0s_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_K0s_TU_FM  = new THnSparseD("hpeakd1_K0s_K0s_TU_FM","hpeakd1_K0s_K0s_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_K0s_TU_FM  = new THnSparseD("hpeakd2_K0s_K0s_TU_FM","hpeakd2_K0s_K0s_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_K0s_TU_FM  = new THnSparseD("hpeakd12_K0s_K0s_TU_FM","hpeakd12_K0s_K0s_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_TU_FM_bias  = new THnSparseD("hpeak_K0s_K0s_TU_FM_bias","hpeak_K0s_K0s_TU_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TU_FM_biasd1  = new THnSparseD("hpeak_K0s_K0s_TU_FM_biasd1","hpeak_K0s_K0s_TU_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TU_FM_biasd2  = new THnSparseD("hpeak_K0s_K0s_TU_FM_biasd2","hpeak_K0s_K0s_TU_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_K0s_TU_FM_biasd12  = new THnSparseD("hpeak_K0s_K0s_TU_FM_biasd12","hpeak_K0s_K0s_TU_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_FM_K0s=new TH1D("V0chi2d1_diff_TU_FM_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_FM_K0s=new TH1D("V0chi2d2_diff_TU_FM_K0s","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_FM_K0s=new TH1D("V0chi2d12_diff_TU_FM_K0s","",10000,-0.01,0.11);


//Lam
//True

THnSparseD *hMass_Lam_Lam_T  = new THnSparseD("hMass_Lam_Lam_T","hMass_Lam_Lam_T",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_Lam_1_T  = new THnSparseD("hMass_Lam_1_T","hMass_Lam_1_T",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_Lam_2_T  = new THnSparseD("hMass_Lam_2_T","hMass_Lam_2_T",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   

 
THnSparseD *hpeak_Lam_Lam_T  = new THnSparseD("hpeak_Lam_Lam_T","hpeak_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakcos_Lam_Lam_T  = new THnSparseD("hpeakcos_Lam_Lam_T","hpeakcos_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_rot_Lam_Lam_T  = new THnSparseD("hpeak_rot_Lam_Lam_T","hpeak_rot_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_inv_Lam_Lam_T  = new THnSparseD("hpeak_inv_Lam_Lam_T","hpeak_inv_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_Lam_Lam_T  = new THnSparseD("hside_Lam_Lam_T","hside_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideL_Lam_Lam_T  = new THnSparseD("hsideL_Lam_Lam_T","hsideL_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideR_Lam_Lam_T  = new THnSparseD("hsideR_Lam_Lam_T","hsideR_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_rot_Lam_Lam_T  = new THnSparseD("hside_rot_Lam_Lam_T","hside_rot_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_inv_Lam_Lam_T  = new THnSparseD("hside_inv_Lam_Lam_T","hside_inv_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_Lam_Lam_T  = new THnSparseD("hpeakside_Lam_Lam_T","hpeakside_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideL_Lam_Lam_T  = new THnSparseD("hpeaksideL_Lam_Lam_T","hpeaksideL_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideR_Lam_Lam_T  = new THnSparseD("hpeaksideR_Lam_Lam_T","hpeaksideR_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_rot_Lam_Lam_T  = new THnSparseD("hpeakside_rot_Lam_Lam_T","hpeakside_rot_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_inv_Lam_Lam_T  = new THnSparseD("hpeakside_inv_Lam_Lam_T","hpeakside_inv_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_T  = new THnSparseD("hpeakd1_Lam_Lam_T","hpeakd1_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_T  = new THnSparseD("hpeakd2_Lam_Lam_T","hpeakd2_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_T  = new THnSparseD("hpeakd12_Lam_Lam_T","hpeakd12_Lam_Lam_T",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_Lam_Lam_T_mix  = new THnSparseD("hpeak_Lam_Lam_T_mix","hpeak_Lam_Lam_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_Lam_Lam_T_mix  = new THnSparseD("hside_Lam_Lam_T_mix","hside_Lam_Lam_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideL_Lam_Lam_T_mix  = new THnSparseD("hsideL_Lam_Lam_T_mix","hsideL_Lam_Lam_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideR_Lam_Lam_T_mix  = new THnSparseD("hsideR_Lam_Lam_T_mix","hsideR_Lam_Lam_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_Lam_Lam_T_mix  = new THnSparseD("hpeakside_Lam_Lam_T_mix","hpeakside_Lam_Lam_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideL_Lam_Lam_T_mix  = new THnSparseD("hpeaksideL_Lam_Lam_T_mix","hpeaksideL_Lam_Lam_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideR_Lam_Lam_T_mix  = new THnSparseD("hpeaksideR_Lam_Lam_T_mix","hpeaksideR_Lam_Lam_T_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_Lam_Lam_T_etamix  = new THnSparseD("hpeak_Lam_Lam_T_etamix","hpeak_Lam_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_Lam_Lam_T_etamix  = new THnSparseD("hside_Lam_Lam_T_etamix","hside_Lam_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideL_Lam_Lam_T_etamix  = new THnSparseD("hsideL_Lam_Lam_T_etamix","hsideL_Lam_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideR_Lam_Lam_T_etamix  = new THnSparseD("hsideR_Lam_Lam_T_etamix","hsideR_Lam_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_Lam_Lam_T_etamix  = new THnSparseD("hpeakside_Lam_Lam_T_etamix","hpeakside_Lam_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideL_Lam_Lam_T_etamix  = new THnSparseD("hpeaksideL_Lam_Lam_T_etamix","hpeaksideL_Lam_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideR_Lam_Lam_T_etamix  = new THnSparseD("hpeaksideR_Lam_Lam_T_etamix","hpeaksideR_Lam_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);

TH1D *V0chi2d1_diff_T_Lam=new TH1D("V0chi2d1_diff_T_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_T_Lam=new TH1D("V0chi2d2_diff_T_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_T_Lam=new TH1D("V0chi2d12_diff_T_Lam","",10000,-0.01,0.11);
  
//False

THnSparseD *hMass_Lam_Lam_F  = new THnSparseD("hMass_Lam_Lam_F","hMass_Lam_Lam_F",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hpeak_Lam_Lam_F  = new THnSparseD("hpeak_Lam_Lam_F","hpeak_Lam_Lam_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_Lam_Lam_F  = new THnSparseD("hside_Lam_Lam_F","hside_Lam_Lam_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_Lam_Lam_F  = new THnSparseD("hpeakside_Lam_Lam_F","hpeakside_Lam_Lam_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_F  = new THnSparseD("hpeakd1_Lam_Lam_F","hpeakd1_Lam_Lam_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_F  = new THnSparseD("hpeakd2_Lam_Lam_F","hpeakd2_Lam_Lam_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_F  = new THnSparseD("hpeakd12_Lam_Lam_F","hpeakd12_Lam_Lam_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_F_bias  = new THnSparseD("hpeak_Lam_Lam_F_bias","hpeak_Lam_Lam_F_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_F_biasd1  = new THnSparseD("hpeak_Lam_Lam_F_biasd1","hpeak_Lam_Lam_F_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_F_biasd2  = new THnSparseD("hpeak_Lam_Lam_F_biasd2","hpeak_Lam_Lam_F_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_F_biasd12  = new THnSparseD("hpeak_Lam_Lam_F_biasd12","hpeak_Lam_Lam_F_biasd12",3,bins3Db,xmin3Db,xmax3Db);


TH1D *V0chi2d1_diff_F_Lam=new TH1D("V0chi2d1_diff_F_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_F_Lam=new TH1D("V0chi2d2_diff_F_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_F_Lam=new TH1D("V0chi2d12_diff_F_Lam","",10000,-0.01,0.11);
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
THnSparseD *hpeak_Lam_Lam_TF  = new THnSparseD("hpeak_Lam_Lam_TF","hpeak_Lam_Lam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_TF  = new THnSparseD("hpeakd1_Lam_Lam_TF","hpeakd1_Lam_Lam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_TF  = new THnSparseD("hpeakd2_Lam_Lam_TF","hpeakd2_Lam_Lam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_TF  = new THnSparseD("hpeakd12_Lam_Lam_TF","hpeakd12_Lam_Lam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_Lam_Lam_TF  = new THnSparseD("hside_Lam_Lam_TF","hside_Lam_Lam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_Lam_Lam_TF  = new THnSparseD("hpeakside_Lam_Lam_TF","hpeakside_Lam_Lam_TF",3,bins3D,xmin3D,xmax3D);  
THnSparseD *hpeak_Lam_Lam_TF_bias  = new THnSparseD("hpeak_Lam_Lam_TF_bias","hpeak_Lam_Lam_TF_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TF_biasd1  = new THnSparseD("hpeak_Lam_Lam_TF_biasd1","hpeak_Lam_Lam_TF_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TF_biasd2  = new THnSparseD("hpeak_Lam_Lam_TF_biasd2","hpeak_Lam_Lam_TF_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TF_biasd12  = new THnSparseD("hpeak_Lam_Lam_TF_biasd12","hpeak_Lam_Lam_TF_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TF_Lam=new TH1D("V0chi2d1_diff_TF_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TF_Lam=new TH1D("V0chi2d2_diff_TF_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TF_Lam=new TH1D("V0chi2d12_diff_TF_Lam","",10000,-0.01,0.11);
  
  
/////// matching hist ////// 
  
THnSparseD *hMass_Lam_Lam_TM  = new THnSparseD("hMass_Lam_Lam_TM","hMass_Lam_Lam_TM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_Lam_Lam_TU  = new THnSparseD("hMass_Lam_Lam_TU","hMass_Lam_Lam_TU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_Lam_Lam_FM  = new THnSparseD("hMass_Lam_Lam_FM","hMass_Lam_Lam_FM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_Lam_Lam_FU  = new THnSparseD("hMass_Lam_Lam_FU","hMass_Lam_Lam_FU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
 
THnSparseD *hpeak_Lam_Lam_TM  = new THnSparseD("hpeak_Lam_Lam_TM","hpeak_Lam_Lam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_TM  = new THnSparseD("hpeakd1_Lam_Lam_TM","hpeakd1_Lam_Lam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_TM  = new THnSparseD("hpeakd2_Lam_Lam_TM","hpeakd2_Lam_Lam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_TM  = new THnSparseD("hpeakd12_Lam_Lam_TM","hpeakd12_Lam_Lam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_TM_bias  = new THnSparseD("hpeak_Lam_Lam_TM_bias","hpeak_Lam_Lam_TM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TM_biasd1  = new THnSparseD("hpeak_Lam_Lam_TM_biasd1","hpeak_Lam_Lam_TM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TM_biasd2  = new THnSparseD("hpeak_Lam_Lam_TM_biasd2","hpeak_Lam_Lam_TM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TM_biasd12  = new THnSparseD("hpeak_Lam_Lam_TM_biasd12","hpeak_Lam_Lam_TM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_Lam=new TH1D("V0chi2d1_diff_TM_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_Lam=new TH1D("V0chi2d2_diff_TM_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_Lam=new TH1D("V0chi2d12_diff_TM_Lam","",10000,-0.01,0.11);

THnSparseD *hpeak_Lam_Lam_TU  = new THnSparseD("hpeak_Lam_Lam_TU","hpeak_Lam_Lam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_TU  = new THnSparseD("hpeakd1_Lam_Lam_TU","hpeakd1_Lam_Lam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_TU  = new THnSparseD("hpeakd2_Lam_Lam_TU","hpeakd2_Lam_Lam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_TU  = new THnSparseD("hpeakd12_Lam_Lam_TU","hpeakd12_Lam_Lam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_TU_bias  = new THnSparseD("hpeak_Lam_Lam_TU_bias","hpeak_Lam_Lam_TU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TU_biasd1  = new THnSparseD("hpeak_Lam_Lam_TU_biasd1","hpeak_Lam_Lam_TU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TU_biasd2  = new THnSparseD("hpeak_Lam_Lam_TU_biasd2","hpeak_Lam_Lam_TU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TU_biasd12  = new THnSparseD("hpeak_Lam_Lam_TU_biasd12","hpeak_Lam_Lam_TU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_Lam=new TH1D("V0chi2d1_diff_TU_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_Lam=new TH1D("V0chi2d2_diff_TU_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_Lam=new TH1D("V0chi2d12_diff_TU_Lam","",10000,-0.01,0.11);

THnSparseD *hpeak_Lam_Lam_FM  = new THnSparseD("hpeak_Lam_Lam_FM","hpeak_Lam_Lam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_FM  = new THnSparseD("hpeakd1_Lam_Lam_FM","hpeakd1_Lam_Lam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_FM  = new THnSparseD("hpeakd2_Lam_Lam_FM","hpeakd2_Lam_Lam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_FM  = new THnSparseD("hpeakd12_Lam_Lam_FM","hpeakd12_Lam_Lam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_FM_bias  = new THnSparseD("hpeak_Lam_Lam_FM_bias","hpeak_Lam_Lam_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_FM_biasd1  = new THnSparseD("hpeak_Lam_Lam_FM_biasd1","hpeak_Lam_Lam_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_FM_biasd2  = new THnSparseD("hpeak_Lam_Lam_FM_biasd2","hpeak_Lam_Lam_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_FM_biasd12  = new THnSparseD("hpeak_Lam_Lam_FM_biasd12","hpeak_Lam_Lam_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FM_Lam=new TH1D("V0chi2d1_diff_FM_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FM_Lam=new TH1D("V0chi2d2_diff_FM_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FM_Lam=new TH1D("V0chi2d12_diff_FM_Lam","",10000,-0.01,0.11);

THnSparseD *hpeak_Lam_Lam_FU  = new THnSparseD("hpeak_Lam_Lam_FU","hpeak_Lam_Lam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_FU  = new THnSparseD("hpeakd1_Lam_Lam_FU","hpeakd1_Lam_Lam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_FU  = new THnSparseD("hpeakd2_Lam_Lam_FU","hpeakd2_Lam_Lam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_FU  = new THnSparseD("hpeakd12_Lam_Lam_FU","hpeakd12_Lam_Lam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_FU_bias  = new THnSparseD("hpeak_Lam_Lam_FU_bias","hpeak_Lam_Lam_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_FU_biasd1  = new THnSparseD("hpeak_Lam_Lam_FU_biasd1","hpeak_Lam_Lam_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_FU_biasd2  = new THnSparseD("hpeak_Lam_Lam_FU_biasd2","hpeak_Lam_Lam_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_FU_biasd12  = new THnSparseD("hpeak_Lam_Lam_FU_biasd12","hpeak_Lam_Lam_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FU_Lam=new TH1D("V0chi2d1_diff_FU_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FU_Lam=new TH1D("V0chi2d2_diff_FU_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FU_Lam=new TH1D("V0chi2d12_diff_FU_Lam","",10000,-0.01,0.11);

//true matched + true unmatched

THnSparseD *hpeak_Lam_Lam_TM_TU  = new THnSparseD("hpeak_Lam_Lam_TM_TU","hpeak_Lam_Lam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_TM_TU  = new THnSparseD("hpeakd1_Lam_Lam_TM_TU","hpeakd1_Lam_Lam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_TM_TU  = new THnSparseD("hpeakd2_Lam_Lam_TM_TU","hpeakd2_Lam_Lam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_TM_TU  = new THnSparseD("hpeakd12_Lam_Lam_TM_TU","hpeakd12_Lam_Lam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_TM_TU_bias  = new THnSparseD("hpeak_Lam_Lam_TM_TU_bias","hpeak_Lam_Lam_TM_TU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TM_TU_biasd1  = new THnSparseD("hpeak_Lam_Lam_TM_TU_biasd1","hpeak_Lam_Lam_TM_TU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TM_TU_biasd2  = new THnSparseD("hpeak_Lam_Lam_TM_TU_biasd2","hpeak_Lam_Lam_TM_TU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TM_TU_biasd12  = new THnSparseD("hpeak_Lam_Lam_TM_TU_biasd12","hpeak_Lam_Lam_TM_TU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_TU_Lam=new TH1D("V0chi2d1_diff_TM_TU_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_TU_Lam=new TH1D("V0chi2d2_diff_TM_TU_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_TU_Lam=new TH1D("V0chi2d12_diff_TM_TU_Lam","",10000,-0.01,0.11);

//fake matched + fake unmatched

THnSparseD *hpeak_Lam_Lam_FM_FU  = new THnSparseD("hpeak_Lam_Lam_FM_FU","hpeak_Lam_Lam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_FM_FU  = new THnSparseD("hpeakd1_Lam_Lam_FM_FU","hpeakd1_Lam_Lam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_FM_FU  = new THnSparseD("hpeakd2_Lam_Lam_FM_FU","hpeakd2_Lam_Lam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_FM_FU  = new THnSparseD("hpeakd12_Lam_Lam_FM_FU","hpeakd12_Lam_Lam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_FM_FU_bias  = new THnSparseD("hpeak_Lam_Lam_FM_FU_bias","hpeak_Lam_Lam_FM_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_FM_FU_biasd1  = new THnSparseD("hpeak_Lam_Lam_FM_FU_biasd1","hpeak_Lam_Lam_FM_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_FM_FU_biasd2  = new THnSparseD("hpeak_Lam_Lam_FM_FU_biasd2","hpeak_Lam_Lam_FM_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_FM_FU_biasd12  = new THnSparseD("hpeak_Lam_Lam_FM_FU_biasd12","hpeak_Lam_Lam_FM_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FM_FU_Lam=new TH1D("V0chi2d1_diff_FM_FU_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FM_FU_Lam=new TH1D("V0chi2d2_diff_FM_FU_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FM_FU_Lam=new TH1D("V0chi2d12_diff_FM_FU_Lam","",10000,-0.01,0.11);


THnSparseD *hpeak_Lam_Lam_TM_FM  = new THnSparseD("hpeak_Lam_Lam_TM_FM","hpeak_Lam_Lam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_TM_FM  = new THnSparseD("hpeakd1_Lam_Lam_TM_FM","hpeakd1_Lam_Lam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_TM_FM  = new THnSparseD("hpeakd2_Lam_Lam_TM_FM","hpeakd2_Lam_Lam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_TM_FM  = new THnSparseD("hpeakd12_Lam_Lam_TM_FM","hpeakd12_Lam_Lam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_TM_FM_bias  = new THnSparseD("hpeak_Lam_Lam_TM_FM_bias","hpeak_Lam_Lam_TM_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TM_FM_biasd1  = new THnSparseD("hpeak_Lam_Lam_TM_FM_biasd1","hpeak_Lam_Lam_TM_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TM_FM_biasd2  = new THnSparseD("hpeak_Lam_Lam_TM_FM_biasd2","hpeak_Lam_Lam_TM_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TM_FM_biasd12  = new THnSparseD("hpeak_Lam_Lam_TM_FM_biasd12","hpeak_Lam_Lam_TM_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_FM_Lam=new TH1D("V0chi2d1_diff_TM_FM_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_FM_Lam=new TH1D("V0chi2d2_diff_TM_FM_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_FM_Lam=new TH1D("V0chi2d12_diff_TM_FM_Lam","",10000,-0.01,0.11);

THnSparseD *hpeak_Lam_Lam_TM_FU  = new THnSparseD("hpeak_Lam_Lam_TM_FU","hpeak_Lam_Lam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_TM_FU  = new THnSparseD("hpeakd1_Lam_Lam_TM_FU","hpeakd1_Lam_Lam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_TM_FU  = new THnSparseD("hpeakd2_Lam_Lam_TM_FU","hpeakd2_Lam_Lam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_TM_FU  = new THnSparseD("hpeakd12_Lam_Lam_TM_FU","hpeakd12_Lam_Lam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_TM_FU_bias  = new THnSparseD("hpeak_Lam_Lam_TM_FU_bias","hpeak_Lam_Lam_TM_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TM_FU_biasd1  = new THnSparseD("hpeak_Lam_Lam_TM_FU_biasd1","hpeak_Lam_Lam_TM_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TM_FU_biasd2  = new THnSparseD("hpeak_Lam_Lam_TM_FU_biasd2","hpeak_Lam_Lam_TM_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TM_FU_biasd12  = new THnSparseD("hpeak_Lam_Lam_TM_FU_biasd12","hpeak_Lam_Lam_TM_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_FU_Lam=new TH1D("V0chi2d1_diff_TM_FU_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_FU_Lam=new TH1D("V0chi2d2_diff_TM_FU_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_FU_Lam=new TH1D("V0chi2d12_diff_TM_FU_Lam","",10000,-0.01,0.11);

THnSparseD *hpeak_Lam_Lam_TU_FU  = new THnSparseD("hpeak_Lam_Lam_TU_FU","hpeak_Lam_Lam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_TU_FU  = new THnSparseD("hpeakd1_Lam_Lam_TU_FU","hpeakd1_Lam_Lam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_TU_FU  = new THnSparseD("hpeakd2_Lam_Lam_TU_FU","hpeakd2_Lam_Lam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_TU_FU  = new THnSparseD("hpeakd12_Lam_Lam_TU_FU","hpeakd12_Lam_Lam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_TU_FU_bias  = new THnSparseD("hpeak_Lam_Lam_TU_FU_bias","hpeak_Lam_Lam_TU_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TU_FU_biasd1  = new THnSparseD("hpeak_Lam_Lam_TU_FU_biasd1","hpeak_Lam_Lam_TU_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TU_FU_biasd2  = new THnSparseD("hpeak_Lam_Lam_TU_FU_biasd2","hpeak_Lam_Lam_TU_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TU_FU_biasd12  = new THnSparseD("hpeak_Lam_Lam_TU_FU_biasd12","hpeak_Lam_Lam_TU_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_FU_Lam=new TH1D("V0chi2d1_diff_TU_FU_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_FU_Lam=new TH1D("V0chi2d2_diff_TU_FU_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_FU_Lam=new TH1D("V0chi2d12_diff_TU_FU_Lam","",10000,-0.01,0.11);

THnSparseD *hpeak_Lam_Lam_TU_FM  = new THnSparseD("hpeak_Lam_Lam_TU_FM","hpeak_Lam_Lam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_Lam_TU_FM  = new THnSparseD("hpeakd1_Lam_Lam_TU_FM","hpeakd1_Lam_Lam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_Lam_TU_FM  = new THnSparseD("hpeakd2_Lam_Lam_TU_FM","hpeakd2_Lam_Lam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_Lam_TU_FM  = new THnSparseD("hpeakd12_Lam_Lam_TU_FM","hpeakd12_Lam_Lam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_TU_FM_bias  = new THnSparseD("hpeak_Lam_Lam_TU_FM_bias","hpeak_Lam_Lam_TU_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TU_FM_biasd1  = new THnSparseD("hpeak_Lam_Lam_TU_FM_biasd1","hpeak_Lam_Lam_TU_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TU_FM_biasd2  = new THnSparseD("hpeak_Lam_Lam_TU_FM_biasd2","hpeak_Lam_Lam_TU_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_Lam_TU_FM_biasd12  = new THnSparseD("hpeak_Lam_Lam_TU_FM_biasd12","hpeak_Lam_Lam_TU_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_FM_Lam=new TH1D("V0chi2d1_diff_TU_FM_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_FM_Lam=new TH1D("V0chi2d2_diff_TU_FM_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_FM_Lam=new TH1D("V0chi2d12_diff_TU_FM_Lam","",10000,-0.01,0.11);

//ALam
//True
THnSparseD *hMass_ALam_ALam_T  = new THnSparseD("hMass_ALam_ALam_T","hMass_ALam_ALam_T",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_ALam_1_T  = new THnSparseD("hMass_ALam_1_T","hMass_ALam_1_T",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_ALam_2_T  = new THnSparseD("hMass_ALam_2_T","hMass_ALam_2_T",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   


THnSparseD *hpeak_ALam_ALam_T  = new THnSparseD("hpeak_ALam_ALam_T","hpeak_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakcos_ALam_ALam_T  = new THnSparseD("hpeakcos_ALam_ALam_T","hpeakcos_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_rot_ALam_ALam_T  = new THnSparseD("hpeak_rot_ALam_ALam_T","hpeak_rot_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_inv_ALam_ALam_T  = new THnSparseD("hpeak_inv_ALam_ALam_T","hpeak_inv_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_ALam_ALam_T  = new THnSparseD("hside_ALam_ALam_T","hside_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideL_ALam_ALam_T  = new THnSparseD("hsideL_ALam_ALam_T","hsideL_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideR_ALam_ALam_T  = new THnSparseD("hsideR_ALam_ALam_T","hsideR_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_rot_ALam_ALam_T  = new THnSparseD("hside_rot_ALam_ALam_T","hside_rot_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_inv_ALam_ALam_T  = new THnSparseD("hside_inv_ALam_ALam_T","hside_inv_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_ALam_ALam_T  = new THnSparseD("hpeakside_ALam_ALam_T","hpeakside_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideL_ALam_ALam_T  = new THnSparseD("hpeaksideL_ALam_ALam_T","hpeaksideL_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideR_ALam_ALam_T  = new THnSparseD("hpeaksideR_ALam_ALam_T","hpeaksideR_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_rot_ALam_ALam_T  = new THnSparseD("hpeakside_rot_ALam_ALam_T","hpeakside_rot_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_inv_ALam_ALam_T  = new THnSparseD("hpeakside_inv_ALam_ALam_T","hpeakside_inv_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_T  = new THnSparseD("hpeakd1_ALam_ALam_T","hpeakd1_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_T  = new THnSparseD("hpeakd2_ALam_ALam_T","hpeakd2_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_T  = new THnSparseD("hpeakd12_ALam_ALam_T","hpeakd12_ALam_ALam_T",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_ALam_ALam_T_mix  = new THnSparseD("hpeak_ALam_ALam_T_mix","hpeak_ALam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_ALam_ALam_T_mix  = new THnSparseD("hside_ALam_ALam_T_mix","hside_ALam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideL_ALam_ALam_T_mix  = new THnSparseD("hsideL_ALam_ALam_T_mix","hsideL_ALam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideR_ALam_ALam_T_mix  = new THnSparseD("hsideR_ALam_ALam_T_mix","hsideR_ALam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_ALam_ALam_T_mix  = new THnSparseD("hpeakside_ALam_ALam_T_mix","hpeakside_ALam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideL_ALam_ALam_T_mix  = new THnSparseD("hpeaksideL_ALam_ALam_T_mix","hpeaksideL_ALam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideR_ALam_ALam_T_mix  = new THnSparseD("hpeaksideR_ALam_ALam_T_mix","hpeaksideR_ALam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_ALam_ALam_T_etamix  = new THnSparseD("hpeak_ALam_ALam_T_etamix","hpeak_ALam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_ALam_ALam_T_etamix  = new THnSparseD("hside_ALam_ALam_T_etamix","hside_ALam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideL_ALam_ALam_T_etamix  = new THnSparseD("hsideL_ALam_ALam_T_etamix","hsideL_ALam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hsideR_ALam_ALam_T_etamix  = new THnSparseD("hsideR_ALam_ALam_T_etamix","hsideR_ALam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_ALam_ALam_T_etamix  = new THnSparseD("hpeakside_ALam_ALam_T_etamix","hpeakside_ALam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideL_ALam_ALam_T_etamix  = new THnSparseD("hpeaksideL_ALam_ALam_T_etamix","hpeaksideL_ALam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeaksideR_ALam_ALam_T_etamix  = new THnSparseD("hpeaksideR_ALam_ALam_T_etamix","hpeaksideR_ALam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);

TH1D *V0chi2d1_diff_T_ALam=new TH1D("V0chi2d1_diff_T_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_T_ALam=new TH1D("V0chi2d2_diff_T_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_T_ALam=new TH1D("V0chi2d12_diff_T_ALam","",10000,-0.01,0.11);

//False

THnSparseD *hMass_ALam_ALam_F  = new THnSparseD("hMass_ALam_ALam_F","hMass_ALam_ALam_F",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   

THnSparseD *hpeak_ALam_ALam_F  = new THnSparseD("hpeak_ALam_ALam_F","hpeak_ALam_ALam_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_ALam_ALam_F  = new THnSparseD("hside_ALam_ALam_F","hside_ALam_ALam_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_ALam_ALam_F  = new THnSparseD("hpeakside_ALam_ALam_F","hpeakside_ALam_ALam_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_F  = new THnSparseD("hpeakd1_ALam_ALam_F","hpeakd1_ALam_ALam_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_F  = new THnSparseD("hpeakd2_ALam_ALam_F","hpeakd2_ALam_ALam_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_F  = new THnSparseD("hpeakd12_ALam_ALam_F","hpeakd12_ALam_ALam_F",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_F_bias  = new THnSparseD("hpeak_ALam_ALam_F_bias","hpeak_ALam_ALam_F_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_F_biasd1  = new THnSparseD("hpeak_ALam_ALam_F_biasd1","hpeak_ALam_ALam_F_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_F_biasd2  = new THnSparseD("hpeak_ALam_ALam_F_biasd2","hpeak_ALam_ALam_F_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_F_biasd12  = new THnSparseD("hpeak_ALam_ALam_F_biasd12","hpeak_ALam_ALam_F_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_F_ALam=new TH1D("V0chi2d1_diff_F_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_F_ALam=new TH1D("V0chi2d2_diff_F_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_F_ALam=new TH1D("V0chi2d12_diff_F_ALam","",10000,-0.01,0.11);
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
THnSparseD *hpeak_ALam_ALam_TF  = new THnSparseD("hpeak_ALam_ALam_TF","hpeak_ALam_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_TF  = new THnSparseD("hpeakd1_ALam_ALam_TF","hpeakd1_ALam_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_TF  = new THnSparseD("hpeakd2_ALam_ALam_TF","hpeakd2_ALam_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_TF  = new THnSparseD("hpeakd12_ALam_ALam_TF","hpeakd12_ALam_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_ALam_ALam_TF  = new THnSparseD("hside_ALam_ALam_TF","hside_ALam_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_ALam_ALam_TF  = new THnSparseD("hpeakside_ALam_ALam_TF","hpeakside_ALam_ALam_TF",3,bins3D,xmin3D,xmax3D);  
THnSparseD *hpeak_ALam_ALam_TF_bias  = new THnSparseD("hpeak_ALam_ALam_TF_bias","hpeak_ALam_ALam_TF_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TF_biasd1  = new THnSparseD("hpeak_ALam_ALam_TF_biasd1","hpeak_ALam_ALam_TF_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TF_biasd2  = new THnSparseD("hpeak_ALam_ALam_TF_biasd2","hpeak_ALam_ALam_TF_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TF_biasd12  = new THnSparseD("hpeak_ALam_ALam_TF_biasd12","hpeak_ALam_ALam_TF_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TF_ALam=new TH1D("V0chi2d1_diff_TF_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TF_ALam=new TH1D("V0chi2d2_diff_TF_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TF_ALam=new TH1D("V0chi2d12_diff_TF_ALam","",10000,-0.01,0.11);

  
/////// matching hist ////// 
  
THnSparseD *hMass_ALam_ALam_TM  = new THnSparseD("hMass_ALam_ALam_TM","hMass_ALam_ALam_TM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_ALam_ALam_TU  = new THnSparseD("hMass_ALam_ALam_TU","hMass_ALam_ALam_TU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_ALam_ALam_FM  = new THnSparseD("hMass_ALam_ALam_FM","hMass_ALam_ALam_FM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_ALam_ALam_FU  = new THnSparseD("hMass_ALam_ALam_FU","hMass_ALam_ALam_FU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
 
THnSparseD *hpeak_ALam_ALam_TM  = new THnSparseD("hpeak_ALam_ALam_TM","hpeak_ALam_ALam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_TM  = new THnSparseD("hpeakd1_ALam_ALam_TM","hpeakd1_ALam_ALam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_TM  = new THnSparseD("hpeakd2_ALam_ALam_TM","hpeakd2_ALam_ALam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_TM  = new THnSparseD("hpeakd12_ALam_ALam_TM","hpeakd12_ALam_ALam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_TM_bias  = new THnSparseD("hpeak_ALam_ALam_TM_bias","hpeak_ALam_ALam_TM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TM_biasd1  = new THnSparseD("hpeak_ALam_ALam_TM_biasd1","hpeak_ALam_ALam_TM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TM_biasd2  = new THnSparseD("hpeak_ALam_ALam_TM_biasd2","hpeak_ALam_ALam_TM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TM_biasd12  = new THnSparseD("hpeak_ALam_ALam_TM_biasd12","hpeak_ALam_ALam_TM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_ALam=new TH1D("V0chi2d1_diff_TM_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_ALam=new TH1D("V0chi2d2_diff_TM_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_ALam=new TH1D("V0chi2d12_diff_TM_ALam","",10000,-0.01,0.11);

THnSparseD *hpeak_ALam_ALam_TU  = new THnSparseD("hpeak_ALam_ALam_TU","hpeak_ALam_ALam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_TU  = new THnSparseD("hpeakd1_ALam_ALam_TU","hpeakd1_ALam_ALam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_TU  = new THnSparseD("hpeakd2_ALam_ALam_TU","hpeakd2_ALam_ALam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_TU  = new THnSparseD("hpeakd12_ALam_ALam_TU","hpeakd12_ALam_ALam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_TU_bias  = new THnSparseD("hpeak_ALam_ALam_TU_bias","hpeak_ALam_ALam_TU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TU_biasd1  = new THnSparseD("hpeak_ALam_ALam_TU_biasd1","hpeak_ALam_ALam_TU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TU_biasd2  = new THnSparseD("hpeak_ALam_ALam_TU_biasd2","hpeak_ALam_ALam_TU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TU_biasd12  = new THnSparseD("hpeak_ALam_ALam_TU_biasd12","hpeak_ALam_ALam_TU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_ALam=new TH1D("V0chi2d1_diff_TU_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_ALam=new TH1D("V0chi2d2_diff_TU_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_ALam=new TH1D("V0chi2d12_diff_TU_ALam","",10000,-0.01,0.11);

THnSparseD *hpeak_ALam_ALam_FM  = new THnSparseD("hpeak_ALam_ALam_FM","hpeak_ALam_ALam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_FM  = new THnSparseD("hpeakd1_ALam_ALam_FM","hpeakd1_ALam_ALam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_FM  = new THnSparseD("hpeakd2_ALam_ALam_FM","hpeakd2_ALam_ALam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_FM  = new THnSparseD("hpeakd12_ALam_ALam_FM","hpeakd12_ALam_ALam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_FM_bias  = new THnSparseD("hpeak_ALam_ALam_FM_bias","hpeak_ALam_ALam_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_FM_biasd1  = new THnSparseD("hpeak_ALam_ALam_FM_biasd1","hpeak_ALam_ALam_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_FM_biasd2  = new THnSparseD("hpeak_ALam_ALam_FM_biasd2","hpeak_ALam_ALam_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_FM_biasd12  = new THnSparseD("hpeak_ALam_ALam_FM_biasd12","hpeak_ALam_ALam_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FM_ALam=new TH1D("V0chi2d1_diff_FM_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FM_ALam=new TH1D("V0chi2d2_diff_FM_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FM_ALam=new TH1D("V0chi2d12_diff_FM_ALam","",10000,-0.01,0.11);

THnSparseD *hpeak_ALam_ALam_FU  = new THnSparseD("hpeak_ALam_ALam_FU","hpeak_ALam_ALam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_FU  = new THnSparseD("hpeakd1_ALam_ALam_FU","hpeakd1_ALam_ALam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_FU  = new THnSparseD("hpeakd2_ALam_ALam_FU","hpeakd2_ALam_ALam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_FU  = new THnSparseD("hpeakd12_ALam_ALam_FU","hpeakd12_ALam_ALam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_FU_bias  = new THnSparseD("hpeak_ALam_ALam_FU_bias","hpeak_ALam_ALam_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_FU_biasd1  = new THnSparseD("hpeak_ALam_ALam_FU_biasd1","hpeak_ALam_ALam_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_FU_biasd2  = new THnSparseD("hpeak_ALam_ALam_FU_biasd2","hpeak_ALam_ALam_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_FU_biasd12  = new THnSparseD("hpeak_ALam_ALam_FU_biasd12","hpeak_ALam_ALam_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FU_ALam=new TH1D("V0chi2d1_diff_FU_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FU_ALam=new TH1D("V0chi2d2_diff_FU_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FU_ALam=new TH1D("V0chi2d12_diff_FU_ALam","",10000,-0.01,0.11);

//true matched + true unmatched

THnSparseD *hpeak_ALam_ALam_TM_TU  = new THnSparseD("hpeak_ALam_ALam_TM_TU","hpeak_ALam_ALam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_TM_TU  = new THnSparseD("hpeakd1_ALam_ALam_TM_TU","hpeakd1_ALam_ALam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_TM_TU  = new THnSparseD("hpeakd2_ALam_ALam_TM_TU","hpeakd2_ALam_ALam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_TM_TU  = new THnSparseD("hpeakd12_ALam_ALam_TM_TU","hpeakd12_ALam_ALam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_TM_TU_bias  = new THnSparseD("hpeak_ALam_ALam_TM_TU_bias","hpeak_ALam_ALam_TM_TU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TM_TU_biasd1  = new THnSparseD("hpeak_ALam_ALam_TM_TU_biasd1","hpeak_ALam_ALam_TM_TU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TM_TU_biasd2  = new THnSparseD("hpeak_ALam_ALam_TM_TU_biasd2","hpeak_ALam_ALam_TM_TU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TM_TU_biasd12  = new THnSparseD("hpeak_ALam_ALam_TM_TU_biasd12","hpeak_ALam_ALam_TM_TU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_TU_ALam=new TH1D("V0chi2d1_diff_TM_TU_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_TU_ALam=new TH1D("V0chi2d2_diff_TM_TU_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_TU_ALam=new TH1D("V0chi2d12_diff_TM_TU_ALam","",10000,-0.01,0.11);

//fake matched + fake unmatched

THnSparseD *hpeak_ALam_ALam_FM_FU  = new THnSparseD("hpeak_ALam_ALam_FM_FU","hpeak_ALam_ALam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_FM_FU  = new THnSparseD("hpeakd1_ALam_ALam_FM_FU","hpeakd1_ALam_ALam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_FM_FU  = new THnSparseD("hpeakd2_ALam_ALam_FM_FU","hpeakd2_ALam_ALam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_FM_FU  = new THnSparseD("hpeakd12_ALam_ALam_FM_FU","hpeakd12_ALam_ALam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_FM_FU_bias  = new THnSparseD("hpeak_ALam_ALam_FM_FU_bias","hpeak_ALam_ALam_FM_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_FM_FU_biasd1  = new THnSparseD("hpeak_ALam_ALam_FM_FU_biasd1","hpeak_ALam_ALam_FM_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_FM_FU_biasd2  = new THnSparseD("hpeak_ALam_ALam_FM_FU_biasd2","hpeak_ALam_ALam_FM_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_FM_FU_biasd12  = new THnSparseD("hpeak_ALam_ALam_FM_FU_biasd12","hpeak_ALam_ALam_FM_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FM_FU_ALam=new TH1D("V0chi2d1_diff_FM_FU_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FM_FU_ALam=new TH1D("V0chi2d2_diff_FM_FU_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FM_FU_ALam=new TH1D("V0chi2d12_diff_FM_FU_ALam","",10000,-0.01,0.11);


THnSparseD *hpeak_ALam_ALam_TM_FM  = new THnSparseD("hpeak_ALam_ALam_TM_FM","hpeak_ALam_ALam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_TM_FM  = new THnSparseD("hpeakd1_ALam_ALam_TM_FM","hpeakd1_ALam_ALam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_TM_FM  = new THnSparseD("hpeakd2_ALam_ALam_TM_FM","hpeakd2_ALam_ALam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_TM_FM  = new THnSparseD("hpeakd12_ALam_ALam_TM_FM","hpeakd12_ALam_ALam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_TM_FM_bias  = new THnSparseD("hpeak_ALam_ALam_TM_FM_bias","hpeak_ALam_ALam_TM_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TM_FM_biasd1  = new THnSparseD("hpeak_ALam_ALam_TM_FM_biasd1","hpeak_ALam_ALam_TM_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TM_FM_biasd2  = new THnSparseD("hpeak_ALam_ALam_TM_FM_biasd2","hpeak_ALam_ALam_TM_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TM_FM_biasd12  = new THnSparseD("hpeak_ALam_ALam_TM_FM_biasd12","hpeak_ALam_ALam_TM_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_FM_ALam=new TH1D("V0chi2d1_diff_TM_FM_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_FM_ALam=new TH1D("V0chi2d2_diff_TM_FM_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_FM_ALam=new TH1D("V0chi2d12_diff_TM_FM_ALam","",10000,-0.01,0.11);

THnSparseD *hpeak_ALam_ALam_TM_FU  = new THnSparseD("hpeak_ALam_ALam_TM_FU","hpeak_ALam_ALam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_TM_FU  = new THnSparseD("hpeakd1_ALam_ALam_TM_FU","hpeakd1_ALam_ALam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_TM_FU  = new THnSparseD("hpeakd2_ALam_ALam_TM_FU","hpeakd2_ALam_ALam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_TM_FU  = new THnSparseD("hpeakd12_ALam_ALam_TM_FU","hpeakd12_ALam_ALam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_TM_FU_bias  = new THnSparseD("hpeak_ALam_ALam_TM_FU_bias","hpeak_ALam_ALam_TM_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TM_FU_biasd1  = new THnSparseD("hpeak_ALam_ALam_TM_FU_biasd1","hpeak_ALam_ALam_TM_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TM_FU_biasd2  = new THnSparseD("hpeak_ALam_ALam_TM_FU_biasd2","hpeak_ALam_ALam_TM_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TM_FU_biasd12  = new THnSparseD("hpeak_ALam_ALam_TM_FU_biasd12","hpeak_ALam_ALam_TM_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_FU_ALam=new TH1D("V0chi2d1_diff_TM_FU_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_FU_ALam=new TH1D("V0chi2d2_diff_TM_FU_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_FU_ALam=new TH1D("V0chi2d12_diff_TM_FU_ALam","",10000,-0.01,0.11);


THnSparseD *hpeak_ALam_ALam_TU_FU  = new THnSparseD("hpeak_ALam_ALam_TU_FU","hpeak_ALam_ALam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_TU_FU  = new THnSparseD("hpeakd1_ALam_ALam_TU_FU","hpeakd1_ALam_ALam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_TU_FU  = new THnSparseD("hpeakd2_ALam_ALam_TU_FU","hpeakd2_ALam_ALam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_TU_FU  = new THnSparseD("hpeakd12_ALam_ALam_TU_FU","hpeakd12_ALam_ALam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_TU_FU_bias  = new THnSparseD("hpeak_ALam_ALam_TU_FU_bias","hpeak_ALam_ALam_TU_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TU_FU_biasd1  = new THnSparseD("hpeak_ALam_ALam_TU_FU_biasd1","hpeak_ALam_ALam_TU_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TU_FU_biasd2  = new THnSparseD("hpeak_ALam_ALam_TU_FU_biasd2","hpeak_ALam_ALam_TU_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TU_FU_biasd12  = new THnSparseD("hpeak_ALam_ALam_TU_FU_biasd12","hpeak_ALam_ALam_TU_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_FU_ALam=new TH1D("V0chi2d1_diff_TU_FU_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_FU_ALam=new TH1D("V0chi2d2_diff_TU_FU_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_FU_ALam=new TH1D("V0chi2d12_diff_TU_FU_ALam","",10000,-0.01,0.11);

THnSparseD *hpeak_ALam_ALam_TU_FM  = new THnSparseD("hpeak_ALam_ALam_TU_FM","hpeak_ALam_ALam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_ALam_ALam_TU_FM  = new THnSparseD("hpeakd1_ALam_ALam_TU_FM","hpeakd1_ALam_ALam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_ALam_ALam_TU_FM  = new THnSparseD("hpeakd2_ALam_ALam_TU_FM","hpeakd2_ALam_ALam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_ALam_ALam_TU_FM  = new THnSparseD("hpeakd12_ALam_ALam_TU_FM","hpeakd12_ALam_ALam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_TU_FM_bias  = new THnSparseD("hpeak_ALam_ALam_TU_FM_bias","hpeak_ALam_ALam_TU_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TU_FM_biasd1  = new THnSparseD("hpeak_ALam_ALam_TU_FM_biasd1","hpeak_ALam_ALam_TU_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TU_FM_biasd2  = new THnSparseD("hpeak_ALam_ALam_TU_FM_biasd2","hpeak_ALam_ALam_TU_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_ALam_ALam_TU_FM_biasd12  = new THnSparseD("hpeak_ALam_ALam_TU_FM_biasd12","hpeak_ALam_ALam_TU_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_FM_ALam=new TH1D("V0chi2d1_diff_TU_FM_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_FM_ALam=new TH1D("V0chi2d2_diff_TU_FM_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_FM_ALam=new TH1D("V0chi2d12_diff_TU_FM_ALam","",10000,-0.01,0.11);

//cross 
//LAL

  THnSparseD *hMass_Lam_T  = new THnSparseD("hMass_Lam_T","hMass_Lam_T",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
  THnSparseD *hMass_ALam_T  = new THnSparseD("hMass_ALam_T","hMass_ALam_T",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
  THnSparseD *hMass_Lam_ALam_T  = new THnSparseD("hMass_Lam_ALam_T","hMass_Lam_ALam_T",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
  
  THnSparseD *hpeak_Lam_ALam_T  = new THnSparseD("hpeak_Lam_ALam_T","hpeak_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakcos_Lam_ALam_T  = new THnSparseD("hpeakcos_Lam_ALam_T","hpeakcos_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeak_rot_Lam_ALam_T  = new THnSparseD("hpeak_rot_Lam_ALam_T","hpeak_rot_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeak_inv_Lam_ALam_T  = new THnSparseD("hpeak_inv_Lam_ALam_T","hpeak_inv_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_Lam_ALam_T  = new THnSparseD("hside_Lam_ALam_T","hside_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideL_Lam_ALam_T  = new THnSparseD("hsideL_Lam_ALam_T","hsideL_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideR_Lam_ALam_T  = new THnSparseD("hsideR_Lam_ALam_T","hsideR_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_rot_Lam_ALam_T  = new THnSparseD("hside_rot_Lam_ALam_T","hside_rot_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_inv_Lam_ALam_T  = new THnSparseD("hside_inv_Lam_ALam_T","hside_inv_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_Lam_ALam_T  = new THnSparseD("hpeakside_Lam_ALam_T","hpeakside_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideL_Lam_ALam_T  = new THnSparseD("hpeaksideL_Lam_ALam_T","hpeaksideL_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideR_Lam_ALam_T  = new THnSparseD("hpeaksideR_Lam_ALam_T","hpeaksideR_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_rot_Lam_ALam_T  = new THnSparseD("hpeakside_rot_Lam_ALam_T","hpeakside_rot_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_inv_Lam_ALam_T  = new THnSparseD("hpeakside_inv_Lam_ALam_T","hpeakside_inv_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd1_Lam_ALam_T  = new THnSparseD("hpeakd1_Lam_ALam_T","hpeakd1_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd2_Lam_ALam_T  = new THnSparseD("hpeakd2_Lam_ALam_T","hpeakd2_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd12_Lam_ALam_T  = new THnSparseD("hpeakd12_Lam_ALam_T","hpeakd12_Lam_ALam_T",3,bins3D,xmin3D,xmax3D);

  THnSparseD *hpeak_Lam_ALam_T_mix  = new THnSparseD("hpeak_Lam_ALam_T_mix","hpeak_Lam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_Lam_ALam_T_mix  = new THnSparseD("hside_Lam_ALam_T_mix","hside_Lam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideL_Lam_ALam_T_mix  = new THnSparseD("hsideL_Lam_ALam_T_mix","hsideL_Lam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideR_Lam_ALam_T_mix  = new THnSparseD("hsideR_Lam_ALam_mix","hsideR_Lam_ALam_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_Lam_ALam_T_mix  = new THnSparseD("hpeakside_Lam_ALam_T_mix","hpeakside_Lam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideL_Lam_ALam_T_mix  = new THnSparseD("hpeaksideL_Lam_ALam_T_mix","hpeaksideL_Lam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideR_Lam_ALam_T_mix  = new THnSparseD("hpeaksideR_Lam_ALam_T_mix","hpeaksideR_Lam_ALam_T_mix",3,bins3D,xmin3D,xmax3D);

  THnSparseD *hpeak_Lam_ALam_T_etamix  = new THnSparseD("hpeak_Lam_ALam_T_etamix","hpeak_Lam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_Lam_ALam_T_etamix  = new THnSparseD("hside_Lam_ALam_T_etamix","hside_Lam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideL_Lam_ALam_T_etamix  = new THnSparseD("hsideL_Lam_ALam_T_etamix","hsideL_Lam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideR_Lam_ALam_T_etamix  = new THnSparseD("hsideR_Lam_ALam_T_etamix","hsideR_Lam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_Lam_ALam_T_etamix  = new THnSparseD("hpeakside_Lam_ALam_T_etamix","hpeakside_Lam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideL_Lam_ALam_T_etamix  = new THnSparseD("hpeaksideL_Lam_ALam_T_etamix","hpeaksideL_Lam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideR_Lam_ALam_T_etamix  = new THnSparseD("hpeaksideR_Lam_ALam_T_etamix","hpeaksideR_Lam_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);  
  
  
  TH1D *V0chi2d1_diff_T_Lam_ALam=new TH1D("V0chi2d1_diff_T_Lam_ALam","",10000,-0.01,0.11);
  TH1D *V0chi2d2_diff_T_Lam_ALam=new TH1D("V0chi2d2_diff_T_Lam_ALam","",10000,-0.01,0.11);
  TH1D *V0chi2d12_diff_T_Lam_ALam=new TH1D("V0chi2d12_diff_T_Lam_ALam","",10000,-0.01,0.11);
  
//------------------------------------------------------------------------
//fake
//------------------------------------------------------------------------ 

  THnSparseD *hMass_Lam_F  = new THnSparseD("hMass_Lam_F","hMass_Lam_F",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
  THnSparseD *hMass_ALam_F  = new THnSparseD("hMass_ALam_F","hMass_ALam_F",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
  THnSparseD *hMass_Lam_ALam_F  = new THnSparseD("hMass_Lam_ALam_F","hMass_Lam_ALam_F",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
  
  THnSparseD *hpeak_Lam_ALam_F  = new THnSparseD("hpeak_Lam_ALam_F","hpeak_Lam_ALam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd1_Lam_ALam_F  = new THnSparseD("hpeakd1_Lam_ALam_F","hpeakd1_Lam_ALam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd2_Lam_ALam_F  = new THnSparseD("hpeakd2_Lam_ALam_F","hpeakd2_Lam_ALam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd12_Lam_ALam_F  = new THnSparseD("hpeakd12_Lam_ALam_F","hpeakd12_Lam_ALam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_Lam_ALam_F  = new THnSparseD("hside_Lam_ALam_F","hside_Lam_ALam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_Lam_ALam_F  = new THnSparseD("hpeakside_Lam_ALam_F","hpeakside_Lam_ALam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeak_Lam_ALam_F_bias  = new THnSparseD("hpeak_Lam_ALam_F_bias","hpeak_Lam_ALam_F_bias",3,bins3Db,xmin3Db,xmax3Db);
  THnSparseD *hpeak_Lam_ALam_F_biasd1  = new THnSparseD("hpeak_Lam_ALam_F_biasd1","hpeak_Lam_ALam_F_biasd1",3,bins3Db,xmin3Db,xmax3Db);
  THnSparseD *hpeak_Lam_ALam_F_biasd2  = new THnSparseD("hpeak_Lam_ALam_F_biasd2","hpeak_Lam_ALam_F_biasd2",3,bins3Db,xmin3Db,xmax3Db);
  THnSparseD *hpeak_Lam_ALam_F_biasd12  = new THnSparseD("hpeak_Lam_ALam_F_biasd12","hpeak_Lam_ALam_F_biasd12",3,bins3Db,xmin3Db,xmax3Db);

  TH1D *V0chi2d1_diff_F_Lam_ALam=new TH1D("V0chi2d1_diff_F_Lam_ALam","",10000,-0.01,0.11);
  TH1D *V0chi2d2_diff_F_Lam_ALam=new TH1D("V0chi2d2_diff_F_Lam_ALam","",10000,-0.01,0.11);
  TH1D *V0chi2d12_diff_F_Lam_ALam=new TH1D("V0chi2d12_diff_F_Lam_ALam","",10000,-0.01,0.11);
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
THnSparseD *hpeak_Lam_ALam_TF  = new THnSparseD("hpeak_Lam_ALam_TF","hpeak_Lam_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_ALam_TF  = new THnSparseD("hpeakd1_Lam_ALam_TF","hpeakd1_Lam_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_ALam_TF  = new THnSparseD("hpeakd2_Lam_ALam_TF","hpeakd2_Lam_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_ALam_TF  = new THnSparseD("hpeakd12_Lam_ALam_TF","hpeakd12_Lam_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_Lam_ALam_TF  = new THnSparseD("hside_Lam_ALam_TF","hside_Lam_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_Lam_ALam_TF  = new THnSparseD("hpeakside_Lam_ALam_TF","hpeakside_Lam_ALam_TF",3,bins3D,xmin3D,xmax3D);  
THnSparseD *hpeak_Lam_ALam_TF_bias  = new THnSparseD("hpeak_Lam_ALam_TF_bias","hpeak_Lam_ALam_TF_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TF_biasd1  = new THnSparseD("hpeak_Lam_ALam_TF_biasd1","hpeak_Lam_ALam_TF_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TF_biasd2  = new THnSparseD("hpeak_Lam_ALam_TF_biasd2","hpeak_Lam_ALam_TF_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TF_biasd12  = new THnSparseD("hpeak_Lam_ALam_TF_biasd12","hpeak_Lam_ALam_TF_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TF_Lam_ALam=new TH1D("V0chi2d1_diff_TF_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TF_Lam_ALam=new TH1D("V0chi2d2_diff_TF_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TF_Lam_ALam=new TH1D("V0chi2d12_diff_TF_Lam_ALam","",10000,-0.01,0.11);

  
/////// matching hist ////// 

THnSparseD *hMass_Lam_ALam_TM  = new THnSparseD("hMass_Lam_ALam_TM","hMass_Lam_ALam_TM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_Lam_TM  = new THnSparseD("hMass_Lam_TM","hMass_Lam_TM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_ALam_TM  = new THnSparseD("hMass_ALam_TM","hMass_ALam_TM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   

THnSparseD *hMass_Lam_ALam_TU  = new THnSparseD("hMass_Lam_ALam_TU","hMass_Lam_ALam_TU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_Lam_TU  = new THnSparseD("hMass_Lam_TU","hMass_Lam_TU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_ALam_TU  = new THnSparseD("hMass_ALam_TU","hMass_ALam_TU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   

THnSparseD *hMass_Lam_ALam_FM  = new THnSparseD("hMass_Lam_ALam_FM","hMass_Lam_ALam_FM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_Lam_FM  = new THnSparseD("hMass_Lam_FM","hMass_Lam_FM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_ALam_FM  = new THnSparseD("hMass_ALam_FM","hMass_ALam_FM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   

THnSparseD *hMass_Lam_ALam_FU  = new THnSparseD("hMass_Lam_ALam_FU","hMass_Lam_ALam_FU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_Lam_FU  = new THnSparseD("hMass_Lam_FU","hMass_Lam_FU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_ALam_FU  = new THnSparseD("hMass_ALam_FU","hMass_ALam_FU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   


THnSparseD *hpeak_Lam_ALam_TM  = new THnSparseD("hpeak_Lam_ALam_TM","hpeak_Lam_ALam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_ALam_TM  = new THnSparseD("hpeakd1_Lam_ALam_TM","hpeakd1_Lam_ALam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_ALam_TM  = new THnSparseD("hpeakd2_Lam_ALam_TM","hpeakd2_Lam_ALam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_ALam_TM  = new THnSparseD("hpeakd12_Lam_ALam_TM","hpeakd12_Lam_ALam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_ALam_TM_bias  = new THnSparseD("hpeak_Lam_ALam_TM_bias","hpeak_Lam_ALam_TM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TM_biasd1  = new THnSparseD("hpeak_Lam_ALam_TM_biasd1","hpeak_Lam_ALam_TM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TM_biasd2  = new THnSparseD("hpeak_Lam_ALam_TM_biasd2","hpeak_Lam_ALam_TM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TM_biasd12  = new THnSparseD("hpeak_Lam_ALam_TM_biasd12","hpeak_Lam_ALam_TM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_Lam_ALam=new TH1D("V0chi2d1_diff_TM_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_Lam_ALam=new TH1D("V0chi2d2_diff_TM_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_Lam_ALam=new TH1D("V0chi2d12_diff_TM_Lam_ALam","",10000,-0.01,0.11);

THnSparseD *hpeak_Lam_ALam_TU  = new THnSparseD("hpeak_Lam_ALam_TU","hpeak_Lam_ALam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_ALam_TU  = new THnSparseD("hpeakd1_Lam_ALam_TU","hpeakd1_Lam_ALam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_ALam_TU  = new THnSparseD("hpeakd2_Lam_ALam_TU","hpeakd2_Lam_ALam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_ALam_TU  = new THnSparseD("hpeakd12_Lam_ALam_TU","hpeakd12_Lam_ALam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_ALam_TU_bias  = new THnSparseD("hpeak_Lam_ALam_TU_bias","hpeak_Lam_ALam_TU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TU_biasd1  = new THnSparseD("hpeak_Lam_ALam_TU_biasd1","hpeak_Lam_ALam_TU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TU_biasd2  = new THnSparseD("hpeak_Lam_ALam_TU_biasd2","hpeak_Lam_ALam_TU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TU_biasd12  = new THnSparseD("hpeak_Lam_ALam_TU_biasd12","hpeak_Lam_ALam_TU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_Lam_ALam=new TH1D("V0chi2d1_diff_TU_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_Lam_ALam=new TH1D("V0chi2d2_diff_TU_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_Lam_ALam=new TH1D("V0chi2d12_diff_TU_Lam_ALam","",10000,-0.01,0.11);

THnSparseD *hpeak_Lam_ALam_FM  = new THnSparseD("hpeak_Lam_ALam_FM","hpeak_Lam_ALam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_ALam_FM  = new THnSparseD("hpeakd1_Lam_ALam_FM","hpeakd1_Lam_ALam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_ALam_FM  = new THnSparseD("hpeakd2_Lam_ALam_FM","hpeakd2_Lam_ALam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_ALam_FM  = new THnSparseD("hpeakd12_Lam_ALam_FM","hpeakd12_Lam_ALam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_ALam_FM_bias  = new THnSparseD("hpeak_Lam_ALam_FM_bias","hpeak_Lam_ALam_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_FM_biasd1  = new THnSparseD("hpeak_Lam_ALam_FM_biasd1","hpeak_Lam_ALam_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_FM_biasd2  = new THnSparseD("hpeak_Lam_ALam_FM_biasd2","hpeak_Lam_ALam_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_FM_biasd12  = new THnSparseD("hpeak_Lam_ALam_FM_biasd12","hpeak_Lam_ALam_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FM_Lam_ALam=new TH1D("V0chi2d1_diff_FM_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FM_Lam_ALam=new TH1D("V0chi2d2_diff_FM_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FM_Lam_ALam=new TH1D("V0chi2d12_diff_FM_Lam_ALam","",10000,-0.01,0.11);

THnSparseD *hpeak_Lam_ALam_FU  = new THnSparseD("hpeak_Lam_ALam_FU","hpeak_Lam_ALam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_ALam_FU  = new THnSparseD("hpeakd1_Lam_ALam_FU","hpeakd1_Lam_ALam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_ALam_FU  = new THnSparseD("hpeakd2_Lam_ALam_FU","hpeakd2_Lam_ALam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_ALam_FU  = new THnSparseD("hpeakd12_Lam_ALam_FU","hpeakd12_Lam_ALam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_ALam_FU_bias  = new THnSparseD("hpeak_Lam_ALam_FU_bias","hpeak_Lam_ALam_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_FU_biasd1  = new THnSparseD("hpeak_Lam_ALam_FU_biasd1","hpeak_Lam_ALam_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_FU_biasd2  = new THnSparseD("hpeak_Lam_ALam_FU_biasd2","hpeak_Lam_ALam_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_FU_biasd12  = new THnSparseD("hpeak_Lam_ALam_FU_biasd12","hpeak_Lam_ALam_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FU_Lam_ALam=new TH1D("V0chi2d1_diff_FU_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FU_Lam_ALam=new TH1D("V0chi2d2_diff_FU_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FU_Lam_ALam=new TH1D("V0chi2d12_diff_FU_Lam_ALam","",10000,-0.01,0.11);

//true matched + true unmatched

THnSparseD *hpeak_Lam_ALam_TM_TU  = new THnSparseD("hpeak_Lam_ALam_TM_TU","hpeak_Lam_ALam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_ALam_TM_TU  = new THnSparseD("hpeakd1_Lam_ALam_TM_TU","hpeakd1_Lam_ALam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_ALam_TM_TU  = new THnSparseD("hpeakd2_Lam_ALam_TM_TU","hpeakd2_Lam_ALam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_ALam_TM_TU  = new THnSparseD("hpeakd12_Lam_ALam_TM_TU","hpeakd12_Lam_ALam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_ALam_TM_TU_bias  = new THnSparseD("hpeak_Lam_ALam_TM_TU_bias","hpeak_Lam_ALam_TM_TU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TM_TU_biasd1  = new THnSparseD("hpeak_Lam_ALam_TM_TU_biasd1","hpeak_Lam_ALam_TM_TU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TM_TU_biasd2  = new THnSparseD("hpeak_Lam_ALam_TM_TU_biasd2","hpeak_Lam_ALam_TM_TU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TM_TU_biasd12  = new THnSparseD("hpeak_Lam_ALam_TM_TU_biasd12","hpeak_Lam_ALam_TM_TU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_TU_Lam_ALam=new TH1D("V0chi2d1_diff_TM_TU_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_TU_Lam_ALam=new TH1D("V0chi2d2_diff_TM_TU_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_TU_Lam_ALam=new TH1D("V0chi2d12_diff_TM_TU_Lam_ALam","",10000,-0.01,0.11);

//fake matched + fake unmatched

THnSparseD *hpeak_Lam_ALam_FM_FU  = new THnSparseD("hpeak_Lam_ALam_FM_FU","hpeak_Lam_ALam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_ALam_FM_FU  = new THnSparseD("hpeakd1_Lam_ALam_FM_FU","hpeakd1_Lam_ALam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_ALam_FM_FU  = new THnSparseD("hpeakd2_Lam_ALam_FM_FU","hpeakd2_Lam_ALam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_ALam_FM_FU  = new THnSparseD("hpeakd12_Lam_ALam_FM_FU","hpeakd12_Lam_ALam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_ALam_FM_FU_bias  = new THnSparseD("hpeak_Lam_ALam_FM_FU_bias","hpeak_Lam_ALam_FM_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_FM_FU_biasd1  = new THnSparseD("hpeak_Lam_ALam_FM_FU_biasd1","hpeak_Lam_ALam_FM_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_FM_FU_biasd2  = new THnSparseD("hpeak_Lam_ALam_FM_FU_biasd2","hpeak_Lam_ALam_FM_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_FM_FU_biasd12  = new THnSparseD("hpeak_Lam_ALam_FM_FU_biasd12","hpeak_Lam_ALam_FM_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FM_FU_Lam_ALam=new TH1D("V0chi2d1_diff_FM_FU_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FM_FU_Lam_ALam=new TH1D("V0chi2d2_diff_FM_FU_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FM_FU_Lam_ALam=new TH1D("V0chi2d12_diff_FM_FU_Lam_ALam","",10000,-0.01,0.11);

//true matched + fake matched

THnSparseD *hpeak_Lam_ALam_TM_FM  = new THnSparseD("hpeak_Lam_ALam_TM_FM","hpeak_Lam_ALam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_ALam_TM_FM  = new THnSparseD("hpeakd1_Lam_ALam_TM_FM","hpeakd1_Lam_ALam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_ALam_TM_FM  = new THnSparseD("hpeakd2_Lam_ALam_TM_FM","hpeakd2_Lam_ALam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_ALam_TM_FM  = new THnSparseD("hpeakd12_Lam_ALam_TM_FM","hpeakd12_Lam_ALam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_ALam_TM_FM_bias  = new THnSparseD("hpeak_Lam_ALam_TM_FM_bias","hpeak_Lam_ALam_TM_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TM_FM_biasd1  = new THnSparseD("hpeak_Lam_ALam_TM_FM_biasd1","hpeak_Lam_ALam_TM_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TM_FM_biasd2  = new THnSparseD("hpeak_Lam_ALam_TM_FM_biasd2","hpeak_Lam_ALam_TM_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TM_FM_biasd12  = new THnSparseD("hpeak_Lam_ALam_TM_FM_biasd12","hpeak_Lam_ALam_TM_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_FM_Lam_ALam=new TH1D("V0chi2d1_diff_TM_FM_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_FM_Lam_ALam=new TH1D("V0chi2d2_diff_TM_FM_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_FM_Lam_ALam=new TH1D("V0chi2d12_diff_TM_FM_Lam_ALam","",10000,-0.01,0.11);

//true matched + fake unmatched

THnSparseD *hpeak_Lam_ALam_TM_FU  = new THnSparseD("hpeak_Lam_ALam_TM_FU","hpeak_Lam_ALam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_ALam_TM_FU  = new THnSparseD("hpeakd1_Lam_ALam_TM_FU","hpeakd1_Lam_ALam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_ALam_TM_FU  = new THnSparseD("hpeakd2_Lam_ALam_TM_FU","hpeakd2_Lam_ALam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_ALam_TM_FU  = new THnSparseD("hpeakd12_Lam_ALam_TM_FU","hpeakd12_Lam_ALam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_ALam_TM_FU_bias  = new THnSparseD("hpeak_Lam_ALam_TM_FU_bias","hpeak_Lam_ALam_TM_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TM_FU_biasd1  = new THnSparseD("hpeak_Lam_ALam_TM_FU_biasd1","hpeak_Lam_ALam_TM_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TM_FU_biasd2  = new THnSparseD("hpeak_Lam_ALam_TM_FU_biasd2","hpeak_Lam_ALam_TM_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TM_FU_biasd12  = new THnSparseD("hpeak_Lam_ALam_TM_FU_biasd12","hpeak_Lam_ALam_TM_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_FU_Lam_ALam=new TH1D("V0chi2d1_diff_TM_FU_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_FU_Lam_ALam=new TH1D("V0chi2d2_diff_TM_FU_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_FU_Lam_ALam=new TH1D("V0chi2d12_diff_TM_FU_Lam_ALam","",10000,-0.01,0.11);

//true unmatched + fake unmatched

THnSparseD *hpeak_Lam_ALam_TU_FU  = new THnSparseD("hpeak_Lam_ALam_TU_FU","hpeak_Lam_ALam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_ALam_TU_FU  = new THnSparseD("hpeakd1_Lam_ALam_TU_FU","hpeakd1_Lam_ALam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_ALam_TU_FU  = new THnSparseD("hpeakd2_Lam_ALam_TU_FU","hpeakd2_Lam_ALam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_ALam_TU_FU  = new THnSparseD("hpeakd12_Lam_ALam_TU_FU","hpeakd12_Lam_ALam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_ALam_TU_FU_bias  = new THnSparseD("hpeak_Lam_ALam_TU_FU_bias","hpeak_Lam_ALam_TU_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TU_FU_biasd1  = new THnSparseD("hpeak_Lam_ALam_TU_FU_biasd1","hpeak_Lam_ALam_TU_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TU_FU_biasd2  = new THnSparseD("hpeak_Lam_ALam_TU_FU_biasd2","hpeak_Lam_ALam_TU_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TU_FU_biasd12  = new THnSparseD("hpeak_Lam_ALam_TU_FU_biasd12","hpeak_Lam_ALam_TU_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_FU_Lam_ALam=new TH1D("V0chi2d1_diff_TU_FU_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_FU_Lam_ALam=new TH1D("V0chi2d2_diff_TU_FU_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_FU_Lam_ALam=new TH1D("V0chi2d12_diff_TU_FU_Lam_ALam","",10000,-0.01,0.11);

//true unmatched + fake matched

THnSparseD *hpeak_Lam_ALam_TU_FM  = new THnSparseD("hpeak_Lam_ALam_TU_FM","hpeak_Lam_ALam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_Lam_ALam_TU_FM  = new THnSparseD("hpeakd1_Lam_ALam_TU_FM","hpeakd1_Lam_ALam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_Lam_ALam_TU_FM  = new THnSparseD("hpeakd2_Lam_ALam_TU_FM","hpeakd2_Lam_ALam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_Lam_ALam_TU_FM  = new THnSparseD("hpeakd12_Lam_ALam_TU_FM","hpeakd12_Lam_ALam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_ALam_TU_FM_bias  = new THnSparseD("hpeak_Lam_ALam_TU_FM_bias","hpeak_Lam_ALam_TU_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TU_FM_biasd1  = new THnSparseD("hpeak_Lam_ALam_TU_FM_biasd1","hpeak_Lam_ALam_TU_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TU_FM_biasd2  = new THnSparseD("hpeak_Lam_ALam_TU_FM_biasd2","hpeak_Lam_ALam_TU_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_Lam_ALam_TU_FM_biasd12  = new THnSparseD("hpeak_Lam_ALam_TU_FM_biasd12","hpeak_Lam_ALam_TU_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_FM_Lam_ALam=new TH1D("V0chi2d1_diff_TU_FM_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_FM_Lam_ALam=new TH1D("V0chi2d2_diff_TU_FM_Lam_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_FM_Lam_ALam=new TH1D("V0chi2d12_diff_TU_FM_Lam_ALam","",10000,-0.01,0.11);


//K0sLam

  THnSparseD *hMass_K0sL_T  = new THnSparseD("hMass_K0sL_T","hMass_K0sL_T",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
  THnSparseD *hMass_LamK_T  = new THnSparseD("hMass_LamK_T","hMass_LamK_T",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
  THnSparseD *hMass_K0s_Lam_T  = new THnSparseD("hMass_K0s_Lam_T","hMass_K0s_Lam_T",4,bins3DM_KL,xmin3DM_KL,xmax3DM_KL);   
  
  THnSparseD *hpeak_K0s_Lam_T  = new THnSparseD("hpeak_K0s_Lam_T","hpeak_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakcos_K0s_Lam_T  = new THnSparseD("hpeakcos_K0s_Lam_T","hpeakcos_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeak_rot_K0s_Lam_T  = new THnSparseD("hpeak_rot_K0s_Lam_T","hpeak_rot_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeak_inv_K0s_Lam_T  = new THnSparseD("hpeak_inv_K0s_Lam_T","hpeak_inv_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_K0s_Lam_T  = new THnSparseD("hside_K0s_Lam_T","hside_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideL_K0s_Lam_T  = new THnSparseD("hsideL_K0s_Lam_T","hsideL_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideR_K0s_Lam_T  = new THnSparseD("hsideR_K0s_Lam_T","hsideR_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_rot_K0s_Lam_T  = new THnSparseD("hside_rot_K0s_Lam_T","hside_rot_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_inv_K0s_Lam_T  = new THnSparseD("hside_inv_K0s_Lam_T","hside_inv_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_K0s_Lam_T  = new THnSparseD("hpeakside_K0s_Lam_T","hpeakside_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideL_K0s_Lam_T  = new THnSparseD("hpeaksideL_K0s_Lam_T","hpeaksideL_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideR_K0s_Lam_T  = new THnSparseD("hpeaksideR_K0s_Lam_T","hpeaksideR_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_rot_K0s_Lam_T  = new THnSparseD("hpeakside_rot_K0s_Lam_T","hpeakside_rot_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_inv_K0s_Lam_T  = new THnSparseD("hpeakside_inv_K0s_Lam_T","hpeakside_inv_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd1_K0s_Lam_T  = new THnSparseD("hpeakd1_K0s_Lam_T","hpeakd1_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd2_K0s_Lam_T  = new THnSparseD("hpeakd2_K0s_Lam_T","hpeakd2_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd12_K0s_Lam_T  = new THnSparseD("hpeakd12_K0s_Lam_T","hpeakd12_K0s_Lam_T",3,bins3D,xmin3D,xmax3D);

  THnSparseD *hpeak_K0s_Lam_T_mix  = new THnSparseD("hpeak_K0s_Lam_T_mix","hpeak_K0s_Lam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_K0s_Lam_T_mix  = new THnSparseD("hside_K0s_Lam_T_mix","hside_K0s_Lam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideL_K0s_Lam_T_mix  = new THnSparseD("hsideL_K0s_Lam_T_mix","hsideL_K0s_Lam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideR_K0s_Lam_T_mix  = new THnSparseD("hsideR_K0s_Lam_mix","hsideR_K0s_Lam_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_K0s_Lam_T_mix  = new THnSparseD("hpeakside_K0s_Lam_T_mix","hpeakside_K0s_Lam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideL_K0s_Lam_T_mix  = new THnSparseD("hpeaksideL_K0s_Lam_T_mix","hpeaksideL_K0s_Lam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideR_K0s_Lam_T_mix  = new THnSparseD("hpeaksideR_K0s_Lam_T_mix","hpeaksideR_K0s_Lam_T_mix",3,bins3D,xmin3D,xmax3D);

  THnSparseD *hpeak_K0s_Lam_T_etamix  = new THnSparseD("hpeak_K0s_Lam_T_etamix","hpeak_K0s_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_K0s_Lam_T_etamix  = new THnSparseD("hside_K0s_Lam_T_etamix","hside_K0s_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideL_K0s_Lam_T_etamix  = new THnSparseD("hsideL_K0s_Lam_T_etamix","hsideL_K0s_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideR_K0s_Lam_T_etamix  = new THnSparseD("hsideR_K0s_Lam_T_etamix","hsideR_K0s_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_K0s_Lam_T_etamix  = new THnSparseD("hpeakside_K0s_Lam_T_etamix","hpeakside_K0s_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideL_K0s_Lam_T_etamix  = new THnSparseD("hpeaksideL_K0s_Lam_T_etamix","hpeaksideL_K0s_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideR_K0s_Lam_T_etamix  = new THnSparseD("hpeaksideR_K0s_Lam_T_etamix","hpeaksideR_K0s_Lam_T_etamix",3,bins3D,xmin3D,xmax3D);  
  
  
  TH1D *V0chi2d1_diff_T_K0s_Lam=new TH1D("V0chi2d1_diff_T_K0s_Lam","",10000,-0.01,0.11);
  TH1D *V0chi2d2_diff_T_K0s_Lam=new TH1D("V0chi2d2_diff_T_K0s_Lam","",10000,-0.01,0.11);
  TH1D *V0chi2d12_diff_T_K0s_Lam=new TH1D("V0chi2d12_diff_T_K0s_Lam","",10000,-0.01,0.11);
  
//------------------------------------------------------------------------
//fake
//------------------------------------------------------------------------ 

  THnSparseD *hMass_K0sL_F  = new THnSparseD("hMass_K0sL_F","hMass_K0sL_F",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
  THnSparseD *hMass_LamK_F  = new THnSparseD("hMass_LamK_F","hMass_LamK_F",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
  THnSparseD *hMass_K0s_Lam_F  = new THnSparseD("hMass_K0s_Lam_F","hMass_K0s_Lam_F",4,bins3DM_KL,xmin3DM_KL,xmax3DM_KL);   
  
  THnSparseD *hpeak_K0s_Lam_F  = new THnSparseD("hpeak_K0s_Lam_F","hpeak_K0s_Lam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd1_K0s_Lam_F  = new THnSparseD("hpeakd1_K0s_Lam_F","hpeakd1_K0s_Lam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd2_K0s_Lam_F  = new THnSparseD("hpeakd2_K0s_Lam_F","hpeakd2_K0s_Lam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd12_K0s_Lam_F  = new THnSparseD("hpeakd12_K0s_Lam_F","hpeakd12_K0s_Lam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_K0s_Lam_F  = new THnSparseD("hside_K0s_Lam_F","hside_K0s_Lam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_K0s_Lam_F  = new THnSparseD("hpeakside_K0s_Lam_F","hpeakside_K0s_Lam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeak_K0s_Lam_F_bias  = new THnSparseD("hpeak_K0s_Lam_F_bias","hpeak_K0s_Lam_F_bias",3,bins3Db,xmin3Db,xmax3Db);
  THnSparseD *hpeak_K0s_Lam_F_biasd1  = new THnSparseD("hpeak_K0s_Lam_F_biasd1","hpeak_K0s_Lam_F_biasd1",3,bins3Db,xmin3Db,xmax3Db);
  THnSparseD *hpeak_K0s_Lam_F_biasd2  = new THnSparseD("hpeak_K0s_Lam_F_biasd2","hpeak_K0s_Lam_F_biasd2",3,bins3Db,xmin3Db,xmax3Db);
  THnSparseD *hpeak_K0s_Lam_F_biasd12  = new THnSparseD("hpeak_K0s_Lam_F_biasd12","hpeak_K0s_Lam_F_biasd12",3,bins3Db,xmin3Db,xmax3Db);

  TH1D *V0chi2d1_diff_F_K0s_Lam=new TH1D("V0chi2d1_diff_F_K0s_Lam","",10000,-0.01,0.11);
  TH1D *V0chi2d2_diff_F_K0s_Lam=new TH1D("V0chi2d2_diff_F_K0s_Lam","",10000,-0.01,0.11);
  TH1D *V0chi2d12_diff_F_K0s_Lam=new TH1D("V0chi2d12_diff_F_K0s_Lam","",10000,-0.01,0.11);
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
THnSparseD *hpeak_K0s_Lam_TF  = new THnSparseD("hpeak_K0s_Lam_TF","hpeak_K0s_Lam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_Lam_TF  = new THnSparseD("hpeakd1_K0s_Lam_TF","hpeakd1_K0s_Lam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_Lam_TF  = new THnSparseD("hpeakd2_K0s_Lam_TF","hpeakd2_K0s_Lam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_Lam_TF  = new THnSparseD("hpeakd12_K0s_Lam_TF","hpeakd12_K0s_Lam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_K0s_Lam_TF  = new THnSparseD("hside_K0s_Lam_TF","hside_K0s_Lam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_K0s_Lam_TF  = new THnSparseD("hpeakside_K0s_Lam_TF","hpeakside_K0s_Lam_TF",3,bins3D,xmin3D,xmax3D);  
THnSparseD *hpeak_K0s_Lam_TF_bias  = new THnSparseD("hpeak_K0s_Lam_TF_bias","hpeak_K0s_Lam_TF_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TF_biasd1  = new THnSparseD("hpeak_K0s_Lam_TF_biasd1","hpeak_K0s_Lam_TF_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TF_biasd2  = new THnSparseD("hpeak_K0s_Lam_TF_biasd2","hpeak_K0s_Lam_TF_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TF_biasd12  = new THnSparseD("hpeak_K0s_Lam_TF_biasd12","hpeak_K0s_Lam_TF_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TF_K0s_Lam=new TH1D("V0chi2d1_diff_TF_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TF_K0s_Lam=new TH1D("V0chi2d2_diff_TF_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TF_K0s_Lam=new TH1D("V0chi2d12_diff_TF_K0s_Lam","",10000,-0.01,0.11);

  
/////// matching hist ////// 

THnSparseD *hMass_K0sL_TM  = new THnSparseD("hMass_K0sL_TM","hMass_K0sL_TM",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_LamK_TM  = new THnSparseD("hMass_LamK_TM","hMass_LamK_TM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_K0s_Lam_TM  = new THnSparseD("hMass_K0s_Lam_TM","hMass_K0s_Lam_TM",4,bins3DM_KL,xmin3DM_KL,xmax3DM_KL);   

THnSparseD *hMass_K0sL_TU  = new THnSparseD("hMass_K0sL_TU","hMass_K0sL_TU",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_LamK_TU  = new THnSparseD("hMass_LamK_TU","hMass_LamK_TU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_K0s_Lam_TU  = new THnSparseD("hMass_K0s_Lam_TU","hMass_K0s_Lam_TU",4,bins3DM_KL,xmin3DM_KL,xmax3DM_KL);   

THnSparseD *hMass_K0sL_FM  = new THnSparseD("hMass_K0sL_FM","hMass_K0sL_FM",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_LamK_FM  = new THnSparseD("hMass_LamK_FM","hMass_LamK_FM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_K0s_Lam_FM  = new THnSparseD("hMass_K0s_Lam_FM","hMass_K0s_Lam_FM",4,bins3DM_KL,xmin3DM_KL,xmax3DM_KL);   

THnSparseD *hMass_K0sL_FU  = new THnSparseD("hMass_K0sL_FU","hMass_K0sL_FU",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_LamK_FU  = new THnSparseD("hMass_LamK_FU","hMass_LamK_FU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_K0s_Lam_FU  = new THnSparseD("hMass_K0s_Lam_FU","hMass_K0s_Lam_FU",4,bins3DM_KL,xmin3DM_KL,xmax3DM_KL);   


THnSparseD *hpeak_K0s_Lam_TM  = new THnSparseD("hpeak_K0s_Lam_TM","hpeak_K0s_Lam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_Lam_TM  = new THnSparseD("hpeakd1_K0s_Lam_TM","hpeakd1_K0s_Lam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_Lam_TM  = new THnSparseD("hpeakd2_K0s_Lam_TM","hpeakd2_K0s_Lam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_Lam_TM  = new THnSparseD("hpeakd12_K0s_Lam_TM","hpeakd12_K0s_Lam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_Lam_TM_bias  = new THnSparseD("hpeak_K0s_Lam_TM_bias","hpeak_K0s_Lam_TM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TM_biasd1  = new THnSparseD("hpeak_K0s_Lam_TM_biasd1","hpeak_K0s_Lam_TM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TM_biasd2  = new THnSparseD("hpeak_K0s_Lam_TM_biasd2","hpeak_K0s_Lam_TM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TM_biasd12  = new THnSparseD("hpeak_K0s_Lam_TM_biasd12","hpeak_K0s_Lam_TM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_K0s_Lam=new TH1D("V0chi2d1_diff_TM_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_K0s_Lam=new TH1D("V0chi2d2_diff_TM_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_K0s_Lam=new TH1D("V0chi2d12_diff_TM_K0s_Lam","",10000,-0.01,0.11);

THnSparseD *hpeak_K0s_Lam_TU  = new THnSparseD("hpeak_K0s_Lam_TU","hpeak_K0s_Lam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_Lam_TU  = new THnSparseD("hpeakd1_K0s_Lam_TU","hpeakd1_K0s_Lam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_Lam_TU  = new THnSparseD("hpeakd2_K0s_Lam_TU","hpeakd2_K0s_Lam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_Lam_TU  = new THnSparseD("hpeakd12_K0s_Lam_TU","hpeakd12_K0s_Lam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_Lam_TU_bias  = new THnSparseD("hpeak_K0s_Lam_TU_bias","hpeak_K0s_Lam_TU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TU_biasd1  = new THnSparseD("hpeak_K0s_Lam_TU_biasd1","hpeak_K0s_Lam_TU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TU_biasd2  = new THnSparseD("hpeak_K0s_Lam_TU_biasd2","hpeak_K0s_Lam_TU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TU_biasd12  = new THnSparseD("hpeak_K0s_Lam_TU_biasd12","hpeak_K0s_Lam_TU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_K0s_Lam=new TH1D("V0chi2d1_diff_TU_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_K0s_Lam=new TH1D("V0chi2d2_diff_TU_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_K0s_Lam=new TH1D("V0chi2d12_diff_TU_K0s_Lam","",10000,-0.01,0.11);

THnSparseD *hpeak_K0s_Lam_FM  = new THnSparseD("hpeak_K0s_Lam_FM","hpeak_K0s_Lam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_Lam_FM  = new THnSparseD("hpeakd1_K0s_Lam_FM","hpeakd1_K0s_Lam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_Lam_FM  = new THnSparseD("hpeakd2_K0s_Lam_FM","hpeakd2_K0s_Lam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_Lam_FM  = new THnSparseD("hpeakd12_K0s_Lam_FM","hpeakd12_K0s_Lam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_Lam_FM_bias  = new THnSparseD("hpeak_K0s_Lam_FM_bias","hpeak_K0s_Lam_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_FM_biasd1  = new THnSparseD("hpeak_K0s_Lam_FM_biasd1","hpeak_K0s_Lam_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_FM_biasd2  = new THnSparseD("hpeak_K0s_Lam_FM_biasd2","hpeak_K0s_Lam_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_FM_biasd12  = new THnSparseD("hpeak_K0s_Lam_FM_biasd12","hpeak_K0s_Lam_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FM_K0s_Lam=new TH1D("V0chi2d1_diff_FM_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FM_K0s_Lam=new TH1D("V0chi2d2_diff_FM_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FM_K0s_Lam=new TH1D("V0chi2d12_diff_FM_K0s_Lam","",10000,-0.01,0.11);

THnSparseD *hpeak_K0s_Lam_FU  = new THnSparseD("hpeak_K0s_Lam_FU","hpeak_K0s_Lam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_Lam_FU  = new THnSparseD("hpeakd1_K0s_Lam_FU","hpeakd1_K0s_Lam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_Lam_FU  = new THnSparseD("hpeakd2_K0s_Lam_FU","hpeakd2_K0s_Lam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_Lam_FU  = new THnSparseD("hpeakd12_K0s_Lam_FU","hpeakd12_K0s_Lam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_Lam_FU_bias  = new THnSparseD("hpeak_K0s_Lam_FU_bias","hpeak_K0s_Lam_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_FU_biasd1  = new THnSparseD("hpeak_K0s_Lam_FU_biasd1","hpeak_K0s_Lam_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_FU_biasd2  = new THnSparseD("hpeak_K0s_Lam_FU_biasd2","hpeak_K0s_Lam_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_FU_biasd12  = new THnSparseD("hpeak_K0s_Lam_FU_biasd12","hpeak_K0s_Lam_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FU_K0s_Lam=new TH1D("V0chi2d1_diff_FU_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FU_K0s_Lam=new TH1D("V0chi2d2_diff_FU_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FU_K0s_Lam=new TH1D("V0chi2d12_diff_FU_K0s_Lam","",10000,-0.01,0.11);

//true matched + true unmatched

THnSparseD *hpeak_K0s_Lam_TM_TU  = new THnSparseD("hpeak_K0s_Lam_TM_TU","hpeak_K0s_Lam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_Lam_TM_TU  = new THnSparseD("hpeakd1_K0s_Lam_TM_TU","hpeakd1_K0s_Lam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_Lam_TM_TU  = new THnSparseD("hpeakd2_K0s_Lam_TM_TU","hpeakd2_K0s_Lam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_Lam_TM_TU  = new THnSparseD("hpeakd12_K0s_Lam_TM_TU","hpeakd12_K0s_Lam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_Lam_TM_TU_bias  = new THnSparseD("hpeak_K0s_Lam_TM_TU_bias","hpeak_K0s_Lam_TM_TU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TM_TU_biasd1  = new THnSparseD("hpeak_K0s_Lam_TM_TU_biasd1","hpeak_K0s_Lam_TM_TU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TM_TU_biasd2  = new THnSparseD("hpeak_K0s_Lam_TM_TU_biasd2","hpeak_K0s_Lam_TM_TU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TM_TU_biasd12  = new THnSparseD("hpeak_K0s_Lam_TM_TU_biasd12","hpeak_K0s_Lam_TM_TU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_TU_K0s_Lam=new TH1D("V0chi2d1_diff_TM_TU_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_TU_K0s_Lam=new TH1D("V0chi2d2_diff_TM_TU_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_TU_K0s_Lam=new TH1D("V0chi2d12_diff_TM_TU_K0s_Lam","",10000,-0.01,0.11);

//fake matched + fake unmatched

THnSparseD *hpeak_K0s_Lam_FM_FU  = new THnSparseD("hpeak_K0s_Lam_FM_FU","hpeak_K0s_Lam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_Lam_FM_FU  = new THnSparseD("hpeakd1_K0s_Lam_FM_FU","hpeakd1_K0s_Lam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_Lam_FM_FU  = new THnSparseD("hpeakd2_K0s_Lam_FM_FU","hpeakd2_K0s_Lam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_Lam_FM_FU  = new THnSparseD("hpeakd12_K0s_Lam_FM_FU","hpeakd12_K0s_Lam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_Lam_FM_FU_bias  = new THnSparseD("hpeak_K0s_Lam_FM_FU_bias","hpeak_K0s_Lam_FM_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_FM_FU_biasd1  = new THnSparseD("hpeak_K0s_Lam_FM_FU_biasd1","hpeak_K0s_Lam_FM_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_FM_FU_biasd2  = new THnSparseD("hpeak_K0s_Lam_FM_FU_biasd2","hpeak_K0s_Lam_FM_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_FM_FU_biasd12  = new THnSparseD("hpeak_K0s_Lam_FM_FU_biasd12","hpeak_K0s_Lam_FM_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FM_FU_K0s_Lam=new TH1D("V0chi2d1_diff_FM_FU_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FM_FU_K0s_Lam=new TH1D("V0chi2d2_diff_FM_FU_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FM_FU_K0s_Lam=new TH1D("V0chi2d12_diff_FM_FU_K0s_Lam","",10000,-0.01,0.11);

//true matched + fake matched

THnSparseD *hpeak_K0s_Lam_TM_FM  = new THnSparseD("hpeak_K0s_Lam_TM_FM","hpeak_K0s_Lam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_Lam_TM_FM  = new THnSparseD("hpeakd1_K0s_Lam_TM_FM","hpeakd1_K0s_Lam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_Lam_TM_FM  = new THnSparseD("hpeakd2_K0s_Lam_TM_FM","hpeakd2_K0s_Lam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_Lam_TM_FM  = new THnSparseD("hpeakd12_K0s_Lam_TM_FM","hpeakd12_K0s_Lam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_Lam_TM_FM_bias  = new THnSparseD("hpeak_K0s_Lam_TM_FM_bias","hpeak_K0s_Lam_TM_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TM_FM_biasd1  = new THnSparseD("hpeak_K0s_Lam_TM_FM_biasd1","hpeak_K0s_Lam_TM_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TM_FM_biasd2  = new THnSparseD("hpeak_K0s_Lam_TM_FM_biasd2","hpeak_K0s_Lam_TM_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TM_FM_biasd12  = new THnSparseD("hpeak_K0s_Lam_TM_FM_biasd12","hpeak_K0s_Lam_TM_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_FM_K0s_Lam=new TH1D("V0chi2d1_diff_TM_FM_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_FM_K0s_Lam=new TH1D("V0chi2d2_diff_TM_FM_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_FM_K0s_Lam=new TH1D("V0chi2d12_diff_TM_FM_K0s_Lam","",10000,-0.01,0.11);

//true matched + fake unmatched

THnSparseD *hpeak_K0s_Lam_TM_FU  = new THnSparseD("hpeak_K0s_Lam_TM_FU","hpeak_K0s_Lam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_Lam_TM_FU  = new THnSparseD("hpeakd1_K0s_Lam_TM_FU","hpeakd1_K0s_Lam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_Lam_TM_FU  = new THnSparseD("hpeakd2_K0s_Lam_TM_FU","hpeakd2_K0s_Lam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_Lam_TM_FU  = new THnSparseD("hpeakd12_K0s_Lam_TM_FU","hpeakd12_K0s_Lam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_Lam_TM_FU_bias  = new THnSparseD("hpeak_K0s_Lam_TM_FU_bias","hpeak_K0s_Lam_TM_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TM_FU_biasd1  = new THnSparseD("hpeak_K0s_Lam_TM_FU_biasd1","hpeak_K0s_Lam_TM_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TM_FU_biasd2  = new THnSparseD("hpeak_K0s_Lam_TM_FU_biasd2","hpeak_K0s_Lam_TM_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TM_FU_biasd12  = new THnSparseD("hpeak_K0s_Lam_TM_FU_biasd12","hpeak_K0s_Lam_TM_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_FU_K0s_Lam=new TH1D("V0chi2d1_diff_TM_FU_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_FU_K0s_Lam=new TH1D("V0chi2d2_diff_TM_FU_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_FU_K0s_Lam=new TH1D("V0chi2d12_diff_TM_FU_K0s_Lam","",10000,-0.01,0.11);

//true unmatched + fake unmatched

THnSparseD *hpeak_K0s_Lam_TU_FU  = new THnSparseD("hpeak_K0s_Lam_TU_FU","hpeak_K0s_Lam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_Lam_TU_FU  = new THnSparseD("hpeakd1_K0s_Lam_TU_FU","hpeakd1_K0s_Lam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_Lam_TU_FU  = new THnSparseD("hpeakd2_K0s_Lam_TU_FU","hpeakd2_K0s_Lam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_Lam_TU_FU  = new THnSparseD("hpeakd12_K0s_Lam_TU_FU","hpeakd12_K0s_Lam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_Lam_TU_FU_bias  = new THnSparseD("hpeak_K0s_Lam_TU_FU_bias","hpeak_K0s_Lam_TU_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TU_FU_biasd1  = new THnSparseD("hpeak_K0s_Lam_TU_FU_biasd1","hpeak_K0s_Lam_TU_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TU_FU_biasd2  = new THnSparseD("hpeak_K0s_Lam_TU_FU_biasd2","hpeak_K0s_Lam_TU_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TU_FU_biasd12  = new THnSparseD("hpeak_K0s_Lam_TU_FU_biasd12","hpeak_K0s_Lam_TU_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_FU_K0s_Lam=new TH1D("V0chi2d1_diff_TU_FU_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_FU_K0s_Lam=new TH1D("V0chi2d2_diff_TU_FU_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_FU_K0s_Lam=new TH1D("V0chi2d12_diff_TU_FU_K0s_Lam","",10000,-0.01,0.11);

//true unmatched + fake matched

THnSparseD *hpeak_K0s_Lam_TU_FM  = new THnSparseD("hpeak_K0s_Lam_TU_FM","hpeak_K0s_Lam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_Lam_TU_FM  = new THnSparseD("hpeakd1_K0s_Lam_TU_FM","hpeakd1_K0s_Lam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_Lam_TU_FM  = new THnSparseD("hpeakd2_K0s_Lam_TU_FM","hpeakd2_K0s_Lam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_Lam_TU_FM  = new THnSparseD("hpeakd12_K0s_Lam_TU_FM","hpeakd12_K0s_Lam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_Lam_TU_FM_bias  = new THnSparseD("hpeak_K0s_Lam_TU_FM_bias","hpeak_K0s_Lam_TU_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TU_FM_biasd1  = new THnSparseD("hpeak_K0s_Lam_TU_FM_biasd1","hpeak_K0s_Lam_TU_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TU_FM_biasd2  = new THnSparseD("hpeak_K0s_Lam_TU_FM_biasd2","hpeak_K0s_Lam_TU_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_Lam_TU_FM_biasd12  = new THnSparseD("hpeak_K0s_Lam_TU_FM_biasd12","hpeak_K0s_Lam_TU_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_FM_K0s_Lam=new TH1D("V0chi2d1_diff_TU_FM_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_FM_K0s_Lam=new TH1D("V0chi2d2_diff_TU_FM_K0s_Lam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_FM_K0s_Lam=new TH1D("V0chi2d12_diff_TU_FM_K0s_Lam","",10000,-0.01,0.11);


//K0sALam

  THnSparseD *hMass_K0sAL_T  = new THnSparseD("hMass_K0sAL_T","hMass_K0sAL_T",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
  THnSparseD *hMass_ALamK_T  = new THnSparseD("hMass_ALamK_T","hMass_ALamK_T",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
  THnSparseD *hMass_K0s_ALam_T  = new THnSparseD("hMass_K0s_ALam_T","hMass_K0s_ALam_T",4,bins3DM_KL,xmin3DM_KL,xmax3DM_KL);   
  
  THnSparseD *hpeak_K0s_ALam_T  = new THnSparseD("hpeak_K0s_ALam_T","hpeak_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakcos_K0s_ALam_T  = new THnSparseD("hpeakcos_K0s_ALam_T","hpeakcos_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeak_rot_K0s_ALam_T  = new THnSparseD("hpeak_rot_K0s_ALam_T","hpeak_rot_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeak_inv_K0s_ALam_T  = new THnSparseD("hpeak_inv_K0s_ALam_T","hpeak_inv_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_K0s_ALam_T  = new THnSparseD("hside_K0s_ALam_T","hside_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideL_K0s_ALam_T  = new THnSparseD("hsideL_K0s_ALam_T","hsideL_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideR_K0s_ALam_T  = new THnSparseD("hsideR_K0s_ALam_T","hsideR_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_rot_K0s_ALam_T  = new THnSparseD("hside_rot_K0s_ALam_T","hside_rot_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_inv_K0s_ALam_T  = new THnSparseD("hside_inv_K0s_ALam_T","hside_inv_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_K0s_ALam_T  = new THnSparseD("hpeakside_K0s_ALam_T","hpeakside_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideL_K0s_ALam_T  = new THnSparseD("hpeaksideL_K0s_ALam_T","hpeaksideL_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideR_K0s_ALam_T  = new THnSparseD("hpeaksideR_K0s_ALam_T","hpeaksideR_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_rot_K0s_ALam_T  = new THnSparseD("hpeakside_rot_K0s_ALam_T","hpeakside_rot_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_inv_K0s_ALam_T  = new THnSparseD("hpeakside_inv_K0s_ALam_T","hpeakside_inv_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd1_K0s_ALam_T  = new THnSparseD("hpeakd1_K0s_ALam_T","hpeakd1_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd2_K0s_ALam_T  = new THnSparseD("hpeakd2_K0s_ALam_T","hpeakd2_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd12_K0s_ALam_T  = new THnSparseD("hpeakd12_K0s_ALam_T","hpeakd12_K0s_ALam_T",3,bins3D,xmin3D,xmax3D);

  THnSparseD *hpeak_K0s_ALam_T_mix  = new THnSparseD("hpeak_K0s_ALam_T_mix","hpeak_K0s_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_K0s_ALam_T_mix  = new THnSparseD("hside_K0s_ALam_T_mix","hside_K0s_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideL_K0s_ALam_T_mix  = new THnSparseD("hsideL_K0s_ALam_T_mix","hsideL_K0s_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideR_K0s_ALam_T_mix  = new THnSparseD("hsideR_K0s_ALam_mix","hsideR_K0s_ALam_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_K0s_ALam_T_mix  = new THnSparseD("hpeakside_K0s_ALam_T_mix","hpeakside_K0s_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideL_K0s_ALam_T_mix  = new THnSparseD("hpeaksideL_K0s_ALam_T_mix","hpeaksideL_K0s_ALam_T_mix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideR_K0s_ALam_T_mix  = new THnSparseD("hpeaksideR_K0s_ALam_T_mix","hpeaksideR_K0s_ALam_T_mix",3,bins3D,xmin3D,xmax3D);

  THnSparseD *hpeak_K0s_ALam_T_etamix  = new THnSparseD("hpeak_K0s_ALam_T_etamix","hpeak_K0s_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_K0s_ALam_T_etamix  = new THnSparseD("hside_K0s_ALam_T_etamix","hside_K0s_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideL_K0s_ALam_T_etamix  = new THnSparseD("hsideL_K0s_ALam_T_etamix","hsideL_K0s_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hsideR_K0s_ALam_T_etamix  = new THnSparseD("hsideR_K0s_ALam_T_etamix","hsideR_K0s_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_K0s_ALam_T_etamix  = new THnSparseD("hpeakside_K0s_ALam_T_etamix","hpeakside_K0s_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideL_K0s_ALam_T_etamix  = new THnSparseD("hpeaksideL_K0s_ALam_T_etamix","hpeaksideL_K0s_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeaksideR_K0s_ALam_T_etamix  = new THnSparseD("hpeaksideR_K0s_ALam_T_etamix","hpeaksideR_K0s_ALam_T_etamix",3,bins3D,xmin3D,xmax3D);  
  
  
  TH1D *V0chi2d1_diff_T_K0s_ALam=new TH1D("V0chi2d1_diff_T_K0s_ALam","",10000,-0.01,0.11);
  TH1D *V0chi2d2_diff_T_K0s_ALam=new TH1D("V0chi2d2_diff_T_K0s_ALam","",10000,-0.01,0.11);
  TH1D *V0chi2d12_diff_T_K0s_ALam=new TH1D("V0chi2d12_diff_T_K0s_ALam","",10000,-0.01,0.11);
  
//------------------------------------------------------------------------
//fake
//------------------------------------------------------------------------ 

  THnSparseD *hMass_K0sAL_F  = new THnSparseD("hMass_K0sAL_F","hMass_K0sAL_F",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
  THnSparseD *hMass_ALamK_F  = new THnSparseD("hMass_ALamK_F","hMass_ALamK_F",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
  THnSparseD *hMass_K0s_ALam_F  = new THnSparseD("hMass_K0s_ALam_F","hMass_K0s_ALam_F",4,bins3DM_KL,xmin3DM_KL,xmax3DM_KL);   
  
  THnSparseD *hpeak_K0s_ALam_F  = new THnSparseD("hpeak_K0s_ALam_F","hpeak_K0s_ALam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd1_K0s_ALam_F  = new THnSparseD("hpeakd1_K0s_ALam_F","hpeakd1_K0s_ALam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd2_K0s_ALam_F  = new THnSparseD("hpeakd2_K0s_ALam_F","hpeakd2_K0s_ALam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakd12_K0s_ALam_F  = new THnSparseD("hpeakd12_K0s_ALam_F","hpeakd12_K0s_ALam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hside_K0s_ALam_F  = new THnSparseD("hside_K0s_ALam_F","hside_K0s_ALam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeakside_K0s_ALam_F  = new THnSparseD("hpeakside_K0s_ALam_F","hpeakside_K0s_ALam_F",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hpeak_K0s_ALam_F_bias  = new THnSparseD("hpeak_K0s_ALam_F_bias","hpeak_K0s_ALam_F_bias",3,bins3Db,xmin3Db,xmax3Db);
  THnSparseD *hpeak_K0s_ALam_F_biasd1  = new THnSparseD("hpeak_K0s_ALam_F_biasd1","hpeak_K0s_ALam_F_biasd1",3,bins3Db,xmin3Db,xmax3Db);
  THnSparseD *hpeak_K0s_ALam_F_biasd2  = new THnSparseD("hpeak_K0s_ALam_F_biasd2","hpeak_K0s_ALam_F_biasd2",3,bins3Db,xmin3Db,xmax3Db);
  THnSparseD *hpeak_K0s_ALam_F_biasd12  = new THnSparseD("hpeak_K0s_ALam_F_biasd12","hpeak_K0s_ALam_F_biasd12",3,bins3Db,xmin3Db,xmax3Db);

  TH1D *V0chi2d1_diff_F_K0s_ALam=new TH1D("V0chi2d1_diff_F_K0s_ALam","",10000,-0.01,0.11);
  TH1D *V0chi2d2_diff_F_K0s_ALam=new TH1D("V0chi2d2_diff_F_K0s_ALam","",10000,-0.01,0.11);
  TH1D *V0chi2d12_diff_F_K0s_ALam=new TH1D("V0chi2d12_diff_F_K0s_ALam","",10000,-0.01,0.11);
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
THnSparseD *hpeak_K0s_ALam_TF  = new THnSparseD("hpeak_K0s_ALam_TF","hpeak_K0s_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_ALam_TF  = new THnSparseD("hpeakd1_K0s_ALam_TF","hpeakd1_K0s_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_ALam_TF  = new THnSparseD("hpeakd2_K0s_ALam_TF","hpeakd2_K0s_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_ALam_TF  = new THnSparseD("hpeakd12_K0s_ALam_TF","hpeakd12_K0s_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hside_K0s_ALam_TF  = new THnSparseD("hside_K0s_ALam_TF","hside_K0s_ALam_TF",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakside_K0s_ALam_TF  = new THnSparseD("hpeakside_K0s_ALam_TF","hpeakside_K0s_ALam_TF",3,bins3D,xmin3D,xmax3D);  
THnSparseD *hpeak_K0s_ALam_TF_bias  = new THnSparseD("hpeak_K0s_ALam_TF_bias","hpeak_K0s_ALam_TF_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TF_biasd1  = new THnSparseD("hpeak_K0s_ALam_TF_biasd1","hpeak_K0s_ALam_TF_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TF_biasd2  = new THnSparseD("hpeak_K0s_ALam_TF_biasd2","hpeak_K0s_ALam_TF_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TF_biasd12  = new THnSparseD("hpeak_K0s_ALam_TF_biasd12","hpeak_K0s_ALam_TF_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TF_K0s_ALam=new TH1D("V0chi2d1_diff_TF_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TF_K0s_ALam=new TH1D("V0chi2d2_diff_TF_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TF_K0s_ALam=new TH1D("V0chi2d12_diff_TF_K0s_ALam","",10000,-0.01,0.11);

  
/////// matching hist ////// 

THnSparseD *hMass_K0sAL_TM  = new THnSparseD("hMass_K0sAL_TM","hMass_K0sAL_TM",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_ALamK_TM  = new THnSparseD("hMass_ALamK_TM","hMass_ALamK_TM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_K0s_ALam_TM  = new THnSparseD("hMass_K0s_ALam_TM","hMass_K0s_ALam_TM",4,bins3DM_KL,xmin3DM_KL,xmax3DM_KL);   

THnSparseD *hMass_K0sAL_TU  = new THnSparseD("hMass_K0sAL_TU","hMass_K0sAL_TU",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_ALamK_TU  = new THnSparseD("hMass_ALamK_TU","hMass_ALamK_TU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_K0s_ALam_TU  = new THnSparseD("hMass_K0s_ALam_TU","hMass_K0s_ALam_TU",4,bins3DM_KL,xmin3DM_KL,xmax3DM_KL);   

THnSparseD *hMass_K0sAL_FM  = new THnSparseD("hMass_K0sAL_FM","hMass_K0sAL_FM",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_ALamK_FM  = new THnSparseD("hMass_ALamK_FM","hMass_ALamK_FM",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_K0s_ALam_FM  = new THnSparseD("hMass_K0s_ALam_FM","hMass_K0s_ALam_FM",4,bins3DM_KL,xmin3DM_KL,xmax3DM_KL);   

THnSparseD *hMass_K0sAL_FU  = new THnSparseD("hMass_K0sAL_FU","hMass_K0sAL_FU",4,bins3DM_K0s,xmin3DM_K0s,xmax3DM_K0s);   
THnSparseD *hMass_ALamK_FU  = new THnSparseD("hMass_ALamK_FU","hMass_ALamK_FU",4,bins3DM_LAL,xmin3DM_LAL,xmax3DM_LAL);   
THnSparseD *hMass_K0s_ALam_FU  = new THnSparseD("hMass_K0s_ALam_FU","hMass_K0s_ALam_FU",4,bins3DM_KL,xmin3DM_KL,xmax3DM_KL);   


THnSparseD *hpeak_K0s_ALam_TM  = new THnSparseD("hpeak_K0s_ALam_TM","hpeak_K0s_ALam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_ALam_TM  = new THnSparseD("hpeakd1_K0s_ALam_TM","hpeakd1_K0s_ALam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_ALam_TM  = new THnSparseD("hpeakd2_K0s_ALam_TM","hpeakd2_K0s_ALam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_ALam_TM  = new THnSparseD("hpeakd12_K0s_ALam_TM","hpeakd12_K0s_ALam_TM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_ALam_TM_bias  = new THnSparseD("hpeak_K0s_ALam_TM_bias","hpeak_K0s_ALam_TM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TM_biasd1  = new THnSparseD("hpeak_K0s_ALam_TM_biasd1","hpeak_K0s_ALam_TM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TM_biasd2  = new THnSparseD("hpeak_K0s_ALam_TM_biasd2","hpeak_K0s_ALam_TM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TM_biasd12  = new THnSparseD("hpeak_K0s_ALam_TM_biasd12","hpeak_K0s_ALam_TM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_K0s_ALam=new TH1D("V0chi2d1_diff_TM_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_K0s_ALam=new TH1D("V0chi2d2_diff_TM_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_K0s_ALam=new TH1D("V0chi2d12_diff_TM_K0s_ALam","",10000,-0.01,0.11);

THnSparseD *hpeak_K0s_ALam_TU  = new THnSparseD("hpeak_K0s_ALam_TU","hpeak_K0s_ALam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_ALam_TU  = new THnSparseD("hpeakd1_K0s_ALam_TU","hpeakd1_K0s_ALam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_ALam_TU  = new THnSparseD("hpeakd2_K0s_ALam_TU","hpeakd2_K0s_ALam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_ALam_TU  = new THnSparseD("hpeakd12_K0s_ALam_TU","hpeakd12_K0s_ALam_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_ALam_TU_bias  = new THnSparseD("hpeak_K0s_ALam_TU_bias","hpeak_K0s_ALam_TU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TU_biasd1  = new THnSparseD("hpeak_K0s_ALam_TU_biasd1","hpeak_K0s_ALam_TU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TU_biasd2  = new THnSparseD("hpeak_K0s_ALam_TU_biasd2","hpeak_K0s_ALam_TU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TU_biasd12  = new THnSparseD("hpeak_K0s_ALam_TU_biasd12","hpeak_K0s_ALam_TU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_K0s_ALam=new TH1D("V0chi2d1_diff_TU_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_K0s_ALam=new TH1D("V0chi2d2_diff_TU_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_K0s_ALam=new TH1D("V0chi2d12_diff_TU_K0s_ALam","",10000,-0.01,0.11);

THnSparseD *hpeak_K0s_ALam_FM  = new THnSparseD("hpeak_K0s_ALam_FM","hpeak_K0s_ALam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_ALam_FM  = new THnSparseD("hpeakd1_K0s_ALam_FM","hpeakd1_K0s_ALam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_ALam_FM  = new THnSparseD("hpeakd2_K0s_ALam_FM","hpeakd2_K0s_ALam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_ALam_FM  = new THnSparseD("hpeakd12_K0s_ALam_FM","hpeakd12_K0s_ALam_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_ALam_FM_bias  = new THnSparseD("hpeak_K0s_ALam_FM_bias","hpeak_K0s_ALam_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_FM_biasd1  = new THnSparseD("hpeak_K0s_ALam_FM_biasd1","hpeak_K0s_ALam_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_FM_biasd2  = new THnSparseD("hpeak_K0s_ALam_FM_biasd2","hpeak_K0s_ALam_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_FM_biasd12  = new THnSparseD("hpeak_K0s_ALam_FM_biasd12","hpeak_K0s_ALam_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FM_K0s_ALam=new TH1D("V0chi2d1_diff_FM_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FM_K0s_ALam=new TH1D("V0chi2d2_diff_FM_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FM_K0s_ALam=new TH1D("V0chi2d12_diff_FM_K0s_ALam","",10000,-0.01,0.11);

THnSparseD *hpeak_K0s_ALam_FU  = new THnSparseD("hpeak_K0s_ALam_FU","hpeak_K0s_ALam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_ALam_FU  = new THnSparseD("hpeakd1_K0s_ALam_FU","hpeakd1_K0s_ALam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_ALam_FU  = new THnSparseD("hpeakd2_K0s_ALam_FU","hpeakd2_K0s_ALam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_ALam_FU  = new THnSparseD("hpeakd12_K0s_ALam_FU","hpeakd12_K0s_ALam_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_ALam_FU_bias  = new THnSparseD("hpeak_K0s_ALam_FU_bias","hpeak_K0s_ALam_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_FU_biasd1  = new THnSparseD("hpeak_K0s_ALam_FU_biasd1","hpeak_K0s_ALam_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_FU_biasd2  = new THnSparseD("hpeak_K0s_ALam_FU_biasd2","hpeak_K0s_ALam_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_FU_biasd12  = new THnSparseD("hpeak_K0s_ALam_FU_biasd12","hpeak_K0s_ALam_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FU_K0s_ALam=new TH1D("V0chi2d1_diff_FU_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FU_K0s_ALam=new TH1D("V0chi2d2_diff_FU_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FU_K0s_ALam=new TH1D("V0chi2d12_diff_FU_K0s_ALam","",10000,-0.01,0.11);

//true matched + true unmatched

THnSparseD *hpeak_K0s_ALam_TM_TU  = new THnSparseD("hpeak_K0s_ALam_TM_TU","hpeak_K0s_ALam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_ALam_TM_TU  = new THnSparseD("hpeakd1_K0s_ALam_TM_TU","hpeakd1_K0s_ALam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_ALam_TM_TU  = new THnSparseD("hpeakd2_K0s_ALam_TM_TU","hpeakd2_K0s_ALam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_ALam_TM_TU  = new THnSparseD("hpeakd12_K0s_ALam_TM_TU","hpeakd12_K0s_ALam_TM_TU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_ALam_TM_TU_bias  = new THnSparseD("hpeak_K0s_ALam_TM_TU_bias","hpeak_K0s_ALam_TM_TU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TM_TU_biasd1  = new THnSparseD("hpeak_K0s_ALam_TM_TU_biasd1","hpeak_K0s_ALam_TM_TU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TM_TU_biasd2  = new THnSparseD("hpeak_K0s_ALam_TM_TU_biasd2","hpeak_K0s_ALam_TM_TU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TM_TU_biasd12  = new THnSparseD("hpeak_K0s_ALam_TM_TU_biasd12","hpeak_K0s_ALam_TM_TU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_TU_K0s_ALam=new TH1D("V0chi2d1_diff_TM_TU_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_TU_K0s_ALam=new TH1D("V0chi2d2_diff_TM_TU_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_TU_K0s_ALam=new TH1D("V0chi2d12_diff_TM_TU_K0s_ALam","",10000,-0.01,0.11);

//fake matched + fake unmatched

THnSparseD *hpeak_K0s_ALam_FM_FU  = new THnSparseD("hpeak_K0s_ALam_FM_FU","hpeak_K0s_ALam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_ALam_FM_FU  = new THnSparseD("hpeakd1_K0s_ALam_FM_FU","hpeakd1_K0s_ALam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_ALam_FM_FU  = new THnSparseD("hpeakd2_K0s_ALam_FM_FU","hpeakd2_K0s_ALam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_ALam_FM_FU  = new THnSparseD("hpeakd12_K0s_ALam_FM_FU","hpeakd12_K0s_ALam_FM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_ALam_FM_FU_bias  = new THnSparseD("hpeak_K0s_ALam_FM_FU_bias","hpeak_K0s_ALam_FM_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_FM_FU_biasd1  = new THnSparseD("hpeak_K0s_ALam_FM_FU_biasd1","hpeak_K0s_ALam_FM_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_FM_FU_biasd2  = new THnSparseD("hpeak_K0s_ALam_FM_FU_biasd2","hpeak_K0s_ALam_FM_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_FM_FU_biasd12  = new THnSparseD("hpeak_K0s_ALam_FM_FU_biasd12","hpeak_K0s_ALam_FM_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_FM_FU_K0s_ALam=new TH1D("V0chi2d1_diff_FM_FU_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_FM_FU_K0s_ALam=new TH1D("V0chi2d2_diff_FM_FU_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_FM_FU_K0s_ALam=new TH1D("V0chi2d12_diff_FM_FU_K0s_ALam","",10000,-0.01,0.11);

//true matched + fake matched

THnSparseD *hpeak_K0s_ALam_TM_FM  = new THnSparseD("hpeak_K0s_ALam_TM_FM","hpeak_K0s_ALam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_ALam_TM_FM  = new THnSparseD("hpeakd1_K0s_ALam_TM_FM","hpeakd1_K0s_ALam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_ALam_TM_FM  = new THnSparseD("hpeakd2_K0s_ALam_TM_FM","hpeakd2_K0s_ALam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_ALam_TM_FM  = new THnSparseD("hpeakd12_K0s_ALam_TM_FM","hpeakd12_K0s_ALam_TM_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_ALam_TM_FM_bias  = new THnSparseD("hpeak_K0s_ALam_TM_FM_bias","hpeak_K0s_ALam_TM_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TM_FM_biasd1  = new THnSparseD("hpeak_K0s_ALam_TM_FM_biasd1","hpeak_K0s_ALam_TM_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TM_FM_biasd2  = new THnSparseD("hpeak_K0s_ALam_TM_FM_biasd2","hpeak_K0s_ALam_TM_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TM_FM_biasd12  = new THnSparseD("hpeak_K0s_ALam_TM_FM_biasd12","hpeak_K0s_ALam_TM_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_FM_K0s_ALam=new TH1D("V0chi2d1_diff_TM_FM_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_FM_K0s_ALam=new TH1D("V0chi2d2_diff_TM_FM_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_FM_K0s_ALam=new TH1D("V0chi2d12_diff_TM_FM_K0s_ALam","",10000,-0.01,0.11);

//true matched + fake unmatched

THnSparseD *hpeak_K0s_ALam_TM_FU  = new THnSparseD("hpeak_K0s_ALam_TM_FU","hpeak_K0s_ALam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_ALam_TM_FU  = new THnSparseD("hpeakd1_K0s_ALam_TM_FU","hpeakd1_K0s_ALam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_ALam_TM_FU  = new THnSparseD("hpeakd2_K0s_ALam_TM_FU","hpeakd2_K0s_ALam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_ALam_TM_FU  = new THnSparseD("hpeakd12_K0s_ALam_TM_FU","hpeakd12_K0s_ALam_TM_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_ALam_TM_FU_bias  = new THnSparseD("hpeak_K0s_ALam_TM_FU_bias","hpeak_K0s_ALam_TM_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TM_FU_biasd1  = new THnSparseD("hpeak_K0s_ALam_TM_FU_biasd1","hpeak_K0s_ALam_TM_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TM_FU_biasd2  = new THnSparseD("hpeak_K0s_ALam_TM_FU_biasd2","hpeak_K0s_ALam_TM_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TM_FU_biasd12  = new THnSparseD("hpeak_K0s_ALam_TM_FU_biasd12","hpeak_K0s_ALam_TM_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TM_FU_K0s_ALam=new TH1D("V0chi2d1_diff_TM_FU_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TM_FU_K0s_ALam=new TH1D("V0chi2d2_diff_TM_FU_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TM_FU_K0s_ALam=new TH1D("V0chi2d12_diff_TM_FU_K0s_ALam","",10000,-0.01,0.11);

//true unmatched + fake unmatched

THnSparseD *hpeak_K0s_ALam_TU_FU  = new THnSparseD("hpeak_K0s_ALam_TU_FU","hpeak_K0s_ALam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_ALam_TU_FU  = new THnSparseD("hpeakd1_K0s_ALam_TU_FU","hpeakd1_K0s_ALam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_ALam_TU_FU  = new THnSparseD("hpeakd2_K0s_ALam_TU_FU","hpeakd2_K0s_ALam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_ALam_TU_FU  = new THnSparseD("hpeakd12_K0s_ALam_TU_FU","hpeakd12_K0s_ALam_TU_FU",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_ALam_TU_FU_bias  = new THnSparseD("hpeak_K0s_ALam_TU_FU_bias","hpeak_K0s_ALam_TU_FU_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TU_FU_biasd1  = new THnSparseD("hpeak_K0s_ALam_TU_FU_biasd1","hpeak_K0s_ALam_TU_FU_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TU_FU_biasd2  = new THnSparseD("hpeak_K0s_ALam_TU_FU_biasd2","hpeak_K0s_ALam_TU_FU_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TU_FU_biasd12  = new THnSparseD("hpeak_K0s_ALam_TU_FU_biasd12","hpeak_K0s_ALam_TU_FU_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_FU_K0s_ALam=new TH1D("V0chi2d1_diff_TU_FU_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_FU_K0s_ALam=new TH1D("V0chi2d2_diff_TU_FU_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_FU_K0s_ALam=new TH1D("V0chi2d12_diff_TU_FU_K0s_ALam","",10000,-0.01,0.11);

//true unmatched + fake matched

THnSparseD *hpeak_K0s_ALam_TU_FM  = new THnSparseD("hpeak_K0s_ALam_TU_FM","hpeak_K0s_ALam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd1_K0s_ALam_TU_FM  = new THnSparseD("hpeakd1_K0s_ALam_TU_FM","hpeakd1_K0s_ALam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd2_K0s_ALam_TU_FM  = new THnSparseD("hpeakd2_K0s_ALam_TU_FM","hpeakd2_K0s_ALam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeakd12_K0s_ALam_TU_FM  = new THnSparseD("hpeakd12_K0s_ALam_TU_FM","hpeakd12_K0s_ALam_TU_FM",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_ALam_TU_FM_bias  = new THnSparseD("hpeak_K0s_ALam_TU_FM_bias","hpeak_K0s_ALam_TU_FM_bias",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TU_FM_biasd1  = new THnSparseD("hpeak_K0s_ALam_TU_FM_biasd1","hpeak_K0s_ALam_TU_FM_biasd1",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TU_FM_biasd2  = new THnSparseD("hpeak_K0s_ALam_TU_FM_biasd2","hpeak_K0s_ALam_TU_FM_biasd2",3,bins3Db,xmin3Db,xmax3Db);
THnSparseD *hpeak_K0s_ALam_TU_FM_biasd12  = new THnSparseD("hpeak_K0s_ALam_TU_FM_biasd12","hpeak_K0s_ALam_TU_FM_biasd12",3,bins3Db,xmin3Db,xmax3Db);

TH1D *V0chi2d1_diff_TU_FM_K0s_ALam=new TH1D("V0chi2d1_diff_TU_FM_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d2_diff_TU_FM_K0s_ALam=new TH1D("V0chi2d2_diff_TU_FM_K0s_ALam","",10000,-0.01,0.11);
TH1D *V0chi2d12_diff_TU_FM_K0s_ALam=new TH1D("V0chi2d12_diff_TU_FM_K0s_ALam","",10000,-0.01,0.11);

  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//for MC shape

THnSparseD *hpeak_K0s_K0s_TM_mix  = new THnSparseD("hpeak_K0s_K0s_TM_mix","hpeak_K0s_K0s_TM_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_TU_mix  = new THnSparseD("hpeak_K0s_K0s_TU_mix","hpeak_K0s_K0s_TU_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_K0s_TM_TU_mix  = new THnSparseD("hpeak_K0s_K0s_TM_TU_mix","hpeak_K0s_K0s_TM_TU_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_Lam_Lam_TM_mix  = new THnSparseD("hpeak_Lam_Lam_TM_mix","hpeak_Lam_Lam_TM_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_TU_mix  = new THnSparseD("hpeak_Lam_Lam_TU_mix","hpeak_Lam_Lam_TU_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_Lam_TM_TU_mix  = new THnSparseD("hpeak_Lam_Lam_TM_TU_mix","hpeak_Lam_Lam_TM_TU_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_ALam_ALam_TM_mix  = new THnSparseD("hpeak_ALam_ALam_TM_mix","hpeak_ALam_ALam_TM_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_TU_mix  = new THnSparseD("hpeak_ALam_ALam_TU_mix","hpeak_ALam_ALam_TU_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALam_TM_TU_mix  = new THnSparseD("hpeak_ALam_ALam_TM_TU_mix","hpeak_ALam_ALam_TM_TU_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_Lam_ALam_TM_mix  = new THnSparseD("hpeak_Lam_ALam_TM_mix","hpeak_Lam_ALam_TM_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_ALam_TU_mix  = new THnSparseD("hpeak_Lam_ALam_TU_mix","hpeak_Lam_ALam_TU_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_ALam_TM_TU_mix  = new THnSparseD("hpeak_Lam_ALam_TM_TU_mix","hpeak_Lam_ALam_TM_TU_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_K0s_Lam_TM_mix  = new THnSparseD("hpeak_K0s_Lam_TM_mix","hpeak_K0s_Lam_TM_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_Lam_TU_mix  = new THnSparseD("hpeak_K0s_Lam_TU_mix","hpeak_K0s_Lam_TU_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_Lam_TM_TU_mix  = new THnSparseD("hpeak_K0s_Lam_TM_TU_mix","hpeak_K0s_Lam_TM_TU_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_K0s_ALam_TM_mix  = new THnSparseD("hpeak_K0s_ALam_TM_mix","hpeak_K0s_ALam_TM_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_ALam_TU_mix  = new THnSparseD("hpeak_K0s_ALam_TU_mix","hpeak_K0s_ALam_TU_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_K0s_ALam_TM_TU_mix  = new THnSparseD("hpeak_K0s_ALam_TM_TU_mix","hpeak_K0s_ALam_TM_TU_mix",3,bins3D,xmin3D,xmax3D);


//this are the histograms from feed down

THnSparseD *hpeak_Lam_LamFD  = new THnSparseD("hpeak_Lam_LamFD","hpeak_Lam_LamFD",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_LamFD_LamFD  = new THnSparseD("hpeak_LamFD_LamFD","hpeak_LamFD_LamFD",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_Lam_LamFD_mix  = new THnSparseD("hpeak_Lam_LamFD_mix","hpeak_Lam_LamFD_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_LamFD_LamFD_mix  = new THnSparseD("hpeak_LamFD_LamFD_mix","hpeak_LamFD_LamFD_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_ALam_ALamFD  = new THnSparseD("hpeak_ALam_ALamFD","hpeak_ALam_ALamFD",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALamFD_ALamFD  = new THnSparseD("hpeak_ALamFD_ALamFD","hpeak_ALamFD_ALamFD",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALam_ALamFD_mix  = new THnSparseD("hpeak_ALam_ALamFD_mix","hpeak_ALam_ALamFD_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_ALamFD_ALamFD_mix  = new THnSparseD("hpeak_ALamFD_ALamFD_mix","hpeak_ALamFD_ALamFD_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_LALFD  = new THnSparseD("hpeak_LALFD","hpeak_LALFD",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_LALFD_mix  = new THnSparseD("hpeak_LALFD_mix","hpeak_LALFD_mix",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_LFDALFD  = new THnSparseD("hpeak_LFDALFD","hpeak_LFDALFD",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_LFDALFD_mix  = new THnSparseD("hpeak_LFDALFD_mix","hpeak_LFDALFD_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_KLFD  = new THnSparseD("hpeak_KLFD","hpeak_KLFD",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_KLFD_mix  = new THnSparseD("hpeak_KLFD_mix","hpeak_KLFD_mix",3,bins3D,xmin3D,xmax3D);

THnSparseD *hpeak_KALFD  = new THnSparseD("hpeak_KALFD","hpeak_KALFD",3,bins3D,xmin3D,xmax3D);
THnSparseD *hpeak_KALFD_mix  = new THnSparseD("hpeak_KALFD_mix","hpeak_KALFD_mix",3,bins3D,xmin3D,xmax3D);



TH1D *V0_chi2_SS_LLFD_XI=new TH1D("V0_chi2_SS_LLFD_XI","",10000,-0.01,0.11);
TH1D *V0_chi2_OS_LLFD_XI=new TH1D("V0_chi2_OS_LLFD_XI","",10000,-0.01,0.11);
TH1D *V0_chi2_SS_ALFD_XI=new TH1D("V0_chi2_SS_ALFD_XI","",10000,-0.01,0.11);
TH1D *V0_chi2_OS_ALFD_XI=new TH1D("V0_chi2_OS_ALFD_XI","",10000,-0.01,0.11);
TH1D *V0_chi2_SS_KLFD_XI=new TH1D("V0_chi2_SS_KLFD_XI","",10000,-0.01,0.11);
TH1D *V0_chi2_OS_KLFD_XI=new TH1D("V0_chi2_OS_KLFD_XI","",10000,-0.01,0.11);

TH1D *V0_chi2_SS_LLFD_AXI=new TH1D("V0_chi2_SS_LLFD_AXI","",10000,-0.01,0.11);
TH1D *V0_chi2_OS_LLFD_AXI=new TH1D("V0_chi2_OS_LLFD_AXI","",10000,-0.01,0.11);
TH1D *V0_chi2_SS_ALFD_AXI=new TH1D("V0_chi2_SS_ALFD_AXI","",10000,-0.01,0.11);
TH1D *V0_chi2_OS_ALFD_AXI=new TH1D("V0_chi2_OS_ALFD_AXI","",10000,-0.01,0.11);
TH1D *V0_chi2_SS_KLFD_AXI=new TH1D("V0_chi2_SS_KLFD_AXI","",10000,-0.01,0.11);
TH1D *V0_chi2_OS_KLFD_AXI=new TH1D("V0_chi2_OS_KLFD_AXI","",10000,-0.01,0.11);

TH1D *V0_chi2_SS_LLFD_OM=new TH1D("V0_chi2_SS_LLFD_OM","",10000,-0.01,0.11);
TH1D *V0_chi2_OS_LLFD_OM=new TH1D("V0_chi2_OS_LLFD_OM","",10000,-0.01,0.11);
TH1D *V0_chi2_SS_ALFD_OM=new TH1D("V0_chi2_SS_ALFD_OM","",10000,-0.01,0.11);
TH1D *V0_chi2_OS_ALFD_OM=new TH1D("V0_chi2_OS_ALFD_OM","",10000,-0.01,0.11);
TH1D *V0_chi2_SS_KLFD_OM=new TH1D("V0_chi2_SS_KLFD_OM","",10000,-0.01,0.11);
TH1D *V0_chi2_OS_KLFD_OM=new TH1D("V0_chi2_OS_KLFD_OM","",10000,-0.01,0.11);

TH1D *V0_chi2_SS_LLFD_AOM=new TH1D("V0_chi2_SS_LLFD_AOM","",10000,-0.01,0.11);
TH1D *V0_chi2_OS_LLFD_AOM=new TH1D("V0_chi2_OS_LLFD_AOM","",10000,-0.01,0.11);
TH1D *V0_chi2_SS_ALFD_AOM=new TH1D("V0_chi2_SS_ALFD_AOM","",10000,-0.01,0.11);
TH1D *V0_chi2_OS_ALFD_AOM=new TH1D("V0_chi2_OS_ALFD_AOM","",10000,-0.01,0.11);
TH1D *V0_chi2_SS_KLFD_AOM=new TH1D("V0_chi2_SS_KLFD_AOM","",10000,-0.01,0.11);
TH1D *V0_chi2_OS_KLFD_AOM=new TH1D("V0_chi2_OS_KLFD_AOM","",10000,-0.01,0.11);


//jets

  THnSparseD *h_jet_jet  = new THnSparseD("h_jet_jet","h_jet_jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_jet_jet  = new THnSparseD("h_inv_jet_jet","h_inv_jet_jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_jet_jet  = new THnSparseD("h_rot_jet_jet","h_rot_jet_jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_jet_jet  = new THnSparseD("hMix_jet_jet","hMix_jet_jet",3,bins3D,xmin3D,xmax3D);

//K0sK0s

  THnSparseD *h_K0s_K0s_samejet  = new THnSparseD("h_K0s_K0s_samejet","h_K0s_K0s_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_K0s_K0s_samejet  = new THnSparseD("h_inv_K0s_K0s_samejet","h_inv_K0s_K0s_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_K0s_K0s_samejet  = new THnSparseD("h_rot_K0s_K0s_samejet","h_rot_K0s_K0s_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_K0s_samejet  = new THnSparseD("hMix_K0s_K0s_samejet","hMix_K0s_K0s_samejet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_K0s_K0s_diffjet  = new THnSparseD("h_K0s_K0s_diffjet","h_K0s_K0s_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_K0s_K0s_diffjet  = new THnSparseD("h_inv_K0s_K0s_diffjet","h_inv_K0s_K0s_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_K0s_K0s_diffjet  = new THnSparseD("h_rot_K0s_K0s_diffjet","h_rot_K0s_K0s_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_K0s_diffjet  = new THnSparseD("hMix_K0s_K0s_diffjet","hMix_K0s_K0s_diffjet",3,bins3D,xmin3D,xmax3D);
  
  THnSparseD *h_K0s_K0s_nojet  = new THnSparseD("h_K0s_K0s_nojet","h_K0s_K0s_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_K0s_K0s_nojet  = new THnSparseD("h_inv_K0s_K0s_nojet","h_inv_K0s_K0s_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_K0s_K0s_nojet  = new THnSparseD("h_rot_K0s_K0s_nojet","h_rot_K0s_K0s_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_K0s_nojet  = new THnSparseD("hMix_K0s_K0s_nojet","hMix_K0s_K0s_nojet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_K0s_K0s_only1jet  = new THnSparseD("h_K0s_K0s_only1jet","h_K0s_K0s_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_K0s_K0s_only1jet  = new THnSparseD("h_inv_K0s_K0s_only1jet","h_inv_K0s_K0s_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_K0s_K0s_only1jet  = new THnSparseD("h_rot_K0s_K0s_only1jet","h_rot_K0s_K0s_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_K0s_only1jet  = new THnSparseD("hMix_K0s_K0s_only1jet","hMix_K0s_K0s_only1jet",3,bins3D,xmin3D,xmax3D);

//LamLam

  THnSparseD *h_Lam_Lam_samejet  = new THnSparseD("h_Lam_Lam_samejet","h_Lam_Lam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_Lam_Lam_samejet  = new THnSparseD("h_inv_Lam_Lam_samejet","h_inv_Lam_Lam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_Lam_Lam_samejet  = new THnSparseD("h_rot_Lam_Lam_samejet","h_rot_Lam_Lam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_Lam_Lam_samejet  = new THnSparseD("hMix_Lam_Lam_samejet","hMix_Lam_Lam_samejet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_Lam_Lam_diffjet  = new THnSparseD("h_Lam_Lam_diffjet","h_Lam_Lam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_Lam_Lam_diffjet  = new THnSparseD("h_inv_Lam_Lam_diffjet","h_inv_Lam_Lam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_Lam_Lam_diffjet  = new THnSparseD("h_rot_Lam_Lam_diffjet","h_rot_Lam_Lam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_Lam_Lam_diffjet  = new THnSparseD("hMix_Lam_Lam_diffjet","hMix_Lam_Lam_diffjet",3,bins3D,xmin3D,xmax3D);
  
  THnSparseD *h_Lam_Lam_nojet  = new THnSparseD("h_Lam_Lam_nojet","h_Lam_Lam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_Lam_Lam_nojet  = new THnSparseD("h_inv_Lam_Lam_nojet","h_inv_Lam_Lam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_Lam_Lam_nojet  = new THnSparseD("h_rot_Lam_Lam_nojet","h_rot_Lam_Lam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_Lam_Lam_nojet  = new THnSparseD("hMix_Lam_Lam_nojet","hMix_Lam_Lam_nojet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_Lam_Lam_only1jet  = new THnSparseD("h_Lam_Lam_only1jet","h_Lam_Lam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_Lam_Lam_only1jet  = new THnSparseD("h_inv_Lam_Lam_only1jet","h_inv_Lam_Lam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_Lam_Lam_only1jet  = new THnSparseD("h_rot_Lam_Lam_only1jet","h_rot_Lam_Lam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_Lam_Lam_only1jet  = new THnSparseD("hMix_Lam_Lam_only1jet","hMix_Lam_Lam_only1jet",3,bins3D,xmin3D,xmax3D);

//ALamALam

  THnSparseD *h_ALam_ALam_samejet  = new THnSparseD("h_ALam_ALam_samejet","h_ALam_ALam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_ALam_ALam_samejet  = new THnSparseD("h_inv_ALam_ALam_samejet","h_inv_ALam_ALam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_ALam_ALam_samejet  = new THnSparseD("h_rot_ALam_ALam_samejet","h_rot_ALam_ALam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_ALam_ALam_samejet  = new THnSparseD("hMix_ALam_ALam_samejet","hMix_ALam_ALam_samejet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_ALam_ALam_diffjet  = new THnSparseD("h_ALam_ALam_diffjet","h_ALam_ALam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_ALam_ALam_diffjet  = new THnSparseD("h_inv_ALam_ALam_diffjet","h_inv_ALam_ALam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_ALam_ALam_diffjet  = new THnSparseD("h_rot_ALam_ALam_diffjet","h_rot_ALam_ALam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_ALam_ALam_diffjet  = new THnSparseD("hMix_ALam_ALam_diffjet","hMix_ALam_ALam_diffjet",3,bins3D,xmin3D,xmax3D);
  
  THnSparseD *h_ALam_ALam_nojet  = new THnSparseD("h_ALam_ALam_nojet","h_ALam_ALam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_ALam_ALam_nojet  = new THnSparseD("h_inv_ALam_ALam_nojet","h_inv_ALam_ALam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_ALam_ALam_nojet  = new THnSparseD("h_rot_ALam_ALam_nojet","h_rot_ALam_ALam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_ALam_ALam_nojet  = new THnSparseD("hMix_ALam_ALam_nojet","hMix_ALam_ALam_nojet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_ALam_ALam_only1jet  = new THnSparseD("h_ALam_ALam_only1jet","h_ALam_ALam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_ALam_ALam_only1jet  = new THnSparseD("h_inv_ALam_ALam_only1jet","h_inv_ALam_ALam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_ALam_ALam_only1jet  = new THnSparseD("h_rot_ALam_ALam_only1jet","h_rot_ALam_ALam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_ALam_ALam_only1jet  = new THnSparseD("hMix_ALam_ALam_only1jet","hMix_ALam_ALam_only1jet",3,bins3D,xmin3D,xmax3D);
  
//LAL

  THnSparseD *h_Lam_ALam_samejet  = new THnSparseD("h_Lam_ALam_samejet","h_Lam_ALam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_Lam_ALam_samejet  = new THnSparseD("h_inv_Lam_ALam_samejet","h_inv_Lam_ALam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_Lam_ALam_samejet  = new THnSparseD("h_rot_Lam_ALam_samejet","h_rot_Lam_ALam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_Lam_ALam_samejet  = new THnSparseD("hMix_Lam_ALam_samejet","hMix_Lam_ALam_samejet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_Lam_ALam_diffjet  = new THnSparseD("h_Lam_ALam_diffjet","h_Lam_ALam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_Lam_ALam_diffjet  = new THnSparseD("h_inv_Lam_ALam_diffjet","h_inv_Lam_ALam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_Lam_ALam_diffjet  = new THnSparseD("h_rot_Lam_ALam_diffjet","h_rot_Lam_ALam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_Lam_ALam_diffjet  = new THnSparseD("hMix_Lam_ALam_diffjet","hMix_Lam_ALam_diffjet",3,bins3D,xmin3D,xmax3D);
  
  THnSparseD *h_Lam_ALam_nojet  = new THnSparseD("h_Lam_ALam_nojet","h_Lam_ALam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_Lam_ALam_nojet  = new THnSparseD("h_inv_Lam_ALam_nojet","h_inv_Lam_ALam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_Lam_ALam_nojet  = new THnSparseD("h_rot_Lam_ALam_nojet","h_rot_Lam_ALam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_Lam_ALam_nojet  = new THnSparseD("hMix_Lam_ALam_nojet","hMix_Lam_ALam_nojet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_Lam_ALam_only1jet  = new THnSparseD("h_Lam_ALam_only1jet","h_Lam_ALam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_Lam_ALam_only1jet  = new THnSparseD("h_inv_Lam_ALam_only1jet","h_inv_Lam_ALam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_Lam_ALam_only1jet  = new THnSparseD("h_rot_Lam_ALam_only1jet","h_rot_Lam_ALam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_Lam_ALam_only1jet  = new THnSparseD("hMix_Lam_ALam_only1jet","hMix_Lam_ALam_only1jet",3,bins3D,xmin3D,xmax3D);
  
//KL

  THnSparseD *h_K0s_Lam_samejet  = new THnSparseD("h_K0s_Lam_samejet","h_K0s_Lam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_K0s_Lam_samejet  = new THnSparseD("h_inv_K0s_Lam_samejet","h_inv_K0s_Lam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_K0s_Lam_samejet  = new THnSparseD("h_rot_K0s_Lam_samejet","h_rot_K0s_Lam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_Lam_samejet  = new THnSparseD("hMix_K0s_Lam_samejet","hMix_K0s_Lam_samejet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_K0s_Lam_diffjet  = new THnSparseD("h_K0s_Lam_diffjet","h_K0s_Lam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_K0s_Lam_diffjet  = new THnSparseD("h_inv_K0s_Lam_diffjet","h_inv_K0s_Lam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_K0s_Lam_diffjet  = new THnSparseD("h_rot_K0s_Lam_diffjet","h_rot_K0s_Lam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_Lam_diffjet  = new THnSparseD("hMix_K0s_Lam_diffjet","hMix_K0s_Lam_diffjet",3,bins3D,xmin3D,xmax3D);
  
  THnSparseD *h_K0s_Lam_nojet  = new THnSparseD("h_K0s_Lam_nojet","h_K0s_Lam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_K0s_Lam_nojet  = new THnSparseD("h_inv_K0s_Lam_nojet","h_inv_K0s_Lam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_K0s_Lam_nojet  = new THnSparseD("h_rot_K0s_Lam_nojet","h_rot_K0s_Lam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_Lam_nojet  = new THnSparseD("hMix_K0s_Lam_nojet","hMix_K0s_Lam_nojet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_K0s_Lam_only1jet  = new THnSparseD("h_K0s_Lam_only1jet","h_K0s_Lam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_K0s_Lam_only1jet  = new THnSparseD("h_inv_K0s_Lam_only1jet","h_inv_K0s_Lam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_K0s_Lam_only1jet  = new THnSparseD("h_rot_K0s_Lam_only1jet","h_rot_K0s_Lam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_Lam_only1jet  = new THnSparseD("hMix_K0s_Lam_only1jet","hMix_K0s_Lam_only1jet",3,bins3D,xmin3D,xmax3D);

//KAL

  THnSparseD *h_K0s_ALam_samejet  = new THnSparseD("h_K0s_ALam_samejet","h_K0s_ALam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_K0s_ALam_samejet  = new THnSparseD("h_inv_K0s_ALam_samejet","h_inv_K0s_ALam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_K0s_ALam_samejet  = new THnSparseD("h_rot_K0s_ALam_samejet","h_rot_K0s_ALam_samejet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_ALam_samejet  = new THnSparseD("hMix_K0s_ALam_samejet","hMix_K0s_ALam_samejet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_K0s_ALam_diffjet  = new THnSparseD("h_K0s_ALam_diffjet","h_K0s_ALam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_K0s_ALam_diffjet  = new THnSparseD("h_inv_K0s_ALam_diffjet","h_inv_K0s_ALam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_K0s_ALam_diffjet  = new THnSparseD("h_rot_K0s_ALam_diffjet","h_rot_K0s_ALam_diffjet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_ALam_diffjet  = new THnSparseD("hMix_K0s_ALam_diffjet","hMix_K0s_ALam_diffjet",3,bins3D,xmin3D,xmax3D);
  
  THnSparseD *h_K0s_ALam_nojet  = new THnSparseD("h_K0s_ALam_nojet","h_K0s_ALam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_K0s_ALam_nojet  = new THnSparseD("h_inv_K0s_ALam_nojet","h_inv_K0s_ALam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_K0s_ALam_nojet  = new THnSparseD("h_rot_K0s_ALam_nojet","h_rot_K0s_ALam_nojet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_ALam_nojet  = new THnSparseD("hMix_K0s_ALam_nojet","hMix_K0s_ALam_nojet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_K0s_ALam_only1jet  = new THnSparseD("h_K0s_ALam_only1jet","h_K0s_ALam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_inv_K0s_ALam_only1jet  = new THnSparseD("h_inv_K0s_ALam_only1jet","h_inv_K0s_ALam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_rot_K0s_ALam_only1jet  = new THnSparseD("h_rot_K0s_ALam_only1jet","h_rot_K0s_ALam_only1jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_ALam_only1jet  = new THnSparseD("hMix_K0s_ALam_only1jet","hMix_K0s_ALam_only1jet",3,bins3D,xmin3D,xmax3D);
  
  
// hist for V0 + jet  

  THnSparseD *h_K0snojet_Jet  = new THnSparseD("h_K0snojet_Jet","h_K0snojet_Jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_K0snojet_Jet_rot  = new THnSparseD("h_K0snojet_Jet_rot","h_K0snojet_Jet_rot",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_K0snojet_Jet_inv  = new THnSparseD("h_K0snojet_Jet_inv","h_K0snojet_Jet_inv",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0snojet_Jet  = new THnSparseD("hMix_K0snojet_Jet","hMix_K0snojet_Jet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_Lamnojet_Jet  = new THnSparseD("h_Lamnojet_Jet","h_Lamnojet_Jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_Lamnojet_Jet_rot  = new THnSparseD("h_Lamnojet_Jet_rot","h_Lamnojet_Jet_rot",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_Lamnojet_Jet_inv  = new THnSparseD("h_Lamnojet_Jet_inv","h_Lamnojet_Jet_inv",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_Lamnojet_Jet  = new THnSparseD("hMix_Lamnojet_Jet","hMix_Lamnojet_Jet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_ALamnojet_Jet  = new THnSparseD("h_ALamnojet_Jet","h_ALamnojet_Jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_ALamnojet_Jet_rot  = new THnSparseD("h_ALamnojet_Jet_rot","h_ALamnojet_Jet_rot",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_ALamnojet_Jet_inv  = new THnSparseD("h_ALamnojet_Jet_inv","h_ALamnojet_Jet_inv",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_ALamnojet_Jet  = new THnSparseD("hMix_ALamnojet_Jet","hMix_ALamnojet_Jet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_K0sfromjet_Jet  = new THnSparseD("h_K0sfromjet_Jet","h_K0sfromjet_Jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_K0sfromjet_Jet_rot  = new THnSparseD("h_K0sfromjet_Jet_rot","h_K0sfromjet_Jet_rot",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_K0sfromjet_Jet_inv  = new THnSparseD("h_K0sfromjet_Jet_inv","h_K0sfromjet_Jet_inv",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0sfromjet_Jet  = new THnSparseD("hMix_K0sfromjet_Jet","hMix_K0sfromjet_Jet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_Lamfromjet_Jet  = new THnSparseD("h_Lamfromjet_Jet","h_Lamfromjet_Jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_Lamfromjet_Jet_rot  = new THnSparseD("h_Lamfromjet_Jet_rot","h_Lamfromjet_Jet_rot",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_Lamfromjet_Jet_inv  = new THnSparseD("h_Lamfromjet_Jet_inv","h_Lamfromjet_Jet_inv",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_Lamfromjet_Jet  = new THnSparseD("hMix_Lamfromjet_Jet","hMix_Lamfromjet_Jet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_ALamfromjet_Jet  = new THnSparseD("h_ALamfromjet_Jet","h_ALamfromjet_Jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_ALamfromjet_Jet_rot  = new THnSparseD("h_ALamfromjet_Jet_rot","h_ALamfromjet_Jet_rot",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_ALamfromjet_Jet_inv  = new THnSparseD("h_ALamfromjet_Jet_inv","h_ALamfromjet_Jet_inv",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_ALamfromjet_Jet  = new THnSparseD("hMix_ALamfromjet_Jet","hMix_ALamfromjet_Jet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_K0s_Jet  = new THnSparseD("h_K0s_Jet","h_K0s_Jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_K0s_Jet_rot  = new THnSparseD("h_K0s_Jet_rot","h_K0s_Jet_rot",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_K0s_Jet_inv  = new THnSparseD("h_K0s_Jet_inv","h_K0s_Jet_inv",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_K0s_Jet  = new THnSparseD("hMix_K0s_Jet","hMix_K0s_Jet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_Lam_Jet  = new THnSparseD("h_Lam_Jet","h_Lam_Jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_Lam_Jet_rot  = new THnSparseD("h_Lam_Jet_rot","h_Lam_Jet_rot",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_Lam_Jet_inv  = new THnSparseD("h_Lam_Jet_inv","h_Lam_Jet_inv",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_Lam_Jet  = new THnSparseD("hMix_Lam_Jet","hMix_Lam_Jet",3,bins3D,xmin3D,xmax3D);

  THnSparseD *h_ALam_Jet  = new THnSparseD("h_ALam_Jet","h_ALam_Jet",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_ALam_Jet_rot  = new THnSparseD("h_ALam_Jet_rot","h_ALam_Jet_rot",3,bins3D,xmin3D,xmax3D);
  THnSparseD *h_ALam_Jet_inv  = new THnSparseD("h_ALam_Jet_inv","h_ALam_Jet_inv",3,bins3D,xmin3D,xmax3D);
  THnSparseD *hMix_ALam_Jet  = new THnSparseD("hMix_ALam_Jet","hMix_ALam_Jet",3,bins3D,xmin3D,xmax3D);
  
  
  //HDibaryon
  Int_t bins3D_HDibaryon[3]=   {100,75,25};
  Double_t xmin3D_HDibaryon[3]={2.2,0.,0.};
  Double_t xmax3D_HDibaryon[3]={2.4,15.0,500.};

  THnSparseD *h_HDibaryon_Lam_Lam  = new THnSparseD("h_HDibaryon_Lam_Lam","h_HDibaryon_Lam_Lam",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_rot_Lam_Lam  = new THnSparseD("h_HDibaryon_rot_Lam_Lam","h_HDibaryon_rot_Lam_Lam",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_inv_Lam_Lam  = new THnSparseD("h_HDibaryon_inv_Lam_Lam","h_HDibaryon_inv_Lam_Lam",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_Lam_Lam_mix  = new THnSparseD("h_HDibaryon_Lam_Lam_mix","h_HDibaryon_Lam_Lam_mix",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);

  THnSparseD *h_HDibaryon_Lam_Lam_side  = new THnSparseD("h_HDibaryon_Lam_Lam_side","h_HDibaryon_Lam_Lam_side",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_rot_Lam_Lam_side  = new THnSparseD("h_HDibaryon_rot_Lam_Lam_side","h_HDibaryon_rot_Lam_Lam_side",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_inv_Lam_Lam_side  = new THnSparseD("h_HDibaryon_inv_Lam_Lam_side","h_HDibaryon_inv_Lam_Lam_side",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_Lam_Lam_mix_side  = new THnSparseD("h_HDibaryon_Lam_Lam_mix_side","h_HDibaryon_Lam_Lam_mix_side",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);

  THnSparseD *h_HDibaryon_Lam_Lam_peakside  = new THnSparseD("h_HDibaryon_Lam_Lam_peakside","h_HDibaryon_Lam_Lam_peakside",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_rot_Lam_Lam_peakside  = new THnSparseD("h_HDibaryon_rot_Lam_Lam_peakside","h_HDibaryon_rot_Lam_Lam_peakside",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_inv_Lam_Lam_peakside  = new THnSparseD("h_HDibaryon_inv_Lam_Lam_peakside","h_HDibaryon_inv_Lam_Lam_peakside",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_Lam_Lam_mix_peakside  = new THnSparseD("h_HDibaryon_Lam_Lam_mix_peakside","h_HDibaryon_Lam_Lam_mix_peakside",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);


  THnSparseD *h_HDibaryon_ALam_ALam  = new THnSparseD("h_HDibaryon_ALam_ALam","h_HDibaryon_ALam_ALam",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_rot_ALam_ALam  = new THnSparseD("h_HDibaryon_rot_ALam_ALam","h_HDibaryon_rot_ALam_ALam",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_inv_ALam_ALam  = new THnSparseD("h_HDibaryon_inv_ALam_ALam","h_HDibaryon_inv_ALam_ALam",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_ALam_ALam_mix  = new THnSparseD("h_HDibaryon_ALam_ALam_mix","h_HDibaryon_ALam_ALam_mix",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);

  THnSparseD *h_HDibaryon_ALam_ALam_side  = new THnSparseD("h_HDibaryon_ALam_ALam_side","h_HDibaryon_ALam_ALam_side",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_rot_ALam_ALam_side  = new THnSparseD("h_HDibaryon_rot_ALam_ALam_side","h_HDibaryon_rot_ALam_ALam_side",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_inv_ALam_ALam_side  = new THnSparseD("h_HDibaryon_inv_ALam_ALam_side","h_HDibaryon_inv_ALam_ALam_side",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_ALam_ALam_mix_side  = new THnSparseD("h_HDibaryon_ALam_ALam_mix_side","h_HDibaryon_ALam_ALam_mix_side",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);

  THnSparseD *h_HDibaryon_ALam_ALam_peakside  = new THnSparseD("h_HDibaryon_ALam_ALam_peakside","h_HDibaryon_ALam_ALam_peakside",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_rot_ALam_ALam_peakside  = new THnSparseD("h_HDibaryon_rot_ALam_ALam_peakside","h_HDibaryon_rot_ALam_ALam_peakside",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_inv_ALam_ALam_peakside  = new THnSparseD("h_HDibaryon_inv_ALam_ALam_peakside","h_HDibaryon_inv_ALam_ALam_peakside",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_ALam_ALam_mix_peakside  = new THnSparseD("h_HDibaryon_ALam_ALam_mix_peakside","h_HDibaryon_ALam_ALam_mix_peakside",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);


  THnSparseD *h_HDibaryon_Lam_ALam  = new THnSparseD("h_HDibaryon_Lam_ALam","h_HDibaryon_Lam_ALam",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_Lam_ALam_side  = new THnSparseD("h_HDibaryon_Lam_ALam_side","h_HDibaryon_Lam_ALam_side",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);
  THnSparseD *h_HDibaryon_Lam_ALam_peakside  = new THnSparseD("h_HDibaryon_Lam_ALam_peakside","h_HDibaryon_Lam_ALam_peakside",3,bins3D_HDibaryon,xmin3D_HDibaryon,xmax3D_HDibaryon);

  //Hf2
  Int_t bins3D_Hf2[3]=   {200,75,25};
  Double_t xmin3D_Hf2[3]={1.3,0.,0.};
  Double_t xmax3D_Hf2[3]={1.7,15.0,500.};

  THnSparseD *h_Hf2_K0s_K0s  = new THnSparseD("h_Hf2_K0s_K0s","h_Hf2_K0s_K0s",3,bins3D_Hf2,xmin3D_Hf2,xmax3D_Hf2);
  THnSparseD *h_Hf2_rot_K0s_K0s  = new THnSparseD("h_Hf2_rot_K0s_K0s","h_Hf2_rot_K0s_K0s",3,bins3D_Hf2,xmin3D_Hf2,xmax3D_Hf2);
  THnSparseD *h_Hf2_inv_K0s_K0s  = new THnSparseD("h_Hf2_inv_K0s_K0s","h_Hf2_inv_K0s_K0s",3,bins3D_Hf2,xmin3D_Hf2,xmax3D_Hf2);
  THnSparseD *h_Hf2_K0s_K0s_mix  = new THnSparseD("h_Hf2_K0s_K0s_mix","h_Hf2_K0s_K0s_mix",3,bins3D_Hf2,xmin3D_Hf2,xmax3D_Hf2);

  THnSparseD *h_Hf2_K0s_K0s_side  = new THnSparseD("h_Hf2_K0s_K0s_side","h_Hf2_K0s_K0s_side",3,bins3D_Hf2,xmin3D_Hf2,xmax3D_Hf2);
  THnSparseD *h_Hf2_rot_K0s_K0s_side  = new THnSparseD("h_Hf2_rot_K0s_K0s_side","h_Hf2_rot_K0s_K0s_side",3,bins3D_Hf2,xmin3D_Hf2,xmax3D_Hf2);
  THnSparseD *h_Hf2_inv_K0s_K0s_side  = new THnSparseD("h_Hf2_inv_K0s_K0s_side","h_Hf2_inv_K0s_K0s_side",3,bins3D_Hf2,xmin3D_Hf2,xmax3D_Hf2);
  THnSparseD *h_Hf2_K0s_K0s_mix_side  = new THnSparseD("h_Hf2_K0s_K0s_mix_side","h_Hf2_K0s_K0s_mix_side",3,bins3D_Hf2,xmin3D_Hf2,xmax3D_Hf2);

  THnSparseD *h_Hf2_K0s_K0s_peakside  = new THnSparseD("h_Hf2_K0s_K0s_peakside","h_Hf2_K0s_K0s_peakside",3,bins3D_Hf2,xmin3D_Hf2,xmax3D_Hf2);
  THnSparseD *h_Hf2_rot_K0s_K0s_peakside  = new THnSparseD("h_Hf2_rot_K0s_K0s_peakside","h_Hf2_rot_K0s_K0s_peakside",3,bins3D_Hf2,xmin3D_Hf2,xmax3D_Hf2);
  THnSparseD *h_Hf2_inv_K0s_K0s_peakside  = new THnSparseD("h_Hf2_inv_K0s_K0s_peakside","h_Hf2_inv_K0s_K0s_peakside",3,bins3D_Hf2,xmin3D_Hf2,xmax3D_Hf2);
  THnSparseD *h_Hf2_K0s_K0s_mix_peakside  = new THnSparseD("h_Hf2_K0s_K0s_mix_peakside","h_Hf2_K0s_K0s_mix_peakside",3,bins3D_Hf2,xmin3D_Hf2,xmax3D_Hf2);

  //Hcasc1820
  Int_t bins3D_Hcasc1820[3]=   {200,75,25};
  Double_t xmin3D_Hcasc1820[3]={1.5,0.,0.};
  Double_t xmax3D_Hcasc1820[3]={1.9,15.0,500.};

  THnSparseD *h_Hcasc1820_K0s_Lam  = new THnSparseD("h_Hcasc1820_K0s_Lam","h_Hcasc1820_K0s_Lam",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_rot_K0s_Lam  = new THnSparseD("h_Hcasc1820_rot_K0s_Lam","h_Hcasc1820_rot_K0s_Lam",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_inv_K0s_Lam  = new THnSparseD("h_Hcasc1820_inv_K0s_Lam","h_Hcasc1820_inv_K0s_Lam",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_K0s_Lam_mix  = new THnSparseD("h_Hcasc1820_K0s_Lam_mix","h_Hcasc1820_K0s_Lam_mix",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);

  THnSparseD *h_Hcasc1820_K0s_Lam_side  = new THnSparseD("h_Hcasc1820_K0s_Lam_side","h_Hcasc1820_K0s_Lam_side",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_rot_K0s_Lam_side  = new THnSparseD("h_Hcasc1820_rot_K0s_Lam_side","h_Hcasc1820_rot_K0s_Lam_side",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_inv_K0s_Lam_side  = new THnSparseD("h_Hcasc1820_inv_K0s_Lam_side","h_Hcasc1820_inv_K0s_Lam_side",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_K0s_Lam_mix_side  = new THnSparseD("h_Hcasc1820_K0s_Lam_mix_side","h_Hcasc1820_K0s_Lam_mix_side",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);

  THnSparseD *h_Hcasc1820_K0s_Lam_peakside  = new THnSparseD("h_Hcasc1820_K0s_Lam_peakside","h_Hcasc1820_K0s_Lam_peakside",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_rot_K0s_Lam_peakside  = new THnSparseD("h_Hcasc1820_rot_K0s_Lam_peakside","h_Hcasc1820_rot_K0s_Lam_peakside",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_inv_K0s_Lam_peakside  = new THnSparseD("h_Hcasc1820_inv_K0s_Lam_peakside","h_Hcasc1820_inv_K0s_Lam_peakside",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_K0s_Lam_mix_peakside  = new THnSparseD("h_Hcasc1820_K0s_Lam_mix_peakside","h_Hcasc1820_K0s_Lam_mix_peakside",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);

  THnSparseD *h_Hcasc1820_K0s_ALam  = new THnSparseD("h_Hcasc1820_K0s_ALam","h_Hcasc1820_K0s_ALam",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_rot_K0s_ALam  = new THnSparseD("h_Hcasc1820_rot_K0s_ALam","h_Hcasc1820_rot_K0s_ALam",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_inv_K0s_ALam  = new THnSparseD("h_Hcasc1820_inv_K0s_ALam","h_Hcasc1820_inv_K0s_ALam",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_K0s_ALam_mix  = new THnSparseD("h_Hcasc1820_K0s_ALam_mix","h_Hcasc1820_K0s_ALam_mix",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);

  THnSparseD *h_Hcasc1820_K0s_ALam_side  = new THnSparseD("h_Hcasc1820_K0s_ALam_side","h_Hcasc1820_K0s_ALam_side",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_rot_K0s_ALam_side  = new THnSparseD("h_Hcasc1820_rot_K0s_ALam_side","h_Hcasc1820_rot_K0s_ALam_side",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_inv_K0s_ALam_side  = new THnSparseD("h_Hcasc1820_inv_K0s_ALam_side","h_Hcasc1820_inv_K0s_ALam_side",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_K0s_ALam_mix_side  = new THnSparseD("h_Hcasc1820_K0s_ALam_mix_side","h_Hcasc1820_K0s_ALam_mix_side",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);

  THnSparseD *h_Hcasc1820_K0s_ALam_peakside  = new THnSparseD("h_Hcasc1820_K0s_ALam_peakside","h_Hcasc1820_K0s_ALam_peakside",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_rot_K0s_ALam_peakside  = new THnSparseD("h_Hcasc1820_rot_K0s_ALam_peakside","h_Hcasc1820_rot_K0s_ALam_peakside",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_inv_K0s_ALam_peakside  = new THnSparseD("h_Hcasc1820_inv_K0s_ALam_peakside","h_Hcasc1820_inv_K0s_ALam_peakside",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  THnSparseD *h_Hcasc1820_K0s_ALam_mix_peakside  = new THnSparseD("h_Hcasc1820_K0s_ALam_mix_peakside","h_Hcasc1820_K0s_ALam_mix_peakside",3,bins3D_Hcasc1820,xmin3D_Hcasc1820,xmax3D_Hcasc1820);
  
 TH2D *Ntrk_VZ=new TH2D("Ntrk_VZ","",400,-20.,20.,500,0.,500.);
 
//sizes
 TH1D *K_size=new TH1D("K_size","",100,0.,100.);
 TH1D *L_size=new TH1D("L_size","",100,0.,100.);
 TH1D *AL_size=new TH1D("AL_size","",100,0.,100.);

 TH1D *K_size_nocut=new TH1D("K_size_nocut","",100,0.,100.);
 TH1D *L_size_nocut=new TH1D("L_size_nocut","",100,0.,100.);
 TH1D *AL_size_nocut=new TH1D("AL_size_nocut","",100,0.,100.);
 
// Armeteros Podolanski 
 TH2D *APwocut_K0s=new TH2D("APwocut_K0s","APwocut_K0s",400,-1.,1.,400,0.,1.);
 TH2D *APwcut_K0s=new TH2D("APwcut_K0s","APwcut_K0s",400,-1.,1.,400,0.,1.);
 TH2D *APwocut_LAL=new TH2D("APwocut_LAL","APwocut_LAL",400,-1.,1.,400,0.,1.);
 TH2D *APwcut_LAL=new TH2D("APwcut_LAL","APwcut_LAL",400,-1.,1.,400,0.,1.);
 
 TH1D *K0smiss_mass_hyp=new TH1D("K0smiss_hyp","",400,1.0,1.2);
 TH1D *K0smiss_mass_ee=new TH1D("K0smissee","",500,0.,1.);
 TH1D *LALmiss_mass_hyp=new TH1D("LALmiss_mass_hyp","",400,0.4,0.6);
 TH1D *LALmiss_mass_ee=new TH1D("LALmiss_mass_ee","",500,0.,1.);
 
//  control histograms (single particles)

//pt vs mass
 TH2D *pt_mass_K=new TH2D("pt_mass_K","",100,0.0,10.0,400,0.4,0.6);
 TH2D *pt_mass_L=new TH2D("pt_mass_L","",100,0.0,10.0,400,1.0,1.2);
 TH2D *pt_mass_AL=new TH2D("pt_mass_AL","",100,0.0,10.0,400,1.0,1.2);

//data and MC
TH1D *K0s_pt_reco=new TH1D("K0s_pt_reco","K0s_pt_reco",100,0.0,10.);
TH1D *K0s_eta_reco=new TH1D("K0s_eta_reco","K0s_eta_reco",50,-2.5,2.5);
TH1D *K0s_phi_reco=new TH1D("K0s_phi_reco","K0s_phi_reco",64,-3.2,3.2);
TH1D *K0s_mass_reco=new TH1D("K0s_mass_reco","K0s_mass_reco",400,0.4,0.6);

TH1D *L_pt_reco=new TH1D("L_pt_reco","L_pt_reco",100,0.0,10.);
TH1D *L_eta_reco=new TH1D("L_eta_reco","L_eta_reco",50,-2.5,2.5);
TH1D *L_phi_reco=new TH1D("L_phi_reco","L_phi_reco",64,-3.2,3.2);
TH1D *L_mass_reco=new TH1D("L_mass_reco","L_mass_reco",400,1.0,1.2);

TH1D *AL_pt_reco=new TH1D("AL_pt_reco","AL_pt_reco",100,0.0,10.);
TH1D *AL_eta_reco=new TH1D("AL_eta_reco","AL_eta_reco",50,-2.5,2.5);
TH1D *AL_phi_reco=new TH1D("AL_phi_reco","AL_phi_reco",64,-3.2,3.2);
TH1D *AL_mass_reco=new TH1D("AL_mass_reco","AL_mass_reco",400,1.0,1.2);

//daughter 1
TH1D *K0s_pt_reco_d1=new TH1D("K0s_pt_reco_d1","K0s_pt_reco_d1",100,0.0,10.);
TH1D *K0s_eta_reco_d1=new TH1D("K0s_eta_reco_d1","K0s_eta_reco_d1",50,-2.5,2.5);
TH1D *K0s_phi_reco_d1=new TH1D("K0s_phi_reco_d1","K0s_phi_reco_d1",64,-3.2,3.2);
TH1D *K0s_mass_reco_d1=new TH1D("K0s_mass_reco_d1","K0s_mass_reco_d1",500,0.135,0.140);

TH1D *L_pt_reco_d1=new TH1D("L_pt_reco_d1","L_pt_reco_d1",100,0.0,10.);
TH1D *L_eta_reco_d1=new TH1D("L_eta_reco_d1","L_eta_reco_d1",50,-2.5,2.5);
TH1D *L_phi_reco_d1=new TH1D("L_phi_reco_d1","L_phi_reco_d1",64,-3.2,3.2);
TH1D *L_mass_reco_d1=new TH1D("L_mass_reco_d1","L_mass_reco_d1",500,0.935,0.940);

TH1D *AL_pt_reco_d1=new TH1D("AL_pt_reco_d1","AL_pt_reco_d1",100,0.0,10.);
TH1D *AL_eta_reco_d1=new TH1D("AL_eta_reco_d1","AL_eta_reco_d1",50,-2.5,2.5);
TH1D *AL_phi_reco_d1=new TH1D("AL_phi_reco_d1","AL_phi_reco_d1",64,-3.2,3.2);
TH1D *AL_mass_reco_d1=new TH1D("AL_mass_reco_d1","AL_mass_reco_d1",500,0.935,0.940);

//daughter 2

TH1D *K0s_pt_reco_d2=new TH1D("K0s_pt_reco_d2","K0s_pt_reco_d2",100,0.0,10.);
TH1D *K0s_eta_reco_d2=new TH1D("K0s_eta_reco_d2","K0s_eta_reco_d2",50,-2.5,2.5);
TH1D *K0s_phi_reco_d2=new TH1D("K0s_phi_reco_d2","K0s_phi_reco_d2",64,-3.2,3.2);
TH1D *K0s_mass_reco_d2=new TH1D("K0s_mass_reco_d2","K0s_mass_reco_d2",500,0.135,0.140);

TH1D *L_pt_reco_d2=new TH1D("L_pt_reco_d2","L_pt_reco_d2",100,0.0,10.);
TH1D *L_eta_reco_d2=new TH1D("L_eta_reco_d2","L_eta_reco_d2",50,-2.5,2.5);
TH1D *L_phi_reco_d2=new TH1D("L_phi_reco_d2","L_phi_reco_d2",64,-3.2,3.2);
TH1D *L_mass_reco_d2=new TH1D("L_mass_reco_d2","L_mass_reco_d2",500,0.135,0.140);

TH1D *AL_pt_reco_d2=new TH1D("AL_pt_reco_d2","AL_pt_reco_d2",100,0.0,10.);
TH1D *AL_eta_reco_d2=new TH1D("AL_eta_reco_d2","AL_eta_reco_d2",50,-2.5,2.5);
TH1D *AL_phi_reco_d2=new TH1D("AL_phi_reco_d2","AL_phi_reco_d2",64,-3.2,3.2);
TH1D *AL_mass_reco_d2=new TH1D("AL_mass_reco_d2","AL_mass_reco_d2",500,0.135,0.140);

 
// MC only 
//Gen level

TH1D *K0s_pt_gen=new TH1D("K0s_pt_gen","K0s_pt_gen",100,0.0,10.);
TH1D *K0s_eta_gen=new TH1D("K0s_eta_gen","K0s_eta_gen",50,-2.5,2.5);
TH1D *K0s_phi_gen=new TH1D("K0s_phi_gen","K0s_phi_gen",64,-3.2,3.2);
TH1D *K0s_mass_gen=new TH1D("K0s_mass_gen","K0s_mass_gen",400,0.4,0.6);

TH1D *L_pt_gen=new TH1D("L_pt_gen","L_pt_gen",100,0.0,10.);
TH1D *L_eta_gen=new TH1D("L_eta_gen","L_eta_gen",50,-2.5,2.5);
TH1D *L_phi_gen=new TH1D("L_phi_gen","L_phi_gen",64,-3.2,3.2);
TH1D *L_mass_gen=new TH1D("L_mass_gen","L_mass_gen",400,0.4,0.6);

TH1D *AL_pt_gen=new TH1D("AL_pt_gen","AL_pt_gen",100,0.0,10.);
TH1D *AL_eta_gen=new TH1D("AL_eta_gen","AL_eta_gen",50,-2.5,2.5);
TH1D *AL_phi_gen=new TH1D("AL_phi_gen","AL_phi_gen",64,-3.2,3.2);
TH1D *AL_mass_gen=new TH1D("AL_mass_gen","AL_mass_gen",400,0.4,0.6);

//Match level
TH1D *K0s_pt_mat=new TH1D("K0s_pt_mat","K0s_pt_mat",100,0.0,10.);
TH1D *K0s_eta_mat=new TH1D("K0s_eta_mat","K0s_eta_mat",50,-2.5,2.5);
TH1D *K0s_phi_mat=new TH1D("K0s_phi_mat","K0s_phi_mat",64,-3.2,3.2);
TH1D *K0s_mass_mat=new TH1D("K0s_mass_mat","K0s_mass_mat",400,0.4,0.6);

TH1D *L_pt_mat=new TH1D("L_pt_mat","L_pt_mat",100,0.0,10.);
TH1D *L_eta_mat=new TH1D("L_eta_mat","L_eta_mat",50,-2.5,2.5);
TH1D *L_phi_mat=new TH1D("L_phi_mat","L_phi_mat",64,-3.2,3.2);
TH1D *L_mass_mat=new TH1D("L_mass_mat","L_mass_mat",400,1.0,1.2);

TH1D *AL_pt_mat=new TH1D("AL_pt_mat","AL_pt_mat",100,0.0,10.);
TH1D *AL_eta_mat=new TH1D("AL_eta_mat","AL_eta_mat",50,-2.5,2.5);
TH1D *AL_phi_mat=new TH1D("AL_phi_mat","AL_phi_mat",64,-3.2,3.2);
TH1D *AL_mass_mat=new TH1D("AL_mass_mat","AL_mass_mat",400,1.0,1.2);


//cascade and omega mass
TH1D *XXi_mass=new TH1D("XXi_mass","XXi_mass",400,1.2,1.4);
TH1D *AXXi_mass=new TH1D("AXXi_mass","AXXi_mass",400,1.2,1.4);
TH1D *OOm_mass=new TH1D("OOm_mass","OOm_mass",400,1.6,1.8);
TH1D *AOOm_mass=new TH1D("AOOm_mass","AOOm_mass",400,1.6,1.8);

TH1D *XXi_mass_LambMom=new TH1D("XXi_mass_LambMom","XXi_mass_LambMom",400,1.2,1.4);
TH1D *AXXi_mass_LambMom=new TH1D("AXXi_mass_LambMom","AXXi_mass_LambMom",400,1.2,1.4);
TH1D *OOm_mass_LambMom=new TH1D("OOm_mass_LambMom","OOm_mass_LambMom",400,1.6,1.8);
TH1D *AOOm_mass_LambMom=new TH1D("AOOm_mass_LambMom","AOOm_mass_LambMom",400,1.6,1.8);


void sw2(){
    
pT_Xi->Sumw2();
pT_Om->Sumw2();
pT_Lamb->Sumw2();

XXi_mass->Sumw2();
AXXi_mass->Sumw2();
OOm_mass->Sumw2();
AOOm_mass->Sumw2();

XXi_mass_LambMom->Sumw2();
AXXi_mass_LambMom->Sumw2();
OOm_mass_LambMom->Sumw2();
AOOm_mass_LambMom->Sumw2();

nocutev->Sumw2(); 

nocutPV->Sumw2();  
nocutPU->Sumw2();  
nocutHF->Sumw2(); 
nocutSC->Sumw2(); 
nocutALL->Sumw2();
nocutALLX->Sumw2();

nev_K0s_ini->Sumw2(); 
nev_Lam_ini->Sumw2(); 
nev_ALam_ini->Sumw2(); 

nev_K0s_AS->Sumw2(); 
nev_Lam_AS->Sumw2(); 
nev_ALam_AS->Sumw2(); 
  
//K0s

nev_K0sT->Sumw2();
nev_K0sF->Sumw2();
nev_K0sTF->Sumw2();
nev_K0s_ssT->Sumw2(); 
nev_K0s_bbT->Sumw2();
nev_K0s_sbT->Sumw2();
nev_K0s_ssT_mix->Sumw2();
nev_K0s_bbT_mix->Sumw2();
nev_K0s_sbT_mix->Sumw2();
nev_K0s_ssT_etamix->Sumw2();
nev_K0s_bbT_etamix->Sumw2();
nev_K0s_sbT_etamix->Sumw2();

K0s_d1_dxy->Sumw2();
K0s_d2_dxy->Sumw2();
K0s_d1_dz->Sumw2();
K0s_d2_dz->Sumw2();
Lam_d1_dxy->Sumw2();
Lam_d2_dxy->Sumw2();
Lam_d1_dz->Sumw2();
Lam_d2_dz->Sumw2();
ALam_d1_dxy->Sumw2();
ALam_d2_dxy->Sumw2();
ALam_d1_dz->Sumw2();
ALam_d2_dz->Sumw2();

K0s_d1_Nhits->Sumw2();
K0s_d2_Nhits->Sumw2();
Lam_d1_Nhits->Sumw2();
Lam_d2_Nhits->Sumw2();
ALam_d1_Nhits->Sumw2();
ALam_d2_Nhits->Sumw2();
  
//Lam

nev_LamT->Sumw2();
nev_LamF->Sumw2();
nev_LamTF->Sumw2();
nev_Lam_ssT->Sumw2();
nev_Lam_bbT->Sumw2();
nev_Lam_sbT->Sumw2();
nev_Lam_ssT_mix->Sumw2();
nev_Lam_bbT_mix->Sumw2();
nev_Lam_sbT_mix->Sumw2();
nev_Lam_ssT_etamix->Sumw2();
nev_Lam_bbT_etamix->Sumw2();
nev_Lam_sbT_etamix->Sumw2();

  
  
//ALam

nev_ALamT->Sumw2();
nev_ALamF->Sumw2();
nev_ALamTF->Sumw2();
nev_ALam_ssT->Sumw2();
nev_ALam_bbT->Sumw2();
nev_ALam_sbT->Sumw2();
nev_ALam_ssT_mix->Sumw2();
nev_ALam_bbT_mix->Sumw2();
nev_ALam_sbT_mix->Sumw2();
nev_ALam_ssT_etamix->Sumw2();
nev_ALam_bbT_etamix->Sumw2();
nev_ALam_sbT_etamix->Sumw2();

  
  
//LAL
nev_LALT->Sumw2();
nev_LALF->Sumw2();
nev_LALTF->Sumw2();
nev_LAL_ssT->Sumw2();
nev_LAL_bbT->Sumw2();
nev_LAL_sbT->Sumw2();
nev_LAL_ssT_mix->Sumw2();
nev_LAL_bbT_mix->Sumw2();
nev_LAL_sbT_mix->Sumw2();
nev_LAL_ssT_etamix->Sumw2();
nev_LAL_bbT_etamix->Sumw2();
nev_LAL_sbT_etamix->Sumw2();

//KL

nev_KLT->Sumw2();
nev_KLF->Sumw2(); 
nev_KLTF->Sumw2();
nev_KL_ssT->Sumw2(); 
nev_KL_bbT->Sumw2();
nev_KL_sbT->Sumw2(); 
nev_KL_ssT_mix->Sumw2();
nev_KL_bbT_mix->Sumw2();
nev_KL_sbT_mix->Sumw2();
nev_KL_ssT_etamix->Sumw2();
nev_KL_bbT_etamix->Sumw2();
nev_KL_sbT_etamix->Sumw2();

//KAL

nev_KALT->Sumw2();
nev_KALF->Sumw2();
nev_KALTF->Sumw2();
nev_KAL_ssT->Sumw2();
nev_KAL_bbT->Sumw2();
nev_KAL_sbT->Sumw2();
nev_KAL_ssT_mix->Sumw2();
nev_KAL_bbT_mix->Sumw2();
nev_KAL_sbT_mix->Sumw2();
nev_KAL_ssT_etamix->Sumw2();
nev_KAL_bbT_etamix->Sumw2();
nev_KAL_sbT_etamix->Sumw2();

cone_K0s->Sumw2();
cone_K0sT->Sumw2();
cone_K0sD1->Sumw2();
cone_K0sTD1->Sumw2();
cone_K0sD2->Sumw2();
cone_K0sTD2->Sumw2();

cone_Lam->Sumw2();
cone_LamT->Sumw2();
cone_LamD1->Sumw2();
cone_LamTD1->Sumw2();
cone_LamD2->Sumw2();
cone_LamTD2->Sumw2();

cone_Jet->Sumw2();

  //vertex / ntrk / centrality  
h_vtx_z->Sumw2();  
h_vtx_rho->Sumw2();
h_ntrk_cent->Sumw2();
Ntrk_VZ->Sumw2();

h_ntrk_cent_1Ks->Sumw2();
h_ntrk_cent_1Lam->Sumw2();
h_ntrk_cent_1ALam->Sumw2();

//mass
K0s_Mass->Sumw2();
Lam_Mass->Sumw2();
ALam_Mass->Sumw2();
LAL_Mass->Sumw2();

K0sctau->Sumw2();
K0sdca3D->Sumw2();
K0sDL->Sumw2();
K0strkdca->Sumw2();
K0scostheta->Sumw2();
K0s_Vtx->Sumw2();


Lamctau->Sumw2();
Lamdca3D->Sumw2();
LamDL->Sumw2();
Lamtrkdca->Sumw2();
Lamcostheta->Sumw2();
Lam_Vtx->Sumw2();


ALamctau->Sumw2();
ALamdca3D->Sumw2();
ALamDL->Sumw2();
ALamtrkdca->Sumw2();
ALamcostheta->Sumw2();
ALam_Vtx->Sumw2();

DR_K0s->Sumw2();
DR_Lam->Sumw2();
DR_ALam->Sumw2();

DPT_K0s->Sumw2();
DPT_Lam->Sumw2();
DPT_ALam->Sumw2();

pT_Lamb_from_Xi_T->Sumw2();
pT_Lamb_from_Om_T->Sumw2();
pT_Lamb_from_prompt_T->Sumw2();

pT_Lamb_from_Xi_TF->Sumw2();
pT_Lamb_from_Om_TF->Sumw2();
pT_Lamb_from_prompt_TF->Sumw2();

pT_Lamb_daughter_Xi->Sumw2();
pT_Lamb_daughter_Om->Sumw2();
pT_Lamb_daughter_prompt->Sumw2();


hGen_K0s_K0s_Matched->Sumw2();
hGen_K0s_K0s_Matched_mix->Sumw2();

hGen_Lam_Lam_Matched->Sumw2();
hGen_Lam_Lam_Matched_mix->Sumw2();

hGen_ALam_ALam_Matched->Sumw2();
hGen_ALam_ALam_Matched_mix->Sumw2();

hGen_Lam_ALam_Matched->Sumw2();
hGen_Lam_ALam_Matched_mix->Sumw2();

hGen_K0s_Lam_Matched->Sumw2();
hGen_K0s_Lam_Matched_mix->Sumw2();

hGen_K0s_ALam_Matched->Sumw2();
hGen_K0s_ALam_Matched_mix->Sumw2();


hGen_K0s_K0s_MatchedT->Sumw2();
hGen_K0s_K0s_MatchedT_mix->Sumw2();

hGen_Lam_Lam_MatchedT->Sumw2();
hGen_Lam_Lam_MatchedT_mix->Sumw2();

hGen_ALam_ALam_MatchedT->Sumw2();
hGen_ALam_ALam_MatchedT_mix->Sumw2();

hGen_Lam_ALam_MatchedT->Sumw2();
hGen_Lam_ALam_MatchedT_mix->Sumw2();

hGen_K0s_Lam_MatchedT->Sumw2();
hGen_K0s_Lam_MatchedT_mix->Sumw2();

hGen_K0s_ALam_MatchedT->Sumw2();
hGen_K0s_ALam_MatchedT_mix->Sumw2();

//K0s
//True

hMass_K0s_K0s_T->Sumw2();   
hMass_K0s_1_T->Sumw2();   
hMass_K0s_2_T->Sumw2();   
  
hpeak_K0s_K0s_T->Sumw2();
hpeakcos_K0s_K0s_T->Sumw2();
hpeak_rot_K0s_K0s_T->Sumw2();
hpeak_inv_K0s_K0s_T->Sumw2();
hside_K0s_K0s_T->Sumw2();
hsideL_K0s_K0s_T->Sumw2();
hsideR_K0s_K0s_T->Sumw2();
hside_rot_K0s_K0s_T->Sumw2();
hside_inv_K0s_K0s_T->Sumw2();
hpeakside_K0s_K0s_T->Sumw2();
hpeaksideL_K0s_K0s_T->Sumw2();
hpeaksideR_K0s_K0s_T->Sumw2();
hpeakside_rot_K0s_K0s_T->Sumw2();
hpeakside_inv_K0s_K0s_T->Sumw2();
hpeakd1_K0s_K0s_T->Sumw2();
hpeakd2_K0s_K0s_T->Sumw2();
hpeakd12_K0s_K0s_T->Sumw2();

hpeak_K0s_K0s_T_mix->Sumw2();
hside_K0s_K0s_T_mix->Sumw2();
hsideL_K0s_K0s_T_mix->Sumw2();
hsideR_K0s_K0s_T_mix->Sumw2();
hpeakside_K0s_K0s_T_mix->Sumw2();
hpeaksideL_K0s_K0s_T_mix->Sumw2();
hpeaksideR_K0s_K0s_T_mix->Sumw2();

hpeak_K0s_K0s_T_etamix->Sumw2();
hside_K0s_K0s_T_etamix->Sumw2();
hsideL_K0s_K0s_T_etamix->Sumw2();
hsideR_K0s_K0s_T_etamix->Sumw2();
hpeakside_K0s_K0s_T_etamix->Sumw2();
hpeaksideL_K0s_K0s_T_etamix->Sumw2();
hpeaksideR_K0s_K0s_T_etamix->Sumw2();


V0chi2d1_diff_T_K0s->Sumw2();
V0chi2d2_diff_T_K0s->Sumw2();
V0chi2d12_diff_T_K0s->Sumw2();

//False

hMass_K0s_K0s_F->Sumw2();
  
hpeak_K0s_K0s_F->Sumw2();
hside_K0s_K0s_F->Sumw2();
hpeakside_K0s_K0s_F->Sumw2();
hpeakd1_K0s_K0s_F->Sumw2();
hpeakd2_K0s_K0s_F->Sumw2();
hpeakd12_K0s_K0s_F->Sumw2();
hpeak_K0s_K0s_F_bias->Sumw2();
hpeak_K0s_K0s_F_biasd1->Sumw2();
hpeak_K0s_K0s_F_biasd2->Sumw2();
hpeak_K0s_K0s_F_biasd12->Sumw2();

V0chi2d1_diff_F_K0s->Sumw2();
V0chi2d2_diff_F_K0s->Sumw2();
V0chi2d12_diff_F_K0s->Sumw2();
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
hpeak_K0s_K0s_TF->Sumw2();
hpeakd1_K0s_K0s_TF->Sumw2();
hpeakd2_K0s_K0s_TF->Sumw2();
hpeakd12_K0s_K0s_TF->Sumw2();
hside_K0s_K0s_TF->Sumw2();
hpeakside_K0s_K0s_TF->Sumw2();
hpeak_K0s_K0s_TF_bias->Sumw2();
hpeak_K0s_K0s_TF_biasd1->Sumw2();
hpeak_K0s_K0s_TF_biasd2->Sumw2();
hpeak_K0s_K0s_TF_biasd12->Sumw2();

V0chi2d1_diff_TF_K0s->Sumw2();
V0chi2d2_diff_TF_K0s->Sumw2();
V0chi2d12_diff_TF_K0s->Sumw2();

  
/////// matching hist ////// 
  
hMass_K0s_K0s_TM->Sumw2();
hMass_K0s_K0s_TU->Sumw2();
hMass_K0s_K0s_FM->Sumw2();
hMass_K0s_K0s_FU->Sumw2();
 
hpeak_K0s_K0s_TM->Sumw2();
hpeakd1_K0s_K0s_TM->Sumw2();
hpeakd2_K0s_K0s_TM->Sumw2();
hpeakd12_K0s_K0s_TM->Sumw2();
hpeak_K0s_K0s_TM_bias->Sumw2();
hpeak_K0s_K0s_TM_biasd1->Sumw2();
hpeak_K0s_K0s_TM_biasd2->Sumw2();
hpeak_K0s_K0s_TM_biasd12->Sumw2();

V0chi2d1_diff_TM_K0s->Sumw2();
V0chi2d2_diff_TM_K0s->Sumw2();
V0chi2d12_diff_TM_K0s->Sumw2();

hpeak_K0s_K0s_TU->Sumw2();
hpeakd1_K0s_K0s_TU->Sumw2();
hpeakd2_K0s_K0s_TU->Sumw2();
hpeakd12_K0s_K0s_TU->Sumw2();
hpeak_K0s_K0s_TU_bias->Sumw2();
hpeak_K0s_K0s_TU_biasd1->Sumw2();
hpeak_K0s_K0s_TU_biasd2->Sumw2();
hpeak_K0s_K0s_TU_biasd12->Sumw2();

V0chi2d1_diff_TU_K0s->Sumw2();
V0chi2d2_diff_TU_K0s->Sumw2();
V0chi2d12_diff_TU_K0s->Sumw2();

hpeak_K0s_K0s_FM->Sumw2();
hpeakd1_K0s_K0s_FM->Sumw2();
hpeakd2_K0s_K0s_FM->Sumw2();
hpeakd12_K0s_K0s_FM->Sumw2();
hpeak_K0s_K0s_FM_bias->Sumw2();
hpeak_K0s_K0s_FM_biasd1->Sumw2();
hpeak_K0s_K0s_FM_biasd2->Sumw2();
hpeak_K0s_K0s_FM_biasd12->Sumw2();

V0chi2d1_diff_FM_K0s->Sumw2();
V0chi2d2_diff_FM_K0s->Sumw2();
V0chi2d12_diff_FM_K0s->Sumw2();

hpeak_K0s_K0s_FU->Sumw2();
hpeakd1_K0s_K0s_FU->Sumw2();
hpeakd2_K0s_K0s_FU->Sumw2();
hpeakd12_K0s_K0s_FU->Sumw2();
hpeak_K0s_K0s_FU_bias->Sumw2();
hpeak_K0s_K0s_FU_biasd1->Sumw2();
hpeak_K0s_K0s_FU_biasd2->Sumw2();
hpeak_K0s_K0s_FU_biasd12->Sumw2();

V0chi2d1_diff_FU_K0s->Sumw2();
V0chi2d2_diff_FU_K0s->Sumw2();
V0chi2d12_diff_FU_K0s->Sumw2();


//true matched + true unmatched

hpeak_K0s_K0s_TM_TU->Sumw2();
hpeakd1_K0s_K0s_TM_TU->Sumw2();
hpeakd2_K0s_K0s_TM_TU->Sumw2();
hpeakd12_K0s_K0s_TM_TU->Sumw2();
hpeak_K0s_K0s_TM_TU_bias->Sumw2();
hpeak_K0s_K0s_TM_TU_biasd1->Sumw2();
hpeak_K0s_K0s_TM_TU_biasd2->Sumw2();
hpeak_K0s_K0s_TM_TU_biasd12->Sumw2();

V0chi2d1_diff_TM_TU_K0s->Sumw2();
V0chi2d2_diff_TM_TU_K0s->Sumw2();
V0chi2d12_diff_TM_TU_K0s->Sumw2();

//fake matched + fake unmatched

hpeak_K0s_K0s_FM_FU->Sumw2();
hpeakd1_K0s_K0s_FM_FU->Sumw2();
hpeakd2_K0s_K0s_FM_FU->Sumw2();
hpeakd12_K0s_K0s_FM_FU->Sumw2();
hpeak_K0s_K0s_FM_FU_bias->Sumw2();
hpeak_K0s_K0s_FM_FU_biasd1->Sumw2();
hpeak_K0s_K0s_FM_FU_biasd2->Sumw2();
hpeak_K0s_K0s_FM_FU_biasd12->Sumw2();

V0chi2d1_diff_FM_FU_K0s->Sumw2();
V0chi2d2_diff_FM_FU_K0s->Sumw2();
V0chi2d12_diff_FM_FU_K0s->Sumw2();


hpeak_K0s_K0s_TM_FM->Sumw2();
hpeakd1_K0s_K0s_TM_FM->Sumw2();
hpeakd2_K0s_K0s_TM_FM->Sumw2();
hpeakd12_K0s_K0s_TM_FM->Sumw2();
hpeak_K0s_K0s_TM_FM_bias->Sumw2();
hpeak_K0s_K0s_TM_FM_biasd1->Sumw2();
hpeak_K0s_K0s_TM_FM_biasd2->Sumw2();
hpeak_K0s_K0s_TM_FM_biasd12->Sumw2();

V0chi2d1_diff_TM_FM_K0s->Sumw2();
V0chi2d2_diff_TM_FM_K0s->Sumw2();
V0chi2d12_diff_TM_FM_K0s->Sumw2();

hpeak_K0s_K0s_TM_FU->Sumw2();
hpeakd1_K0s_K0s_TM_FU->Sumw2();
hpeakd2_K0s_K0s_TM_FU->Sumw2();
hpeakd12_K0s_K0s_TM_FU->Sumw2();
hpeak_K0s_K0s_TM_FU_bias->Sumw2();
hpeak_K0s_K0s_TM_FU_biasd1->Sumw2();
hpeak_K0s_K0s_TM_FU_biasd2->Sumw2();
hpeak_K0s_K0s_TM_FU_biasd12->Sumw2();

V0chi2d1_diff_TM_FU_K0s->Sumw2();
V0chi2d2_diff_TM_FU_K0s->Sumw2();
V0chi2d12_diff_TM_FU_K0s->Sumw2();

hpeak_K0s_K0s_TU_FU->Sumw2();
hpeakd1_K0s_K0s_TU_FU->Sumw2();
hpeakd2_K0s_K0s_TU_FU->Sumw2();
hpeakd12_K0s_K0s_TU_FU->Sumw2();
hpeak_K0s_K0s_TU_FU_bias->Sumw2();
hpeak_K0s_K0s_TU_FU_biasd1->Sumw2();
hpeak_K0s_K0s_TU_FU_biasd2->Sumw2();
hpeak_K0s_K0s_TU_FU_biasd12->Sumw2();

V0chi2d1_diff_TU_FU_K0s->Sumw2();
V0chi2d2_diff_TU_FU_K0s->Sumw2();
V0chi2d12_diff_TU_FU_K0s->Sumw2();

hpeak_K0s_K0s_TU_FM->Sumw2();
hpeakd1_K0s_K0s_TU_FM->Sumw2();
hpeakd2_K0s_K0s_TU_FM->Sumw2();
hpeakd12_K0s_K0s_TU_FM->Sumw2();
hpeak_K0s_K0s_TU_FM_bias->Sumw2();
hpeak_K0s_K0s_TU_FM_biasd1->Sumw2();
hpeak_K0s_K0s_TU_FM_biasd2->Sumw2();
hpeak_K0s_K0s_TU_FM_biasd12->Sumw2();

V0chi2d1_diff_TU_FM_K0s->Sumw2();
V0chi2d2_diff_TU_FM_K0s->Sumw2();
V0chi2d12_diff_TU_FM_K0s->Sumw2();
    
//Lam
//True

hMass_Lam_Lam_T->Sumw2();   
hMass_Lam_1_T->Sumw2();   
hMass_Lam_2_T->Sumw2();   

  
hpeak_Lam_Lam_T->Sumw2();
hpeakcos_Lam_Lam_T->Sumw2();
hpeak_rot_Lam_Lam_T->Sumw2();
hpeak_inv_Lam_Lam_T->Sumw2();
hside_Lam_Lam_T->Sumw2();
hsideL_Lam_Lam_T->Sumw2();
hsideR_Lam_Lam_T->Sumw2();
hside_rot_Lam_Lam_T->Sumw2();
hside_inv_Lam_Lam_T->Sumw2();
hpeakside_Lam_Lam_T->Sumw2();
hpeaksideL_Lam_Lam_T->Sumw2();
hpeaksideR_Lam_Lam_T->Sumw2();
hpeakside_rot_Lam_Lam_T->Sumw2();
hpeakside_inv_Lam_Lam_T->Sumw2();
hpeakd1_Lam_Lam_T->Sumw2();
hpeakd2_Lam_Lam_T->Sumw2();
hpeakd12_Lam_Lam_T->Sumw2();

hpeak_Lam_Lam_T_mix->Sumw2();
hside_Lam_Lam_T_mix->Sumw2();
hsideL_Lam_Lam_T_mix->Sumw2();
hsideR_Lam_Lam_T_mix->Sumw2();
hpeakside_Lam_Lam_T_mix->Sumw2();
hpeaksideL_Lam_Lam_T_mix->Sumw2();
hpeaksideR_Lam_Lam_T_mix->Sumw2();

hpeak_Lam_Lam_T_etamix->Sumw2();
hside_Lam_Lam_T_etamix->Sumw2();
hsideL_Lam_Lam_T_etamix->Sumw2();
hsideR_Lam_Lam_T_etamix->Sumw2();
hpeakside_Lam_Lam_T_etamix->Sumw2();
hpeaksideL_Lam_Lam_T_etamix->Sumw2();
hpeaksideR_Lam_Lam_T_etamix->Sumw2();


V0chi2d1_diff_T_Lam->Sumw2();
V0chi2d2_diff_T_Lam->Sumw2();
V0chi2d12_diff_T_Lam->Sumw2();

//False

hMass_Lam_Lam_F->Sumw2();
  
hpeak_Lam_Lam_F->Sumw2();
hside_Lam_Lam_F->Sumw2();
hpeakside_Lam_Lam_F->Sumw2();
hpeakd1_Lam_Lam_F->Sumw2();
hpeakd2_Lam_Lam_F->Sumw2();
hpeakd12_Lam_Lam_F->Sumw2();
hpeak_Lam_Lam_F_bias->Sumw2();
hpeak_Lam_Lam_F_biasd1->Sumw2();
hpeak_Lam_Lam_F_biasd2->Sumw2();
hpeak_Lam_Lam_F_biasd12->Sumw2();

V0chi2d1_diff_F_Lam->Sumw2();
V0chi2d2_diff_F_Lam->Sumw2();
V0chi2d12_diff_F_Lam->Sumw2();
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
hpeak_Lam_Lam_TF->Sumw2();
hpeakd1_Lam_Lam_TF->Sumw2();
hpeakd2_Lam_Lam_TF->Sumw2();
hpeakd12_Lam_Lam_TF->Sumw2();
hside_Lam_Lam_TF->Sumw2();
hpeakside_Lam_Lam_TF->Sumw2();
hpeak_Lam_Lam_TF_bias->Sumw2();
hpeak_Lam_Lam_TF_biasd1->Sumw2();
hpeak_Lam_Lam_TF_biasd2->Sumw2();
hpeak_Lam_Lam_TF_biasd12->Sumw2();

V0chi2d1_diff_TF_Lam->Sumw2();
V0chi2d2_diff_TF_Lam->Sumw2();
V0chi2d12_diff_TF_Lam->Sumw2();

/////// matching hist ////// 
  
hMass_Lam_Lam_TM->Sumw2();
hMass_Lam_Lam_TU->Sumw2();
hMass_Lam_Lam_FM->Sumw2();
hMass_Lam_Lam_FU->Sumw2();
 
hpeak_Lam_Lam_TM->Sumw2();
hpeakd1_Lam_Lam_TM->Sumw2();
hpeakd2_Lam_Lam_TM->Sumw2();
hpeakd12_Lam_Lam_TM->Sumw2();
hpeak_Lam_Lam_TM_bias->Sumw2();
hpeak_Lam_Lam_TM_biasd1->Sumw2();
hpeak_Lam_Lam_TM_biasd2->Sumw2();
hpeak_Lam_Lam_TM_biasd12->Sumw2();

V0chi2d1_diff_TM_Lam->Sumw2();
V0chi2d2_diff_TM_Lam->Sumw2();
V0chi2d12_diff_TM_Lam->Sumw2();

hpeak_Lam_Lam_TU->Sumw2();
hpeakd1_Lam_Lam_TU->Sumw2();
hpeakd2_Lam_Lam_TU->Sumw2();
hpeakd12_Lam_Lam_TU->Sumw2();
hpeak_Lam_Lam_TU_bias->Sumw2();
hpeak_Lam_Lam_TU_biasd1->Sumw2();
hpeak_Lam_Lam_TU_biasd2->Sumw2();
hpeak_Lam_Lam_TU_biasd12->Sumw2();

V0chi2d1_diff_TU_Lam->Sumw2();
V0chi2d2_diff_TU_Lam->Sumw2();
V0chi2d12_diff_TU_Lam->Sumw2();

hpeak_Lam_Lam_FM->Sumw2();
hpeakd1_Lam_Lam_FM->Sumw2();
hpeakd2_Lam_Lam_FM->Sumw2();
hpeakd12_Lam_Lam_FM->Sumw2();
hpeak_Lam_Lam_FM_bias->Sumw2();
hpeak_Lam_Lam_FM_biasd1->Sumw2();
hpeak_Lam_Lam_FM_biasd2->Sumw2();
hpeak_Lam_Lam_FM_biasd12->Sumw2();

V0chi2d1_diff_FM_Lam->Sumw2();
V0chi2d2_diff_FM_Lam->Sumw2();
V0chi2d12_diff_FM_Lam->Sumw2();

hpeak_Lam_Lam_FU->Sumw2();
hpeakd1_Lam_Lam_FU->Sumw2();
hpeakd2_Lam_Lam_FU->Sumw2();
hpeakd12_Lam_Lam_FU->Sumw2();
hpeak_Lam_Lam_FU_bias->Sumw2();
hpeak_Lam_Lam_FU_biasd1->Sumw2();
hpeak_Lam_Lam_FU_biasd2->Sumw2();
hpeak_Lam_Lam_FU_biasd12->Sumw2();

V0chi2d1_diff_FU_Lam->Sumw2();
V0chi2d2_diff_FU_Lam->Sumw2();
V0chi2d12_diff_FU_Lam->Sumw2();


//true matched + true unmatched

hpeak_Lam_Lam_TM_TU->Sumw2();
hpeakd1_Lam_Lam_TM_TU->Sumw2();
hpeakd2_Lam_Lam_TM_TU->Sumw2();
hpeakd12_Lam_Lam_TM_TU->Sumw2();
hpeak_Lam_Lam_TM_TU_bias->Sumw2();
hpeak_Lam_Lam_TM_TU_biasd1->Sumw2();
hpeak_Lam_Lam_TM_TU_biasd2->Sumw2();
hpeak_Lam_Lam_TM_TU_biasd12->Sumw2();

V0chi2d1_diff_TM_TU_Lam->Sumw2();
V0chi2d2_diff_TM_TU_Lam->Sumw2();
V0chi2d12_diff_TM_TU_Lam->Sumw2();

//fake matched + fake unmatched

hpeak_Lam_Lam_FM_FU->Sumw2();
hpeakd1_Lam_Lam_FM_FU->Sumw2();
hpeakd2_Lam_Lam_FM_FU->Sumw2();
hpeakd12_Lam_Lam_FM_FU->Sumw2();
hpeak_Lam_Lam_FM_FU_bias->Sumw2();
hpeak_Lam_Lam_FM_FU_biasd1->Sumw2();
hpeak_Lam_Lam_FM_FU_biasd2->Sumw2();
hpeak_Lam_Lam_FM_FU_biasd12->Sumw2();

V0chi2d1_diff_FM_FU_Lam->Sumw2();
V0chi2d2_diff_FM_FU_Lam->Sumw2();
V0chi2d12_diff_FM_FU_Lam->Sumw2();



hpeak_Lam_Lam_TM_FM->Sumw2();
hpeakd1_Lam_Lam_TM_FM->Sumw2();
hpeakd2_Lam_Lam_TM_FM->Sumw2();
hpeakd12_Lam_Lam_TM_FM->Sumw2();
hpeak_Lam_Lam_TM_FM_bias->Sumw2();
hpeak_Lam_Lam_TM_FM_biasd1->Sumw2();
hpeak_Lam_Lam_TM_FM_biasd2->Sumw2();
hpeak_Lam_Lam_TM_FM_biasd12->Sumw2();

V0chi2d1_diff_TM_FM_Lam->Sumw2();
V0chi2d2_diff_TM_FM_Lam->Sumw2();
V0chi2d12_diff_TM_FM_Lam->Sumw2();

hpeak_Lam_Lam_TM_FU->Sumw2();
hpeakd1_Lam_Lam_TM_FU->Sumw2();
hpeakd2_Lam_Lam_TM_FU->Sumw2();
hpeakd12_Lam_Lam_TM_FU->Sumw2();
hpeak_Lam_Lam_TM_FU_bias->Sumw2();
hpeak_Lam_Lam_TM_FU_biasd1->Sumw2();
hpeak_Lam_Lam_TM_FU_biasd2->Sumw2();
hpeak_Lam_Lam_TM_FU_biasd12->Sumw2();

V0chi2d1_diff_TM_FU_Lam->Sumw2();
V0chi2d2_diff_TM_FU_Lam->Sumw2();
V0chi2d12_diff_TM_FU_Lam->Sumw2();

hpeak_Lam_Lam_TU_FU->Sumw2();
hpeakd1_Lam_Lam_TU_FU->Sumw2();
hpeakd2_Lam_Lam_TU_FU->Sumw2();
hpeakd12_Lam_Lam_TU_FU->Sumw2();
hpeak_Lam_Lam_TU_FU_bias->Sumw2();
hpeak_Lam_Lam_TU_FU_biasd1->Sumw2();
hpeak_Lam_Lam_TU_FU_biasd2->Sumw2();
hpeak_Lam_Lam_TU_FU_biasd12->Sumw2();

V0chi2d1_diff_TU_FU_Lam->Sumw2();
V0chi2d2_diff_TU_FU_Lam->Sumw2();
V0chi2d12_diff_TU_FU_Lam->Sumw2();

hpeak_Lam_Lam_TU_FM->Sumw2();
hpeakd1_Lam_Lam_TU_FM->Sumw2();
hpeakd2_Lam_Lam_TU_FM->Sumw2();
hpeakd12_Lam_Lam_TU_FM->Sumw2();
hpeak_Lam_Lam_TU_FM_bias->Sumw2();
hpeak_Lam_Lam_TU_FM_biasd1->Sumw2();
hpeak_Lam_Lam_TU_FM_biasd2->Sumw2();
hpeak_Lam_Lam_TU_FM_biasd12->Sumw2();

V0chi2d1_diff_TU_FM_Lam->Sumw2();
V0chi2d2_diff_TU_FM_Lam->Sumw2();
V0chi2d12_diff_TU_FM_Lam->Sumw2();
    
    
//ALam
//True

hMass_ALam_ALam_T->Sumw2();   
hMass_ALam_1_T->Sumw2();   
hMass_ALam_2_T->Sumw2();   
  
  
hpeak_ALam_ALam_T->Sumw2();
hpeakcos_ALam_ALam_T->Sumw2();
hpeak_rot_ALam_ALam_T->Sumw2();
hpeak_inv_ALam_ALam_T->Sumw2();
hside_ALam_ALam_T->Sumw2();
hsideL_ALam_ALam_T->Sumw2();
hsideR_ALam_ALam_T->Sumw2();
hside_rot_ALam_ALam_T->Sumw2();
hside_inv_ALam_ALam_T->Sumw2();
hpeakside_ALam_ALam_T->Sumw2();
hpeaksideL_ALam_ALam_T->Sumw2();
hpeaksideR_ALam_ALam_T->Sumw2();
hpeakside_rot_ALam_ALam_T->Sumw2();
hpeakside_inv_ALam_ALam_T->Sumw2();
hpeakd1_ALam_ALam_T->Sumw2();
hpeakd2_ALam_ALam_T->Sumw2();
hpeakd12_ALam_ALam_T->Sumw2();

hpeak_ALam_ALam_T_mix->Sumw2();
hside_ALam_ALam_T_mix->Sumw2();
hsideL_ALam_ALam_T_mix->Sumw2();
hsideR_ALam_ALam_T_mix->Sumw2();
hpeakside_ALam_ALam_T_mix->Sumw2();
hpeaksideL_ALam_ALam_T_mix->Sumw2();
hpeaksideR_ALam_ALam_T_mix->Sumw2();

hpeak_ALam_ALam_T_etamix->Sumw2();
hside_ALam_ALam_T_etamix->Sumw2();
hsideL_ALam_ALam_T_etamix->Sumw2();
hsideR_ALam_ALam_T_etamix->Sumw2();
hpeakside_ALam_ALam_T_etamix->Sumw2();
hpeaksideL_ALam_ALam_T_etamix->Sumw2();
hpeaksideR_ALam_ALam_T_etamix->Sumw2();


V0chi2d1_diff_T_ALam->Sumw2();
V0chi2d2_diff_T_ALam->Sumw2();
V0chi2d12_diff_T_ALam->Sumw2();
//False

hMass_ALam_ALam_F->Sumw2();
  
hpeak_ALam_ALam_F->Sumw2();
hside_ALam_ALam_F->Sumw2();
hpeakside_ALam_ALam_F->Sumw2();
hpeakd1_ALam_ALam_F->Sumw2();
hpeakd2_ALam_ALam_F->Sumw2();
hpeakd12_ALam_ALam_F->Sumw2();
hpeak_ALam_ALam_F_bias->Sumw2();
hpeak_ALam_ALam_F_biasd1->Sumw2();
hpeak_ALam_ALam_F_biasd2->Sumw2();
hpeak_ALam_ALam_F_biasd12->Sumw2();

V0chi2d1_diff_F_ALam->Sumw2();
V0chi2d2_diff_F_ALam->Sumw2();
V0chi2d12_diff_F_ALam->Sumw2();
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
hpeak_ALam_ALam_TF->Sumw2();
hpeakd1_ALam_ALam_TF->Sumw2();
hpeakd2_ALam_ALam_TF->Sumw2();
hpeakd12_ALam_ALam_TF->Sumw2();
hside_ALam_ALam_TF->Sumw2();
hpeakside_ALam_ALam_TF->Sumw2();
hpeak_ALam_ALam_TF_bias->Sumw2();
hpeak_ALam_ALam_TF_biasd1->Sumw2();
hpeak_ALam_ALam_TF_biasd2->Sumw2();
hpeak_ALam_ALam_TF_biasd12->Sumw2();

V0chi2d1_diff_TF_ALam->Sumw2();
V0chi2d2_diff_TF_ALam->Sumw2();
V0chi2d12_diff_TF_ALam->Sumw2();
  
/////// matching hist ////// 
  
hMass_ALam_ALam_TM->Sumw2();
hMass_ALam_ALam_TU->Sumw2();
hMass_ALam_ALam_FM->Sumw2();
hMass_ALam_ALam_FU->Sumw2();
 
hpeak_ALam_ALam_TM->Sumw2();
hpeakd1_ALam_ALam_TM->Sumw2();
hpeakd2_ALam_ALam_TM->Sumw2();
hpeakd12_ALam_ALam_TM->Sumw2();
hpeak_ALam_ALam_TM_bias->Sumw2();
hpeak_ALam_ALam_TM_biasd1->Sumw2();
hpeak_ALam_ALam_TM_biasd2->Sumw2();
hpeak_ALam_ALam_TM_biasd12->Sumw2();

V0chi2d1_diff_TM_ALam->Sumw2();
V0chi2d2_diff_TM_ALam->Sumw2();
V0chi2d12_diff_TM_ALam->Sumw2();

hpeak_ALam_ALam_TU->Sumw2();
hpeakd1_ALam_ALam_TU->Sumw2();
hpeakd2_ALam_ALam_TU->Sumw2();
hpeakd12_ALam_ALam_TU->Sumw2();
hpeak_ALam_ALam_TU_bias->Sumw2();
hpeak_ALam_ALam_TU_biasd1->Sumw2();
hpeak_ALam_ALam_TU_biasd2->Sumw2();
hpeak_ALam_ALam_TU_biasd12->Sumw2();

V0chi2d1_diff_TU_ALam->Sumw2();
V0chi2d2_diff_TU_ALam->Sumw2();
V0chi2d12_diff_TU_ALam->Sumw2();

hpeak_ALam_ALam_FM->Sumw2();
hpeakd1_ALam_ALam_FM->Sumw2();
hpeakd2_ALam_ALam_FM->Sumw2();
hpeakd12_ALam_ALam_FM->Sumw2();
hpeak_ALam_ALam_FM_bias->Sumw2();
hpeak_ALam_ALam_FM_biasd1->Sumw2();
hpeak_ALam_ALam_FM_biasd2->Sumw2();
hpeak_ALam_ALam_FM_biasd12->Sumw2();

V0chi2d1_diff_FM_ALam->Sumw2();
V0chi2d2_diff_FM_ALam->Sumw2();
V0chi2d12_diff_FM_ALam->Sumw2();

hpeak_ALam_ALam_FU->Sumw2();
hpeakd1_ALam_ALam_FU->Sumw2();
hpeakd2_ALam_ALam_FU->Sumw2();
hpeakd12_ALam_ALam_FU->Sumw2();
hpeak_ALam_ALam_FU_bias->Sumw2();
hpeak_ALam_ALam_FU_biasd1->Sumw2();
hpeak_ALam_ALam_FU_biasd2->Sumw2();
hpeak_ALam_ALam_FU_biasd12->Sumw2();

V0chi2d1_diff_FU_ALam->Sumw2();
V0chi2d2_diff_FU_ALam->Sumw2();
V0chi2d12_diff_FU_ALam->Sumw2();


//true matched + true unmatched

hpeak_ALam_ALam_TM_TU->Sumw2();
hpeakd1_ALam_ALam_TM_TU->Sumw2();
hpeakd2_ALam_ALam_TM_TU->Sumw2();
hpeakd12_ALam_ALam_TM_TU->Sumw2();
hpeak_ALam_ALam_TM_TU_bias->Sumw2();
hpeak_ALam_ALam_TM_TU_biasd1->Sumw2();
hpeak_ALam_ALam_TM_TU_biasd2->Sumw2();
hpeak_ALam_ALam_TM_TU_biasd12->Sumw2();

V0chi2d1_diff_TM_TU_ALam->Sumw2();
V0chi2d2_diff_TM_TU_ALam->Sumw2();
V0chi2d12_diff_TM_TU_ALam->Sumw2();

//fake matched + fake unmatched

hpeak_ALam_ALam_FM_FU->Sumw2();
hpeakd1_ALam_ALam_FM_FU->Sumw2();
hpeakd2_ALam_ALam_FM_FU->Sumw2();
hpeakd12_ALam_ALam_FM_FU->Sumw2();
hpeak_ALam_ALam_FM_FU_bias->Sumw2();
hpeak_ALam_ALam_FM_FU_biasd1->Sumw2();
hpeak_ALam_ALam_FM_FU_biasd2->Sumw2();
hpeak_ALam_ALam_FM_FU_biasd12->Sumw2();

V0chi2d1_diff_FM_FU_ALam->Sumw2();
V0chi2d2_diff_FM_FU_ALam->Sumw2();
V0chi2d12_diff_FM_FU_ALam->Sumw2();

hpeak_ALam_ALam_TM_FM->Sumw2();
hpeakd1_ALam_ALam_TM_FM->Sumw2();
hpeakd2_ALam_ALam_TM_FM->Sumw2();
hpeakd12_ALam_ALam_TM_FM->Sumw2();
hpeak_ALam_ALam_TM_FM_bias->Sumw2();
hpeak_ALam_ALam_TM_FM_biasd1->Sumw2();
hpeak_ALam_ALam_TM_FM_biasd2->Sumw2();
hpeak_ALam_ALam_TM_FM_biasd12->Sumw2();

V0chi2d1_diff_TM_FM_ALam->Sumw2();
V0chi2d2_diff_TM_FM_ALam->Sumw2();
V0chi2d12_diff_TM_FM_ALam->Sumw2();

hpeak_ALam_ALam_TM_FU->Sumw2();
hpeakd1_ALam_ALam_TM_FU->Sumw2();
hpeakd2_ALam_ALam_TM_FU->Sumw2();
hpeakd12_ALam_ALam_TM_FU->Sumw2();
hpeak_ALam_ALam_TM_FU_bias->Sumw2();
hpeak_ALam_ALam_TM_FU_biasd1->Sumw2();
hpeak_ALam_ALam_TM_FU_biasd2->Sumw2();
hpeak_ALam_ALam_TM_FU_biasd12->Sumw2();

V0chi2d1_diff_TM_FU_ALam->Sumw2();
V0chi2d2_diff_TM_FU_ALam->Sumw2();
V0chi2d12_diff_TM_FU_ALam->Sumw2();

hpeak_ALam_ALam_TU_FU->Sumw2();
hpeakd1_ALam_ALam_TU_FU->Sumw2();
hpeakd2_ALam_ALam_TU_FU->Sumw2();
hpeakd12_ALam_ALam_TU_FU->Sumw2();
hpeak_ALam_ALam_TU_FU_bias->Sumw2();
hpeak_ALam_ALam_TU_FU_biasd1->Sumw2();
hpeak_ALam_ALam_TU_FU_biasd2->Sumw2();
hpeak_ALam_ALam_TU_FU_biasd12->Sumw2();

V0chi2d1_diff_TU_FU_ALam->Sumw2();
V0chi2d2_diff_TU_FU_ALam->Sumw2();
V0chi2d12_diff_TU_FU_ALam->Sumw2();

hpeak_ALam_ALam_TU_FM->Sumw2();
hpeakd1_ALam_ALam_TU_FM->Sumw2();
hpeakd2_ALam_ALam_TU_FM->Sumw2();
hpeakd12_ALam_ALam_TU_FM->Sumw2();
hpeak_ALam_ALam_TU_FM_bias->Sumw2();
hpeak_ALam_ALam_TU_FM_biasd1->Sumw2();
hpeak_ALam_ALam_TU_FM_biasd2->Sumw2();
hpeak_ALam_ALam_TU_FM_biasd12->Sumw2();

V0chi2d1_diff_TU_FM_ALam->Sumw2();
V0chi2d2_diff_TU_FM_ALam->Sumw2();
V0chi2d12_diff_TU_FM_ALam->Sumw2();
    
//cross
 
//LAL
hMass_Lam_T->Sumw2(); 
hMass_ALam_T->Sumw2(); 
hMass_Lam_ALam_T->Sumw2(); 

hpeak_Lam_ALam_T->Sumw2(); 
hpeakcos_Lam_ALam_T->Sumw2();
hpeak_rot_Lam_ALam_T->Sumw2(); 
hpeak_inv_Lam_ALam_T->Sumw2(); 
hside_Lam_ALam_T->Sumw2(); 
hsideL_Lam_ALam_T->Sumw2(); 
hsideR_Lam_ALam_T->Sumw2(); 
hside_rot_Lam_ALam_T->Sumw2(); 
hside_inv_Lam_ALam_T->Sumw2(); 
hpeakside_Lam_ALam_T->Sumw2(); 
hpeaksideL_Lam_ALam_T->Sumw2(); 
hpeaksideR_Lam_ALam_T->Sumw2(); 
hpeakside_rot_Lam_ALam_T->Sumw2(); 
hpeakside_inv_Lam_ALam_T->Sumw2(); 
hpeakd1_Lam_ALam_T->Sumw2(); 
hpeakd2_Lam_ALam_T->Sumw2(); 
hpeakd12_Lam_ALam_T->Sumw2();
   

hpeak_Lam_ALam_T_mix->Sumw2(); 
hside_Lam_ALam_T_mix->Sumw2(); 
hsideL_Lam_ALam_T_mix->Sumw2(); 
hsideR_Lam_ALam_T_mix->Sumw2(); 
hpeakside_Lam_ALam_T_mix->Sumw2(); 
hpeaksideL_Lam_ALam_T_mix->Sumw2(); 
hpeaksideR_Lam_ALam_T_mix->Sumw2(); 

hpeak_Lam_ALam_T_etamix->Sumw2(); 
hside_Lam_ALam_T_etamix->Sumw2(); 
hsideL_Lam_ALam_T_etamix->Sumw2(); 
hsideR_Lam_ALam_T_etamix->Sumw2(); 
hpeakside_Lam_ALam_T_etamix->Sumw2(); 
hpeaksideL_Lam_ALam_T_etamix->Sumw2(); 
hpeaksideR_Lam_ALam_T_etamix->Sumw2(); 


V0chi2d1_diff_T_Lam_ALam->Sumw2(); 
V0chi2d2_diff_T_Lam_ALam->Sumw2(); 
V0chi2d12_diff_T_Lam_ALam->Sumw2();     


//------------------------------------------------------------------------
//fake
//------------------------------------------------------------------------ 

hMass_Lam_F->Sumw2();
hMass_ALam_F->Sumw2();
hMass_Lam_ALam_F->Sumw2();
  
hpeak_Lam_ALam_F->Sumw2();
hpeakd1_Lam_ALam_F->Sumw2();
hpeakd2_Lam_ALam_F->Sumw2();
hpeakd12_Lam_ALam_F->Sumw2();
hside_Lam_ALam_F->Sumw2();
hpeakside_Lam_ALam_F->Sumw2();
hpeak_Lam_ALam_F_bias->Sumw2();
hpeak_Lam_ALam_F_biasd1->Sumw2();
hpeak_Lam_ALam_F_biasd2->Sumw2();
hpeak_Lam_ALam_F_biasd12->Sumw2();

V0chi2d1_diff_F_Lam_ALam->Sumw2();
V0chi2d2_diff_F_Lam_ALam->Sumw2();
V0chi2d12_diff_F_Lam_ALam->Sumw2();


//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
hpeak_Lam_ALam_TF->Sumw2();
hpeakd1_Lam_ALam_TF->Sumw2();
hpeakd2_Lam_ALam_TF->Sumw2();
hpeakd12_Lam_ALam_TF->Sumw2();
hside_Lam_ALam_TF->Sumw2();
hpeakside_Lam_ALam_TF->Sumw2();
hpeak_Lam_ALam_TF_bias->Sumw2();
hpeak_Lam_ALam_TF_biasd1->Sumw2();
hpeak_Lam_ALam_TF_biasd2->Sumw2();
hpeak_Lam_ALam_TF_biasd12->Sumw2();

V0chi2d1_diff_TF_Lam_ALam->Sumw2();
V0chi2d2_diff_TF_Lam_ALam->Sumw2();
V0chi2d12_diff_TF_Lam_ALam->Sumw2();

  
/////// matching hist ////// 


hMass_Lam_ALam_TM->Sumw2();
hMass_Lam_TM->Sumw2();
hMass_ALam_TM->Sumw2();

hMass_Lam_ALam_TU->Sumw2();
hMass_Lam_TU->Sumw2();
hMass_ALam_TU->Sumw2();

hMass_Lam_ALam_FM->Sumw2();
hMass_Lam_FM->Sumw2();
hMass_ALam_FM->Sumw2();

hMass_Lam_ALam_FU->Sumw2();
hMass_Lam_FU->Sumw2();
hMass_ALam_FU->Sumw2();
 
hpeak_Lam_ALam_TM->Sumw2();
hpeakd1_Lam_ALam_TM->Sumw2();
hpeakd2_Lam_ALam_TM->Sumw2();
hpeakd12_Lam_ALam_TM->Sumw2();
hpeak_Lam_ALam_TM_bias->Sumw2();
hpeak_Lam_ALam_TM_biasd1->Sumw2();
hpeak_Lam_ALam_TM_biasd2->Sumw2();
hpeak_Lam_ALam_TM_biasd12->Sumw2();

V0chi2d1_diff_TM_Lam_ALam->Sumw2();
V0chi2d2_diff_TM_Lam_ALam->Sumw2();
V0chi2d12_diff_TM_Lam_ALam->Sumw2();

hpeak_Lam_ALam_TU->Sumw2();
hpeakd1_Lam_ALam_TU->Sumw2();
hpeakd2_Lam_ALam_TU->Sumw2();
hpeakd12_Lam_ALam_TU->Sumw2();
hpeak_Lam_ALam_TU_bias->Sumw2();
hpeak_Lam_ALam_TU_biasd1->Sumw2();
hpeak_Lam_ALam_TU_biasd2->Sumw2();
hpeak_Lam_ALam_TU_biasd12->Sumw2();

V0chi2d1_diff_TU_Lam_ALam->Sumw2();
V0chi2d2_diff_TU_Lam_ALam->Sumw2();
V0chi2d12_diff_TU_Lam_ALam->Sumw2();

hpeak_Lam_ALam_FM->Sumw2();
hpeakd1_Lam_ALam_FM->Sumw2();
hpeakd2_Lam_ALam_FM->Sumw2();
hpeakd12_Lam_ALam_FM->Sumw2();
hpeak_Lam_ALam_FM_bias->Sumw2();
hpeak_Lam_ALam_FM_biasd1->Sumw2();
hpeak_Lam_ALam_FM_biasd2->Sumw2();
hpeak_Lam_ALam_FM_biasd12->Sumw2();

V0chi2d1_diff_FM_Lam_ALam->Sumw2();
V0chi2d2_diff_FM_Lam_ALam->Sumw2();
V0chi2d12_diff_FM_Lam_ALam->Sumw2();

hpeak_Lam_ALam_FU->Sumw2();
hpeakd1_Lam_ALam_FU->Sumw2();
hpeakd2_Lam_ALam_FU->Sumw2();
hpeakd12_Lam_ALam_FU->Sumw2();
hpeak_Lam_ALam_FU_bias->Sumw2();
hpeak_Lam_ALam_FU_biasd1->Sumw2();
hpeak_Lam_ALam_FU_biasd2->Sumw2();
hpeak_Lam_ALam_FU_biasd12->Sumw2();

V0chi2d1_diff_FU_Lam_ALam->Sumw2();
V0chi2d2_diff_FU_Lam_ALam->Sumw2();
V0chi2d12_diff_FU_Lam_ALam->Sumw2();

//true matched + true unmatched

hpeak_Lam_ALam_TM_TU->Sumw2();
hpeakd1_Lam_ALam_TM_TU->Sumw2();
hpeakd2_Lam_ALam_TM_TU->Sumw2();
hpeakd12_Lam_ALam_TM_TU->Sumw2();
hpeak_Lam_ALam_TM_TU_bias->Sumw2();
hpeak_Lam_ALam_TM_TU_biasd1->Sumw2();
hpeak_Lam_ALam_TM_TU_biasd2->Sumw2();
hpeak_Lam_ALam_TM_TU_biasd12->Sumw2();

V0chi2d1_diff_TM_TU_Lam_ALam->Sumw2();
V0chi2d2_diff_TM_TU_Lam_ALam->Sumw2();
V0chi2d12_diff_TM_TU_Lam_ALam->Sumw2();

//fake matched + fake unmatched

hpeak_Lam_ALam_FM_FU->Sumw2();
hpeakd1_Lam_ALam_FM_FU->Sumw2();
hpeakd2_Lam_ALam_FM_FU->Sumw2();
hpeakd12_Lam_ALam_FM_FU->Sumw2();
hpeak_Lam_ALam_FM_FU_bias->Sumw2();
hpeak_Lam_ALam_FM_FU_biasd1->Sumw2();
hpeak_Lam_ALam_FM_FU_biasd2->Sumw2();
hpeak_Lam_ALam_FM_FU_biasd12->Sumw2();

V0chi2d1_diff_FM_FU_Lam_ALam->Sumw2();
V0chi2d2_diff_FM_FU_Lam_ALam->Sumw2();
V0chi2d12_diff_FM_FU_Lam_ALam->Sumw2();

//true matched + fake matched

hpeak_Lam_ALam_TM_FM->Sumw2();
hpeakd1_Lam_ALam_TM_FM->Sumw2();
hpeakd2_Lam_ALam_TM_FM->Sumw2();
hpeakd12_Lam_ALam_TM_FM->Sumw2();
hpeak_Lam_ALam_TM_FM_bias->Sumw2();
hpeak_Lam_ALam_TM_FM_biasd1->Sumw2();
hpeak_Lam_ALam_TM_FM_biasd2->Sumw2();
hpeak_Lam_ALam_TM_FM_biasd12->Sumw2();

V0chi2d1_diff_TM_FM_Lam_ALam->Sumw2();
V0chi2d2_diff_TM_FM_Lam_ALam->Sumw2();
V0chi2d12_diff_TM_FM_Lam_ALam->Sumw2();

//true matched + fake unmatched

hpeak_Lam_ALam_TM_FU->Sumw2();
hpeakd1_Lam_ALam_TM_FU->Sumw2();
hpeakd2_Lam_ALam_TM_FU->Sumw2();
hpeakd12_Lam_ALam_TM_FU->Sumw2();
hpeak_Lam_ALam_TM_FU_bias->Sumw2();
hpeak_Lam_ALam_TM_FU_biasd1->Sumw2();
hpeak_Lam_ALam_TM_FU_biasd2->Sumw2();
hpeak_Lam_ALam_TM_FU_biasd12->Sumw2();

V0chi2d1_diff_TM_FU_Lam_ALam->Sumw2();
V0chi2d2_diff_TM_FU_Lam_ALam->Sumw2();
V0chi2d12_diff_TM_FU_Lam_ALam->Sumw2();

//true unmatched + fake unmatched

hpeak_Lam_ALam_TU_FU->Sumw2();
hpeakd1_Lam_ALam_TU_FU->Sumw2();
hpeakd2_Lam_ALam_TU_FU->Sumw2();
hpeakd12_Lam_ALam_TU_FU->Sumw2();
hpeak_Lam_ALam_TU_FU_bias->Sumw2();
hpeak_Lam_ALam_TU_FU_biasd1->Sumw2();
hpeak_Lam_ALam_TU_FU_biasd2->Sumw2();
hpeak_Lam_ALam_TU_FU_biasd12->Sumw2();

V0chi2d1_diff_TU_FU_Lam_ALam->Sumw2();
V0chi2d2_diff_TU_FU_Lam_ALam->Sumw2();
V0chi2d12_diff_TU_FU_Lam_ALam->Sumw2();

//true unmatched + fake matched

hpeak_Lam_ALam_TU_FM->Sumw2();
hpeakd1_Lam_ALam_TU_FM->Sumw2();
hpeakd2_Lam_ALam_TU_FM->Sumw2();
hpeakd12_Lam_ALam_TU_FM->Sumw2();
hpeak_Lam_ALam_TU_FM_bias->Sumw2();
hpeak_Lam_ALam_TU_FM_biasd1->Sumw2();
hpeak_Lam_ALam_TU_FM_biasd2->Sumw2();
hpeak_Lam_ALam_TU_FM_biasd12->Sumw2();

V0chi2d1_diff_TU_FM_Lam_ALam->Sumw2();
V0chi2d2_diff_TU_FM_Lam_ALam->Sumw2();
V0chi2d12_diff_TU_FM_Lam_ALam->Sumw2();
    
  

hMass_K0sL_T->Sumw2(); 
hMass_LamK_T->Sumw2(); 
hMass_K0s_Lam_T->Sumw2();
hpeak_K0s_Lam_T->Sumw2();
hpeakcos_K0s_Lam_T->Sumw2();
hpeak_rot_K0s_Lam_T->Sumw2();
hpeak_inv_K0s_Lam_T->Sumw2();
hside_K0s_Lam_T->Sumw2();
hsideL_K0s_Lam_T->Sumw2();
hsideR_K0s_Lam_T->Sumw2();
hside_rot_K0s_Lam_T->Sumw2();
hside_inv_K0s_Lam_T->Sumw2();
hpeakside_K0s_Lam_T->Sumw2();
hpeaksideL_K0s_Lam_T->Sumw2();
hpeaksideR_K0s_Lam_T->Sumw2();
hpeakside_rot_K0s_Lam_T->Sumw2();
hpeakside_inv_K0s_Lam_T->Sumw2();
hpeakd1_K0s_Lam_T->Sumw2();
hpeakd2_K0s_Lam_T->Sumw2();
hpeakd12_K0s_Lam_T->Sumw2();
hpeak_K0s_Lam_T_mix->Sumw2();
hside_K0s_Lam_T_mix->Sumw2();
hsideL_K0s_Lam_T_mix->Sumw2();
hsideR_K0s_Lam_T_mix->Sumw2();
hpeakside_K0s_Lam_T_mix->Sumw2();
hpeaksideL_K0s_Lam_T_mix->Sumw2();
hpeaksideR_K0s_Lam_T_mix->Sumw2();
hpeak_K0s_Lam_T_etamix->Sumw2();
hside_K0s_Lam_T_etamix->Sumw2();
hsideL_K0s_Lam_T_etamix->Sumw2();
hsideR_K0s_Lam_T_etamix->Sumw2();
hpeakside_K0s_Lam_T_etamix->Sumw2();
hpeaksideL_K0s_Lam_T_etamix->Sumw2();
hpeaksideR_K0s_Lam_T_etamix->Sumw2();
V0chi2d1_diff_T_K0s_Lam->Sumw2();
V0chi2d2_diff_T_K0s_Lam->Sumw2();
V0chi2d12_diff_T_K0s_Lam->Sumw2();
  
//------------------------------------------------------------------------
//fake
//------------------------------------------------------------------------ 

hMass_K0sL_F->Sumw2();   
hMass_LamK_F->Sumw2();
hMass_K0s_Lam_F->Sumw2();
  
hpeak_K0s_Lam_F->Sumw2();
hpeakd1_K0s_Lam_F->Sumw2();
hpeakd2_K0s_Lam_F->Sumw2();
hpeakd12_K0s_Lam_F->Sumw2();
hside_K0s_Lam_F->Sumw2();
hpeakside_K0s_Lam_F->Sumw2();
hpeak_K0s_Lam_F_bias->Sumw2();
hpeak_K0s_Lam_F_biasd1->Sumw2();
hpeak_K0s_Lam_F_biasd2->Sumw2();
hpeak_K0s_Lam_F_biasd12->Sumw2();

V0chi2d1_diff_F_K0s_Lam->Sumw2();
V0chi2d2_diff_F_K0s_Lam->Sumw2();
V0chi2d12_diff_F_K0s_Lam->Sumw2();
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
hpeak_K0s_Lam_TF->Sumw2();
hpeakd1_K0s_Lam_TF->Sumw2();
hpeakd2_K0s_Lam_TF->Sumw2();
hpeakd12_K0s_Lam_TF->Sumw2();
hside_K0s_Lam_TF->Sumw2();
hpeakside_K0s_Lam_TF->Sumw2(); 
hpeak_K0s_Lam_TF_bias->Sumw2();
hpeak_K0s_Lam_TF_biasd1->Sumw2();
hpeak_K0s_Lam_TF_biasd2->Sumw2();
hpeak_K0s_Lam_TF_biasd12->Sumw2();

V0chi2d1_diff_TF_K0s_Lam->Sumw2();
V0chi2d2_diff_TF_K0s_Lam->Sumw2();
V0chi2d12_diff_TF_K0s_Lam->Sumw2();

  
/////// matching hist ////// 

hMass_K0sL_TM->Sumw2();  
hMass_LamK_TM->Sumw2();
hMass_K0s_Lam_TM->Sumw2();

hMass_K0sL_TU->Sumw2();
hMass_LamK_TU->Sumw2();
hMass_K0s_Lam_TU->Sumw2();

hMass_K0sL_FM->Sumw2();
hMass_LamK_FM->Sumw2();
hMass_K0s_Lam_FM->Sumw2();

hMass_K0sL_FU->Sumw2();
hMass_LamK_FU->Sumw2();
hMass_K0s_Lam_FU->Sumw2();


hpeak_K0s_Lam_TM->Sumw2();
hpeakd1_K0s_Lam_TM->Sumw2();
hpeakd2_K0s_Lam_TM->Sumw2();
hpeakd12_K0s_Lam_TM->Sumw2();
hpeak_K0s_Lam_TM_bias->Sumw2();
hpeak_K0s_Lam_TM_biasd1->Sumw2();
hpeak_K0s_Lam_TM_biasd2->Sumw2();
hpeak_K0s_Lam_TM_biasd12->Sumw2();

V0chi2d1_diff_TM_K0s_Lam->Sumw2();
V0chi2d2_diff_TM_K0s_Lam->Sumw2();
V0chi2d12_diff_TM_K0s_Lam->Sumw2();

hpeak_K0s_Lam_TU->Sumw2();
hpeakd1_K0s_Lam_TU->Sumw2();
hpeakd2_K0s_Lam_TU->Sumw2();
hpeakd12_K0s_Lam_TU->Sumw2();
hpeak_K0s_Lam_TU_bias->Sumw2();
hpeak_K0s_Lam_TU_biasd1->Sumw2();
hpeak_K0s_Lam_TU_biasd2->Sumw2();
hpeak_K0s_Lam_TU_biasd12->Sumw2();

V0chi2d1_diff_TU_K0s_Lam->Sumw2();
V0chi2d2_diff_TU_K0s_Lam->Sumw2();
V0chi2d12_diff_TU_K0s_Lam->Sumw2();

hpeak_K0s_Lam_FM->Sumw2();
hpeakd1_K0s_Lam_FM->Sumw2();
hpeakd2_K0s_Lam_FM->Sumw2();
hpeakd12_K0s_Lam_FM->Sumw2();
hpeak_K0s_Lam_FM_bias->Sumw2();
hpeak_K0s_Lam_FM_biasd1->Sumw2();
hpeak_K0s_Lam_FM_biasd2->Sumw2();
hpeak_K0s_Lam_FM_biasd12->Sumw2();

V0chi2d1_diff_FM_K0s_Lam->Sumw2();
V0chi2d2_diff_FM_K0s_Lam->Sumw2();
V0chi2d12_diff_FM_K0s_Lam->Sumw2();

hpeak_K0s_Lam_FU->Sumw2();
hpeakd1_K0s_Lam_FU->Sumw2();
hpeakd2_K0s_Lam_FU->Sumw2();
hpeakd12_K0s_Lam_FU->Sumw2();
hpeak_K0s_Lam_FU_bias->Sumw2();
hpeak_K0s_Lam_FU_biasd1->Sumw2();
hpeak_K0s_Lam_FU_biasd2->Sumw2();
hpeak_K0s_Lam_FU_biasd12->Sumw2();

V0chi2d1_diff_FU_K0s_Lam->Sumw2();
V0chi2d2_diff_FU_K0s_Lam->Sumw2();
V0chi2d12_diff_FU_K0s_Lam->Sumw2();

//true matched + true unmatched

hpeak_K0s_Lam_TM_TU->Sumw2();
hpeakd1_K0s_Lam_TM_TU->Sumw2();
hpeakd2_K0s_Lam_TM_TU->Sumw2();
hpeakd12_K0s_Lam_TM_TU->Sumw2();
hpeak_K0s_Lam_TM_TU_bias->Sumw2();
hpeak_K0s_Lam_TM_TU_biasd1->Sumw2();
hpeak_K0s_Lam_TM_TU_biasd2->Sumw2();
hpeak_K0s_Lam_TM_TU_biasd12->Sumw2();

V0chi2d1_diff_TM_TU_K0s_Lam->Sumw2();
V0chi2d2_diff_TM_TU_K0s_Lam->Sumw2();
V0chi2d12_diff_TM_TU_K0s_Lam->Sumw2();

//fake matched + fake unmatched

hpeak_K0s_Lam_FM_FU->Sumw2();
hpeakd1_K0s_Lam_FM_FU->Sumw2();
hpeakd2_K0s_Lam_FM_FU->Sumw2();
hpeakd12_K0s_Lam_FM_FU->Sumw2();
hpeak_K0s_Lam_FM_FU_bias->Sumw2();
hpeak_K0s_Lam_FM_FU_biasd1->Sumw2();
hpeak_K0s_Lam_FM_FU_biasd2->Sumw2();
hpeak_K0s_Lam_FM_FU_biasd12->Sumw2();

V0chi2d1_diff_FM_FU_K0s_Lam->Sumw2();
V0chi2d2_diff_FM_FU_K0s_Lam->Sumw2();
V0chi2d12_diff_FM_FU_K0s_Lam->Sumw2();

//true matched + fake matched

hpeak_K0s_Lam_TM_FM->Sumw2();
hpeakd1_K0s_Lam_TM_FM->Sumw2();
hpeakd2_K0s_Lam_TM_FM->Sumw2();
hpeakd12_K0s_Lam_TM_FM->Sumw2();
hpeak_K0s_Lam_TM_FM_bias->Sumw2();
hpeak_K0s_Lam_TM_FM_biasd1->Sumw2();
hpeak_K0s_Lam_TM_FM_biasd2->Sumw2();
hpeak_K0s_Lam_TM_FM_biasd12->Sumw2();

V0chi2d1_diff_TM_FM_K0s_Lam->Sumw2();
V0chi2d2_diff_TM_FM_K0s_Lam->Sumw2();
V0chi2d12_diff_TM_FM_K0s_Lam->Sumw2();

//true matched + fake unmatched

hpeak_K0s_Lam_TM_FU->Sumw2();
hpeakd1_K0s_Lam_TM_FU->Sumw2();
hpeakd2_K0s_Lam_TM_FU->Sumw2();
hpeakd12_K0s_Lam_TM_FU->Sumw2();
hpeak_K0s_Lam_TM_FU_bias->Sumw2();
hpeak_K0s_Lam_TM_FU_biasd1->Sumw2();
hpeak_K0s_Lam_TM_FU_biasd2->Sumw2();
hpeak_K0s_Lam_TM_FU_biasd12->Sumw2();

V0chi2d1_diff_TM_FU_K0s_Lam->Sumw2();
V0chi2d2_diff_TM_FU_K0s_Lam->Sumw2();
V0chi2d12_diff_TM_FU_K0s_Lam->Sumw2();

//true unmatched + fake unmatched

hpeak_K0s_Lam_TU_FU->Sumw2();
hpeakd1_K0s_Lam_TU_FU->Sumw2();
hpeakd2_K0s_Lam_TU_FU->Sumw2();
hpeakd12_K0s_Lam_TU_FU->Sumw2();
hpeak_K0s_Lam_TU_FU_bias->Sumw2();
hpeak_K0s_Lam_TU_FU_biasd1->Sumw2();
hpeak_K0s_Lam_TU_FU_biasd2->Sumw2();
hpeak_K0s_Lam_TU_FU_biasd12->Sumw2();

V0chi2d1_diff_TU_FU_K0s_Lam->Sumw2();
V0chi2d2_diff_TU_FU_K0s_Lam->Sumw2();
V0chi2d12_diff_TU_FU_K0s_Lam->Sumw2();

//true unmatched + fake matched

hpeak_K0s_Lam_TU_FM->Sumw2();
hpeakd1_K0s_Lam_TU_FM->Sumw2();
hpeakd2_K0s_Lam_TU_FM->Sumw2();
hpeakd12_K0s_Lam_TU_FM->Sumw2();
hpeak_K0s_Lam_TU_FM_bias->Sumw2();
hpeak_K0s_Lam_TU_FM_biasd1->Sumw2();
hpeak_K0s_Lam_TU_FM_biasd2->Sumw2();
hpeak_K0s_Lam_TU_FM_biasd12->Sumw2();

V0chi2d1_diff_TU_FM_K0s_Lam->Sumw2();
V0chi2d2_diff_TU_FM_K0s_Lam->Sumw2();
V0chi2d12_diff_TU_FM_K0s_Lam->Sumw2();

hMass_K0sAL_T->Sumw2(); 
hMass_ALamK_T->Sumw2(); 
hMass_K0s_ALam_T->Sumw2();
hpeak_K0s_ALam_T->Sumw2();
hpeakcos_K0s_ALam_T->Sumw2();
hpeak_rot_K0s_ALam_T->Sumw2();
hpeak_inv_K0s_ALam_T->Sumw2();
hside_K0s_ALam_T->Sumw2();
hsideL_K0s_ALam_T->Sumw2();
hsideR_K0s_ALam_T->Sumw2();
hside_rot_K0s_ALam_T->Sumw2();
hside_inv_K0s_ALam_T->Sumw2();
hpeakside_K0s_ALam_T->Sumw2();
hpeaksideL_K0s_ALam_T->Sumw2();
hpeaksideR_K0s_ALam_T->Sumw2();
hpeakside_rot_K0s_ALam_T->Sumw2();
hpeakside_inv_K0s_ALam_T->Sumw2();
hpeakd1_K0s_ALam_T->Sumw2();
hpeakd2_K0s_ALam_T->Sumw2();
hpeakd12_K0s_ALam_T->Sumw2();
hpeak_K0s_ALam_T_mix->Sumw2();
hside_K0s_ALam_T_mix->Sumw2();
hsideL_K0s_ALam_T_mix->Sumw2();
hsideR_K0s_ALam_T_mix->Sumw2();
hpeakside_K0s_ALam_T_mix->Sumw2();
hpeaksideL_K0s_ALam_T_mix->Sumw2();
hpeaksideR_K0s_ALam_T_mix->Sumw2();
hpeak_K0s_ALam_T_etamix->Sumw2();
hside_K0s_ALam_T_etamix->Sumw2();
hsideL_K0s_ALam_T_etamix->Sumw2();
hsideR_K0s_ALam_T_etamix->Sumw2();
hpeakside_K0s_ALam_T_etamix->Sumw2();
hpeaksideL_K0s_ALam_T_etamix->Sumw2();
hpeaksideR_K0s_ALam_T_etamix->Sumw2();
V0chi2d1_diff_T_K0s_ALam->Sumw2();
V0chi2d2_diff_T_K0s_ALam->Sumw2();
V0chi2d12_diff_T_K0s_ALam->Sumw2();
  
//------------------------------------------------------------------------
//fake
//------------------------------------------------------------------------ 

hMass_K0sAL_F->Sumw2();   
hMass_ALamK_F->Sumw2();
hMass_K0s_ALam_F->Sumw2();
hpeak_K0s_ALam_F->Sumw2();
hpeakd1_K0s_ALam_F->Sumw2();
hpeakd2_K0s_ALam_F->Sumw2();
hpeakd12_K0s_ALam_F->Sumw2();
hside_K0s_ALam_F->Sumw2();
hpeakside_K0s_ALam_F->Sumw2();
hpeak_K0s_ALam_F_bias->Sumw2();
hpeak_K0s_ALam_F_biasd1->Sumw2();
hpeak_K0s_ALam_F_biasd2->Sumw2();
hpeak_K0s_ALam_F_biasd12->Sumw2();
V0chi2d1_diff_F_K0s_ALam->Sumw2();
V0chi2d2_diff_F_K0s_ALam->Sumw2();
V0chi2d12_diff_F_K0s_ALam->Sumw2();
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
hpeak_K0s_ALam_TF->Sumw2();
hpeakd1_K0s_ALam_TF->Sumw2();
hpeakd2_K0s_ALam_TF->Sumw2();
hpeakd12_K0s_ALam_TF->Sumw2();
hside_K0s_ALam_TF->Sumw2();
hpeakside_K0s_ALam_TF->Sumw2(); 
hpeak_K0s_ALam_TF_bias->Sumw2();
hpeak_K0s_ALam_TF_biasd1->Sumw2();
hpeak_K0s_ALam_TF_biasd2->Sumw2();
hpeak_K0s_ALam_TF_biasd12->Sumw2();

V0chi2d1_diff_TF_K0s_ALam->Sumw2();
V0chi2d2_diff_TF_K0s_ALam->Sumw2();
V0chi2d12_diff_TF_K0s_ALam->Sumw2();

  
/////// matching hist ////// 

hMass_K0sAL_TM->Sumw2();  
hMass_ALamK_TM->Sumw2();
hMass_K0s_ALam_TM->Sumw2();

hMass_K0sAL_TU->Sumw2();
hMass_ALamK_TU->Sumw2();
hMass_K0s_ALam_TU->Sumw2();

hMass_K0sAL_FM->Sumw2();
hMass_ALamK_FM->Sumw2();
hMass_K0s_ALam_FM->Sumw2();

hMass_K0sAL_FU->Sumw2();
hMass_ALamK_FU->Sumw2();
hMass_K0s_ALam_FU->Sumw2();


hpeak_K0s_ALam_TM->Sumw2();
hpeakd1_K0s_ALam_TM->Sumw2();
hpeakd2_K0s_ALam_TM->Sumw2();
hpeakd12_K0s_ALam_TM->Sumw2();
hpeak_K0s_ALam_TM_bias->Sumw2();
hpeak_K0s_ALam_TM_biasd1->Sumw2();
hpeak_K0s_ALam_TM_biasd2->Sumw2();
hpeak_K0s_ALam_TM_biasd12->Sumw2();

V0chi2d1_diff_TM_K0s_ALam->Sumw2();
V0chi2d2_diff_TM_K0s_ALam->Sumw2();
V0chi2d12_diff_TM_K0s_ALam->Sumw2();

hpeak_K0s_ALam_TU->Sumw2();
hpeakd1_K0s_ALam_TU->Sumw2();
hpeakd2_K0s_ALam_TU->Sumw2();
hpeakd12_K0s_ALam_TU->Sumw2();
hpeak_K0s_ALam_TU_bias->Sumw2();
hpeak_K0s_ALam_TU_biasd1->Sumw2();
hpeak_K0s_ALam_TU_biasd2->Sumw2();
hpeak_K0s_ALam_TU_biasd12->Sumw2();

V0chi2d1_diff_TU_K0s_ALam->Sumw2();
V0chi2d2_diff_TU_K0s_ALam->Sumw2();
V0chi2d12_diff_TU_K0s_ALam->Sumw2();

hpeak_K0s_ALam_FM->Sumw2();
hpeakd1_K0s_ALam_FM->Sumw2();
hpeakd2_K0s_ALam_FM->Sumw2();
hpeakd12_K0s_ALam_FM->Sumw2();
hpeak_K0s_ALam_FM_bias->Sumw2();
hpeak_K0s_ALam_FM_biasd1->Sumw2();
hpeak_K0s_ALam_FM_biasd2->Sumw2();
hpeak_K0s_ALam_FM_biasd12->Sumw2();

V0chi2d1_diff_FM_K0s_ALam->Sumw2();
V0chi2d2_diff_FM_K0s_ALam->Sumw2();
V0chi2d12_diff_FM_K0s_ALam->Sumw2();

hpeak_K0s_ALam_FU->Sumw2();
hpeakd1_K0s_ALam_FU->Sumw2();
hpeakd2_K0s_ALam_FU->Sumw2();
hpeakd12_K0s_ALam_FU->Sumw2();
hpeak_K0s_ALam_FU_bias->Sumw2();
hpeak_K0s_ALam_FU_biasd1->Sumw2();
hpeak_K0s_ALam_FU_biasd2->Sumw2();
hpeak_K0s_ALam_FU_biasd12->Sumw2();

V0chi2d1_diff_FU_K0s_ALam->Sumw2();
V0chi2d2_diff_FU_K0s_ALam->Sumw2();
V0chi2d12_diff_FU_K0s_ALam->Sumw2();

//true matched + true unmatched

hpeak_K0s_ALam_TM_TU->Sumw2();
hpeakd1_K0s_ALam_TM_TU->Sumw2();
hpeakd2_K0s_ALam_TM_TU->Sumw2();
hpeakd12_K0s_ALam_TM_TU->Sumw2();
hpeak_K0s_ALam_TM_TU_bias->Sumw2();
hpeak_K0s_ALam_TM_TU_biasd1->Sumw2();
hpeak_K0s_ALam_TM_TU_biasd2->Sumw2();
hpeak_K0s_ALam_TM_TU_biasd12->Sumw2();

V0chi2d1_diff_TM_TU_K0s_ALam->Sumw2();
V0chi2d2_diff_TM_TU_K0s_ALam->Sumw2();
V0chi2d12_diff_TM_TU_K0s_ALam->Sumw2();

//fake matched + fake unmatched

hpeak_K0s_ALam_FM_FU->Sumw2();
hpeakd1_K0s_ALam_FM_FU->Sumw2();
hpeakd2_K0s_ALam_FM_FU->Sumw2();
hpeakd12_K0s_ALam_FM_FU->Sumw2();
hpeak_K0s_ALam_FM_FU_bias->Sumw2();
hpeak_K0s_ALam_FM_FU_biasd1->Sumw2();
hpeak_K0s_ALam_FM_FU_biasd2->Sumw2();
hpeak_K0s_ALam_FM_FU_biasd12->Sumw2();

V0chi2d1_diff_FM_FU_K0s_ALam->Sumw2();
V0chi2d2_diff_FM_FU_K0s_ALam->Sumw2();
V0chi2d12_diff_FM_FU_K0s_ALam->Sumw2();

//true matched + fake matched

hpeak_K0s_ALam_TM_FM->Sumw2();
hpeakd1_K0s_ALam_TM_FM->Sumw2();
hpeakd2_K0s_ALam_TM_FM->Sumw2();
hpeakd12_K0s_ALam_TM_FM->Sumw2();
hpeak_K0s_ALam_TM_FM_bias->Sumw2();
hpeak_K0s_ALam_TM_FM_biasd1->Sumw2();
hpeak_K0s_ALam_TM_FM_biasd2->Sumw2();
hpeak_K0s_ALam_TM_FM_biasd12->Sumw2();

V0chi2d1_diff_TM_FM_K0s_ALam->Sumw2();
V0chi2d2_diff_TM_FM_K0s_ALam->Sumw2();
V0chi2d12_diff_TM_FM_K0s_ALam->Sumw2();

//true matched + fake unmatched

hpeak_K0s_ALam_TM_FU->Sumw2();
hpeakd1_K0s_ALam_TM_FU->Sumw2();
hpeakd2_K0s_ALam_TM_FU->Sumw2();
hpeakd12_K0s_ALam_TM_FU->Sumw2();
hpeak_K0s_ALam_TM_FU_bias->Sumw2();
hpeak_K0s_ALam_TM_FU_biasd1->Sumw2();
hpeak_K0s_ALam_TM_FU_biasd2->Sumw2();
hpeak_K0s_ALam_TM_FU_biasd12->Sumw2();

V0chi2d1_diff_TM_FU_K0s_ALam->Sumw2();
V0chi2d2_diff_TM_FU_K0s_ALam->Sumw2();
V0chi2d12_diff_TM_FU_K0s_ALam->Sumw2();

//true unmatched + fake unmatched

hpeak_K0s_ALam_TU_FU->Sumw2();
hpeakd1_K0s_ALam_TU_FU->Sumw2();
hpeakd2_K0s_ALam_TU_FU->Sumw2();
hpeakd12_K0s_ALam_TU_FU->Sumw2();
hpeak_K0s_ALam_TU_FU_bias->Sumw2();
hpeak_K0s_ALam_TU_FU_biasd1->Sumw2();
hpeak_K0s_ALam_TU_FU_biasd2->Sumw2();
hpeak_K0s_ALam_TU_FU_biasd12->Sumw2();

V0chi2d1_diff_TU_FU_K0s_ALam->Sumw2();
V0chi2d2_diff_TU_FU_K0s_ALam->Sumw2();
V0chi2d12_diff_TU_FU_K0s_ALam->Sumw2();

//true unmatched + fake matched

hpeak_K0s_ALam_TU_FM->Sumw2();
hpeakd1_K0s_ALam_TU_FM->Sumw2();
hpeakd2_K0s_ALam_TU_FM->Sumw2();
hpeakd12_K0s_ALam_TU_FM->Sumw2();
hpeak_K0s_ALam_TU_FM_bias->Sumw2();
hpeak_K0s_ALam_TU_FM_biasd1->Sumw2();
hpeak_K0s_ALam_TU_FM_biasd2->Sumw2();
hpeak_K0s_ALam_TU_FM_biasd12->Sumw2();

V0chi2d1_diff_TU_FM_K0s_ALam->Sumw2();
V0chi2d2_diff_TU_FM_K0s_ALam->Sumw2();
V0chi2d12_diff_TU_FM_K0s_ALam->Sumw2();

hpeak_K0s_K0s_TM_mix->Sumw2(); 
hpeak_K0s_K0s_TU_mix->Sumw2(); 
hpeak_K0s_K0s_TM_TU_mix->Sumw2(); 

hpeak_Lam_Lam_TM_mix->Sumw2(); 
hpeak_Lam_Lam_TU_mix->Sumw2(); 
hpeak_Lam_Lam_TM_TU_mix->Sumw2(); 

hpeak_ALam_ALam_TM_mix->Sumw2(); 
hpeak_ALam_ALam_TU_mix->Sumw2(); 
hpeak_ALam_ALam_TM_TU_mix->Sumw2(); 

hpeak_Lam_ALam_TM_mix->Sumw2(); 
hpeak_Lam_ALam_TU_mix->Sumw2(); 
hpeak_Lam_ALam_TM_TU_mix->Sumw2(); 

hpeak_K0s_Lam_TM_mix->Sumw2(); 
hpeak_K0s_Lam_TU_mix->Sumw2(); 
hpeak_K0s_Lam_TM_TU_mix->Sumw2(); 

hpeak_K0s_ALam_TM_mix->Sumw2(); 
hpeak_K0s_ALam_TU_mix->Sumw2(); 
hpeak_K0s_ALam_TM_TU_mix->Sumw2(); 


V0_chi2_SS_LLFD_XI->Sumw2(); 
V0_chi2_OS_LLFD_XI->Sumw2(); 
V0_chi2_SS_ALFD_XI->Sumw2(); 
V0_chi2_OS_ALFD_XI->Sumw2(); 
V0_chi2_SS_KLFD_XI->Sumw2(); 
V0_chi2_OS_KLFD_XI->Sumw2(); 

V0_chi2_SS_LLFD_AXI->Sumw2(); 
V0_chi2_OS_LLFD_AXI->Sumw2(); 
V0_chi2_SS_ALFD_AXI->Sumw2(); 
V0_chi2_OS_ALFD_AXI->Sumw2(); 
V0_chi2_SS_KLFD_AXI->Sumw2(); 
V0_chi2_OS_KLFD_AXI->Sumw2(); 

V0_chi2_SS_LLFD_OM->Sumw2(); 
V0_chi2_OS_LLFD_OM->Sumw2(); 
V0_chi2_SS_ALFD_OM->Sumw2(); 
V0_chi2_OS_ALFD_OM->Sumw2(); 
V0_chi2_SS_KLFD_OM->Sumw2(); 
V0_chi2_OS_KLFD_OM->Sumw2(); 

V0_chi2_SS_LLFD_AOM->Sumw2(); 
V0_chi2_OS_LLFD_AOM->Sumw2(); 
V0_chi2_SS_ALFD_AOM->Sumw2(); 
V0_chi2_OS_ALFD_AOM->Sumw2(); 
V0_chi2_SS_KLFD_AOM->Sumw2(); 
V0_chi2_OS_KLFD_AOM->Sumw2(); 



hpeak_Lam_LamFD->Sumw2(); 
hpeak_LamFD_LamFD->Sumw2(); 
hpeak_Lam_LamFD_mix->Sumw2(); 
hpeak_LamFD_LamFD_mix->Sumw2(); 

hpeak_ALam_ALamFD->Sumw2(); 
hpeak_ALamFD_ALamFD->Sumw2(); 
hpeak_ALam_ALamFD_mix->Sumw2(); 
hpeak_ALamFD_ALamFD_mix->Sumw2(); 

hpeak_LALFD->Sumw2(); 
hpeak_LALFD_mix->Sumw2(); 
hpeak_LFDALFD->Sumw2(); 
hpeak_LFDALFD_mix->Sumw2(); 


hpeak_KLFD->Sumw2(); 
hpeak_KLFD_mix->Sumw2(); 

hpeak_KALFD->Sumw2(); 
hpeak_KALFD_mix->Sumw2(); 


h_HDibaryon_Lam_Lam->Sumw2(); 
h_HDibaryon_rot_Lam_Lam->Sumw2(); 
h_HDibaryon_inv_Lam_Lam->Sumw2(); 
h_HDibaryon_Lam_Lam_mix->Sumw2(); 

h_HDibaryon_Lam_Lam_side->Sumw2(); 
h_HDibaryon_rot_Lam_Lam_side->Sumw2(); 
h_HDibaryon_inv_Lam_Lam_side->Sumw2(); 
h_HDibaryon_Lam_Lam_mix_side->Sumw2(); 

h_HDibaryon_Lam_Lam_peakside->Sumw2(); 
h_HDibaryon_rot_Lam_Lam_peakside->Sumw2(); 
h_HDibaryon_inv_Lam_Lam_peakside->Sumw2(); 
h_HDibaryon_Lam_Lam_mix_peakside->Sumw2(); 


h_HDibaryon_ALam_ALam->Sumw2(); 
h_HDibaryon_rot_ALam_ALam->Sumw2(); 
h_HDibaryon_inv_ALam_ALam->Sumw2(); 
h_HDibaryon_ALam_ALam_mix->Sumw2(); 

h_HDibaryon_ALam_ALam_peakside->Sumw2(); 
h_HDibaryon_rot_ALam_ALam_peakside->Sumw2(); 
h_HDibaryon_inv_ALam_ALam_peakside->Sumw2(); 
h_HDibaryon_ALam_ALam_mix_peakside->Sumw2(); 

h_HDibaryon_Lam_ALam->Sumw2(); 
h_HDibaryon_Lam_ALam_side->Sumw2(); 
h_HDibaryon_Lam_ALam_peakside->Sumw2(); 


h_Hf2_K0s_K0s->Sumw2(); 
h_Hf2_rot_K0s_K0s->Sumw2(); 
h_Hf2_inv_K0s_K0s->Sumw2(); 
h_Hf2_K0s_K0s_mix->Sumw2(); 

h_Hf2_K0s_K0s_side->Sumw2(); 
h_Hf2_rot_K0s_K0s_side->Sumw2(); 
h_Hf2_inv_K0s_K0s_side->Sumw2(); 
h_Hf2_K0s_K0s_mix_side->Sumw2(); 

h_Hf2_K0s_K0s_peakside->Sumw2(); 
h_Hf2_rot_K0s_K0s_peakside->Sumw2(); 
h_Hf2_inv_K0s_K0s_peakside->Sumw2(); 
h_Hf2_K0s_K0s_mix_peakside->Sumw2(); 


h_Hcasc1820_K0s_Lam->Sumw2(); 
h_Hcasc1820_rot_K0s_Lam->Sumw2(); 
h_Hcasc1820_inv_K0s_Lam->Sumw2(); 
h_Hcasc1820_K0s_Lam_mix->Sumw2(); 

h_Hcasc1820_K0s_Lam_side->Sumw2(); 
h_Hcasc1820_rot_K0s_Lam_side->Sumw2(); 
h_Hcasc1820_inv_K0s_Lam_side->Sumw2(); 
h_Hcasc1820_K0s_Lam_mix_side->Sumw2(); 

h_Hcasc1820_K0s_Lam_peakside->Sumw2(); 
h_Hcasc1820_rot_K0s_Lam_peakside->Sumw2(); 
h_Hcasc1820_inv_K0s_Lam_peakside->Sumw2(); 
h_Hcasc1820_K0s_Lam_mix_peakside->Sumw2(); 

h_Hcasc1820_K0s_ALam->Sumw2(); 
h_Hcasc1820_rot_K0s_ALam->Sumw2(); 
h_Hcasc1820_inv_K0s_ALam->Sumw2(); 
h_Hcasc1820_K0s_ALam_mix->Sumw2(); 

h_Hcasc1820_K0s_ALam_side->Sumw2(); 
h_Hcasc1820_rot_K0s_ALam_side->Sumw2(); 
h_Hcasc1820_inv_K0s_ALam_side->Sumw2(); 
h_Hcasc1820_K0s_ALam_mix_side->Sumw2(); 

h_Hcasc1820_K0s_ALam_peakside->Sumw2(); 
h_Hcasc1820_rot_K0s_ALam_peakside->Sumw2(); 
h_Hcasc1820_inv_K0s_ALam_peakside->Sumw2(); 
h_Hcasc1820_K0s_ALam_mix_peakside->Sumw2(); 



h_jet_jet->Sumw2(); 
h_inv_jet_jet->Sumw2(); 
h_rot_jet_jet->Sumw2(); 
hMix_jet_jet->Sumw2(); 

//K0sK0s

h_K0s_K0s_samejet->Sumw2(); 
h_inv_K0s_K0s_samejet->Sumw2(); 
h_rot_K0s_K0s_samejet->Sumw2(); 
hMix_K0s_K0s_samejet->Sumw2(); 

h_K0s_K0s_diffjet->Sumw2(); 
h_inv_K0s_K0s_diffjet->Sumw2(); 
h_rot_K0s_K0s_diffjet->Sumw2(); 
hMix_K0s_K0s_diffjet->Sumw2(); 
  
h_K0s_K0s_nojet->Sumw2(); 
h_inv_K0s_K0s_nojet->Sumw2(); 
h_rot_K0s_K0s_nojet->Sumw2(); 
hMix_K0s_K0s_nojet->Sumw2(); 

h_K0s_K0s_only1jet->Sumw2(); 
h_inv_K0s_K0s_only1jet->Sumw2(); 
h_rot_K0s_K0s_only1jet->Sumw2(); 
hMix_K0s_K0s_only1jet->Sumw2(); 

//LamLam

h_Lam_Lam_samejet->Sumw2(); 
h_inv_Lam_Lam_samejet->Sumw2(); 
h_rot_Lam_Lam_samejet->Sumw2(); 
hMix_Lam_Lam_samejet->Sumw2(); 

h_Lam_Lam_diffjet->Sumw2(); 
h_inv_Lam_Lam_diffjet->Sumw2(); 
h_rot_Lam_Lam_diffjet->Sumw2(); 
hMix_Lam_Lam_diffjet->Sumw2(); 
  
h_Lam_Lam_nojet->Sumw2(); 
h_inv_Lam_Lam_nojet->Sumw2(); 
h_rot_Lam_Lam_nojet->Sumw2(); 
hMix_Lam_Lam_nojet->Sumw2(); 

h_Lam_Lam_only1jet->Sumw2();
h_inv_Lam_Lam_only1jet->Sumw2();
h_rot_Lam_Lam_only1jet->Sumw2();
hMix_Lam_Lam_only1jet->Sumw2();

//ALamALam

h_ALam_ALam_samejet->Sumw2();
h_inv_ALam_ALam_samejet->Sumw2();
h_rot_ALam_ALam_samejet->Sumw2();
hMix_ALam_ALam_samejet->Sumw2();

h_ALam_ALam_diffjet->Sumw2();
h_inv_ALam_ALam_diffjet->Sumw2();
h_rot_ALam_ALam_diffjet->Sumw2();
hMix_ALam_ALam_diffjet->Sumw2();
  
h_ALam_ALam_nojet->Sumw2();
h_inv_ALam_ALam_nojet->Sumw2();
h_rot_ALam_ALam_nojet->Sumw2();
hMix_ALam_ALam_nojet->Sumw2();

h_ALam_ALam_only1jet->Sumw2();
h_inv_ALam_ALam_only1jet->Sumw2();
h_rot_ALam_ALam_only1jet->Sumw2();
hMix_ALam_ALam_only1jet->Sumw2();
  
//LAL

h_Lam_ALam_samejet->Sumw2();
h_inv_Lam_ALam_samejet->Sumw2();
h_rot_Lam_ALam_samejet->Sumw2();
hMix_Lam_ALam_samejet->Sumw2();

h_Lam_ALam_diffjet->Sumw2();
h_inv_Lam_ALam_diffjet->Sumw2();
h_rot_Lam_ALam_diffjet->Sumw2();
hMix_Lam_ALam_diffjet->Sumw2();
  
h_Lam_ALam_nojet->Sumw2();
h_inv_Lam_ALam_nojet->Sumw2();
h_rot_Lam_ALam_nojet->Sumw2();
hMix_Lam_ALam_nojet->Sumw2(); 

h_Lam_ALam_only1jet->Sumw2();
h_inv_Lam_ALam_only1jet->Sumw2();
h_rot_Lam_ALam_only1jet->Sumw2();
hMix_Lam_ALam_only1jet->Sumw2(); 
  
//KL

h_K0s_Lam_samejet->Sumw2();
h_inv_K0s_Lam_samejet->Sumw2();
h_rot_K0s_Lam_samejet->Sumw2();
hMix_K0s_Lam_samejet->Sumw2();

h_K0s_Lam_diffjet->Sumw2();
h_inv_K0s_Lam_diffjet->Sumw2();
h_rot_K0s_Lam_diffjet->Sumw2();
hMix_K0s_Lam_diffjet->Sumw2();
  
h_K0s_Lam_nojet->Sumw2();
h_inv_K0s_Lam_nojet->Sumw2();
h_rot_K0s_Lam_nojet->Sumw2();
hMix_K0s_Lam_nojet->Sumw2(); 

h_K0s_Lam_only1jet->Sumw2();
h_inv_K0s_Lam_only1jet->Sumw2();
h_rot_K0s_Lam_only1jet->Sumw2();
hMix_K0s_Lam_only1jet->Sumw2();

//KAL

h_K0s_ALam_samejet->Sumw2();
h_inv_K0s_ALam_samejet->Sumw2();
h_rot_K0s_ALam_samejet->Sumw2();
hMix_K0s_ALam_samejet->Sumw2(); 

h_K0s_ALam_diffjet->Sumw2();
h_inv_K0s_ALam_diffjet->Sumw2();
h_rot_K0s_ALam_diffjet->Sumw2();
hMix_K0s_ALam_diffjet->Sumw2();
  
h_K0s_ALam_nojet->Sumw2();
h_inv_K0s_ALam_nojet->Sumw2();
h_rot_K0s_ALam_nojet->Sumw2();
hMix_K0s_ALam_nojet->Sumw2(); 

h_K0s_ALam_only1jet->Sumw2();
h_inv_K0s_ALam_only1jet->Sumw2();
h_rot_K0s_ALam_only1jet->Sumw2();
hMix_K0s_ALam_only1jet->Sumw2(); 

// hist for V0 + jet  

h_K0snojet_Jet->Sumw2();
h_K0snojet_Jet_rot->Sumw2();
h_K0snojet_Jet_inv->Sumw2();
hMix_K0snojet_Jet->Sumw2(); 

h_Lamnojet_Jet->Sumw2();
h_Lamnojet_Jet_rot->Sumw2();
h_Lamnojet_Jet_inv->Sumw2();
hMix_Lamnojet_Jet->Sumw2(); 

h_ALamnojet_Jet->Sumw2();
h_ALamnojet_Jet_rot->Sumw2();
h_ALamnojet_Jet_inv->Sumw2();
hMix_ALamnojet_Jet->Sumw2();

h_K0sfromjet_Jet->Sumw2(); 
h_K0sfromjet_Jet_rot->Sumw2();
h_K0sfromjet_Jet_inv->Sumw2();
hMix_K0sfromjet_Jet->Sumw2(); 

h_Lamfromjet_Jet->Sumw2();
h_Lamfromjet_Jet_rot->Sumw2();
h_Lamfromjet_Jet_inv->Sumw2();
hMix_Lamfromjet_Jet->Sumw2();

h_ALamfromjet_Jet->Sumw2();
h_ALamfromjet_Jet_rot->Sumw2();
h_ALamfromjet_Jet_inv->Sumw2();
hMix_ALamfromjet_Jet->Sumw2();

h_K0s_Jet->Sumw2();
h_K0s_Jet_rot->Sumw2();
h_K0s_Jet_inv->Sumw2();
hMix_K0s_Jet->Sumw2(); 

h_Lam_Jet->Sumw2(); 
h_Lam_Jet_rot->Sumw2();
h_Lam_Jet_inv->Sumw2();
hMix_Lam_Jet->Sumw2();

h_ALam_Jet->Sumw2();
h_ALam_Jet_rot->Sumw2();
h_ALam_Jet_inv->Sumw2();
hMix_ALam_Jet->Sumw2();


K_size->Sumw2();
L_size->Sumw2();
AL_size->Sumw2();

K_size_nocut->Sumw2();
L_size_nocut->Sumw2();
AL_size_nocut->Sumw2();
 
APwocut_K0s->Sumw2();
APwcut_K0s->Sumw2();
APwocut_LAL->Sumw2();
APwcut_LAL->Sumw2();
 
K0smiss_mass_hyp->Sumw2();
K0smiss_mass_ee->Sumw2();
LALmiss_mass_hyp->Sumw2();
LALmiss_mass_ee->Sumw2();
 
pt_mass_K->Sumw2();
pt_mass_L->Sumw2();
pt_mass_AL->Sumw2();

K0s_pt_reco->Sumw2();
K0s_eta_reco->Sumw2();
K0s_phi_reco->Sumw2();
K0s_mass_reco->Sumw2();

L_pt_reco->Sumw2();
L_eta_reco->Sumw2();
L_phi_reco->Sumw2();
L_mass_reco->Sumw2();

AL_pt_reco->Sumw2();
AL_eta_reco->Sumw2();
AL_phi_reco->Sumw2();
AL_mass_reco->Sumw2();

K0s_pt_reco_d1->Sumw2();
K0s_eta_reco_d1->Sumw2();
K0s_phi_reco_d1->Sumw2();
K0s_mass_reco_d1->Sumw2();

L_pt_reco_d1->Sumw2();
L_eta_reco_d1->Sumw2();
L_phi_reco_d1->Sumw2();
L_mass_reco_d1->Sumw2();

AL_pt_reco_d1->Sumw2();
AL_eta_reco_d1->Sumw2();
AL_phi_reco_d1->Sumw2();
AL_mass_reco_d1->Sumw2();

K0s_pt_reco_d2->Sumw2();
K0s_eta_reco_d2->Sumw2();
K0s_phi_reco_d2->Sumw2();
K0s_mass_reco_d2->Sumw2();

L_pt_reco_d2->Sumw2();
L_eta_reco_d2->Sumw2();
L_phi_reco_d2->Sumw2();
L_mass_reco_d2->Sumw2();

AL_pt_reco_d2->Sumw2();
AL_eta_reco_d2->Sumw2();
AL_phi_reco_d2->Sumw2();
AL_mass_reco_d2->Sumw2();

K0s_pt_gen->Sumw2();
K0s_eta_gen->Sumw2();
K0s_phi_gen->Sumw2();
K0s_mass_gen->Sumw2();

L_pt_gen->Sumw2();
L_eta_gen->Sumw2();
L_phi_gen->Sumw2();
L_mass_gen->Sumw2();

AL_pt_gen->Sumw2();
AL_eta_gen->Sumw2();
AL_phi_gen->Sumw2();
AL_mass_gen->Sumw2();

K0s_pt_mat->Sumw2();
K0s_eta_mat->Sumw2();
K0s_phi_mat->Sumw2();
K0s_mass_mat->Sumw2();

L_pt_mat->Sumw2();
L_eta_mat->Sumw2();
L_phi_mat->Sumw2();
L_mass_mat->Sumw2();

AL_pt_mat->Sumw2();
AL_eta_mat->Sumw2();
AL_phi_mat->Sumw2();
AL_mass_mat->Sumw2();
 


}

void write_info_histos(){

//events
  
nocutev->Write(); 

nocutPV->Write();  
nocutPU->Write();  
nocutHF->Write(); 
nocutSC->Write(); 
nocutALL->Write();
nocutALLX->Write();


nev_K0s_ini->Write(); 
nev_Lam_ini->Write(); 
nev_ALam_ini->Write(); 

nev_K0s_AS->Write(); 
nev_Lam_AS->Write(); 
nev_ALam_AS->Write(); 
  
//K0s

nev_K0sT->Write();
nev_K0sF->Write();
nev_K0sTF->Write();
nev_K0s_ssT->Write(); 
nev_K0s_bbT->Write();
nev_K0s_sbT->Write();
nev_K0s_ssT_mix->Write();
nev_K0s_bbT_mix->Write();
nev_K0s_sbT_mix->Write();
nev_K0s_ssT_etamix->Write();
nev_K0s_bbT_etamix->Write();
nev_K0s_sbT_etamix->Write();

  
//Lam

nev_LamT->Write();
nev_LamF->Write();
nev_LamTF->Write();
nev_Lam_ssT->Write();
nev_Lam_bbT->Write();
nev_Lam_sbT->Write();
nev_Lam_ssT_mix->Write();
nev_Lam_bbT_mix->Write();
nev_Lam_sbT_mix->Write();
nev_Lam_ssT_etamix->Write();
nev_Lam_bbT_etamix->Write();
nev_Lam_sbT_etamix->Write();

  
  
//ALam

nev_ALamT->Write();
nev_ALamF->Write();
nev_ALamTF->Write();
nev_ALam_ssT->Write();
nev_ALam_bbT->Write();
nev_ALam_sbT->Write();
nev_ALam_ssT_mix->Write();
nev_ALam_bbT_mix->Write();
nev_ALam_sbT_mix->Write();
nev_ALam_ssT_etamix->Write();
nev_ALam_bbT_etamix->Write();
nev_ALam_sbT_etamix->Write();

  
  
//LAL
nev_LALT->Write();
nev_LALF->Write();
nev_LALTF->Write();
nev_LAL_ssT->Write();
nev_LAL_bbT->Write();
nev_LAL_sbT->Write();
nev_LAL_ssT_mix->Write();
nev_LAL_bbT_mix->Write();
nev_LAL_sbT_mix->Write();
nev_LAL_ssT_etamix->Write();
nev_LAL_bbT_etamix->Write();
nev_LAL_sbT_etamix->Write();

//KL

nev_KLT->Write();
nev_KLF->Write(); 
nev_KLTF->Write();
nev_KL_ssT->Write(); 
nev_KL_bbT->Write();
nev_KL_sbT->Write(); 
nev_KL_ssT_mix->Write();
nev_KL_bbT_mix->Write();
nev_KL_sbT_mix->Write();
nev_KL_ssT_etamix->Write();
nev_KL_bbT_etamix->Write();
nev_KL_sbT_etamix->Write();

//KAL

nev_KALT->Write();
nev_KALF->Write();
nev_KALTF->Write();
nev_KAL_ssT->Write();
nev_KAL_bbT->Write();
nev_KAL_sbT->Write();
nev_KAL_ssT_mix->Write();
nev_KAL_bbT_mix->Write();
nev_KAL_sbT_mix->Write();
nev_KAL_ssT_etamix->Write();
nev_KAL_bbT_etamix->Write();
nev_KAL_sbT_etamix->Write();

cone_K0s->Write();
cone_K0sT->Write();
cone_K0sD1->Write();
cone_K0sTD1->Write();
cone_K0sD2->Write();
cone_K0sTD2->Write();

cone_Lam->Write();
cone_LamT->Write();
cone_LamD1->Write();
cone_LamTD1->Write();
cone_LamD2->Write();
cone_LamTD2->Write();

cone_Jet->Write();

  //vertex / ntrk / centrality  
h_vtx_z->Write();  
h_vtx_rho->Write();
h_ntrk_cent->Write();
Ntrk_VZ->Write();

h_ntrk_cent_1Ks->Write();
h_ntrk_cent_1Lam->Write();
h_ntrk_cent_1ALam->Write();

//mass
K0s_Mass->Write();
Lam_Mass->Write();
ALam_Mass->Write();
LAL_Mass->Write();
XXi_mass->Write();
AXXi_mass->Write();
OOm_mass->Write();
AOOm_mass->Write();
XXi_mass_LambMom->Write();
AXXi_mass_LambMom->Write();
OOm_mass_LambMom->Write();
AOOm_mass_LambMom->Write();

K0sctau->Write();
K0sdca3D->Write();
K0sDL->Write();
K0strkdca->Write();
K0scostheta->Write();
K0s_Vtx->Write();


Lamctau->Write();
Lamdca3D->Write();
LamDL->Write();
Lamtrkdca->Write();
Lamcostheta->Write();
Lam_Vtx->Write();


ALamctau->Write();
ALamdca3D->Write();
ALamDL->Write();
ALamtrkdca->Write();
ALamcostheta->Write();
ALam_Vtx->Write();


K0s_d1_dxy->Write();
K0s_d2_dxy->Write();
K0s_d1_dz->Write();
K0s_d2_dz->Write();
Lam_d1_dxy->Write();
Lam_d2_dxy->Write();
Lam_d1_dz->Write();
Lam_d2_dz->Write();
ALam_d1_dxy->Write();
ALam_d2_dxy->Write();
ALam_d1_dz->Write();
ALam_d2_dz->Write();

K0s_d1_Nhits->Write();
K0s_d2_Nhits->Write();
Lam_d1_Nhits->Write();
Lam_d2_Nhits->Write();
ALam_d1_Nhits->Write();
ALam_d2_Nhits->Write();

DR_K0s->Write();
DR_Lam->Write();
DR_ALam->Write();

DPT_K0s->Write();
DPT_Lam->Write();
DPT_ALam->Write();

pT_Lamb_from_Xi_T->Write();
pT_Lamb_from_Om_T->Write();
pT_Lamb_from_prompt_T->Write();

pT_Lamb_from_Xi_TF->Write();
pT_Lamb_from_Om_TF->Write();
pT_Lamb_from_prompt_TF->Write();

pT_Lamb_daughter_Xi->Write();
pT_Lamb_daughter_Om->Write();
pT_Lamb_daughter_prompt->Write();


hGen_K0s_K0s_Matched->Write();
hGen_K0s_K0s_Matched_mix->Write();

hGen_Lam_Lam_Matched->Write();
hGen_Lam_Lam_Matched_mix->Write();

hGen_ALam_ALam_Matched->Write();
hGen_ALam_ALam_Matched_mix->Write();

hGen_Lam_ALam_Matched->Write();
hGen_Lam_ALam_Matched_mix->Write();

hGen_K0s_Lam_Matched->Write();
hGen_K0s_Lam_Matched_mix->Write();

hGen_K0s_ALam_Matched->Write();
hGen_K0s_ALam_Matched_mix->Write();


hGen_K0s_K0s_MatchedT->Write();
hGen_K0s_K0s_MatchedT_mix->Write();

hGen_Lam_Lam_MatchedT->Write();
hGen_Lam_Lam_MatchedT_mix->Write();

hGen_ALam_ALam_MatchedT->Write();
hGen_ALam_ALam_MatchedT_mix->Write();

hGen_Lam_ALam_MatchedT->Write();
hGen_Lam_ALam_MatchedT_mix->Write();

hGen_K0s_Lam_MatchedT->Write();
hGen_K0s_Lam_MatchedT_mix->Write();

hGen_K0s_ALam_MatchedT->Write();
hGen_K0s_ALam_MatchedT_mix->Write();

pT_Xi->Write();
pT_Om->Write();
pT_Lamb->Write();

} //DONE

void write_K0s_histos(){
    
//K0s
//True

hMass_K0s_K0s_T->Write();   
hMass_K0s_1_T->Write();   
hMass_K0s_2_T->Write();   
  
hpeak_K0s_K0s_T->Write();
hpeakcos_K0s_K0s_T->Write();
hpeak_rot_K0s_K0s_T->Write();
hpeak_inv_K0s_K0s_T->Write();
hside_K0s_K0s_T->Write();
hsideL_K0s_K0s_T->Write();
hsideR_K0s_K0s_T->Write();
hside_rot_K0s_K0s_T->Write();
hside_inv_K0s_K0s_T->Write();
hpeakside_K0s_K0s_T->Write();
hpeaksideL_K0s_K0s_T->Write();
hpeaksideR_K0s_K0s_T->Write();
hpeakside_rot_K0s_K0s_T->Write();
hpeakside_inv_K0s_K0s_T->Write();
hpeakd1_K0s_K0s_T->Write();
hpeakd2_K0s_K0s_T->Write();
hpeakd12_K0s_K0s_T->Write();

hpeak_K0s_K0s_T_mix->Write();
hside_K0s_K0s_T_mix->Write();
hsideL_K0s_K0s_T_mix->Write();
hsideR_K0s_K0s_T_mix->Write();
hpeakside_K0s_K0s_T_mix->Write();
hpeaksideL_K0s_K0s_T_mix->Write();
hpeaksideR_K0s_K0s_T_mix->Write();

hpeak_K0s_K0s_T_etamix->Write();
hside_K0s_K0s_T_etamix->Write();
hsideL_K0s_K0s_T_etamix->Write();
hsideR_K0s_K0s_T_etamix->Write();
hpeakside_K0s_K0s_T_etamix->Write();
hpeaksideL_K0s_K0s_T_etamix->Write();
hpeaksideR_K0s_K0s_T_etamix->Write();


V0chi2d1_diff_T_K0s->Write();
V0chi2d2_diff_T_K0s->Write();
V0chi2d12_diff_T_K0s->Write();

//False

hMass_K0s_K0s_F->Write();
  
hpeak_K0s_K0s_F->Write();
hside_K0s_K0s_F->Write();
hpeakside_K0s_K0s_F->Write();
hpeakd1_K0s_K0s_F->Write();
hpeakd2_K0s_K0s_F->Write();
hpeakd12_K0s_K0s_F->Write();
hpeak_K0s_K0s_F_bias->Write();
hpeak_K0s_K0s_F_biasd1->Write();
hpeak_K0s_K0s_F_biasd2->Write();
hpeak_K0s_K0s_F_biasd12->Write();

V0chi2d1_diff_F_K0s->Write();
V0chi2d2_diff_F_K0s->Write();
V0chi2d12_diff_F_K0s->Write();
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
hpeak_K0s_K0s_TF->Write();
hpeakd1_K0s_K0s_TF->Write();
hpeakd2_K0s_K0s_TF->Write();
hpeakd12_K0s_K0s_TF->Write();
hside_K0s_K0s_TF->Write();
hpeakside_K0s_K0s_TF->Write();
hpeak_K0s_K0s_TF_bias->Write();
hpeak_K0s_K0s_TF_biasd1->Write();
hpeak_K0s_K0s_TF_biasd2->Write();
hpeak_K0s_K0s_TF_biasd12->Write();

V0chi2d1_diff_TF_K0s->Write();
V0chi2d2_diff_TF_K0s->Write();
V0chi2d12_diff_TF_K0s->Write();

  
/////// matching hist ////// 
  
hMass_K0s_K0s_TM->Write();
hMass_K0s_K0s_TU->Write();
hMass_K0s_K0s_FM->Write();
hMass_K0s_K0s_FU->Write();
 
hpeak_K0s_K0s_TM->Write();
hpeakd1_K0s_K0s_TM->Write();
hpeakd2_K0s_K0s_TM->Write();
hpeakd12_K0s_K0s_TM->Write();
hpeak_K0s_K0s_TM_bias->Write();
hpeak_K0s_K0s_TM_biasd1->Write();
hpeak_K0s_K0s_TM_biasd2->Write();
hpeak_K0s_K0s_TM_biasd12->Write();

V0chi2d1_diff_TM_K0s->Write();
V0chi2d2_diff_TM_K0s->Write();
V0chi2d12_diff_TM_K0s->Write();

hpeak_K0s_K0s_TU->Write();
hpeakd1_K0s_K0s_TU->Write();
hpeakd2_K0s_K0s_TU->Write();
hpeakd12_K0s_K0s_TU->Write();
hpeak_K0s_K0s_TU_bias->Write();
hpeak_K0s_K0s_TU_biasd1->Write();
hpeak_K0s_K0s_TU_biasd2->Write();
hpeak_K0s_K0s_TU_biasd12->Write();

V0chi2d1_diff_TU_K0s->Write();
V0chi2d2_diff_TU_K0s->Write();
V0chi2d12_diff_TU_K0s->Write();

hpeak_K0s_K0s_FM->Write();
hpeakd1_K0s_K0s_FM->Write();
hpeakd2_K0s_K0s_FM->Write();
hpeakd12_K0s_K0s_FM->Write();
hpeak_K0s_K0s_FM_bias->Write();
hpeak_K0s_K0s_FM_biasd1->Write();
hpeak_K0s_K0s_FM_biasd2->Write();
hpeak_K0s_K0s_FM_biasd12->Write();

V0chi2d1_diff_FM_K0s->Write();
V0chi2d2_diff_FM_K0s->Write();
V0chi2d12_diff_FM_K0s->Write();

hpeak_K0s_K0s_FU->Write();
hpeakd1_K0s_K0s_FU->Write();
hpeakd2_K0s_K0s_FU->Write();
hpeakd12_K0s_K0s_FU->Write();
hpeak_K0s_K0s_FU_bias->Write();
hpeak_K0s_K0s_FU_biasd1->Write();
hpeak_K0s_K0s_FU_biasd2->Write();
hpeak_K0s_K0s_FU_biasd12->Write();

V0chi2d1_diff_FU_K0s->Write();
V0chi2d2_diff_FU_K0s->Write();
V0chi2d12_diff_FU_K0s->Write();


//true matched + true unmatched

hpeak_K0s_K0s_TM_TU->Write();
hpeakd1_K0s_K0s_TM_TU->Write();
hpeakd2_K0s_K0s_TM_TU->Write();
hpeakd12_K0s_K0s_TM_TU->Write();
hpeak_K0s_K0s_TM_TU_bias->Write();
hpeak_K0s_K0s_TM_TU_biasd1->Write();
hpeak_K0s_K0s_TM_TU_biasd2->Write();
hpeak_K0s_K0s_TM_TU_biasd12->Write();

V0chi2d1_diff_TM_TU_K0s->Write();
V0chi2d2_diff_TM_TU_K0s->Write();
V0chi2d12_diff_TM_TU_K0s->Write();

//fake matched + fake unmatched

hpeak_K0s_K0s_FM_FU->Write();
hpeakd1_K0s_K0s_FM_FU->Write();
hpeakd2_K0s_K0s_FM_FU->Write();
hpeakd12_K0s_K0s_FM_FU->Write();
hpeak_K0s_K0s_FM_FU_bias->Write();
hpeak_K0s_K0s_FM_FU_biasd1->Write();
hpeak_K0s_K0s_FM_FU_biasd2->Write();
hpeak_K0s_K0s_FM_FU_biasd12->Write();

V0chi2d1_diff_FM_FU_K0s->Write();
V0chi2d2_diff_FM_FU_K0s->Write();
V0chi2d12_diff_FM_FU_K0s->Write();


hpeak_K0s_K0s_TM_FM->Write();
hpeakd1_K0s_K0s_TM_FM->Write();
hpeakd2_K0s_K0s_TM_FM->Write();
hpeakd12_K0s_K0s_TM_FM->Write();
hpeak_K0s_K0s_TM_FM_bias->Write();
hpeak_K0s_K0s_TM_FM_biasd1->Write();
hpeak_K0s_K0s_TM_FM_biasd2->Write();
hpeak_K0s_K0s_TM_FM_biasd12->Write();

V0chi2d1_diff_TM_FM_K0s->Write();
V0chi2d2_diff_TM_FM_K0s->Write();
V0chi2d12_diff_TM_FM_K0s->Write();

hpeak_K0s_K0s_TM_FU->Write();
hpeakd1_K0s_K0s_TM_FU->Write();
hpeakd2_K0s_K0s_TM_FU->Write();
hpeakd12_K0s_K0s_TM_FU->Write();
hpeak_K0s_K0s_TM_FU_bias->Write();
hpeak_K0s_K0s_TM_FU_biasd1->Write();
hpeak_K0s_K0s_TM_FU_biasd2->Write();
hpeak_K0s_K0s_TM_FU_biasd12->Write();

V0chi2d1_diff_TM_FU_K0s->Write();
V0chi2d2_diff_TM_FU_K0s->Write();
V0chi2d12_diff_TM_FU_K0s->Write();

hpeak_K0s_K0s_TU_FU->Write();
hpeakd1_K0s_K0s_TU_FU->Write();
hpeakd2_K0s_K0s_TU_FU->Write();
hpeakd12_K0s_K0s_TU_FU->Write();
hpeak_K0s_K0s_TU_FU_bias->Write();
hpeak_K0s_K0s_TU_FU_biasd1->Write();
hpeak_K0s_K0s_TU_FU_biasd2->Write();
hpeak_K0s_K0s_TU_FU_biasd12->Write();

V0chi2d1_diff_TU_FU_K0s->Write();
V0chi2d2_diff_TU_FU_K0s->Write();
V0chi2d12_diff_TU_FU_K0s->Write();

hpeak_K0s_K0s_TU_FM->Write();
hpeakd1_K0s_K0s_TU_FM->Write();
hpeakd2_K0s_K0s_TU_FM->Write();
hpeakd12_K0s_K0s_TU_FM->Write();
hpeak_K0s_K0s_TU_FM_bias->Write();
hpeak_K0s_K0s_TU_FM_biasd1->Write();
hpeak_K0s_K0s_TU_FM_biasd2->Write();
hpeak_K0s_K0s_TU_FM_biasd12->Write();

V0chi2d1_diff_TU_FM_K0s->Write();
V0chi2d2_diff_TU_FM_K0s->Write();
V0chi2d12_diff_TU_FM_K0s->Write();
    
    
} //DONE

void write_Lam_histos(){

//Lam
//True

hMass_Lam_Lam_T->Write();   
hMass_Lam_1_T->Write();   
hMass_Lam_2_T->Write();   

  
hpeak_Lam_Lam_T->Write();
hpeakcos_Lam_Lam_T->Write();
hpeak_rot_Lam_Lam_T->Write();
hpeak_inv_Lam_Lam_T->Write();
hside_Lam_Lam_T->Write();
hsideL_Lam_Lam_T->Write();
hsideR_Lam_Lam_T->Write();
hside_rot_Lam_Lam_T->Write();
hside_inv_Lam_Lam_T->Write();
hpeakside_Lam_Lam_T->Write();
hpeaksideL_Lam_Lam_T->Write();
hpeaksideR_Lam_Lam_T->Write();
hpeakside_rot_Lam_Lam_T->Write();
hpeakside_inv_Lam_Lam_T->Write();
hpeakd1_Lam_Lam_T->Write();
hpeakd2_Lam_Lam_T->Write();
hpeakd12_Lam_Lam_T->Write();

hpeak_Lam_Lam_T_mix->Write();
hside_Lam_Lam_T_mix->Write();
hsideL_Lam_Lam_T_mix->Write();
hsideR_Lam_Lam_T_mix->Write();
hpeakside_Lam_Lam_T_mix->Write();
hpeaksideL_Lam_Lam_T_mix->Write();
hpeaksideR_Lam_Lam_T_mix->Write();

hpeak_Lam_Lam_T_etamix->Write();
hside_Lam_Lam_T_etamix->Write();
hsideL_Lam_Lam_T_etamix->Write();
hsideR_Lam_Lam_T_etamix->Write();
hpeakside_Lam_Lam_T_etamix->Write();
hpeaksideL_Lam_Lam_T_etamix->Write();
hpeaksideR_Lam_Lam_T_etamix->Write();


V0chi2d1_diff_T_Lam->Write();
V0chi2d2_diff_T_Lam->Write();
V0chi2d12_diff_T_Lam->Write();

//False

hMass_Lam_Lam_F->Write();
  
hpeak_Lam_Lam_F->Write();
hside_Lam_Lam_F->Write();
hpeakside_Lam_Lam_F->Write();
hpeakd1_Lam_Lam_F->Write();
hpeakd2_Lam_Lam_F->Write();
hpeakd12_Lam_Lam_F->Write();
hpeak_Lam_Lam_F_bias->Write();
hpeak_Lam_Lam_F_biasd1->Write();
hpeak_Lam_Lam_F_biasd2->Write();
hpeak_Lam_Lam_F_biasd12->Write();

V0chi2d1_diff_F_Lam->Write();
V0chi2d2_diff_F_Lam->Write();
V0chi2d12_diff_F_Lam->Write();
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
hpeak_Lam_Lam_TF->Write();
hpeakd1_Lam_Lam_TF->Write();
hpeakd2_Lam_Lam_TF->Write();
hpeakd12_Lam_Lam_TF->Write();
hside_Lam_Lam_TF->Write();
hpeakside_Lam_Lam_TF->Write();
hpeak_Lam_Lam_TF_bias->Write();
hpeak_Lam_Lam_TF_biasd1->Write();
hpeak_Lam_Lam_TF_biasd2->Write();
hpeak_Lam_Lam_TF_biasd12->Write();

V0chi2d1_diff_TF_Lam->Write();
V0chi2d2_diff_TF_Lam->Write();
V0chi2d12_diff_TF_Lam->Write();

/////// matching hist ////// 
  
hMass_Lam_Lam_TM->Write();
hMass_Lam_Lam_TU->Write();
hMass_Lam_Lam_FM->Write();
hMass_Lam_Lam_FU->Write();
 
hpeak_Lam_Lam_TM->Write();
hpeakd1_Lam_Lam_TM->Write();
hpeakd2_Lam_Lam_TM->Write();
hpeakd12_Lam_Lam_TM->Write();
hpeak_Lam_Lam_TM_bias->Write();
hpeak_Lam_Lam_TM_biasd1->Write();
hpeak_Lam_Lam_TM_biasd2->Write();
hpeak_Lam_Lam_TM_biasd12->Write();

V0chi2d1_diff_TM_Lam->Write();
V0chi2d2_diff_TM_Lam->Write();
V0chi2d12_diff_TM_Lam->Write();

hpeak_Lam_Lam_TU->Write();
hpeakd1_Lam_Lam_TU->Write();
hpeakd2_Lam_Lam_TU->Write();
hpeakd12_Lam_Lam_TU->Write();
hpeak_Lam_Lam_TU_bias->Write();
hpeak_Lam_Lam_TU_biasd1->Write();
hpeak_Lam_Lam_TU_biasd2->Write();
hpeak_Lam_Lam_TU_biasd12->Write();

V0chi2d1_diff_TU_Lam->Write();
V0chi2d2_diff_TU_Lam->Write();
V0chi2d12_diff_TU_Lam->Write();

hpeak_Lam_Lam_FM->Write();
hpeakd1_Lam_Lam_FM->Write();
hpeakd2_Lam_Lam_FM->Write();
hpeakd12_Lam_Lam_FM->Write();
hpeak_Lam_Lam_FM_bias->Write();
hpeak_Lam_Lam_FM_biasd1->Write();
hpeak_Lam_Lam_FM_biasd2->Write();
hpeak_Lam_Lam_FM_biasd12->Write();

V0chi2d1_diff_FM_Lam->Write();
V0chi2d2_diff_FM_Lam->Write();
V0chi2d12_diff_FM_Lam->Write();

hpeak_Lam_Lam_FU->Write();
hpeakd1_Lam_Lam_FU->Write();
hpeakd2_Lam_Lam_FU->Write();
hpeakd12_Lam_Lam_FU->Write();
hpeak_Lam_Lam_FU_bias->Write();
hpeak_Lam_Lam_FU_biasd1->Write();
hpeak_Lam_Lam_FU_biasd2->Write();
hpeak_Lam_Lam_FU_biasd12->Write();

V0chi2d1_diff_FU_Lam->Write();
V0chi2d2_diff_FU_Lam->Write();
V0chi2d12_diff_FU_Lam->Write();


//true matched + true unmatched

hpeak_Lam_Lam_TM_TU->Write();
hpeakd1_Lam_Lam_TM_TU->Write();
hpeakd2_Lam_Lam_TM_TU->Write();
hpeakd12_Lam_Lam_TM_TU->Write();
hpeak_Lam_Lam_TM_TU_bias->Write();
hpeak_Lam_Lam_TM_TU_biasd1->Write();
hpeak_Lam_Lam_TM_TU_biasd2->Write();
hpeak_Lam_Lam_TM_TU_biasd12->Write();

V0chi2d1_diff_TM_TU_Lam->Write();
V0chi2d2_diff_TM_TU_Lam->Write();
V0chi2d12_diff_TM_TU_Lam->Write();

//fake matched + fake unmatched

hpeak_Lam_Lam_FM_FU->Write();
hpeakd1_Lam_Lam_FM_FU->Write();
hpeakd2_Lam_Lam_FM_FU->Write();
hpeakd12_Lam_Lam_FM_FU->Write();
hpeak_Lam_Lam_FM_FU_bias->Write();
hpeak_Lam_Lam_FM_FU_biasd1->Write();
hpeak_Lam_Lam_FM_FU_biasd2->Write();
hpeak_Lam_Lam_FM_FU_biasd12->Write();

V0chi2d1_diff_FM_FU_Lam->Write();
V0chi2d2_diff_FM_FU_Lam->Write();
V0chi2d12_diff_FM_FU_Lam->Write();



hpeak_Lam_Lam_TM_FM->Write();
hpeakd1_Lam_Lam_TM_FM->Write();
hpeakd2_Lam_Lam_TM_FM->Write();
hpeakd12_Lam_Lam_TM_FM->Write();
hpeak_Lam_Lam_TM_FM_bias->Write();
hpeak_Lam_Lam_TM_FM_biasd1->Write();
hpeak_Lam_Lam_TM_FM_biasd2->Write();
hpeak_Lam_Lam_TM_FM_biasd12->Write();

V0chi2d1_diff_TM_FM_Lam->Write();
V0chi2d2_diff_TM_FM_Lam->Write();
V0chi2d12_diff_TM_FM_Lam->Write();

hpeak_Lam_Lam_TM_FU->Write();
hpeakd1_Lam_Lam_TM_FU->Write();
hpeakd2_Lam_Lam_TM_FU->Write();
hpeakd12_Lam_Lam_TM_FU->Write();
hpeak_Lam_Lam_TM_FU_bias->Write();
hpeak_Lam_Lam_TM_FU_biasd1->Write();
hpeak_Lam_Lam_TM_FU_biasd2->Write();
hpeak_Lam_Lam_TM_FU_biasd12->Write();

V0chi2d1_diff_TM_FU_Lam->Write();
V0chi2d2_diff_TM_FU_Lam->Write();
V0chi2d12_diff_TM_FU_Lam->Write();

hpeak_Lam_Lam_TU_FU->Write();
hpeakd1_Lam_Lam_TU_FU->Write();
hpeakd2_Lam_Lam_TU_FU->Write();
hpeakd12_Lam_Lam_TU_FU->Write();
hpeak_Lam_Lam_TU_FU_bias->Write();
hpeak_Lam_Lam_TU_FU_biasd1->Write();
hpeak_Lam_Lam_TU_FU_biasd2->Write();
hpeak_Lam_Lam_TU_FU_biasd12->Write();

V0chi2d1_diff_TU_FU_Lam->Write();
V0chi2d2_diff_TU_FU_Lam->Write();
V0chi2d12_diff_TU_FU_Lam->Write();

hpeak_Lam_Lam_TU_FM->Write();
hpeakd1_Lam_Lam_TU_FM->Write();
hpeakd2_Lam_Lam_TU_FM->Write();
hpeakd12_Lam_Lam_TU_FM->Write();
hpeak_Lam_Lam_TU_FM_bias->Write();
hpeak_Lam_Lam_TU_FM_biasd1->Write();
hpeak_Lam_Lam_TU_FM_biasd2->Write();
hpeak_Lam_Lam_TU_FM_biasd12->Write();

V0chi2d1_diff_TU_FM_Lam->Write();
V0chi2d2_diff_TU_FM_Lam->Write();
V0chi2d12_diff_TU_FM_Lam->Write();
    
    
} //DONE

void write_ALam_histos(){
    
//ALam
//True

hMass_ALam_ALam_T->Write();   
hMass_ALam_1_T->Write();   
hMass_ALam_2_T->Write();   
  
  
hpeak_ALam_ALam_T->Write();
hpeakcos_ALam_ALam_T->Write();
hpeak_rot_ALam_ALam_T->Write();
hpeak_inv_ALam_ALam_T->Write();
hside_ALam_ALam_T->Write();
hsideL_ALam_ALam_T->Write();
hsideR_ALam_ALam_T->Write();
hside_rot_ALam_ALam_T->Write();
hside_inv_ALam_ALam_T->Write();
hpeakside_ALam_ALam_T->Write();
hpeaksideL_ALam_ALam_T->Write();
hpeaksideR_ALam_ALam_T->Write();
hpeakside_rot_ALam_ALam_T->Write();
hpeakside_inv_ALam_ALam_T->Write();
hpeakd1_ALam_ALam_T->Write();
hpeakd2_ALam_ALam_T->Write();
hpeakd12_ALam_ALam_T->Write();

hpeak_ALam_ALam_T_mix->Write();
hside_ALam_ALam_T_mix->Write();
hsideL_ALam_ALam_T_mix->Write();
hsideR_ALam_ALam_T_mix->Write();
hpeakside_ALam_ALam_T_mix->Write();
hpeaksideL_ALam_ALam_T_mix->Write();
hpeaksideR_ALam_ALam_T_mix->Write();

hpeak_ALam_ALam_T_etamix->Write();
hside_ALam_ALam_T_etamix->Write();
hsideL_ALam_ALam_T_etamix->Write();
hsideR_ALam_ALam_T_etamix->Write();
hpeakside_ALam_ALam_T_etamix->Write();
hpeaksideL_ALam_ALam_T_etamix->Write();
hpeaksideR_ALam_ALam_T_etamix->Write();


V0chi2d1_diff_T_ALam->Write();
V0chi2d2_diff_T_ALam->Write();
V0chi2d12_diff_T_ALam->Write();
//False

hMass_ALam_ALam_F->Write();
  
hpeak_ALam_ALam_F->Write();
hside_ALam_ALam_F->Write();
hpeakside_ALam_ALam_F->Write();
hpeakd1_ALam_ALam_F->Write();
hpeakd2_ALam_ALam_F->Write();
hpeakd12_ALam_ALam_F->Write();
hpeak_ALam_ALam_F_bias->Write();
hpeak_ALam_ALam_F_biasd1->Write();
hpeak_ALam_ALam_F_biasd2->Write();
hpeak_ALam_ALam_F_biasd12->Write();

V0chi2d1_diff_F_ALam->Write();
V0chi2d2_diff_F_ALam->Write();
V0chi2d12_diff_F_ALam->Write();
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
hpeak_ALam_ALam_TF->Write();
hpeakd1_ALam_ALam_TF->Write();
hpeakd2_ALam_ALam_TF->Write();
hpeakd12_ALam_ALam_TF->Write();
hside_ALam_ALam_TF->Write();
hpeakside_ALam_ALam_TF->Write();
hpeak_ALam_ALam_TF_bias->Write();
hpeak_ALam_ALam_TF_biasd1->Write();
hpeak_ALam_ALam_TF_biasd2->Write();
hpeak_ALam_ALam_TF_biasd12->Write();

V0chi2d1_diff_TF_ALam->Write();
V0chi2d2_diff_TF_ALam->Write();
V0chi2d12_diff_TF_ALam->Write();
  
/////// matching hist ////// 
  
hMass_ALam_ALam_TM->Write();
hMass_ALam_ALam_TU->Write();
hMass_ALam_ALam_FM->Write();
hMass_ALam_ALam_FU->Write();
 
hpeak_ALam_ALam_TM->Write();
hpeakd1_ALam_ALam_TM->Write();
hpeakd2_ALam_ALam_TM->Write();
hpeakd12_ALam_ALam_TM->Write();
hpeak_ALam_ALam_TM_bias->Write();
hpeak_ALam_ALam_TM_biasd1->Write();
hpeak_ALam_ALam_TM_biasd2->Write();
hpeak_ALam_ALam_TM_biasd12->Write();

V0chi2d1_diff_TM_ALam->Write();
V0chi2d2_diff_TM_ALam->Write();
V0chi2d12_diff_TM_ALam->Write();

hpeak_ALam_ALam_TU->Write();
hpeakd1_ALam_ALam_TU->Write();
hpeakd2_ALam_ALam_TU->Write();
hpeakd12_ALam_ALam_TU->Write();
hpeak_ALam_ALam_TU_bias->Write();
hpeak_ALam_ALam_TU_biasd1->Write();
hpeak_ALam_ALam_TU_biasd2->Write();
hpeak_ALam_ALam_TU_biasd12->Write();

V0chi2d1_diff_TU_ALam->Write();
V0chi2d2_diff_TU_ALam->Write();
V0chi2d12_diff_TU_ALam->Write();

hpeak_ALam_ALam_FM->Write();
hpeakd1_ALam_ALam_FM->Write();
hpeakd2_ALam_ALam_FM->Write();
hpeakd12_ALam_ALam_FM->Write();
hpeak_ALam_ALam_FM_bias->Write();
hpeak_ALam_ALam_FM_biasd1->Write();
hpeak_ALam_ALam_FM_biasd2->Write();
hpeak_ALam_ALam_FM_biasd12->Write();

V0chi2d1_diff_FM_ALam->Write();
V0chi2d2_diff_FM_ALam->Write();
V0chi2d12_diff_FM_ALam->Write();

hpeak_ALam_ALam_FU->Write();
hpeakd1_ALam_ALam_FU->Write();
hpeakd2_ALam_ALam_FU->Write();
hpeakd12_ALam_ALam_FU->Write();
hpeak_ALam_ALam_FU_bias->Write();
hpeak_ALam_ALam_FU_biasd1->Write();
hpeak_ALam_ALam_FU_biasd2->Write();
hpeak_ALam_ALam_FU_biasd12->Write();

V0chi2d1_diff_FU_ALam->Write();
V0chi2d2_diff_FU_ALam->Write();
V0chi2d12_diff_FU_ALam->Write();


//true matched + true unmatched

hpeak_ALam_ALam_TM_TU->Write();
hpeakd1_ALam_ALam_TM_TU->Write();
hpeakd2_ALam_ALam_TM_TU->Write();
hpeakd12_ALam_ALam_TM_TU->Write();
hpeak_ALam_ALam_TM_TU_bias->Write();
hpeak_ALam_ALam_TM_TU_biasd1->Write();
hpeak_ALam_ALam_TM_TU_biasd2->Write();
hpeak_ALam_ALam_TM_TU_biasd12->Write();

V0chi2d1_diff_TM_TU_ALam->Write();
V0chi2d2_diff_TM_TU_ALam->Write();
V0chi2d12_diff_TM_TU_ALam->Write();

//fake matched + fake unmatched

hpeak_ALam_ALam_FM_FU->Write();
hpeakd1_ALam_ALam_FM_FU->Write();
hpeakd2_ALam_ALam_FM_FU->Write();
hpeakd12_ALam_ALam_FM_FU->Write();
hpeak_ALam_ALam_FM_FU_bias->Write();
hpeak_ALam_ALam_FM_FU_biasd1->Write();
hpeak_ALam_ALam_FM_FU_biasd2->Write();
hpeak_ALam_ALam_FM_FU_biasd12->Write();

V0chi2d1_diff_FM_FU_ALam->Write();
V0chi2d2_diff_FM_FU_ALam->Write();
V0chi2d12_diff_FM_FU_ALam->Write();

hpeak_ALam_ALam_TM_FM->Write();
hpeakd1_ALam_ALam_TM_FM->Write();
hpeakd2_ALam_ALam_TM_FM->Write();
hpeakd12_ALam_ALam_TM_FM->Write();
hpeak_ALam_ALam_TM_FM_bias->Write();
hpeak_ALam_ALam_TM_FM_biasd1->Write();
hpeak_ALam_ALam_TM_FM_biasd2->Write();
hpeak_ALam_ALam_TM_FM_biasd12->Write();

V0chi2d1_diff_TM_FM_ALam->Write();
V0chi2d2_diff_TM_FM_ALam->Write();
V0chi2d12_diff_TM_FM_ALam->Write();

hpeak_ALam_ALam_TM_FU->Write();
hpeakd1_ALam_ALam_TM_FU->Write();
hpeakd2_ALam_ALam_TM_FU->Write();
hpeakd12_ALam_ALam_TM_FU->Write();
hpeak_ALam_ALam_TM_FU_bias->Write();
hpeak_ALam_ALam_TM_FU_biasd1->Write();
hpeak_ALam_ALam_TM_FU_biasd2->Write();
hpeak_ALam_ALam_TM_FU_biasd12->Write();

V0chi2d1_diff_TM_FU_ALam->Write();
V0chi2d2_diff_TM_FU_ALam->Write();
V0chi2d12_diff_TM_FU_ALam->Write();

hpeak_ALam_ALam_TU_FU->Write();
hpeakd1_ALam_ALam_TU_FU->Write();
hpeakd2_ALam_ALam_TU_FU->Write();
hpeakd12_ALam_ALam_TU_FU->Write();
hpeak_ALam_ALam_TU_FU_bias->Write();
hpeak_ALam_ALam_TU_FU_biasd1->Write();
hpeak_ALam_ALam_TU_FU_biasd2->Write();
hpeak_ALam_ALam_TU_FU_biasd12->Write();

V0chi2d1_diff_TU_FU_ALam->Write();
V0chi2d2_diff_TU_FU_ALam->Write();
V0chi2d12_diff_TU_FU_ALam->Write();

hpeak_ALam_ALam_TU_FM->Write();
hpeakd1_ALam_ALam_TU_FM->Write();
hpeakd2_ALam_ALam_TU_FM->Write();
hpeakd12_ALam_ALam_TU_FM->Write();
hpeak_ALam_ALam_TU_FM_bias->Write();
hpeak_ALam_ALam_TU_FM_biasd1->Write();
hpeak_ALam_ALam_TU_FM_biasd2->Write();
hpeak_ALam_ALam_TU_FM_biasd12->Write();

V0chi2d1_diff_TU_FM_ALam->Write();
V0chi2d2_diff_TU_FM_ALam->Write();
V0chi2d12_diff_TU_FM_ALam->Write();
    
} //DONE

void write_LAL_histos(){

//cross
 
//LAL
hMass_Lam_T->Write(); 
hMass_ALam_T->Write(); 
hMass_Lam_ALam_T->Write(); 

hpeak_Lam_ALam_T->Write(); 
hpeakcos_Lam_ALam_T->Write();
hpeak_rot_Lam_ALam_T->Write(); 
hpeak_inv_Lam_ALam_T->Write(); 
hside_Lam_ALam_T->Write(); 
hsideL_Lam_ALam_T->Write(); 
hsideR_Lam_ALam_T->Write(); 
hside_rot_Lam_ALam_T->Write(); 
hside_inv_Lam_ALam_T->Write(); 
hpeakside_Lam_ALam_T->Write(); 
hpeaksideL_Lam_ALam_T->Write(); 
hpeaksideR_Lam_ALam_T->Write(); 
hpeakside_rot_Lam_ALam_T->Write(); 
hpeakside_inv_Lam_ALam_T->Write(); 
hpeakd1_Lam_ALam_T->Write(); 
hpeakd2_Lam_ALam_T->Write(); 
hpeakd12_Lam_ALam_T->Write();
   

hpeak_Lam_ALam_T_mix->Write(); 
hside_Lam_ALam_T_mix->Write(); 
hsideL_Lam_ALam_T_mix->Write(); 
hsideR_Lam_ALam_T_mix->Write(); 
hpeakside_Lam_ALam_T_mix->Write(); 
hpeaksideL_Lam_ALam_T_mix->Write(); 
hpeaksideR_Lam_ALam_T_mix->Write(); 

hpeak_Lam_ALam_T_etamix->Write(); 
hside_Lam_ALam_T_etamix->Write(); 
hsideL_Lam_ALam_T_etamix->Write(); 
hsideR_Lam_ALam_T_etamix->Write(); 
hpeakside_Lam_ALam_T_etamix->Write(); 
hpeaksideL_Lam_ALam_T_etamix->Write(); 
hpeaksideR_Lam_ALam_T_etamix->Write(); 


V0chi2d1_diff_T_Lam_ALam->Write(); 
V0chi2d2_diff_T_Lam_ALam->Write(); 
V0chi2d12_diff_T_Lam_ALam->Write();     


//------------------------------------------------------------------------
//fake
//------------------------------------------------------------------------ 

hMass_Lam_F->Write();
hMass_ALam_F->Write();
hMass_Lam_ALam_F->Write();
  
hpeak_Lam_ALam_F->Write();
hpeakd1_Lam_ALam_F->Write();
hpeakd2_Lam_ALam_F->Write();
hpeakd12_Lam_ALam_F->Write();
hside_Lam_ALam_F->Write();
hpeakside_Lam_ALam_F->Write();
hpeak_Lam_ALam_F_bias->Write();
hpeak_Lam_ALam_F_biasd1->Write();
hpeak_Lam_ALam_F_biasd2->Write();
hpeak_Lam_ALam_F_biasd12->Write();

V0chi2d1_diff_F_Lam_ALam->Write();
V0chi2d2_diff_F_Lam_ALam->Write();
V0chi2d12_diff_F_Lam_ALam->Write();


//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
hpeak_Lam_ALam_TF->Write();
hpeakd1_Lam_ALam_TF->Write();
hpeakd2_Lam_ALam_TF->Write();
hpeakd12_Lam_ALam_TF->Write();
hside_Lam_ALam_TF->Write();
hpeakside_Lam_ALam_TF->Write();
hpeak_Lam_ALam_TF_bias->Write();
hpeak_Lam_ALam_TF_biasd1->Write();
hpeak_Lam_ALam_TF_biasd2->Write();
hpeak_Lam_ALam_TF_biasd12->Write();

V0chi2d1_diff_TF_Lam_ALam->Write();
V0chi2d2_diff_TF_Lam_ALam->Write();
V0chi2d12_diff_TF_Lam_ALam->Write();

  
/////// matching hist ////// 


hMass_Lam_ALam_TM->Write();
hMass_Lam_TM->Write();
hMass_ALam_TM->Write();

hMass_Lam_ALam_TU->Write();
hMass_Lam_TU->Write();
hMass_ALam_TU->Write();

hMass_Lam_ALam_FM->Write();
hMass_Lam_FM->Write();
hMass_ALam_FM->Write();

hMass_Lam_ALam_FU->Write();
hMass_Lam_FU->Write();
hMass_ALam_FU->Write();
 
hpeak_Lam_ALam_TM->Write();
hpeakd1_Lam_ALam_TM->Write();
hpeakd2_Lam_ALam_TM->Write();
hpeakd12_Lam_ALam_TM->Write();
hpeak_Lam_ALam_TM_bias->Write();
hpeak_Lam_ALam_TM_biasd1->Write();
hpeak_Lam_ALam_TM_biasd2->Write();
hpeak_Lam_ALam_TM_biasd12->Write();

V0chi2d1_diff_TM_Lam_ALam->Write();
V0chi2d2_diff_TM_Lam_ALam->Write();
V0chi2d12_diff_TM_Lam_ALam->Write();

hpeak_Lam_ALam_TU->Write();
hpeakd1_Lam_ALam_TU->Write();
hpeakd2_Lam_ALam_TU->Write();
hpeakd12_Lam_ALam_TU->Write();
hpeak_Lam_ALam_TU_bias->Write();
hpeak_Lam_ALam_TU_biasd1->Write();
hpeak_Lam_ALam_TU_biasd2->Write();
hpeak_Lam_ALam_TU_biasd12->Write();

V0chi2d1_diff_TU_Lam_ALam->Write();
V0chi2d2_diff_TU_Lam_ALam->Write();
V0chi2d12_diff_TU_Lam_ALam->Write();

hpeak_Lam_ALam_FM->Write();
hpeakd1_Lam_ALam_FM->Write();
hpeakd2_Lam_ALam_FM->Write();
hpeakd12_Lam_ALam_FM->Write();
hpeak_Lam_ALam_FM_bias->Write();
hpeak_Lam_ALam_FM_biasd1->Write();
hpeak_Lam_ALam_FM_biasd2->Write();
hpeak_Lam_ALam_FM_biasd12->Write();

V0chi2d1_diff_FM_Lam_ALam->Write();
V0chi2d2_diff_FM_Lam_ALam->Write();
V0chi2d12_diff_FM_Lam_ALam->Write();

hpeak_Lam_ALam_FU->Write();
hpeakd1_Lam_ALam_FU->Write();
hpeakd2_Lam_ALam_FU->Write();
hpeakd12_Lam_ALam_FU->Write();
hpeak_Lam_ALam_FU_bias->Write();
hpeak_Lam_ALam_FU_biasd1->Write();
hpeak_Lam_ALam_FU_biasd2->Write();
hpeak_Lam_ALam_FU_biasd12->Write();

V0chi2d1_diff_FU_Lam_ALam->Write();
V0chi2d2_diff_FU_Lam_ALam->Write();
V0chi2d12_diff_FU_Lam_ALam->Write();

//true matched + true unmatched

hpeak_Lam_ALam_TM_TU->Write();
hpeakd1_Lam_ALam_TM_TU->Write();
hpeakd2_Lam_ALam_TM_TU->Write();
hpeakd12_Lam_ALam_TM_TU->Write();
hpeak_Lam_ALam_TM_TU_bias->Write();
hpeak_Lam_ALam_TM_TU_biasd1->Write();
hpeak_Lam_ALam_TM_TU_biasd2->Write();
hpeak_Lam_ALam_TM_TU_biasd12->Write();

V0chi2d1_diff_TM_TU_Lam_ALam->Write();
V0chi2d2_diff_TM_TU_Lam_ALam->Write();
V0chi2d12_diff_TM_TU_Lam_ALam->Write();

//fake matched + fake unmatched

hpeak_Lam_ALam_FM_FU->Write();
hpeakd1_Lam_ALam_FM_FU->Write();
hpeakd2_Lam_ALam_FM_FU->Write();
hpeakd12_Lam_ALam_FM_FU->Write();
hpeak_Lam_ALam_FM_FU_bias->Write();
hpeak_Lam_ALam_FM_FU_biasd1->Write();
hpeak_Lam_ALam_FM_FU_biasd2->Write();
hpeak_Lam_ALam_FM_FU_biasd12->Write();

V0chi2d1_diff_FM_FU_Lam_ALam->Write();
V0chi2d2_diff_FM_FU_Lam_ALam->Write();
V0chi2d12_diff_FM_FU_Lam_ALam->Write();

//true matched + fake matched

hpeak_Lam_ALam_TM_FM->Write();
hpeakd1_Lam_ALam_TM_FM->Write();
hpeakd2_Lam_ALam_TM_FM->Write();
hpeakd12_Lam_ALam_TM_FM->Write();
hpeak_Lam_ALam_TM_FM_bias->Write();
hpeak_Lam_ALam_TM_FM_biasd1->Write();
hpeak_Lam_ALam_TM_FM_biasd2->Write();
hpeak_Lam_ALam_TM_FM_biasd12->Write();

V0chi2d1_diff_TM_FM_Lam_ALam->Write();
V0chi2d2_diff_TM_FM_Lam_ALam->Write();
V0chi2d12_diff_TM_FM_Lam_ALam->Write();

//true matched + fake unmatched

hpeak_Lam_ALam_TM_FU->Write();
hpeakd1_Lam_ALam_TM_FU->Write();
hpeakd2_Lam_ALam_TM_FU->Write();
hpeakd12_Lam_ALam_TM_FU->Write();
hpeak_Lam_ALam_TM_FU_bias->Write();
hpeak_Lam_ALam_TM_FU_biasd1->Write();
hpeak_Lam_ALam_TM_FU_biasd2->Write();
hpeak_Lam_ALam_TM_FU_biasd12->Write();

V0chi2d1_diff_TM_FU_Lam_ALam->Write();
V0chi2d2_diff_TM_FU_Lam_ALam->Write();
V0chi2d12_diff_TM_FU_Lam_ALam->Write();

//true unmatched + fake unmatched

hpeak_Lam_ALam_TU_FU->Write();
hpeakd1_Lam_ALam_TU_FU->Write();
hpeakd2_Lam_ALam_TU_FU->Write();
hpeakd12_Lam_ALam_TU_FU->Write();
hpeak_Lam_ALam_TU_FU_bias->Write();
hpeak_Lam_ALam_TU_FU_biasd1->Write();
hpeak_Lam_ALam_TU_FU_biasd2->Write();
hpeak_Lam_ALam_TU_FU_biasd12->Write();

V0chi2d1_diff_TU_FU_Lam_ALam->Write();
V0chi2d2_diff_TU_FU_Lam_ALam->Write();
V0chi2d12_diff_TU_FU_Lam_ALam->Write();

//true unmatched + fake matched

hpeak_Lam_ALam_TU_FM->Write();
hpeakd1_Lam_ALam_TU_FM->Write();
hpeakd2_Lam_ALam_TU_FM->Write();
hpeakd12_Lam_ALam_TU_FM->Write();
hpeak_Lam_ALam_TU_FM_bias->Write();
hpeak_Lam_ALam_TU_FM_biasd1->Write();
hpeak_Lam_ALam_TU_FM_biasd2->Write();
hpeak_Lam_ALam_TU_FM_biasd12->Write();

V0chi2d1_diff_TU_FM_Lam_ALam->Write();
V0chi2d2_diff_TU_FM_Lam_ALam->Write();
V0chi2d12_diff_TU_FM_Lam_ALam->Write();
    
  

} //DONE

void write_KL_histos(){
hMass_K0sL_T->Write(); 
hMass_LamK_T->Write(); 
hMass_K0s_Lam_T->Write();
hpeak_K0s_Lam_T->Write();
hpeakcos_K0s_Lam_T->Write();
hpeak_rot_K0s_Lam_T->Write();
hpeak_inv_K0s_Lam_T->Write();
hside_K0s_Lam_T->Write();
hsideL_K0s_Lam_T->Write();
hsideR_K0s_Lam_T->Write();
hside_rot_K0s_Lam_T->Write();
hside_inv_K0s_Lam_T->Write();
hpeakside_K0s_Lam_T->Write();
hpeaksideL_K0s_Lam_T->Write();
hpeaksideR_K0s_Lam_T->Write();
hpeakside_rot_K0s_Lam_T->Write();
hpeakside_inv_K0s_Lam_T->Write();
hpeakd1_K0s_Lam_T->Write();
hpeakd2_K0s_Lam_T->Write();
hpeakd12_K0s_Lam_T->Write();
hpeak_K0s_Lam_T_mix->Write();
hside_K0s_Lam_T_mix->Write();
hsideL_K0s_Lam_T_mix->Write();
hsideR_K0s_Lam_T_mix->Write();
hpeakside_K0s_Lam_T_mix->Write();
hpeaksideL_K0s_Lam_T_mix->Write();
hpeaksideR_K0s_Lam_T_mix->Write();
hpeak_K0s_Lam_T_etamix->Write();
hside_K0s_Lam_T_etamix->Write();
hsideL_K0s_Lam_T_etamix->Write();
hsideR_K0s_Lam_T_etamix->Write();
hpeakside_K0s_Lam_T_etamix->Write();
hpeaksideL_K0s_Lam_T_etamix->Write();
hpeaksideR_K0s_Lam_T_etamix->Write();
V0chi2d1_diff_T_K0s_Lam->Write();
V0chi2d2_diff_T_K0s_Lam->Write();
V0chi2d12_diff_T_K0s_Lam->Write();
  
//------------------------------------------------------------------------
//fake
//------------------------------------------------------------------------ 

hMass_K0sL_F->Write();   
hMass_LamK_F->Write();
hMass_K0s_Lam_F->Write();
  
hpeak_K0s_Lam_F->Write();
hpeakd1_K0s_Lam_F->Write();
hpeakd2_K0s_Lam_F->Write();
hpeakd12_K0s_Lam_F->Write();
hside_K0s_Lam_F->Write();
hpeakside_K0s_Lam_F->Write();
hpeak_K0s_Lam_F_bias->Write();
hpeak_K0s_Lam_F_biasd1->Write();
hpeak_K0s_Lam_F_biasd2->Write();
hpeak_K0s_Lam_F_biasd12->Write();

V0chi2d1_diff_F_K0s_Lam->Write();
V0chi2d2_diff_F_K0s_Lam->Write();
V0chi2d12_diff_F_K0s_Lam->Write();
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
hpeak_K0s_Lam_TF->Write();
hpeakd1_K0s_Lam_TF->Write();
hpeakd2_K0s_Lam_TF->Write();
hpeakd12_K0s_Lam_TF->Write();
hside_K0s_Lam_TF->Write();
hpeakside_K0s_Lam_TF->Write(); 
hpeak_K0s_Lam_TF_bias->Write();
hpeak_K0s_Lam_TF_biasd1->Write();
hpeak_K0s_Lam_TF_biasd2->Write();
hpeak_K0s_Lam_TF_biasd12->Write();

V0chi2d1_diff_TF_K0s_Lam->Write();
V0chi2d2_diff_TF_K0s_Lam->Write();
V0chi2d12_diff_TF_K0s_Lam->Write();

  
/////// matching hist ////// 

hMass_K0sL_TM->Write();  
hMass_LamK_TM->Write();
hMass_K0s_Lam_TM->Write();

hMass_K0sL_TU->Write();
hMass_LamK_TU->Write();
hMass_K0s_Lam_TU->Write();

hMass_K0sL_FM->Write();
hMass_LamK_FM->Write();
hMass_K0s_Lam_FM->Write();

hMass_K0sL_FU->Write();
hMass_LamK_FU->Write();
hMass_K0s_Lam_FU->Write();


hpeak_K0s_Lam_TM->Write();
hpeakd1_K0s_Lam_TM->Write();
hpeakd2_K0s_Lam_TM->Write();
hpeakd12_K0s_Lam_TM->Write();
hpeak_K0s_Lam_TM_bias->Write();
hpeak_K0s_Lam_TM_biasd1->Write();
hpeak_K0s_Lam_TM_biasd2->Write();
hpeak_K0s_Lam_TM_biasd12->Write();

V0chi2d1_diff_TM_K0s_Lam->Write();
V0chi2d2_diff_TM_K0s_Lam->Write();
V0chi2d12_diff_TM_K0s_Lam->Write();

hpeak_K0s_Lam_TU->Write();
hpeakd1_K0s_Lam_TU->Write();
hpeakd2_K0s_Lam_TU->Write();
hpeakd12_K0s_Lam_TU->Write();
hpeak_K0s_Lam_TU_bias->Write();
hpeak_K0s_Lam_TU_biasd1->Write();
hpeak_K0s_Lam_TU_biasd2->Write();
hpeak_K0s_Lam_TU_biasd12->Write();

V0chi2d1_diff_TU_K0s_Lam->Write();
V0chi2d2_diff_TU_K0s_Lam->Write();
V0chi2d12_diff_TU_K0s_Lam->Write();

hpeak_K0s_Lam_FM->Write();
hpeakd1_K0s_Lam_FM->Write();
hpeakd2_K0s_Lam_FM->Write();
hpeakd12_K0s_Lam_FM->Write();
hpeak_K0s_Lam_FM_bias->Write();
hpeak_K0s_Lam_FM_biasd1->Write();
hpeak_K0s_Lam_FM_biasd2->Write();
hpeak_K0s_Lam_FM_biasd12->Write();

V0chi2d1_diff_FM_K0s_Lam->Write();
V0chi2d2_diff_FM_K0s_Lam->Write();
V0chi2d12_diff_FM_K0s_Lam->Write();

hpeak_K0s_Lam_FU->Write();
hpeakd1_K0s_Lam_FU->Write();
hpeakd2_K0s_Lam_FU->Write();
hpeakd12_K0s_Lam_FU->Write();
hpeak_K0s_Lam_FU_bias->Write();
hpeak_K0s_Lam_FU_biasd1->Write();
hpeak_K0s_Lam_FU_biasd2->Write();
hpeak_K0s_Lam_FU_biasd12->Write();

V0chi2d1_diff_FU_K0s_Lam->Write();
V0chi2d2_diff_FU_K0s_Lam->Write();
V0chi2d12_diff_FU_K0s_Lam->Write();

//true matched + true unmatched

hpeak_K0s_Lam_TM_TU->Write();
hpeakd1_K0s_Lam_TM_TU->Write();
hpeakd2_K0s_Lam_TM_TU->Write();
hpeakd12_K0s_Lam_TM_TU->Write();
hpeak_K0s_Lam_TM_TU_bias->Write();
hpeak_K0s_Lam_TM_TU_biasd1->Write();
hpeak_K0s_Lam_TM_TU_biasd2->Write();
hpeak_K0s_Lam_TM_TU_biasd12->Write();

V0chi2d1_diff_TM_TU_K0s_Lam->Write();
V0chi2d2_diff_TM_TU_K0s_Lam->Write();
V0chi2d12_diff_TM_TU_K0s_Lam->Write();

//fake matched + fake unmatched

hpeak_K0s_Lam_FM_FU->Write();
hpeakd1_K0s_Lam_FM_FU->Write();
hpeakd2_K0s_Lam_FM_FU->Write();
hpeakd12_K0s_Lam_FM_FU->Write();
hpeak_K0s_Lam_FM_FU_bias->Write();
hpeak_K0s_Lam_FM_FU_biasd1->Write();
hpeak_K0s_Lam_FM_FU_biasd2->Write();
hpeak_K0s_Lam_FM_FU_biasd12->Write();

V0chi2d1_diff_FM_FU_K0s_Lam->Write();
V0chi2d2_diff_FM_FU_K0s_Lam->Write();
V0chi2d12_diff_FM_FU_K0s_Lam->Write();

//true matched + fake matched

hpeak_K0s_Lam_TM_FM->Write();
hpeakd1_K0s_Lam_TM_FM->Write();
hpeakd2_K0s_Lam_TM_FM->Write();
hpeakd12_K0s_Lam_TM_FM->Write();
hpeak_K0s_Lam_TM_FM_bias->Write();
hpeak_K0s_Lam_TM_FM_biasd1->Write();
hpeak_K0s_Lam_TM_FM_biasd2->Write();
hpeak_K0s_Lam_TM_FM_biasd12->Write();

V0chi2d1_diff_TM_FM_K0s_Lam->Write();
V0chi2d2_diff_TM_FM_K0s_Lam->Write();
V0chi2d12_diff_TM_FM_K0s_Lam->Write();

//true matched + fake unmatched

hpeak_K0s_Lam_TM_FU->Write();
hpeakd1_K0s_Lam_TM_FU->Write();
hpeakd2_K0s_Lam_TM_FU->Write();
hpeakd12_K0s_Lam_TM_FU->Write();
hpeak_K0s_Lam_TM_FU_bias->Write();
hpeak_K0s_Lam_TM_FU_biasd1->Write();
hpeak_K0s_Lam_TM_FU_biasd2->Write();
hpeak_K0s_Lam_TM_FU_biasd12->Write();

V0chi2d1_diff_TM_FU_K0s_Lam->Write();
V0chi2d2_diff_TM_FU_K0s_Lam->Write();
V0chi2d12_diff_TM_FU_K0s_Lam->Write();

//true unmatched + fake unmatched

hpeak_K0s_Lam_TU_FU->Write();
hpeakd1_K0s_Lam_TU_FU->Write();
hpeakd2_K0s_Lam_TU_FU->Write();
hpeakd12_K0s_Lam_TU_FU->Write();
hpeak_K0s_Lam_TU_FU_bias->Write();
hpeak_K0s_Lam_TU_FU_biasd1->Write();
hpeak_K0s_Lam_TU_FU_biasd2->Write();
hpeak_K0s_Lam_TU_FU_biasd12->Write();

V0chi2d1_diff_TU_FU_K0s_Lam->Write();
V0chi2d2_diff_TU_FU_K0s_Lam->Write();
V0chi2d12_diff_TU_FU_K0s_Lam->Write();

//true unmatched + fake matched

hpeak_K0s_Lam_TU_FM->Write();
hpeakd1_K0s_Lam_TU_FM->Write();
hpeakd2_K0s_Lam_TU_FM->Write();
hpeakd12_K0s_Lam_TU_FM->Write();
hpeak_K0s_Lam_TU_FM_bias->Write();
hpeak_K0s_Lam_TU_FM_biasd1->Write();
hpeak_K0s_Lam_TU_FM_biasd2->Write();
hpeak_K0s_Lam_TU_FM_biasd12->Write();

V0chi2d1_diff_TU_FM_K0s_Lam->Write();
V0chi2d2_diff_TU_FM_K0s_Lam->Write();
V0chi2d12_diff_TU_FM_K0s_Lam->Write();

} //DONE

void write_KAL_histos(){

hMass_K0sAL_T->Write(); 
hMass_ALamK_T->Write(); 
hMass_K0s_ALam_T->Write();
hpeak_K0s_ALam_T->Write();
hpeakcos_K0s_ALam_T->Write();
hpeak_rot_K0s_ALam_T->Write();
hpeak_inv_K0s_ALam_T->Write();
hside_K0s_ALam_T->Write();
hsideL_K0s_ALam_T->Write();
hsideR_K0s_ALam_T->Write();
hside_rot_K0s_ALam_T->Write();
hside_inv_K0s_ALam_T->Write();
hpeakside_K0s_ALam_T->Write();
hpeaksideL_K0s_ALam_T->Write();
hpeaksideR_K0s_ALam_T->Write();
hpeakside_rot_K0s_ALam_T->Write();
hpeakside_inv_K0s_ALam_T->Write();
hpeakd1_K0s_ALam_T->Write();
hpeakd2_K0s_ALam_T->Write();
hpeakd12_K0s_ALam_T->Write();
hpeak_K0s_ALam_T_mix->Write();
hside_K0s_ALam_T_mix->Write();
hsideL_K0s_ALam_T_mix->Write();
hsideR_K0s_ALam_T_mix->Write();
hpeakside_K0s_ALam_T_mix->Write();
hpeaksideL_K0s_ALam_T_mix->Write();
hpeaksideR_K0s_ALam_T_mix->Write();
hpeak_K0s_ALam_T_etamix->Write();
hside_K0s_ALam_T_etamix->Write();
hsideL_K0s_ALam_T_etamix->Write();
hsideR_K0s_ALam_T_etamix->Write();
hpeakside_K0s_ALam_T_etamix->Write();
hpeaksideL_K0s_ALam_T_etamix->Write();
hpeaksideR_K0s_ALam_T_etamix->Write();
V0chi2d1_diff_T_K0s_ALam->Write();
V0chi2d2_diff_T_K0s_ALam->Write();
V0chi2d12_diff_T_K0s_ALam->Write();
  
//------------------------------------------------------------------------
//fake
//------------------------------------------------------------------------ 

hMass_K0sAL_F->Write();   
hMass_ALamK_F->Write();
hMass_K0s_ALam_F->Write();
hpeak_K0s_ALam_F->Write();
hpeakd1_K0s_ALam_F->Write();
hpeakd2_K0s_ALam_F->Write();
hpeakd12_K0s_ALam_F->Write();
hside_K0s_ALam_F->Write();
hpeakside_K0s_ALam_F->Write();
hpeak_K0s_ALam_F_bias->Write();
hpeak_K0s_ALam_F_biasd1->Write();
hpeak_K0s_ALam_F_biasd2->Write();
hpeak_K0s_ALam_F_biasd12->Write();
V0chi2d1_diff_F_K0s_ALam->Write();
V0chi2d2_diff_F_K0s_ALam->Write();
V0chi2d12_diff_F_K0s_ALam->Write();
  
//------------------------------------------------------------------------
//true + fake
//------------------------------------------------------------------------ 
  
hpeak_K0s_ALam_TF->Write();
hpeakd1_K0s_ALam_TF->Write();
hpeakd2_K0s_ALam_TF->Write();
hpeakd12_K0s_ALam_TF->Write();
hside_K0s_ALam_TF->Write();
hpeakside_K0s_ALam_TF->Write(); 
hpeak_K0s_ALam_TF_bias->Write();
hpeak_K0s_ALam_TF_biasd1->Write();
hpeak_K0s_ALam_TF_biasd2->Write();
hpeak_K0s_ALam_TF_biasd12->Write();

V0chi2d1_diff_TF_K0s_ALam->Write();
V0chi2d2_diff_TF_K0s_ALam->Write();
V0chi2d12_diff_TF_K0s_ALam->Write();

  
/////// matching hist ////// 

hMass_K0sAL_TM->Write();  
hMass_ALamK_TM->Write();
hMass_K0s_ALam_TM->Write();

hMass_K0sAL_TU->Write();
hMass_ALamK_TU->Write();
hMass_K0s_ALam_TU->Write();

hMass_K0sAL_FM->Write();
hMass_ALamK_FM->Write();
hMass_K0s_ALam_FM->Write();

hMass_K0sAL_FU->Write();
hMass_ALamK_FU->Write();
hMass_K0s_ALam_FU->Write();


hpeak_K0s_ALam_TM->Write();
hpeakd1_K0s_ALam_TM->Write();
hpeakd2_K0s_ALam_TM->Write();
hpeakd12_K0s_ALam_TM->Write();
hpeak_K0s_ALam_TM_bias->Write();
hpeak_K0s_ALam_TM_biasd1->Write();
hpeak_K0s_ALam_TM_biasd2->Write();
hpeak_K0s_ALam_TM_biasd12->Write();

V0chi2d1_diff_TM_K0s_ALam->Write();
V0chi2d2_diff_TM_K0s_ALam->Write();
V0chi2d12_diff_TM_K0s_ALam->Write();

hpeak_K0s_ALam_TU->Write();
hpeakd1_K0s_ALam_TU->Write();
hpeakd2_K0s_ALam_TU->Write();
hpeakd12_K0s_ALam_TU->Write();
hpeak_K0s_ALam_TU_bias->Write();
hpeak_K0s_ALam_TU_biasd1->Write();
hpeak_K0s_ALam_TU_biasd2->Write();
hpeak_K0s_ALam_TU_biasd12->Write();

V0chi2d1_diff_TU_K0s_ALam->Write();
V0chi2d2_diff_TU_K0s_ALam->Write();
V0chi2d12_diff_TU_K0s_ALam->Write();

hpeak_K0s_ALam_FM->Write();
hpeakd1_K0s_ALam_FM->Write();
hpeakd2_K0s_ALam_FM->Write();
hpeakd12_K0s_ALam_FM->Write();
hpeak_K0s_ALam_FM_bias->Write();
hpeak_K0s_ALam_FM_biasd1->Write();
hpeak_K0s_ALam_FM_biasd2->Write();
hpeak_K0s_ALam_FM_biasd12->Write();

V0chi2d1_diff_FM_K0s_ALam->Write();
V0chi2d2_diff_FM_K0s_ALam->Write();
V0chi2d12_diff_FM_K0s_ALam->Write();

hpeak_K0s_ALam_FU->Write();
hpeakd1_K0s_ALam_FU->Write();
hpeakd2_K0s_ALam_FU->Write();
hpeakd12_K0s_ALam_FU->Write();
hpeak_K0s_ALam_FU_bias->Write();
hpeak_K0s_ALam_FU_biasd1->Write();
hpeak_K0s_ALam_FU_biasd2->Write();
hpeak_K0s_ALam_FU_biasd12->Write();

V0chi2d1_diff_FU_K0s_ALam->Write();
V0chi2d2_diff_FU_K0s_ALam->Write();
V0chi2d12_diff_FU_K0s_ALam->Write();

//true matched + true unmatched

hpeak_K0s_ALam_TM_TU->Write();
hpeakd1_K0s_ALam_TM_TU->Write();
hpeakd2_K0s_ALam_TM_TU->Write();
hpeakd12_K0s_ALam_TM_TU->Write();
hpeak_K0s_ALam_TM_TU_bias->Write();
hpeak_K0s_ALam_TM_TU_biasd1->Write();
hpeak_K0s_ALam_TM_TU_biasd2->Write();
hpeak_K0s_ALam_TM_TU_biasd12->Write();

V0chi2d1_diff_TM_TU_K0s_ALam->Write();
V0chi2d2_diff_TM_TU_K0s_ALam->Write();
V0chi2d12_diff_TM_TU_K0s_ALam->Write();

//fake matched + fake unmatched

hpeak_K0s_ALam_FM_FU->Write();
hpeakd1_K0s_ALam_FM_FU->Write();
hpeakd2_K0s_ALam_FM_FU->Write();
hpeakd12_K0s_ALam_FM_FU->Write();
hpeak_K0s_ALam_FM_FU_bias->Write();
hpeak_K0s_ALam_FM_FU_biasd1->Write();
hpeak_K0s_ALam_FM_FU_biasd2->Write();
hpeak_K0s_ALam_FM_FU_biasd12->Write();

V0chi2d1_diff_FM_FU_K0s_ALam->Write();
V0chi2d2_diff_FM_FU_K0s_ALam->Write();
V0chi2d12_diff_FM_FU_K0s_ALam->Write();

//true matched + fake matched

hpeak_K0s_ALam_TM_FM->Write();
hpeakd1_K0s_ALam_TM_FM->Write();
hpeakd2_K0s_ALam_TM_FM->Write();
hpeakd12_K0s_ALam_TM_FM->Write();
hpeak_K0s_ALam_TM_FM_bias->Write();
hpeak_K0s_ALam_TM_FM_biasd1->Write();
hpeak_K0s_ALam_TM_FM_biasd2->Write();
hpeak_K0s_ALam_TM_FM_biasd12->Write();

V0chi2d1_diff_TM_FM_K0s_ALam->Write();
V0chi2d2_diff_TM_FM_K0s_ALam->Write();
V0chi2d12_diff_TM_FM_K0s_ALam->Write();

//true matched + fake unmatched

hpeak_K0s_ALam_TM_FU->Write();
hpeakd1_K0s_ALam_TM_FU->Write();
hpeakd2_K0s_ALam_TM_FU->Write();
hpeakd12_K0s_ALam_TM_FU->Write();
hpeak_K0s_ALam_TM_FU_bias->Write();
hpeak_K0s_ALam_TM_FU_biasd1->Write();
hpeak_K0s_ALam_TM_FU_biasd2->Write();
hpeak_K0s_ALam_TM_FU_biasd12->Write();

V0chi2d1_diff_TM_FU_K0s_ALam->Write();
V0chi2d2_diff_TM_FU_K0s_ALam->Write();
V0chi2d12_diff_TM_FU_K0s_ALam->Write();

//true unmatched + fake unmatched

hpeak_K0s_ALam_TU_FU->Write();
hpeakd1_K0s_ALam_TU_FU->Write();
hpeakd2_K0s_ALam_TU_FU->Write();
hpeakd12_K0s_ALam_TU_FU->Write();
hpeak_K0s_ALam_TU_FU_bias->Write();
hpeak_K0s_ALam_TU_FU_biasd1->Write();
hpeak_K0s_ALam_TU_FU_biasd2->Write();
hpeak_K0s_ALam_TU_FU_biasd12->Write();

V0chi2d1_diff_TU_FU_K0s_ALam->Write();
V0chi2d2_diff_TU_FU_K0s_ALam->Write();
V0chi2d12_diff_TU_FU_K0s_ALam->Write();

//true unmatched + fake matched

hpeak_K0s_ALam_TU_FM->Write();
hpeakd1_K0s_ALam_TU_FM->Write();
hpeakd2_K0s_ALam_TU_FM->Write();
hpeakd12_K0s_ALam_TU_FM->Write();
hpeak_K0s_ALam_TU_FM_bias->Write();
hpeak_K0s_ALam_TU_FM_biasd1->Write();
hpeak_K0s_ALam_TU_FM_biasd2->Write();
hpeak_K0s_ALam_TU_FM_biasd12->Write();

V0chi2d1_diff_TU_FM_K0s_ALam->Write();
V0chi2d2_diff_TU_FM_K0s_ALam->Write();
V0chi2d12_diff_TU_FM_K0s_ALam->Write();

} //DONE

void write_shapes(){

hpeak_K0s_K0s_TM_mix->Write(); 
hpeak_K0s_K0s_TU_mix->Write(); 
hpeak_K0s_K0s_TM_TU_mix->Write(); 

hpeak_Lam_Lam_TM_mix->Write(); 
hpeak_Lam_Lam_TU_mix->Write(); 
hpeak_Lam_Lam_TM_TU_mix->Write(); 

hpeak_ALam_ALam_TM_mix->Write(); 
hpeak_ALam_ALam_TU_mix->Write(); 
hpeak_ALam_ALam_TM_TU_mix->Write(); 

hpeak_Lam_ALam_TM_mix->Write(); 
hpeak_Lam_ALam_TU_mix->Write(); 
hpeak_Lam_ALam_TM_TU_mix->Write(); 

hpeak_K0s_Lam_TM_mix->Write(); 
hpeak_K0s_Lam_TU_mix->Write(); 
hpeak_K0s_Lam_TM_TU_mix->Write(); 

hpeak_K0s_ALam_TM_mix->Write(); 
hpeak_K0s_ALam_TU_mix->Write(); 
hpeak_K0s_ALam_TM_TU_mix->Write(); 


} //DONE

void write_feeddown(){

V0_chi2_SS_LLFD_XI->Write(); 
V0_chi2_OS_LLFD_XI->Write(); 
V0_chi2_SS_ALFD_XI->Write(); 
V0_chi2_OS_ALFD_XI->Write(); 
V0_chi2_SS_KLFD_XI->Write(); 
V0_chi2_OS_KLFD_XI->Write(); 

V0_chi2_SS_LLFD_AXI->Write(); 
V0_chi2_OS_LLFD_AXI->Write(); 
V0_chi2_SS_ALFD_AXI->Write(); 
V0_chi2_OS_ALFD_AXI->Write(); 
V0_chi2_SS_KLFD_AXI->Write(); 
V0_chi2_OS_KLFD_AXI->Write(); 

V0_chi2_SS_LLFD_OM->Write(); 
V0_chi2_OS_LLFD_OM->Write(); 
V0_chi2_SS_ALFD_OM->Write(); 
V0_chi2_OS_ALFD_OM->Write(); 
V0_chi2_SS_KLFD_OM->Write(); 
V0_chi2_OS_KLFD_OM->Write(); 

V0_chi2_SS_LLFD_AOM->Write(); 
V0_chi2_OS_LLFD_AOM->Write(); 
V0_chi2_SS_ALFD_AOM->Write(); 
V0_chi2_OS_ALFD_AOM->Write(); 
V0_chi2_SS_KLFD_AOM->Write(); 
V0_chi2_OS_KLFD_AOM->Write(); 



hpeak_Lam_LamFD->Write(); 
hpeak_LamFD_LamFD->Write(); 
hpeak_Lam_LamFD_mix->Write(); 
hpeak_LamFD_LamFD_mix->Write(); 

hpeak_ALam_ALamFD->Write(); 
hpeak_ALamFD_ALamFD->Write(); 
hpeak_ALam_ALamFD_mix->Write(); 
hpeak_ALamFD_ALamFD_mix->Write(); 

hpeak_LALFD->Write(); 
hpeak_LALFD_mix->Write(); 
hpeak_LFDALFD->Write(); 
hpeak_LFDALFD_mix->Write(); 


hpeak_KLFD->Write(); 
hpeak_KLFD_mix->Write(); 

hpeak_KALFD->Write(); 
hpeak_KALFD_mix->Write(); 


} //DONE

void write_hdibaryon(){

h_HDibaryon_Lam_Lam->Write(); 
h_HDibaryon_rot_Lam_Lam->Write(); 
h_HDibaryon_inv_Lam_Lam->Write(); 
h_HDibaryon_Lam_Lam_mix->Write(); 

h_HDibaryon_Lam_Lam_side->Write(); 
h_HDibaryon_rot_Lam_Lam_side->Write(); 
h_HDibaryon_inv_Lam_Lam_side->Write(); 
h_HDibaryon_Lam_Lam_mix_side->Write(); 

h_HDibaryon_Lam_Lam_peakside->Write(); 
h_HDibaryon_rot_Lam_Lam_peakside->Write(); 
h_HDibaryon_inv_Lam_Lam_peakside->Write(); 
h_HDibaryon_Lam_Lam_mix_peakside->Write(); 


h_HDibaryon_ALam_ALam->Write(); 
h_HDibaryon_rot_ALam_ALam->Write(); 
h_HDibaryon_inv_ALam_ALam->Write(); 
h_HDibaryon_ALam_ALam_mix->Write(); 

h_HDibaryon_ALam_ALam_peakside->Write(); 
h_HDibaryon_rot_ALam_ALam_peakside->Write(); 
h_HDibaryon_inv_ALam_ALam_peakside->Write(); 
h_HDibaryon_ALam_ALam_mix_peakside->Write(); 

h_HDibaryon_Lam_ALam->Write(); 
h_HDibaryon_Lam_ALam_side->Write(); 
h_HDibaryon_Lam_ALam_peakside->Write(); 


}

void write_f2(){

h_Hf2_K0s_K0s->Write(); 
h_Hf2_rot_K0s_K0s->Write(); 
h_Hf2_inv_K0s_K0s->Write(); 
h_Hf2_K0s_K0s_mix->Write(); 

h_Hf2_K0s_K0s_side->Write(); 
h_Hf2_rot_K0s_K0s_side->Write(); 
h_Hf2_inv_K0s_K0s_side->Write(); 
h_Hf2_K0s_K0s_mix_side->Write(); 

h_Hf2_K0s_K0s_peakside->Write(); 
h_Hf2_rot_K0s_K0s_peakside->Write(); 
h_Hf2_inv_K0s_K0s_peakside->Write(); 
h_Hf2_K0s_K0s_mix_peakside->Write(); 


}

void write_cas1820(){

h_Hcasc1820_K0s_Lam->Write(); 
h_Hcasc1820_rot_K0s_Lam->Write(); 
h_Hcasc1820_inv_K0s_Lam->Write(); 
h_Hcasc1820_K0s_Lam_mix->Write(); 

h_Hcasc1820_K0s_Lam_side->Write(); 
h_Hcasc1820_rot_K0s_Lam_side->Write(); 
h_Hcasc1820_inv_K0s_Lam_side->Write(); 
h_Hcasc1820_K0s_Lam_mix_side->Write(); 

h_Hcasc1820_K0s_Lam_peakside->Write(); 
h_Hcasc1820_rot_K0s_Lam_peakside->Write(); 
h_Hcasc1820_inv_K0s_Lam_peakside->Write(); 
h_Hcasc1820_K0s_Lam_mix_peakside->Write(); 

h_Hcasc1820_K0s_ALam->Write(); 
h_Hcasc1820_rot_K0s_ALam->Write(); 
h_Hcasc1820_inv_K0s_ALam->Write(); 
h_Hcasc1820_K0s_ALam_mix->Write(); 

h_Hcasc1820_K0s_ALam_side->Write(); 
h_Hcasc1820_rot_K0s_ALam_side->Write(); 
h_Hcasc1820_inv_K0s_ALam_side->Write(); 
h_Hcasc1820_K0s_ALam_mix_side->Write(); 

h_Hcasc1820_K0s_ALam_peakside->Write(); 
h_Hcasc1820_rot_K0s_ALam_peakside->Write(); 
h_Hcasc1820_inv_K0s_ALam_peakside->Write(); 
h_Hcasc1820_K0s_ALam_mix_peakside->Write(); 


}

void write_jets(){

//jets

h_jet_jet->Write(); 
h_inv_jet_jet->Write(); 
h_rot_jet_jet->Write(); 
hMix_jet_jet->Write(); 

//K0sK0s

h_K0s_K0s_samejet->Write(); 
h_inv_K0s_K0s_samejet->Write(); 
h_rot_K0s_K0s_samejet->Write(); 
hMix_K0s_K0s_samejet->Write(); 

h_K0s_K0s_diffjet->Write(); 
h_inv_K0s_K0s_diffjet->Write(); 
h_rot_K0s_K0s_diffjet->Write(); 
hMix_K0s_K0s_diffjet->Write(); 
  
h_K0s_K0s_nojet->Write(); 
h_inv_K0s_K0s_nojet->Write(); 
h_rot_K0s_K0s_nojet->Write(); 
hMix_K0s_K0s_nojet->Write(); 

h_K0s_K0s_only1jet->Write(); 
h_inv_K0s_K0s_only1jet->Write(); 
h_rot_K0s_K0s_only1jet->Write(); 
hMix_K0s_K0s_only1jet->Write(); 

//LamLam

h_Lam_Lam_samejet->Write(); 
h_inv_Lam_Lam_samejet->Write(); 
h_rot_Lam_Lam_samejet->Write(); 
hMix_Lam_Lam_samejet->Write(); 

h_Lam_Lam_diffjet->Write(); 
h_inv_Lam_Lam_diffjet->Write(); 
h_rot_Lam_Lam_diffjet->Write(); 
hMix_Lam_Lam_diffjet->Write(); 
  
h_Lam_Lam_nojet->Write(); 
h_inv_Lam_Lam_nojet->Write(); 
h_rot_Lam_Lam_nojet->Write(); 
hMix_Lam_Lam_nojet->Write(); 

h_Lam_Lam_only1jet->Write();
h_inv_Lam_Lam_only1jet->Write();
h_rot_Lam_Lam_only1jet->Write();
hMix_Lam_Lam_only1jet->Write();

//ALamALam

h_ALam_ALam_samejet->Write();
h_inv_ALam_ALam_samejet->Write();
h_rot_ALam_ALam_samejet->Write();
hMix_ALam_ALam_samejet->Write();

h_ALam_ALam_diffjet->Write();
h_inv_ALam_ALam_diffjet->Write();
h_rot_ALam_ALam_diffjet->Write();
hMix_ALam_ALam_diffjet->Write();
  
h_ALam_ALam_nojet->Write();
h_inv_ALam_ALam_nojet->Write();
h_rot_ALam_ALam_nojet->Write();
hMix_ALam_ALam_nojet->Write();

h_ALam_ALam_only1jet->Write();
h_inv_ALam_ALam_only1jet->Write();
h_rot_ALam_ALam_only1jet->Write();
hMix_ALam_ALam_only1jet->Write();
  
//LAL

h_Lam_ALam_samejet->Write();
h_inv_Lam_ALam_samejet->Write();
h_rot_Lam_ALam_samejet->Write();
hMix_Lam_ALam_samejet->Write();

h_Lam_ALam_diffjet->Write();
h_inv_Lam_ALam_diffjet->Write();
h_rot_Lam_ALam_diffjet->Write();
hMix_Lam_ALam_diffjet->Write();
  
h_Lam_ALam_nojet->Write();
h_inv_Lam_ALam_nojet->Write();
h_rot_Lam_ALam_nojet->Write();
hMix_Lam_ALam_nojet->Write(); 

h_Lam_ALam_only1jet->Write();
h_inv_Lam_ALam_only1jet->Write();
h_rot_Lam_ALam_only1jet->Write();
hMix_Lam_ALam_only1jet->Write(); 
  
//KL

h_K0s_Lam_samejet->Write();
h_inv_K0s_Lam_samejet->Write();
h_rot_K0s_Lam_samejet->Write();
hMix_K0s_Lam_samejet->Write();

h_K0s_Lam_diffjet->Write();
h_inv_K0s_Lam_diffjet->Write();
h_rot_K0s_Lam_diffjet->Write();
hMix_K0s_Lam_diffjet->Write();
  
h_K0s_Lam_nojet->Write();
h_inv_K0s_Lam_nojet->Write();
h_rot_K0s_Lam_nojet->Write();
hMix_K0s_Lam_nojet->Write(); 

h_K0s_Lam_only1jet->Write();
h_inv_K0s_Lam_only1jet->Write();
h_rot_K0s_Lam_only1jet->Write();
hMix_K0s_Lam_only1jet->Write();

//KAL

h_K0s_ALam_samejet->Write();
h_inv_K0s_ALam_samejet->Write();
h_rot_K0s_ALam_samejet->Write();
hMix_K0s_ALam_samejet->Write(); 

h_K0s_ALam_diffjet->Write();
h_inv_K0s_ALam_diffjet->Write();
h_rot_K0s_ALam_diffjet->Write();
hMix_K0s_ALam_diffjet->Write();
  
h_K0s_ALam_nojet->Write();
h_inv_K0s_ALam_nojet->Write();
h_rot_K0s_ALam_nojet->Write();
hMix_K0s_ALam_nojet->Write(); 

h_K0s_ALam_only1jet->Write();
h_inv_K0s_ALam_only1jet->Write();
h_rot_K0s_ALam_only1jet->Write();
hMix_K0s_ALam_only1jet->Write(); 

// hist for V0 + jet  

h_K0snojet_Jet->Write();
h_K0snojet_Jet_rot->Write();
h_K0snojet_Jet_inv->Write();
hMix_K0snojet_Jet->Write(); 

h_Lamnojet_Jet->Write();
h_Lamnojet_Jet_rot->Write();
h_Lamnojet_Jet_inv->Write();
hMix_Lamnojet_Jet->Write(); 

h_ALamnojet_Jet->Write();
h_ALamnojet_Jet_rot->Write();
h_ALamnojet_Jet_inv->Write();
hMix_ALamnojet_Jet->Write();

h_K0sfromjet_Jet->Write(); 
h_K0sfromjet_Jet_rot->Write();
h_K0sfromjet_Jet_inv->Write();
hMix_K0sfromjet_Jet->Write(); 

h_Lamfromjet_Jet->Write();
h_Lamfromjet_Jet_rot->Write();
h_Lamfromjet_Jet_inv->Write();
hMix_Lamfromjet_Jet->Write();

h_ALamfromjet_Jet->Write();
h_ALamfromjet_Jet_rot->Write();
h_ALamfromjet_Jet_inv->Write();
hMix_ALamfromjet_Jet->Write();

h_K0s_Jet->Write();
h_K0s_Jet_rot->Write();
h_K0s_Jet_inv->Write();
hMix_K0s_Jet->Write(); 

h_Lam_Jet->Write(); 
h_Lam_Jet_rot->Write();
h_Lam_Jet_inv->Write();
hMix_Lam_Jet->Write();

h_ALam_Jet->Write();
h_ALam_Jet_rot->Write();
h_ALam_Jet_inv->Write();
hMix_ALam_Jet->Write();


}

void write_cp(){


K_size->Write();
L_size->Write();
AL_size->Write();

K_size_nocut->Write();
L_size_nocut->Write();
AL_size_nocut->Write();
 
K0smiss_mass_hyp->Write();
K0smiss_mass_ee->Write();
LALmiss_mass_hyp->Write();
LALmiss_mass_ee->Write();

APwocut_K0s->Write();
APwcut_K0s->Write();
APwocut_LAL->Write();
APwcut_LAL->Write();
 
pt_mass_K->Write();
pt_mass_L->Write();
pt_mass_AL->Write();

K0s_pt_reco->Write();
K0s_eta_reco->Write();
K0s_phi_reco->Write();
K0s_mass_reco->Write();

L_pt_reco->Write();
L_eta_reco->Write();
L_phi_reco->Write();
L_mass_reco->Write();

AL_pt_reco->Write();
AL_eta_reco->Write();
AL_phi_reco->Write();
AL_mass_reco->Write();

K0s_pt_reco_d1->Write();
K0s_eta_reco_d1->Write();
K0s_phi_reco_d1->Write();
K0s_mass_reco_d1->Write();

L_pt_reco_d1->Write();
L_eta_reco_d1->Write();
L_phi_reco_d1->Write();
L_mass_reco_d1->Write();

AL_pt_reco_d1->Write();
AL_eta_reco_d1->Write();
AL_phi_reco_d1->Write();
AL_mass_reco_d1->Write();

K0s_pt_reco_d2->Write();
K0s_eta_reco_d2->Write();
K0s_phi_reco_d2->Write();
K0s_mass_reco_d2->Write();

L_pt_reco_d2->Write();
L_eta_reco_d2->Write();
L_phi_reco_d2->Write();
L_mass_reco_d2->Write();

AL_pt_reco_d2->Write();
AL_eta_reco_d2->Write();
AL_phi_reco_d2->Write();
AL_mass_reco_d2->Write();

K0s_pt_gen->Write();
K0s_eta_gen->Write();
K0s_phi_gen->Write();
K0s_mass_gen->Write();

L_pt_gen->Write();
L_eta_gen->Write();
L_phi_gen->Write();
L_mass_gen->Write();

AL_pt_gen->Write();
AL_eta_gen->Write();
AL_phi_gen->Write();
AL_mass_gen->Write();

K0s_pt_mat->Write();
K0s_eta_mat->Write();
K0s_phi_mat->Write();
K0s_mass_mat->Write();

L_pt_mat->Write();
L_eta_mat->Write();
L_phi_mat->Write();
L_mass_mat->Write();

AL_pt_mat->Write();
AL_eta_mat->Write();
AL_phi_mat->Write();
AL_mass_mat->Write();
 
    
}
