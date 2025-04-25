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
#include <stdlib.h>
#include "hbt.h"

using namespace std;

void hbtV0(TString filein, int Evt_mix, float vtxcut, int meth, int multiplicity, int input, int systemat, int mcc, int beamside){

cout << "Welcome to Femtoscopy Analysis Code for Data/MC" << endl;

cout << "This code measure the follow correlations: " << endl;

cout << "K0sK0s, K0sLam, K0sALam, LamLam, ALamALam and LamALam"<< endl;

clock_t sec0, sec1;
sec0 = clock();

sw2();

//how much mothers to remove: method = 0 remove both mothers, method = 1 remove one mother
int method = meth;
string mth;
if(method==0){mth = "remboth";
}else if(method==1){mth = "remonerand";
}else if(method==2){mth = "remonemass";
}else if(method==3){mth = "remoneV0chi2";
}else if(method>3){mth = "donotremove";}

//file X    
string fileX = to_string(input);

// system: 0 = pp, 1 = pPb, 2 >= PbPb
int beam = beamside;
string sss;
if(beam==0){sss = "pPb";}else if(beam==1){sss = "Pbp";}else{sss = "No_trig_weight";}
cout << "Beam: "+sss << endl;
//mult = 0: [0-80], mult = 1: [80,104], mult = 2: [105,129] and mult = 3: [130,...] for pp at 13 TeV
//mult = 0: [0-120], mult = 1: [120,185], mult = 2: [185,250] and mult = 3: [250,400] for pPb ar 8.16 TeV
//mult = 0: [0-50000] for PbPb at 5.02 TeV

int mult = multiplicity;
string mmm = to_string(mult);

//Events to mixing
string mix = to_string(Evt_mix);
string mix_vtx_z = to_string(vtxcut);


string mc_or_data;
bool isMC;
int mc = mcc;
if(mc == 0){isMC = false;}else{isMC = true;}
if(isMC){mc_or_data = "MC";}else{mc_or_data = "Data";}
bool DoFeedDown_Xi = false;
bool DoFeedDown_Om = false;
bool DoFeedDown=false;
if(DoFeedDown_Xi==true || DoFeedDown_Om==true)DoFeedDown=true;

bool remove_all_mis = false;

bool removedoublejet = false;
bool removeV0from2jets = false;
bool no_jets_et_all = false;
bool doveto20gev = false;

bool trig_rew_usebits = true;

bool usematrix = false;

//Read the file

TFile *file = TFile::Open(Form("%s",filein.Data()));
TTree *tree = (TTree *)file->Get("K0SAnalysis/my_tree");
Int_t n_ev = tree->GetEntriesFast(); //total number of events

//Access the TTree
read_tree(tree, isMC);

int p = systemat;
string syst;

syst = get_sys(p);


TString filename;
filename = "./"+syst+"/hist_strange_"+sss+"_"+mc_or_data+"_Mult_"+mmm+"_Mix_"+mix+"_vtx_cut_mix_"+mix_vtx_z+"_Method_"+mth+"_file_"+fileX+".root";


//open output file
TFile *MyFile = new TFile(Form("%s",filename.Data()),"RECREATE");
if ( MyFile->IsOpen() ) cout << "File opened successfully" << endl;
MyFile->cd();

cout << "Total number of events: " << n_ev << endl;
cout << "Systematics: "+syst << endl;
cout << "Number of events to mix: " << Evt_mix << endl;
cout << "Delta Vz: "+mix_vtx_z << endl;

if(isMC){cout << "Running with MC" << endl;}else{cout << "Running with Data" << endl;}

TFile *file_reweight = TFile::Open("/eos/cms/store/group/phys_heavyions/ddesouza/DATA/Vtx_reweight.root");

if(isMC && !usematrix){

if(mc!=3 && mc!=4){
    Vz_w_hist = (TH1D *)file_reweight->Get("Vtx_reweight_epos"); 
    Ntk_w_hist = (TH1D *)file_reweight->Get("Ntk_reweight_epos");
}
if(mc==3){
    Vz_w_hist = (TH1D *)file_reweight->Get("Vtx_reweight_ampt"); 
    Ntk_w_hist = (TH1D *)file_reweight->Get("Ntk_reweight_ampt");    
}
if(mc==4){
    Vz_w_hist = (TH1D *)file_reweight->Get("Vtx_reweight_hijing"); 
    Ntk_w_hist = (TH1D *)file_reweight->Get("Ntk_reweight_hijing");
}

}

if(isMC && usematrix){

if(mc!=3 && mc!=4){
    NtkVz_w_hist = (TH2D *)file_reweight->Get("NtkVz_reweight_epos"); 
}
if(mc==3){
    NtkVz_w_hist = (TH2D *)file_reweight->Get("NtkVz_reweight_ampt"); 
}
if(mc==4){
    NtkVz_w_hist = (TH2D *)file_reweight->Get("NtkVz_reweight_hijing"); 
}    
    
}


double nev = (double) n_ev;
for(Int_t i=0; i<n_ev; i++){//loop in events

tree->GetEntry(i);

nocutev->Fill(1);

if (i!=0 && (i%100000)==0 ){
double alpha = (double) i;
cout << " Running --> i = " << i << ", this is: " << ((alpha/nev) * 100)  << "%" << endl;
}

//apply filters

bool noscramping = false;
bool primaryvertexfilter = false;
bool hffilter = false;
bool PUFilter = false;

if(p < 23){
	if(i==0){cout << "PU filter: " << "olvFilter_pPb8TeV_dz1p0" << endl;}
	if(olvFilter_pPb8TeV_dz1p0 > 0) { PUFilter = true; }
}else if(p == 23){
	if(i==0){cout << "PU filter: " << "pileUpFilter_pPb8TeV_Gplus" << endl;}
	if(pileUpFilter_pPb8TeV_Gplus > 0) { PUFilter = true; }
}else if(p == 24){
	if(i==0){cout << "PU filter: " << "pileUpFilter_pPb8TeV_vtx1" << endl;}
	if(pileUpFilter_pPb8TeV_vtx1 > 0) { PUFilter = true; }
}else if(p == 25){
	if(i==0){cout << "PU filter: " << "No PU filter applied" << endl;}
	PUFilter = true;
}

if(minnTowersTh3HF >= 1){hffilter = true;}
if(validPV > 0 && tracksizePV >= 2 && abs(vtx_z) < 25 && vtx_rho <= 2.0){primaryvertexfilter = true;}
if(trkColl_noscr > 10){
   if(fraction_noscr > 0.25){ noscramping = true; }
}else if(trkColl_noscr <= 10){ noscramping = true; }

if(PUFilter){nocutPU->Fill(1);}
if(primaryvertexfilter){nocutPV->Fill(1);}
if(hffilter){nocutHF->Fill(1);}
if(noscramping){nocutSC->Fill(1);}

if(phffilter_1tw > 0 && pprimaryvertexfilter > 0 && pnoscrampingfilter > 0 && olvFilter_pPb8TeV_dz1p0 > 0){nocutALLX->Fill(1);}

if(mult==0){if(N_tkoff >= 120)continue;}//for pPb
if(mult==1){if(N_tkoff < 120 || N_tkoff >= 150)continue;}//for pPb
if(mult==2){if(N_tkoff < 150 || N_tkoff >= 185)continue;}//for pPb
if(mult==3){if(N_tkoff < 185 || N_tkoff >= 250)continue;}//for pPb
if(mult==4){if(N_tkoff < 250 || N_tkoff >= 500)continue;}//for pPb

if(noscramping && primaryvertexfilter && hffilter && PUFilter){
nocutALL->Fill(1);
//Vertex and Ntrk weight

double Vz_w = 1.0;
double Ntk_w = 1.0;
if(isMC && mc > 1 && !usematrix){
Vz_w = Vz_w_hist->GetBinContent(Vz_w_hist->GetXaxis()->FindBin(vtx_z));
Ntk_w = Ntk_w_hist->GetBinContent(Ntk_w_hist->GetXaxis()->FindBin(N_tkoff));
}

if(isMC && mc > 1 && usematrix){
Ntk_w = Ntk_w_hist->GetBinContent(Ntk_w_hist->GetXaxis()->FindBin(vtx_z),Ntk_w_hist->GetYaxis()->FindBin(N_tkoff));
}

//trigger weight

double SF_lumi = 1.0;
double SF_bits = 1.0;

double factlumibeam0 = 1.0;
double factlumibeam1 = 1.0;
double factbitsbeam0 = 1.0;
double factbitsbeam1 = 1.0;


if(!isMC && beam==0){
if(mult==0){SF_lumi=56.3291/factlumibeam0; SF_bits=59.3841/factbitsbeam0;}//MB1to8 -> divided by 8 because # of samples
if(mult==1){SF_lumi=66.8549/factlumibeam0; SF_bits=70.4712/factbitsbeam0;}//HM0120
if(mult==2){SF_lumi=24.8170/factlumibeam0; SF_bits=25.5365/factbitsbeam0;}//HM0150
if(mult==3){SF_lumi=2.18724/factlumibeam0; SF_bits=2.22348/factbitsbeam0;}//HM1to6 -> divided by 6 because # of samples
if(mult==4){SF_lumi=1.0/factlumibeam0; SF_bits=1.0/factbitsbeam0;}//HM7 -> unprescaled
}else if(!isMC && beam==1){
if(mult==0){SF_lumi=37.1378/factlumibeam1; SF_bits=35.222/factbitsbeam1;}//MB1to20 -> divided by 20 because # of samples
if(mult==1){SF_lumi=65.7231/factlumibeam1; SF_bits=70.0289/factbitsbeam1;}//HM0120
if(mult==2){SF_lumi=34.6225/factlumibeam1; SF_bits=35.0635/factbitsbeam1;}//HM0150
if(mult==3){SF_lumi=1.64455/factlumibeam1; SF_bits=1.712/factbitsbeam1;}//HM1to6 -> divided by 6 because # of samples
if(mult==4){SF_lumi=1.0/factlumibeam1; SF_bits=1.0/factbitsbeam1;}//HM7 -> unprescaled
}else if(!isMC && beam>1){SF_lumi=1.0; SF_bits=1.0;}//cout << "No trigger rewight" << endl;}

if(isMC){SF_lumi=1.0; SF_bits=1.0;}

double Ntk_Vz_weight;
if(trig_rew_usebits){Ntk_Vz_weight = Vz_w*Ntk_w*SF_bits;}else{Ntk_Vz_weight = Vz_w*Ntk_w*SF_lumi;}

Int_t gk0size;
Int_t gLamsize;


//general cuts
if(vtx_z >= vtxzcutmax ||  vtx_z <= vtxzcutmin)continue;
if(p==14){if(vtx_z <= 3.0 && vtx_z >= -3.0)continue;}
if(!isMC){if(vtx_rho > vtxrhocut)continue;}

h_vtx_z->Fill(vtx_z,Vz_w);
h_vtx_rho->Fill(vtx_rho,Vz_w);
h_ntrk_cent->Fill(N_tkoff,Ntk_Vz_weight);

Ntrk_VZ->Fill(vtx_z,N_tkoff);

Int_t k0size = K0s_pt->size();
Int_t Lamsize = Lam_pt->size();
Int_t Xisize = Xi_pt->size();
Int_t Omsize = Om_pt->size();

Int_t jetsize;
if(isMC){

jetsize = Jet_pt->size();
if(no_jets_et_all){if(jetsize!=0)continue;} //events with

//veto
bool vetob = false;
if(doveto20gev){
    if(jetsize >= 1){
        for(Int_t aa=0; aa<jetsize; aa++){ 
            if(Jet_pt->at(aa)>=20.){
                vetob=true; 
                break;
            }
        }
    }
    if(vetob)continue;
}

}

//cout << "Xisize: " << Xisize << endl;
//cout << "Omsize: " << Omsize << endl;

Int_t aux_N_tk_offline = N_tkoff;
Double_t aux_vtxz = vtx_z;

std::vector<TLorentzVector> map_GenK0s; 
std::vector<TLorentzVector> map_GenK0s_matched; 
std::vector<TLorentzVector> map_GenK0s_matchedT; 
std::vector<TLorentzVector> map_GenLam; 
std::vector<TLorentzVector> map_GenLam_matched; 
std::vector<TLorentzVector> map_GenLam_matchedT; 
std::vector<TLorentzVector> map_GenALam; 
std::vector<TLorentzVector> map_GenALam_matched; 
std::vector<TLorentzVector> map_GenALam_matchedT; 

if(isMC){ //DONE

gk0size = gK0s_pt->size();
gLamsize = gLam_pt->size();

//Gen stuff ==========================

if(gk0size>=1){
for(Int_t aj=0; aj<gk0size; aj++){ 
if(gK0s_stat->at(aj) > 1.5 || gK0s_stat->at(aj) < 0.5){continue;}
if(gK0s_pt->at(aj) <= ptminK0s || gK0s_pt->at(aj) >= ptmaxK0s){continue;}
if(fabs(gK0s_eta->at(aj)) > 2.4){continue;}
double mid1 = gK0s_mom1->at(aj);
double mid2 = gK0s_mom2->at(aj);
if(gK0s_mass->at(aj) < 0.496){continue;}
if(gK0s_mass->at(aj) > 0.498){continue;}
TLorentzVector pvectors;
pvectors.SetPtEtaPhiM(gK0s_pt->at(aj),gK0s_eta->at(aj),gK0s_phi->at(aj),gK0s_mass->at(aj));
if(fabs(pvectors.Rapidity()) >= rapvar)continue;
map_GenK0s.push_back(pvectors); 
}}

if(gLamsize>=1){
for(Int_t ak=0; ak<gLamsize; ak++){ 
if(gLam_stat->at(ak) > 1.5 || gLam_stat->at(ak) < 0.5){continue;}
if(gLam_pt->at(ak) <= ptminLAL || gLam_pt->at(ak) >= ptmaxLAL)continue;
if(fabs(gLam_eta->at(ak)) > 2.4)continue;
if(gLam_mass->at(ak) < massLAL - nsigmapeak*sigmaLAL)continue;
if(gLam_mass->at(ak) > massLAL + nsigmapeak*sigmaLAL)continue;
double mid1 = gLam_mom1->at(ak);
double mid2 = gLam_mom2->at(ak);
TLorentzVector pvectors;
pvectors.SetPtEtaPhiM(gLam_pt->at(ak),gLam_eta->at(ak),gLam_phi->at(ak),gLam_mass->at(ak));
if(fabs(pvectors.Rapidity()) >= rapvar)continue;
if(gLam_id->at(ak) > 0){
map_GenLam.push_back(pvectors);
}else if(gLam_id->at(ak) < 0){
map_GenALam.push_back(pvectors); 
}}}
//====================================

cplots(map_GenK0s,K0s_pt_gen,K0s_eta_gen,K0s_phi_gen,K0s_mass_gen,Ntk_Vz_weight);
cplots(map_GenLam,L_pt_gen,L_eta_gen,L_phi_gen,L_mass_gen,Ntk_Vz_weight);
cplots(map_GenALam,AL_pt_gen,AL_eta_gen,AL_phi_gen,AL_mass_gen,Ntk_Vz_weight);

}


if(i==0)cout << "Start to work with RECO"  << endl;  
//Reco stuff ========================== 

////////////////////////////////////K0s/////////////////////////////////////////////////////////////////////////////////////

//4-vectors for no chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectord1_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2_K0s;
std::vector<TLorentzVector> GoodTrackFourVector_K0s;
std::vector<double> chi2_K0s;
std::vector<double> chi21_K0s;
std::vector<double> chi22_K0s;
std::vector<Bool_t> map_K0s_chi2;
std::vector<Bool_t> map_K0s_match;

//4-vectors for no chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectord1_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2_Lam;
std::vector<TLorentzVector> GoodTrackFourVector_Lam;
std::vector<double> chi2_Lam;
std::vector<double> chi21_Lam;
std::vector<double> chi22_Lam;
std::vector<Bool_t> map_Lam_chi2;
std::vector<Bool_t> map_Lam_match;


std::vector<TLorentzVector> GoodTrackFourVectord1_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2_ALam;
std::vector<TLorentzVector> GoodTrackFourVector_ALam;
std::vector<double> chi2_ALam;
std::vector<double> chi21_ALam;
std::vector<double> chi22_ALam;
std::vector<Bool_t> map_ALam_chi2;
std::vector<Bool_t> map_ALam_match;

std::vector<TLorentzVector> GoodTrackFourVector_Jet;

//for non-identical

std::vector<Bool_t> map_Lam_chi2_LAL;
std::vector<Bool_t> map_ALam_chi2_LAL;

std::vector<Bool_t> map_K0s_chi2_KL;
std::vector<Bool_t> map_Lam_chi2_KL;

std::vector<Bool_t> map_K0s_chi2_KAL;
std::vector<Bool_t> map_ALam_chi2_KAL;

//for feed down
std::vector<Bool_t> lamb_from_feed_xi_K;
std::vector<Bool_t> lamb_from_feed_axi_K;
std::vector<Bool_t> alamb_from_feed_xi_K;
std::vector<Bool_t> alamb_from_feed_axi_K;

std::vector<Bool_t> lamb_from_feed_om_K;
std::vector<Bool_t> lamb_from_feed_aom_K;
std::vector<Bool_t> alamb_from_feed_om_K;
std::vector<Bool_t> alamb_from_feed_aom_K;

std::vector<Bool_t> lamb_from_feed_xi_L;
std::vector<Bool_t> lamb_from_feed_axi_L;
std::vector<Bool_t> alamb_from_feed_xi_L;
std::vector<Bool_t> alamb_from_feed_axi_L;

std::vector<Bool_t> lamb_from_feed_om_L;
std::vector<Bool_t> lamb_from_feed_aom_L;
std::vector<Bool_t> alamb_from_feed_om_L;
std::vector<Bool_t> alamb_from_feed_aom_L;

std::vector<Bool_t> lamb_from_feed_xi_AL;
std::vector<Bool_t> lamb_from_feed_axi_AL;
std::vector<Bool_t> alamb_from_feed_xi_AL;
std::vector<Bool_t> alamb_from_feed_axi_AL;

std::vector<Bool_t> lamb_from_feed_om_AL;
std::vector<Bool_t> lamb_from_feed_aom_AL;
std::vector<Bool_t> alamb_from_feed_om_AL;
std::vector<Bool_t> alamb_from_feed_aom_AL;

std::vector<TLorentzVector> GoodTrackFourVector_K0s_Jet;
std::vector<Int_t> GoodTrackFourVector_K0s_NJet;
std::vector<TLorentzVector> GoodTrackFourVector_K0s_NotaJet;

std::vector<TLorentzVector> GoodTrackFourVector_Lam_Jet;
std::vector<Int_t> GoodTrackFourVector_Lam_NJet;
std::vector<TLorentzVector> GoodTrackFourVector_Lam_NotaJet;

std::vector<TLorentzVector> GoodTrackFourVector_ALam_Jet;
std::vector<Int_t> GoodTrackFourVector_ALam_NJet;
std::vector<TLorentzVector> GoodTrackFourVector_ALam_NotaJet;

int contK = 0;

for(Int_t j=0; j<k0size; j++){ 

if(contK==0){nev_K0s_ini->Fill(1);}
contK=contK+1;
       
if(K0s_pt->at(j) <= ptminK0s || K0s_pt->at(j) >= ptmaxK0s)continue;
if(fabs(K0s_eta->at(j)) > 2.4)continue;

if(fabs(K0s_dxy1->at(j)) < dxyz)continue;
if(fabs(K0s_dxy2->at(j)) < dxyz)continue;
if(fabs(K0s_dz1->at(j)) < dxyz)continue;
if(fabs(K0s_dz2->at(j)) < dxyz)continue;
if(K0s_d1Nhit->at(j) < nhitss)continue;
if(K0s_d2Nhit->at(j) < nhitss)continue;
if(K0s_d1Pix->at(j) < pixelhits)continue;
if(K0s_d2Pix->at(j) < pixelhits)continue;
if(K0s_3Dagl->at(j) <= cospoint)continue;
if(K0s_3Ddl->at(j) <= decaylength)continue;
if(K0s_dca->at(j) > trkdca)continue;

if(K0s_ctau->at(j) <= ctaucutmin || K0s_ctau->at(j) >= ctaucutmax)continue;

double Pxp = K0s_d1px->at(j);
double Pyp = K0s_d1py->at(j);
double Pzp = K0s_d1pz->at(j);
double Pmp = K0s_d1M->at(j);

double Pxn = K0s_d2px->at(j);
double Pyn = K0s_d2py->at(j);
double Pzn = K0s_d2pz->at(j);
double Pmn = K0s_d2M->at(j);

TLorentzVector pvector1;
pvector1.SetXYZM(Pxp,Pyp,Pzp,Pmp);
TLorentzVector pvector2;
pvector2.SetXYZM(Pxn,Pyn,Pzn,Pmn);

TVector3 dauvec1(pvector1.Px(),pvector1.Py(),pvector1.Pz());
TVector3 dauvec2(pvector2.Px(),pvector2.Py(),pvector2.Pz());
TVector3 dauvecsum(dauvec1+dauvec2);

// Armenteros-Podolanski
 double Pp=sqrt(Pxp*Pxp+Pyp*Pyp+Pzp*Pzp);
 double Pn=sqrt(Pxn*Pxn+Pyn*Pyn+Pzn*Pzn);
 double PN=sqrt((Pxp+Pxn)*(Pxp+Pxn)+(Pyp+Pyn)*(Pyp+Pyn)+(Pzp+Pzn)*(Pzp+Pzn));
 double cosx1=(Pxp*(Pxp+Pxn)+Pyp*(Pyp+Pyn)+Pzp*(Pzp+Pzn))/(Pp*PN);
 double cosx2=(Pxn*(Pxp+Pxn)+Pyn*(Pyp+Pyn)+Pzn*(Pzp+Pzn))/(Pn*PN);
 double QT = Pp*sqrt(1-cosx1*cosx1);
 double Alpha;
 Alpha = (Pp*cosx1-Pn*cosx2)/(Pp*cosx1+Pn*cosx2);
 APwocut_K0s->Fill(Alpha,QT,Ntk_Vz_weight);
 
double energyd1e = sqrt(electronMass*electronMass+pvector1.P()*pvector1.P());
double energyd2e = sqrt(electronMass*electronMass+pvector2.P()*pvector2.P());
double invmass_ee = sqrt((energyd1e+energyd2e)*(energyd1e+energyd2e)-dauvecsum.Mag2());

K0smiss_mass_ee->Fill(invmass_ee,Ntk_Vz_weight);

if(invmass_ee<misee)continue;// do not change anymore

double massd1, massd2, energyd1, energyd2, invmass;

massd1=pimass;
massd2=pro_mass;
energyd1 = sqrt(massd1*massd1+pvector1.P()*pvector1.P());
energyd2 = sqrt(massd2*massd2+pvector2.P()*pvector2.P());
invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());

K0smiss_mass_hyp->Fill(invmass,Ntk_Vz_weight);

if(fabs(invmass-la_mass)<misK0s) continue;//ok

massd1=pro_mass;
massd2=pimass;
energyd1 = sqrt(massd1*massd1+pvector1.P()*pvector1.P());
energyd2 = sqrt(massd2*massd2+pvector2.P()*pvector2.P());
invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
if(fabs(invmass-la_mass)<misK0s) continue;// do not change anymore

TLorentzVector pvector;
pvector.SetPtEtaPhiM(K0s_pt->at(j),K0s_eta->at(j),K0s_phi->at(j),K0s_mass->at(j));
if(fabs(pvector.Rapidity()) >= rapvar)continue;

APwcut_K0s->Fill(Alpha,QT,Ntk_Vz_weight);

chi2_K0s.push_back(K0s_vtx->at(j));
GoodTrackFourVector_K0s.push_back(pvector);
GoodTrackFourVectord1_K0s.push_back(pvector1);
GoodTrackFourVectord2_K0s.push_back(pvector2);
chi21_K0s.push_back(K0s_chi21->at(j));
chi22_K0s.push_back(K0s_chi22->at(j));
map_K0s_chi2.push_back(kTRUE);
map_K0s_match.push_back(kFALSE);
map_K0s_chi2_KL.push_back(kTRUE);
map_K0s_chi2_KAL.push_back(kTRUE);

K0sctau->Fill(K0s_ctau->at(j),Ntk_Vz_weight);
K0sdca3D->Fill(K0s_3Ddca->at(j),Ntk_Vz_weight);
K0sDL->Fill(K0s_3Ddl->at(j),Ntk_Vz_weight);
K0strkdca->Fill(K0s_dca->at(j),Ntk_Vz_weight);
K0scostheta->Fill(K0s_3Dagl->at(j),Ntk_Vz_weight);
K0s_Vtx->Fill(K0s_vtx->at(j),Ntk_Vz_weight);

K0s_d1_dxy->Fill(fabs(K0s_dxy1->at(j)),Ntk_Vz_weight);
K0s_d1_dz->Fill(fabs(K0s_dz1->at(j)),Ntk_Vz_weight);
K0s_d2_dxy->Fill(fabs(K0s_dxy2->at(j)),Ntk_Vz_weight);
K0s_d2_dz->Fill(fabs(K0s_dz2->at(j)),Ntk_Vz_weight);
K0s_d1_Nhits->Fill(fabs(K0s_d1Nhit->at(j)),Ntk_Vz_weight);
K0s_d2_Nhits->Fill(fabs(K0s_d2Nhit->at(j)),Ntk_Vz_weight);
pt_mass_K->Fill(pvector.Pt(), pvector.M(),Ntk_Vz_weight);


}//K0s selection

K_size_nocut->Fill(contK,Ntk_Vz_weight);


int contL = 0;
int contAL = 0;

for(Int_t k=0; k<Lamsize; k++){ 

if(Lam_id->at(k)>0 && contL==0){nev_Lam_ini->Fill(1);}
contL=contL+1;
if(Lam_id->at(k)<0 && contAL==0){nev_ALam_ini->Fill(1);}
contAL=contAL+1;
    
if(Lam_pt->at(k) <= ptminLAL || Lam_pt->at(k) >= ptmaxLAL)continue;
if(fabs(Lam_eta->at(k)) > 2.4)continue;

if(fabs(Lam_dxy1->at(k)) < dxyz)continue;
if(fabs(Lam_dxy2->at(k)) < dxyz)continue;
if(fabs(Lam_dz1->at(k)) < dxyz)continue;
if(fabs(Lam_dz2->at(k)) < dxyz)continue;
if(Lam_d1Nhit->at(k) < nhitss)continue;
if(Lam_d2Nhit->at(k) < nhitss)continue;
if(Lam_d1Pix->at(k) < pixelhits)continue;
if(Lam_d2Pix->at(k) < pixelhits)continue;
if(Lam_3Dagl->at(k) <= cospoint)continue;
if(Lam_3Ddl->at(k) <= decaylength)continue;
if(Lam_dca->at(k) > trkdca)continue;

if(Lam_ctau->at(k) <= ctaucutmin || Lam_ctau->at(k) >= ctaucutmax)continue;

double Pxp = Lam_d1px->at(k);
double Pyp = Lam_d1py->at(k);
double Pzp = Lam_d1pz->at(k);
double Pmp = Lam_d1M->at(k);

double Pxn = Lam_d2px->at(k);
double Pyn = Lam_d2py->at(k);
double Pzn = Lam_d2pz->at(k);
double Pmn = Lam_d2M->at(k);

TLorentzVector pvector1;
pvector1.SetXYZM(Pxp,Pyp,Pzp,Pmp);
TLorentzVector pvector2;
pvector2.SetXYZM(Pxn,Pyn,Pzn,Pmn);

TVector3 dauvec1(pvector1.Px(),pvector1.Py(),pvector1.Pz());
TVector3 dauvec2(pvector2.Px(),pvector2.Py(),pvector2.Pz());
TVector3 dauvecsum(dauvec1+dauvec2);

// Armenteros-Podolanski
 double Pp=sqrt(Pxp*Pxp+Pyp*Pyp+Pzp*Pzp);
 double Pn=sqrt(Pxn*Pxn+Pyn*Pyn+Pzn*Pzn);
 double PN=sqrt((Pxp+Pxn)*(Pxp+Pxn)+(Pyp+Pyn)*(Pyp+Pyn)+(Pzp+Pzn)*(Pzp+Pzn));
 double cosx1=(Pxp*(Pxp+Pxn)+Pyp*(Pyp+Pyn)+Pzp*(Pzp+Pzn))/(Pp*PN);
 double cosx2=(Pxn*(Pxp+Pxn)+Pyn*(Pyp+Pyn)+Pzn*(Pzp+Pzn))/(Pn*PN);
 double QT = Pp*sqrt(1-cosx1*cosx1);
 double Alpha;
 if(Lam_id->at(k) < 0){Alpha = (-Pp*cosx1+Pn*cosx2)/(Pp*cosx1+Pn*cosx2);}
 if(Lam_id->at(k) > 0){Alpha = (Pp*cosx1-Pn*cosx2)/(Pp*cosx1+Pn*cosx2);}
 APwocut_LAL->Fill(Alpha,QT,Ntk_Vz_weight);

double energyd1e = sqrt(electronMass*electronMass+pvector1.P()*pvector1.P());
double energyd2e = sqrt(electronMass*electronMass+pvector2.P()*pvector2.P());
double invmass_ee = sqrt((energyd1e+energyd2e)*(energyd1e+energyd2e)-dauvecsum.Mag2());
LALmiss_mass_ee->Fill(invmass_ee,Ntk_Vz_weight);

if(invmass_ee<misee)continue;// do not change anymore

double massd1, massd2, energyd1, energyd2, invmass;

massd1=pimass;
massd2=pimass;
energyd1 = sqrt(massd1*massd1+pvector1.P()*pvector1.P());
energyd2 = sqrt(massd2*massd2+pvector2.P()*pvector2.P());
invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
LALmiss_mass_hyp->Fill(invmass,Ntk_Vz_weight);

if(fabs(invmass-k0s_mass)<misLam) continue;//ok

TLorentzVector pvector;
pvector.SetPtEtaPhiM(Lam_pt->at(k),Lam_eta->at(k),Lam_phi->at(k),Lam_mass->at(k));
if(fabs(pvector.Rapidity()) >= rapvar)continue;

APwcut_LAL->Fill(Alpha,QT,Ntk_Vz_weight);

pT_Lamb_from_prompt_TF->Fill(Lam_pt->at(k));
pT_Lamb->Fill(Lam_pt->at(k),Lam_mass->at(k));

if(Lam_id->at(k) > 0){
chi2_Lam.push_back(Lam_vtx->at(k));
GoodTrackFourVector_Lam.push_back(pvector);
GoodTrackFourVectord1_Lam.push_back(pvector1);
GoodTrackFourVectord2_Lam.push_back(pvector2);
chi21_Lam.push_back(Lam_chi21->at(k));
chi22_Lam.push_back(Lam_chi22->at(k));     
map_Lam_chi2.push_back(kTRUE);
map_Lam_match.push_back(kFALSE);
map_Lam_chi2_LAL.push_back(kTRUE);
map_Lam_chi2_KL.push_back(kTRUE);


Lamctau->Fill(Lam_ctau->at(k),Ntk_Vz_weight);
Lamdca3D->Fill(Lam_3Ddca->at(k),Ntk_Vz_weight);
LamDL->Fill(Lam_3Ddl->at(k),Ntk_Vz_weight);
Lamtrkdca->Fill(Lam_dca->at(k),Ntk_Vz_weight);
Lamcostheta->Fill(Lam_3Dagl->at(k),Ntk_Vz_weight);
Lam_Vtx->Fill(Lam_vtx->at(k),Ntk_Vz_weight);

Lam_d1_dxy->Fill(fabs(Lam_dxy1->at(k)),Ntk_Vz_weight);
Lam_d1_dz->Fill(fabs(Lam_dz1->at(k)),Ntk_Vz_weight);
Lam_d2_dxy->Fill(fabs(Lam_dxy2->at(k)),Ntk_Vz_weight);
Lam_d2_dz->Fill(fabs(Lam_dz2->at(k)),Ntk_Vz_weight);
Lam_d1_Nhits->Fill(fabs(Lam_d1Nhit->at(k)),Ntk_Vz_weight);
Lam_d2_Nhits->Fill(fabs(Lam_d2Nhit->at(k)),Ntk_Vz_weight);

lamb_from_feed_xi_K.push_back(kFALSE);
lamb_from_feed_axi_K.push_back(kFALSE);
lamb_from_feed_om_K.push_back(kFALSE);
lamb_from_feed_aom_K.push_back(kFALSE);

lamb_from_feed_xi_L.push_back(kFALSE);
lamb_from_feed_axi_L.push_back(kFALSE);
lamb_from_feed_om_L.push_back(kFALSE);
lamb_from_feed_aom_L.push_back(kFALSE);

lamb_from_feed_xi_AL.push_back(kFALSE);
lamb_from_feed_axi_AL.push_back(kFALSE);
lamb_from_feed_om_AL.push_back(kFALSE);
lamb_from_feed_aom_AL.push_back(kFALSE);

pt_mass_L->Fill(pvector.Pt(), pvector.M(),Ntk_Vz_weight);


}else if(Lam_id->at(k) < 0){
chi2_ALam.push_back(Lam_vtx->at(k));
GoodTrackFourVector_ALam.push_back(pvector);
GoodTrackFourVectord1_ALam.push_back(pvector1);
GoodTrackFourVectord2_ALam.push_back(pvector2);
chi21_ALam.push_back(Lam_chi21->at(k));
chi22_ALam.push_back(Lam_chi22->at(k));
map_ALam_chi2.push_back(kTRUE);
map_ALam_match.push_back(kFALSE);
map_ALam_chi2_LAL.push_back(kTRUE);
map_ALam_chi2_KAL.push_back(kTRUE);

ALamctau->Fill(Lam_ctau->at(k),Ntk_Vz_weight);
ALamdca3D->Fill(Lam_3Ddca->at(k),Ntk_Vz_weight);
ALamDL->Fill(Lam_3Ddl->at(k),Ntk_Vz_weight);
ALamtrkdca->Fill(Lam_dca->at(k),Ntk_Vz_weight);
ALamcostheta->Fill(Lam_3Dagl->at(k),Ntk_Vz_weight);
ALam_Vtx->Fill(Lam_vtx->at(k),Ntk_Vz_weight);

ALam_d1_dxy->Fill(fabs(Lam_dxy1->at(k)),Ntk_Vz_weight);
ALam_d1_dz->Fill(fabs(Lam_dz1->at(k)),Ntk_Vz_weight);
ALam_d2_dxy->Fill(fabs(Lam_dxy2->at(k)),Ntk_Vz_weight);
ALam_d2_dz->Fill(fabs(Lam_dz2->at(k)),Ntk_Vz_weight);
ALam_d1_Nhits->Fill(fabs(Lam_d1Nhit->at(k)),Ntk_Vz_weight);
ALam_d2_Nhits->Fill(fabs(Lam_d2Nhit->at(k)),Ntk_Vz_weight);

alamb_from_feed_xi_K.push_back(kFALSE);
alamb_from_feed_axi_K.push_back(kFALSE);
alamb_from_feed_om_K.push_back(kFALSE);
alamb_from_feed_aom_K.push_back(kFALSE);

alamb_from_feed_xi_L.push_back(kFALSE);
alamb_from_feed_axi_L.push_back(kFALSE);
alamb_from_feed_om_L.push_back(kFALSE);
alamb_from_feed_aom_L.push_back(kFALSE);

alamb_from_feed_xi_AL.push_back(kFALSE);
alamb_from_feed_axi_AL.push_back(kFALSE);
alamb_from_feed_om_AL.push_back(kFALSE);
alamb_from_feed_aom_AL.push_back(kFALSE);

pt_mass_AL->Fill(pvector.Pt(), pvector.M(),Ntk_Vz_weight);


}

}//Lam/ALam selection

L_size_nocut->Fill(contL,Ntk_Vz_weight);
AL_size_nocut->Fill(contAL,Ntk_Vz_weight);

//call Xi and Om for feed down studies

std::vector<TLorentzVector> GoodTrackFourVector_Xi;
std::vector<TLorentzVector> GoodTrackFourVector_AXi;
std::vector<double> GoodTrackFourVector_Xi_pion_chi2;
std::vector<double> GoodTrackFourVector_AXi_pion_chi2;
std::vector<Bool_t> GoodTrackFourVector_Xi_bool;
std::vector<Bool_t> GoodTrackFourVector_AXi_bool;


std::vector<TLorentzVector> GoodTrackFourVector_Xi_Lamb;
std::vector<TLorentzVector> GoodTrackFourVector_Xi_Lamb_proton;
std::vector<TLorentzVector> GoodTrackFourVector_Xi_Lamb_pion;

std::vector<TLorentzVector> GoodTrackFourVector_AXi_Lamb;
std::vector<TLorentzVector> GoodTrackFourVector_AXi_Lamb_proton;
std::vector<TLorentzVector> GoodTrackFourVector_AXi_Lamb_pion;

std::vector<double> GoodTrackFourVector_Xi_Lamb_proton_chi2;
std::vector<double> GoodTrackFourVector_Xi_Lamb_pion_chi2;
std::vector<double> GoodTrackFourVector_AXi_Lamb_proton_chi2;
std::vector<double> GoodTrackFourVector_AXi_Lamb_pion_chi2;

std::vector<Bool_t> GoodTrackFourVector_Xi_Lamb_bool_K;
std::vector<Bool_t> GoodTrackFourVector_Xi_Lamb_bool_L;
std::vector<Bool_t> GoodTrackFourVector_Xi_Lamb_bool_AL;
std::vector<Bool_t> GoodTrackFourVector_AXi_Lamb_bool_K;
std::vector<Bool_t> GoodTrackFourVector_AXi_Lamb_bool_L;
std::vector<Bool_t> GoodTrackFourVector_AXi_Lamb_bool_AL;


if(DoFeedDown_Xi){


//Xi + AXi ------------------------
for(Int_t qq=0; qq<Xisize; qq++){ 
if(Xi_pt->at(qq) <= ptminXAX || Xi_pt->at(qq) >= ptmaxXAX)continue;
if(fabs(Xi_eta->at(qq)) > 2.4)continue;
if(Xi_VTrkP3DIpSigValue->at(qq) <= Xi_VTrkP3DIpSigValue_val)continue;
if(Xi_VTrkPi3DIpSigValue->at(qq) <= Xi_VTrkPi3DIpSigValue_val)continue;
if(Xi_casPi3DIpSigValue->at(qq) <= Xi_casPi3DIpSigValue_val)continue;
if(Xi_cas3DIpSigValue->at(qq) >= Xi_cas3DIpSigValue_val)continue;
if(Xi_distanceSigValue->at(qq) <= Xi_distanceSigValue_val)continue;
if(Xi_casFlightSigValue->at(qq) <= Xi_casFlightSigValue_val)continue;
//daughter 2: pion
if(Xi_d2Pix->at(qq) < pixelhits)continue;
//daughter 1: Lambda daughters
//  -> proton cuts 
if(Xi_d1Pix_1->at(qq) < pixelhits)continue;
//  -> pion cuts 
if(Xi_d1Pix_2->at(qq) < pixelhits)continue; 

TLorentzVector pvector1; //lambda
pvector1.SetPtEtaPhiM(Xi_d1pt->at(qq),Xi_d1eta->at(qq),Xi_d1phi->at(qq),Xi_d1mass->at(qq));
TLorentzVector pvector2; //pion
pvector2.SetPtEtaPhiM(Xi_d2pt->at(qq),Xi_d2eta->at(qq),Xi_d2phi->at(qq),Xi_d2mass->at(qq));

double Pxp = pvector1.Px();
double Pyp = pvector1.Py();
double Pzp = pvector1.Pz();
double Pmp = pvector1.M();
double Pxn = pvector2.Px();
double Pyn = pvector2.Py();
double Pzn = pvector2.Pz();
double Pmn = pvector2.M();

TVector3 dauvec1(pvector1.Px(),pvector1.Py(),pvector1.Pz());
TVector3 dauvec2(pvector2.Px(),pvector2.Py(),pvector2.Pz());
TVector3 dauvecsum(dauvec1+dauvec2);


// Armenteros-Podolanski
 double Pp=sqrt(Pxp*Pxp+Pyp*Pyp+Pzp*Pzp);
 double Pn=sqrt(Pxn*Pxn+Pyn*Pyn+Pzn*Pzn);
 double PN=sqrt((Pxp+Pxn)*(Pxp+Pxn)+(Pyp+Pyn)*(Pyp+Pyn)+(Pzp+Pzn)*(Pzp+Pzn));
 double cosx1=(Pxp*(Pxp+Pxn)+Pyp*(Pyp+Pyn)+Pzp*(Pzp+Pzn))/(Pp*PN);
 double cosx2=(Pxn*(Pxp+Pxn)+Pyn*(Pyp+Pyn)+Pzn*(Pzp+Pzn))/(Pn*PN);
 double QT = Pp*sqrt(1-cosx1*cosx1);
 double Alpha;
 if(Xi_id->at(qq) < 0){Alpha = (-Pp*cosx1+Pn*cosx2)/(Pp*cosx1+Pn*cosx2);}
 if(Xi_id->at(qq) > 0){Alpha = (Pp*cosx1-Pn*cosx2)/(Pp*cosx1+Pn*cosx2);}
 
double energyd1e = sqrt(electronMass*electronMass+pvector1.P()*pvector1.P());
double energyd2e = sqrt(electronMass*electronMass+pvector2.P()*pvector2.P());
double invmass_ee = sqrt((energyd1e+energyd2e)*(energyd1e+energyd2e)-dauvecsum.Mag2());
if(invmass_ee<misee)continue;// do not change anymore

double massd1, massd2, energyd1, energyd2, invmass;
massd1=kaon_mass;
massd2=la_mass;
energyd1 = sqrt(massd1*massd1+pvector1.P()*pvector1.P());
energyd2 = sqrt(massd2*massd2+pvector2.P()*pvector2.P());
invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
if(fabs(invmass-om_mass)<misXi) continue;//ok
massd2=kaon_mass;
massd1=la_mass;
energyd1 = sqrt(massd1*massd1+pvector1.P()*pvector1.P());
energyd2 = sqrt(massd2*massd2+pvector2.P()*pvector2.P());
invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
if(fabs(invmass-om_mass)<misXi) continue;// do not change anymore

TLorentzVector pvector;
pvector.SetPtEtaPhiM(Xi_pt->at(qq),Xi_eta->at(qq),Xi_phi->at(qq),Xi_mass->at(qq));
if(fabs(pvector.Rapidity()) >= rapvar)continue;


//AP for lambda daughter
TLorentzVector pvector1_l; //proton
pvector1_l.SetPtEtaPhiM(Xi_d1pt_1->at(qq),Xi_d1eta_1->at(qq),Xi_d1phi_1->at(qq),Xi_d1mass_1->at(qq));
TLorentzVector pvector2_l; //pion
pvector2_l.SetPtEtaPhiM(Xi_d1pt_2->at(qq),Xi_d1eta_2->at(qq),Xi_d1phi_2->at(qq),Xi_d1mass_2->at(qq));

double Pxp_l = pvector1_l.Px();
double Pyp_l = pvector1_l.Py();
double Pzp_l = pvector1_l.Pz();
double Pmp_l = pvector1_l.M();

double Pxn_l = pvector2_l.Px();
double Pyn_l = pvector2_l.Py();
double Pzn_l = pvector2_l.Pz();
double Pmn_l = pvector2_l.M();

TVector3 dauvec1_l(pvector1_l.Px(),pvector1_l.Py(),pvector1_l.Pz());
TVector3 dauvec2_l(pvector2_l.Px(),pvector2_l.Py(),pvector2_l.Pz());
TVector3 dauvecsum_l(dauvec1_l+dauvec2_l);
 
double energyd1e_l = sqrt(electronMass*electronMass+pvector1_l.P()*pvector1_l.P());
double energyd2e_l = sqrt(electronMass*electronMass+pvector2_l.P()*pvector2_l.P());
double invmass_ee_l = sqrt((energyd1e_l+energyd2e_l)*(energyd1e_l+energyd2e_l)-dauvecsum_l.Mag2());
if(invmass_ee_l<misee)continue;// do not change anymore

double massd1_l, massd2_l, energyd1_l, energyd2_l, invmass_l;
massd1_l=pimass;
massd2_l=pimass;
energyd1_l = sqrt(massd1_l*massd1_l+pvector1_l.P()*pvector1_l.P());
energyd2_l = sqrt(massd2_l*massd2_l+pvector2_l.P()*pvector2_l.P());
invmass_l = sqrt((energyd1_l+energyd2_l)*(energyd1_l+energyd2_l)-dauvecsum_l.Mag2());
if(fabs(invmass_l-k0s_mass)<misLam) continue;//ok


//lambda daughters
//proton
TLorentzVector lvector1;
lvector1.SetPtEtaPhiM(Xi_d1pt_1->at(qq),Xi_d1eta_1->at(qq),Xi_d1phi_1->at(qq),Xi_d1mass_1->at(qq));
//pion
TLorentzVector lvector2;
lvector2.SetPtEtaPhiM(Xi_d1pt_2->at(qq),Xi_d1eta_2->at(qq),Xi_d1phi_2->at(qq),Xi_d1mass_2->at(qq));

pT_Xi->Fill(Xi_pt->at(qq),Xi_mass->at(qq));

if(Xi_id->at(qq) < 0){XXi_mass->Fill(Xi_mass->at(qq));}else{AXXi_mass->Fill(Xi_mass->at(qq));}

if(Xi_mass->at(qq) < massXAX - nsigmapeak*sigmaXAX || Xi_mass->at(qq) > massXAX + nsigmapeak*sigmaXAX)continue;

pT_Lamb_from_Xi_TF->Fill(Xi_pt->at(qq));

if(Xi_id->at(qq) < 0){

GoodTrackFourVector_AXi.push_back(pvector);
GoodTrackFourVector_AXi_pion_chi2.push_back(Xi_chi22->at(qq));
GoodTrackFourVector_AXi_bool.push_back(kTRUE);

GoodTrackFourVector_AXi_Lamb_bool_K.push_back(kTRUE);
GoodTrackFourVector_AXi_Lamb_bool_L.push_back(kTRUE);
GoodTrackFourVector_AXi_Lamb_bool_AL.push_back(kTRUE);
GoodTrackFourVector_AXi_Lamb.push_back(pvector1);
GoodTrackFourVector_AXi_Lamb_proton.push_back(lvector1);
GoodTrackFourVector_AXi_Lamb_pion.push_back(lvector2);
GoodTrackFourVector_AXi_Lamb_proton_chi2.push_back(Xi_chi21_1->at(qq));
GoodTrackFourVector_AXi_Lamb_pion_chi2.push_back(Xi_chi21_2->at(qq));
}else if(Xi_id->at(qq) > 0){

GoodTrackFourVector_Xi.push_back(pvector);
GoodTrackFourVector_Xi_pion_chi2.push_back(Xi_chi22->at(qq));
GoodTrackFourVector_Xi_bool.push_back(kTRUE);

GoodTrackFourVector_Xi_Lamb_bool_K.push_back(kTRUE);
GoodTrackFourVector_Xi_Lamb_bool_L.push_back(kTRUE);
GoodTrackFourVector_Xi_Lamb_bool_AL.push_back(kTRUE);
GoodTrackFourVector_Xi_Lamb.push_back(pvector1);
GoodTrackFourVector_Xi_Lamb_proton.push_back(lvector1);
GoodTrackFourVector_Xi_Lamb_pion.push_back(lvector2);
GoodTrackFourVector_Xi_Lamb_proton_chi2.push_back(Xi_chi21_1->at(qq));
GoodTrackFourVector_Xi_Lamb_pion_chi2.push_back(Xi_chi21_2->at(qq));
}
}
}

std::vector<TLorentzVector> GoodTrackFourVector_Om;
std::vector<TLorentzVector> GoodTrackFourVector_AOm;
std::vector<double> GoodTrackFourVector_Om_kaon_chi2;
std::vector<double> GoodTrackFourVector_AOm_kaon_chi2;
std::vector<Bool_t> GoodTrackFourVector_Om_bool;
std::vector<Bool_t> GoodTrackFourVector_AOm_bool;

std::vector<TLorentzVector> GoodTrackFourVector_Om_Lamb;
std::vector<TLorentzVector> GoodTrackFourVector_Om_Lamb_proton;
std::vector<TLorentzVector> GoodTrackFourVector_Om_Lamb_kaon;

std::vector<TLorentzVector> GoodTrackFourVector_AOm_Lamb;
std::vector<TLorentzVector> GoodTrackFourVector_AOm_Lamb_proton;
std::vector<TLorentzVector> GoodTrackFourVector_AOm_Lamb_kaon;

std::vector<double> GoodTrackFourVector_Om_Lamb_proton_chi2;
std::vector<double> GoodTrackFourVector_Om_Lamb_kaon_chi2;
std::vector<double> GoodTrackFourVector_AOm_Lamb_proton_chi2;
std::vector<double> GoodTrackFourVector_AOm_Lamb_kaon_chi2;

std::vector<Bool_t> GoodTrackFourVector_Om_Lamb_bool_K;
std::vector<Bool_t> GoodTrackFourVector_Om_Lamb_bool_L;
std::vector<Bool_t> GoodTrackFourVector_Om_Lamb_bool_AL;
std::vector<Bool_t> GoodTrackFourVector_AOm_Lamb_bool_K;
std::vector<Bool_t> GoodTrackFourVector_AOm_Lamb_bool_L;
std::vector<Bool_t> GoodTrackFourVector_AOm_Lamb_bool_AL;


if(DoFeedDown_Om){

//Omega + anti Omega
for(Int_t qw=0; qw<Omsize; qw++){ 
//cascade cuts
if(Om_pt->at(qw) <= ptminOAO || Om_pt->at(qw) >= ptmaxOAO)continue;
if(fabs(Om_eta->at(qw)) > 2.4)continue;
if(Om_VTrkP3DIpSigValue->at(qw) <= Om_VTrkP3DIpSigValue_val)continue;
if(Om_VTrkPi3DIpSigValue->at(qw) <= Om_VTrkPi3DIpSigValue_val)continue;
if(Om_casPi3DIpSigValue->at(qw) <= Om_casPi3DIpSigValue_val)continue;
if(Om_cas3DIpSigValue->at(qw) >= Om_cas3DIpSigValue_val)continue;
if(Om_distanceSigValue->at(qw) <= Om_distanceSigValue_val)continue;
if(Om_casFlightSigValue->at(qw) <= Om_casFlightSigValue_val)continue;

//daughter 2: pion
if(Om_d2Pix->at(qw) < pixelhits)continue;
//daughter 1: Lambda
//  -> proton cuts
if(Om_d1Pix_1->at(qw) < pixelhits)continue;
//  -> pion cuts 
if(Om_d1Pix_2->at(qw) < pixelhits)continue;

TLorentzVector pvector1;
pvector1.SetPtEtaPhiM(Om_d1pt->at(qw),Om_d1eta->at(qw),Om_d1phi->at(qw),Om_d1mass->at(qw));
TLorentzVector pvector2;
pvector2.SetPtEtaPhiM(Om_d2pt->at(qw),Om_d2eta->at(qw),Om_d2phi->at(qw),Om_d2mass->at(qw));

double Pxp = pvector1.Px();
double Pyp = pvector1.Py();
double Pzp = pvector1.Pz();
double Pmp = pvector1.M();

double Pxn = pvector2.Px();
double Pyn = pvector2.Py();
double Pzn = pvector2.Pz();
double Pmn = pvector2.M();


TVector3 dauvec1(pvector1.Px(),pvector1.Py(),pvector1.Pz());
TVector3 dauvec2(pvector2.Px(),pvector2.Py(),pvector2.Pz());
TVector3 dauvecsum(dauvec1+dauvec2);

// Armenteros-Podolanski
 double Pp=sqrt(Pxp*Pxp+Pyp*Pyp+Pzp*Pzp);
 double Pn=sqrt(Pxn*Pxn+Pyn*Pyn+Pzn*Pzn);
 double PN=sqrt((Pxp+Pxn)*(Pxp+Pxn)+(Pyp+Pyn)*(Pyp+Pyn)+(Pzp+Pzn)*(Pzp+Pzn));
 double cosx1=(Pxp*(Pxp+Pxn)+Pyp*(Pyp+Pyn)+Pzp*(Pzp+Pzn))/(Pp*PN);
 double cosx2=(Pxn*(Pxp+Pxn)+Pyn*(Pyp+Pyn)+Pzn*(Pzp+Pzn))/(Pn*PN);
 double QT = Pp*sqrt(1-cosx1*cosx1);
 double Alpha;
 if(Om_id->at(qw) < 0){Alpha = (-Pp*cosx1+Pn*cosx2)/(Pp*cosx1+Pn*cosx2);}
 if(Om_id->at(qw) > 0){Alpha = (Pp*cosx1-Pn*cosx2)/(Pp*cosx1+Pn*cosx2);}
 
double energyd1e = sqrt(electronMass*electronMass+pvector1.P()*pvector1.P());
double energyd2e = sqrt(electronMass*electronMass+pvector2.P()*pvector2.P());
double invmass_ee = sqrt((energyd1e+energyd2e)*(energyd1e+energyd2e)-dauvecsum.Mag2());

if(invmass_ee<misee)continue;// do not change anymore

double massd1, massd2, energyd1, energyd2, invmass;

massd1=pimass;
massd2=la_mass;
energyd1 = sqrt(massd1*massd1+pvector1.P()*pvector1.P());
energyd2 = sqrt(massd2*massd2+pvector2.P()*pvector2.P());
invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
if(fabs(invmass-xi_mass)<misOm) continue;//ok

massd1=la_mass;
massd2=pimass;
energyd1 = sqrt(massd1*massd1+pvector1.P()*pvector1.P());
energyd2 = sqrt(massd2*massd2+pvector2.P()*pvector2.P());
invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
if(fabs(invmass-xi_mass)<misOm) continue;// do not change anymore

TLorentzVector pvector;
pvector.SetPtEtaPhiM(Om_pt->at(qw),Om_eta->at(qw),Om_phi->at(qw),Om_mass->at(qw));
if(fabs(pvector.Rapidity()) >= rapvar)continue;

//AP for lambda daughter
TLorentzVector pvector1_l; //proton
pvector1_l.SetPtEtaPhiM(Om_d1pt_1->at(qw),Om_d1eta_1->at(qw),Om_d1phi_1->at(qw),Om_d1mass_1->at(qw));
TLorentzVector pvector2_l; //pion
pvector2_l.SetPtEtaPhiM(Om_d1pt_2->at(qw),Om_d1eta_2->at(qw),Om_d1phi_2->at(qw),Om_d1mass_2->at(qw));

double Pxp_l = pvector1_l.Px();
double Pyp_l = pvector1_l.Py();
double Pzp_l = pvector1_l.Pz();
double Pmp_l = pvector1_l.M();
double Pxn_l = pvector2_l.Px();
double Pyn_l = pvector2_l.Py();
double Pzn_l = pvector2_l.Pz();
double Pmn_l = pvector2_l.M();

TVector3 dauvec1_l(pvector1_l.Px(),pvector1_l.Py(),pvector1_l.Pz());
TVector3 dauvec2_l(pvector2_l.Px(),pvector2_l.Py(),pvector2_l.Pz());
TVector3 dauvecsum_l(dauvec1_l+dauvec2_l);
 
double energyd1e_l = sqrt(electronMass*electronMass+pvector1_l.P()*pvector1_l.P());
double energyd2e_l = sqrt(electronMass*electronMass+pvector2_l.P()*pvector2_l.P());
double invmass_ee_l = sqrt((energyd1e_l+energyd2e_l)*(energyd1e_l+energyd2e_l)-dauvecsum_l.Mag2());
if(invmass_ee_l<misee)continue;// do not change anymore

double massd1_l, massd2_l, energyd1_l, energyd2_l, invmass_l;
massd1_l=pimass;
massd2_l=pimass;
energyd1_l = sqrt(massd1_l*massd1_l+pvector1_l.P()*pvector1_l.P());
energyd2_l = sqrt(massd2_l*massd2_l+pvector2_l.P()*pvector2_l.P());
invmass_l = sqrt((energyd1_l+energyd2_l)*(energyd1_l+energyd2_l)-dauvecsum_l.Mag2());
if(fabs(invmass_l-k0s_mass)<misLam) continue;//ok

//lambda daughters
//proton
TLorentzVector lvector1;
lvector1.SetPtEtaPhiM(Om_d1pt_1->at(qw),Om_d1eta_1->at(qw),Om_d1phi_1->at(qw),Om_d1mass_1->at(qw));
//kaon
TLorentzVector lvector2;
lvector2.SetPtEtaPhiM(Om_d1pt_2->at(qw),Om_d1eta_2->at(qw),Om_d1phi_2->at(qw),Om_d1mass_2->at(qw));

pT_Om->Fill(Om_pt->at(qw),Om_mass->at(qw));
if(Om_id->at(qw) < 0){OOm_mass->Fill(Om_mass->at(qw));}else{AOOm_mass->Fill(Om_mass->at(qw));}


if(Om_mass->at(qw) < massOAO - nsigmapeak*sigmaOAO || Om_mass->at(qw) > massOAO + nsigmapeak*sigmaOAO)continue;

pT_Lamb_from_Om_TF->Fill(Om_pt->at(qw));

if(Om_id->at(qw) < 0){
 
GoodTrackFourVector_AOm.push_back(pvector);
GoodTrackFourVector_AOm_kaon_chi2.push_back(Om_chi22->at(qw));
GoodTrackFourVector_AOm_bool.push_back(kTRUE);

GoodTrackFourVector_AOm_Lamb_bool_K.push_back(kTRUE);
GoodTrackFourVector_AOm_Lamb_bool_L.push_back(kTRUE);
GoodTrackFourVector_AOm_Lamb_bool_AL.push_back(kTRUE);
GoodTrackFourVector_AOm_Lamb.push_back(pvector1); //lambda
GoodTrackFourVector_AOm_Lamb_proton.push_back(lvector1);
GoodTrackFourVector_AOm_Lamb_kaon.push_back(lvector2);
GoodTrackFourVector_AOm_Lamb_proton_chi2.push_back(Om_chi21_1->at(qw));
GoodTrackFourVector_AOm_Lamb_kaon_chi2.push_back(Om_chi21_2->at(qw));
}else if(Om_id->at(qw) > 0){

GoodTrackFourVector_Om.push_back(pvector);
GoodTrackFourVector_Om_kaon_chi2.push_back(Om_chi22->at(qw));
GoodTrackFourVector_Om_bool.push_back(kTRUE);

GoodTrackFourVector_Om_Lamb_bool_K.push_back(kTRUE);
GoodTrackFourVector_Om_Lamb_bool_L.push_back(kTRUE);
GoodTrackFourVector_Om_Lamb_bool_AL.push_back(kTRUE);
GoodTrackFourVector_Om_Lamb.push_back(pvector1); //lambda
GoodTrackFourVector_Om_Lamb_proton.push_back(lvector1);
GoodTrackFourVector_Om_Lamb_kaon.push_back(lvector2);
GoodTrackFourVector_Om_Lamb_proton_chi2.push_back(Om_chi21_1->at(qw));
GoodTrackFourVector_Om_Lamb_kaon_chi2.push_back(Om_chi21_2->at(qw));
}

}

}


//Remove particles sharing the same daughter

if(GoodTrackFourVector_K0s.size()>=2){
for(Int_t ik=0; ik<GoodTrackFourVector_K0s.size(); ik++){
double chi2V0 = chi2_K0s[ik];
double chi2d1 = chi21_K0s[ik];
double chi2d2 = chi22_K0s[ik];
double mass1_dif = fabs(GoodTrackFourVector_K0s[ik].M() - k0s_mass);
for(Int_t iik=ik+1; iik<GoodTrackFourVector_K0s.size(); iik++){
double chi2V0a = chi2_K0s[iik];
double chi2d1a = chi21_K0s[iik];
double chi2d2a = chi22_K0s[iik];
double mass2_dif = fabs(GoodTrackFourVector_K0s[iik].M() - k0s_mass);
if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){
if(method==0){map_K0s_chi2[ik] = kFALSE;map_K0s_chi2[iik] = kFALSE;}
if(method==1){if(i % 2 == 0){map_K0s_chi2[ik] = kFALSE;}else{map_K0s_chi2[iik] = kFALSE;}}
if(method==2){if(mass1_dif >= mass2_dif){map_K0s_chi2[ik] = kFALSE;}else{map_K0s_chi2[iik] = kFALSE;}}
if(method==3){if(chi2V0 >= chi2V0a){map_K0s_chi2[ik] = kFALSE;}else{map_K0s_chi2[iik] = kFALSE;}}
if(method>3){if ( i==0 )cout << "do not remove mothers" << endl;}
}}}}

if(GoodTrackFourVector_Lam.size()>=2){
for(Int_t ik=0; ik<GoodTrackFourVector_Lam.size(); ik++){
double chi2V0 = chi2_Lam[ik];
double chi2d1 = chi21_Lam[ik];
double chi2d2 = chi22_Lam[ik];
double mass1_dif = fabs(GoodTrackFourVector_Lam[ik].M() - la_mass);
for(Int_t iik=ik+1; iik<GoodTrackFourVector_Lam.size(); iik++){
double chi2V0a = chi2_Lam[iik];
double chi2d1a = chi21_Lam[iik];
double chi2d2a = chi22_Lam[iik];
double mass2_dif = fabs(GoodTrackFourVector_Lam[iik].M() - la_mass);
if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){
if(method==0){map_Lam_chi2[ik] = kFALSE;map_Lam_chi2[iik] = kFALSE;}
if(method==1){if(i % 2 == 0){map_Lam_chi2[ik] = kFALSE;}else{map_Lam_chi2[iik] = kFALSE;}}
if(method==2){if(mass1_dif >= mass2_dif){map_Lam_chi2[ik] = kFALSE;}else{map_Lam_chi2[iik] = kFALSE;}}
if(method==3){if(chi2V0 >= chi2V0a){map_Lam_chi2[ik] = kFALSE;}else{map_Lam_chi2[iik] = kFALSE;}}
if(method>3){if ( i==0 )cout << "do not remove mothers" << endl;}
}}}}

if(GoodTrackFourVector_ALam.size()>=2){
for(Int_t ik=0; ik<GoodTrackFourVector_ALam.size(); ik++){
double chi2V0 = chi2_ALam[ik];
double chi2d1 = chi21_ALam[ik];
double chi2d2 = chi22_ALam[ik];
double mass1_dif = fabs(GoodTrackFourVector_ALam[ik].M() - la_mass);
for(Int_t iik=ik+1; iik<GoodTrackFourVector_ALam.size(); iik++){
double chi2V0a = chi2_ALam[iik];
double chi2d1a = chi21_ALam[iik];
double chi2d2a = chi22_ALam[iik];
double mass2_dif = fabs(GoodTrackFourVector_ALam[iik].M() - la_mass);
if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){
if(method==0){map_ALam_chi2[ik] = kFALSE;map_ALam_chi2[iik] = kFALSE;}
if(method==1){if(i % 2 == 0){map_ALam_chi2[ik] = kFALSE;}else{map_ALam_chi2[iik] = kFALSE;}}
if(method==2){if(mass1_dif >= mass2_dif){map_ALam_chi2[ik] = kFALSE;}else{map_ALam_chi2[iik] = kFALSE;}}
if(method==3){if(chi2V0 >= chi2V0a){map_ALam_chi2[ik] = kFALSE;}else{map_ALam_chi2[iik] = kFALSE;}}
if(method>3){if ( i==0 )cout << "do not remove mothers" << endl;}
}}}}

if(GoodTrackFourVector_Lam.size()>=1 && GoodTrackFourVector_ALam.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_Lam.size(); ik++){
double chi2V0 = chi2_Lam[ik];
double chi2d1 = chi21_Lam[ik];
double chi2d2 = chi22_Lam[ik];
double mass1_dif = fabs(GoodTrackFourVector_Lam[ik].M() - la_mass);
for(Int_t iik=0; iik<GoodTrackFourVector_ALam.size(); iik++){
double chi2V0a = chi2_ALam[iik];
double chi2d1a = chi21_ALam[iik];
double chi2d2a = chi22_ALam[iik];
double mass2_dif = fabs(GoodTrackFourVector_ALam[iik].M() - la_mass);
if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){
if(method==0){map_Lam_chi2_LAL[ik] = kFALSE;map_ALam_chi2_LAL[iik] = kFALSE;}
if(method==1){if(i % 2 == 0){map_Lam_chi2_LAL[ik] = kFALSE;}else{map_ALam_chi2_LAL[iik] = kFALSE;}}
if(method==2){if(mass1_dif >= mass2_dif){map_Lam_chi2_LAL[ik] = kFALSE;}else{map_ALam_chi2_LAL[iik] = kFALSE;}}
if(method==3){if(chi2V0 >= chi2V0a){map_Lam_chi2_LAL[ik] = kFALSE;}else{map_ALam_chi2_LAL[iik] = kFALSE;}}
if(method>3){if ( i==0 )cout << "do not remove mothers" << endl;}
}}}}

if(GoodTrackFourVector_K0s.size()>=1 && GoodTrackFourVector_ALam.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_K0s.size(); ik++){
double chi2V0 = chi2_K0s[ik];
double chi2d1 = chi21_K0s[ik];
double chi2d2 = chi22_K0s[ik];
double mass1_dif = fabs(GoodTrackFourVector_K0s[ik].M() - k0s_mass);
for(Int_t iik=0; iik<GoodTrackFourVector_ALam.size(); iik++){
double chi2V0a = chi2_ALam[iik];
double chi2d1a = chi21_ALam[iik];
double chi2d2a = chi22_ALam[iik];
double mass2_dif = fabs(GoodTrackFourVector_ALam[iik].M() - la_mass);
if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){
if(method==0){map_K0s_chi2_KAL[ik] = kFALSE;map_ALam_chi2_KAL[iik] = kFALSE;}
if(method==1){if(i % 2 == 0){map_K0s_chi2_KAL[ik] = kFALSE;}else{map_ALam_chi2_KAL[iik] = kFALSE;}}
if(method==2){if(mass1_dif >= mass2_dif){map_K0s_chi2_KAL[ik] = kFALSE;}else{map_ALam_chi2_KAL[iik] = kFALSE;}}
if(method==3){if(chi2V0 >= chi2V0a){map_K0s_chi2_KAL[ik] = kFALSE;}else{map_ALam_chi2_KAL[iik] = kFALSE;}}
if(method>3){if ( i==0 )cout << "do not remove mothers" << endl;}
}}}}

if(GoodTrackFourVector_K0s.size()>=1 && GoodTrackFourVector_Lam.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_K0s.size(); ik++){
double chi2V0 = chi2_K0s[ik];
double chi2d1 = chi21_K0s[ik];
double chi2d2 = chi22_K0s[ik];
double mass1_dif = fabs(GoodTrackFourVector_K0s[ik].M() - k0s_mass);
for(Int_t iik=0; iik<GoodTrackFourVector_Lam.size(); iik++){
double chi2V0a = chi2_Lam[iik];
double chi2d1a = chi21_Lam[iik];
double chi2d2a = chi22_Lam[iik];
double mass2_dif = fabs(GoodTrackFourVector_Lam[iik].M() - la_mass);
if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){
if(method==0){map_K0s_chi2_KL[ik] = kFALSE;map_Lam_chi2_KL[iik] = kFALSE;}
if(method==1){if(i % 2 == 0){map_K0s_chi2_KL[ik] = kFALSE;}else{map_Lam_chi2_KL[iik] = kFALSE;}}
if(method==2){if(mass1_dif >= mass2_dif){map_K0s_chi2_KL[ik] = kFALSE;}else{map_Lam_chi2_KL[iik] = kFALSE;}}
if(method==3){if(chi2V0 >= chi2V0a){map_K0s_chi2_KL[ik] = kFALSE;}else{map_Lam_chi2_KL[iik] = kFALSE;}}
if(method>3){if ( i==0 )cout << "do not remove mothers" << endl;}
}}}}

if(DoFeedDown_Xi){

if(GoodTrackFourVector_K0s.size()>=1 && GoodTrackFourVector_Xi_Lamb.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_K0s.size(); ik++){
double chi2d1 = chi21_K0s[ik];
double chi2d2 = chi22_K0s[ik];
for(Int_t iik=0; iik<GoodTrackFourVector_Xi_Lamb.size(); iik++){
double chi2d1a = GoodTrackFourVector_Xi_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_Xi_Lamb_pion_chi2[iik];
if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){
GoodTrackFourVector_Xi_Lamb_bool_K[iik] = kFALSE;
}}}}

if(GoodTrackFourVector_K0s.size()>=1 && GoodTrackFourVector_AXi_Lamb.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_K0s.size(); ik++){
double chi2d1 = chi21_K0s[ik];
double chi2d2 = chi22_K0s[ik];
for(Int_t iik=0; iik<GoodTrackFourVector_AXi_Lamb.size(); iik++){
double chi2d1a = GoodTrackFourVector_AXi_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_AXi_Lamb_pion_chi2[iik];
if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){
GoodTrackFourVector_AXi_Lamb_bool_K[iik] = kFALSE;
}}}}

if(GoodTrackFourVector_Lam.size()>=1 && GoodTrackFourVector_Xi_Lamb.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_Lam.size(); ik++){
double chi2d1 = chi21_Lam[ik];
double chi2d2 = chi22_Lam[ik];
double pt_lam = GoodTrackFourVector_Lam[ik].Pt();
double m_lam = GoodTrackFourVector_Lam[ik].M();
double phi_lam = GoodTrackFourVector_Lam[ik].Phi();
double y_lam = GoodTrackFourVector_Lam[ik].Rapidity();

for(Int_t iik=0; iik<GoodTrackFourVector_Xi_Lamb.size(); iik++){
double chi2d1a = GoodTrackFourVector_Xi_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_Xi_Lamb_pion_chi2[iik];
double pt_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].Rapidity();

if(pt_lam == pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed){continue;}

if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){GoodTrackFourVector_Xi_Lamb_bool_L[iik] = kFALSE;}
}}}

if(GoodTrackFourVector_Lam.size()>=1 && GoodTrackFourVector_AXi_Lamb.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_Lam.size(); ik++){
double chi2d1 = chi21_Lam[ik];
double chi2d2 = chi22_Lam[ik];
double pt_lam = GoodTrackFourVector_Lam[ik].Pt();
double m_lam = GoodTrackFourVector_Lam[ik].M();
double phi_lam = GoodTrackFourVector_Lam[ik].Phi();
double y_lam = GoodTrackFourVector_Lam[ik].Rapidity();

for(Int_t iik=0; iik<GoodTrackFourVector_AXi_Lamb.size(); iik++){
double chi2d1a = GoodTrackFourVector_AXi_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_AXi_Lamb_pion_chi2[iik];
double pt_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].Rapidity();

if(pt_lam == pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed){continue;}

if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){GoodTrackFourVector_AXi_Lamb_bool_L[iik] = kFALSE;}
}}}

if(GoodTrackFourVector_ALam.size()>=1 && GoodTrackFourVector_Xi_Lamb.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_ALam.size(); ik++){
double chi2d1 = chi21_ALam[ik];
double chi2d2 = chi22_ALam[ik];
double pt_lam = GoodTrackFourVector_ALam[ik].Pt();
double m_lam = GoodTrackFourVector_ALam[ik].M();
double phi_lam = GoodTrackFourVector_ALam[ik].Phi();
double y_lam = GoodTrackFourVector_ALam[ik].Rapidity();

for(Int_t iik=0; iik<GoodTrackFourVector_Xi_Lamb.size(); iik++){
double chi2d1a = GoodTrackFourVector_Xi_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_Xi_Lamb_pion_chi2[iik];
double pt_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].Rapidity();

if(pt_lam == pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed){continue;}

if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){GoodTrackFourVector_Xi_Lamb_bool_AL[iik] = kFALSE;}
}}}

if(GoodTrackFourVector_ALam.size()>=1 && GoodTrackFourVector_AXi_Lamb.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_ALam.size(); ik++){
double chi2d1 = chi21_ALam[ik];
double chi2d2 = chi22_ALam[ik];
double pt_lam = GoodTrackFourVector_ALam[ik].Pt();
double m_lam = GoodTrackFourVector_ALam[ik].M();
double phi_lam = GoodTrackFourVector_ALam[ik].Phi();
double y_lam = GoodTrackFourVector_ALam[ik].Rapidity();

for(Int_t iik=0; iik<GoodTrackFourVector_AXi_Lamb.size(); iik++){
double chi2d1a = GoodTrackFourVector_AXi_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_AXi_Lamb_pion_chi2[iik];
double pt_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].Rapidity();

if(pt_lam == pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed){continue;}

if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){GoodTrackFourVector_AXi_Lamb_bool_AL[iik] = kFALSE;}
}}}

if(GoodTrackFourVector_Xi.size()>=2){
for(Int_t ik=0; ik<GoodTrackFourVector_Xi.size(); ik++){
double chi2d1 = GoodTrackFourVector_Xi_Lamb_proton_chi2[ik];
double chi2d2 = GoodTrackFourVector_Xi_Lamb_pion_chi2[ik];
double chi2dX = GoodTrackFourVector_Xi_pion_chi2[ik];
for(Int_t iik=ik+1; iik<GoodTrackFourVector_Xi.size(); iik++){
double chi2d1a = GoodTrackFourVector_Xi_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_Xi_Lamb_pion_chi2[iik];
double chi2dXa = GoodTrackFourVector_Xi_pion_chi2[iik];
if((chi2d1 == chi2d1a) || (chi2d1 == chi2d2a) || (chi2d1 == chi2dXa) || (chi2d2 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d2 == chi2dXa) || (chi2dX == chi2d1a) || (chi2dX == chi2d2a) || (chi2dX == chi2dXa) ){GoodTrackFourVector_Xi_bool[ik] = kFALSE; GoodTrackFourVector_Xi_bool[iik] = kFALSE;}
}}}

if(GoodTrackFourVector_AXi.size()>=2){
for(Int_t ik=0; ik<GoodTrackFourVector_AXi.size(); ik++){
double chi2d1 = GoodTrackFourVector_AXi_Lamb_proton_chi2[ik];
double chi2d2 = GoodTrackFourVector_AXi_Lamb_pion_chi2[ik];
double chi2dX = GoodTrackFourVector_AXi_pion_chi2[ik];
for(Int_t iik=ik+1; iik<GoodTrackFourVector_AXi.size(); iik++){
double chi2d1a = GoodTrackFourVector_AXi_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_AXi_Lamb_pion_chi2[iik];
double chi2dXa = GoodTrackFourVector_AXi_pion_chi2[iik];
if((chi2d1 == chi2d1a) || (chi2d1 == chi2d2a) || (chi2d1 == chi2dXa) || (chi2d2 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d2 == chi2dXa) || (chi2dX == chi2d1a) || (chi2dX == chi2d2a) || (chi2dX == chi2dXa) ){GoodTrackFourVector_AXi_bool[ik] = kFALSE; GoodTrackFourVector_AXi_bool[iik] = kFALSE;}
}}}

if(GoodTrackFourVector_Xi.size()>=1 && GoodTrackFourVector_AXi.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_Xi.size(); ik++){
double chi2d1 = GoodTrackFourVector_Xi_Lamb_proton_chi2[ik];
double chi2d2 = GoodTrackFourVector_Xi_Lamb_pion_chi2[ik];
double chi2dX = GoodTrackFourVector_Xi_pion_chi2[ik];
for(Int_t iik=0; iik<GoodTrackFourVector_AXi.size(); iik++){
double chi2d1a = GoodTrackFourVector_AXi_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_AXi_Lamb_pion_chi2[iik];
double chi2dXa = GoodTrackFourVector_AXi_pion_chi2[iik];
if((chi2d1 == chi2d1a) || (chi2d1 == chi2d2a) || (chi2d1 == chi2dXa) || (chi2d2 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d2 == chi2dXa) || (chi2dX == chi2d1a) || (chi2dX == chi2d2a) || (chi2dX == chi2dXa) ){GoodTrackFourVector_Xi_bool[ik] = kFALSE; GoodTrackFourVector_AXi_bool[iik] = kFALSE;}
}}}


}

if(DoFeedDown_Om){

if(GoodTrackFourVector_K0s.size()>=1 && GoodTrackFourVector_Om_Lamb.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_K0s.size(); ik++){
double chi2d1 = chi21_K0s[ik];
double chi2d2 = chi22_K0s[ik];
for(Int_t iik=0; iik<GoodTrackFourVector_Om_Lamb.size(); iik++){
double chi2d1a = GoodTrackFourVector_Om_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_Om_Lamb_kaon_chi2[iik];
if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){
GoodTrackFourVector_Om_Lamb_bool_K[iik] = kFALSE;
}}}}

if(GoodTrackFourVector_K0s.size()>=1 && GoodTrackFourVector_AOm_Lamb.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_K0s.size(); ik++){
double chi2d1 = chi21_K0s[ik];
double chi2d2 = chi22_K0s[ik];
for(Int_t iik=0; iik<GoodTrackFourVector_AOm_Lamb.size(); iik++){
double chi2d1a = GoodTrackFourVector_AOm_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_AOm_Lamb_kaon_chi2[iik];
if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){
GoodTrackFourVector_AOm_Lamb_bool_K[iik] = kFALSE;
}}}}

if(GoodTrackFourVector_Lam.size()>=1 && GoodTrackFourVector_Om_Lamb.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_Lam.size(); ik++){
double chi2d1 = chi21_Lam[ik];
double chi2d2 = chi22_Lam[ik];
double pt_lam = GoodTrackFourVector_Lam[ik].Pt();
double m_lam = GoodTrackFourVector_Lam[ik].M();
double phi_lam = GoodTrackFourVector_Lam[ik].Phi();
double y_lam = GoodTrackFourVector_Lam[ik].Rapidity();

for(Int_t iik=0; iik<GoodTrackFourVector_Om_Lamb.size(); iik++){
double chi2d1a = GoodTrackFourVector_Om_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_Om_Lamb_kaon_chi2[iik];
double pt_lam_feed = GoodTrackFourVector_Om_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_Om_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_Om_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_Om_Lamb[iik].Rapidity();

if(pt_lam == pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed){continue;}

if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){GoodTrackFourVector_Om_Lamb_bool_L[iik] = kFALSE;}
}}}

if(GoodTrackFourVector_Lam.size()>=1 && GoodTrackFourVector_AOm_Lamb.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_Lam.size(); ik++){
double chi2d1 = chi21_Lam[ik];
double chi2d2 = chi22_Lam[ik];
double pt_lam = GoodTrackFourVector_Lam[ik].Pt();
double m_lam = GoodTrackFourVector_Lam[ik].M();
double phi_lam = GoodTrackFourVector_Lam[ik].Phi();
double y_lam = GoodTrackFourVector_Lam[ik].Rapidity();

for(Int_t iik=0; iik<GoodTrackFourVector_AOm_Lamb.size(); iik++){
double chi2d1a = GoodTrackFourVector_AOm_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_AOm_Lamb_kaon_chi2[iik];
double pt_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].Rapidity();

if(pt_lam == pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed){continue;}

if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){GoodTrackFourVector_AOm_Lamb_bool_L[iik] = kFALSE;}
}}}

if(GoodTrackFourVector_ALam.size()>=1 && GoodTrackFourVector_Om_Lamb.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_ALam.size(); ik++){
double chi2d1 = chi21_ALam[ik];
double chi2d2 = chi22_ALam[ik];
double pt_lam = GoodTrackFourVector_ALam[ik].Pt();
double m_lam = GoodTrackFourVector_ALam[ik].M();
double phi_lam = GoodTrackFourVector_ALam[ik].Phi();
double y_lam = GoodTrackFourVector_ALam[ik].Rapidity();

for(Int_t iik=0; iik<GoodTrackFourVector_Om_Lamb.size(); iik++){
double chi2d1a = GoodTrackFourVector_Om_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_Om_Lamb_kaon_chi2[iik];
double pt_lam_feed = GoodTrackFourVector_Om_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_Om_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_Om_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_Om_Lamb[iik].Rapidity();

if(pt_lam == pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed){continue;}

if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){GoodTrackFourVector_Om_Lamb_bool_AL[iik] = kFALSE;}
}}}

if(GoodTrackFourVector_ALam.size()>=1 && GoodTrackFourVector_AOm_Lamb.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_ALam.size(); ik++){
double chi2d1 = chi21_ALam[ik];
double chi2d2 = chi22_ALam[ik];
double pt_lam = GoodTrackFourVector_ALam[ik].Pt();
double m_lam = GoodTrackFourVector_ALam[ik].M();
double phi_lam = GoodTrackFourVector_ALam[ik].Phi();
double y_lam = GoodTrackFourVector_ALam[ik].Rapidity();

for(Int_t iik=0; iik<GoodTrackFourVector_AOm_Lamb.size(); iik++){
double chi2d1a = GoodTrackFourVector_AOm_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_AOm_Lamb_kaon_chi2[iik];
double pt_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].Rapidity();

if(pt_lam == pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed){continue;}

if((chi2d1 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d1 == chi2d2a) || (chi2d2 == chi2d1a)){GoodTrackFourVector_AOm_Lamb_bool_AL[iik] = kFALSE;}
}}}

if(GoodTrackFourVector_Om.size()>=2){
for(Int_t ik=0; ik<GoodTrackFourVector_Om.size(); ik++){
double chi2d1 = GoodTrackFourVector_Om_Lamb_proton_chi2[ik];
double chi2d2 = GoodTrackFourVector_Om_Lamb_kaon_chi2[ik];
double chi2dX = GoodTrackFourVector_Om_kaon_chi2[ik];
for(Int_t iik=ik+1; iik<GoodTrackFourVector_Om.size(); iik++){
double chi2d1a = GoodTrackFourVector_Om_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_Om_Lamb_kaon_chi2[iik];
double chi2dXa = GoodTrackFourVector_Om_kaon_chi2[iik];
if((chi2d1 == chi2d1a) || (chi2d1 == chi2d2a) || (chi2d1 == chi2dXa) || (chi2d2 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d2 == chi2dXa) || (chi2dX == chi2d1a) || (chi2dX == chi2d2a) || (chi2dX == chi2dXa) ){GoodTrackFourVector_Om_bool[ik] = kFALSE;GoodTrackFourVector_Om_bool[iik] = kFALSE;}
}}}

if(GoodTrackFourVector_AOm.size()>=2){
for(Int_t ik=0; ik<GoodTrackFourVector_AOm.size(); ik++){
double chi2d1 = GoodTrackFourVector_AOm_Lamb_proton_chi2[ik];
double chi2d2 = GoodTrackFourVector_AOm_Lamb_kaon_chi2[ik];
double chi2dX = GoodTrackFourVector_AOm_kaon_chi2[ik];
for(Int_t iik=ik+1; iik<GoodTrackFourVector_AOm.size(); iik++){
double chi2d1a = GoodTrackFourVector_AOm_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_AOm_Lamb_kaon_chi2[iik];
double chi2dXa = GoodTrackFourVector_AOm_kaon_chi2[iik];
if((chi2d1 == chi2d1a) || (chi2d1 == chi2d2a) || (chi2d1 == chi2dXa) || (chi2d2 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d2 == chi2dXa) || (chi2dX == chi2d1a) || (chi2dX == chi2d2a) || (chi2dX == chi2dXa) ){GoodTrackFourVector_AOm_bool[ik] = kFALSE;GoodTrackFourVector_AOm_bool[iik] = kFALSE;}
}}}

if(GoodTrackFourVector_Om.size()>=1 && GoodTrackFourVector_AOm.size()>=1){
for(Int_t ik=0; ik<GoodTrackFourVector_Om.size(); ik++){
double chi2d1 = GoodTrackFourVector_Om_Lamb_proton_chi2[ik];
double chi2d2 = GoodTrackFourVector_Om_Lamb_kaon_chi2[ik];
double chi2dX = GoodTrackFourVector_Om_kaon_chi2[ik];
for(Int_t iik=0; iik<GoodTrackFourVector_AOm.size(); iik++){
double chi2d1a = GoodTrackFourVector_AOm_Lamb_proton_chi2[iik];
double chi2d2a = GoodTrackFourVector_AOm_Lamb_kaon_chi2[iik];
double chi2dXa = GoodTrackFourVector_AOm_kaon_chi2[iik];
if((chi2d1 == chi2d1a) || (chi2d1 == chi2d2a) || (chi2d1 == chi2dXa) || (chi2d2 == chi2d1a) || (chi2d2 == chi2d2a) || (chi2d2 == chi2dXa) || (chi2dX == chi2d1a) || (chi2dX == chi2d2a) || (chi2dX == chi2dXa) ){GoodTrackFourVector_Om_bool[ik] = kFALSE; GoodTrackFourVector_AOm_bool[iik] = kFALSE;}
}}}

}

if(isMC){

const int a1 = GoodTrackFourVector_K0s.size();
Int_t NmatchK0s_Gen[a1];
const int a2 = GoodTrackFourVector_Lam.size();
Int_t NmatchLam_Gen[a2];
const int a3 = GoodTrackFourVector_ALam.size();
Int_t NmatchALam_Gen[a3];

//matching
if(map_GenK0s.size() >= 1 && GoodTrackFourVector_K0s.size() >= 1){

for(Int_t ik1x=0; ik1x<GoodTrackFourVector_K0s.size(); ++ik1x){
	double      pTR = GoodTrackFourVector_K0s[ik1x].Pt();
	double      etaR = GoodTrackFourVector_K0s[ik1x].Eta();
	double      phiR = GoodTrackFourVector_K0s[ik1x].Phi();
    NmatchK0s_Gen[ik1x] = -999;
for(Int_t ik2x=0; ik2x<map_GenK0s.size(); ++ik2x){
	double pTG = map_GenK0s[ik2x].Pt();
	double etaG = map_GenK0s[ik2x].Eta();
	double phiG = map_GenK0s[ik2x].Phi();
	double dphi = GoodTrackFourVector_K0s[ik1x].DeltaPhi(map_GenK0s[ik2x]);//phiR - phiG;
    double deta = etaR - etaG;
    double dR = sqrt(dphi*dphi + deta*deta);
    double dpt = (pTG - pTR)/(pTG);
	DR_K0s->Fill(dR);
	DPT_K0s->Fill(dpt);
    if(dR < drdpt && fabs(dpt) < drdpt){
        NmatchK0s_Gen[ik1x] = ik2x;
        break;
}else{}}}

for(Int_t ik1x=0; ik1x<GoodTrackFourVector_K0s.size(); ++ik1x){
	double      pTR = GoodTrackFourVector_K0s[ik1x].Pt();
	double      etaR = GoodTrackFourVector_K0s[ik1x].Eta();
	double      phiR = GoodTrackFourVector_K0s[ik1x].Phi();
	bool samegen = false;
	for(Int_t ik1xx=ik1x+1; ik1xx<GoodTrackFourVector_K0s.size(); ++ik1xx){if(NmatchK0s_Gen[ik1x] >= 0 && NmatchK0s_Gen[ik1xx] >= 0 && NmatchK0s_Gen[ik1x] == NmatchK0s_Gen[ik1xx]){samegen=true;}}
	if(samegen)continue;
for(Int_t ik2x=0; ik2x<map_GenK0s.size(); ++ik2x){
	double pTG = map_GenK0s[ik2x].Pt();
	double etaG = map_GenK0s[ik2x].Eta();
	double phiG = map_GenK0s[ik2x].Phi();
	double dphi = GoodTrackFourVector_K0s[ik1x].DeltaPhi(map_GenK0s[ik2x]);//phiR - phiG;
    double deta = etaR - etaG;
    double dR = sqrt(dphi*dphi + deta*deta);
    double dpt = (pTG - pTR)/(pTG);
    if(dR < drdpt && fabs(dpt) < drdpt){
    map_K0s_match[ik1x] = kTRUE;
    if(GoodTrackFourVector_K0s[ik1x].M() > massK0s - nsigmapeak*sigmaK0s && GoodTrackFourVector_K0s[ik1x].M() < massK0s + nsigmapeak*sigmaK0s){map_GenK0s_matched.push_back(map_GenK0s[ik2x]);}
    break;
}else{}}}

}

if(map_GenLam.size() >= 1 && GoodTrackFourVector_Lam.size() >= 1){

for(Int_t ik1x=0; ik1x<GoodTrackFourVector_Lam.size(); ++ik1x){
	double      pTR = GoodTrackFourVector_Lam[ik1x].Pt();
	double      etaR = GoodTrackFourVector_Lam[ik1x].Eta();
	double      phiR = GoodTrackFourVector_Lam[ik1x].Phi();
    NmatchLam_Gen[ik1x] = -999;
for(Int_t ik2x=0; ik2x<map_GenLam.size(); ++ik2x){
	double pTG = map_GenLam[ik2x].Pt();
	double etaG = map_GenLam[ik2x].Eta();
	double phiG = map_GenLam[ik2x].Phi();
	double dphi = GoodTrackFourVector_Lam[ik1x].DeltaPhi(map_GenLam[ik2x]);//phiR - phiG;
    double deta = etaR - etaG;
    double dR = sqrt(dphi*dphi + deta*deta);
    double dpt = (pTG - pTR)/(pTG);
	DR_Lam->Fill(dR);
	DPT_Lam->Fill(dpt);
    if(dR < drdpt && fabs(dpt) < drdpt){
        NmatchLam_Gen[ik1x] = ik2x;
        break;
}else{}}}

for(Int_t ik1x=0; ik1x<GoodTrackFourVector_Lam.size(); ++ik1x){
	double      pTR = GoodTrackFourVector_Lam[ik1x].Pt();
	double      etaR = GoodTrackFourVector_Lam[ik1x].Eta();
	double      phiR = GoodTrackFourVector_Lam[ik1x].Phi();
	bool samegen = false;
	for(Int_t ik1xx=ik1x+1; ik1xx<GoodTrackFourVector_Lam.size(); ++ik1xx){if(NmatchLam_Gen[ik1x] >= 0 && NmatchLam_Gen[ik1xx] >= 0 && NmatchLam_Gen[ik1x] == NmatchLam_Gen[ik1xx]){samegen=true;}}
	if(samegen)continue;
for(Int_t ik2x=0; ik2x<map_GenLam.size(); ++ik2x){
	double pTG = map_GenLam[ik2x].Pt();
	double etaG = map_GenLam[ik2x].Eta();
	double phiG = map_GenLam[ik2x].Phi();
	double dphi = GoodTrackFourVector_Lam[ik1x].DeltaPhi(map_GenLam[ik2x]);//phiR - phiG;
    double deta = etaR - etaG;
    double dR = sqrt(dphi*dphi + deta*deta);
    double dpt = (pTG - pTR)/(pTG);
    if(dR < drdpt && fabs(dpt) < drdpt){
    map_Lam_match[ik1x] = kTRUE;
    if(GoodTrackFourVector_Lam[ik1x].M() > massLAL - nsigmapeak*sigmaLAL && GoodTrackFourVector_Lam[ik1x].M() < massLAL + nsigmapeak*sigmaLAL){map_GenLam_matched.push_back(map_GenLam[ik2x]);}
    break;
}else{}}}

}

if(map_GenALam.size() >= 1 && GoodTrackFourVector_ALam.size() >= 1){

for(Int_t ik1x=0; ik1x<GoodTrackFourVector_ALam.size(); ++ik1x){
	double      pTR = GoodTrackFourVector_ALam[ik1x].Pt();
	double      etaR = GoodTrackFourVector_ALam[ik1x].Eta();
	double      phiR = GoodTrackFourVector_ALam[ik1x].Phi();
    NmatchALam_Gen[ik1x] = -999;
for(Int_t ik2x=0; ik2x<map_GenALam.size(); ++ik2x){
	double pTG = map_GenALam[ik2x].Pt();
	double etaG = map_GenALam[ik2x].Eta();
	double phiG = map_GenALam[ik2x].Phi();
	double dphi = GoodTrackFourVector_ALam[ik1x].DeltaPhi(map_GenALam[ik2x]);//phiR - phiG;
    double deta = etaR - etaG;
    double dR = sqrt(dphi*dphi + deta*deta);
    double dpt = (pTG - pTR)/(pTG);
	DR_ALam->Fill(dR);
	DPT_ALam->Fill(dpt);
    if(dR < drdpt && fabs(dpt) < drdpt){
        NmatchALam_Gen[ik1x] = ik2x;
        break;
}else{}}}

for(Int_t ik1x=0; ik1x<GoodTrackFourVector_ALam.size(); ++ik1x){
	double      pTR = GoodTrackFourVector_ALam[ik1x].Pt();
	double      etaR = GoodTrackFourVector_ALam[ik1x].Eta();
	double      phiR = GoodTrackFourVector_ALam[ik1x].Phi();
	bool samegen = false;
	for(Int_t ik1xx=ik1x+1; ik1xx<GoodTrackFourVector_ALam.size(); ++ik1xx){if(NmatchALam_Gen[ik1x] >= 0 && NmatchALam_Gen[ik1xx] >= 0 && NmatchALam_Gen[ik1x] == NmatchALam_Gen[ik1xx]){samegen=true;}}
	if(samegen)continue;
for(Int_t ik2x=0; ik2x<map_GenALam.size(); ++ik2x){
	double pTG = map_GenALam[ik2x].Pt();
	double etaG = map_GenALam[ik2x].Eta();
	double phiG = map_GenALam[ik2x].Phi();
	double dphi = GoodTrackFourVector_ALam[ik1x].DeltaPhi(map_GenALam[ik2x]);//phiR - phiG;
    double deta = etaR - etaG;
    double dR = sqrt(dphi*dphi + deta*deta);
    double dpt = (pTG - pTR)/(pTG);
    if(dR < drdpt && fabs(dpt) < drdpt){
    map_ALam_match[ik1x] = kTRUE;
    if(GoodTrackFourVector_ALam[ik1x].M() > massLAL - nsigmapeak*sigmaLAL && GoodTrackFourVector_ALam[ik1x].M() < massLAL + nsigmapeak*sigmaLAL){map_GenALam_matched.push_back(map_GenALam[ik2x]);}
    break;
}else{}}}

}


if(map_GenK0s_matched.size()>=2){

hbt(map_GenK0s_matched, (Double_t) aux_N_tk_offline, hGen_K0s_K0s_Matched);

ev_z_vtxK0sGen.push_back(aux_vtxz);
ev_ntrkoff_vecGenK0s.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecGenK0s.push_back(map_GenK0s_matched); 

}

if(map_GenLam_matched.size()>=2){

hbt(map_GenLam_matched, (Double_t) aux_N_tk_offline, hGen_Lam_Lam_Matched);

ev_z_vtxLamGen.push_back(aux_vtxz);
ev_ntrkoff_vecGenLam.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecGenLam.push_back(map_GenLam_matched); 

}

if(map_GenALam_matched.size()>=2){

hbt(map_GenALam_matched, (Double_t) aux_N_tk_offline, hGen_ALam_ALam_Matched);

ev_z_vtxALamGen.push_back(aux_vtxz);
ev_ntrkoff_vecGenALam.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecGenALam.push_back(map_GenALam_matched); 

}

if(map_GenLam_matched.size()>=1 && map_GenALam_matched.size()>=1){

hbtOS(map_GenLam_matched, map_GenALam_matched, (Double_t) aux_N_tk_offline, hGen_Lam_ALam_Matched);

ev_z_vtxLamALamGen.push_back(aux_vtxz);
ev_ntrkoff_vecGenLamALam.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecGenLamALam.push_back(map_GenLam_matched); 
ev_GoodTrackFourVector_vecGenALamLam.push_back(map_GenALam_matched); 

}

if(map_GenK0s_matched.size()>=1 && map_GenLam_matched.size()>=1){

hbtOS(map_GenK0s_matched, map_GenLam_matched, (Double_t) aux_N_tk_offline, hGen_K0s_Lam_Matched);

ev_z_vtxK0sLamGen.push_back(aux_vtxz);
ev_ntrkoff_vecGenK0sLam.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecGenK0sLam.push_back(map_GenK0s_matched); 
ev_GoodTrackFourVector_vecGenLamK0s.push_back(map_GenLam_matched); 

}

if(map_GenK0s_matched.size()>=1 && map_GenALam_matched.size()>=1){

hbtOS(map_GenK0s_matched, map_GenALam_matched, (Double_t) aux_N_tk_offline, hGen_K0s_ALam_Matched);

ev_z_vtxK0sALamGen.push_back(aux_vtxz);
ev_ntrkoff_vecGenK0sALam.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecGenK0sALam.push_back(map_GenK0s_matched); 
ev_GoodTrackFourVector_vecGenALamK0s.push_back(map_GenALam_matched); 

}

}

map_GenK0s.clear();
map_GenLam.clear();
map_GenALam.clear();


if(DoFeedDown_Xi){

//Xi-Lam
if(GoodTrackFourVector_Xi_Lamb.size()>=1 && GoodTrackFourVector_Lam.size()>=1){

for(Int_t ik=0; ik<GoodTrackFourVector_Lam.size(); ik++){

double pt_lam = GoodTrackFourVector_Lam[ik].Pt();
double m_lam = GoodTrackFourVector_Lam[ik].M();
double phi_lam = GoodTrackFourVector_Lam[ik].Phi();
double y_lam = GoodTrackFourVector_Lam[ik].Rapidity();

double pt_lam_proton = GoodTrackFourVectord1_Lam[ik].Pt();
double phi_lam_proton = GoodTrackFourVectord1_Lam[ik].Phi();
double eta_lam_proton = GoodTrackFourVectord1_Lam[ik].Eta();

double pt_lam_pion = GoodTrackFourVectord2_Lam[ik].Pt();
double phi_lam_pion = GoodTrackFourVectord2_Lam[ik].Phi();
double eta_lam_pion = GoodTrackFourVectord2_Lam[ik].Eta();


for(Int_t iik=0; iik<GoodTrackFourVector_Xi_Lamb.size(); iik++){

double pt_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].Rapidity();

double pt_lam_proton_feed = GoodTrackFourVector_Xi_Lamb_proton[iik].Pt();
double phi_lam_proton_feed = GoodTrackFourVector_Xi_Lamb_proton[iik].Phi();
double eta_lam_proton_feed = GoodTrackFourVector_Xi_Lamb_proton[iik].Eta();

double pt_lam_pion_feed = GoodTrackFourVector_Xi_Lamb_pion[iik].Pt();
double phi_lam_pion_feed = GoodTrackFourVector_Xi_Lamb_pion[iik].Phi();
double eta_lam_pion_feed = GoodTrackFourVector_Xi_Lamb_pion[iik].Eta();

/*
cout << "====================================" << endl;
cout << "---------------Xi+Lam---------------" << endl;
cout << "====================================" << endl;

cout << "pT(lamb): " << pt_lam << " ; pT(lamb->casc): " << pt_lam_feed << endl;
cout << "M(lamb): " << m_lam << " ; M(lamb->casc): " << m_lam_feed << endl;
cout << "phi(lamb): " << phi_lam << " ; phi(lamb->casc): " << phi_lam_feed << endl;
cout << "y(lamb): " << y_lam << " ; y(lamb->casc): " << y_lam_feed << endl;

cout << "pT(lamb_d1): " << pt_lam_proton << " ; pT(lamb_d1->casc): " << pt_lam_proton_feed << endl;
cout << "phi(lamb_d1): " << phi_lam_proton << " ; phi(lamb_d1->casc): " << phi_lam_proton_feed << endl;
cout << "eta(lamb_d1): " << eta_lam_proton << " ; eta(lamb_d1->casc): " << eta_lam_proton_feed << endl;

cout << "pT(lamb_d2): " << pt_lam_pion << " ; pT(lamb_d2->casc): " << pt_lam_pion_feed << endl;
cout << "phi(lamb_d2): " << phi_lam_pion << " ; phi(lamb_d2->casc): " << phi_lam_pion_feed << endl;
cout << "eta(lamb_d2): " << eta_lam_pion << " ; eta(lamb_d2->casc): " << eta_lam_pion_feed << endl;
*/

if(pt_lam==pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed 
&& pt_lam_proton==pt_lam_proton_feed && phi_lam_proton==phi_lam_proton_feed && eta_lam_proton==eta_lam_proton_feed
&& pt_lam_pion==pt_lam_pion_feed && phi_lam_pion==phi_lam_pion_feed && eta_lam_pion==eta_lam_pion_feed){
cout << "non-prompt Lambda" << endl;
if(GoodTrackFourVector_Xi_Lamb_bool_K[iik])lamb_from_feed_xi_K[ik] = kTRUE;
if(GoodTrackFourVector_Xi_Lamb_bool_L[iik])lamb_from_feed_xi_L[ik] = kTRUE;
if(GoodTrackFourVector_Xi_Lamb_bool_AL[iik])lamb_from_feed_xi_AL[ik] = kTRUE;
}}}

}

//Xi-ALam
if(GoodTrackFourVector_Xi_Lamb.size()>=1 && GoodTrackFourVector_ALam.size()>=1){

for(Int_t ik=0; ik<GoodTrackFourVector_ALam.size(); ik++){

double pt_lam = GoodTrackFourVector_ALam[ik].Pt();
double m_lam = GoodTrackFourVector_ALam[ik].M();
double phi_lam = GoodTrackFourVector_ALam[ik].Phi();
double y_lam = GoodTrackFourVector_ALam[ik].Rapidity();

double pt_lam_proton = GoodTrackFourVectord1_ALam[ik].Pt();
double phi_lam_proton = GoodTrackFourVectord1_ALam[ik].Phi();
double eta_lam_proton = GoodTrackFourVectord1_ALam[ik].Eta();

double pt_lam_pion = GoodTrackFourVectord2_ALam[ik].Pt();
double phi_lam_pion = GoodTrackFourVectord2_ALam[ik].Phi();
double eta_lam_pion = GoodTrackFourVectord2_ALam[ik].Eta();

for(Int_t iik=0; iik<GoodTrackFourVector_Xi_Lamb.size(); iik++){

double pt_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_Xi_Lamb[iik].Rapidity();

double pt_lam_proton_feed = GoodTrackFourVector_Xi_Lamb_proton[iik].Pt();
double phi_lam_proton_feed = GoodTrackFourVector_Xi_Lamb_proton[iik].Phi();
double eta_lam_proton_feed = GoodTrackFourVector_Xi_Lamb_proton[iik].Eta();

double pt_lam_pion_feed = GoodTrackFourVector_Xi_Lamb_pion[iik].Pt();
double phi_lam_pion_feed = GoodTrackFourVector_Xi_Lamb_pion[iik].Phi();
double eta_lam_pion_feed = GoodTrackFourVector_Xi_Lamb_pion[iik].Eta();


/*

cout << "====================================" << endl;
cout << "---------------Xi+ALam---------------" << endl;
cout << "====================================" << endl;

cout << "pT(lamb): " << pt_lam << " ; pT(lamb->casc): " << pt_lam_feed << endl;
cout << "M(lamb): " << m_lam << " ; M(lamb->casc): " << m_lam_feed << endl;
cout << "phi(lamb): " << phi_lam << " ; phi(lamb->casc): " << phi_lam_feed << endl;
cout << "y(lamb): " << y_lam << " ; y(lamb->casc): " << y_lam_feed << endl;

cout << "pT(lamb_d1): " << pt_lam_proton << " ; pT(lamb_d1->casc): " << pt_lam_proton_feed << endl;
cout << "phi(lamb_d1): " << phi_lam_proton << " ; phi(lamb_d1->casc): " << phi_lam_proton_feed << endl;
cout << "eta(lamb_d1): " << eta_lam_proton << " ; eta(lamb_d1->casc): " << eta_lam_proton_feed << endl;

cout << "pT(lamb_d2): " << pt_lam_pion << " ; pT(lamb_d2->casc): " << pt_lam_pion_feed << endl;
cout << "phi(lamb_d2): " << phi_lam_pion << " ; phi(lamb_d2->casc): " << phi_lam_pion_feed << endl;
cout << "eta(lamb_d2): " << eta_lam_pion << " ; eta(lamb_d2->casc): " << eta_lam_pion_feed << endl;

*/

if(pt_lam==pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed 
&& pt_lam_proton==pt_lam_proton_feed && phi_lam_proton==phi_lam_proton_feed && eta_lam_proton==eta_lam_proton_feed
&& pt_lam_pion==pt_lam_pion_feed && phi_lam_pion==phi_lam_pion_feed && eta_lam_pion==eta_lam_pion_feed){
cout << "non-prompt Lambda" << endl;
if(GoodTrackFourVector_Xi_Lamb_bool_K[iik])alamb_from_feed_xi_K[ik] = kTRUE;
if(GoodTrackFourVector_Xi_Lamb_bool_L[iik])alamb_from_feed_xi_L[ik] = kTRUE;
if(GoodTrackFourVector_Xi_Lamb_bool_AL[iik])alamb_from_feed_xi_AL[ik] = kTRUE;
}}}

}

//AXi-Lam
if(GoodTrackFourVector_AXi_Lamb.size()>=1 && GoodTrackFourVector_Lam.size()>=1){

for(Int_t ik=0; ik<GoodTrackFourVector_Lam.size(); ik++){

double pt_lam = GoodTrackFourVector_Lam[ik].Pt();
double m_lam = GoodTrackFourVector_Lam[ik].M();
double phi_lam = GoodTrackFourVector_Lam[ik].Phi();
double y_lam = GoodTrackFourVector_Lam[ik].Rapidity();

double pt_lam_proton = GoodTrackFourVectord1_Lam[ik].Pt();
double phi_lam_proton = GoodTrackFourVectord1_Lam[ik].Phi();
double eta_lam_proton = GoodTrackFourVectord1_Lam[ik].Eta();

double pt_lam_pion = GoodTrackFourVectord2_Lam[ik].Pt();
double phi_lam_pion = GoodTrackFourVectord2_Lam[ik].Phi();
double eta_lam_pion = GoodTrackFourVectord2_Lam[ik].Eta();


for(Int_t iik=0; iik<GoodTrackFourVector_AXi_Lamb.size(); iik++){

double pt_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].Rapidity();

double pt_lam_proton_feed = GoodTrackFourVector_AXi_Lamb_proton[iik].Pt();
double phi_lam_proton_feed = GoodTrackFourVector_AXi_Lamb_proton[iik].Phi();
double eta_lam_proton_feed = GoodTrackFourVector_AXi_Lamb_proton[iik].Eta();

double pt_lam_pion_feed = GoodTrackFourVector_AXi_Lamb_pion[iik].Pt();
double phi_lam_pion_feed = GoodTrackFourVector_AXi_Lamb_pion[iik].Phi();
double eta_lam_pion_feed = GoodTrackFourVector_AXi_Lamb_pion[iik].Eta();

/*
cout << "====================================" << endl;
cout << "---------------AXi+Lam---------------" << endl;
cout << "====================================" << endl;

cout << "pT(lamb): " << pt_lam << " ; pT(lamb->casc): " << pt_lam_feed << endl;
cout << "M(lamb): " << m_lam << " ; M(lamb->casc): " << m_lam_feed << endl;
cout << "phi(lamb): " << phi_lam << " ; phi(lamb->casc): " << phi_lam_feed << endl;
cout << "y(lamb): " << y_lam << " ; y(lamb->casc): " << y_lam_feed << endl;

cout << "pT(lamb_d1): " << pt_lam_proton << " ; pT(lamb_d1->casc): " << pt_lam_proton_feed << endl;
cout << "phi(lamb_d1): " << phi_lam_proton << " ; phi(lamb_d1->casc): " << phi_lam_proton_feed << endl;
cout << "eta(lamb_d1): " << eta_lam_proton << " ; eta(lamb_d1->casc): " << eta_lam_proton_feed << endl;

cout << "pT(lamb_d2): " << pt_lam_pion << " ; pT(lamb_d2->casc): " << pt_lam_pion_feed << endl;
cout << "phi(lamb_d2): " << phi_lam_pion << " ; phi(lamb_d2->casc): " << phi_lam_pion_feed << endl;
cout << "eta(lamb_d2): " << eta_lam_pion << " ; eta(lamb_d2->casc): " << eta_lam_pion_feed << endl;
*/

if(pt_lam==pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed 
&& pt_lam_proton==pt_lam_proton_feed && phi_lam_proton==phi_lam_proton_feed && eta_lam_proton==eta_lam_proton_feed
&& pt_lam_pion==pt_lam_pion_feed && phi_lam_pion==phi_lam_pion_feed && eta_lam_pion==eta_lam_pion_feed){
cout << "non-prompt Lambda" << endl;
if(GoodTrackFourVector_AXi_Lamb_bool_K[iik])lamb_from_feed_axi_K[ik] = kTRUE;
if(GoodTrackFourVector_AXi_Lamb_bool_L[iik])lamb_from_feed_axi_L[ik] = kTRUE;
if(GoodTrackFourVector_AXi_Lamb_bool_AL[iik])lamb_from_feed_axi_AL[ik] = kTRUE;
}}}

}

//AXi-ALam
if(GoodTrackFourVector_AXi_Lamb.size()>=1 && GoodTrackFourVector_ALam.size()>=1){

for(Int_t ik=0; ik<GoodTrackFourVector_ALam.size(); ik++){

double pt_lam = GoodTrackFourVector_ALam[ik].Pt();
double m_lam = GoodTrackFourVector_ALam[ik].M();
double phi_lam = GoodTrackFourVector_ALam[ik].Phi();
double y_lam = GoodTrackFourVector_ALam[ik].Rapidity();

double pt_lam_proton = GoodTrackFourVectord1_ALam[ik].Pt();
double phi_lam_proton = GoodTrackFourVectord1_ALam[ik].Phi();
double eta_lam_proton = GoodTrackFourVectord1_ALam[ik].Eta();

double pt_lam_pion = GoodTrackFourVectord2_ALam[ik].Pt();
double phi_lam_pion = GoodTrackFourVectord2_ALam[ik].Phi();
double eta_lam_pion = GoodTrackFourVectord2_ALam[ik].Eta();


for(Int_t iik=0; iik<GoodTrackFourVector_AXi_Lamb.size(); iik++){

double pt_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_AXi_Lamb[iik].Rapidity();

double pt_lam_proton_feed = GoodTrackFourVector_AXi_Lamb_proton[iik].Pt();
double phi_lam_proton_feed = GoodTrackFourVector_AXi_Lamb_proton[iik].Phi();
double eta_lam_proton_feed = GoodTrackFourVector_AXi_Lamb_proton[iik].Eta();

double pt_lam_pion_feed = GoodTrackFourVector_AXi_Lamb_pion[iik].Pt();
double phi_lam_pion_feed = GoodTrackFourVector_AXi_Lamb_pion[iik].Phi();
double eta_lam_pion_feed = GoodTrackFourVector_AXi_Lamb_pion[iik].Eta();

/*
cout << "====================================" << endl;
cout << "---------------AXi+ALam---------------" << endl;
cout << "====================================" << endl;

cout << "pT(lamb): " << pt_lam << " ; pT(lamb->casc): " << pt_lam_feed << endl;
cout << "M(lamb): " << m_lam << " ; M(lamb->casc): " << m_lam_feed << endl;
cout << "phi(lamb): " << phi_lam << " ; phi(lamb->casc): " << phi_lam_feed << endl;
cout << "y(lamb): " << y_lam << " ; y(lamb->casc): " << y_lam_feed << endl;

cout << "pT(lamb_d1): " << pt_lam_proton << " ; pT(lamb_d1->casc): " << pt_lam_proton_feed << endl;
cout << "phi(lamb_d1): " << phi_lam_proton << " ; phi(lamb_d1->casc): " << phi_lam_proton_feed << endl;
cout << "eta(lamb_d1): " << eta_lam_proton << " ; eta(lamb_d1->casc): " << eta_lam_proton_feed << endl;

cout << "pT(lamb_d2): " << pt_lam_pion << " ; pT(lamb_d2->casc): " << pt_lam_pion_feed << endl;
cout << "phi(lamb_d2): " << phi_lam_pion << " ; phi(lamb_d2->casc): " << phi_lam_pion_feed << endl;
cout << "eta(lamb_d2): " << eta_lam_pion << " ; eta(lamb_d2->casc): " << eta_lam_pion_feed << endl;
*/

if(pt_lam==pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed 
&& pt_lam_proton==pt_lam_proton_feed && phi_lam_proton==phi_lam_proton_feed && eta_lam_proton==eta_lam_proton_feed
&& pt_lam_pion==pt_lam_pion_feed && phi_lam_pion==phi_lam_pion_feed && eta_lam_pion==eta_lam_pion_feed){
cout << "non-prompt Lambda" << endl;
if(GoodTrackFourVector_AXi_Lamb_bool_K[iik])alamb_from_feed_axi_K[ik] = kTRUE;
if(GoodTrackFourVector_AXi_Lamb_bool_L[iik])alamb_from_feed_axi_L[ik] = kTRUE;
if(GoodTrackFourVector_AXi_Lamb_bool_AL[iik])alamb_from_feed_axi_AL[ik] = kTRUE;
}}}

}


}

if(DoFeedDown_Om){

//Om-Lam
if(GoodTrackFourVector_Om_Lamb.size()>=1 && GoodTrackFourVector_Lam.size()>=1){

for(Int_t ik=0; ik<GoodTrackFourVector_Lam.size(); ik++){

double pt_lam = GoodTrackFourVector_Lam[ik].Pt();
double m_lam = GoodTrackFourVector_Lam[ik].M();
double phi_lam = GoodTrackFourVector_Lam[ik].Phi();
double y_lam = GoodTrackFourVector_Lam[ik].Rapidity();

double pt_lam_proton = GoodTrackFourVectord1_Lam[ik].Pt();
double phi_lam_proton = GoodTrackFourVectord1_Lam[ik].Phi();
double eta_lam_proton = GoodTrackFourVectord1_Lam[ik].Eta();

double pt_lam_pion = GoodTrackFourVectord2_Lam[ik].Pt();
double phi_lam_pion = GoodTrackFourVectord2_Lam[ik].Phi();
double eta_lam_pion = GoodTrackFourVectord2_Lam[ik].Eta();


for(Int_t iik=0; iik<GoodTrackFourVector_Om_Lamb.size(); iik++){

double pt_lam_feed = GoodTrackFourVector_Om_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_Om_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_Om_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_Om_Lamb[iik].Rapidity();

double pt_lam_proton_feed = GoodTrackFourVector_Om_Lamb_proton[iik].Pt();
double phi_lam_proton_feed = GoodTrackFourVector_Om_Lamb_proton[iik].Phi();
double eta_lam_proton_feed = GoodTrackFourVector_Om_Lamb_proton[iik].Eta();

double pt_lam_pion_feed = GoodTrackFourVector_Om_Lamb_kaon[iik].Pt();
double phi_lam_pion_feed = GoodTrackFourVector_Om_Lamb_kaon[iik].Phi();
double eta_lam_pion_feed = GoodTrackFourVector_Om_Lamb_kaon[iik].Eta();

/*
cout << "====================================" << endl;
cout << "---------------Om+Lam---------------" << endl;
cout << "====================================" << endl;

cout << "pT(lamb): " << pt_lam << " ; pT(lamb->omega): " << pt_lam_feed << endl;
cout << "M(lamb): " << m_lam << " ; M(lamb->omega): " << m_lam_feed << endl;
cout << "phi(lamb): " << phi_lam << " ; phi(lamb->omega): " << phi_lam_feed << endl;
cout << "y(lamb): " << y_lam << " ; y(lamb->omega): " << y_lam_feed << endl;

cout << "pT(lamb_d1): " << pt_lam_proton << " ; pT(lamb_d1->omega): " << pt_lam_proton_feed << endl;
cout << "phi(lamb_d1): " << phi_lam_proton << " ; phi(lamb_d1->omega): " << phi_lam_proton_feed << endl;
cout << "eta(lamb_d1): " << eta_lam_proton << " ; eta(lamb_d1->omega): " << eta_lam_proton_feed << endl;

cout << "pT(lamb_d2): " << pt_lam_pion << " ; pT(lamb_d2->omega): " << pt_lam_pion_feed << endl;
cout << "phi(lamb_d2): " << phi_lam_pion << " ; phi(lamb_d2->omega): " << phi_lam_pion_feed << endl;
cout << "eta(lamb_d2): " << eta_lam_pion << " ; eta(lamb_d2->omega): " << eta_lam_pion_feed << endl;
*/

if(pt_lam==pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed 
&& pt_lam_proton==pt_lam_proton_feed && phi_lam_proton==phi_lam_proton_feed && eta_lam_proton==eta_lam_proton_feed
&& pt_lam_pion==pt_lam_pion_feed && phi_lam_pion==phi_lam_pion_feed && eta_lam_pion==eta_lam_pion_feed){
cout << "non-prompt Lambda" << endl;
if(GoodTrackFourVector_Om_Lamb_bool_K[iik])lamb_from_feed_om_K[ik] = kTRUE;
if(GoodTrackFourVector_Om_Lamb_bool_L[iik])lamb_from_feed_om_L[ik] = kTRUE;
if(GoodTrackFourVector_Om_Lamb_bool_AL[iik])lamb_from_feed_om_AL[ik] = kTRUE;
}}}

}

//Om-ALam
if(GoodTrackFourVector_Om_Lamb.size()>=1 && GoodTrackFourVector_ALam.size()>=1){

for(Int_t ik=0; ik<GoodTrackFourVector_ALam.size(); ik++){

double pt_lam = GoodTrackFourVector_ALam[ik].Pt();
double m_lam = GoodTrackFourVector_ALam[ik].M();
double phi_lam = GoodTrackFourVector_ALam[ik].Phi();
double y_lam = GoodTrackFourVector_ALam[ik].Rapidity();

double pt_lam_proton = GoodTrackFourVectord1_ALam[ik].Pt();
double phi_lam_proton = GoodTrackFourVectord1_ALam[ik].Phi();
double eta_lam_proton = GoodTrackFourVectord1_ALam[ik].Eta();

double pt_lam_pion = GoodTrackFourVectord2_ALam[ik].Pt();
double phi_lam_pion = GoodTrackFourVectord2_ALam[ik].Phi();
double eta_lam_pion = GoodTrackFourVectord2_ALam[ik].Eta();


for(Int_t iik=0; iik<GoodTrackFourVector_Om_Lamb.size(); iik++){

double pt_lam_feed = GoodTrackFourVector_Om_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_Om_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_Om_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_Om_Lamb[iik].Rapidity();

double pt_lam_proton_feed = GoodTrackFourVector_Om_Lamb_proton[iik].Pt();
double phi_lam_proton_feed = GoodTrackFourVector_Om_Lamb_proton[iik].Phi();
double eta_lam_proton_feed = GoodTrackFourVector_Om_Lamb_proton[iik].Eta();

double pt_lam_pion_feed = GoodTrackFourVector_Om_Lamb_kaon[iik].Pt();
double phi_lam_pion_feed = GoodTrackFourVector_Om_Lamb_kaon[iik].Phi();
double eta_lam_pion_feed = GoodTrackFourVector_Om_Lamb_kaon[iik].Eta();

/*
cout << "====================================" << endl;
cout << "---------------Om+ALam---------------" << endl;
cout << "====================================" << endl;

cout << "pT(lamb): " << pt_lam << " ; pT(lamb->omega): " << pt_lam_feed << endl;
cout << "M(lamb): " << m_lam << " ; M(lamb->omega): " << m_lam_feed << endl;
cout << "phi(lamb): " << phi_lam << " ; phi(lamb->omega): " << phi_lam_feed << endl;
cout << "y(lamb): " << y_lam << " ; y(lamb->omega): " << y_lam_feed << endl;

cout << "pT(lamb_d1): " << pt_lam_proton << " ; pT(lamb_d1->omega): " << pt_lam_proton_feed << endl;
cout << "phi(lamb_d1): " << phi_lam_proton << " ; phi(lamb_d1->omega): " << phi_lam_proton_feed << endl;
cout << "eta(lamb_d1): " << eta_lam_proton << " ; eta(lamb_d1->omega): " << eta_lam_proton_feed << endl;

cout << "pT(lamb_d2): " << pt_lam_pion << " ; pT(lamb_d2->omega): " << pt_lam_pion_feed << endl;
cout << "phi(lamb_d2): " << phi_lam_pion << " ; phi(lamb_d2->omega): " << phi_lam_pion_feed << endl;
cout << "eta(lamb_d2): " << eta_lam_pion << " ; eta(lamb_d2->omega): " << eta_lam_pion_feed << endl;
*/

if(pt_lam==pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed 
&& pt_lam_proton==pt_lam_proton_feed && phi_lam_proton==phi_lam_proton_feed && eta_lam_proton==eta_lam_proton_feed
&& pt_lam_pion==pt_lam_pion_feed && phi_lam_pion==phi_lam_pion_feed && eta_lam_pion==eta_lam_pion_feed){
cout << "non-prompt Lambda" << endl;
if(GoodTrackFourVector_Om_Lamb_bool_K[iik])alamb_from_feed_om_K[ik] = kTRUE;
if(GoodTrackFourVector_Om_Lamb_bool_L[iik])alamb_from_feed_om_L[ik] = kTRUE;
if(GoodTrackFourVector_Om_Lamb_bool_AL[iik])alamb_from_feed_om_AL[ik] = kTRUE;
}}}

}

//AOm-Lam
if(GoodTrackFourVector_AOm_Lamb.size()>=1 && GoodTrackFourVector_Lam.size()>=1){

for(Int_t ik=0; ik<GoodTrackFourVector_Lam.size(); ik++){

double pt_lam = GoodTrackFourVector_Lam[ik].Pt();
double m_lam = GoodTrackFourVector_Lam[ik].M();
double phi_lam = GoodTrackFourVector_Lam[ik].Phi();
double y_lam = GoodTrackFourVector_Lam[ik].Rapidity();

double pt_lam_proton = GoodTrackFourVectord1_Lam[ik].Pt();
double phi_lam_proton = GoodTrackFourVectord1_Lam[ik].Phi();
double eta_lam_proton = GoodTrackFourVectord1_Lam[ik].Eta();

double pt_lam_pion = GoodTrackFourVectord2_Lam[ik].Pt();
double phi_lam_pion = GoodTrackFourVectord2_Lam[ik].Phi();
double eta_lam_pion = GoodTrackFourVectord2_Lam[ik].Eta();


for(Int_t iik=0; iik<GoodTrackFourVector_AOm_Lamb.size(); iik++){

double pt_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].Rapidity();

double pt_lam_proton_feed = GoodTrackFourVector_AOm_Lamb_proton[iik].Pt();
double phi_lam_proton_feed = GoodTrackFourVector_AOm_Lamb_proton[iik].Phi();
double eta_lam_proton_feed = GoodTrackFourVector_AOm_Lamb_proton[iik].Eta();

double pt_lam_pion_feed = GoodTrackFourVector_AOm_Lamb_kaon[iik].Pt();
double phi_lam_pion_feed = GoodTrackFourVector_AOm_Lamb_kaon[iik].Phi();
double eta_lam_pion_feed = GoodTrackFourVector_AOm_Lamb_kaon[iik].Eta();

/*
cout << "====================================" << endl;
cout << "---------------AOm+Lam---------------" << endl;
cout << "====================================" << endl;

cout << "pT(lamb): " << pt_lam << " ; pT(lamb->omega): " << pt_lam_feed << endl;
cout << "M(lamb): " << m_lam << " ; M(lamb->omega): " << m_lam_feed << endl;
cout << "phi(lamb): " << phi_lam << " ; phi(lamb->omega): " << phi_lam_feed << endl;
cout << "y(lamb): " << y_lam << " ; y(lamb->omega): " << y_lam_feed << endl;

cout << "pT(lamb_d1): " << pt_lam_proton << " ; pT(lamb_d1->omega): " << pt_lam_proton_feed << endl;
cout << "phi(lamb_d1): " << phi_lam_proton << " ; phi(lamb_d1->omega): " << phi_lam_proton_feed << endl;
cout << "eta(lamb_d1): " << eta_lam_proton << " ; eta(lamb_d1->omega): " << eta_lam_proton_feed << endl;

cout << "pT(lamb_d2): " << pt_lam_pion << " ; pT(lamb_d2->omega): " << pt_lam_pion_feed << endl;
cout << "phi(lamb_d2): " << phi_lam_pion << " ; phi(lamb_d2->omega): " << phi_lam_pion_feed << endl;
cout << "eta(lamb_d2): " << eta_lam_pion << " ; eta(lamb_d2->omega): " << eta_lam_pion_feed << endl;
*/

if(pt_lam==pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed 
&& pt_lam_proton==pt_lam_proton_feed && phi_lam_proton==phi_lam_proton_feed && eta_lam_proton==eta_lam_proton_feed
&& pt_lam_pion==pt_lam_pion_feed && phi_lam_pion==phi_lam_pion_feed && eta_lam_pion==eta_lam_pion_feed){
cout << "non-prompt Lambda" << endl;
if(GoodTrackFourVector_AOm_Lamb_bool_K[iik])lamb_from_feed_aom_K[ik] = kTRUE;
if(GoodTrackFourVector_AOm_Lamb_bool_L[iik])lamb_from_feed_aom_L[ik] = kTRUE;
if(GoodTrackFourVector_AOm_Lamb_bool_AL[iik])lamb_from_feed_aom_AL[ik] = kTRUE;
}}}

}

//AOm-ALam
if(GoodTrackFourVector_AOm_Lamb.size()>=1 && GoodTrackFourVector_ALam.size()>=1){

for(Int_t ik=0; ik<GoodTrackFourVector_ALam.size(); ik++){

double pt_lam = GoodTrackFourVector_ALam[ik].Pt();
double m_lam = GoodTrackFourVector_ALam[ik].M();
double phi_lam = GoodTrackFourVector_ALam[ik].Phi();
double y_lam = GoodTrackFourVector_ALam[ik].Rapidity();

double pt_lam_proton = GoodTrackFourVectord1_ALam[ik].Pt();
double phi_lam_proton = GoodTrackFourVectord1_ALam[ik].Phi();
double eta_lam_proton = GoodTrackFourVectord1_ALam[ik].Eta();

double pt_lam_pion = GoodTrackFourVectord2_ALam[ik].Pt();
double phi_lam_pion = GoodTrackFourVectord2_ALam[ik].Phi();
double eta_lam_pion = GoodTrackFourVectord2_ALam[ik].Eta();


for(Int_t iik=0; iik<GoodTrackFourVector_AOm_Lamb.size(); iik++){

double pt_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].Pt();
double m_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].M();
double phi_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].Phi();
double y_lam_feed = GoodTrackFourVector_AOm_Lamb[iik].Rapidity();

double pt_lam_proton_feed = GoodTrackFourVector_AOm_Lamb_proton[iik].Pt();
double phi_lam_proton_feed = GoodTrackFourVector_AOm_Lamb_proton[iik].Phi();
double eta_lam_proton_feed = GoodTrackFourVector_AOm_Lamb_proton[iik].Eta();

double pt_lam_pion_feed = GoodTrackFourVector_AOm_Lamb_kaon[iik].Pt();
double phi_lam_pion_feed = GoodTrackFourVector_AOm_Lamb_kaon[iik].Phi();
double eta_lam_pion_feed = GoodTrackFourVector_AOm_Lamb_kaon[iik].Eta();

/*
cout << "====================================" << endl;
cout << "---------------AOm+ALam---------------" << endl;
cout << "====================================" << endl;

cout << "pT(lamb): " << pt_lam << " ; pT(lamb->omega): " << pt_lam_feed << endl;
cout << "M(lamb): " << m_lam << " ; M(lamb->omega): " << m_lam_feed << endl;
cout << "phi(lamb): " << phi_lam << " ; phi(lamb->omega): " << phi_lam_feed << endl;
cout << "y(lamb): " << y_lam << " ; y(lamb->omega): " << y_lam_feed << endl;

cout << "pT(lamb_d1): " << pt_lam_proton << " ; pT(lamb_d1->omega): " << pt_lam_proton_feed << endl;
cout << "phi(lamb_d1): " << phi_lam_proton << " ; phi(lamb_d1->omega): " << phi_lam_proton_feed << endl;
cout << "eta(lamb_d1): " << eta_lam_proton << " ; eta(lamb_d1->omega): " << eta_lam_proton_feed << endl;

cout << "pT(lamb_d2): " << pt_lam_pion << " ; pT(lamb_d2->omega): " << pt_lam_pion_feed << endl;
cout << "phi(lamb_d2): " << phi_lam_pion << " ; phi(lamb_d2->omega): " << phi_lam_pion_feed << endl;
cout << "eta(lamb_d2): " << eta_lam_pion << " ; eta(lamb_d2->omega): " << eta_lam_pion_feed << endl;
*/

if(pt_lam==pt_lam_feed && m_lam == m_lam_feed &&  phi_lam == phi_lam_feed && y_lam == y_lam_feed 
&& pt_lam_proton==pt_lam_proton_feed && phi_lam_proton==phi_lam_proton_feed && eta_lam_proton==eta_lam_proton_feed
&& pt_lam_pion==pt_lam_pion_feed && phi_lam_pion==phi_lam_pion_feed && eta_lam_pion==eta_lam_pion_feed){
cout << "non-prompt Lambda" << endl;
if(GoodTrackFourVector_AOm_Lamb_bool_K[iik])alamb_from_feed_aom_K[ik] = kTRUE;
if(GoodTrackFourVector_AOm_Lamb_bool_L[iik])alamb_from_feed_aom_L[ik] = kTRUE;
if(GoodTrackFourVector_AOm_Lamb_bool_AL[iik])alamb_from_feed_aom_AL[ik] = kTRUE;
}}}

}


}

std::vector<double> GoodTrackFourVector_Xi_Lamb_proton_chi2_T_K;
std::vector<double> GoodTrackFourVector_Xi_Lamb_pion_chi2_T_K;
std::vector<double> GoodTrackFourVector_AXi_Lamb_proton_chi2_T_K;
std::vector<double> GoodTrackFourVector_AXi_Lamb_pion_chi2_T_K;

std::vector<double> GoodTrackFourVector_Om_Lamb_proton_chi2_T_K;
std::vector<double> GoodTrackFourVector_Om_Lamb_kaon_chi2_T_K;
std::vector<double> GoodTrackFourVector_AOm_Lamb_proton_chi2_T_K;
std::vector<double> GoodTrackFourVector_AOm_Lamb_kaon_chi2_T_K;

std::vector<double> GoodTrackFourVector_Xi_Lamb_proton_chi2_T_L;
std::vector<double> GoodTrackFourVector_Xi_Lamb_pion_chi2_T_L;
std::vector<double> GoodTrackFourVector_AXi_Lamb_proton_chi2_T_L;
std::vector<double> GoodTrackFourVector_AXi_Lamb_pion_chi2_T_L;

std::vector<double> GoodTrackFourVector_Om_Lamb_proton_chi2_T_L;
std::vector<double> GoodTrackFourVector_Om_Lamb_kaon_chi2_T_L;
std::vector<double> GoodTrackFourVector_AOm_Lamb_proton_chi2_T_L;
std::vector<double> GoodTrackFourVector_AOm_Lamb_kaon_chi2_T_L;

std::vector<double> GoodTrackFourVector_Xi_Lamb_proton_chi2_T_AL;
std::vector<double> GoodTrackFourVector_Xi_Lamb_pion_chi2_T_AL;
std::vector<double> GoodTrackFourVector_AXi_Lamb_proton_chi2_T_AL;
std::vector<double> GoodTrackFourVector_AXi_Lamb_pion_chi2_T_AL;

std::vector<double> GoodTrackFourVector_Om_Lamb_proton_chi2_T_AL;
std::vector<double> GoodTrackFourVector_Om_Lamb_kaon_chi2_T_AL;
std::vector<double> GoodTrackFourVector_AOm_Lamb_proton_chi2_T_AL;
std::vector<double> GoodTrackFourVector_AOm_Lamb_kaon_chi2_T_AL;


if(DoFeedDown_Xi){

for(Int_t ii=0; ii<GoodTrackFourVector_Xi_Lamb.size(); ii++){
double chi2_1 = GoodTrackFourVector_Xi_Lamb_proton_chi2[ii];
double chi2_2 = GoodTrackFourVector_Xi_Lamb_pion_chi2[ii];

if(GoodTrackFourVector_Xi_Lamb_bool_K[ii] == kTRUE && GoodTrackFourVector_Xi_Lamb_bool_L[ii] == kTRUE && GoodTrackFourVector_Xi_Lamb_bool_AL[ii] == kTRUE && GoodTrackFourVector_Xi_bool[ii] == kTRUE){pT_Lamb_from_Xi_T->Fill(GoodTrackFourVector_Xi[ii].Pt());XXi_mass_LambMom->Fill(GoodTrackFourVector_Xi[ii].M());}

if(GoodTrackFourVector_Xi_Lamb_bool_K[ii]){
GoodTrackFourVector_Xi_Lamb_proton_chi2_T_K.push_back(chi2_1);
GoodTrackFourVector_Xi_Lamb_pion_chi2_T_K.push_back(chi2_2);
}

if(GoodTrackFourVector_Xi_Lamb_bool_L[ii]){
GoodTrackFourVector_Xi_Lamb_proton_chi2_T_L.push_back(chi2_1);
GoodTrackFourVector_Xi_Lamb_pion_chi2_T_L.push_back(chi2_2);
}

if(GoodTrackFourVector_Xi_Lamb_bool_AL[ii]){
GoodTrackFourVector_Xi_Lamb_proton_chi2_T_AL.push_back(chi2_1);
GoodTrackFourVector_Xi_Lamb_pion_chi2_T_AL.push_back(chi2_2);
}

}

for(Int_t ii=0; ii<GoodTrackFourVector_AXi_Lamb.size(); ii++){
double chi2_1 = GoodTrackFourVector_AXi_Lamb_proton_chi2[ii];
double chi2_2 = GoodTrackFourVector_AXi_Lamb_pion_chi2[ii];

if(GoodTrackFourVector_AXi_Lamb_bool_K[ii] == kTRUE && GoodTrackFourVector_AXi_Lamb_bool_L[ii] == kTRUE && GoodTrackFourVector_AXi_Lamb_bool_AL[ii] && GoodTrackFourVector_AXi_bool[ii]){pT_Lamb_from_Xi_T->Fill(GoodTrackFourVector_AXi[ii].Pt());AXXi_mass_LambMom->Fill(GoodTrackFourVector_AXi[ii].M());}

if(GoodTrackFourVector_AXi_Lamb_bool_K[ii]){
GoodTrackFourVector_AXi_Lamb_proton_chi2_T_K.push_back(chi2_1);
GoodTrackFourVector_AXi_Lamb_pion_chi2_T_K.push_back(chi2_2);
}

if(GoodTrackFourVector_AXi_Lamb_bool_L[ii]){
GoodTrackFourVector_AXi_Lamb_proton_chi2_T_L.push_back(chi2_1);
GoodTrackFourVector_AXi_Lamb_pion_chi2_T_L.push_back(chi2_2);
}

if(GoodTrackFourVector_AXi_Lamb_bool_AL[ii]){
GoodTrackFourVector_AXi_Lamb_proton_chi2_T_AL.push_back(chi2_1);
GoodTrackFourVector_AXi_Lamb_pion_chi2_T_AL.push_back(chi2_2);
}
}

}

if(DoFeedDown_Om){

for(Int_t ii=0; ii<GoodTrackFourVector_Om_Lamb.size(); ii++){
double chi2_1 = GoodTrackFourVector_Om_Lamb_proton_chi2[ii];
double chi2_2 = GoodTrackFourVector_Om_Lamb_kaon_chi2[ii];

if(GoodTrackFourVector_Om_Lamb_bool_K[ii] == kTRUE && GoodTrackFourVector_Om_Lamb_bool_L[ii] == kTRUE && GoodTrackFourVector_Om_Lamb_bool_AL[ii] == kTRUE && GoodTrackFourVector_Om_bool[ii] == kTRUE){pT_Lamb_from_Om_T->Fill(GoodTrackFourVector_Om[ii].Pt());OOm_mass_LambMom->Fill(GoodTrackFourVector_Om[ii].M());}

if(GoodTrackFourVector_Om_Lamb_bool_K[ii]){
GoodTrackFourVector_Om_Lamb_proton_chi2_T_K.push_back(chi2_1);
GoodTrackFourVector_Om_Lamb_kaon_chi2_T_K.push_back(chi2_2);
}

if(GoodTrackFourVector_Om_Lamb_bool_L[ii]){
GoodTrackFourVector_Om_Lamb_proton_chi2_T_L.push_back(chi2_1);
GoodTrackFourVector_Om_Lamb_kaon_chi2_T_L.push_back(chi2_2);
}

if(GoodTrackFourVector_Om_Lamb_bool_AL[ii]){
GoodTrackFourVector_Om_Lamb_proton_chi2_T_AL.push_back(chi2_1);
GoodTrackFourVector_Om_Lamb_kaon_chi2_T_AL.push_back(chi2_2);
}

}

for(Int_t ii=0; ii<GoodTrackFourVector_AOm_Lamb.size(); ii++){
double chi2_1 = GoodTrackFourVector_AOm_Lamb_proton_chi2[ii];
double chi2_2 = GoodTrackFourVector_AOm_Lamb_kaon_chi2[ii];

if(GoodTrackFourVector_AOm_Lamb_bool_K[ii] == kTRUE && GoodTrackFourVector_AOm_Lamb_bool_L[ii] == kTRUE && GoodTrackFourVector_AOm_Lamb_bool_AL[ii] == kTRUE && GoodTrackFourVector_AOm_bool[ii] == kTRUE){pT_Lamb_from_Om_T->Fill(GoodTrackFourVector_AOm[ii].Pt());AOOm_mass_LambMom->Fill(GoodTrackFourVector_AOm[ii].M());}

if(GoodTrackFourVector_AOm_Lamb_bool_K[ii]){
GoodTrackFourVector_AOm_Lamb_proton_chi2_T_K.push_back(chi2_1);
GoodTrackFourVector_AOm_Lamb_kaon_chi2_T_K.push_back(chi2_2);
}

if(GoodTrackFourVector_AOm_Lamb_bool_L[ii]){
GoodTrackFourVector_AOm_Lamb_proton_chi2_T_L.push_back(chi2_1);
GoodTrackFourVector_AOm_Lamb_kaon_chi2_T_L.push_back(chi2_2);
}

if(GoodTrackFourVector_AOm_Lamb_bool_AL[ii]){
GoodTrackFourVector_AOm_Lamb_proton_chi2_T_AL.push_back(chi2_1);
GoodTrackFourVector_AOm_Lamb_kaon_chi2_T_AL.push_back(chi2_2);
}
}

}


/*========================= K0sK0s =========================*/

//4-vectors for true particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorT_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1T_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2T_K0s;
std::vector<double> chi21T_K0s;
std::vector<double> chi22T_K0s;
//4-vectors for true particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorF_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1F_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2F_K0s;
std::vector<double> chi21F_K0s;
std::vector<double> chi22F_K0s;
//4-vectors for true matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTM_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1TM_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2TM_K0s;
std::vector<double> chi21TM_K0s;
std::vector<double> chi22TM_K0s;
//4-vectors for true unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTU_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1TU_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2TU_K0s;
std::vector<double> chi21TU_K0s;
std::vector<double> chi22TU_K0s;
//4-vectors for fake matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFM_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1FM_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2FM_K0s;
std::vector<double> chi21FM_K0s;
std::vector<double> chi22FM_K0s;
//4-vectors for fake unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFU_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1FU_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2FU_K0s;
std::vector<double> chi21FU_K0s;
std::vector<double> chi22FU_K0s;

/*========================= K0sLam =========================*/

//K0sLam
//4-vectors for true particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorT_K0s_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1T_K0s_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2T_K0s_Lam;
std::vector<double> chi21T_K0s_Lam;
std::vector<double> chi22T_K0s_Lam;
//4-vectors for true particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorF_K0s_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1F_K0s_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2F_K0s_Lam;
std::vector<double> chi21F_K0s_Lam;
std::vector<double> chi22F_K0s_Lam;
//4-vectors for true matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTM_K0s_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1TM_K0s_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2TM_K0s_Lam;
std::vector<double> chi21TM_K0s_Lam;
std::vector<double> chi22TM_K0s_Lam;
//4-vectors for true unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTU_K0s_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1TU_K0s_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2TU_K0s_Lam;
std::vector<double> chi21TU_K0s_Lam;
std::vector<double> chi22TU_K0s_Lam;
//4-vectors for fake matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFM_K0s_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1FM_K0s_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2FM_K0s_Lam;
std::vector<double> chi21FM_K0s_Lam;
std::vector<double> chi22FM_K0s_Lam;
//4-vectors for fake unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFU_K0s_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1FU_K0s_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2FU_K0s_Lam;
std::vector<double> chi21FU_K0s_Lam;
std::vector<double> chi22FU_K0s_Lam;


/*========================= K0sALam =========================*/
//K0sALam
//4-vectors for true particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorT_K0s_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1T_K0s_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2T_K0s_ALam;
std::vector<double> chi21T_K0s_ALam;
std::vector<double> chi22T_K0s_ALam;
//4-vectors for true particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorF_K0s_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1F_K0s_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2F_K0s_ALam;
std::vector<double> chi21F_K0s_ALam;
std::vector<double> chi22F_K0s_ALam;
//4-vectors for true matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTM_K0s_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1TM_K0s_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2TM_K0s_ALam;
std::vector<double> chi21TM_K0s_ALam;
std::vector<double> chi22TM_K0s_ALam;
//4-vectors for true unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTU_K0s_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1TU_K0s_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2TU_K0s_ALam;
std::vector<double> chi21TU_K0s_ALam;
std::vector<double> chi22TU_K0s_ALam;
//4-vectors for fake matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFM_K0s_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1FM_K0s_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2FM_K0s_ALam;
std::vector<double> chi21FM_K0s_ALam;
std::vector<double> chi22FM_K0s_ALam;
//4-vectors for fake unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFU_K0s_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1FU_K0s_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2FU_K0s_ALam;
std::vector<double> chi21FU_K0s_ALam;
std::vector<double> chi22FU_K0s_ALam;

for(Int_t jj=0; jj<GoodTrackFourVector_K0s.size(); jj++){

double pt = GoodTrackFourVector_K0s[jj].Pt();
double eta = GoodTrackFourVector_K0s[jj].Eta();
double phi = GoodTrackFourVector_K0s[jj].Phi();
double m = GoodTrackFourVector_K0s[jj].M();

double ptd1 = GoodTrackFourVectord1_K0s[jj].Pt();
double etad1 = GoodTrackFourVectord1_K0s[jj].Eta();
double phid1 = GoodTrackFourVectord1_K0s[jj].Phi();
double md1 = GoodTrackFourVectord1_K0s[jj].M();

double ptd2 = GoodTrackFourVectord2_K0s[jj].Pt();
double etad2 = GoodTrackFourVectord2_K0s[jj].Eta();
double phid2 = GoodTrackFourVectord2_K0s[jj].Phi();
double md2 = GoodTrackFourVectord2_K0s[jj].M();

double chi21X = chi21_K0s[jj];
double chi22X = chi22_K0s[jj];

TLorentzVector pvectorM;
pvectorM.SetPtEtaPhiM(pt,eta,phi,m);
TLorentzVector pvectorD1;
pvectorD1.SetPtEtaPhiM(ptd1,etad1,phid1,md1);
TLorentzVector pvectorD2;
pvectorD2.SetPtEtaPhiM(ptd2,etad2,phid2,md2);


if(isMC){

bool V0from2 = false;
Int_t njetmatch = 0;
Int_t NJet = -999;

for(Int_t ajj=0; ajj<jetsize; ajj++){ 
TLorentzVector pvectorj;
pvectorj.SetPtEtaPhiM(Jet_pt->at(ajj),Jet_eta->at(ajj),Jet_phi->at(ajj),Jet_mass->at(ajj));

double deltaPhi=pvectorM.DeltaPhi(pvectorj);
double deltaEta=pvectorM.Eta() -  pvectorj.Eta();
double DR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

double deltaPhid1=pvectorD1.DeltaPhi(pvectorj);
double deltaEtad1=pvectorD1.Eta() -  pvectorj.Eta();
double DRd1 = sqrt(deltaPhid1*deltaPhid1 + deltaEtad1*deltaEtad1);

double deltaPhid2=pvectorD2.DeltaPhi(pvectorj);
double deltaEtad2=pvectorD2.Eta() -  pvectorj.Eta();
double DRd2 = sqrt(deltaPhid2*deltaPhid2 + deltaEtad2*deltaEtad2);

cone_K0s->Fill(DR); cone_K0sD1->Fill(DRd1); cone_K0sD2->Fill(DRd2);
if(map_K0s_chi2[jj]==kTRUE){cone_K0sT->Fill(DR); cone_K0sTD1->Fill(DRd1); cone_K0sTD2->Fill(DRd2);
if(DR <= 0.4){njetmatch=njetmatch+1; NJet = ajj; if(!removeV0from2jets){break;}}}
}

if(njetmatch>1){V0from2=true;}
if(V0from2)continue;

if(map_K0s_chi2[jj]==kTRUE){

if(njetmatch==1){GoodTrackFourVector_K0s_Jet.push_back(pvectorM); GoodTrackFourVector_K0s_NJet.push_back(NJet);}
if(njetmatch==0){GoodTrackFourVector_K0s_NotaJet.push_back(pvectorM);}

}

}

//K0sK0s
if(remove_all_mis){

if(map_K0s_chi2[jj]==kTRUE && map_K0s_chi2_KL[jj]==kTRUE && map_K0s_chi2_KAL[jj]==kTRUE){
GoodTrackFourVectorT_K0s.push_back(pvectorM);
GoodTrackFourVectord1T_K0s.push_back(pvectorD1);
GoodTrackFourVectord2T_K0s.push_back(pvectorD2);
chi21T_K0s.push_back(chi21X);
chi22T_K0s.push_back(chi22X);
GoodTrackFourVectorT_K0s_Lam.push_back(pvectorM);
GoodTrackFourVectord1T_K0s_Lam.push_back(pvectorD1);
GoodTrackFourVectord2T_K0s_Lam.push_back(pvectorD2);
chi21T_K0s_Lam.push_back(chi21X);
chi22T_K0s_Lam.push_back(chi22X);
GoodTrackFourVectorT_K0s_ALam.push_back(pvectorM);
GoodTrackFourVectord1T_K0s_ALam.push_back(pvectorD1);
GoodTrackFourVectord2T_K0s_ALam.push_back(pvectorD2);
chi21T_K0s_ALam.push_back(chi21X);
chi22T_K0s_ALam.push_back(chi22X);
if(map_K0s_match[jj]==kTRUE){
GoodTrackFourVectorTM_K0s.push_back(pvectorM);
GoodTrackFourVectord1TM_K0s.push_back(pvectorD1);
GoodTrackFourVectord2TM_K0s.push_back(pvectorD2);
chi21TM_K0s.push_back(chi21X);
chi22TM_K0s.push_back(chi22X);
GoodTrackFourVectorTM_K0s_Lam.push_back(pvectorM);
GoodTrackFourVectord1TM_K0s_Lam.push_back(pvectorD1);
GoodTrackFourVectord2TM_K0s_Lam.push_back(pvectorD2);
chi21TM_K0s_Lam.push_back(chi21X);
chi22TM_K0s_Lam.push_back(chi22X);
GoodTrackFourVectorTM_K0s_ALam.push_back(pvectorM);
GoodTrackFourVectord1TM_K0s_ALam.push_back(pvectorD1);
GoodTrackFourVectord2TM_K0s_ALam.push_back(pvectorD2);
chi21TM_K0s_ALam.push_back(chi21X);
chi22TM_K0s_ALam.push_back(chi22X);
}else{
GoodTrackFourVectorTU_K0s.push_back(pvectorM);
GoodTrackFourVectord1TU_K0s.push_back(pvectorD1);
GoodTrackFourVectord2TU_K0s.push_back(pvectorD2);
chi21TU_K0s.push_back(chi21X);
chi22TU_K0s.push_back(chi22X);
GoodTrackFourVectorTU_K0s_Lam.push_back(pvectorM);
GoodTrackFourVectord1TU_K0s_Lam.push_back(pvectorD1);
GoodTrackFourVectord2TU_K0s_Lam.push_back(pvectorD2);
chi21TU_K0s_Lam.push_back(chi21X);
chi22TU_K0s_Lam.push_back(chi22X);
GoodTrackFourVectorTU_K0s_ALam.push_back(pvectorM);
GoodTrackFourVectord1TU_K0s_ALam.push_back(pvectorD1);
GoodTrackFourVectord2TU_K0s_ALam.push_back(pvectorD2);
chi21TU_K0s_ALam.push_back(chi21X);
chi22TU_K0s_ALam.push_back(chi22X);
}}else{
GoodTrackFourVectorF_K0s.push_back(pvectorM);
GoodTrackFourVectord1F_K0s.push_back(pvectorD1);
GoodTrackFourVectord2F_K0s.push_back(pvectorD2);
chi21F_K0s.push_back(chi21X);
chi22F_K0s.push_back(chi22X);
GoodTrackFourVectorF_K0s_Lam.push_back(pvectorM);
GoodTrackFourVectord1F_K0s_Lam.push_back(pvectorD1);
GoodTrackFourVectord2F_K0s_Lam.push_back(pvectorD2);
chi21F_K0s_Lam.push_back(chi21X);
chi22F_K0s_Lam.push_back(chi22X);
GoodTrackFourVectorF_K0s_ALam.push_back(pvectorM);
GoodTrackFourVectord1F_K0s_ALam.push_back(pvectorD1);
GoodTrackFourVectord2F_K0s_ALam.push_back(pvectorD2);
chi21F_K0s_ALam.push_back(chi21X);
chi22F_K0s_ALam.push_back(chi22X);
if(map_K0s_match[jj]==kTRUE){
GoodTrackFourVectorFM_K0s.push_back(pvectorM);
GoodTrackFourVectord1FM_K0s.push_back(pvectorD1);
GoodTrackFourVectord2FM_K0s.push_back(pvectorD2);
chi21FM_K0s.push_back(chi21X);
chi22FM_K0s.push_back(chi22X);
GoodTrackFourVectorFM_K0s_Lam.push_back(pvectorM);
GoodTrackFourVectord1FM_K0s_Lam.push_back(pvectorD1);
GoodTrackFourVectord2FM_K0s_Lam.push_back(pvectorD2);
chi21FM_K0s_Lam.push_back(chi21X);
chi22FM_K0s_Lam.push_back(chi22X);
GoodTrackFourVectorFM_K0s_ALam.push_back(pvectorM);
GoodTrackFourVectord1FM_K0s_ALam.push_back(pvectorD1);
GoodTrackFourVectord2FM_K0s_ALam.push_back(pvectorD2);
chi21FM_K0s_ALam.push_back(chi21X);
chi22FM_K0s_ALam.push_back(chi22X);
}else{
GoodTrackFourVectorFU_K0s.push_back(pvectorM);
GoodTrackFourVectord1FU_K0s.push_back(pvectorD1);
GoodTrackFourVectord2FU_K0s.push_back(pvectorD2);
chi21FU_K0s.push_back(chi21X);
chi22FU_K0s.push_back(chi22X);
GoodTrackFourVectorFU_K0s_Lam.push_back(pvectorM);
GoodTrackFourVectord1FU_K0s_Lam.push_back(pvectorD1);
GoodTrackFourVectord2FU_K0s_Lam.push_back(pvectorD2);
chi21FU_K0s_Lam.push_back(chi21X);
chi22FU_K0s_Lam.push_back(chi22X);
GoodTrackFourVectorFU_K0s_ALam.push_back(pvectorM);
GoodTrackFourVectord1FU_K0s_ALam.push_back(pvectorD1);
GoodTrackFourVectord2FU_K0s_ALam.push_back(pvectorD2);
chi21FU_K0s_ALam.push_back(chi21X);
chi22FU_K0s_ALam.push_back(chi22X);
}}
}else{

if(map_K0s_chi2[jj]==kTRUE){
GoodTrackFourVectorT_K0s.push_back(pvectorM);
GoodTrackFourVectord1T_K0s.push_back(pvectorD1);
GoodTrackFourVectord2T_K0s.push_back(pvectorD2);
chi21T_K0s.push_back(chi21X);
chi22T_K0s.push_back(chi22X);
if(map_K0s_match[jj]==kTRUE){
GoodTrackFourVectorTM_K0s.push_back(pvectorM);
GoodTrackFourVectord1TM_K0s.push_back(pvectorD1);
GoodTrackFourVectord2TM_K0s.push_back(pvectorD2);
chi21TM_K0s.push_back(chi21X);
chi22TM_K0s.push_back(chi22X);
}else{
GoodTrackFourVectorTU_K0s.push_back(pvectorM);
GoodTrackFourVectord1TU_K0s.push_back(pvectorD1);
GoodTrackFourVectord2TU_K0s.push_back(pvectorD2);
chi21TU_K0s.push_back(chi21X);
chi22TU_K0s.push_back(chi22X);
}}else{
GoodTrackFourVectorF_K0s.push_back(pvectorM);
GoodTrackFourVectord1F_K0s.push_back(pvectorD1);
GoodTrackFourVectord2F_K0s.push_back(pvectorD2);
chi21F_K0s.push_back(chi21X);
chi22F_K0s.push_back(chi22X);
if(map_K0s_match[jj]==kTRUE){
GoodTrackFourVectorFM_K0s.push_back(pvectorM);
GoodTrackFourVectord1FM_K0s.push_back(pvectorD1);
GoodTrackFourVectord2FM_K0s.push_back(pvectorD2);
chi21FM_K0s.push_back(chi21X);
chi22FM_K0s.push_back(chi22X);
}else{
GoodTrackFourVectorFU_K0s.push_back(pvectorM);
GoodTrackFourVectord1FU_K0s.push_back(pvectorD1);
GoodTrackFourVectord2FU_K0s.push_back(pvectorD2);
chi21FU_K0s.push_back(chi21X);
chi22FU_K0s.push_back(chi22X);
}}

//K0sLam
if(map_K0s_chi2_KL[jj]==kTRUE){
GoodTrackFourVectorT_K0s_Lam.push_back(pvectorM);
GoodTrackFourVectord1T_K0s_Lam.push_back(pvectorD1);
GoodTrackFourVectord2T_K0s_Lam.push_back(pvectorD2);
chi21T_K0s_Lam.push_back(chi21X);
chi22T_K0s_Lam.push_back(chi22X);
if(map_K0s_match[jj]==kTRUE){
GoodTrackFourVectorTM_K0s_Lam.push_back(pvectorM);
GoodTrackFourVectord1TM_K0s_Lam.push_back(pvectorD1);
GoodTrackFourVectord2TM_K0s_Lam.push_back(pvectorD2);
chi21TM_K0s_Lam.push_back(chi21X);
chi22TM_K0s_Lam.push_back(chi22X);
}else{
GoodTrackFourVectorTU_K0s_Lam.push_back(pvectorM);
GoodTrackFourVectord1TU_K0s_Lam.push_back(pvectorD1);
GoodTrackFourVectord2TU_K0s_Lam.push_back(pvectorD2);
chi21TU_K0s_Lam.push_back(chi21X);
chi22TU_K0s_Lam.push_back(chi22X);
}}else{
GoodTrackFourVectorF_K0s_Lam.push_back(pvectorM);
GoodTrackFourVectord1F_K0s_Lam.push_back(pvectorD1);
GoodTrackFourVectord2F_K0s_Lam.push_back(pvectorD2);
chi21F_K0s_Lam.push_back(chi21X);
chi22F_K0s_Lam.push_back(chi22X);
if(map_K0s_match[jj]==kTRUE){
GoodTrackFourVectorFM_K0s_Lam.push_back(pvectorM);
GoodTrackFourVectord1FM_K0s_Lam.push_back(pvectorD1);
GoodTrackFourVectord2FM_K0s_Lam.push_back(pvectorD2);
chi21FM_K0s_Lam.push_back(chi21X);
chi22FM_K0s_Lam.push_back(chi22X);
}else{
GoodTrackFourVectorFU_K0s_Lam.push_back(pvectorM);
GoodTrackFourVectord1FU_K0s_Lam.push_back(pvectorD1);
GoodTrackFourVectord2FU_K0s_Lam.push_back(pvectorD2);
chi21FU_K0s_Lam.push_back(chi21X);
chi22FU_K0s_Lam.push_back(chi22X);
}}

//K0sALam
if(map_K0s_chi2_KAL[jj]==kTRUE){
GoodTrackFourVectorT_K0s_ALam.push_back(pvectorM);
GoodTrackFourVectord1T_K0s_ALam.push_back(pvectorD1);
GoodTrackFourVectord2T_K0s_ALam.push_back(pvectorD2);
chi21T_K0s_ALam.push_back(chi21X);
chi22T_K0s_ALam.push_back(chi22X);
if(map_K0s_match[jj]==kTRUE){
GoodTrackFourVectorTM_K0s_ALam.push_back(pvectorM);
GoodTrackFourVectord1TM_K0s_ALam.push_back(pvectorD1);
GoodTrackFourVectord2TM_K0s_ALam.push_back(pvectorD2);
chi21TM_K0s_ALam.push_back(chi21X);
chi22TM_K0s_ALam.push_back(chi22X);
}else{
GoodTrackFourVectorTU_K0s_ALam.push_back(pvectorM);
GoodTrackFourVectord1TU_K0s_ALam.push_back(pvectorD1);
GoodTrackFourVectord2TU_K0s_ALam.push_back(pvectorD2);
chi21TU_K0s_ALam.push_back(chi21X);
chi22TU_K0s_ALam.push_back(chi22X);
}}else{
GoodTrackFourVectorF_K0s_ALam.push_back(pvectorM);
GoodTrackFourVectord1F_K0s_ALam.push_back(pvectorD1);
GoodTrackFourVectord2F_K0s_ALam.push_back(pvectorD2);
chi21F_K0s_ALam.push_back(chi21X);
chi22F_K0s_ALam.push_back(chi22X);
if(map_K0s_match[jj]==kTRUE){
GoodTrackFourVectorFM_K0s_ALam.push_back(pvectorM);
GoodTrackFourVectord1FM_K0s_ALam.push_back(pvectorD1);
GoodTrackFourVectord2FM_K0s_ALam.push_back(pvectorD2);
chi21FM_K0s_ALam.push_back(chi21X);
chi22FM_K0s_ALam.push_back(chi22X);
}else{
GoodTrackFourVectorFU_K0s_ALam.push_back(pvectorM);
GoodTrackFourVectord1FU_K0s_ALam.push_back(pvectorD1);
GoodTrackFourVectord2FU_K0s_ALam.push_back(pvectorD2);
chi21FU_K0s_ALam.push_back(chi21X);
chi22FU_K0s_ALam.push_back(chi22X);
}}

}

}

//--------------------true--------------------

if(GoodTrackFourVectorT_K0s.size()>=2){

nev_K0sT->Fill(1); 

for(unsigned int iy=0; iy<GoodTrackFourVectorT_K0s.size();iy++){K0s_Mass->Fill(GoodTrackFourVectorT_K0s[iy].M(),Ntk_Vz_weight);}

Mass_dep(GoodTrackFourVectorT_K0s, (Double_t) aux_N_tk_offline, hMass_K0s_K0s_T,Ntk_Vz_weight);
Mass_dep_sep(GoodTrackFourVectorT_K0s, (Double_t) aux_N_tk_offline, hMass_K0s_1_T, hMass_K0s_2_T, Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorT_K0s,chi21T_K0s,chi22T_K0s,V0chi2d1_diff_T_K0s,V0chi2d2_diff_T_K0s,V0chi2d12_diff_T_K0s,Ntk_Vz_weight);
    
hbt_reco(GoodTrackFourVectorT_K0s, GoodTrackFourVectord1T_K0s, GoodTrackFourVectord2T_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_T, hpeakcos_K0s_K0s_T, hpeakd1_K0s_K0s_T, hpeakd2_K0s_K0s_T, hpeakd12_K0s_K0s_T, hpeak_rot_K0s_K0s_T, hpeak_inv_K0s_K0s_T, hside_K0s_K0s_T, hside_rot_K0s_K0s_T, hside_inv_K0s_K0s_T, hpeakside_K0s_K0s_T, hpeakside_rot_K0s_K0s_T, hpeakside_inv_K0s_K0s_T, hsideL_K0s_K0s_T, hsideR_K0s_K0s_T, hpeaksideL_K0s_K0s_T, hpeaksideR_K0s_K0s_T, p, effhisto_K0s, nev_K0s_ssT, nev_K0s_bbT, nev_K0s_sbT,Ntk_Vz_weight);

hdibaryon(GoodTrackFourVectorT_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, nsigmapeak, nsigmaside, h_Hf2_K0s_K0s, h_Hf2_rot_K0s_K0s, h_Hf2_inv_K0s_K0s, h_Hf2_K0s_K0s_side, h_Hf2_rot_K0s_K0s_side, h_Hf2_inv_K0s_K0s_side, h_Hf2_K0s_K0s_peakside, h_Hf2_rot_K0s_K0s_peakside, h_Hf2_inv_K0s_K0s_peakside, Ntk_Vz_weight);

//weight
weight_K0sK0s.push_back(Ntk_Vz_weight);

//Mix Random

ev_z_vtxK0s.push_back(aux_vtxz);
ev_ntrkoff_vecK0s.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0s.push_back(GoodTrackFourVectorT_K0s); 

//Eta Mix

Double_t aux_etaMix_w = Eta_weight;
std::pair<Double_t, std::vector<TLorentzVector>> aux_pair_GoodTrackFourVector_etaMixWeight = make_pair(aux_etaMix_w, GoodTrackFourVectorT_K0s);
ev_GoodTrackFourVector_etaMixWeight_vecK0s.push_back(aux_pair_GoodTrackFourVector_etaMixWeight);
std::pair<Double_t,int> aux_pair_ntrkoff_etaMixWeight = make_pair(aux_etaMix_w,aux_N_tk_offline);
ev_ntrkoff_etaMixWeight_vecK0s.push_back(aux_pair_ntrkoff_etaMixWeight); 

}//K0sK0s

//--------------------false--------------------

if(GoodTrackFourVectorF_K0s.size()>=2){

nev_K0sF->Fill(1); 

Mass_dep(GoodTrackFourVectorF_K0s, (Double_t) aux_N_tk_offline, hMass_K0s_K0s_F,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorF_K0s,chi21F_K0s,chi22F_K0s,V0chi2d1_diff_F_K0s,V0chi2d2_diff_F_K0s,V0chi2d12_diff_F_K0s,Ntk_Vz_weight);
    
hbt_reco_fake(GoodTrackFourVectorF_K0s, GoodTrackFourVectord1F_K0s, GoodTrackFourVectord2F_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_F, hpeakd1_K0s_K0s_F, hpeakd2_K0s_K0s_F, hpeakd12_K0s_K0s_F, hside_K0s_K0s_F, hpeakside_K0s_K0s_F, hpeak_K0s_K0s_F_bias, hpeak_K0s_K0s_F_biasd1, hpeak_K0s_K0s_F_biasd2, hpeak_K0s_K0s_F_biasd12,Ntk_Vz_weight);

}

//--------------------TF-------------------

if(GoodTrackFourVectorT_K0s.size()>=1 && GoodTrackFourVectorF_K0s.size()>=1){

nev_K0sTF->Fill(1);

get_chi2plot_TF(GoodTrackFourVectorT_K0s, chi21T_K0s, chi22T_K0s, GoodTrackFourVectorF_K0s, chi21F_K0s, chi22F_K0s, V0chi2d1_diff_TF_K0s,V0chi2d2_diff_TF_K0s,V0chi2d12_diff_TF_K0s,Ntk_Vz_weight);

hbt_reco_truefake(GoodTrackFourVectorT_K0s, GoodTrackFourVectord1T_K0s, GoodTrackFourVectord2T_K0s, GoodTrackFourVectorF_K0s, GoodTrackFourVectord1F_K0s, GoodTrackFourVectord2F_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_TF, hpeakd1_K0s_K0s_TF, hpeakd2_K0s_K0s_TF, hpeakd12_K0s_K0s_TF, hside_K0s_K0s_TF, hpeakside_K0s_K0s_TF, hpeak_K0s_K0s_TF_bias, hpeak_K0s_K0s_TF_biasd1, hpeak_K0s_K0s_TF_biasd2, hpeak_K0s_K0s_TF_biasd12,Ntk_Vz_weight);
    
}

if(isMC){

//-------------------True Matched-----------------------

if(GoodTrackFourVectorTM_K0s.size()>=2){

Mass_dep(GoodTrackFourVectorTM_K0s, (Double_t) aux_N_tk_offline, hMass_K0s_K0s_TM,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorTM_K0s,chi21TM_K0s,chi22TM_K0s,V0chi2d1_diff_TM_K0s,V0chi2d2_diff_TM_K0s,V0chi2d12_diff_TM_K0s,Ntk_Vz_weight);
    
hbt_MM_UU(GoodTrackFourVectorTM_K0s, GoodTrackFourVectord1TM_K0s, GoodTrackFourVectord2TM_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_TM, hpeakd1_K0s_K0s_TM, hpeakd2_K0s_K0s_TM, hpeakd12_K0s_K0s_TM, hpeak_K0s_K0s_TM_bias, hpeak_K0s_K0s_TM_biasd1, hpeak_K0s_K0s_TM_biasd2, hpeak_K0s_K0s_TM_biasd12,Ntk_Vz_weight);

ev_z_vtxK0s_TM.push_back(aux_vtxz);
ev_ntrkoff_vecK0s_TM.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0s_TM.push_back(GoodTrackFourVectorTM_K0s); 

}

//-------------------True Unmatched-----------------------

if(GoodTrackFourVectorTU_K0s.size()>=2){

Mass_dep(GoodTrackFourVectorTU_K0s, (Double_t) aux_N_tk_offline, hMass_K0s_K0s_TU,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorTU_K0s,chi21TU_K0s,chi22TU_K0s,V0chi2d1_diff_TU_K0s,V0chi2d2_diff_TU_K0s,V0chi2d12_diff_TU_K0s,Ntk_Vz_weight);
    
hbt_MM_UU(GoodTrackFourVectorTU_K0s, GoodTrackFourVectord1TU_K0s, GoodTrackFourVectord2TU_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_TU, hpeakd1_K0s_K0s_TU, hpeakd2_K0s_K0s_TU, hpeakd12_K0s_K0s_TU, hpeak_K0s_K0s_TU_bias, hpeak_K0s_K0s_TU_biasd1, hpeak_K0s_K0s_TU_biasd2, hpeak_K0s_K0s_TU_biasd12,Ntk_Vz_weight);

ev_z_vtxK0s_TU.push_back(aux_vtxz);
ev_ntrkoff_vecK0s_TU.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0s_TU.push_back(GoodTrackFourVectorTU_K0s); 

}

//-------------------Fake Matched-----------------------

if(GoodTrackFourVectorFM_K0s.size()>=2){

Mass_dep(GoodTrackFourVectorFM_K0s, (Double_t) aux_N_tk_offline, hMass_K0s_K0s_FM,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorFM_K0s,chi21FM_K0s,chi22FM_K0s,V0chi2d1_diff_FM_K0s,V0chi2d2_diff_FM_K0s,V0chi2d12_diff_FM_K0s,Ntk_Vz_weight);
    
hbt_MM_UU(GoodTrackFourVectorFM_K0s, GoodTrackFourVectord1FM_K0s, GoodTrackFourVectord2FM_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_FM, hpeakd1_K0s_K0s_FM, hpeakd2_K0s_K0s_FM, hpeakd12_K0s_K0s_FM, hpeak_K0s_K0s_FM_bias, hpeak_K0s_K0s_FM_biasd1, hpeak_K0s_K0s_FM_biasd2, hpeak_K0s_K0s_FM_biasd12,Ntk_Vz_weight);

}

//-------------------Fake Unmatched-----------------------

if(GoodTrackFourVectorFU_K0s.size()>=2){

Mass_dep(GoodTrackFourVectorFU_K0s, (Double_t) aux_N_tk_offline, hMass_K0s_K0s_FU,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorFU_K0s,chi21FU_K0s,chi22FU_K0s,V0chi2d1_diff_FU_K0s,V0chi2d2_diff_FU_K0s,V0chi2d12_diff_FU_K0s,Ntk_Vz_weight);
    
hbt_MM_UU(GoodTrackFourVectorFU_K0s, GoodTrackFourVectord1FU_K0s, GoodTrackFourVectord2FU_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_FU, hpeakd1_K0s_K0s_FU, hpeakd2_K0s_K0s_FU, hpeakd12_K0s_K0s_FU, hpeak_K0s_K0s_FU_bias, hpeak_K0s_K0s_FU_biasd1, hpeak_K0s_K0s_FU_biasd2, hpeak_K0s_K0s_FU_biasd12,Ntk_Vz_weight);

}


//-------------------True Matched + True Unmatched-----------------------

if(GoodTrackFourVectorTM_K0s.size()>=1 && GoodTrackFourVectorTU_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_K0s, chi21TM_K0s, chi22TM_K0s, GoodTrackFourVectorTU_K0s, chi21TU_K0s, chi22TU_K0s, V0chi2d1_diff_TM_TU_K0s, V0chi2d2_diff_TM_TU_K0s, V0chi2d12_diff_TM_TU_K0s,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTM_K0s, GoodTrackFourVectord1TM_K0s, GoodTrackFourVectord2TM_K0s, GoodTrackFourVectorTU_K0s, GoodTrackFourVectord1TU_K0s, GoodTrackFourVectord2TU_K0s,(Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_TM_TU, hpeakd1_K0s_K0s_TM_TU, hpeakd2_K0s_K0s_TM_TU, hpeakd12_K0s_K0s_TM_TU, hpeak_K0s_K0s_TM_TU_bias, hpeak_K0s_K0s_TM_TU_biasd1, hpeak_K0s_K0s_TM_TU_biasd2, hpeak_K0s_K0s_TM_TU_biasd12,Ntk_Vz_weight);

ev_z_vtxK0s_TM_TU.push_back(aux_vtxz);
ev_ntrkoff_vecK0s_TM_TU.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0s_TM_TU.push_back(GoodTrackFourVectorTM_K0s); 
ev_GoodTrackFourVector_vecK0s_TU_TM.push_back(GoodTrackFourVectorTU_K0s); 

}

//-------------------True Matched + Fake Matched-----------------------

if(GoodTrackFourVectorTM_K0s.size()>=1 && GoodTrackFourVectorFM_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_K0s, chi21TM_K0s, chi22TM_K0s, GoodTrackFourVectorFM_K0s, chi21FM_K0s, chi22FM_K0s, V0chi2d1_diff_TM_FM_K0s, V0chi2d2_diff_TM_FM_K0s, V0chi2d12_diff_TM_FM_K0s,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTM_K0s, GoodTrackFourVectord1TM_K0s, GoodTrackFourVectord2TM_K0s, GoodTrackFourVectorFM_K0s, GoodTrackFourVectord1FM_K0s, GoodTrackFourVectord2FM_K0s,(Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_TM_FM, hpeakd1_K0s_K0s_TM_FM, hpeakd2_K0s_K0s_TM_FM, hpeakd12_K0s_K0s_TM_FM, hpeak_K0s_K0s_TM_FM_bias, hpeak_K0s_K0s_TM_FM_biasd1, hpeak_K0s_K0s_TM_FM_biasd2, hpeak_K0s_K0s_TM_FM_biasd12,Ntk_Vz_weight);

}

//-------------------True Matched + Fake Unmatched-----------------------

if(GoodTrackFourVectorTM_K0s.size()>=1 && GoodTrackFourVectorFU_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_K0s, chi21TM_K0s, chi22TM_K0s, GoodTrackFourVectorFU_K0s, chi21FU_K0s, chi22FU_K0s, V0chi2d1_diff_TM_FU_K0s, V0chi2d2_diff_TM_FU_K0s, V0chi2d12_diff_TM_FU_K0s,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTM_K0s, GoodTrackFourVectord1TM_K0s, GoodTrackFourVectord2TM_K0s, GoodTrackFourVectorFU_K0s, GoodTrackFourVectord1FU_K0s, GoodTrackFourVectord2FU_K0s,(Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_TM_FU, hpeakd1_K0s_K0s_TM_FU, hpeakd2_K0s_K0s_TM_FU, hpeakd12_K0s_K0s_TM_FU, hpeak_K0s_K0s_TM_FU_bias, hpeak_K0s_K0s_TM_FU_biasd1, hpeak_K0s_K0s_TM_FU_biasd2, hpeak_K0s_K0s_TM_FU_biasd12,Ntk_Vz_weight);

}


//-------------------True Unmatched + Fake Matched-----------------------

if(GoodTrackFourVectorTU_K0s.size()>=1 && GoodTrackFourVectorFM_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_K0s, chi21TU_K0s, chi22TU_K0s, GoodTrackFourVectorFM_K0s, chi21FM_K0s, chi22FM_K0s, V0chi2d1_diff_TU_FM_K0s, V0chi2d2_diff_TU_FM_K0s, V0chi2d12_diff_TU_FM_K0s,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTU_K0s, GoodTrackFourVectord1TU_K0s, GoodTrackFourVectord2TU_K0s, GoodTrackFourVectorFM_K0s, GoodTrackFourVectord1FM_K0s, GoodTrackFourVectord2FM_K0s,(Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_TU_FM, hpeakd1_K0s_K0s_TU_FM, hpeakd2_K0s_K0s_TU_FM, hpeakd12_K0s_K0s_TU_FM, hpeak_K0s_K0s_TU_FM_bias, hpeak_K0s_K0s_TU_FM_biasd1, hpeak_K0s_K0s_TU_FM_biasd2, hpeak_K0s_K0s_TU_FM_biasd12,Ntk_Vz_weight);

}

//-------------------True Unmatched + Fake Unmatched-----------------------

if(GoodTrackFourVectorTU_K0s.size()>=1 && GoodTrackFourVectorFU_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_K0s, chi21TU_K0s, chi22TU_K0s, GoodTrackFourVectorFU_K0s, chi21FU_K0s, chi22FU_K0s, V0chi2d1_diff_TU_FU_K0s, V0chi2d2_diff_TU_FU_K0s, V0chi2d12_diff_TU_FU_K0s,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTU_K0s, GoodTrackFourVectord1TU_K0s, GoodTrackFourVectord2TU_K0s, GoodTrackFourVectorFU_K0s, GoodTrackFourVectord1FU_K0s, GoodTrackFourVectord2FU_K0s,(Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_TU_FU, hpeakd1_K0s_K0s_TU_FU, hpeakd2_K0s_K0s_TU_FU, hpeakd12_K0s_K0s_TU_FU, hpeak_K0s_K0s_TU_FU_bias, hpeak_K0s_K0s_TU_FU_biasd1, hpeak_K0s_K0s_TU_FU_biasd2, hpeak_K0s_K0s_TU_FU_biasd12,Ntk_Vz_weight);

}


//-------------------Fake Matched + Fake Unmatched-----------------------

if(GoodTrackFourVectorFM_K0s.size()>=1 && GoodTrackFourVectorFU_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFM_K0s, chi21FM_K0s, chi22FM_K0s, GoodTrackFourVectorFU_K0s, chi21FU_K0s, chi22FU_K0s, V0chi2d1_diff_FM_FU_K0s, V0chi2d2_diff_FM_FU_K0s, V0chi2d12_diff_FM_FU_K0s,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorFM_K0s, GoodTrackFourVectord1FM_K0s, GoodTrackFourVectord2FM_K0s, GoodTrackFourVectorFU_K0s, GoodTrackFourVectord1FU_K0s, GoodTrackFourVectord2FU_K0s,(Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_FM_FU, hpeakd1_K0s_K0s_FM_FU, hpeakd2_K0s_K0s_FM_FU, hpeakd12_K0s_K0s_FM_FU, hpeak_K0s_K0s_FM_FU_bias, hpeak_K0s_K0s_FM_FU_biasd1, hpeak_K0s_K0s_FM_FU_biasd2, hpeak_K0s_K0s_FM_FU_biasd12,Ntk_Vz_weight);

}

}

/*========================= LambdaLambda =========================*/

//4-vectors for true particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectord1T_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2T_Lam;
std::vector<TLorentzVector> GoodTrackFourVectorT_Lam;
std::vector<double> chi21T_Lam;
std::vector<double> chi22T_Lam;
//4-vectors for fake particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectord1F_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2F_Lam;
std::vector<TLorentzVector> GoodTrackFourVectorF_Lam;
std::vector<double> chi21F_Lam;
std::vector<double> chi22F_Lam;
//4-vectors for true matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTM_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1TM_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2TM_Lam;
std::vector<double> chi21TM_Lam;
std::vector<double> chi22TM_Lam;
//4-vectors for true unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTU_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1TU_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2TU_Lam;
std::vector<double> chi21TU_Lam;
std::vector<double> chi22TU_Lam;
//4-vectors for fake matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFM_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1FM_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2FM_Lam;
std::vector<double> chi21FM_Lam;
std::vector<double> chi22FM_Lam;
//4-vectors for fake unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFU_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1FU_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2FU_Lam;
std::vector<double> chi21FU_Lam;
std::vector<double> chi22FU_Lam;

/*========================= LambdaAntiLambda =========================*/

//4-vectors for true particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectord1T_Lam_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2T_Lam_ALam;
std::vector<TLorentzVector> GoodTrackFourVectorT_Lam_ALam;
std::vector<double> chi21T_Lam_ALam;
std::vector<double> chi22T_Lam_ALam;
//4-vectors for fake particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectord1F_Lam_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2F_Lam_ALam;
std::vector<TLorentzVector> GoodTrackFourVectorF_Lam_ALam;
std::vector<double> chi21F_Lam_ALam;
std::vector<double> chi22F_Lam_ALam;
//4-vectors for true matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTM_Lam_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1TM_Lam_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2TM_Lam_ALam;
std::vector<double> chi21TM_Lam_ALam;
std::vector<double> chi22TM_Lam_ALam;
//4-vectors for true unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTU_Lam_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1TU_Lam_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2TU_Lam_ALam;
std::vector<double> chi21TU_Lam_ALam;
std::vector<double> chi22TU_Lam_ALam;
//4-vectors for fake matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFM_Lam_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1FM_Lam_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2FM_Lam_ALam;
std::vector<double> chi21FM_Lam_ALam;
std::vector<double> chi22FM_Lam_ALam;
//4-vectors for fake unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFU_Lam_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1FU_Lam_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2FU_Lam_ALam;
std::vector<double> chi21FU_Lam_ALam;
std::vector<double> chi22FU_Lam_ALam;

/*========================= LambdaK0s =========================*/

//4-vectors for true particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectord1T_Lam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2T_Lam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectorT_Lam_K0s;
std::vector<double> chi21T_Lam_K0s;
std::vector<double> chi22T_Lam_K0s;
//4-vectors for fake particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectord1F_Lam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2F_Lam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectorF_Lam_K0s;
std::vector<double> chi21F_Lam_K0s;
std::vector<double> chi22F_Lam_K0s;
//4-vectors for true matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTM_Lam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1TM_Lam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2TM_Lam_K0s;
std::vector<double> chi21TM_Lam_K0s;
std::vector<double> chi22TM_Lam_K0s;
//4-vectors for true unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTU_Lam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1TU_Lam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2TU_Lam_K0s;
std::vector<double> chi21TU_Lam_K0s;
std::vector<double> chi22TU_Lam_K0s;
//4-vectors for fake matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFM_Lam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1FM_Lam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2FM_Lam_K0s;
std::vector<double> chi21FM_Lam_K0s;
std::vector<double> chi22FM_Lam_K0s;
//4-vectors for fake unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFU_Lam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1FU_Lam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2FU_Lam_K0s;
std::vector<double> chi21FU_Lam_K0s;
std::vector<double> chi22FU_Lam_K0s;

std::vector<TLorentzVector> GoodTrackFourVectorT_Lam_FD_K;
std::vector<TLorentzVector> GoodTrackFourVectorT_Lam_FD_L;
std::vector<TLorentzVector> GoodTrackFourVectorT_Lam_FD_AL;


for(Int_t jj=0; jj<GoodTrackFourVector_Lam.size(); jj++){

double pt = GoodTrackFourVector_Lam[jj].Pt();
double eta = GoodTrackFourVector_Lam[jj].Eta();
double phi = GoodTrackFourVector_Lam[jj].Phi();
double m = GoodTrackFourVector_Lam[jj].M();

double ptd1 = GoodTrackFourVectord1_Lam[jj].Pt();
double etad1 = GoodTrackFourVectord1_Lam[jj].Eta();
double phid1 = GoodTrackFourVectord1_Lam[jj].Phi();
double md1 = GoodTrackFourVectord1_Lam[jj].M();

double ptd2 = GoodTrackFourVectord2_Lam[jj].Pt();
double etad2 = GoodTrackFourVectord2_Lam[jj].Eta();
double phid2 = GoodTrackFourVectord2_Lam[jj].Phi();
double md2 = GoodTrackFourVectord2_Lam[jj].M();

double chi21X = chi21_Lam[jj];
double chi22X = chi22_Lam[jj];

TLorentzVector pvectorM;
pvectorM.SetPtEtaPhiM(pt,eta,phi,m);
TLorentzVector pvectorD1;
pvectorD1.SetPtEtaPhiM(ptd1,etad1,phid1,md1);
TLorentzVector pvectorD2;
pvectorD2.SetPtEtaPhiM(ptd2,etad2,phid2,md2);

//feed down
if(DoFeedDown_Xi){if(lamb_from_feed_xi_K[jj]==kTRUE || lamb_from_feed_axi_K[jj]==kTRUE || lamb_from_feed_xi_L[jj]==kTRUE || lamb_from_feed_axi_L[jj]==kTRUE || lamb_from_feed_xi_AL[jj]==kTRUE || lamb_from_feed_axi_AL[jj]==kTRUE){pT_Lamb_daughter_Xi->Fill(pt);}}
if(DoFeedDown_Om){if(lamb_from_feed_om_K[jj]==kTRUE || lamb_from_feed_aom_K[jj]==kTRUE || lamb_from_feed_om_L[jj]==kTRUE || lamb_from_feed_aom_L[jj]==kTRUE || lamb_from_feed_om_AL[jj]==kTRUE || lamb_from_feed_aom_AL[jj]==kTRUE){pT_Lamb_daughter_Om->Fill(pt);}}
if(DoFeedDown_Xi){if(lamb_from_feed_xi_K[jj]==kTRUE || lamb_from_feed_axi_K[jj]==kTRUE){GoodTrackFourVectorT_Lam_FD_K.push_back(pvectorM);}}
if(DoFeedDown_Om){if(lamb_from_feed_om_K[jj]==kTRUE || lamb_from_feed_aom_K[jj]==kTRUE){GoodTrackFourVectorT_Lam_FD_K.push_back(pvectorM);}}
if(DoFeedDown_Xi){if(lamb_from_feed_xi_L[jj]==kTRUE || lamb_from_feed_axi_L[jj]==kTRUE){GoodTrackFourVectorT_Lam_FD_L.push_back(pvectorM);}}
if(DoFeedDown_Om){if(lamb_from_feed_om_L[jj]==kTRUE || lamb_from_feed_aom_L[jj]==kTRUE){GoodTrackFourVectorT_Lam_FD_L.push_back(pvectorM);}}
if(DoFeedDown_Xi){if(lamb_from_feed_xi_AL[jj]==kTRUE || lamb_from_feed_axi_AL[jj]==kTRUE){GoodTrackFourVectorT_Lam_FD_AL.push_back(pvectorM);}}
if(DoFeedDown_Om){if(lamb_from_feed_om_AL[jj]==kTRUE || lamb_from_feed_aom_AL[jj]==kTRUE){GoodTrackFourVectorT_Lam_FD_AL.push_back(pvectorM);}}

if(DoFeedDown_Xi){if(lamb_from_feed_xi_K[jj]==kTRUE || lamb_from_feed_axi_K[jj]==kTRUE)continue;}
if(DoFeedDown_Om){if(lamb_from_feed_om_K[jj]==kTRUE || lamb_from_feed_aom_K[jj]==kTRUE)continue;}
if(DoFeedDown_Xi){if(lamb_from_feed_xi_L[jj]==kTRUE || lamb_from_feed_axi_L[jj]==kTRUE)continue;}
if(DoFeedDown_Om){if(lamb_from_feed_om_L[jj]==kTRUE || lamb_from_feed_aom_L[jj]==kTRUE)continue;}
if(DoFeedDown_Xi){if(lamb_from_feed_xi_AL[jj]==kTRUE || lamb_from_feed_axi_AL[jj]==kTRUE)continue;}
if(DoFeedDown_Om){if(lamb_from_feed_om_AL[jj]==kTRUE || lamb_from_feed_aom_AL[jj]==kTRUE)continue;}

//LamLam
if(map_Lam_chi2[jj]==kTRUE && map_Lam_chi2_LAL[jj]==kTRUE && map_Lam_chi2_KL[jj]==kTRUE){pT_Lamb_daughter_prompt->Fill(pt);pT_Lamb_from_prompt_T->Fill(pt);}


if(isMC){

bool V0from2 = false;
Int_t njetmatch = 0;
Int_t NJet = -999;

for(Int_t ajj=0; ajj<jetsize; ajj++){ 
TLorentzVector pvectorj;
pvectorj.SetPtEtaPhiM(Jet_pt->at(ajj),Jet_eta->at(ajj),Jet_phi->at(ajj),Jet_mass->at(ajj));

double deltaPhi=pvectorM.DeltaPhi(pvectorj);
double deltaEta=pvectorM.Eta() -  pvectorj.Eta();
double DR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

double deltaPhid1=pvectorD1.DeltaPhi(pvectorj);
double deltaEtad1=pvectorD1.Eta() -  pvectorj.Eta();
double DRd1 = sqrt(deltaPhid1*deltaPhid1 + deltaEtad1*deltaEtad1);

double deltaPhid2=pvectorD2.DeltaPhi(pvectorj);
double deltaEtad2=pvectorD2.Eta() -  pvectorj.Eta();
double DRd2 = sqrt(deltaPhid2*deltaPhid2 + deltaEtad2*deltaEtad2);

cone_Lam->Fill(DR); cone_LamD1->Fill(DRd1); cone_LamD2->Fill(DRd2);
if(map_Lam_chi2[jj]==kTRUE){cone_LamT->Fill(DR); cone_LamTD1->Fill(DRd1); cone_LamTD2->Fill(DRd2);
if(DR <= 0.4){njetmatch=njetmatch+1; NJet = ajj; if(!removeV0from2jets){break;}}}
}

if(njetmatch>1){V0from2=true;}
if(V0from2)continue;

if(map_Lam_chi2[jj]==kTRUE){

if(njetmatch==1){GoodTrackFourVector_Lam_Jet.push_back(pvectorM); GoodTrackFourVector_Lam_NJet.push_back(NJet);}
if(njetmatch==0){GoodTrackFourVector_Lam_NotaJet.push_back(pvectorM);}

}

}


if(remove_all_mis){

if(map_Lam_chi2[jj]==kTRUE && map_Lam_chi2_LAL[jj]==kTRUE && map_Lam_chi2_KL[jj]==kTRUE){
GoodTrackFourVectorT_Lam.push_back(pvectorM);
GoodTrackFourVectord1T_Lam.push_back(pvectorD1);
GoodTrackFourVectord2T_Lam.push_back(pvectorD2);
chi21T_Lam.push_back(chi21X);
chi22T_Lam.push_back(chi22X);
GoodTrackFourVectorT_Lam_ALam.push_back(pvectorM);
GoodTrackFourVectord1T_Lam_ALam.push_back(pvectorD1);
GoodTrackFourVectord2T_Lam_ALam.push_back(pvectorD2);
chi21T_Lam_ALam.push_back(chi21X);
chi22T_Lam_ALam.push_back(chi22X);
GoodTrackFourVectorT_Lam_K0s.push_back(pvectorM);
GoodTrackFourVectord1T_Lam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2T_Lam_K0s.push_back(pvectorD2);
chi21T_Lam_K0s.push_back(chi21X);
chi22T_Lam_K0s.push_back(chi22X);
if(map_Lam_match[jj]==kTRUE){
GoodTrackFourVectorTM_Lam.push_back(pvectorM);
GoodTrackFourVectord1TM_Lam.push_back(pvectorD1);
GoodTrackFourVectord2TM_Lam.push_back(pvectorD2);
chi21TM_Lam.push_back(chi21X);
chi22TM_Lam.push_back(chi22X);
GoodTrackFourVectorTM_Lam_ALam.push_back(pvectorM);
GoodTrackFourVectord1TM_Lam_ALam.push_back(pvectorD1);
GoodTrackFourVectord2TM_Lam_ALam.push_back(pvectorD2);
chi21TM_Lam_ALam.push_back(chi21X);
chi22TM_Lam_ALam.push_back(chi22X);
GoodTrackFourVectorTM_Lam_K0s.push_back(pvectorM);
GoodTrackFourVectord1TM_Lam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2TM_Lam_K0s.push_back(pvectorD2);
chi21TM_Lam_K0s.push_back(chi21X);
chi22TM_Lam_K0s.push_back(chi22X);
}else{
GoodTrackFourVectorTU_Lam.push_back(pvectorM);
GoodTrackFourVectord1TU_Lam.push_back(pvectorD1);
GoodTrackFourVectord2TU_Lam.push_back(pvectorD2);
chi21TU_Lam.push_back(chi21X);
chi22TU_Lam.push_back(chi22X);
GoodTrackFourVectorTU_Lam_ALam.push_back(pvectorM);
GoodTrackFourVectord1TU_Lam_ALam.push_back(pvectorD1);
GoodTrackFourVectord2TU_Lam_ALam.push_back(pvectorD2);
chi21TU_Lam_ALam.push_back(chi21X);
chi22TU_Lam_ALam.push_back(chi22X);
GoodTrackFourVectorTU_Lam_K0s.push_back(pvectorM);
GoodTrackFourVectord1TU_Lam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2TU_Lam_K0s.push_back(pvectorD2);
chi21TU_Lam_K0s.push_back(chi21X);
chi22TU_Lam_K0s.push_back(chi22X);
}}else{
GoodTrackFourVectorF_Lam.push_back(pvectorM);
GoodTrackFourVectord1F_Lam.push_back(pvectorD1);
GoodTrackFourVectord2F_Lam.push_back(pvectorD2);
chi21F_Lam.push_back(chi21X);
chi22F_Lam.push_back(chi22X);
GoodTrackFourVectorF_Lam_ALam.push_back(pvectorM);
GoodTrackFourVectord1F_Lam_ALam.push_back(pvectorD1);
GoodTrackFourVectord2F_Lam_ALam.push_back(pvectorD2);
chi21F_Lam_ALam.push_back(chi21X);
chi22F_Lam_ALam.push_back(chi22X);
GoodTrackFourVectorF_Lam_K0s.push_back(pvectorM);
GoodTrackFourVectord1F_Lam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2F_Lam_K0s.push_back(pvectorD2);
chi21F_Lam_K0s.push_back(chi21X);
chi22F_Lam_K0s.push_back(chi22X);
if(map_Lam_match[jj]==kTRUE){
GoodTrackFourVectorFM_Lam.push_back(pvectorM);
GoodTrackFourVectord1FM_Lam.push_back(pvectorD1);
GoodTrackFourVectord2FM_Lam.push_back(pvectorD2);
chi21FM_Lam.push_back(chi21X);
chi22FM_Lam.push_back(chi22X);
GoodTrackFourVectorFM_Lam_ALam.push_back(pvectorM);
GoodTrackFourVectord1FM_Lam_ALam.push_back(pvectorD1);
GoodTrackFourVectord2FM_Lam_ALam.push_back(pvectorD2);
chi21FM_Lam_ALam.push_back(chi21X);
chi22FM_Lam_ALam.push_back(chi22X);
GoodTrackFourVectorFM_Lam_K0s.push_back(pvectorM);
GoodTrackFourVectord1FM_Lam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2FM_Lam_K0s.push_back(pvectorD2);
chi21FM_Lam_K0s.push_back(chi21X);
chi22FM_Lam_K0s.push_back(chi22X);
}else{
GoodTrackFourVectorFU_Lam.push_back(pvectorM);
GoodTrackFourVectord1FU_Lam.push_back(pvectorD1);
GoodTrackFourVectord2FU_Lam.push_back(pvectorD2);
chi21FU_Lam.push_back(chi21X);
chi22FU_Lam.push_back(chi22X);
GoodTrackFourVectorFU_Lam_ALam.push_back(pvectorM);
GoodTrackFourVectord1FU_Lam_ALam.push_back(pvectorD1);
GoodTrackFourVectord2FU_Lam_ALam.push_back(pvectorD2);
chi21FU_Lam_ALam.push_back(chi21X);
chi22FU_Lam_ALam.push_back(chi22X);
GoodTrackFourVectorFU_Lam_K0s.push_back(pvectorM);
GoodTrackFourVectord1FU_Lam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2FU_Lam_K0s.push_back(pvectorD2);
chi21FU_Lam_K0s.push_back(chi21X);
chi22FU_Lam_K0s.push_back(chi22X);
}}

}else{

if(map_Lam_chi2[jj]==kTRUE){
GoodTrackFourVectorT_Lam.push_back(pvectorM);
GoodTrackFourVectord1T_Lam.push_back(pvectorD1);
GoodTrackFourVectord2T_Lam.push_back(pvectorD2);
chi21T_Lam.push_back(chi21X);
chi22T_Lam.push_back(chi22X);
if(map_Lam_match[jj]==kTRUE){
GoodTrackFourVectorTM_Lam.push_back(pvectorM);
GoodTrackFourVectord1TM_Lam.push_back(pvectorD1);
GoodTrackFourVectord2TM_Lam.push_back(pvectorD2);
chi21TM_Lam.push_back(chi21X);
chi22TM_Lam.push_back(chi22X);
}else{
GoodTrackFourVectorTU_Lam.push_back(pvectorM);
GoodTrackFourVectord1TU_Lam.push_back(pvectorD1);
GoodTrackFourVectord2TU_Lam.push_back(pvectorD2);
chi21TU_Lam.push_back(chi21X);
chi22TU_Lam.push_back(chi22X);
}}else{
GoodTrackFourVectorF_Lam.push_back(pvectorM);
GoodTrackFourVectord1F_Lam.push_back(pvectorD1);
GoodTrackFourVectord2F_Lam.push_back(pvectorD2);
chi21F_Lam.push_back(chi21X);
chi22F_Lam.push_back(chi22X);
if(map_Lam_match[jj]==kTRUE){
GoodTrackFourVectorFM_Lam.push_back(pvectorM);
GoodTrackFourVectord1FM_Lam.push_back(pvectorD1);
GoodTrackFourVectord2FM_Lam.push_back(pvectorD2);
chi21FM_Lam.push_back(chi21X);
chi22FM_Lam.push_back(chi22X);
}else{
GoodTrackFourVectorFU_Lam.push_back(pvectorM);
GoodTrackFourVectord1FU_Lam.push_back(pvectorD1);
GoodTrackFourVectord2FU_Lam.push_back(pvectorD2);
chi21FU_Lam.push_back(chi21X);
chi22FU_Lam.push_back(chi22X);
}}

//LamALam
if(map_Lam_chi2_LAL[jj]==kTRUE){
GoodTrackFourVectorT_Lam_ALam.push_back(pvectorM);
GoodTrackFourVectord1T_Lam_ALam.push_back(pvectorD1);
GoodTrackFourVectord2T_Lam_ALam.push_back(pvectorD2);
chi21T_Lam_ALam.push_back(chi21X);
chi22T_Lam_ALam.push_back(chi22X);
if(map_Lam_match[jj]==kTRUE){
GoodTrackFourVectorTM_Lam_ALam.push_back(pvectorM);
GoodTrackFourVectord1TM_Lam_ALam.push_back(pvectorD1);
GoodTrackFourVectord2TM_Lam_ALam.push_back(pvectorD2);
chi21TM_Lam_ALam.push_back(chi21X);
chi22TM_Lam_ALam.push_back(chi22X);
}else{
GoodTrackFourVectorTU_Lam_ALam.push_back(pvectorM);
GoodTrackFourVectord1TU_Lam_ALam.push_back(pvectorD1);
GoodTrackFourVectord2TU_Lam_ALam.push_back(pvectorD2);
chi21TU_Lam_ALam.push_back(chi21X);
chi22TU_Lam_ALam.push_back(chi22X);
}}else{
GoodTrackFourVectorF_Lam_ALam.push_back(pvectorM);
GoodTrackFourVectord1F_Lam_ALam.push_back(pvectorD1);
GoodTrackFourVectord2F_Lam_ALam.push_back(pvectorD2);
chi21F_Lam_ALam.push_back(chi21X);
chi22F_Lam_ALam.push_back(chi22X);
if(map_Lam_match[jj]==kTRUE){
GoodTrackFourVectorFM_Lam_ALam.push_back(pvectorM);
GoodTrackFourVectord1FM_Lam_ALam.push_back(pvectorD1);
GoodTrackFourVectord2FM_Lam_ALam.push_back(pvectorD2);
chi21FM_Lam_ALam.push_back(chi21X);
chi22FM_Lam_ALam.push_back(chi22X);
}else{
GoodTrackFourVectorFU_Lam_ALam.push_back(pvectorM);
GoodTrackFourVectord1FU_Lam_ALam.push_back(pvectorD1);
GoodTrackFourVectord2FU_Lam_ALam.push_back(pvectorD2);
chi21FU_Lam_ALam.push_back(chi21X);
chi22FU_Lam_ALam.push_back(chi22X);
}}

//LamK0s
if(map_Lam_chi2_KL[jj]==kTRUE){
GoodTrackFourVectorT_Lam_K0s.push_back(pvectorM);
GoodTrackFourVectord1T_Lam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2T_Lam_K0s.push_back(pvectorD2);
chi21T_Lam_K0s.push_back(chi21X);
chi22T_Lam_K0s.push_back(chi22X);
if(map_Lam_match[jj]==kTRUE){
GoodTrackFourVectorTM_Lam_K0s.push_back(pvectorM);
GoodTrackFourVectord1TM_Lam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2TM_Lam_K0s.push_back(pvectorD2);
chi21TM_Lam_K0s.push_back(chi21X);
chi22TM_Lam_K0s.push_back(chi22X);
}else{
GoodTrackFourVectorTU_Lam_K0s.push_back(pvectorM);
GoodTrackFourVectord1TU_Lam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2TU_Lam_K0s.push_back(pvectorD2);
chi21TU_Lam_K0s.push_back(chi21X);
chi22TU_Lam_K0s.push_back(chi22X);
}}else{
GoodTrackFourVectorF_Lam_K0s.push_back(pvectorM);
GoodTrackFourVectord1F_Lam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2F_Lam_K0s.push_back(pvectorD2);
chi21F_Lam_K0s.push_back(chi21X);
chi22F_Lam_K0s.push_back(chi22X);
if(map_Lam_match[jj]==kTRUE){
GoodTrackFourVectorFM_Lam_K0s.push_back(pvectorM);
GoodTrackFourVectord1FM_Lam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2FM_Lam_K0s.push_back(pvectorD2);
chi21FM_Lam_K0s.push_back(chi21X);
chi22FM_Lam_K0s.push_back(chi22X);
}else{
GoodTrackFourVectorFU_Lam_K0s.push_back(pvectorM);
GoodTrackFourVectord1FU_Lam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2FU_Lam_K0s.push_back(pvectorD2);
chi21FU_Lam_K0s.push_back(chi21X);
chi22FU_Lam_K0s.push_back(chi22X);
}}

}

}

//--------------------true--------------------

if(GoodTrackFourVectorT_Lam.size()>=2){

nev_LamT->Fill(1); 

for(unsigned int iy=0; iy<GoodTrackFourVectorT_Lam.size();iy++){Lam_Mass->Fill(GoodTrackFourVectorT_Lam[iy].M(),Ntk_Vz_weight);}

Mass_dep(GoodTrackFourVectorT_Lam, (Double_t) aux_N_tk_offline, hMass_Lam_Lam_T,Ntk_Vz_weight);
Mass_dep_sep(GoodTrackFourVectorT_Lam, (Double_t) aux_N_tk_offline, hMass_Lam_1_T, hMass_Lam_2_T, Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorT_Lam,chi21T_Lam,chi22T_Lam,V0chi2d1_diff_T_Lam,V0chi2d2_diff_T_Lam,V0chi2d12_diff_T_Lam,Ntk_Vz_weight);
    
hbt_reco(GoodTrackFourVectorT_Lam, GoodTrackFourVectord1T_Lam, GoodTrackFourVectord2T_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_T, hpeakcos_Lam_Lam_T, hpeakd1_Lam_Lam_T, hpeakd2_Lam_Lam_T, hpeakd12_Lam_Lam_T, hpeak_rot_Lam_Lam_T, hpeak_inv_Lam_Lam_T, hside_Lam_Lam_T, hside_rot_Lam_Lam_T, hside_inv_Lam_Lam_T, hpeakside_Lam_Lam_T, hpeakside_rot_Lam_Lam_T, hpeakside_inv_Lam_Lam_T, hsideL_Lam_Lam_T, hsideR_Lam_Lam_T, hpeaksideL_Lam_Lam_T, hpeaksideR_Lam_Lam_T, p, effhisto_Lam, nev_Lam_ssT, nev_Lam_bbT, nev_Lam_sbT,Ntk_Vz_weight);

hdibaryon(GoodTrackFourVectorT_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, h_HDibaryon_Lam_Lam, h_HDibaryon_rot_Lam_Lam, h_HDibaryon_inv_Lam_Lam, h_HDibaryon_Lam_Lam_side, h_HDibaryon_rot_Lam_Lam_side, h_HDibaryon_inv_Lam_Lam_side, h_HDibaryon_Lam_Lam_peakside, h_HDibaryon_rot_Lam_Lam_peakside, h_HDibaryon_inv_Lam_Lam_peakside, Ntk_Vz_weight);

//weight
weight_LL.push_back(Ntk_Vz_weight);


//Mix Random

ev_z_vtxLam.push_back(aux_vtxz);
ev_ntrkoff_vecLam.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLam.push_back(GoodTrackFourVectorT_Lam); 

//Eta Mix

Double_t aux_etaMix_w = Eta_weight;
std::pair<Double_t, std::vector<TLorentzVector>> aux_pair_GoodTrackFourVector_etaMixWeight = make_pair(aux_etaMix_w, GoodTrackFourVectorT_Lam);
ev_GoodTrackFourVector_etaMixWeight_vecLam.push_back(aux_pair_GoodTrackFourVector_etaMixWeight);
std::pair<Double_t,int> aux_pair_ntrkoff_etaMixWeight = make_pair(aux_etaMix_w,aux_N_tk_offline);
ev_ntrkoff_etaMixWeight_vecLam.push_back(aux_pair_ntrkoff_etaMixWeight); 

}//LamLam

//--------------------false--------------------

if(GoodTrackFourVectorF_Lam.size()>=2){

nev_LamF->Fill(1); 

Mass_dep(GoodTrackFourVectorF_Lam, (Double_t) aux_N_tk_offline, hMass_Lam_Lam_F,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorF_Lam,chi21F_Lam,chi22F_Lam,V0chi2d1_diff_F_Lam,V0chi2d2_diff_F_Lam,V0chi2d12_diff_F_Lam,Ntk_Vz_weight);
    
hbt_reco_fake(GoodTrackFourVectorF_Lam, GoodTrackFourVectord1F_Lam, GoodTrackFourVectord2F_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_F, hpeakd1_Lam_Lam_F, hpeakd2_Lam_Lam_F, hpeakd12_Lam_Lam_F, hside_Lam_Lam_F, hpeakside_Lam_Lam_F, hpeak_Lam_Lam_F_bias, hpeak_Lam_Lam_F_biasd1, hpeak_Lam_Lam_F_biasd2, hpeak_Lam_Lam_F_biasd12,Ntk_Vz_weight);

}

//--------------------TF-------------------

if(GoodTrackFourVectorT_Lam.size()>=1 && GoodTrackFourVectorF_Lam.size()>=1){

nev_LamTF->Fill(1);

get_chi2plot_TF(GoodTrackFourVectorT_Lam, chi21T_Lam, chi22T_Lam, GoodTrackFourVectorF_Lam, chi21F_Lam, chi22F_Lam, V0chi2d1_diff_TF_Lam,V0chi2d2_diff_TF_Lam,V0chi2d12_diff_TF_Lam,Ntk_Vz_weight);

hbt_reco_truefake(GoodTrackFourVectorT_Lam, GoodTrackFourVectord1T_Lam, GoodTrackFourVectord2T_Lam, GoodTrackFourVectorF_Lam, GoodTrackFourVectord1F_Lam, GoodTrackFourVectord2F_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_TF, hpeakd1_Lam_Lam_TF, hpeakd2_Lam_Lam_TF, hpeakd12_Lam_Lam_TF, hside_Lam_Lam_TF, hpeakside_Lam_Lam_TF, hpeak_Lam_Lam_TF_bias, hpeak_Lam_Lam_TF_biasd1, hpeak_Lam_Lam_TF_biasd2, hpeak_Lam_Lam_TF_biasd12,Ntk_Vz_weight);
    
}

if(isMC){

//-------------------True Matched-----------------------

if(GoodTrackFourVectorTM_Lam.size()>=2){

Mass_dep(GoodTrackFourVectorTM_Lam, (Double_t) aux_N_tk_offline, hMass_Lam_Lam_TM,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorTM_Lam,chi21TM_Lam,chi22TM_Lam,V0chi2d1_diff_TM_Lam,V0chi2d2_diff_TM_Lam,V0chi2d12_diff_TM_Lam,Ntk_Vz_weight);
    
hbt_MM_UU(GoodTrackFourVectorTM_Lam, GoodTrackFourVectord1TM_Lam, GoodTrackFourVectord2TM_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_TM, hpeakd1_Lam_Lam_TM, hpeakd2_Lam_Lam_TM, hpeakd12_Lam_Lam_TM, hpeak_Lam_Lam_TM_bias, hpeak_Lam_Lam_TM_biasd1, hpeak_Lam_Lam_TM_biasd2, hpeak_Lam_Lam_TM_biasd12,Ntk_Vz_weight);

ev_z_vtxLam_TM.push_back(aux_vtxz);
ev_ntrkoff_vecLam_TM.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLam_TM.push_back(GoodTrackFourVectorTM_Lam); 


}

//-------------------True Unmatched-----------------------

if(GoodTrackFourVectorTU_Lam.size()>=2){

Mass_dep(GoodTrackFourVectorTU_Lam, (Double_t) aux_N_tk_offline, hMass_Lam_Lam_TU,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorTU_Lam,chi21TU_Lam,chi22TU_Lam,V0chi2d1_diff_TU_Lam,V0chi2d2_diff_TU_Lam,V0chi2d12_diff_TU_Lam,Ntk_Vz_weight);
    
hbt_MM_UU(GoodTrackFourVectorTU_Lam, GoodTrackFourVectord1TU_Lam, GoodTrackFourVectord2TU_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_TU, hpeakd1_Lam_Lam_TU, hpeakd2_Lam_Lam_TU, hpeakd12_Lam_Lam_TU, hpeak_Lam_Lam_TU_bias, hpeak_Lam_Lam_TU_biasd1, hpeak_Lam_Lam_TU_biasd2, hpeak_Lam_Lam_TU_biasd12,Ntk_Vz_weight);

ev_z_vtxLam_TU.push_back(aux_vtxz);
ev_ntrkoff_vecLam_TU.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLam_TU.push_back(GoodTrackFourVectorTU_Lam); 

}

//-------------------Fake Matched-----------------------

if(GoodTrackFourVectorFM_Lam.size()>=2){

Mass_dep(GoodTrackFourVectorFM_Lam, (Double_t) aux_N_tk_offline, hMass_Lam_Lam_FM,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorFM_Lam,chi21FM_Lam,chi22FM_Lam,V0chi2d1_diff_FM_Lam,V0chi2d2_diff_FM_Lam,V0chi2d12_diff_FM_Lam,Ntk_Vz_weight);
    
hbt_MM_UU(GoodTrackFourVectorFM_Lam, GoodTrackFourVectord1FM_Lam, GoodTrackFourVectord2FM_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_FM, hpeakd1_Lam_Lam_FM, hpeakd2_Lam_Lam_FM, hpeakd12_Lam_Lam_FM, hpeak_Lam_Lam_FM_bias, hpeak_Lam_Lam_FM_biasd1, hpeak_Lam_Lam_FM_biasd2, hpeak_Lam_Lam_FM_biasd12,Ntk_Vz_weight);

}

//-------------------Fake Unmatched-----------------------

if(GoodTrackFourVectorFU_Lam.size()>=2){

Mass_dep(GoodTrackFourVectorFU_Lam, (Double_t) aux_N_tk_offline, hMass_Lam_Lam_FU,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorFU_Lam,chi21FU_Lam,chi22FU_Lam,V0chi2d1_diff_FU_Lam,V0chi2d2_diff_FU_Lam,V0chi2d12_diff_FU_Lam,Ntk_Vz_weight);
    
hbt_MM_UU(GoodTrackFourVectorFU_Lam, GoodTrackFourVectord1FU_Lam, GoodTrackFourVectord2FU_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_FU, hpeakd1_Lam_Lam_FU, hpeakd2_Lam_Lam_FU, hpeakd12_Lam_Lam_FU, hpeak_Lam_Lam_FU_bias, hpeak_Lam_Lam_FU_biasd1, hpeak_Lam_Lam_FU_biasd2, hpeak_Lam_Lam_FU_biasd12,Ntk_Vz_weight);

}


//-------------------True Matched + True Unmatched-----------------------

if(GoodTrackFourVectorTM_Lam.size()>=1 && GoodTrackFourVectorTU_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_Lam, chi21TM_Lam, chi22TM_Lam, GoodTrackFourVectorTU_Lam, chi21TU_Lam, chi22TU_Lam, V0chi2d1_diff_TM_TU_Lam, V0chi2d2_diff_TM_TU_Lam, V0chi2d12_diff_TM_TU_Lam,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTM_Lam, GoodTrackFourVectord1TM_Lam, GoodTrackFourVectord2TM_Lam, GoodTrackFourVectorTU_Lam, GoodTrackFourVectord1TU_Lam, GoodTrackFourVectord2TU_Lam,(Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_TM_TU, hpeakd1_Lam_Lam_TM_TU, hpeakd2_Lam_Lam_TM_TU, hpeakd12_Lam_Lam_TM_TU, hpeak_Lam_Lam_TM_TU_bias, hpeak_Lam_Lam_TM_TU_biasd1, hpeak_Lam_Lam_TM_TU_biasd2, hpeak_Lam_Lam_TM_TU_biasd12,Ntk_Vz_weight);

ev_z_vtxLam_TM_TU.push_back(aux_vtxz);
ev_ntrkoff_vecLam_TM_TU.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLam_TM_TU.push_back(GoodTrackFourVectorTM_Lam); 
ev_GoodTrackFourVector_vecLam_TU_TM.push_back(GoodTrackFourVectorTU_Lam);

}

//-------------------True Matched + Fake Matched-----------------------

if(GoodTrackFourVectorTM_Lam.size()>=1 && GoodTrackFourVectorFM_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_Lam, chi21TM_Lam, chi22TM_Lam, GoodTrackFourVectorFM_Lam, chi21FM_Lam, chi22FM_Lam, V0chi2d1_diff_TM_FM_Lam, V0chi2d2_diff_TM_FM_Lam, V0chi2d12_diff_TM_FM_Lam,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTM_Lam, GoodTrackFourVectord1TM_Lam, GoodTrackFourVectord2TM_Lam, GoodTrackFourVectorFM_Lam, GoodTrackFourVectord1FM_Lam, GoodTrackFourVectord2FM_Lam,(Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  hpeak_Lam_Lam_TM_FM, hpeakd1_Lam_Lam_TM_FM, hpeakd2_Lam_Lam_TM_FM, hpeakd12_Lam_Lam_TM_FM, hpeak_Lam_Lam_TM_FM_bias, hpeak_Lam_Lam_TM_FM_biasd1, hpeak_Lam_Lam_TM_FM_biasd2, hpeak_Lam_Lam_TM_FM_biasd12,Ntk_Vz_weight);

}

//-------------------True Matched + Fake Unmatched-----------------------

if(GoodTrackFourVectorTM_Lam.size()>=1 && GoodTrackFourVectorFU_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_Lam, chi21TM_Lam, chi22TM_Lam, GoodTrackFourVectorFU_Lam, chi21FU_Lam, chi22FU_Lam, V0chi2d1_diff_TM_FU_Lam, V0chi2d2_diff_TM_FU_Lam, V0chi2d12_diff_TM_FU_Lam,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTM_Lam, GoodTrackFourVectord1TM_Lam, GoodTrackFourVectord2TM_Lam, GoodTrackFourVectorFU_Lam, GoodTrackFourVectord1FU_Lam, GoodTrackFourVectord2FU_Lam,(Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  hpeak_Lam_Lam_TM_FU, hpeakd1_Lam_Lam_TM_FU, hpeakd2_Lam_Lam_TM_FU, hpeakd12_Lam_Lam_TM_FU, hpeak_Lam_Lam_TM_FU_bias, hpeak_Lam_Lam_TM_FU_biasd1, hpeak_Lam_Lam_TM_FU_biasd2, hpeak_Lam_Lam_TM_FU_biasd12,Ntk_Vz_weight);

}


//-------------------True Unmatched + Fake Matched-----------------------

if(GoodTrackFourVectorTU_Lam.size()>=1 && GoodTrackFourVectorFM_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_Lam, chi21TU_Lam, chi22TU_Lam, GoodTrackFourVectorFM_Lam, chi21FM_Lam, chi22FM_Lam, V0chi2d1_diff_TU_FM_Lam, V0chi2d2_diff_TU_FM_Lam, V0chi2d12_diff_TU_FM_Lam,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTU_Lam, GoodTrackFourVectord1TU_Lam, GoodTrackFourVectord2TU_Lam, GoodTrackFourVectorFM_Lam, GoodTrackFourVectord1FM_Lam, GoodTrackFourVectord2FM_Lam,(Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  hpeak_Lam_Lam_TU_FM, hpeakd1_Lam_Lam_TU_FM, hpeakd2_Lam_Lam_TU_FM, hpeakd12_Lam_Lam_TU_FM, hpeak_Lam_Lam_TU_FM_bias, hpeak_Lam_Lam_TU_FM_biasd1, hpeak_Lam_Lam_TU_FM_biasd2, hpeak_Lam_Lam_TU_FM_biasd12,Ntk_Vz_weight);

}

//-------------------True Unmatched + Fake Unmatched-----------------------

if(GoodTrackFourVectorTU_Lam.size()>=1 && GoodTrackFourVectorFU_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_Lam, chi21TU_Lam, chi22TU_Lam, GoodTrackFourVectorFU_Lam, chi21FU_Lam, chi22FU_Lam, V0chi2d1_diff_TU_FU_Lam, V0chi2d2_diff_TU_FU_Lam, V0chi2d12_diff_TU_FU_Lam,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTU_Lam, GoodTrackFourVectord1TU_Lam, GoodTrackFourVectord2TU_Lam, GoodTrackFourVectorFU_Lam, GoodTrackFourVectord1FU_Lam, GoodTrackFourVectord2FU_Lam,(Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  hpeak_Lam_Lam_TU_FU, hpeakd1_Lam_Lam_TU_FU, hpeakd2_Lam_Lam_TU_FU, hpeakd12_Lam_Lam_TU_FU, hpeak_Lam_Lam_TU_FU_bias, hpeak_Lam_Lam_TU_FU_biasd1, hpeak_Lam_Lam_TU_FU_biasd2, hpeak_Lam_Lam_TU_FU_biasd12,Ntk_Vz_weight);

}


//-------------------Fake Matched + Fake Unmatched-----------------------

if(GoodTrackFourVectorFM_Lam.size()>=1 && GoodTrackFourVectorFU_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFM_Lam, chi21FM_Lam, chi22FM_Lam, GoodTrackFourVectorFU_Lam, chi21FU_Lam, chi22FU_Lam, V0chi2d1_diff_FM_FU_Lam, V0chi2d2_diff_FM_FU_Lam, V0chi2d12_diff_FM_FU_Lam,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorFM_Lam, GoodTrackFourVectord1FM_Lam, GoodTrackFourVectord2FM_Lam, GoodTrackFourVectorFU_Lam, GoodTrackFourVectord1FU_Lam, GoodTrackFourVectord2FU_Lam,(Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  hpeak_Lam_Lam_FM_FU, hpeakd1_Lam_Lam_FM_FU, hpeakd2_Lam_Lam_FM_FU, hpeakd12_Lam_Lam_FM_FU, hpeak_Lam_Lam_FM_FU_bias, hpeak_Lam_Lam_FM_FU_biasd1, hpeak_Lam_Lam_FM_FU_biasd2, hpeak_Lam_Lam_FM_FU_biasd12,Ntk_Vz_weight);

}

}


/*========================= AntiLambdaAntiLambda =========================*/

std::vector<TLorentzVector> GoodTrackFourVectord1T_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2T_ALam;
std::vector<TLorentzVector> GoodTrackFourVectorT_ALam;
std::vector<double> chi21T_ALam;
std::vector<double> chi22T_ALam;

std::vector<TLorentzVector> GoodTrackFourVectord1F_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2F_ALam;
std::vector<TLorentzVector> GoodTrackFourVectorF_ALam;
std::vector<double> chi21F_ALam;
std::vector<double> chi22F_ALam;

//4-vectors for true matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTM_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1TM_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2TM_ALam;
std::vector<double> chi21TM_ALam;
std::vector<double> chi22TM_ALam;

//4-vectors for true unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTU_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1TU_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2TU_ALam;
std::vector<double> chi21TU_ALam;
std::vector<double> chi22TU_ALam;

//4-vectors for fake matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFM_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1FM_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2FM_ALam;
std::vector<double> chi21FM_ALam;
std::vector<double> chi22FM_ALam;

//4-vectors for fake unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFU_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord1FU_ALam;
std::vector<TLorentzVector> GoodTrackFourVectord2FU_ALam;
std::vector<double> chi21FU_ALam;
std::vector<double> chi22FU_ALam;

/*========================= AntiLambdaLambda =========================*/

//4-vectors for true particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectord1T_ALam_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2T_ALam_Lam;
std::vector<TLorentzVector> GoodTrackFourVectorT_ALam_Lam;
std::vector<double> chi21T_ALam_Lam;
std::vector<double> chi22T_ALam_Lam;
//4-vectors for fake particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectord1F_ALam_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2F_ALam_Lam;
std::vector<TLorentzVector> GoodTrackFourVectorF_ALam_Lam;
std::vector<double> chi21F_ALam_Lam;
std::vector<double> chi22F_ALam_Lam;
//4-vectors for true matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTM_ALam_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1TM_ALam_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2TM_ALam_Lam;
std::vector<double> chi21TM_ALam_Lam;
std::vector<double> chi22TM_ALam_Lam;
//4-vectors for true unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTU_ALam_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1TU_ALam_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2TU_ALam_Lam;
std::vector<double> chi21TU_ALam_Lam;
std::vector<double> chi22TU_ALam_Lam;
//4-vectors for fake matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFM_ALam_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1FM_ALam_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2FM_ALam_Lam;
std::vector<double> chi21FM_ALam_Lam;
std::vector<double> chi22FM_ALam_Lam;
//4-vectors for fake unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFU_ALam_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord1FU_ALam_Lam;
std::vector<TLorentzVector> GoodTrackFourVectord2FU_ALam_Lam;
std::vector<double> chi21FU_ALam_Lam;
std::vector<double> chi22FU_ALam_Lam;

/*========================= AntiLambdaK0s =========================*/

//4-vectors for true particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectord1T_ALam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2T_ALam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectorT_ALam_K0s;
std::vector<double> chi21T_ALam_K0s;
std::vector<double> chi22T_ALam_K0s;
//4-vectors for fake particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectord1F_ALam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2F_ALam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectorF_ALam_K0s;
std::vector<double> chi21F_ALam_K0s;
std::vector<double> chi22F_ALam_K0s;
//4-vectors for true matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTM_ALam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1TM_ALam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2TM_ALam_K0s;
std::vector<double> chi21TM_ALam_K0s;
std::vector<double> chi22TM_ALam_K0s;
//4-vectors for true unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorTU_ALam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1TU_ALam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2TU_ALam_K0s;
std::vector<double> chi21TU_ALam_K0s;
std::vector<double> chi22TU_ALam_K0s;
//4-vectors for fake matched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFM_ALam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1FM_ALam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2FM_ALam_K0s;
std::vector<double> chi21FM_ALam_K0s;
std::vector<double> chi22FM_ALam_K0s;
//4-vectors for fake unmatched particles after chi2 cut
std::vector<TLorentzVector> GoodTrackFourVectorFU_ALam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord1FU_ALam_K0s;
std::vector<TLorentzVector> GoodTrackFourVectord2FU_ALam_K0s;
std::vector<double> chi21FU_ALam_K0s;
std::vector<double> chi22FU_ALam_K0s;

//for feed down

std::vector<TLorentzVector> GoodTrackFourVectorT_ALam_FD_K;
std::vector<TLorentzVector> GoodTrackFourVectorT_ALam_FD_L;
std::vector<TLorentzVector> GoodTrackFourVectorT_ALam_FD_AL;


for(Int_t jj=0; jj<GoodTrackFourVector_ALam.size(); jj++){

double pt = GoodTrackFourVector_ALam[jj].Pt();
double eta = GoodTrackFourVector_ALam[jj].Eta();
double phi = GoodTrackFourVector_ALam[jj].Phi();
double m = GoodTrackFourVector_ALam[jj].M();

double ptd1 = GoodTrackFourVectord1_ALam[jj].Pt();
double etad1 = GoodTrackFourVectord1_ALam[jj].Eta();
double phid1 = GoodTrackFourVectord1_ALam[jj].Phi();
double md1 = GoodTrackFourVectord1_ALam[jj].M();

double ptd2 = GoodTrackFourVectord2_ALam[jj].Pt();
double etad2 = GoodTrackFourVectord2_ALam[jj].Eta();
double phid2 = GoodTrackFourVectord2_ALam[jj].Phi();
double md2 = GoodTrackFourVectord2_ALam[jj].M();

double chi21X = chi21_ALam[jj];
double chi22X = chi22_ALam[jj];

TLorentzVector pvectorM;
pvectorM.SetPtEtaPhiM(pt,eta,phi,m);
TLorentzVector pvectorD1;
pvectorD1.SetPtEtaPhiM(ptd1,etad1,phid1,md1);
TLorentzVector pvectorD2;
pvectorD2.SetPtEtaPhiM(ptd2,etad2,phid2,md2);


//feed down
if(DoFeedDown_Xi){if(alamb_from_feed_xi_K[jj]==kTRUE || alamb_from_feed_axi_K[jj]==kTRUE || alamb_from_feed_xi_L[jj]==kTRUE || alamb_from_feed_axi_L[jj]==kTRUE || alamb_from_feed_xi_AL[jj]==kTRUE || alamb_from_feed_axi_AL[jj]==kTRUE){pT_Lamb_daughter_Xi->Fill(pt);}}
if(DoFeedDown_Om){if(alamb_from_feed_om_K[jj]==kTRUE || alamb_from_feed_aom_K[jj]==kTRUE || alamb_from_feed_om_L[jj]==kTRUE || alamb_from_feed_aom_L[jj]==kTRUE || alamb_from_feed_om_AL[jj]==kTRUE || alamb_from_feed_aom_AL[jj]==kTRUE){pT_Lamb_daughter_Om->Fill(pt);}}
if(DoFeedDown_Xi){if(alamb_from_feed_xi_K[jj]==kTRUE || alamb_from_feed_axi_K[jj]==kTRUE){GoodTrackFourVectorT_ALam_FD_K.push_back(pvectorM);}}
if(DoFeedDown_Om){if(alamb_from_feed_om_K[jj]==kTRUE || alamb_from_feed_aom_K[jj]==kTRUE){GoodTrackFourVectorT_ALam_FD_K.push_back(pvectorM);}}
if(DoFeedDown_Xi){if(alamb_from_feed_xi_L[jj]==kTRUE || alamb_from_feed_axi_L[jj]==kTRUE){GoodTrackFourVectorT_ALam_FD_L.push_back(pvectorM);}}
if(DoFeedDown_Om){if(alamb_from_feed_om_L[jj]==kTRUE || alamb_from_feed_aom_L[jj]==kTRUE){GoodTrackFourVectorT_ALam_FD_L.push_back(pvectorM);}}
if(DoFeedDown_Xi){if(alamb_from_feed_xi_AL[jj]==kTRUE || alamb_from_feed_axi_AL[jj]==kTRUE){GoodTrackFourVectorT_ALam_FD_AL.push_back(pvectorM);}}
if(DoFeedDown_Om){if(alamb_from_feed_om_AL[jj]==kTRUE || alamb_from_feed_aom_AL[jj]==kTRUE){GoodTrackFourVectorT_ALam_FD_AL.push_back(pvectorM);}}

if(DoFeedDown_Xi){if(alamb_from_feed_xi_K[jj]==kTRUE || alamb_from_feed_axi_K[jj]==kTRUE)continue;}
if(DoFeedDown_Om){if(alamb_from_feed_om_K[jj]==kTRUE || alamb_from_feed_aom_K[jj]==kTRUE)continue;}
if(DoFeedDown_Xi){if(alamb_from_feed_xi_L[jj]==kTRUE || alamb_from_feed_axi_L[jj]==kTRUE)continue;}
if(DoFeedDown_Om){if(alamb_from_feed_om_L[jj]==kTRUE || alamb_from_feed_aom_L[jj]==kTRUE)continue;}
if(DoFeedDown_Xi){if(alamb_from_feed_xi_AL[jj]==kTRUE || alamb_from_feed_axi_AL[jj]==kTRUE)continue;}
if(DoFeedDown_Om){if(alamb_from_feed_om_AL[jj]==kTRUE || alamb_from_feed_aom_AL[jj]==kTRUE)continue;}


//ALamALam
if(map_ALam_chi2[jj]==kTRUE && map_ALam_chi2_LAL[jj]==kTRUE && map_ALam_chi2_KAL[jj]==kTRUE){pT_Lamb_daughter_prompt->Fill(pt);pT_Lamb_from_prompt_T->Fill(pt);}


if(isMC){

bool V0from2 = false;
Int_t njetmatch = 0;
Int_t NJet = -999;

for(Int_t ajj=0; ajj<jetsize; ajj++){ 
TLorentzVector pvectorj;
pvectorj.SetPtEtaPhiM(Jet_pt->at(ajj),Jet_eta->at(ajj),Jet_phi->at(ajj),Jet_mass->at(ajj));

double deltaPhi=pvectorM.DeltaPhi(pvectorj);
double deltaEta=pvectorM.Eta() -  pvectorj.Eta();
double DR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

double deltaPhid1=pvectorD1.DeltaPhi(pvectorj);
double deltaEtad1=pvectorD1.Eta() -  pvectorj.Eta();
double DRd1 = sqrt(deltaPhid1*deltaPhid1 + deltaEtad1*deltaEtad1);

double deltaPhid2=pvectorD2.DeltaPhi(pvectorj);
double deltaEtad2=pvectorD2.Eta() -  pvectorj.Eta();
double DRd2 = sqrt(deltaPhid2*deltaPhid2 + deltaEtad2*deltaEtad2);

cone_Lam->Fill(DR); cone_LamD1->Fill(DRd1); cone_LamD2->Fill(DRd2);
if(map_ALam_chi2[jj]==kTRUE){cone_LamT->Fill(DR); cone_LamTD1->Fill(DRd1); cone_LamTD2->Fill(DRd2);
if(DR <= 0.4){njetmatch=njetmatch+1; NJet = ajj; if(!removeV0from2jets){break;}}}
}

if(njetmatch>1){V0from2=true;}
if(V0from2)continue;

if(map_ALam_chi2[jj]==kTRUE){

if(njetmatch==1){GoodTrackFourVector_ALam_Jet.push_back(pvectorM); GoodTrackFourVector_ALam_NJet.push_back(NJet);}
if(njetmatch==0){GoodTrackFourVector_ALam_NotaJet.push_back(pvectorM);}

}

}


if(remove_all_mis){

if(map_ALam_chi2[jj]==kTRUE){
GoodTrackFourVectorT_ALam.push_back(pvectorM);
GoodTrackFourVectord1T_ALam.push_back(pvectorD1);
GoodTrackFourVectord2T_ALam.push_back(pvectorD2);
chi21T_ALam.push_back(chi21X);
chi22T_ALam.push_back(chi22X);
GoodTrackFourVectorT_ALam_Lam.push_back(pvectorM);
GoodTrackFourVectord1T_ALam_Lam.push_back(pvectorD1);
GoodTrackFourVectord2T_ALam_Lam.push_back(pvectorD2);
chi21T_ALam_Lam.push_back(chi21X);
chi22T_ALam_Lam.push_back(chi22X);
GoodTrackFourVectorT_ALam_K0s.push_back(pvectorM);
GoodTrackFourVectord1T_ALam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2T_ALam_K0s.push_back(pvectorD2);
chi21T_ALam_K0s.push_back(chi21X);
chi22T_ALam_K0s.push_back(chi22X);
if(map_ALam_match[jj]==kTRUE){
GoodTrackFourVectorTM_ALam.push_back(pvectorM);
GoodTrackFourVectord1TM_ALam.push_back(pvectorD1);
GoodTrackFourVectord2TM_ALam.push_back(pvectorD2);
chi21TM_ALam.push_back(chi21X);
chi22TM_ALam.push_back(chi22X);
GoodTrackFourVectorTM_ALam_Lam.push_back(pvectorM);
GoodTrackFourVectord1TM_ALam_Lam.push_back(pvectorD1);
GoodTrackFourVectord2TM_ALam_Lam.push_back(pvectorD2);
chi21TM_ALam_Lam.push_back(chi21X);
chi22TM_ALam_Lam.push_back(chi22X);
GoodTrackFourVectorTM_ALam_K0s.push_back(pvectorM);
GoodTrackFourVectord1TM_ALam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2TM_ALam_K0s.push_back(pvectorD2);
chi21TM_ALam_K0s.push_back(chi21X);
chi22TM_ALam_K0s.push_back(chi22X);
}else{
GoodTrackFourVectorTU_ALam.push_back(pvectorM);
GoodTrackFourVectord1TU_ALam.push_back(pvectorD1);
GoodTrackFourVectord2TU_ALam.push_back(pvectorD2);
chi21TU_ALam.push_back(chi21X);
chi22TU_ALam.push_back(chi22X);
GoodTrackFourVectorTU_ALam_Lam.push_back(pvectorM);
GoodTrackFourVectord1TU_ALam_Lam.push_back(pvectorD1);
GoodTrackFourVectord2TU_ALam_Lam.push_back(pvectorD2);
chi21TU_ALam_Lam.push_back(chi21X);
chi22TU_ALam_Lam.push_back(chi22X);
GoodTrackFourVectorTU_ALam_K0s.push_back(pvectorM);
GoodTrackFourVectord1TU_ALam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2TU_ALam_K0s.push_back(pvectorD2);
chi21TU_ALam_K0s.push_back(chi21X);
chi22TU_ALam_K0s.push_back(chi22X);
}}else{
GoodTrackFourVectorF_ALam.push_back(pvectorM);
GoodTrackFourVectord1F_ALam.push_back(pvectorD1);
GoodTrackFourVectord2F_ALam.push_back(pvectorD2);
chi21F_ALam.push_back(chi21X);
chi22F_ALam.push_back(chi22X);
GoodTrackFourVectorF_ALam_Lam.push_back(pvectorM);
GoodTrackFourVectord1F_ALam_Lam.push_back(pvectorD1);
GoodTrackFourVectord2F_ALam_Lam.push_back(pvectorD2);
chi21F_ALam_Lam.push_back(chi21X);
chi22F_ALam_Lam.push_back(chi22X);
GoodTrackFourVectorF_ALam_K0s.push_back(pvectorM);
GoodTrackFourVectord1F_ALam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2F_ALam_K0s.push_back(pvectorD2);
chi21F_ALam_K0s.push_back(chi21X);
chi22F_ALam_K0s.push_back(chi22X);
if(map_ALam_match[jj]==kTRUE){
GoodTrackFourVectorFM_ALam.push_back(pvectorM);
GoodTrackFourVectord1FM_ALam.push_back(pvectorD1);
GoodTrackFourVectord2FM_ALam.push_back(pvectorD2);
chi21FM_ALam.push_back(chi21X);
chi22FM_ALam.push_back(chi22X);
GoodTrackFourVectorFM_ALam_Lam.push_back(pvectorM);
GoodTrackFourVectord1FM_ALam_Lam.push_back(pvectorD1);
GoodTrackFourVectord2FM_ALam_Lam.push_back(pvectorD2);
chi21FM_ALam_Lam.push_back(chi21X);
chi22FM_ALam_Lam.push_back(chi22X);
GoodTrackFourVectorFM_ALam_K0s.push_back(pvectorM);
GoodTrackFourVectord1FM_ALam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2FM_ALam_K0s.push_back(pvectorD2);
chi21FM_ALam_K0s.push_back(chi21X);
chi22FM_ALam_K0s.push_back(chi22X);
}else{
GoodTrackFourVectorFU_ALam.push_back(pvectorM);
GoodTrackFourVectord1FU_ALam.push_back(pvectorD1);
GoodTrackFourVectord2FU_ALam.push_back(pvectorD2);
chi21FU_ALam.push_back(chi21X);
chi22FU_ALam.push_back(chi22X);
GoodTrackFourVectorFU_ALam_Lam.push_back(pvectorM);
GoodTrackFourVectord1FU_ALam_Lam.push_back(pvectorD1);
GoodTrackFourVectord2FU_ALam_Lam.push_back(pvectorD2);
chi21FU_ALam_Lam.push_back(chi21X);
chi22FU_ALam_Lam.push_back(chi22X);
GoodTrackFourVectorFU_ALam_K0s.push_back(pvectorM);
GoodTrackFourVectord1FU_ALam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2FU_ALam_K0s.push_back(pvectorD2);
chi21FU_ALam_K0s.push_back(chi21X);
chi22FU_ALam_K0s.push_back(chi22X);
}}

}else{

if(map_ALam_chi2[jj]==kTRUE){
GoodTrackFourVectorT_ALam.push_back(pvectorM);
GoodTrackFourVectord1T_ALam.push_back(pvectorD1);
GoodTrackFourVectord2T_ALam.push_back(pvectorD2);
chi21T_ALam.push_back(chi21X);
chi22T_ALam.push_back(chi22X);
if(map_ALam_match[jj]==kTRUE){
GoodTrackFourVectorTM_ALam.push_back(pvectorM);
GoodTrackFourVectord1TM_ALam.push_back(pvectorD1);
GoodTrackFourVectord2TM_ALam.push_back(pvectorD2);
chi21TM_ALam.push_back(chi21X);
chi22TM_ALam.push_back(chi22X);
}else{
GoodTrackFourVectorTU_ALam.push_back(pvectorM);
GoodTrackFourVectord1TU_ALam.push_back(pvectorD1);
GoodTrackFourVectord2TU_ALam.push_back(pvectorD2);
chi21TU_ALam.push_back(chi21X);
chi22TU_ALam.push_back(chi22X);
}}else{
GoodTrackFourVectorF_ALam.push_back(pvectorM);
GoodTrackFourVectord1F_ALam.push_back(pvectorD1);
GoodTrackFourVectord2F_ALam.push_back(pvectorD2);
chi21F_ALam.push_back(chi21X);
chi22F_ALam.push_back(chi22X);
if(map_ALam_match[jj]==kTRUE){
GoodTrackFourVectorFM_ALam.push_back(pvectorM);
GoodTrackFourVectord1FM_ALam.push_back(pvectorD1);
GoodTrackFourVectord2FM_ALam.push_back(pvectorD2);
chi21FM_ALam.push_back(chi21X);
chi22FM_ALam.push_back(chi22X);
}else{
GoodTrackFourVectorFU_ALam.push_back(pvectorM);
GoodTrackFourVectord1FU_ALam.push_back(pvectorD1);
GoodTrackFourVectord2FU_ALam.push_back(pvectorD2);
chi21FU_ALam.push_back(chi21X);
chi22FU_ALam.push_back(chi22X);
}}

//LamALam
if(map_ALam_chi2_LAL[jj]==kTRUE){
GoodTrackFourVectorT_ALam_Lam.push_back(pvectorM);
GoodTrackFourVectord1T_ALam_Lam.push_back(pvectorD1);
GoodTrackFourVectord2T_ALam_Lam.push_back(pvectorD2);
chi21T_ALam_Lam.push_back(chi21X);
chi22T_ALam_Lam.push_back(chi22X);
if(map_ALam_match[jj]==kTRUE){
GoodTrackFourVectorTM_ALam_Lam.push_back(pvectorM);
GoodTrackFourVectord1TM_ALam_Lam.push_back(pvectorD1);
GoodTrackFourVectord2TM_ALam_Lam.push_back(pvectorD2);
chi21TM_ALam_Lam.push_back(chi21X);
chi22TM_ALam_Lam.push_back(chi22X);
}else{
GoodTrackFourVectorTU_ALam_Lam.push_back(pvectorM);
GoodTrackFourVectord1TU_ALam_Lam.push_back(pvectorD1);
GoodTrackFourVectord2TU_ALam_Lam.push_back(pvectorD2);
chi21TU_ALam_Lam.push_back(chi21X);
chi22TU_ALam_Lam.push_back(chi22X);
}}else{
GoodTrackFourVectorF_ALam_Lam.push_back(pvectorM);
GoodTrackFourVectord1F_ALam_Lam.push_back(pvectorD1);
GoodTrackFourVectord2F_ALam_Lam.push_back(pvectorD2);
chi21F_ALam_Lam.push_back(chi21X);
chi22F_ALam_Lam.push_back(chi22X);
if(map_ALam_match[jj]==kTRUE){
GoodTrackFourVectorFM_ALam_Lam.push_back(pvectorM);
GoodTrackFourVectord1FM_ALam_Lam.push_back(pvectorD1);
GoodTrackFourVectord2FM_ALam_Lam.push_back(pvectorD2);
chi21FM_ALam_Lam.push_back(chi21X);
chi22FM_ALam_Lam.push_back(chi22X);
}else{
GoodTrackFourVectorFU_ALam_Lam.push_back(pvectorM);
GoodTrackFourVectord1FU_ALam_Lam.push_back(pvectorD1);
GoodTrackFourVectord2FU_ALam_Lam.push_back(pvectorD2);
chi21FU_ALam_Lam.push_back(chi21X);
chi22FU_ALam_Lam.push_back(chi22X);
}}

//LamK0s
if(map_ALam_chi2_KAL[jj]==kTRUE){
GoodTrackFourVectorT_ALam_K0s.push_back(pvectorM);
GoodTrackFourVectord1T_ALam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2T_ALam_K0s.push_back(pvectorD2);
chi21T_ALam_K0s.push_back(chi21X);
chi22T_ALam_K0s.push_back(chi22X);
if(map_ALam_match[jj]==kTRUE){
GoodTrackFourVectorTM_ALam_K0s.push_back(pvectorM);
GoodTrackFourVectord1TM_ALam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2TM_ALam_K0s.push_back(pvectorD2);
chi21TM_ALam_K0s.push_back(chi21X);
chi22TM_ALam_K0s.push_back(chi22X);
}else{
GoodTrackFourVectorTU_ALam_K0s.push_back(pvectorM);
GoodTrackFourVectord1TU_ALam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2TU_ALam_K0s.push_back(pvectorD2);
chi21TU_ALam_K0s.push_back(chi21X);
chi22TU_ALam_K0s.push_back(chi22X);
}}else{
GoodTrackFourVectorF_ALam_K0s.push_back(pvectorM);
GoodTrackFourVectord1F_ALam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2F_ALam_K0s.push_back(pvectorD2);
chi21F_ALam_K0s.push_back(chi21X);
chi22F_ALam_K0s.push_back(chi22X);
if(map_ALam_match[jj]==kTRUE){
GoodTrackFourVectorFM_ALam_K0s.push_back(pvectorM);
GoodTrackFourVectord1FM_ALam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2FM_ALam_K0s.push_back(pvectorD2);
chi21FM_ALam_K0s.push_back(chi21X);
chi22FM_ALam_K0s.push_back(chi22X);
}else{
GoodTrackFourVectorFU_ALam_K0s.push_back(pvectorM);
GoodTrackFourVectord1FU_ALam_K0s.push_back(pvectorD1);
GoodTrackFourVectord2FU_ALam_K0s.push_back(pvectorD2);
chi21FU_ALam_K0s.push_back(chi21X);
chi22FU_ALam_K0s.push_back(chi22X);
}}

}

}

//--------------------true--------------------

if(GoodTrackFourVectorT_ALam.size()>=2){

nev_ALamT->Fill(1); 

for(unsigned int iy=0; iy<GoodTrackFourVectorT_ALam.size();iy++){ALam_Mass->Fill(GoodTrackFourVectorT_ALam[iy].M(),Ntk_Vz_weight);}

Mass_dep(GoodTrackFourVectorT_ALam, (Double_t) aux_N_tk_offline, hMass_ALam_ALam_T,Ntk_Vz_weight);
Mass_dep_sep(GoodTrackFourVectorT_ALam, (Double_t) aux_N_tk_offline, hMass_ALam_1_T, hMass_ALam_2_T, Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorT_ALam,chi21T_ALam,chi22T_ALam,V0chi2d1_diff_T_ALam,V0chi2d2_diff_T_ALam,V0chi2d12_diff_T_ALam,Ntk_Vz_weight);
    
hbt_reco(GoodTrackFourVectorT_ALam, GoodTrackFourVectord1T_ALam, GoodTrackFourVectord2T_ALam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_T, hpeakcos_ALam_ALam_T, hpeakd1_ALam_ALam_T, hpeakd2_ALam_ALam_T, hpeakd12_ALam_ALam_T, hpeak_rot_ALam_ALam_T, hpeak_inv_ALam_ALam_T, hside_ALam_ALam_T, hside_rot_ALam_ALam_T, hside_inv_ALam_ALam_T, hpeakside_ALam_ALam_T, hpeakside_rot_ALam_ALam_T, hpeakside_inv_ALam_ALam_T, hsideL_ALam_ALam_T, hsideR_ALam_ALam_T, hpeaksideL_ALam_ALam_T, hpeaksideR_ALam_ALam_T, p, effhisto_ALam, nev_ALam_ssT, nev_ALam_bbT, nev_ALam_sbT,Ntk_Vz_weight);

hdibaryon(GoodTrackFourVectorT_ALam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, h_HDibaryon_ALam_ALam, h_HDibaryon_rot_ALam_ALam, h_HDibaryon_inv_ALam_ALam, h_HDibaryon_ALam_ALam_side, h_HDibaryon_rot_ALam_ALam_side, h_HDibaryon_inv_ALam_ALam_side, h_HDibaryon_ALam_ALam_peakside, h_HDibaryon_rot_ALam_ALam_peakside, h_HDibaryon_inv_ALam_ALam_peakside, Ntk_Vz_weight);

//weight
weight_ALAL.push_back(Ntk_Vz_weight);


//Mix Random

ev_z_vtxALam.push_back(aux_vtxz);
ev_ntrkoff_vecALam.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecALam.push_back(GoodTrackFourVectorT_ALam); 

//Eta Mix

Double_t aux_etaMix_w = Eta_weight;
std::pair<Double_t, std::vector<TLorentzVector>> aux_pair_GoodTrackFourVector_etaMixWeight = make_pair(aux_etaMix_w, GoodTrackFourVectorT_ALam);
ev_GoodTrackFourVector_etaMixWeight_vecALam.push_back(aux_pair_GoodTrackFourVector_etaMixWeight);
std::pair<Double_t,int> aux_pair_ntrkoff_etaMixWeight = make_pair(aux_etaMix_w,aux_N_tk_offline);
ev_ntrkoff_etaMixWeight_vecALam.push_back(aux_pair_ntrkoff_etaMixWeight); 

}//ALamALam

//--------------------false--------------------

if(GoodTrackFourVectorF_ALam.size()>=2){

nev_ALamF->Fill(1); 

Mass_dep(GoodTrackFourVectorF_ALam, (Double_t) aux_N_tk_offline, hMass_ALam_ALam_F,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorF_ALam,chi21F_ALam,chi22F_ALam,V0chi2d1_diff_F_ALam,V0chi2d2_diff_F_ALam,V0chi2d12_diff_F_ALam,Ntk_Vz_weight);
    
hbt_reco_fake(GoodTrackFourVectorF_ALam, GoodTrackFourVectord1F_ALam, GoodTrackFourVectord2F_ALam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_F, hpeakd1_ALam_ALam_F, hpeakd2_ALam_ALam_F, hpeakd12_ALam_ALam_F, hside_ALam_ALam_F, hpeakside_ALam_ALam_F, hpeak_ALam_ALam_F_bias, hpeak_ALam_ALam_F_biasd1, hpeak_ALam_ALam_F_biasd2, hpeak_ALam_ALam_F_biasd12,Ntk_Vz_weight);

}

//--------------------TF-------------------

if(GoodTrackFourVectorT_ALam.size()>=1 && GoodTrackFourVectorF_ALam.size()>=1){

nev_ALamTF->Fill(1);

get_chi2plot_TF(GoodTrackFourVectorT_ALam, chi21T_ALam, chi22T_ALam, GoodTrackFourVectorF_ALam, chi21F_ALam, chi22F_ALam, V0chi2d1_diff_TF_ALam,V0chi2d2_diff_TF_ALam,V0chi2d12_diff_TF_ALam,Ntk_Vz_weight);

hbt_reco_truefake(GoodTrackFourVectorT_ALam, GoodTrackFourVectord1T_ALam, GoodTrackFourVectord2T_ALam, GoodTrackFourVectorF_ALam, GoodTrackFourVectord1F_ALam, GoodTrackFourVectord2F_ALam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_TF, hpeakd1_ALam_ALam_TF, hpeakd2_ALam_ALam_TF, hpeakd12_ALam_ALam_TF, hside_ALam_ALam_TF, hpeakside_ALam_ALam_TF, hpeak_ALam_ALam_TF_bias, hpeak_ALam_ALam_TF_biasd1, hpeak_ALam_ALam_TF_biasd2, hpeak_ALam_ALam_TF_biasd12,Ntk_Vz_weight);
    
}

if(isMC){

//-------------------True Matched-----------------------

if(GoodTrackFourVectorTM_ALam.size()>=2){

Mass_dep(GoodTrackFourVectorTM_ALam, (Double_t) aux_N_tk_offline, hMass_ALam_ALam_TM, Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorTM_ALam,chi21TM_ALam,chi22TM_ALam,V0chi2d1_diff_TM_ALam,V0chi2d2_diff_TM_ALam,V0chi2d12_diff_TM_ALam,Ntk_Vz_weight);
    
hbt_MM_UU(GoodTrackFourVectorTM_ALam, GoodTrackFourVectord1TM_ALam, GoodTrackFourVectord2TM_ALam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_TM, hpeakd1_ALam_ALam_TM, hpeakd2_ALam_ALam_TM, hpeakd12_ALam_ALam_TM, hpeak_ALam_ALam_TM_bias, hpeak_ALam_ALam_TM_biasd1, hpeak_ALam_ALam_TM_biasd2, hpeak_ALam_ALam_TM_biasd12,Ntk_Vz_weight);

ev_z_vtxALam_TM.push_back(aux_vtxz);
ev_ntrkoff_vecALam_TM.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecALam_TM.push_back(GoodTrackFourVectorTM_ALam); 


}

//-------------------True Unmatched-----------------------

if(GoodTrackFourVectorTU_ALam.size()>=2){

Mass_dep(GoodTrackFourVectorTU_ALam, (Double_t) aux_N_tk_offline, hMass_ALam_ALam_TU,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorTU_ALam,chi21TU_ALam,chi22TU_ALam,V0chi2d1_diff_TU_ALam,V0chi2d2_diff_TU_ALam,V0chi2d12_diff_TU_ALam,Ntk_Vz_weight);
    
hbt_MM_UU(GoodTrackFourVectorTU_ALam, GoodTrackFourVectord1TU_ALam, GoodTrackFourVectord2TU_ALam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_TU, hpeakd1_ALam_ALam_TU, hpeakd2_ALam_ALam_TU, hpeakd12_ALam_ALam_TU, hpeak_ALam_ALam_TU_bias, hpeak_ALam_ALam_TU_biasd1, hpeak_ALam_ALam_TU_biasd2, hpeak_ALam_ALam_TU_biasd12,Ntk_Vz_weight);

ev_z_vtxALam_TU.push_back(aux_vtxz);
ev_ntrkoff_vecALam_TU.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecALam_TU.push_back(GoodTrackFourVectorTU_ALam); 


}

//-------------------Fake Matched-----------------------

if(GoodTrackFourVectorFM_ALam.size()>=2){

Mass_dep(GoodTrackFourVectorFM_ALam, (Double_t) aux_N_tk_offline, hMass_ALam_ALam_FM,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorFM_ALam,chi21FM_ALam,chi22FM_ALam,V0chi2d1_diff_FM_ALam,V0chi2d2_diff_FM_ALam,V0chi2d12_diff_FM_ALam,Ntk_Vz_weight);
    
hbt_MM_UU(GoodTrackFourVectorFM_ALam, GoodTrackFourVectord1FM_ALam, GoodTrackFourVectord2FM_ALam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_FM, hpeakd1_ALam_ALam_FM, hpeakd2_ALam_ALam_FM, hpeakd12_ALam_ALam_FM, hpeak_ALam_ALam_FM_bias, hpeak_ALam_ALam_FM_biasd1, hpeak_ALam_ALam_FM_biasd2, hpeak_ALam_ALam_FM_biasd12,Ntk_Vz_weight);

}

//-------------------Fake Unmatched-----------------------

if(GoodTrackFourVectorFU_ALam.size()>=2){

Mass_dep(GoodTrackFourVectorFU_ALam, (Double_t) aux_N_tk_offline, hMass_ALam_ALam_FU,Ntk_Vz_weight);

get_chi2plot(GoodTrackFourVectorFU_ALam,chi21FU_ALam,chi22FU_ALam,V0chi2d1_diff_FU_ALam,V0chi2d2_diff_FU_ALam,V0chi2d12_diff_FU_ALam,Ntk_Vz_weight);
    
hbt_MM_UU(GoodTrackFourVectorFU_ALam, GoodTrackFourVectord1FU_ALam, GoodTrackFourVectord2FU_ALam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_FU, hpeakd1_ALam_ALam_FU, hpeakd2_ALam_ALam_FU, hpeakd12_ALam_ALam_FU, hpeak_ALam_ALam_FU_bias, hpeak_ALam_ALam_FU_biasd1, hpeak_ALam_ALam_FU_biasd2, hpeak_ALam_ALam_FU_biasd12,Ntk_Vz_weight);

}


//-------------------True Matched + True Unmatched-----------------------

if(GoodTrackFourVectorTM_ALam.size()>=1 && GoodTrackFourVectorTU_ALam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_ALam, chi21TM_ALam, chi22TM_ALam, GoodTrackFourVectorTU_ALam, chi21TU_ALam, chi22TU_ALam, V0chi2d1_diff_TM_TU_ALam, V0chi2d2_diff_TM_TU_ALam, V0chi2d12_diff_TM_TU_ALam,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTM_ALam, GoodTrackFourVectord1TM_ALam, GoodTrackFourVectord2TM_ALam, GoodTrackFourVectorTU_ALam, GoodTrackFourVectord1TU_ALam, GoodTrackFourVectord2TU_ALam,(Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_TM_TU, hpeakd1_ALam_ALam_TM_TU, hpeakd2_ALam_ALam_TM_TU, hpeakd12_ALam_ALam_TM_TU, hpeak_ALam_ALam_TM_TU_bias, hpeak_ALam_ALam_TM_TU_biasd1, hpeak_ALam_ALam_TM_TU_biasd2, hpeak_ALam_ALam_TM_TU_biasd12,Ntk_Vz_weight);

ev_z_vtxALam_TM_TU.push_back(aux_vtxz);
ev_ntrkoff_vecALam_TM_TU.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecALam_TM_TU.push_back(GoodTrackFourVectorTM_ALam); 
ev_GoodTrackFourVector_vecALam_TU_TM.push_back(GoodTrackFourVectorTU_ALam);

}

//-------------------True Matched + Fake Matched-----------------------

if(GoodTrackFourVectorTM_ALam.size()>=1 && GoodTrackFourVectorFM_ALam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_ALam, chi21TM_ALam, chi22TM_ALam, GoodTrackFourVectorFM_ALam, chi21FM_ALam, chi22FM_ALam, V0chi2d1_diff_TM_FM_ALam, V0chi2d2_diff_TM_FM_ALam, V0chi2d12_diff_TM_FM_ALam,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTM_ALam, GoodTrackFourVectord1TM_ALam, GoodTrackFourVectord2TM_ALam, GoodTrackFourVectorFM_ALam, GoodTrackFourVectord1FM_ALam, GoodTrackFourVectord2FM_ALam,(Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_TM_FM, hpeakd1_ALam_ALam_TM_FM, hpeakd2_ALam_ALam_TM_FM, hpeakd12_ALam_ALam_TM_FM, hpeak_ALam_ALam_TM_FM_bias, hpeak_ALam_ALam_TM_FM_biasd1, hpeak_ALam_ALam_TM_FM_biasd2, hpeak_ALam_ALam_TM_FM_biasd12,Ntk_Vz_weight);

}

//-------------------True Matched + Fake Unmatched-----------------------

if(GoodTrackFourVectorTM_ALam.size()>=1 && GoodTrackFourVectorFU_ALam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_ALam, chi21TM_ALam, chi22TM_ALam, GoodTrackFourVectorFU_ALam, chi21FU_ALam, chi22FU_ALam, V0chi2d1_diff_TM_FU_ALam, V0chi2d2_diff_TM_FU_ALam, V0chi2d12_diff_TM_FU_ALam,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTM_ALam, GoodTrackFourVectord1TM_ALam, GoodTrackFourVectord2TM_ALam, GoodTrackFourVectorFU_ALam, GoodTrackFourVectord1FU_ALam, GoodTrackFourVectord2FU_ALam,(Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  hpeak_ALam_ALam_TM_FU, hpeakd1_ALam_ALam_TM_FU, hpeakd2_ALam_ALam_TM_FU, hpeakd12_ALam_ALam_TM_FU, hpeak_ALam_ALam_TM_FU_bias, hpeak_ALam_ALam_TM_FU_biasd1, hpeak_ALam_ALam_TM_FU_biasd2, hpeak_ALam_ALam_TM_FU_biasd12,Ntk_Vz_weight);

}


//-------------------True Unmatched + Fake Matched-----------------------

if(GoodTrackFourVectorTU_ALam.size()>=1 && GoodTrackFourVectorFM_ALam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_ALam, chi21TU_ALam, chi22TU_ALam, GoodTrackFourVectorFM_ALam, chi21FM_ALam, chi22FM_ALam, V0chi2d1_diff_TU_FM_ALam, V0chi2d2_diff_TU_FM_ALam, V0chi2d12_diff_TU_FM_ALam,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTU_ALam, GoodTrackFourVectord1TU_ALam, GoodTrackFourVectord2TU_ALam, GoodTrackFourVectorFM_ALam, GoodTrackFourVectord1FM_ALam, GoodTrackFourVectord2FM_ALam,(Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  hpeak_ALam_ALam_TU_FM, hpeakd1_ALam_ALam_TU_FM, hpeakd2_ALam_ALam_TU_FM, hpeakd12_ALam_ALam_TU_FM, hpeak_ALam_ALam_TU_FM_bias, hpeak_ALam_ALam_TU_FM_biasd1, hpeak_ALam_ALam_TU_FM_biasd2, hpeak_ALam_ALam_TU_FM_biasd12,Ntk_Vz_weight);

}

//-------------------True Unmatched + Fake Unmatched-----------------------

if(GoodTrackFourVectorTU_ALam.size()>=1 && GoodTrackFourVectorFU_ALam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_ALam, chi21TU_ALam, chi22TU_ALam, GoodTrackFourVectorFU_ALam, chi21FU_ALam, chi22FU_ALam, V0chi2d1_diff_TU_FU_ALam, V0chi2d2_diff_TU_FU_ALam, V0chi2d12_diff_TU_FU_ALam,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorTU_ALam, GoodTrackFourVectord1TU_ALam, GoodTrackFourVectord2TU_ALam, GoodTrackFourVectorFU_ALam, GoodTrackFourVectord1FU_ALam, GoodTrackFourVectord2FU_ALam,(Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  hpeak_ALam_ALam_TU_FU, hpeakd1_ALam_ALam_TU_FU, hpeakd2_ALam_ALam_TU_FU, hpeakd12_ALam_ALam_TU_FU, hpeak_ALam_ALam_TU_FU_bias, hpeak_ALam_ALam_TU_FU_biasd1, hpeak_ALam_ALam_TU_FU_biasd2, hpeak_ALam_ALam_TU_FU_biasd12,Ntk_Vz_weight);

}


//-------------------Fake Matched + Fake Unmatched-----------------------

if(GoodTrackFourVectorFM_ALam.size()>=1 && GoodTrackFourVectorFU_ALam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFM_ALam, chi21FM_ALam, chi22FM_ALam, GoodTrackFourVectorFU_ALam, chi21FU_ALam, chi22FU_ALam, V0chi2d1_diff_FM_FU_ALam, V0chi2d2_diff_FM_FU_ALam, V0chi2d12_diff_FM_FU_ALam,Ntk_Vz_weight);

hbt_MM_UU_cross(GoodTrackFourVectorFM_ALam, GoodTrackFourVectord1FM_ALam, GoodTrackFourVectord2FM_ALam, GoodTrackFourVectorFU_ALam, GoodTrackFourVectord1FU_ALam, GoodTrackFourVectord2FU_ALam,(Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_FM_FU, hpeakd1_ALam_ALam_FM_FU, hpeakd2_ALam_ALam_FM_FU, hpeakd12_ALam_ALam_FM_FU, hpeak_ALam_ALam_FM_FU_bias, hpeak_ALam_ALam_FM_FU_biasd1, hpeak_ALam_ALam_FM_FU_biasd2, hpeak_ALam_ALam_FM_FU_biasd12,Ntk_Vz_weight);

}

}

/*========================= LambdaAntiLambda =========================*/

if(GoodTrackFourVectorT_Lam_ALam.size()>=1 && GoodTrackFourVectorT_ALam_Lam.size()>=1){

nev_LALT->Fill(1);

for(unsigned int iy=0; iy<GoodTrackFourVectorT_Lam_ALam.size();iy++){LAL_Mass->Fill(GoodTrackFourVectorT_Lam_ALam[iy].M(),Ntk_Vz_weight);}

for(unsigned int iy=0; iy<GoodTrackFourVectorT_ALam_Lam.size();iy++){LAL_Mass->Fill(GoodTrackFourVectorT_ALam_Lam[iy].M(),Ntk_Vz_weight);}

Mass_dep_OS(GoodTrackFourVectorT_Lam_ALam,GoodTrackFourVectorT_ALam_Lam,(Double_t) aux_N_tk_offline,hMass_Lam_T,hMass_ALam_T,hMass_Lam_ALam_T,Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorT_Lam_ALam,chi21T_Lam_ALam,chi22T_Lam_ALam,GoodTrackFourVectorT_ALam_Lam,chi21T_ALam_Lam,chi22T_ALam_Lam,V0chi2d1_diff_T_Lam_ALam,V0chi2d2_diff_T_Lam_ALam,V0chi2d12_diff_T_Lam_ALam,Ntk_Vz_weight);

hbt_reco_cross(GoodTrackFourVectorT_Lam_ALam, GoodTrackFourVectord1T_Lam_ALam, GoodTrackFourVectord2T_Lam_ALam, GoodTrackFourVectorT_ALam_Lam, GoodTrackFourVectord1T_ALam_Lam, GoodTrackFourVectord2T_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_T, hpeakcos_Lam_ALam_T, hpeakd1_Lam_ALam_T, hpeakd2_Lam_ALam_T, hpeakd12_Lam_ALam_T, hpeak_rot_Lam_ALam_T, hpeak_inv_Lam_ALam_T, hside_Lam_ALam_T, hside_rot_Lam_ALam_T, hside_inv_Lam_ALam_T, hpeakside_Lam_ALam_T, hpeakside_rot_Lam_ALam_T, hpeakside_inv_Lam_ALam_T, hsideL_Lam_ALam_T, hsideR_Lam_ALam_T, hpeaksideL_Lam_ALam_T, hpeaksideR_Lam_ALam_T, p, effhisto_Lam, effhisto_ALam, nev_LAL_ssT, nev_LAL_bbT, nev_LAL_sbT,Ntk_Vz_weight);

hdibaryon_OS(GoodTrackFourVectorT_Lam_ALam, GoodTrackFourVectorT_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, nsigmaside, h_HDibaryon_Lam_ALam, h_HDibaryon_Lam_ALam_side, h_HDibaryon_Lam_ALam_peakside, Ntk_Vz_weight);

//weight
weight_LAL.push_back(Ntk_Vz_weight);


ev_z_vtxLamALam.push_back(aux_vtxz);
ev_ntrkoff_vecLamALam.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLamALam.push_back(GoodTrackFourVectorT_Lam_ALam); 
ev_GoodTrackFourVector_vecALamLam.push_back(GoodTrackFourVectorT_ALam_Lam); 

Double_t aux_etaMix_w = Eta_weight;
std::pair<Double_t, std::vector<TLorentzVector>> aux_pair_GoodTrackFourVector_etaMixWeight = make_pair(aux_etaMix_w, GoodTrackFourVectorT_Lam_ALam);
ev_GoodTrackFourVector_etaMixWeight_vecLamALam.push_back(aux_pair_GoodTrackFourVector_etaMixWeight);
std::pair<Double_t, std::vector<TLorentzVector>> aux_pair_GoodTrackFourVector_etaMixWeight2 = make_pair(aux_etaMix_w, GoodTrackFourVectorT_ALam_Lam);
ev_GoodTrackFourVector_etaMixWeight_vecALamLam.push_back(aux_pair_GoodTrackFourVector_etaMixWeight2);
std::pair<Double_t,int> aux_pair_ntrkoff_etaMixWeight = make_pair(aux_etaMix_w,aux_N_tk_offline);
ev_ntrkoff_etaMixWeight_vecLamALam.push_back(aux_pair_ntrkoff_etaMixWeight); 


}

if(GoodTrackFourVectorF_Lam_ALam.size()>=1 && GoodTrackFourVectorF_ALam_Lam.size()>=1){

nev_LALF->Fill(1);

Mass_dep_OS(GoodTrackFourVectorF_Lam_ALam,GoodTrackFourVectorF_ALam_Lam,(Double_t) aux_N_tk_offline,hMass_Lam_F,hMass_ALam_F,hMass_Lam_ALam_F,Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorF_Lam_ALam, chi21F_Lam_ALam, chi22F_Lam_ALam, GoodTrackFourVectorF_ALam_Lam, chi21F_ALam_Lam, chi22F_ALam_Lam, V0chi2d1_diff_F_Lam_ALam, V0chi2d2_diff_F_Lam_ALam, V0chi2d12_diff_F_Lam_ALam, Ntk_Vz_weight);

hbt_reco_truefake(GoodTrackFourVectorF_Lam_ALam, GoodTrackFourVectord1F_Lam_ALam, GoodTrackFourVectord2F_Lam_ALam, GoodTrackFourVectorF_ALam_Lam, GoodTrackFourVectord1F_ALam_Lam, GoodTrackFourVectord2F_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_F, hpeakd1_Lam_ALam_F, hpeakd2_Lam_ALam_F, hpeakd12_Lam_ALam_F, hside_Lam_ALam_F, hpeakside_Lam_ALam_F, hpeak_Lam_ALam_F_bias, hpeak_Lam_ALam_F_biasd1, hpeak_Lam_ALam_F_biasd2, hpeak_Lam_ALam_F_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorF_Lam_ALam.size()>=1 && GoodTrackFourVectorT_ALam_Lam.size()>=1){

nev_LALTF->Fill(1);

get_chi2plot_TF(GoodTrackFourVectorF_Lam_ALam, chi21F_Lam_ALam, chi22F_Lam_ALam, GoodTrackFourVectorT_ALam_Lam, chi21T_ALam_Lam, chi22T_ALam_Lam, V0chi2d1_diff_TF_Lam_ALam, V0chi2d2_diff_TF_Lam_ALam, V0chi2d12_diff_TF_Lam_ALam, Ntk_Vz_weight);

hbt_reco_truefake(GoodTrackFourVectorF_Lam_ALam, GoodTrackFourVectord1F_Lam_ALam, GoodTrackFourVectord2F_Lam_ALam, GoodTrackFourVectorT_ALam_Lam, GoodTrackFourVectord1T_ALam_Lam, GoodTrackFourVectord2T_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TF, hpeakd1_Lam_ALam_TF, hpeakd2_Lam_ALam_TF, hpeakd12_Lam_ALam_TF, hside_Lam_ALam_TF, hpeakside_Lam_ALam_TF, hpeak_Lam_ALam_TF_bias, hpeak_Lam_ALam_TF_biasd1, hpeak_Lam_ALam_TF_biasd2, hpeak_Lam_ALam_TF_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorT_Lam_ALam.size()>=1 && GoodTrackFourVectorF_ALam_Lam.size()>=1){

nev_LALTF->Fill(1);

get_chi2plot_TF(GoodTrackFourVectorT_Lam_ALam, chi21T_Lam_ALam, chi22T_Lam_ALam, GoodTrackFourVectorF_ALam_Lam, chi21F_ALam_Lam, chi22F_ALam_Lam, V0chi2d1_diff_TF_Lam_ALam, V0chi2d2_diff_TF_Lam_ALam, V0chi2d12_diff_TF_Lam_ALam, Ntk_Vz_weight);

hbt_reco_truefake(GoodTrackFourVectorT_Lam_ALam, GoodTrackFourVectord1T_Lam_ALam, GoodTrackFourVectord2T_Lam_ALam, GoodTrackFourVectorF_ALam_Lam, GoodTrackFourVectord1F_ALam_Lam, GoodTrackFourVectord2F_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TF, hpeakd1_Lam_ALam_TF, hpeakd2_Lam_ALam_TF, hpeakd12_Lam_ALam_TF, hside_Lam_ALam_TF, hpeakside_Lam_ALam_TF, hpeak_Lam_ALam_TF_bias, hpeak_Lam_ALam_TF_biasd1, hpeak_Lam_ALam_TF_biasd2, hpeak_Lam_ALam_TF_biasd12, Ntk_Vz_weight);
    
}

if(isMC){

//-------------------True Matched-----------------------

if(GoodTrackFourVectorTM_Lam_ALam.size()>=1 && GoodTrackFourVectorTM_ALam_Lam.size()>=1){

Mass_dep_OS(GoodTrackFourVectorTM_Lam_ALam,GoodTrackFourVectorTM_ALam_Lam,(Double_t) aux_N_tk_offline,hMass_Lam_TM,hMass_ALam_TM,hMass_Lam_ALam_TM, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorTM_Lam_ALam,chi21TM_Lam_ALam,chi22TM_Lam_ALam,GoodTrackFourVectorTM_ALam_Lam,chi21TM_ALam_Lam,chi22TM_ALam_Lam,V0chi2d1_diff_TM_Lam_ALam,V0chi2d2_diff_TM_Lam_ALam,V0chi2d12_diff_TM_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTM_Lam_ALam, GoodTrackFourVectord1TM_Lam_ALam, GoodTrackFourVectord2TM_Lam_ALam, GoodTrackFourVectorTM_ALam_Lam, GoodTrackFourVectord1TM_ALam_Lam, GoodTrackFourVectord2TM_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TM, hpeakd1_Lam_ALam_TM, hpeakd2_Lam_ALam_TM, hpeakd12_Lam_ALam_TM, hpeak_Lam_ALam_TM_bias, hpeak_Lam_ALam_TM_biasd1, hpeak_Lam_ALam_TM_biasd2, hpeak_Lam_ALam_TM_biasd12, Ntk_Vz_weight);

ev_z_vtxLAL_TM.push_back(aux_vtxz);
ev_ntrkoff_vecLAL_TM.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLAL_TM.push_back(GoodTrackFourVectorTM_Lam_ALam); 
ev_GoodTrackFourVector_vecALL_TM.push_back(GoodTrackFourVectorTM_ALam_Lam); 

}

//-------------------True Unmatched-----------------------

if(GoodTrackFourVectorTU_Lam_ALam.size()>=1 && GoodTrackFourVectorTU_ALam_Lam.size()>=1){

Mass_dep_OS(GoodTrackFourVectorTU_Lam_ALam,GoodTrackFourVectorTU_ALam_Lam,(Double_t) aux_N_tk_offline,hMass_Lam_TU,hMass_ALam_TU,hMass_Lam_ALam_TU, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorTU_Lam_ALam,chi21TU_Lam_ALam,chi22TU_Lam_ALam,GoodTrackFourVectorTU_ALam_Lam,chi21TU_ALam_Lam,chi22TU_ALam_Lam,V0chi2d1_diff_TU_Lam_ALam,V0chi2d2_diff_TU_Lam_ALam,V0chi2d12_diff_TU_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTU_Lam_ALam, GoodTrackFourVectord1TU_Lam_ALam, GoodTrackFourVectord2TU_Lam_ALam, GoodTrackFourVectorTU_ALam_Lam, GoodTrackFourVectord1TU_ALam_Lam, GoodTrackFourVectord2TU_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TU, hpeakd1_Lam_ALam_TU, hpeakd2_Lam_ALam_TU, hpeakd12_Lam_ALam_TU, hpeak_Lam_ALam_TU_bias, hpeak_Lam_ALam_TU_biasd1, hpeak_Lam_ALam_TU_biasd2, hpeak_Lam_ALam_TU_biasd12, Ntk_Vz_weight);

ev_z_vtxLAL_TU.push_back(aux_vtxz);
ev_ntrkoff_vecLAL_TU.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLAL_TU.push_back(GoodTrackFourVectorTU_Lam_ALam); 
ev_GoodTrackFourVector_vecALL_TU.push_back(GoodTrackFourVectorTU_ALam_Lam); 


}

//-------------------Fake Matched-----------------------

if(GoodTrackFourVectorFM_Lam_ALam.size()>=1 && GoodTrackFourVectorFM_ALam_Lam.size()>=1){

Mass_dep_OS(GoodTrackFourVectorFM_Lam_ALam,GoodTrackFourVectorFM_ALam_Lam,(Double_t) aux_N_tk_offline,hMass_Lam_FM,hMass_ALam_FM,hMass_Lam_ALam_FM, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorFM_Lam_ALam,chi21FM_Lam_ALam,chi22FM_Lam_ALam,GoodTrackFourVectorFM_ALam_Lam,chi21FM_ALam_Lam,chi22FM_ALam_Lam,V0chi2d1_diff_FM_Lam_ALam,V0chi2d2_diff_FM_Lam_ALam,V0chi2d12_diff_FM_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFM_Lam_ALam, GoodTrackFourVectord1FM_Lam_ALam, GoodTrackFourVectord2FM_Lam_ALam, GoodTrackFourVectorFM_ALam_Lam, GoodTrackFourVectord1FM_ALam_Lam, GoodTrackFourVectord2FM_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_FM, hpeakd1_Lam_ALam_FM, hpeakd2_Lam_ALam_FM, hpeakd12_Lam_ALam_FM, hpeak_Lam_ALam_FM_bias, hpeak_Lam_ALam_FM_biasd1, hpeak_Lam_ALam_FM_biasd2, hpeak_Lam_ALam_FM_biasd12, Ntk_Vz_weight);


}

//-------------------True Unmatched-----------------------

if(GoodTrackFourVectorFU_Lam_ALam.size()>=1 && GoodTrackFourVectorFU_ALam_Lam.size()>=1){

Mass_dep_OS(GoodTrackFourVectorFU_Lam_ALam,GoodTrackFourVectorFU_ALam_Lam,(Double_t) aux_N_tk_offline,hMass_Lam_FU,hMass_ALam_FU,hMass_Lam_ALam_FU, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorFU_Lam_ALam,chi21FU_Lam_ALam,chi22FU_Lam_ALam,GoodTrackFourVectorFU_ALam_Lam,chi21FU_ALam_Lam,chi22FU_ALam_Lam,V0chi2d1_diff_FU_Lam_ALam,V0chi2d2_diff_FU_Lam_ALam,V0chi2d12_diff_FU_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFU_Lam_ALam, GoodTrackFourVectord1FU_Lam_ALam, GoodTrackFourVectord2FU_Lam_ALam, GoodTrackFourVectorFU_ALam_Lam, GoodTrackFourVectord1FU_ALam_Lam, GoodTrackFourVectord2FU_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_FU, hpeakd1_Lam_ALam_FU, hpeakd2_Lam_ALam_FU, hpeakd12_Lam_ALam_FU, hpeak_Lam_ALam_FU_bias, hpeak_Lam_ALam_FU_biasd1, hpeak_Lam_ALam_FU_biasd2, hpeak_Lam_ALam_FU_biasd12, Ntk_Vz_weight);

}

//-------------------True matched + True unmatched-----------------------

if(GoodTrackFourVectorTM_Lam_ALam.size()>=1 && GoodTrackFourVectorTU_ALam_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_Lam_ALam,chi21TM_Lam_ALam,chi22TM_Lam_ALam,GoodTrackFourVectorTU_ALam_Lam,chi21TU_ALam_Lam,chi22TU_ALam_Lam,V0chi2d1_diff_TM_TU_Lam_ALam,V0chi2d2_diff_TM_TU_Lam_ALam,V0chi2d12_diff_TM_TU_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTM_Lam_ALam, GoodTrackFourVectord1TM_Lam_ALam, GoodTrackFourVectord2TM_Lam_ALam, GoodTrackFourVectorTU_ALam_Lam, GoodTrackFourVectord1TU_ALam_Lam, GoodTrackFourVectord2TU_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TM_TU, hpeakd1_Lam_ALam_TM_TU, hpeakd2_Lam_ALam_TM_TU, hpeakd12_Lam_ALam_TM_TU, hpeak_Lam_ALam_TM_TU_bias, hpeak_Lam_ALam_TM_TU_biasd1, hpeak_Lam_ALam_TM_TU_biasd2, hpeak_Lam_ALam_TM_TU_biasd12, Ntk_Vz_weight);

ev_z_vtxLAL_TM_TU.push_back(aux_vtxz);
ev_ntrkoff_vecLAL_TM_TU.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLAL_TM_TU.push_back(GoodTrackFourVectorTM_Lam_ALam); 
ev_GoodTrackFourVector_vecALL_TM_TU.push_back(GoodTrackFourVectorTU_ALam_Lam); 


}

if(GoodTrackFourVectorTU_Lam_ALam.size()>=1 && GoodTrackFourVectorTM_ALam_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_Lam_ALam,chi21TU_Lam_ALam,chi22TU_Lam_ALam,GoodTrackFourVectorTM_ALam_Lam,chi21TM_ALam_Lam,chi22TM_ALam_Lam,V0chi2d1_diff_TM_TU_Lam_ALam,V0chi2d2_diff_TM_TU_Lam_ALam,V0chi2d12_diff_TM_TU_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTU_Lam_ALam, GoodTrackFourVectord1TU_Lam_ALam, GoodTrackFourVectord2TU_Lam_ALam, GoodTrackFourVectorTM_ALam_Lam, GoodTrackFourVectord1TM_ALam_Lam, GoodTrackFourVectord2TM_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TM_TU, hpeakd1_Lam_ALam_TM_TU, hpeakd2_Lam_ALam_TM_TU, hpeakd12_Lam_ALam_TM_TU, hpeak_Lam_ALam_TM_TU_bias, hpeak_Lam_ALam_TM_TU_biasd1, hpeak_Lam_ALam_TM_TU_biasd2, hpeak_Lam_ALam_TM_TU_biasd12, Ntk_Vz_weight);

ev_z_vtxLAL_TU_TM.push_back(aux_vtxz);
ev_ntrkoff_vecLAL_TU_TM.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLAL_TU_TM.push_back(GoodTrackFourVectorTU_Lam_ALam); 
ev_GoodTrackFourVector_vecALL_TU_TM.push_back(GoodTrackFourVectorTM_ALam_Lam); 


}

//-------------------Fake matched + Fake unmatched-----------------------

if(GoodTrackFourVectorFM_Lam_ALam.size()>=1 && GoodTrackFourVectorFU_ALam_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFM_Lam_ALam,chi21FM_Lam_ALam,chi22FM_Lam_ALam,GoodTrackFourVectorFU_ALam_Lam,chi21FU_ALam_Lam,chi22FU_ALam_Lam,V0chi2d1_diff_FM_FU_Lam_ALam,V0chi2d2_diff_FM_FU_Lam_ALam,V0chi2d12_diff_FM_FU_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFM_Lam_ALam, GoodTrackFourVectord1FM_Lam_ALam, GoodTrackFourVectord2FM_Lam_ALam, GoodTrackFourVectorFU_ALam_Lam, GoodTrackFourVectord1FU_ALam_Lam, GoodTrackFourVectord2FU_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_FM_FU, hpeakd1_Lam_ALam_FM_FU, hpeakd2_Lam_ALam_FM_FU, hpeakd12_Lam_ALam_FM_FU, hpeak_Lam_ALam_FM_FU_bias, hpeak_Lam_ALam_FM_FU_biasd1, hpeak_Lam_ALam_FM_FU_biasd2, hpeak_Lam_ALam_FM_FU_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFU_Lam_ALam.size()>=1 && GoodTrackFourVectorFM_ALam_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFU_Lam_ALam,chi21FU_Lam_ALam,chi22FU_Lam_ALam,GoodTrackFourVectorFM_ALam_Lam,chi21FM_ALam_Lam,chi22FM_ALam_Lam,V0chi2d1_diff_FM_FU_Lam_ALam,V0chi2d2_diff_FM_FU_Lam_ALam,V0chi2d12_diff_FM_FU_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFU_Lam_ALam, GoodTrackFourVectord1FU_Lam_ALam, GoodTrackFourVectord2FU_Lam_ALam, GoodTrackFourVectorFM_ALam_Lam, GoodTrackFourVectord1FM_ALam_Lam, GoodTrackFourVectord2FM_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_FM_FU, hpeakd1_Lam_ALam_FM_FU, hpeakd2_Lam_ALam_FM_FU, hpeakd12_Lam_ALam_FM_FU, hpeak_Lam_ALam_FM_FU_bias, hpeak_Lam_ALam_FM_FU_biasd1, hpeak_Lam_ALam_FM_FU_biasd2, hpeak_Lam_ALam_FM_FU_biasd12, Ntk_Vz_weight);

}

//-------------------True matched + Fake unmatched-----------------------

if(GoodTrackFourVectorTM_Lam_ALam.size()>=1 && GoodTrackFourVectorFU_ALam_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_Lam_ALam,chi21TM_Lam_ALam,chi22TM_Lam_ALam,GoodTrackFourVectorFU_ALam_Lam,chi21FU_ALam_Lam,chi22FU_ALam_Lam,V0chi2d1_diff_TM_FU_Lam_ALam,V0chi2d2_diff_TM_FU_Lam_ALam,V0chi2d12_diff_TM_FU_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTM_Lam_ALam, GoodTrackFourVectord1TM_Lam_ALam, GoodTrackFourVectord2TM_Lam_ALam, GoodTrackFourVectorFU_ALam_Lam, GoodTrackFourVectord1FU_ALam_Lam, GoodTrackFourVectord2FU_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TM_FU, hpeakd1_Lam_ALam_TM_FU, hpeakd2_Lam_ALam_TM_FU, hpeakd12_Lam_ALam_TM_FU, hpeak_Lam_ALam_TM_FU_bias, hpeak_Lam_ALam_TM_FU_biasd1, hpeak_Lam_ALam_TM_FU_biasd2, hpeak_Lam_ALam_TM_FU_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFU_Lam_ALam.size()>=1 && GoodTrackFourVectorTM_ALam_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFU_Lam_ALam,chi21FU_Lam_ALam,chi22FU_Lam_ALam,GoodTrackFourVectorTM_ALam_Lam,chi21TM_ALam_Lam,chi22TM_ALam_Lam,V0chi2d1_diff_TM_FU_Lam_ALam,V0chi2d2_diff_TM_FU_Lam_ALam,V0chi2d12_diff_TM_FU_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFU_Lam_ALam, GoodTrackFourVectord1FU_Lam_ALam, GoodTrackFourVectord2FU_Lam_ALam, GoodTrackFourVectorTM_ALam_Lam, GoodTrackFourVectord1TM_ALam_Lam, GoodTrackFourVectord2TM_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TM_FU, hpeakd1_Lam_ALam_TM_FU, hpeakd2_Lam_ALam_TM_FU, hpeakd12_Lam_ALam_TM_FU, hpeak_Lam_ALam_TM_FU_bias, hpeak_Lam_ALam_TM_FU_biasd1, hpeak_Lam_ALam_TM_FU_biasd2, hpeak_Lam_ALam_TM_FU_biasd12, Ntk_Vz_weight);

}

//-------------------True matched + Fake matched-----------------------

if(GoodTrackFourVectorTM_Lam_ALam.size()>=1 && GoodTrackFourVectorFM_ALam_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_Lam_ALam,chi21TM_Lam_ALam,chi22TM_Lam_ALam,GoodTrackFourVectorFM_ALam_Lam,chi21FM_ALam_Lam,chi22FM_ALam_Lam,V0chi2d1_diff_TM_FM_Lam_ALam,V0chi2d2_diff_TM_FM_Lam_ALam,V0chi2d12_diff_TM_FM_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTM_Lam_ALam, GoodTrackFourVectord1TM_Lam_ALam, GoodTrackFourVectord2TM_Lam_ALam, GoodTrackFourVectorFM_ALam_Lam, GoodTrackFourVectord1FM_ALam_Lam, GoodTrackFourVectord2FM_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TM_FM, hpeakd1_Lam_ALam_TM_FM, hpeakd2_Lam_ALam_TM_FM, hpeakd12_Lam_ALam_TM_FM, hpeak_Lam_ALam_TM_FM_bias, hpeak_Lam_ALam_TM_FM_biasd1, hpeak_Lam_ALam_TM_FM_biasd2, hpeak_Lam_ALam_TM_FM_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFM_Lam_ALam.size()>=1 && GoodTrackFourVectorTM_ALam_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFM_Lam_ALam,chi21FM_Lam_ALam,chi22FM_Lam_ALam,GoodTrackFourVectorTM_ALam_Lam,chi21TM_ALam_Lam,chi22TM_ALam_Lam,V0chi2d1_diff_TM_FM_Lam_ALam,V0chi2d2_diff_TM_FM_Lam_ALam,V0chi2d12_diff_TM_FM_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFM_Lam_ALam, GoodTrackFourVectord1FM_Lam_ALam, GoodTrackFourVectord2FM_Lam_ALam, GoodTrackFourVectorTM_ALam_Lam, GoodTrackFourVectord1TM_ALam_Lam, GoodTrackFourVectord2TM_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TM_FM, hpeakd1_Lam_ALam_TM_FM, hpeakd2_Lam_ALam_TM_FM, hpeakd12_Lam_ALam_TM_FM, hpeak_Lam_ALam_TM_FM_bias, hpeak_Lam_ALam_TM_FM_biasd1, hpeak_Lam_ALam_TM_FM_biasd2, hpeak_Lam_ALam_TM_FM_biasd12, Ntk_Vz_weight);

}

//-------------------True unmatched + Fake unmatched-----------------------

if(GoodTrackFourVectorTU_Lam_ALam.size()>=1 && GoodTrackFourVectorFU_ALam_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_Lam_ALam,chi21TU_Lam_ALam,chi22TU_Lam_ALam,GoodTrackFourVectorFU_ALam_Lam,chi21FU_ALam_Lam,chi22FU_ALam_Lam,V0chi2d1_diff_TU_FU_Lam_ALam,V0chi2d2_diff_TU_FU_Lam_ALam,V0chi2d12_diff_TU_FU_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTU_Lam_ALam, GoodTrackFourVectord1TU_Lam_ALam, GoodTrackFourVectord2TU_Lam_ALam, GoodTrackFourVectorFU_ALam_Lam, GoodTrackFourVectord1FU_ALam_Lam, GoodTrackFourVectord2FU_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TU_FU, hpeakd1_Lam_ALam_TU_FU, hpeakd2_Lam_ALam_TU_FU, hpeakd12_Lam_ALam_TU_FU, hpeak_Lam_ALam_TU_FU_bias, hpeak_Lam_ALam_TU_FU_biasd1, hpeak_Lam_ALam_TU_FU_biasd2, hpeak_Lam_ALam_TU_FU_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFU_Lam_ALam.size()>=1 && GoodTrackFourVectorTU_ALam_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFU_Lam_ALam,chi21FU_Lam_ALam,chi22FU_Lam_ALam,GoodTrackFourVectorTU_ALam_Lam,chi21TU_ALam_Lam,chi22TU_ALam_Lam,V0chi2d1_diff_TU_FU_Lam_ALam,V0chi2d2_diff_TU_FU_Lam_ALam,V0chi2d12_diff_TU_FU_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFU_Lam_ALam, GoodTrackFourVectord1FU_Lam_ALam, GoodTrackFourVectord2FU_Lam_ALam, GoodTrackFourVectorTU_ALam_Lam, GoodTrackFourVectord1TU_ALam_Lam, GoodTrackFourVectord2TU_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TU_FU, hpeakd1_Lam_ALam_TU_FU, hpeakd2_Lam_ALam_TU_FU, hpeakd12_Lam_ALam_TU_FU, hpeak_Lam_ALam_TU_FU_bias, hpeak_Lam_ALam_TU_FU_biasd1, hpeak_Lam_ALam_TU_FU_biasd2, hpeak_Lam_ALam_TU_FU_biasd12, Ntk_Vz_weight);

}

//-------------------True unmatched + Fake matched-----------------------

if(GoodTrackFourVectorTU_Lam_ALam.size()>=1 && GoodTrackFourVectorFM_ALam_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_Lam_ALam,chi21TU_Lam_ALam,chi22TU_Lam_ALam,GoodTrackFourVectorFM_ALam_Lam,chi21FM_ALam_Lam,chi22FM_ALam_Lam,V0chi2d1_diff_TU_FM_Lam_ALam,V0chi2d2_diff_TU_FM_Lam_ALam,V0chi2d12_diff_TU_FM_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTU_Lam_ALam, GoodTrackFourVectord1TU_Lam_ALam, GoodTrackFourVectord2TU_Lam_ALam, GoodTrackFourVectorFM_ALam_Lam, GoodTrackFourVectord1FM_ALam_Lam, GoodTrackFourVectord2FM_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TU_FM, hpeakd1_Lam_ALam_TU_FM, hpeakd2_Lam_ALam_TU_FM, hpeakd12_Lam_ALam_TU_FM, hpeak_Lam_ALam_TU_FM_bias, hpeak_Lam_ALam_TU_FM_biasd1, hpeak_Lam_ALam_TU_FM_biasd2, hpeak_Lam_ALam_TU_FM_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFM_Lam_ALam.size()>=1 && GoodTrackFourVectorTU_ALam_Lam.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFM_Lam_ALam,chi21FM_Lam_ALam,chi22FM_Lam_ALam,GoodTrackFourVectorTU_ALam_Lam,chi21TU_ALam_Lam,chi22TU_ALam_Lam,V0chi2d1_diff_TU_FM_Lam_ALam,V0chi2d2_diff_TU_FM_Lam_ALam,V0chi2d12_diff_TU_FM_Lam_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFM_Lam_ALam, GoodTrackFourVectord1FM_Lam_ALam, GoodTrackFourVectord2FM_Lam_ALam, GoodTrackFourVectorTU_ALam_Lam, GoodTrackFourVectord1TU_ALam_Lam, GoodTrackFourVectord2TU_ALam_Lam, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TU_FM, hpeakd1_Lam_ALam_TU_FM, hpeakd2_Lam_ALam_TU_FM, hpeakd12_Lam_ALam_TU_FM, hpeak_Lam_ALam_TU_FM_bias, hpeak_Lam_ALam_TU_FM_biasd1, hpeak_Lam_ALam_TU_FM_biasd2, hpeak_Lam_ALam_TU_FM_biasd12, Ntk_Vz_weight);

}

}

/*========================= K0sLambda =========================*/

if(GoodTrackFourVectorT_K0s_Lam.size()>=1 && GoodTrackFourVectorT_Lam_K0s.size()>=1){

nev_KLT->Fill(1);

Mass_dep_OS(GoodTrackFourVectorT_K0s_Lam,GoodTrackFourVectorT_Lam_K0s,(Double_t) aux_N_tk_offline,hMass_K0sL_T,hMass_LamK_T,hMass_K0s_Lam_T, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorT_K0s_Lam,chi21T_K0s_Lam,chi22T_K0s_Lam,GoodTrackFourVectorT_Lam_K0s,chi21T_Lam_K0s,chi22T_Lam_K0s,V0chi2d1_diff_T_K0s_Lam,V0chi2d2_diff_T_K0s_Lam,V0chi2d12_diff_T_K0s_Lam, Ntk_Vz_weight);

hbt_reco_cross(GoodTrackFourVectorT_K0s_Lam, GoodTrackFourVectord1T_K0s_Lam, GoodTrackFourVectord2T_K0s_Lam, GoodTrackFourVectorT_Lam_K0s, GoodTrackFourVectord1T_Lam_K0s, GoodTrackFourVectord2T_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_T, hpeakcos_K0s_Lam_T, hpeakd1_K0s_Lam_T, hpeakd2_K0s_Lam_T, hpeakd12_K0s_Lam_T, hpeak_rot_K0s_Lam_T, hpeak_inv_K0s_Lam_T, hside_K0s_Lam_T, hside_rot_K0s_Lam_T, hside_inv_K0s_Lam_T, hpeakside_K0s_Lam_T, hpeakside_rot_K0s_Lam_T, hpeakside_inv_K0s_Lam_T, hsideL_K0s_Lam_T, hsideR_K0s_Lam_T, hpeaksideL_K0s_Lam_T, hpeaksideR_K0s_Lam_T, p, effhisto_K0s, effhisto_Lam, nev_KL_ssT, nev_KL_bbT, nev_KL_sbT, Ntk_Vz_weight);

hdibaryon_OS2(GoodTrackFourVectorT_K0s_Lam, GoodTrackFourVectorT_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, h_Hcasc1820_K0s_Lam, h_Hcasc1820_rot_K0s_Lam, h_Hcasc1820_inv_K0s_Lam, h_Hcasc1820_K0s_Lam_side, h_Hcasc1820_rot_K0s_Lam_side, h_Hcasc1820_inv_K0s_Lam_side, h_Hcasc1820_K0s_Lam_peakside, h_Hcasc1820_rot_K0s_Lam_peakside, h_Hcasc1820_inv_K0s_Lam_peakside, Ntk_Vz_weight);

weight_K0sL.push_back(Ntk_Vz_weight);

ev_z_vtxK0sLam.push_back(aux_vtxz);
ev_ntrkoff_vecK0sLam.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0sLam.push_back(GoodTrackFourVectorT_K0s_Lam); 
ev_GoodTrackFourVector_vecLamK0s.push_back(GoodTrackFourVectorT_Lam_K0s); 

Double_t aux_etaMix_w = Eta_weight;
std::pair<Double_t, std::vector<TLorentzVector>> aux_pair_GoodTrackFourVector_etaMixWeight = make_pair(aux_etaMix_w, GoodTrackFourVectorT_K0s_Lam);
ev_GoodTrackFourVector_etaMixWeight_vecK0sLam.push_back(aux_pair_GoodTrackFourVector_etaMixWeight);
std::pair<Double_t, std::vector<TLorentzVector>> aux_pair_GoodTrackFourVector_etaMixWeight2 = make_pair(aux_etaMix_w, GoodTrackFourVectorT_Lam_K0s);
ev_GoodTrackFourVector_etaMixWeight_vecLamK0s.push_back(aux_pair_GoodTrackFourVector_etaMixWeight2);
std::pair<Double_t,int> aux_pair_ntrkoff_etaMixWeight = make_pair(aux_etaMix_w,aux_N_tk_offline);
ev_ntrkoff_etaMixWeight_vecK0sLam.push_back(aux_pair_ntrkoff_etaMixWeight); 

}

if(GoodTrackFourVectorF_K0s_Lam.size()>=1 && GoodTrackFourVectorF_Lam_K0s.size()>=1){

nev_KLF->Fill(1);

Mass_dep_OS(GoodTrackFourVectorF_K0s_Lam,GoodTrackFourVectorF_Lam_K0s,(Double_t) aux_N_tk_offline,hMass_K0sL_F,hMass_LamK_F,hMass_K0s_Lam_F, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorF_K0s_Lam, chi21F_K0s_Lam, chi22F_K0s_Lam, GoodTrackFourVectorF_Lam_K0s, chi21F_Lam_K0s, chi22F_Lam_K0s, V0chi2d1_diff_F_K0s_Lam, V0chi2d2_diff_F_K0s_Lam, V0chi2d12_diff_F_K0s_Lam, Ntk_Vz_weight);

hbt_reco_truefake(GoodTrackFourVectorF_K0s_Lam, GoodTrackFourVectord1F_K0s_Lam, GoodTrackFourVectord2F_K0s_Lam, GoodTrackFourVectorF_Lam_K0s, GoodTrackFourVectord1F_Lam_K0s, GoodTrackFourVectord2F_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_F, hpeakd1_K0s_Lam_F, hpeakd2_K0s_Lam_F, hpeakd12_K0s_Lam_F, hside_K0s_Lam_F, hpeakside_K0s_Lam_F, hpeak_K0s_Lam_F_bias, hpeak_K0s_Lam_F_biasd1, hpeak_K0s_Lam_F_biasd2, hpeak_K0s_Lam_F_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorF_K0s_Lam.size()>=1 && GoodTrackFourVectorT_Lam_K0s.size()>=1){

nev_KLTF->Fill(1);

get_chi2plot_TF(GoodTrackFourVectorF_K0s_Lam, chi21F_K0s_Lam, chi22F_K0s_Lam, GoodTrackFourVectorT_Lam_K0s, chi21T_Lam_K0s, chi22T_Lam_K0s, V0chi2d1_diff_TF_K0s_Lam, V0chi2d2_diff_TF_K0s_Lam, V0chi2d12_diff_TF_K0s_Lam, Ntk_Vz_weight);

hbt_reco_truefake(GoodTrackFourVectorF_K0s_Lam, GoodTrackFourVectord1F_K0s_Lam, GoodTrackFourVectord2F_K0s_Lam, GoodTrackFourVectorT_Lam_K0s, GoodTrackFourVectord1T_Lam_K0s, GoodTrackFourVectord2T_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TF, hpeakd1_K0s_Lam_TF, hpeakd2_K0s_Lam_TF, hpeakd12_K0s_Lam_TF, hside_K0s_Lam_TF, hpeakside_K0s_Lam_TF, hpeak_K0s_Lam_TF_bias, hpeak_K0s_Lam_TF_biasd1, hpeak_K0s_Lam_TF_biasd2, hpeak_K0s_Lam_TF_biasd12, Ntk_Vz_weight);
    
}

if(GoodTrackFourVectorT_K0s_Lam.size()>=1 && GoodTrackFourVectorF_Lam_K0s.size()>=1){

nev_KLTF->Fill(1);

get_chi2plot_TF(GoodTrackFourVectorT_K0s_Lam, chi21T_K0s_Lam, chi22T_K0s_Lam, GoodTrackFourVectorF_Lam_K0s, chi21F_Lam_K0s, chi22F_Lam_K0s, V0chi2d1_diff_TF_K0s_Lam, V0chi2d2_diff_TF_K0s_Lam, V0chi2d12_diff_TF_K0s_Lam, Ntk_Vz_weight);

hbt_reco_truefake(GoodTrackFourVectorT_K0s_Lam, GoodTrackFourVectord1T_K0s_Lam, GoodTrackFourVectord2T_K0s_Lam, GoodTrackFourVectorF_Lam_K0s, GoodTrackFourVectord1F_Lam_K0s, GoodTrackFourVectord2F_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TF, hpeakd1_K0s_Lam_TF, hpeakd2_K0s_Lam_TF, hpeakd12_K0s_Lam_TF, hside_K0s_Lam_TF, hpeakside_K0s_Lam_TF, hpeak_K0s_Lam_TF_bias, hpeak_K0s_Lam_TF_biasd1, hpeak_K0s_Lam_TF_biasd2, hpeak_K0s_Lam_TF_biasd12, Ntk_Vz_weight);
    
}

if(isMC){

//-------------------True Matched-----------------------

if(GoodTrackFourVectorTM_K0s_Lam.size()>=1 && GoodTrackFourVectorTM_Lam_K0s.size()>=1){

Mass_dep_OS(GoodTrackFourVectorTM_K0s_Lam,GoodTrackFourVectorTM_Lam_K0s,(Double_t) aux_N_tk_offline,hMass_K0sL_TM,hMass_LamK_TM,hMass_K0s_Lam_TM, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorTM_K0s_Lam,chi21TM_K0s_Lam,chi22TM_K0s_Lam,GoodTrackFourVectorTM_Lam_K0s,chi21TM_Lam_K0s,chi22TM_Lam_K0s,V0chi2d1_diff_TM_K0s_Lam,V0chi2d2_diff_TM_K0s_Lam,V0chi2d12_diff_TM_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTM_K0s_Lam, GoodTrackFourVectord1TM_K0s_Lam, GoodTrackFourVectord2TM_K0s_Lam, GoodTrackFourVectorTM_Lam_K0s, GoodTrackFourVectord1TM_Lam_K0s, GoodTrackFourVectord2TM_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TM, hpeakd1_K0s_Lam_TM, hpeakd2_K0s_Lam_TM, hpeakd12_K0s_Lam_TM, hpeak_K0s_Lam_TM_bias, hpeak_K0s_Lam_TM_biasd1, hpeak_K0s_Lam_TM_biasd2, hpeak_K0s_Lam_TM_biasd12, Ntk_Vz_weight);

ev_z_vtxKL_TM.push_back(aux_vtxz);
ev_ntrkoff_vecKL_TM.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecKL_TM.push_back(GoodTrackFourVectorTM_K0s_Lam); 
ev_GoodTrackFourVector_vecLK_TM.push_back(GoodTrackFourVectorTM_Lam_K0s); 


}

//-------------------True Unmatched-----------------------

if(GoodTrackFourVectorTU_K0s_Lam.size()>=1 && GoodTrackFourVectorTU_Lam_K0s.size()>=1){

Mass_dep_OS(GoodTrackFourVectorTU_K0s_Lam,GoodTrackFourVectorTU_Lam_K0s,(Double_t) aux_N_tk_offline,hMass_K0sL_TU,hMass_LamK_TU,hMass_K0s_Lam_TU, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorTU_K0s_Lam,chi21TU_K0s_Lam,chi22TU_K0s_Lam,GoodTrackFourVectorTU_Lam_K0s,chi21TU_Lam_K0s,chi22TU_Lam_K0s,V0chi2d1_diff_TU_K0s_Lam,V0chi2d2_diff_TU_K0s_Lam,V0chi2d12_diff_TU_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTU_K0s_Lam, GoodTrackFourVectord1TU_K0s_Lam, GoodTrackFourVectord2TU_K0s_Lam, GoodTrackFourVectorTU_Lam_K0s, GoodTrackFourVectord1TU_Lam_K0s, GoodTrackFourVectord2TU_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TU, hpeakd1_K0s_Lam_TU, hpeakd2_K0s_Lam_TU, hpeakd12_K0s_Lam_TU, hpeak_K0s_Lam_TU_bias, hpeak_K0s_Lam_TU_biasd1, hpeak_K0s_Lam_TU_biasd2, hpeak_K0s_Lam_TU_biasd12, Ntk_Vz_weight);

ev_z_vtxKL_TU.push_back(aux_vtxz);
ev_ntrkoff_vecKL_TU.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecKL_TU.push_back(GoodTrackFourVectorTU_K0s_Lam); 
ev_GoodTrackFourVector_vecLK_TU.push_back(GoodTrackFourVectorTU_Lam_K0s); 


}

//-------------------Fake Matched-----------------------

if(GoodTrackFourVectorFM_K0s_Lam.size()>=1 && GoodTrackFourVectorFM_Lam_K0s.size()>=1){

Mass_dep_OS(GoodTrackFourVectorFM_K0s_Lam,GoodTrackFourVectorFM_Lam_K0s,(Double_t) aux_N_tk_offline,hMass_K0sL_FM,hMass_LamK_FM,hMass_K0s_Lam_FM, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorFM_K0s_Lam,chi21FM_K0s_Lam,chi22FM_K0s_Lam,GoodTrackFourVectorFM_Lam_K0s,chi21FM_Lam_K0s,chi22FM_Lam_K0s,V0chi2d1_diff_FM_K0s_Lam,V0chi2d2_diff_FM_K0s_Lam,V0chi2d12_diff_FM_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFM_K0s_Lam, GoodTrackFourVectord1FM_K0s_Lam, GoodTrackFourVectord2FM_K0s_Lam, GoodTrackFourVectorFM_Lam_K0s, GoodTrackFourVectord1FM_Lam_K0s, GoodTrackFourVectord2FM_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_FM, hpeakd1_K0s_Lam_FM, hpeakd2_K0s_Lam_FM, hpeakd12_K0s_Lam_FM, hpeak_K0s_Lam_FM_bias, hpeak_K0s_Lam_FM_biasd1, hpeak_K0s_Lam_FM_biasd2, hpeak_K0s_Lam_FM_biasd12, Ntk_Vz_weight);

}

//-------------------True Unmatched-----------------------

if(GoodTrackFourVectorFU_K0s_Lam.size()>=1 && GoodTrackFourVectorFU_Lam_K0s.size()>=1){

Mass_dep_OS(GoodTrackFourVectorFU_K0s_Lam,GoodTrackFourVectorFU_Lam_K0s,(Double_t) aux_N_tk_offline,hMass_K0sL_FU,hMass_LamK_FU,hMass_K0s_Lam_FU, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorFU_K0s_Lam,chi21FU_K0s_Lam,chi22FU_K0s_Lam,GoodTrackFourVectorFU_Lam_K0s,chi21FU_Lam_K0s,chi22FU_Lam_K0s,V0chi2d1_diff_FU_K0s_Lam,V0chi2d2_diff_FU_K0s_Lam,V0chi2d12_diff_FU_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFU_K0s_Lam, GoodTrackFourVectord1FU_K0s_Lam, GoodTrackFourVectord2FU_K0s_Lam, GoodTrackFourVectorFU_Lam_K0s, GoodTrackFourVectord1FU_Lam_K0s, GoodTrackFourVectord2FU_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_FU, hpeakd1_K0s_Lam_FU, hpeakd2_K0s_Lam_FU, hpeakd12_K0s_Lam_FU, hpeak_K0s_Lam_FU_bias, hpeak_K0s_Lam_FU_biasd1, hpeak_K0s_Lam_FU_biasd2, hpeak_K0s_Lam_FU_biasd12, Ntk_Vz_weight);

}

//-------------------True matched + True unmatched-----------------------

if(GoodTrackFourVectorTM_K0s_Lam.size()>=1 && GoodTrackFourVectorTU_Lam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_K0s_Lam,chi21TM_K0s_Lam,chi22TM_K0s_Lam,GoodTrackFourVectorTU_Lam_K0s,chi21TU_Lam_K0s,chi22TU_Lam_K0s,V0chi2d1_diff_TM_TU_K0s_Lam,V0chi2d2_diff_TM_TU_K0s_Lam,V0chi2d12_diff_TM_TU_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTM_K0s_Lam, GoodTrackFourVectord1TM_K0s_Lam, GoodTrackFourVectord2TM_K0s_Lam, GoodTrackFourVectorTU_Lam_K0s, GoodTrackFourVectord1TU_Lam_K0s, GoodTrackFourVectord2TU_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TM_TU, hpeakd1_K0s_Lam_TM_TU, hpeakd2_K0s_Lam_TM_TU, hpeakd12_K0s_Lam_TM_TU, hpeak_K0s_Lam_TM_TU_bias, hpeak_K0s_Lam_TM_TU_biasd1, hpeak_K0s_Lam_TM_TU_biasd2, hpeak_K0s_Lam_TM_TU_biasd12, Ntk_Vz_weight);

ev_z_vtxKL_TM_TU.push_back(aux_vtxz);
ev_ntrkoff_vecKL_TM_TU.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecKL_TM_TU.push_back(GoodTrackFourVectorTM_K0s_Lam); 
ev_GoodTrackFourVector_vecLK_TM_TU.push_back(GoodTrackFourVectorTU_Lam_K0s); 


}

if(GoodTrackFourVectorTU_K0s_Lam.size()>=1 && GoodTrackFourVectorTM_Lam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_K0s_Lam,chi21TU_K0s_Lam,chi22TU_K0s_Lam,GoodTrackFourVectorTM_Lam_K0s,chi21TM_Lam_K0s,chi22TM_Lam_K0s,V0chi2d1_diff_TM_TU_K0s_Lam,V0chi2d2_diff_TM_TU_K0s_Lam,V0chi2d12_diff_TM_TU_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTU_K0s_Lam, GoodTrackFourVectord1TU_K0s_Lam, GoodTrackFourVectord2TU_K0s_Lam, GoodTrackFourVectorTM_Lam_K0s, GoodTrackFourVectord1TM_Lam_K0s, GoodTrackFourVectord2TM_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TM_TU, hpeakd1_K0s_Lam_TM_TU, hpeakd2_K0s_Lam_TM_TU, hpeakd12_K0s_Lam_TM_TU, hpeak_K0s_Lam_TM_TU_bias, hpeak_K0s_Lam_TM_TU_biasd1, hpeak_K0s_Lam_TM_TU_biasd2, hpeak_K0s_Lam_TM_TU_biasd12, Ntk_Vz_weight);

ev_z_vtxKL_TU_TM.push_back(aux_vtxz);
ev_ntrkoff_vecKL_TU_TM.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecKL_TU_TM.push_back(GoodTrackFourVectorTU_K0s_Lam); 
ev_GoodTrackFourVector_vecLK_TU_TM.push_back(GoodTrackFourVectorTM_Lam_K0s); 


}

//-------------------Fake matched + Fake unmatched-----------------------

if(GoodTrackFourVectorFM_K0s_Lam.size()>=1 && GoodTrackFourVectorFU_Lam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFM_K0s_Lam,chi21FM_K0s_Lam,chi22FM_K0s_Lam,GoodTrackFourVectorFU_Lam_K0s,chi21FU_Lam_K0s,chi22FU_Lam_K0s,V0chi2d1_diff_FM_FU_K0s_Lam,V0chi2d2_diff_FM_FU_K0s_Lam,V0chi2d12_diff_FM_FU_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFM_K0s_Lam, GoodTrackFourVectord1FM_K0s_Lam, GoodTrackFourVectord2FM_K0s_Lam, GoodTrackFourVectorFU_Lam_K0s, GoodTrackFourVectord1FU_Lam_K0s, GoodTrackFourVectord2FU_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_FM_FU, hpeakd1_K0s_Lam_FM_FU, hpeakd2_K0s_Lam_FM_FU, hpeakd12_K0s_Lam_FM_FU, hpeak_K0s_Lam_FM_FU_bias, hpeak_K0s_Lam_FM_FU_biasd1, hpeak_K0s_Lam_FM_FU_biasd2, hpeak_K0s_Lam_FM_FU_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFU_K0s_Lam.size()>=1 && GoodTrackFourVectorFM_Lam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFU_K0s_Lam,chi21FU_K0s_Lam,chi22FU_K0s_Lam,GoodTrackFourVectorFM_Lam_K0s,chi21FM_Lam_K0s,chi22FM_Lam_K0s,V0chi2d1_diff_FM_FU_K0s_Lam,V0chi2d2_diff_FM_FU_K0s_Lam,V0chi2d12_diff_FM_FU_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFU_K0s_Lam, GoodTrackFourVectord1FU_K0s_Lam, GoodTrackFourVectord2FU_K0s_Lam, GoodTrackFourVectorFM_Lam_K0s, GoodTrackFourVectord1FM_Lam_K0s, GoodTrackFourVectord2FM_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_FM_FU, hpeakd1_K0s_Lam_FM_FU, hpeakd2_K0s_Lam_FM_FU, hpeakd12_K0s_Lam_FM_FU, hpeak_K0s_Lam_FM_FU_bias, hpeak_K0s_Lam_FM_FU_biasd1, hpeak_K0s_Lam_FM_FU_biasd2, hpeak_K0s_Lam_FM_FU_biasd12, Ntk_Vz_weight);

}

//-------------------True matched + Fake unmatched-----------------------

if(GoodTrackFourVectorTM_K0s_Lam.size()>=1 && GoodTrackFourVectorFU_Lam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_K0s_Lam,chi21TM_K0s_Lam,chi22TM_K0s_Lam,GoodTrackFourVectorFU_Lam_K0s,chi21FU_Lam_K0s,chi22FU_Lam_K0s,V0chi2d1_diff_TM_FU_K0s_Lam,V0chi2d2_diff_TM_FU_K0s_Lam,V0chi2d12_diff_TM_FU_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTM_K0s_Lam, GoodTrackFourVectord1TM_K0s_Lam, GoodTrackFourVectord2TM_K0s_Lam, GoodTrackFourVectorFU_Lam_K0s, GoodTrackFourVectord1FU_Lam_K0s, GoodTrackFourVectord2FU_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TM_FU, hpeakd1_K0s_Lam_TM_FU, hpeakd2_K0s_Lam_TM_FU, hpeakd12_K0s_Lam_TM_FU, hpeak_K0s_Lam_TM_FU_bias, hpeak_K0s_Lam_TM_FU_biasd1, hpeak_K0s_Lam_TM_FU_biasd2, hpeak_K0s_Lam_TM_FU_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFU_K0s_Lam.size()>=1 && GoodTrackFourVectorTM_Lam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFU_K0s_Lam,chi21FU_K0s_Lam,chi22FU_K0s_Lam,GoodTrackFourVectorTM_Lam_K0s,chi21TM_Lam_K0s,chi22TM_Lam_K0s,V0chi2d1_diff_TM_FU_K0s_Lam,V0chi2d2_diff_TM_FU_K0s_Lam,V0chi2d12_diff_TM_FU_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFU_K0s_Lam, GoodTrackFourVectord1FU_K0s_Lam, GoodTrackFourVectord2FU_K0s_Lam, GoodTrackFourVectorTM_Lam_K0s, GoodTrackFourVectord1TM_Lam_K0s, GoodTrackFourVectord2TM_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TM_FU, hpeakd1_K0s_Lam_TM_FU, hpeakd2_K0s_Lam_TM_FU, hpeakd12_K0s_Lam_TM_FU, hpeak_K0s_Lam_TM_FU_bias, hpeak_K0s_Lam_TM_FU_biasd1, hpeak_K0s_Lam_TM_FU_biasd2, hpeak_K0s_Lam_TM_FU_biasd12, Ntk_Vz_weight);

}

//-------------------True matched + Fake matched-----------------------

if(GoodTrackFourVectorTM_K0s_Lam.size()>=1 && GoodTrackFourVectorFM_Lam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_K0s_Lam,chi21TM_K0s_Lam,chi22TM_K0s_Lam,GoodTrackFourVectorFM_Lam_K0s,chi21FM_Lam_K0s,chi22FM_Lam_K0s,V0chi2d1_diff_TM_FM_K0s_Lam,V0chi2d2_diff_TM_FM_K0s_Lam,V0chi2d12_diff_TM_FM_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTM_K0s_Lam, GoodTrackFourVectord1TM_K0s_Lam, GoodTrackFourVectord2TM_K0s_Lam, GoodTrackFourVectorFM_Lam_K0s, GoodTrackFourVectord1FM_Lam_K0s, GoodTrackFourVectord2FM_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TM_FM, hpeakd1_K0s_Lam_TM_FM, hpeakd2_K0s_Lam_TM_FM, hpeakd12_K0s_Lam_TM_FM, hpeak_K0s_Lam_TM_FM_bias, hpeak_K0s_Lam_TM_FM_biasd1, hpeak_K0s_Lam_TM_FM_biasd2, hpeak_K0s_Lam_TM_FM_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFM_K0s_Lam.size()>=1 && GoodTrackFourVectorTM_Lam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFM_K0s_Lam,chi21FM_K0s_Lam,chi22FM_K0s_Lam,GoodTrackFourVectorTM_Lam_K0s,chi21TM_Lam_K0s,chi22TM_Lam_K0s,V0chi2d1_diff_TM_FM_K0s_Lam,V0chi2d2_diff_TM_FM_K0s_Lam,V0chi2d12_diff_TM_FM_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFM_K0s_Lam, GoodTrackFourVectord1FM_K0s_Lam, GoodTrackFourVectord2FM_K0s_Lam, GoodTrackFourVectorTM_Lam_K0s, GoodTrackFourVectord1TM_Lam_K0s, GoodTrackFourVectord2TM_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TM_FM, hpeakd1_K0s_Lam_TM_FM, hpeakd2_K0s_Lam_TM_FM, hpeakd12_K0s_Lam_TM_FM, hpeak_K0s_Lam_TM_FM_bias, hpeak_K0s_Lam_TM_FM_biasd1, hpeak_K0s_Lam_TM_FM_biasd2, hpeak_K0s_Lam_TM_FM_biasd12, Ntk_Vz_weight);

}

//-------------------True unmatched + Fake unmatched-----------------------

if(GoodTrackFourVectorTU_K0s_Lam.size()>=1 && GoodTrackFourVectorFU_Lam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_K0s_Lam,chi21TU_K0s_Lam,chi22TU_K0s_Lam,GoodTrackFourVectorFU_Lam_K0s,chi21FU_Lam_K0s,chi22FU_Lam_K0s,V0chi2d1_diff_TU_FU_K0s_Lam,V0chi2d2_diff_TU_FU_K0s_Lam,V0chi2d12_diff_TU_FU_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTU_K0s_Lam, GoodTrackFourVectord1TU_K0s_Lam, GoodTrackFourVectord2TU_K0s_Lam, GoodTrackFourVectorFU_Lam_K0s, GoodTrackFourVectord1FU_Lam_K0s, GoodTrackFourVectord2FU_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TU_FU, hpeakd1_K0s_Lam_TU_FU, hpeakd2_K0s_Lam_TU_FU, hpeakd12_K0s_Lam_TU_FU, hpeak_K0s_Lam_TU_FU_bias, hpeak_K0s_Lam_TU_FU_biasd1, hpeak_K0s_Lam_TU_FU_biasd2, hpeak_K0s_Lam_TU_FU_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFU_K0s_Lam.size()>=1 && GoodTrackFourVectorTU_Lam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFU_K0s_Lam,chi21FU_K0s_Lam,chi22FU_K0s_Lam,GoodTrackFourVectorTU_Lam_K0s,chi21TU_Lam_K0s,chi22TU_Lam_K0s,V0chi2d1_diff_TU_FU_K0s_Lam,V0chi2d2_diff_TU_FU_K0s_Lam,V0chi2d12_diff_TU_FU_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFU_K0s_Lam, GoodTrackFourVectord1FU_K0s_Lam, GoodTrackFourVectord2FU_K0s_Lam, GoodTrackFourVectorTU_Lam_K0s, GoodTrackFourVectord1TU_Lam_K0s, GoodTrackFourVectord2TU_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TU_FU, hpeakd1_K0s_Lam_TU_FU, hpeakd2_K0s_Lam_TU_FU, hpeakd12_K0s_Lam_TU_FU, hpeak_K0s_Lam_TU_FU_bias, hpeak_K0s_Lam_TU_FU_biasd1, hpeak_K0s_Lam_TU_FU_biasd2, hpeak_K0s_Lam_TU_FU_biasd12, Ntk_Vz_weight);

}

//-------------------True unmatched + Fake matched-----------------------

if(GoodTrackFourVectorTU_K0s_Lam.size()>=1 && GoodTrackFourVectorFM_Lam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_K0s_Lam,chi21TU_K0s_Lam,chi22TU_K0s_Lam,GoodTrackFourVectorFM_Lam_K0s,chi21FM_Lam_K0s,chi22FM_Lam_K0s,V0chi2d1_diff_TU_FM_K0s_Lam,V0chi2d2_diff_TU_FM_K0s_Lam,V0chi2d12_diff_TU_FM_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTU_K0s_Lam, GoodTrackFourVectord1TU_K0s_Lam, GoodTrackFourVectord2TU_K0s_Lam, GoodTrackFourVectorFM_Lam_K0s, GoodTrackFourVectord1FM_Lam_K0s, GoodTrackFourVectord2FM_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TU_FM, hpeakd1_K0s_Lam_TU_FM, hpeakd2_K0s_Lam_TU_FM, hpeakd12_K0s_Lam_TU_FM, hpeak_K0s_Lam_TU_FM_bias, hpeak_K0s_Lam_TU_FM_biasd1, hpeak_K0s_Lam_TU_FM_biasd2, hpeak_K0s_Lam_TU_FM_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFM_K0s_Lam.size()>=1 && GoodTrackFourVectorTU_Lam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFM_K0s_Lam,chi21FM_K0s_Lam,chi22FM_K0s_Lam,GoodTrackFourVectorTU_Lam_K0s,chi21TU_Lam_K0s,chi22TU_Lam_K0s,V0chi2d1_diff_TU_FM_K0s_Lam,V0chi2d2_diff_TU_FM_K0s_Lam,V0chi2d12_diff_TU_FM_K0s_Lam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFM_K0s_Lam, GoodTrackFourVectord1FM_K0s_Lam, GoodTrackFourVectord2FM_K0s_Lam, GoodTrackFourVectorTU_Lam_K0s, GoodTrackFourVectord1TU_Lam_K0s, GoodTrackFourVectord2TU_Lam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TU_FM, hpeakd1_K0s_Lam_TU_FM, hpeakd2_K0s_Lam_TU_FM, hpeakd12_K0s_Lam_TU_FM, hpeak_K0s_Lam_TU_FM_bias, hpeak_K0s_Lam_TU_FM_biasd1, hpeak_K0s_Lam_TU_FM_biasd2, hpeak_K0s_Lam_TU_FM_biasd12, Ntk_Vz_weight);

}

}

/*========================= K0sAntiLambda =========================*/

if(GoodTrackFourVectorT_K0s_ALam.size()>=1 && GoodTrackFourVectorT_ALam_K0s.size()>=1){

nev_KALT->Fill(1);

Mass_dep_OS(GoodTrackFourVectorT_K0s_ALam,GoodTrackFourVectorT_ALam_K0s,(Double_t) aux_N_tk_offline,hMass_K0sAL_T,hMass_ALamK_T,hMass_K0s_ALam_T, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorT_K0s_ALam,chi21T_K0s_ALam,chi22T_K0s_ALam,GoodTrackFourVectorT_ALam_K0s,chi21T_ALam_K0s,chi22T_ALam_K0s,V0chi2d1_diff_T_K0s_ALam,V0chi2d2_diff_T_K0s_ALam,V0chi2d12_diff_T_K0s_ALam, Ntk_Vz_weight);

hbt_reco_cross(GoodTrackFourVectorT_K0s_ALam, GoodTrackFourVectord1T_K0s_ALam, GoodTrackFourVectord2T_K0s_ALam, GoodTrackFourVectorT_ALam_K0s, GoodTrackFourVectord1T_ALam_K0s, GoodTrackFourVectord2T_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_T, hpeakcos_K0s_ALam_T, hpeakd1_K0s_ALam_T, hpeakd2_K0s_ALam_T, hpeakd12_K0s_ALam_T, hpeak_rot_K0s_ALam_T, hpeak_inv_K0s_ALam_T, hside_K0s_ALam_T, hside_rot_K0s_ALam_T, hside_inv_K0s_ALam_T, hpeakside_K0s_ALam_T, hpeakside_rot_K0s_ALam_T, hpeakside_inv_K0s_ALam_T, hsideL_K0s_ALam_T, hsideR_K0s_ALam_T, hpeaksideL_K0s_ALam_T, hpeaksideR_K0s_ALam_T, p, effhisto_K0s, effhisto_ALam, nev_KAL_ssT, nev_KAL_bbT, nev_KAL_sbT, Ntk_Vz_weight);

hdibaryon_OS2(GoodTrackFourVectorT_K0s_ALam, GoodTrackFourVectorT_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, h_Hcasc1820_K0s_ALam, h_Hcasc1820_rot_K0s_ALam, h_Hcasc1820_inv_K0s_ALam, h_Hcasc1820_K0s_ALam_side, h_Hcasc1820_rot_K0s_ALam_side, h_Hcasc1820_inv_K0s_ALam_side, h_Hcasc1820_K0s_ALam_peakside, h_Hcasc1820_rot_K0s_ALam_peakside, h_Hcasc1820_inv_K0s_ALam_peakside, Ntk_Vz_weight);

weight_K0sAL.push_back(Ntk_Vz_weight);

ev_z_vtxK0sALam.push_back(aux_vtxz);
ev_ntrkoff_vecK0sALam.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0sALam.push_back(GoodTrackFourVectorT_K0s_ALam); 
ev_GoodTrackFourVector_vecALamK0s.push_back(GoodTrackFourVectorT_ALam_K0s); 

Double_t aux_etaMix_w = Eta_weight;
std::pair<Double_t, std::vector<TLorentzVector>> aux_pair_GoodTrackFourVector_etaMixWeight = make_pair(aux_etaMix_w, GoodTrackFourVectorT_K0s_ALam);
ev_GoodTrackFourVector_etaMixWeight_vecK0sALam.push_back(aux_pair_GoodTrackFourVector_etaMixWeight);
std::pair<Double_t, std::vector<TLorentzVector>> aux_pair_GoodTrackFourVector_etaMixWeight2 = make_pair(aux_etaMix_w, GoodTrackFourVectorT_ALam_K0s);
ev_GoodTrackFourVector_etaMixWeight_vecALamK0s.push_back(aux_pair_GoodTrackFourVector_etaMixWeight2);
std::pair<Double_t,int> aux_pair_ntrkoff_etaMixWeight = make_pair(aux_etaMix_w,aux_N_tk_offline);
ev_ntrkoff_etaMixWeight_vecK0sALam.push_back(aux_pair_ntrkoff_etaMixWeight); 

}

if(GoodTrackFourVectorF_K0s_ALam.size()>=1 && GoodTrackFourVectorF_ALam_K0s.size()>=1){

nev_KALF->Fill(1);

Mass_dep_OS(GoodTrackFourVectorF_K0s_ALam,GoodTrackFourVectorF_ALam_K0s,(Double_t) aux_N_tk_offline,hMass_K0sAL_F,hMass_ALamK_F,hMass_K0s_ALam_F, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorF_K0s_ALam, chi21F_K0s_ALam, chi22F_K0s_ALam, GoodTrackFourVectorF_ALam_K0s, chi21F_ALam_K0s, chi22F_ALam_K0s, V0chi2d1_diff_F_K0s_ALam, V0chi2d2_diff_F_K0s_ALam, V0chi2d12_diff_F_K0s_ALam, Ntk_Vz_weight);

hbt_reco_truefake(GoodTrackFourVectorF_K0s_ALam, GoodTrackFourVectord1F_K0s_ALam, GoodTrackFourVectord2F_K0s_ALam, GoodTrackFourVectorF_ALam_K0s, GoodTrackFourVectord1F_ALam_K0s, GoodTrackFourVectord2F_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_F, hpeakd1_K0s_ALam_F, hpeakd2_K0s_ALam_F, hpeakd12_K0s_ALam_F, hside_K0s_ALam_F, hpeakside_K0s_ALam_F, hpeak_K0s_ALam_F_bias, hpeak_K0s_ALam_F_biasd1, hpeak_K0s_ALam_F_biasd2, hpeak_K0s_ALam_F_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorF_K0s_ALam.size()>=1 && GoodTrackFourVectorT_ALam_K0s.size()>=1){

nev_KALTF->Fill(1);

get_chi2plot_TF(GoodTrackFourVectorF_K0s_ALam, chi21F_K0s_ALam, chi22F_K0s_ALam, GoodTrackFourVectorT_ALam_K0s, chi21T_ALam_K0s, chi22T_ALam_K0s, V0chi2d1_diff_TF_K0s_ALam, V0chi2d2_diff_TF_K0s_ALam, V0chi2d12_diff_TF_K0s_ALam, Ntk_Vz_weight);

hbt_reco_truefake(GoodTrackFourVectorF_K0s_ALam, GoodTrackFourVectord1F_K0s_ALam, GoodTrackFourVectord2F_K0s_ALam, GoodTrackFourVectorT_ALam_K0s, GoodTrackFourVectord1T_ALam_K0s, GoodTrackFourVectord2T_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TF, hpeakd1_K0s_ALam_TF, hpeakd2_K0s_ALam_TF, hpeakd12_K0s_ALam_TF, hside_K0s_ALam_TF, hpeakside_K0s_ALam_TF, hpeak_K0s_ALam_TF_bias, hpeak_K0s_ALam_TF_biasd1, hpeak_K0s_ALam_TF_biasd2, hpeak_K0s_ALam_TF_biasd12, Ntk_Vz_weight);
    
}

if(GoodTrackFourVectorT_K0s_ALam.size()>=1 && GoodTrackFourVectorF_ALam_K0s.size()>=1){

nev_KALTF->Fill(1);

get_chi2plot_TF(GoodTrackFourVectorT_K0s_ALam, chi21T_K0s_ALam, chi22T_K0s_ALam, GoodTrackFourVectorF_ALam_K0s, chi21F_ALam_K0s, chi22F_ALam_K0s, V0chi2d1_diff_TF_K0s_ALam, V0chi2d2_diff_TF_K0s_ALam, V0chi2d12_diff_TF_K0s_ALam, Ntk_Vz_weight);

hbt_reco_truefake(GoodTrackFourVectorT_K0s_ALam, GoodTrackFourVectord1T_K0s_ALam, GoodTrackFourVectord2T_K0s_ALam, GoodTrackFourVectorF_ALam_K0s, GoodTrackFourVectord1F_ALam_K0s, GoodTrackFourVectord2F_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TF, hpeakd1_K0s_ALam_TF, hpeakd2_K0s_ALam_TF, hpeakd12_K0s_ALam_TF,  hside_K0s_ALam_TF, hpeakside_K0s_ALam_TF, hpeak_K0s_ALam_TF_bias, hpeak_K0s_ALam_TF_biasd1, hpeak_K0s_ALam_TF_biasd2, hpeak_K0s_ALam_TF_biasd12, Ntk_Vz_weight);
    
}

if(isMC){

//-------------------True Matched-----------------------

if(GoodTrackFourVectorTM_K0s_ALam.size()>=1 && GoodTrackFourVectorTM_ALam_K0s.size()>=1){

Mass_dep_OS(GoodTrackFourVectorTM_K0s_ALam,GoodTrackFourVectorTM_ALam_K0s,(Double_t) aux_N_tk_offline,hMass_K0sAL_TM,hMass_ALamK_TM,hMass_K0s_ALam_TM, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorTM_K0s_ALam,chi21TM_K0s_ALam,chi22TM_K0s_ALam,GoodTrackFourVectorTM_ALam_K0s,chi21TM_ALam_K0s,chi22TM_ALam_K0s,V0chi2d1_diff_TM_K0s_ALam,V0chi2d2_diff_TM_K0s_ALam,V0chi2d12_diff_TM_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTM_K0s_ALam, GoodTrackFourVectord1TM_K0s_ALam, GoodTrackFourVectord2TM_K0s_ALam, GoodTrackFourVectorTM_ALam_K0s, GoodTrackFourVectord1TM_ALam_K0s, GoodTrackFourVectord2TM_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TM, hpeakd1_K0s_ALam_TM, hpeakd2_K0s_ALam_TM, hpeakd12_K0s_ALam_TM, hpeak_K0s_ALam_TM_bias, hpeak_K0s_ALam_TM_biasd1, hpeak_K0s_ALam_TM_biasd2, hpeak_K0s_ALam_TM_biasd12, Ntk_Vz_weight);

ev_z_vtxKAL_TM.push_back(aux_vtxz);
ev_ntrkoff_vecKAL_TM.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecKAL_TM.push_back(GoodTrackFourVectorTM_K0s_ALam); 
ev_GoodTrackFourVector_vecALK_TM.push_back(GoodTrackFourVectorTM_ALam_K0s); 


}

//-------------------True Unmatched-----------------------

if(GoodTrackFourVectorTU_K0s_ALam.size()>=1 && GoodTrackFourVectorTU_ALam_K0s.size()>=1){

Mass_dep_OS(GoodTrackFourVectorTU_K0s_ALam,GoodTrackFourVectorTU_ALam_K0s,(Double_t) aux_N_tk_offline,hMass_K0sAL_TU,hMass_ALamK_TU,hMass_K0s_ALam_TU, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorTU_K0s_ALam,chi21TU_K0s_ALam,chi22TU_K0s_ALam,GoodTrackFourVectorTU_ALam_K0s,chi21TU_ALam_K0s,chi22TU_ALam_K0s,V0chi2d1_diff_TU_K0s_ALam,V0chi2d2_diff_TU_K0s_ALam,V0chi2d12_diff_TU_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTU_K0s_ALam, GoodTrackFourVectord1TU_K0s_ALam, GoodTrackFourVectord2TU_K0s_ALam, GoodTrackFourVectorTU_ALam_K0s, GoodTrackFourVectord1TU_ALam_K0s, GoodTrackFourVectord2TU_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TU, hpeakd1_K0s_ALam_TU, hpeakd2_K0s_ALam_TU, hpeakd12_K0s_ALam_TU, hpeak_K0s_ALam_TU_bias, hpeak_K0s_ALam_TU_biasd1, hpeak_K0s_ALam_TU_biasd2, hpeak_K0s_ALam_TU_biasd12, Ntk_Vz_weight);

ev_z_vtxKAL_TU.push_back(aux_vtxz);
ev_ntrkoff_vecKAL_TU.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecKAL_TU.push_back(GoodTrackFourVectorTU_K0s_ALam); 
ev_GoodTrackFourVector_vecALK_TU.push_back(GoodTrackFourVectorTU_ALam_K0s); 


}

//-------------------Fake Matched-----------------------

if(GoodTrackFourVectorFM_K0s_ALam.size()>=1 && GoodTrackFourVectorFM_ALam_K0s.size()>=1){

Mass_dep_OS(GoodTrackFourVectorFM_K0s_ALam,GoodTrackFourVectorFM_ALam_K0s,(Double_t) aux_N_tk_offline,hMass_K0sAL_FM,hMass_ALamK_FM,hMass_K0s_ALam_FM, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorFM_K0s_ALam,chi21FM_K0s_ALam,chi22FM_K0s_ALam,GoodTrackFourVectorFM_ALam_K0s,chi21FM_ALam_K0s,chi22FM_ALam_K0s,V0chi2d1_diff_FM_K0s_ALam,V0chi2d2_diff_FM_K0s_ALam,V0chi2d12_diff_FM_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFM_K0s_ALam, GoodTrackFourVectord1FM_K0s_ALam, GoodTrackFourVectord2FM_K0s_ALam, GoodTrackFourVectorFM_ALam_K0s, GoodTrackFourVectord1FM_ALam_K0s, GoodTrackFourVectord2FM_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_FM, hpeakd1_K0s_ALam_FM, hpeakd2_K0s_ALam_FM, hpeakd12_K0s_ALam_FM, hpeak_K0s_ALam_FM_bias, hpeak_K0s_ALam_FM_biasd1, hpeak_K0s_ALam_FM_biasd2, hpeak_K0s_ALam_FM_biasd12, Ntk_Vz_weight);

}

//-------------------True Unmatched-----------------------

if(GoodTrackFourVectorFU_K0s_ALam.size()>=1 && GoodTrackFourVectorFU_ALam_K0s.size()>=1){

Mass_dep_OS(GoodTrackFourVectorFU_K0s_ALam,GoodTrackFourVectorFU_ALam_K0s,(Double_t) aux_N_tk_offline,hMass_K0sAL_FU,hMass_ALamK_FU,hMass_K0s_ALam_FU, Ntk_Vz_weight);

get_chi2plot_TF(GoodTrackFourVectorFU_K0s_ALam,chi21FU_K0s_ALam,chi22FU_K0s_ALam,GoodTrackFourVectorFU_ALam_K0s,chi21FU_ALam_K0s,chi22FU_ALam_K0s,V0chi2d1_diff_FU_K0s_ALam,V0chi2d2_diff_FU_K0s_ALam,V0chi2d12_diff_FU_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFU_K0s_ALam, GoodTrackFourVectord1FU_K0s_ALam, GoodTrackFourVectord2FU_K0s_ALam, GoodTrackFourVectorFU_ALam_K0s, GoodTrackFourVectord1FU_ALam_K0s, GoodTrackFourVectord2FU_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_FU, hpeakd1_K0s_ALam_FU, hpeakd2_K0s_ALam_FU, hpeakd12_K0s_ALam_FU, hpeak_K0s_ALam_FU_bias, hpeak_K0s_ALam_FU_biasd1, hpeak_K0s_ALam_FU_biasd2, hpeak_K0s_ALam_FU_biasd12, Ntk_Vz_weight);

}

//-------------------True matched + True unmatched-----------------------

if(GoodTrackFourVectorTM_K0s_ALam.size()>=1 && GoodTrackFourVectorTU_ALam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_K0s_ALam,chi21TM_K0s_ALam,chi22TM_K0s_ALam,GoodTrackFourVectorTU_ALam_K0s,chi21TU_ALam_K0s,chi22TU_ALam_K0s,V0chi2d1_diff_TM_TU_K0s_ALam,V0chi2d2_diff_TM_TU_K0s_ALam,V0chi2d12_diff_TM_TU_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTM_K0s_ALam, GoodTrackFourVectord1TM_K0s_ALam, GoodTrackFourVectord2TM_K0s_ALam, GoodTrackFourVectorTU_ALam_K0s, GoodTrackFourVectord1TU_ALam_K0s, GoodTrackFourVectord2TU_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TM_TU, hpeakd1_K0s_ALam_TM_TU, hpeakd2_K0s_ALam_TM_TU, hpeakd12_K0s_ALam_TM_TU, hpeak_K0s_ALam_TM_TU_bias, hpeak_K0s_ALam_TM_TU_biasd1, hpeak_K0s_ALam_TM_TU_biasd2, hpeak_K0s_ALam_TM_TU_biasd12, Ntk_Vz_weight);

ev_z_vtxKAL_TM_TU.push_back(aux_vtxz);
ev_ntrkoff_vecKAL_TM_TU.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecKAL_TM_TU.push_back(GoodTrackFourVectorTM_K0s_ALam); 
ev_GoodTrackFourVector_vecALK_TM_TU.push_back(GoodTrackFourVectorTU_ALam_K0s); 


}

if(GoodTrackFourVectorTU_K0s_ALam.size()>=1 && GoodTrackFourVectorTM_ALam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_K0s_ALam,chi21TU_K0s_ALam,chi22TU_K0s_ALam,GoodTrackFourVectorTM_ALam_K0s,chi21TM_ALam_K0s,chi22TM_ALam_K0s,V0chi2d1_diff_TM_TU_K0s_ALam,V0chi2d2_diff_TM_TU_K0s_ALam,V0chi2d12_diff_TM_TU_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTU_K0s_ALam, GoodTrackFourVectord1TU_K0s_ALam, GoodTrackFourVectord2TU_K0s_ALam, GoodTrackFourVectorTM_ALam_K0s, GoodTrackFourVectord1TM_ALam_K0s, GoodTrackFourVectord2TM_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TM_TU, hpeakd1_K0s_ALam_TM_TU, hpeakd2_K0s_ALam_TM_TU, hpeakd12_K0s_ALam_TM_TU, hpeak_K0s_ALam_TM_TU_bias, hpeak_K0s_ALam_TM_TU_biasd1, hpeak_K0s_ALam_TM_TU_biasd2, hpeak_K0s_ALam_TM_TU_biasd12, Ntk_Vz_weight);

ev_z_vtxKAL_TU_TM.push_back(aux_vtxz);
ev_ntrkoff_vecKAL_TU_TM.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecKAL_TU_TM.push_back(GoodTrackFourVectorTU_K0s_ALam); 
ev_GoodTrackFourVector_vecALK_TU_TM.push_back(GoodTrackFourVectorTM_ALam_K0s); 


}

//-------------------Fake matched + Fake unmatched-----------------------

if(GoodTrackFourVectorFM_K0s_ALam.size()>=1 && GoodTrackFourVectorFU_ALam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFM_K0s_ALam,chi21FM_K0s_ALam,chi22FM_K0s_ALam,GoodTrackFourVectorFU_ALam_K0s,chi21FU_ALam_K0s,chi22FU_ALam_K0s,V0chi2d1_diff_FM_FU_K0s_ALam,V0chi2d2_diff_FM_FU_K0s_ALam,V0chi2d12_diff_FM_FU_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFM_K0s_ALam, GoodTrackFourVectord1FM_K0s_ALam, GoodTrackFourVectord2FM_K0s_ALam, GoodTrackFourVectorFU_ALam_K0s, GoodTrackFourVectord1FU_ALam_K0s, GoodTrackFourVectord2FU_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_FM_FU, hpeakd1_K0s_ALam_FM_FU, hpeakd2_K0s_ALam_FM_FU, hpeakd12_K0s_ALam_FM_FU, hpeak_K0s_ALam_FM_FU_bias, hpeak_K0s_ALam_FM_FU_biasd1, hpeak_K0s_ALam_FM_FU_biasd2, hpeak_K0s_ALam_FM_FU_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFU_K0s_ALam.size()>=1 && GoodTrackFourVectorFM_ALam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFU_K0s_ALam,chi21FU_K0s_ALam,chi22FU_K0s_ALam,GoodTrackFourVectorFM_ALam_K0s,chi21FM_ALam_K0s,chi22FM_ALam_K0s,V0chi2d1_diff_FM_FU_K0s_ALam,V0chi2d2_diff_FM_FU_K0s_ALam,V0chi2d12_diff_FM_FU_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFU_K0s_ALam, GoodTrackFourVectord1FU_K0s_ALam, GoodTrackFourVectord2FU_K0s_ALam, GoodTrackFourVectorFM_ALam_K0s, GoodTrackFourVectord1FM_ALam_K0s, GoodTrackFourVectord2FM_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_FM_FU, hpeakd1_K0s_ALam_FM_FU, hpeakd2_K0s_ALam_FM_FU, hpeakd12_K0s_ALam_FM_FU, hpeak_K0s_ALam_FM_FU_bias, hpeak_K0s_ALam_FM_FU_biasd1, hpeak_K0s_ALam_FM_FU_biasd2, hpeak_K0s_ALam_FM_FU_biasd12, Ntk_Vz_weight);

}

//-------------------True matched + Fake unmatched-----------------------

if(GoodTrackFourVectorTM_K0s_ALam.size()>=1 && GoodTrackFourVectorFU_ALam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_K0s_ALam,chi21TM_K0s_ALam,chi22TM_K0s_ALam,GoodTrackFourVectorFU_ALam_K0s,chi21FU_ALam_K0s,chi22FU_ALam_K0s,V0chi2d1_diff_TM_FU_K0s_ALam,V0chi2d2_diff_TM_FU_K0s_ALam,V0chi2d12_diff_TM_FU_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTM_K0s_ALam, GoodTrackFourVectord1TM_K0s_ALam, GoodTrackFourVectord2TM_K0s_ALam, GoodTrackFourVectorFU_ALam_K0s, GoodTrackFourVectord1FU_ALam_K0s, GoodTrackFourVectord2FU_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TM_FU, hpeakd1_K0s_ALam_TM_FU, hpeakd2_K0s_ALam_TM_FU, hpeakd12_K0s_ALam_TM_FU, hpeak_K0s_ALam_TM_FU_bias, hpeak_K0s_ALam_TM_FU_biasd1, hpeak_K0s_ALam_TM_FU_biasd2, hpeak_K0s_ALam_TM_FU_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFU_K0s_ALam.size()>=1 && GoodTrackFourVectorTM_ALam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFU_K0s_ALam,chi21FU_K0s_ALam,chi22FU_K0s_ALam,GoodTrackFourVectorTM_ALam_K0s,chi21TM_ALam_K0s,chi22TM_ALam_K0s,V0chi2d1_diff_TM_FU_K0s_ALam,V0chi2d2_diff_TM_FU_K0s_ALam,V0chi2d12_diff_TM_FU_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFU_K0s_ALam, GoodTrackFourVectord1FU_K0s_ALam, GoodTrackFourVectord2FU_K0s_ALam, GoodTrackFourVectorTM_ALam_K0s, GoodTrackFourVectord1TM_ALam_K0s, GoodTrackFourVectord2TM_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TM_FU, hpeakd1_K0s_ALam_TM_FU, hpeakd2_K0s_ALam_TM_FU, hpeakd12_K0s_ALam_TM_FU, hpeak_K0s_ALam_TM_FU_bias, hpeak_K0s_ALam_TM_FU_biasd1, hpeak_K0s_ALam_TM_FU_biasd2, hpeak_K0s_ALam_TM_FU_biasd12, Ntk_Vz_weight);

}

//-------------------True matched + Fake matched-----------------------

if(GoodTrackFourVectorTM_K0s_ALam.size()>=1 && GoodTrackFourVectorFM_ALam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTM_K0s_ALam,chi21TM_K0s_ALam,chi22TM_K0s_ALam,GoodTrackFourVectorFM_ALam_K0s,chi21FM_ALam_K0s,chi22FM_ALam_K0s,V0chi2d1_diff_TM_FM_K0s_ALam,V0chi2d2_diff_TM_FM_K0s_ALam,V0chi2d12_diff_TM_FM_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTM_K0s_ALam, GoodTrackFourVectord1TM_K0s_ALam, GoodTrackFourVectord2TM_K0s_ALam, GoodTrackFourVectorFM_ALam_K0s, GoodTrackFourVectord1FM_ALam_K0s, GoodTrackFourVectord2FM_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TM_FM, hpeakd1_K0s_ALam_TM_FM, hpeakd2_K0s_ALam_TM_FM, hpeakd12_K0s_ALam_TM_FM, hpeak_K0s_ALam_TM_FM_bias, hpeak_K0s_ALam_TM_FM_biasd1, hpeak_K0s_ALam_TM_FM_biasd2, hpeak_K0s_ALam_TM_FM_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFM_K0s_ALam.size()>=1 && GoodTrackFourVectorTM_ALam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFM_K0s_ALam,chi21FM_K0s_ALam,chi22FM_K0s_ALam,GoodTrackFourVectorTM_ALam_K0s,chi21TM_ALam_K0s,chi22TM_ALam_K0s,V0chi2d1_diff_TM_FM_K0s_ALam,V0chi2d2_diff_TM_FM_K0s_ALam,V0chi2d12_diff_TM_FM_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFM_K0s_ALam, GoodTrackFourVectord1FM_K0s_ALam, GoodTrackFourVectord2FM_K0s_ALam, GoodTrackFourVectorTM_ALam_K0s, GoodTrackFourVectord1TM_ALam_K0s, GoodTrackFourVectord2TM_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TM_FM, hpeakd1_K0s_ALam_TM_FM, hpeakd2_K0s_ALam_TM_FM, hpeakd12_K0s_ALam_TM_FM, hpeak_K0s_ALam_TM_FM_bias, hpeak_K0s_ALam_TM_FM_biasd1, hpeak_K0s_ALam_TM_FM_biasd2, hpeak_K0s_ALam_TM_FM_biasd12, Ntk_Vz_weight);

}

//-------------------True unmatched + Fake unmatched-----------------------

if(GoodTrackFourVectorTU_K0s_ALam.size()>=1 && GoodTrackFourVectorFU_ALam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_K0s_ALam,chi21TU_K0s_ALam,chi22TU_K0s_ALam,GoodTrackFourVectorFU_ALam_K0s,chi21FU_ALam_K0s,chi22FU_ALam_K0s,V0chi2d1_diff_TU_FU_K0s_ALam,V0chi2d2_diff_TU_FU_K0s_ALam,V0chi2d12_diff_TU_FU_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTU_K0s_ALam, GoodTrackFourVectord1TU_K0s_ALam, GoodTrackFourVectord2TU_K0s_ALam, GoodTrackFourVectorFU_ALam_K0s, GoodTrackFourVectord1FU_ALam_K0s, GoodTrackFourVectord2FU_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TU_FU, hpeakd1_K0s_ALam_TU_FU, hpeakd2_K0s_ALam_TU_FU, hpeakd12_K0s_ALam_TU_FU, hpeak_K0s_ALam_TU_FU_bias, hpeak_K0s_ALam_TU_FU_biasd1, hpeak_K0s_ALam_TU_FU_biasd2, hpeak_K0s_ALam_TU_FU_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFU_K0s_ALam.size()>=1 && GoodTrackFourVectorTU_ALam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFU_K0s_ALam,chi21FU_K0s_ALam,chi22FU_K0s_ALam,GoodTrackFourVectorTU_ALam_K0s,chi21TU_ALam_K0s,chi22TU_ALam_K0s,V0chi2d1_diff_TU_FU_K0s_ALam,V0chi2d2_diff_TU_FU_K0s_ALam,V0chi2d12_diff_TU_FU_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFU_K0s_ALam, GoodTrackFourVectord1FU_K0s_ALam, GoodTrackFourVectord2FU_K0s_ALam, GoodTrackFourVectorTU_ALam_K0s, GoodTrackFourVectord1TU_ALam_K0s, GoodTrackFourVectord2TU_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TU_FU, hpeakd1_K0s_ALam_TU_FU, hpeakd2_K0s_ALam_TU_FU, hpeakd12_K0s_ALam_TU_FU, hpeak_K0s_ALam_TU_FU_bias, hpeak_K0s_ALam_TU_FU_biasd1, hpeak_K0s_ALam_TU_FU_biasd2, hpeak_K0s_ALam_TU_FU_biasd12, Ntk_Vz_weight);

}

//-------------------True unmatched + Fake matched-----------------------

if(GoodTrackFourVectorTU_K0s_ALam.size()>=1 && GoodTrackFourVectorFM_ALam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorTU_K0s_ALam,chi21TU_K0s_ALam,chi22TU_K0s_ALam,GoodTrackFourVectorFM_ALam_K0s,chi21FM_ALam_K0s,chi22FM_ALam_K0s,V0chi2d1_diff_TU_FM_K0s_ALam,V0chi2d2_diff_TU_FM_K0s_ALam,V0chi2d12_diff_TU_FM_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorTU_K0s_ALam, GoodTrackFourVectord1TU_K0s_ALam, GoodTrackFourVectord2TU_K0s_ALam, GoodTrackFourVectorFM_ALam_K0s, GoodTrackFourVectord1FM_ALam_K0s, GoodTrackFourVectord2FM_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TU_FM, hpeakd1_K0s_ALam_TU_FM, hpeakd2_K0s_ALam_TU_FM, hpeakd12_K0s_ALam_TU_FM, hpeak_K0s_ALam_TU_FM_bias, hpeak_K0s_ALam_TU_FM_biasd1, hpeak_K0s_ALam_TU_FM_biasd2, hpeak_K0s_ALam_TU_FM_biasd12, Ntk_Vz_weight);

}

if(GoodTrackFourVectorFM_K0s_ALam.size()>=1 && GoodTrackFourVectorTU_ALam_K0s.size()>=1){

get_chi2plot_TF(GoodTrackFourVectorFM_K0s_ALam,chi21FM_K0s_ALam,chi22FM_K0s_ALam,GoodTrackFourVectorTU_ALam_K0s,chi21TU_ALam_K0s,chi22TU_ALam_K0s,V0chi2d1_diff_TU_FM_K0s_ALam,V0chi2d2_diff_TU_FM_K0s_ALam,V0chi2d12_diff_TU_FM_K0s_ALam, Ntk_Vz_weight);
    
hbt_MM_UU_cross(GoodTrackFourVectorFM_K0s_ALam, GoodTrackFourVectord1FM_K0s_ALam, GoodTrackFourVectord2FM_K0s_ALam, GoodTrackFourVectorTU_ALam_K0s, GoodTrackFourVectord1TU_ALam_K0s, GoodTrackFourVectord2TU_ALam_K0s, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TU_FM, hpeakd1_K0s_ALam_TU_FM, hpeakd2_K0s_ALam_TU_FM, hpeakd12_K0s_ALam_TU_FM, hpeak_K0s_ALam_TU_FM_bias, hpeak_K0s_ALam_TU_FM_biasd1, hpeak_K0s_ALam_TU_FM_biasd2, hpeak_K0s_ALam_TU_FM_biasd12, Ntk_Vz_weight);

}

}

//for Feed down studies

if(DoFeedDown){


//LL
if(GoodTrackFourVectorT_Lam_FD_L.size()>=2){
hbt_simples(GoodTrackFourVectorT_Lam_FD_L, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, hpeak_LamFD_LamFD, Ntk_Vz_weight);
//Mix Random
weight_LFDLFD.push_back(Ntk_Vz_weight);
ev_z_vtx_LamFD_LamFD.push_back(aux_vtxz);
ev_ntrkoff_vec_LamFD_LamFD.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vec_LamFD_LamFD.push_back(GoodTrackFourVectorT_Lam_FD_L); 
}

if(GoodTrackFourVectorT_Lam.size()>=1 && GoodTrackFourVectorT_Lam_FD_L.size()>=1){
hbt_simples_cross(GoodTrackFourVectorT_Lam, GoodTrackFourVectorT_Lam_FD_L, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, hpeak_Lam_LamFD, Ntk_Vz_weight);
weight_LLFD.push_back(Ntk_Vz_weight);
ev_z_vtx_Lam_LamFD.push_back(aux_vtxz);
ev_ntrkoff_vec_Lam_LamFD.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vec_Lam_LamFD.push_back(GoodTrackFourVectorT_Lam);
ev_GoodTrackFourVector_vec_LamFD_Lam.push_back(GoodTrackFourVectorT_Lam_FD_L);
//get_chi2plot_TF(GoodTrackFourVectorT_Lam, chi21T_K0s_ALam, chi22T_K0s_ALam, GoodTrackFourVectorT_Lam_FD, chi21F_ALam_K0s, chi22F_ALam_K0s, V0chi2d1_diff_TF_K0s_ALam, V0chi2d2_diff_TF_K0s_ALam, V0chi2d12_diff_TF_K0s_ALam, Ntk_Vz_weight);
}

//ALAL
if(GoodTrackFourVectorT_ALam_FD_AL.size()>=2){
hbt_simples(GoodTrackFourVectorT_ALam_FD_AL, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, nsigmapeak, hpeak_ALamFD_ALamFD, Ntk_Vz_weight);
//Mix Random
weight_ALFDALFD.push_back(Ntk_Vz_weight);
ev_z_vtx_ALamFD_ALamFD.push_back(aux_vtxz);
ev_ntrkoff_vec_ALamFD_ALamFD.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vec_ALamFD_ALamFD.push_back(GoodTrackFourVectorT_ALam_FD_AL); 
}

if(GoodTrackFourVectorT_ALam.size()>=1 && GoodTrackFourVectorT_ALam_FD_AL.size()>=1){
hbt_simples_cross(GoodTrackFourVectorT_ALam, GoodTrackFourVectorT_ALam_FD_AL, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, hpeak_ALam_ALamFD, Ntk_Vz_weight);
weight_ALALFD.push_back(Ntk_Vz_weight);
ev_z_vtx_ALam_ALamFD.push_back(aux_vtxz);
ev_ntrkoff_vec_ALam_ALamFD.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vec_ALam_ALamFD.push_back(GoodTrackFourVectorT_ALam);
ev_GoodTrackFourVector_vec_ALamFD_ALam.push_back(GoodTrackFourVectorT_ALam_FD_AL);
}


//LAL

if(GoodTrackFourVectorT_Lam.size()>=1 && GoodTrackFourVectorT_ALam_FD_L.size()>=1){
hbt_simples_cross(GoodTrackFourVectorT_Lam, GoodTrackFourVectorT_ALam_FD_L, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, hpeak_LALFD, Ntk_Vz_weight);
weight_LALFD.push_back(Ntk_Vz_weight);
ev_z_vtx_LALFD.push_back(aux_vtxz);
ev_ntrkoff_vec_LALFD.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vec_LALFD.push_back(GoodTrackFourVectorT_Lam);
ev_GoodTrackFourVector_vec_LFDAL.push_back(GoodTrackFourVectorT_ALam_FD_L);
}

if(GoodTrackFourVectorT_Lam_FD_AL.size()>=1 && GoodTrackFourVectorT_ALam_FD_L.size()>=1){
hbt_simples_cross(GoodTrackFourVectorT_Lam_FD_AL, GoodTrackFourVectorT_ALam_FD_L, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, hpeak_LFDALFD, Ntk_Vz_weight);
weight_LFDALFD.push_back(Ntk_Vz_weight);
ev_z_vtx_LFDALFD.push_back(aux_vtxz);
ev_ntrkoff_vec_LFDALFD.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vec_LFDALFD.push_back(GoodTrackFourVectorT_Lam_FD_AL);
ev_GoodTrackFourVector_vec_ALFDLFD.push_back(GoodTrackFourVectorT_ALam_FD_L);
}

if(GoodTrackFourVectorT_ALam.size()>=1 && GoodTrackFourVectorT_Lam_FD_AL.size()>=1){
hbt_simples_cross(GoodTrackFourVectorT_ALam, GoodTrackFourVectorT_Lam_FD_AL, (Double_t) aux_N_tk_offline, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, hpeak_LALFD, Ntk_Vz_weight);
weight_ALLFD.push_back(Ntk_Vz_weight);
ev_z_vtx_ALLFD.push_back(aux_vtxz);
ev_ntrkoff_vec_ALLFD.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vec_ALLFD.push_back(GoodTrackFourVectorT_ALam);
ev_GoodTrackFourVector_vec_ALFDL.push_back(GoodTrackFourVectorT_Lam_FD_AL);
}


//KL

if(GoodTrackFourVectorT_K0s.size()>=1 && GoodTrackFourVectorT_Lam_FD_K.size()>=1){
hbt_simples_cross(GoodTrackFourVectorT_K0s, GoodTrackFourVectorT_Lam_FD_K, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, hpeak_KLFD, Ntk_Vz_weight);
weight_KLFD.push_back(Ntk_Vz_weight);
ev_z_vtx_KLFD.push_back(aux_vtxz);
ev_ntrkoff_vec_KLFD.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vec_KLFD.push_back(GoodTrackFourVectorT_K0s);
ev_GoodTrackFourVector_vec_LFDK.push_back(GoodTrackFourVectorT_Lam_FD_K);
}

//KAL

if(GoodTrackFourVectorT_K0s.size()>=1 && GoodTrackFourVectorT_ALam_FD_K.size()>=1){
hbt_simples_cross(GoodTrackFourVectorT_K0s, GoodTrackFourVectorT_ALam_FD_K, (Double_t) aux_N_tk_offline, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, hpeak_KALFD, Ntk_Vz_weight);
weight_KALFD.push_back(Ntk_Vz_weight);
ev_z_vtx_KALFD.push_back(aux_vtxz);
ev_ntrkoff_vec_KALFD.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vec_KALFD.push_back(GoodTrackFourVectorT_K0s);
ev_GoodTrackFourVector_vec_ALFDK.push_back(GoodTrackFourVectorT_ALam_FD_K);
}

}

if(DoFeedDown_Xi){


get_chi2plot_TF_new(chi21_Lam, chi22_Lam, GoodTrackFourVector_Xi_Lamb_proton_chi2_T_L, GoodTrackFourVector_Xi_Lamb_pion_chi2_T_L,V0_chi2_SS_LLFD_XI,V0_chi2_OS_LLFD_XI,Ntk_Vz_weight);
get_chi2plot_TF_new(chi21_ALam, chi22_ALam, GoodTrackFourVector_Xi_Lamb_proton_chi2_T_AL, GoodTrackFourVector_Xi_Lamb_pion_chi2_T_AL,V0_chi2_SS_ALFD_XI,V0_chi2_OS_ALFD_XI,Ntk_Vz_weight);
get_chi2plot_TF_new(chi21_K0s, chi22_K0s, GoodTrackFourVector_Xi_Lamb_proton_chi2_T_K, GoodTrackFourVector_Xi_Lamb_pion_chi2_T_K,V0_chi2_SS_KLFD_XI,V0_chi2_OS_KLFD_XI,Ntk_Vz_weight);

get_chi2plot_TF_new(chi21_Lam, chi22_Lam, GoodTrackFourVector_AXi_Lamb_proton_chi2_T_L, GoodTrackFourVector_AXi_Lamb_pion_chi2_T_L,V0_chi2_SS_LLFD_AXI,V0_chi2_OS_LLFD_AXI,Ntk_Vz_weight);
get_chi2plot_TF_new(chi21_ALam, chi22_ALam, GoodTrackFourVector_AXi_Lamb_proton_chi2_T_AL, GoodTrackFourVector_AXi_Lamb_pion_chi2_T_AL,V0_chi2_SS_ALFD_AXI,V0_chi2_OS_ALFD_AXI,Ntk_Vz_weight);
get_chi2plot_TF_new(chi21_K0s, chi22_K0s, GoodTrackFourVector_AXi_Lamb_proton_chi2_T_K, GoodTrackFourVector_AXi_Lamb_pion_chi2_T_K,V0_chi2_SS_KLFD_AXI,V0_chi2_OS_KLFD_AXI,Ntk_Vz_weight);


}

if(DoFeedDown_Om){

get_chi2plot_TF_new(chi21_Lam, chi22_Lam, GoodTrackFourVector_Om_Lamb_proton_chi2_T_L, GoodTrackFourVector_Om_Lamb_kaon_chi2_T_L,V0_chi2_SS_LLFD_OM,V0_chi2_OS_LLFD_OM,Ntk_Vz_weight);
get_chi2plot_TF_new(chi21_ALam, chi22_ALam, GoodTrackFourVector_Om_Lamb_proton_chi2_T_AL, GoodTrackFourVector_Om_Lamb_kaon_chi2_T_AL,V0_chi2_SS_ALFD_OM,V0_chi2_OS_ALFD_OM,Ntk_Vz_weight);
get_chi2plot_TF_new(chi21_K0s, chi22_K0s, GoodTrackFourVector_Om_Lamb_proton_chi2_T_K, GoodTrackFourVector_Om_Lamb_kaon_chi2_T_K,V0_chi2_SS_KLFD_OM,V0_chi2_OS_KLFD_OM,Ntk_Vz_weight);

get_chi2plot_TF_new(chi21_Lam, chi22_Lam, GoodTrackFourVector_AOm_Lamb_proton_chi2_T_L, GoodTrackFourVector_AOm_Lamb_kaon_chi2_T_L,V0_chi2_SS_LLFD_AOM,V0_chi2_OS_LLFD_AOM,Ntk_Vz_weight);
get_chi2plot_TF_new(chi21_ALam, chi22_ALam, GoodTrackFourVector_AOm_Lamb_proton_chi2_T_AL, GoodTrackFourVector_AOm_Lamb_kaon_chi2_T_AL,V0_chi2_SS_ALFD_AOM,V0_chi2_OS_ALFD_AOM,Ntk_Vz_weight);
get_chi2plot_TF_new(chi21_K0s, chi22_K0s, GoodTrackFourVector_AOm_Lamb_proton_chi2_T_K, GoodTrackFourVector_AOm_Lamb_kaon_chi2_T_K,V0_chi2_SS_KLFD_AOM,V0_chi2_OS_KLFD_AOM,Ntk_Vz_weight);


}

if(GoodTrackFourVectorT_K0s.size()>0){nev_K0s_AS->Fill(1);h_ntrk_cent_1Ks->Fill((Double_t) aux_N_tk_offline,Ntk_Vz_weight);}
if(GoodTrackFourVectorT_Lam.size()>0){nev_Lam_AS->Fill(1);h_ntrk_cent_1Lam->Fill((Double_t) aux_N_tk_offline,Ntk_Vz_weight);}
if(GoodTrackFourVectorT_ALam.size()>0){nev_ALam_AS->Fill(1);h_ntrk_cent_1ALam->Fill((Double_t) aux_N_tk_offline,Ntk_Vz_weight);}

//////////////////////////////////control - Data//////////////////////////////////

K_size->Fill(GoodTrackFourVectorT_K0s.size(),Ntk_Vz_weight);
L_size->Fill(GoodTrackFourVectorT_Lam.size(),Ntk_Vz_weight);
AL_size->Fill(GoodTrackFourVectorT_ALam.size(),Ntk_Vz_weight);


cplots(GoodTrackFourVectorT_K0s,K0s_pt_reco,K0s_eta_reco,K0s_phi_reco,K0s_mass_reco,Ntk_Vz_weight);
cplots(GoodTrackFourVectorT_Lam,L_pt_reco,L_eta_reco,L_phi_reco,L_mass_reco,Ntk_Vz_weight);
cplots(GoodTrackFourVectorT_ALam,AL_pt_reco,AL_eta_reco,AL_phi_reco,AL_mass_reco,Ntk_Vz_weight);

cplots(GoodTrackFourVectord1T_K0s,K0s_pt_reco_d1,K0s_eta_reco_d1,K0s_phi_reco_d1,K0s_mass_reco_d1,Ntk_Vz_weight);
cplots(GoodTrackFourVectord1T_Lam,L_pt_reco_d1,L_eta_reco_d1,L_phi_reco_d1,L_mass_reco_d1,Ntk_Vz_weight);
cplots(GoodTrackFourVectord1T_ALam,AL_pt_reco_d1,AL_eta_reco_d1,AL_phi_reco_d1,AL_mass_reco_d1,Ntk_Vz_weight);

cplots(GoodTrackFourVectord2T_K0s,K0s_pt_reco_d2,K0s_eta_reco_d2,K0s_phi_reco_d2,K0s_mass_reco_d2,Ntk_Vz_weight);
cplots(GoodTrackFourVectord2T_Lam,L_pt_reco_d2,L_eta_reco_d2,L_phi_reco_d2,L_mass_reco_d2,Ntk_Vz_weight);
cplots(GoodTrackFourVectord2T_ALam,AL_pt_reco_d2,AL_eta_reco_d2,AL_phi_reco_d2,AL_mass_reco_d2,Ntk_Vz_weight);


///////////////////////////////////////////////////////////////////////////////



if(isMC){

//////////////////////////////////control//////////////////////////////////

cplots(GoodTrackFourVectorTM_K0s,K0s_pt_mat,K0s_eta_mat,K0s_phi_mat,K0s_mass_mat,Ntk_Vz_weight);
cplots(GoodTrackFourVectorTM_Lam,L_pt_mat,L_eta_mat,L_phi_mat,L_mass_mat,Ntk_Vz_weight);
cplots(GoodTrackFourVectorTM_ALam,AL_pt_mat,AL_eta_mat,AL_phi_mat,AL_mass_mat,Ntk_Vz_weight);

///////////////////////////////////////////////////////////////////////////////
}


//include jets

if(isMC){

if(jetsize >= 1){
for(Int_t ii=0; ii<jetsize; ii++){ 
bool issamejet = false;
if(Jet_pt->at(ii) <= ptminJet || Jet_pt->at(ii) >= ptmaxJet)continue;
TLorentzVector pvectorj;
pvectorj.SetPtEtaPhiM(Jet_pt->at(ii),Jet_eta->at(ii),Jet_phi->at(ii),Jet_mass->at(ii));
for(Int_t j=ii+1; j<jetsize; j++){
TLorentzVector pvectorj2;
pvectorj2.SetPtEtaPhiM(Jet_pt->at(j),Jet_eta->at(j),Jet_phi->at(j),Jet_mass->at(j));
double deltaPhi=pvectorj.DeltaPhi(pvectorj2);
double deltaEta=pvectorj.Eta() - pvectorj2.Eta();
double DR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);
cone_Jet->Fill(DR);
if(DR<=0.8){issamejet=true;}
}
if(removedoublejet){if(issamejet)continue;}
GoodTrackFourVector_Jet.push_back(pvectorj);
}}


//jet+jet

if(GoodTrackFourVector_Jet.size()>=2){

hbt_gen_jet(GoodTrackFourVector_Jet, (Double_t) aux_N_tk_offline, h_jet_jet, h_inv_jet_jet, h_rot_jet_jet);

ev_z_vtxJet.push_back(aux_vtxz);
ev_ntrkoff_vecJet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecJet.push_back(GoodTrackFourVector_Jet); 

}

//K0s+Jet

if(GoodTrackFourVectorT_K0s.size() >= 1 && GoodTrackFourVector_Jet.size()>=1){

hbt_gen_jet_OS(GoodTrackFourVectorT_K0s, GoodTrackFourVector_Jet, (Double_t) aux_N_tk_offline, h_K0s_Jet, h_K0s_Jet_inv, h_K0s_Jet_rot, massK0s, sigmaK0s, nsigmapeak);

ev_z_vtxK0sJet.push_back(aux_vtxz);
ev_ntrkoff_vecK0sJet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0sJet.push_back(GoodTrackFourVectorT_K0s);
ev_GoodTrackFourVector_vecJetK0s.push_back(GoodTrackFourVector_Jet);

}

//Lam+Jet

if(GoodTrackFourVectorT_Lam.size() >= 1 && GoodTrackFourVector_Jet.size()>=1){

hbt_gen_jet_OS(GoodTrackFourVectorT_Lam, GoodTrackFourVector_Jet, (Double_t) aux_N_tk_offline, h_Lam_Jet, h_Lam_Jet_inv, h_Lam_Jet_rot, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxLamJet.push_back(aux_vtxz);
ev_ntrkoff_vecLamJet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLamJet.push_back(GoodTrackFourVectorT_Lam);
ev_GoodTrackFourVector_vecJetLam.push_back(GoodTrackFourVector_Jet);

}

//ALam+Jet

if(GoodTrackFourVectorT_ALam.size() >= 1 && GoodTrackFourVector_Jet.size()>=1){

hbt_gen_jet_OS(GoodTrackFourVectorT_ALam, GoodTrackFourVector_Jet, (Double_t) aux_N_tk_offline, h_ALam_Jet, h_ALam_Jet_inv, h_ALam_Jet_rot, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxALamJet.push_back(aux_vtxz);
ev_ntrkoff_vecALamJet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecALamJet.push_back(GoodTrackFourVectorT_ALam);
ev_GoodTrackFourVector_vecJetALam.push_back(GoodTrackFourVector_Jet);

}

//K0Jets+Jet

if(GoodTrackFourVector_K0s_Jet.size() >= 1 && GoodTrackFourVector_Jet.size()>=1){

hbt_gen_jet_OS(GoodTrackFourVector_K0s_Jet, GoodTrackFourVector_Jet, (Double_t) aux_N_tk_offline, h_K0sfromjet_Jet, h_K0sfromjet_Jet_inv, h_K0sfromjet_Jet_rot, massK0s, sigmaK0s, nsigmapeak);

ev_z_vtxK0sfromjetJet.push_back(aux_vtxz);
ev_ntrkoff_vecK0sfromjetJet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0sfromjetJet.push_back(GoodTrackFourVector_K0s_Jet);
ev_GoodTrackFourVector_vecJetK0sfromjet.push_back(GoodTrackFourVector_Jet);

}

//LamJet+Jet

if(GoodTrackFourVector_Lam_Jet.size() >= 1 && GoodTrackFourVector_Jet.size()>=1){

hbt_gen_jet_OS(GoodTrackFourVector_Lam_Jet, GoodTrackFourVector_Jet, (Double_t) aux_N_tk_offline, h_Lamfromjet_Jet, h_Lamfromjet_Jet_inv, h_Lamfromjet_Jet_rot, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxLamfromjetJet.push_back(aux_vtxz);
ev_ntrkoff_vecLamfromjetJet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLamfromjetJet.push_back(GoodTrackFourVector_Lam_Jet);
ev_GoodTrackFourVector_vecJetLamfromjet.push_back(GoodTrackFourVector_Jet);

}

//ALamJet+Jet

if(GoodTrackFourVector_ALam_Jet.size() >= 1 && GoodTrackFourVector_Jet.size()>=1){

hbt_gen_jet_OS(GoodTrackFourVector_ALam_Jet, GoodTrackFourVector_Jet, (Double_t) aux_N_tk_offline, h_ALamfromjet_Jet, h_ALamfromjet_Jet_inv, h_ALamfromjet_Jet_rot, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxALamfromjetJet.push_back(aux_vtxz);
ev_ntrkoff_vecALamfromjetJet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecALamfromjetJet.push_back(GoodTrackFourVector_ALam_Jet);
ev_GoodTrackFourVector_vecJetALamfromjet.push_back(GoodTrackFourVector_Jet);

}

//K0sNotaJet+Jet

if(GoodTrackFourVector_K0s_NotaJet.size() >= 1 && GoodTrackFourVector_Jet.size()>=1){

hbt_gen_jet_OS(GoodTrackFourVector_K0s_NotaJet, GoodTrackFourVector_Jet, (Double_t) aux_N_tk_offline, h_K0snojet_Jet, h_K0snojet_Jet_inv, h_K0snojet_Jet_rot, massK0s, sigmaK0s, nsigmapeak);

ev_z_vtxK0snojetJet.push_back(aux_vtxz);
ev_ntrkoff_vecK0snojetJet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0snojetJet.push_back(GoodTrackFourVector_K0s_NotaJet);
ev_GoodTrackFourVector_vecJetK0snojet.push_back(GoodTrackFourVector_Jet);

}

//LamNotaJet+Jet

if(GoodTrackFourVector_Lam_NotaJet.size() >= 1 && GoodTrackFourVector_Jet.size()>=1){

hbt_gen_jet_OS(GoodTrackFourVector_Lam_NotaJet, GoodTrackFourVector_Jet, (Double_t) aux_N_tk_offline, h_Lamnojet_Jet, h_Lamnojet_Jet_inv, h_Lamnojet_Jet_rot, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxLamnojetJet.push_back(aux_vtxz);
ev_ntrkoff_vecLamnojetJet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLamnojetJet.push_back(GoodTrackFourVector_Lam_NotaJet);
ev_GoodTrackFourVector_vecJetLamnojet.push_back(GoodTrackFourVector_Jet);

}

//ALamNotaJet+Jet

if(GoodTrackFourVector_ALam_NotaJet.size() >= 1 && GoodTrackFourVector_Jet.size()>=1){

hbt_gen_jet_OS(GoodTrackFourVector_ALam_NotaJet, GoodTrackFourVector_Jet, (Double_t) aux_N_tk_offline, h_ALamnojet_Jet, h_ALamnojet_Jet_inv, h_ALamnojet_Jet_rot, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxALamnojetJet.push_back(aux_vtxz);
ev_ntrkoff_vecALamnojetJet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecALamnojetJet.push_back(GoodTrackFourVector_ALam_NotaJet);
ev_GoodTrackFourVector_vecJetALamnojet.push_back(GoodTrackFourVector_Jet);

}


//V0 from jets

//nojet+nojet

if(GoodTrackFourVector_K0s_NotaJet.size()>=2){

hbt_gen_V0jet(GoodTrackFourVector_K0s_NotaJet, (Double_t) aux_N_tk_offline, h_K0s_K0s_nojet, h_inv_K0s_K0s_nojet, h_rot_K0s_K0s_nojet, massK0s, sigmaK0s, nsigmapeak);

ev_z_vtxK0s_nojet.push_back(aux_vtxz);
ev_ntrkoff_vecK0s_nojet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0s_nojet.push_back(GoodTrackFourVector_K0s_NotaJet); 

} //K0sK0s

if(GoodTrackFourVector_Lam_NotaJet.size()>=2){

hbt_gen_V0jet(GoodTrackFourVector_Lam_NotaJet, (Double_t) aux_N_tk_offline, h_Lam_Lam_nojet, h_inv_Lam_Lam_nojet, h_rot_Lam_Lam_nojet, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxLam_nojet.push_back(aux_vtxz);
ev_ntrkoff_vecLam_nojet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLam_nojet.push_back(GoodTrackFourVector_Lam_NotaJet); 

} //LamLam

if(GoodTrackFourVector_ALam_NotaJet.size()>=2){

hbt_gen_V0jet(GoodTrackFourVector_ALam_NotaJet, (Double_t) aux_N_tk_offline, h_ALam_ALam_nojet, h_inv_ALam_ALam_nojet, h_rot_ALam_ALam_nojet, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxALam_nojet.push_back(aux_vtxz);
ev_ntrkoff_vecALam_nojet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecALam_nojet.push_back(GoodTrackFourVector_ALam_NotaJet); 

} //ALamALam


if(GoodTrackFourVector_Lam_NotaJet.size()>=1 && GoodTrackFourVector_ALam_NotaJet.size()>=1){

hbt_gen_V0jet_OS(GoodTrackFourVector_Lam_NotaJet,GoodTrackFourVector_ALam_NotaJet,(Double_t) aux_N_tk_offline,h_Lam_ALam_nojet,h_inv_Lam_ALam_nojet,h_rot_Lam_ALam_nojet, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxLamALam_nojet.push_back(aux_vtxz);
ev_ntrkoff_vecLamALam_nojet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLamALam_nojet.push_back(GoodTrackFourVector_Lam_NotaJet); 
ev_GoodTrackFourVector_vecALamLam_nojet.push_back(GoodTrackFourVector_ALam_NotaJet); 

}

if(GoodTrackFourVector_K0s_NotaJet.size()>=1 && GoodTrackFourVector_Lam_NotaJet.size()>=1){

hbt_gen_V0jet_OS(GoodTrackFourVector_K0s_NotaJet, GoodTrackFourVector_Lam_NotaJet,(Double_t) aux_N_tk_offline,h_K0s_Lam_nojet,h_inv_K0s_Lam_nojet,h_rot_K0s_Lam_nojet, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxK0sLam_nojet.push_back(aux_vtxz);
ev_ntrkoff_vecK0sLam_nojet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0sLam_nojet.push_back(GoodTrackFourVector_K0s_NotaJet); 
ev_GoodTrackFourVector_vecLamK0s_nojet.push_back(GoodTrackFourVector_Lam_NotaJet); 


} //K0sLam

if(GoodTrackFourVector_K0s_NotaJet.size()>=1 && GoodTrackFourVector_ALam_NotaJet.size()>=1){

hbt_gen_V0jet_OS(GoodTrackFourVector_K0s_NotaJet,GoodTrackFourVector_ALam_NotaJet,(Double_t) aux_N_tk_offline,h_K0s_ALam_nojet,h_inv_K0s_ALam_nojet,h_rot_K0s_ALam_nojet, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxK0sALam_nojet.push_back(aux_vtxz);
ev_ntrkoff_vecK0sALam_nojet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0sALam_nojet.push_back(GoodTrackFourVector_K0s_NotaJet); 
ev_GoodTrackFourVector_vecALamK0s_nojet.push_back(GoodTrackFourVector_ALam_NotaJet); 

} //K0sALam

//jet+nojet

if(GoodTrackFourVector_K0s_NotaJet.size()>=1 && GoodTrackFourVector_K0s_Jet.size()>=1){

hbt_gen_V0jet_OS(GoodTrackFourVector_K0s_NotaJet, GoodTrackFourVector_K0s_Jet,(Double_t) aux_N_tk_offline,h_K0s_K0s_only1jet,h_inv_K0s_K0s_only1jet,h_rot_K0s_K0s_only1jet, massK0s, sigmaK0s, massK0s, sigmaK0s, nsigmapeak);

ev_z_vtxK0s_only1jet.push_back(aux_vtxz);
ev_ntrkoff_vecK0s_only1jet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0s_only1jet.push_back(GoodTrackFourVector_K0s_NotaJet); 
ev_GoodTrackFourVector_vecK0s_only1jetX.push_back(GoodTrackFourVector_K0s_Jet); 

}//K0sjet+K0snotajet

if(GoodTrackFourVector_Lam_NotaJet.size()>=1 && GoodTrackFourVector_Lam_Jet.size()>=1){

hbt_gen_V0jet_OS(GoodTrackFourVector_Lam_NotaJet, GoodTrackFourVector_Lam_Jet,(Double_t) aux_N_tk_offline,h_Lam_Lam_only1jet,h_inv_Lam_Lam_only1jet,h_rot_Lam_Lam_only1jet, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxLam_only1jet.push_back(aux_vtxz);
ev_ntrkoff_vecLam_only1jet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLam_only1jet.push_back(GoodTrackFourVector_Lam_NotaJet); 
ev_GoodTrackFourVector_vecLam_only1jetX.push_back(GoodTrackFourVector_Lam_Jet); 

}//Lamjet+Lamnotajet

if(GoodTrackFourVector_ALam_NotaJet.size()>=1 && GoodTrackFourVector_ALam_Jet.size()>=1){

hbt_gen_V0jet_OS(GoodTrackFourVector_ALam_NotaJet, GoodTrackFourVector_ALam_Jet,(Double_t) aux_N_tk_offline,h_ALam_ALam_only1jet,h_inv_ALam_ALam_only1jet,h_rot_ALam_ALam_only1jet, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxALam_only1jet.push_back(aux_vtxz);
ev_ntrkoff_vecALam_only1jet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecALam_only1jet.push_back(GoodTrackFourVector_ALam_NotaJet); 
ev_GoodTrackFourVector_vecALam_only1jetX.push_back(GoodTrackFourVector_ALam_Jet); 

}//ALamjet+ALamnotajet

if(GoodTrackFourVector_Lam_NotaJet.size()>=1 && GoodTrackFourVector_ALam_Jet.size()>=1){

hbt_gen_V0jet_OS(GoodTrackFourVector_Lam_NotaJet, GoodTrackFourVector_ALam_Jet,(Double_t) aux_N_tk_offline,h_Lam_ALam_only1jet,h_inv_Lam_ALam_only1jet,h_rot_Lam_ALam_only1jet, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxLamALam_only1jet.push_back(aux_vtxz);
ev_ntrkoff_vecLamALam_only1jet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLamALam_only1jet.push_back(GoodTrackFourVector_Lam_NotaJet); 
ev_GoodTrackFourVector_vecALamLam_only1jet.push_back(GoodTrackFourVector_ALam_Jet); 

}//ALamjet+Lamnotajet

if(GoodTrackFourVector_ALam_NotaJet.size()>=1 && GoodTrackFourVector_Lam_Jet.size()>=1){

hbt_gen_V0jet_OS(GoodTrackFourVector_ALam_NotaJet, GoodTrackFourVector_Lam_Jet,(Double_t) aux_N_tk_offline,h_Lam_ALam_only1jet,h_inv_Lam_ALam_only1jet,h_rot_Lam_ALam_only1jet, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxLamALam_only1jetX.push_back(aux_vtxz);
ev_ntrkoff_vecLamALam_only1jetX.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLamALam_only1jetX.push_back(GoodTrackFourVector_ALam_NotaJet); 
ev_GoodTrackFourVector_vecALamLam_only1jetX.push_back(GoodTrackFourVector_Lam_Jet); 

}//Lamjet+ALamnotajet

if(GoodTrackFourVector_K0s_NotaJet.size()>=1 && GoodTrackFourVector_Lam_Jet.size()>=1){

hbt_gen_V0jet_OS(GoodTrackFourVector_K0s_NotaJet, GoodTrackFourVector_Lam_Jet,(Double_t) aux_N_tk_offline,h_K0s_Lam_only1jet,h_inv_K0s_Lam_only1jet,h_rot_K0s_Lam_only1jet, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxK0sLam_only1jet.push_back(aux_vtxz);
ev_ntrkoff_vecK0sLam_only1jet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0sLam_only1jet.push_back(GoodTrackFourVector_K0s_NotaJet); 
ev_GoodTrackFourVector_vecLamK0s_only1jet.push_back(GoodTrackFourVector_Lam_Jet); 

}

if(GoodTrackFourVector_Lam_NotaJet.size()>=1 && GoodTrackFourVector_K0s_Jet.size()>=1){

hbt_gen_V0jet_OS(GoodTrackFourVector_Lam_NotaJet, GoodTrackFourVector_K0s_Jet, (Double_t) aux_N_tk_offline,h_K0s_Lam_only1jet,h_inv_K0s_Lam_only1jet,h_rot_K0s_Lam_only1jet, massLAL, sigmaLAL, massK0s, sigmaK0s, nsigmapeak);

ev_z_vtxK0sLam_only1jetX.push_back(aux_vtxz);
ev_ntrkoff_vecK0sLam_only1jetX.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0sLam_only1jetX.push_back(GoodTrackFourVector_Lam_NotaJet); 
ev_GoodTrackFourVector_vecLamK0s_only1jetX.push_back(GoodTrackFourVector_K0s_Jet); 

}

if(GoodTrackFourVector_K0s_NotaJet.size()>=1 && GoodTrackFourVector_ALam_Jet.size()>=1){

hbt_gen_V0jet_OS(GoodTrackFourVector_K0s_NotaJet, GoodTrackFourVector_ALam_Jet,(Double_t) aux_N_tk_offline,h_K0s_ALam_only1jet,h_inv_K0s_ALam_only1jet,h_rot_K0s_ALam_only1jet, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxK0sALam_only1jet.push_back(aux_vtxz);
ev_ntrkoff_vecK0sALam_only1jet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0sALam_only1jet.push_back(GoodTrackFourVector_K0s_NotaJet); 
ev_GoodTrackFourVector_vecALamK0s_only1jet.push_back(GoodTrackFourVector_ALam_Jet); 

}

if(GoodTrackFourVector_ALam_NotaJet.size()>=1 && GoodTrackFourVector_K0s_Jet.size()>=1){

hbt_gen_V0jet_OS(GoodTrackFourVector_ALam_NotaJet, GoodTrackFourVector_K0s_Jet, (Double_t) aux_N_tk_offline,h_K0s_ALam_only1jet,h_inv_K0s_ALam_only1jet,h_rot_K0s_ALam_only1jet, massLAL, sigmaLAL, massK0s, sigmaK0s, nsigmapeak);

ev_z_vtxK0sALam_only1jetX.push_back(aux_vtxz);
ev_ntrkoff_vecK0sALam_only1jetX.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0sALam_only1jetX.push_back(GoodTrackFourVector_ALam_NotaJet); 
ev_GoodTrackFourVector_vecALamK0s_only1jetX.push_back(GoodTrackFourVector_K0s_Jet); 

}


//fromjet+fromjet

if(GoodTrackFourVector_K0s_Jet.size()>=2){

hbt_gen_jet_V0(GoodTrackFourVector_K0s_Jet, GoodTrackFourVector_K0s_NJet, (Double_t) aux_N_tk_offline, h_K0s_K0s_samejet, h_inv_K0s_K0s_samejet, h_rot_K0s_K0s_samejet, true, massK0s, sigmaK0s, nsigmapeak);
hbt_gen_jet_V0(GoodTrackFourVector_K0s_Jet, GoodTrackFourVector_K0s_NJet, (Double_t) aux_N_tk_offline, h_K0s_K0s_diffjet, h_inv_K0s_K0s_diffjet, h_rot_K0s_K0s_diffjet, false, massK0s, sigmaK0s, nsigmapeak);

ev_z_vtxK0s_samediffjet.push_back(aux_vtxz);
ev_ntrkoff_vecK0s_samediffjet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0s_samediffjet.push_back(GoodTrackFourVector_K0s_Jet); 
ev_njet_vecK0s_samediffjet.push_back(GoodTrackFourVector_K0s_NJet);

}

if(GoodTrackFourVector_Lam_Jet.size()>=2){

hbt_gen_jet_V0(GoodTrackFourVector_Lam_Jet, GoodTrackFourVector_Lam_NJet, (Double_t) aux_N_tk_offline, h_Lam_Lam_samejet, h_inv_Lam_Lam_samejet, h_rot_Lam_Lam_samejet, true, massLAL, sigmaLAL, nsigmapeak);
hbt_gen_jet_V0(GoodTrackFourVector_Lam_Jet, GoodTrackFourVector_Lam_NJet, (Double_t) aux_N_tk_offline, h_Lam_Lam_diffjet, h_inv_Lam_Lam_diffjet, h_rot_Lam_Lam_diffjet, false, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxLam_samediffjet.push_back(aux_vtxz);
ev_ntrkoff_vecLam_samediffjet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLam_samediffjet.push_back(GoodTrackFourVector_Lam_Jet); 
ev_njet_vecLam_samediffjet.push_back(GoodTrackFourVector_Lam_NJet);

}

if(GoodTrackFourVector_ALam_Jet.size()>=2){

hbt_gen_jet_V0(GoodTrackFourVector_ALam_Jet, GoodTrackFourVector_ALam_NJet, (Double_t) aux_N_tk_offline, h_ALam_ALam_samejet, h_inv_ALam_ALam_samejet, h_rot_ALam_ALam_samejet, true, massLAL, sigmaLAL, nsigmapeak);
hbt_gen_jet_V0(GoodTrackFourVector_ALam_Jet, GoodTrackFourVector_ALam_NJet, (Double_t) aux_N_tk_offline, h_ALam_ALam_diffjet, h_inv_ALam_ALam_diffjet, h_rot_ALam_ALam_diffjet, false, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxALam_samediffjet.push_back(aux_vtxz);
ev_ntrkoff_vecALam_samediffjet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecALam_samediffjet.push_back(GoodTrackFourVector_ALam_Jet); 
ev_njet_vecALam_samediffjet.push_back(GoodTrackFourVector_ALam_NJet);

}

if(GoodTrackFourVector_Lam_Jet.size()>=1 && GoodTrackFourVector_ALam_Jet.size()>=1){

hbt_gen_jet_V0_OS(GoodTrackFourVector_Lam_Jet, GoodTrackFourVector_Lam_NJet, GoodTrackFourVector_ALam_Jet, GoodTrackFourVector_ALam_NJet, (Double_t) aux_N_tk_offline, h_Lam_ALam_samejet, h_inv_Lam_ALam_samejet, h_rot_Lam_ALam_samejet, true, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak);
hbt_gen_jet_V0_OS(GoodTrackFourVector_Lam_Jet, GoodTrackFourVector_Lam_NJet, GoodTrackFourVector_ALam_Jet, GoodTrackFourVector_ALam_NJet, (Double_t) aux_N_tk_offline, h_Lam_ALam_diffjet, h_inv_Lam_ALam_diffjet, h_rot_Lam_ALam_diffjet, false, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxLamALam_samediffjet.push_back(aux_vtxz);
ev_ntrkoff_vecLamALam_samediffjet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecLamALam_samediffjet.push_back(GoodTrackFourVector_Lam_Jet); 
ev_njet_vecLamALam_samediffjet.push_back(GoodTrackFourVector_Lam_NJet); 
ev_GoodTrackFourVector_vecALamLam_samediffjet.push_back(GoodTrackFourVector_ALam_Jet); 
ev_njet_vecALamLam_samediffjet.push_back(GoodTrackFourVector_ALam_NJet); 

}

if(GoodTrackFourVector_K0s_Jet.size()>=1 && GoodTrackFourVector_Lam_Jet.size()>=1){

hbt_gen_jet_V0_OS(GoodTrackFourVector_K0s_Jet, GoodTrackFourVector_K0s_NJet, GoodTrackFourVector_Lam_Jet, GoodTrackFourVector_Lam_NJet, (Double_t) aux_N_tk_offline, h_K0s_Lam_samejet, h_inv_K0s_Lam_samejet, h_rot_K0s_Lam_samejet, true, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak);
hbt_gen_jet_V0_OS(GoodTrackFourVector_K0s_Jet, GoodTrackFourVector_K0s_NJet, GoodTrackFourVector_Lam_Jet, GoodTrackFourVector_Lam_NJet, (Double_t) aux_N_tk_offline, h_K0s_Lam_diffjet, h_inv_K0s_Lam_diffjet, h_rot_K0s_Lam_diffjet, false, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxK0sLam_samediffjet.push_back(aux_vtxz);
ev_ntrkoff_vecK0sLam_samediffjet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0sLam_samediffjet.push_back(GoodTrackFourVector_K0s_Jet); 
ev_njet_vecK0sLam_samediffjet.push_back(GoodTrackFourVector_K0s_NJet); 
ev_GoodTrackFourVector_vecLamK0s_samediffjet.push_back(GoodTrackFourVector_Lam_Jet); 
ev_njet_vecLamK0s_samediffjet.push_back(GoodTrackFourVector_Lam_NJet); 

}

if(GoodTrackFourVector_K0s_Jet.size()>=1 && GoodTrackFourVector_ALam_Jet.size()>=1){

hbt_gen_jet_V0_OS(GoodTrackFourVector_K0s_Jet, GoodTrackFourVector_K0s_NJet, GoodTrackFourVector_ALam_Jet, GoodTrackFourVector_ALam_NJet, (Double_t) aux_N_tk_offline, h_K0s_ALam_samejet, h_inv_K0s_ALam_samejet, h_rot_K0s_ALam_samejet, true, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak);
hbt_gen_jet_V0_OS(GoodTrackFourVector_K0s_Jet, GoodTrackFourVector_K0s_NJet, GoodTrackFourVector_ALam_Jet, GoodTrackFourVector_ALam_NJet, (Double_t) aux_N_tk_offline, h_K0s_ALam_diffjet, h_inv_K0s_ALam_diffjet, h_rot_K0s_ALam_diffjet, false, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak);

ev_z_vtxK0sALam_samediffjet.push_back(aux_vtxz);
ev_ntrkoff_vecK0sALam_samediffjet.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecK0sALam_samediffjet.push_back(GoodTrackFourVector_K0s_Jet); 
ev_njet_vecK0sALam_samediffjet.push_back(GoodTrackFourVector_K0s_NJet); 
ev_GoodTrackFourVector_vecALamK0s_samediffjet.push_back(GoodTrackFourVector_ALam_Jet); 
ev_njet_vecALamK0s_samediffjet.push_back(GoodTrackFourVector_ALam_NJet); 

}


}

GoodTrackFourVector_Jet.clear();

GoodTrackFourVector_K0s_Jet.clear();
GoodTrackFourVector_K0s_NJet.clear();
GoodTrackFourVector_K0s_NotaJet.clear();

GoodTrackFourVector_Lam_Jet.clear();
GoodTrackFourVector_Lam_NJet.clear();
GoodTrackFourVector_Lam_NotaJet.clear();

GoodTrackFourVector_ALam_Jet.clear();
GoodTrackFourVector_ALam_NJet.clear();
GoodTrackFourVector_ALam_NotaJet.clear();



if(isMC){
    

const int a1x = GoodTrackFourVectorT_K0s.size();
Int_t NmatchK0s_GenT[a1x];
const int a2x = GoodTrackFourVectorT_Lam.size();
Int_t NmatchLam_GenT[a2x];
const int a3x = GoodTrackFourVectorT_ALam.size();
Int_t NmatchALam_GenT[a3x];

//matching


if(map_GenK0s_matched.size() >= 1 && GoodTrackFourVectorT_K0s.size() >= 1){

for(Int_t ik1x=0; ik1x<GoodTrackFourVectorT_K0s.size(); ++ik1x){
	double      pTR = GoodTrackFourVectorT_K0s[ik1x].Pt();
	double      etaR = GoodTrackFourVectorT_K0s[ik1x].Eta();
	double      phiR = GoodTrackFourVectorT_K0s[ik1x].Phi();
    NmatchK0s_GenT[ik1x] = -9999;
for(Int_t ik2x=0; ik2x<map_GenK0s_matched.size(); ++ik2x){
	double pTG = map_GenK0s_matched[ik2x].Pt();
	double etaG = map_GenK0s_matched[ik2x].Eta();
	double phiG = map_GenK0s_matched[ik2x].Phi();
	double dphi = GoodTrackFourVectorT_K0s[ik1x].DeltaPhi(map_GenK0s_matched[ik2x]);//phiR - phiG;
    double deta = etaR - etaG;
    double dR = sqrt(dphi*dphi + deta*deta);
    double dpt = (pTG - pTR)/(pTG);
	DR_K0s->Fill(dR);
	DPT_K0s->Fill(dpt);
    if(dR < drdpt && fabs(dpt) < drdpt){
        NmatchK0s_GenT[ik1x] = ik2x;
        break;
}else{}}}

for(Int_t ik1x=0; ik1x<GoodTrackFourVectorT_K0s.size(); ++ik1x){
	double      pTR = GoodTrackFourVectorT_K0s[ik1x].Pt();
	double      etaR = GoodTrackFourVectorT_K0s[ik1x].Eta();
	double      phiR = GoodTrackFourVectorT_K0s[ik1x].Phi();
	bool samegen = false;
	for(Int_t ik1xx=ik1x+1; ik1xx<GoodTrackFourVectorT_K0s.size(); ++ik1xx){if(NmatchK0s_GenT[ik1x] >= 0 && NmatchK0s_GenT[ik1xx] >= 0 && NmatchK0s_GenT[ik1x] == NmatchK0s_GenT[ik1xx]){samegen=true;}}
	if(samegen)continue;
for(Int_t ik2x=0; ik2x<map_GenK0s_matched.size(); ++ik2x){
	double pTG = map_GenK0s_matched[ik2x].Pt();
	double etaG = map_GenK0s_matched[ik2x].Eta();
	double phiG = map_GenK0s_matched[ik2x].Phi();
	double dphi = GoodTrackFourVectorT_K0s[ik1x].DeltaPhi(map_GenK0s_matched[ik2x]);//phiR - phiG;
    double deta = etaR - etaG;
    double dR = sqrt(dphi*dphi + deta*deta);
    double dpt = (pTG - pTR)/(pTG);
    if(dR < drdpt && fabs(dpt) < drdpt){
    if(GoodTrackFourVectorT_K0s[ik1x].M() > massK0s - nsigmapeak*sigmaK0s && GoodTrackFourVectorT_K0s[ik1x].M() < massK0s + nsigmapeak*sigmaK0s){map_GenK0s_matchedT.push_back(map_GenK0s_matched[ik2x]);}
    break;
}else{}}}

}


if(map_GenLam_matched.size() >= 1 && GoodTrackFourVectorT_Lam.size() >= 1){

for(Int_t ik1x=0; ik1x<GoodTrackFourVectorT_Lam.size(); ++ik1x){
	double      pTR = GoodTrackFourVectorT_Lam[ik1x].Pt();
	double      etaR = GoodTrackFourVectorT_Lam[ik1x].Eta();
	double      phiR = GoodTrackFourVectorT_Lam[ik1x].Phi();
    NmatchLam_GenT[ik1x] = -9999;
for(Int_t ik2x=0; ik2x<map_GenLam_matched.size(); ++ik2x){
	double pTG = map_GenLam_matched[ik2x].Pt();
	double etaG = map_GenLam_matched[ik2x].Eta();
	double phiG = map_GenLam_matched[ik2x].Phi();
	double dphi = GoodTrackFourVectorT_Lam[ik1x].DeltaPhi(map_GenLam_matched[ik2x]);//phiR - phiG;
    double deta = etaR - etaG;
    double dR = sqrt(dphi*dphi + deta*deta);
    double dpt = (pTG - pTR)/(pTG);
	DR_Lam->Fill(dR);
	DPT_Lam->Fill(dpt);
    if(dR < drdpt && fabs(dpt) < drdpt){
        NmatchLam_GenT[ik1x] = ik2x;
        break;
}else{}}}

for(Int_t ik1x=0; ik1x<GoodTrackFourVectorT_Lam.size(); ++ik1x){
	double      pTR = GoodTrackFourVectorT_Lam[ik1x].Pt();
	double      etaR = GoodTrackFourVectorT_Lam[ik1x].Eta();
	double      phiR = GoodTrackFourVectorT_Lam[ik1x].Phi();
	bool samegen = false;
	for(Int_t ik1xx=ik1x+1; ik1xx<GoodTrackFourVectorT_Lam.size(); ++ik1xx){if(NmatchLam_GenT[ik1x] >= 0 && NmatchLam_GenT[ik1xx] >= 0 && NmatchLam_GenT[ik1x] == NmatchLam_GenT[ik1xx]){samegen=true;}}
	if(samegen)continue;
for(Int_t ik2x=0; ik2x<map_GenLam_matched.size(); ++ik2x){
	double pTG = map_GenLam_matched[ik2x].Pt();
	double etaG = map_GenLam_matched[ik2x].Eta();
	double phiG = map_GenLam_matched[ik2x].Phi();
	double dphi = GoodTrackFourVectorT_Lam[ik1x].DeltaPhi(map_GenLam_matched[ik2x]);//phiR - phiG;
    double deta = etaR - etaG;
    double dR = sqrt(dphi*dphi + deta*deta);
    double dpt = (pTG - pTR)/(pTG);
    if(dR < drdpt && fabs(dpt) < drdpt){
    if(GoodTrackFourVectorT_Lam[ik1x].M() > massLAL - nsigmapeak*sigmaLAL && GoodTrackFourVectorT_Lam[ik1x].M() < massLAL + nsigmapeak*sigmaLAL){map_GenLam_matchedT.push_back(map_GenLam_matched[ik2x]);}
    break;
}else{}}}

}


if(map_GenALam_matched.size() >= 1 && GoodTrackFourVectorT_ALam.size() >= 1){

for(Int_t ik1x=0; ik1x<GoodTrackFourVectorT_ALam.size(); ++ik1x){
	double      pTR = GoodTrackFourVectorT_ALam[ik1x].Pt();
	double      etaR = GoodTrackFourVectorT_ALam[ik1x].Eta();
	double      phiR = GoodTrackFourVectorT_ALam[ik1x].Phi();
    NmatchALam_GenT[ik1x] = -9999;
for(Int_t ik2x=0; ik2x<map_GenALam_matched.size(); ++ik2x){
	double pTG = map_GenALam_matched[ik2x].Pt();
	double etaG = map_GenALam_matched[ik2x].Eta();
	double phiG = map_GenALam_matched[ik2x].Phi();
	double dphi = GoodTrackFourVectorT_ALam[ik1x].DeltaPhi(map_GenALam_matched[ik2x]);//phiR - phiG;
    double deta = etaR - etaG;
    double dR = sqrt(dphi*dphi + deta*deta);
    double dpt = (pTG - pTR)/(pTG);
	DR_ALam->Fill(dR);
	DPT_ALam->Fill(dpt);
    if(dR < drdpt && fabs(dpt) < drdpt){
        NmatchALam_GenT[ik1x] = ik2x;
        break;
}else{}}}

for(Int_t ik1x=0; ik1x<GoodTrackFourVectorT_ALam.size(); ++ik1x){
	double      pTR = GoodTrackFourVectorT_ALam[ik1x].Pt();
	double      etaR = GoodTrackFourVectorT_ALam[ik1x].Eta();
	double      phiR = GoodTrackFourVectorT_ALam[ik1x].Phi();
	bool samegen = false;
	for(Int_t ik1xx=ik1x+1; ik1xx<GoodTrackFourVectorT_ALam.size(); ++ik1xx){if(NmatchALam_GenT[ik1x] >= 0 && NmatchALam_GenT[ik1xx] >= 0 && NmatchALam_GenT[ik1x] == NmatchALam_GenT[ik1xx]){samegen=true;}}
	if(samegen)continue;
for(Int_t ik2x=0; ik2x<map_GenALam_matched.size(); ++ik2x){
	double pTG = map_GenALam_matched[ik2x].Pt();
	double etaG = map_GenALam_matched[ik2x].Eta();
	double phiG = map_GenALam_matched[ik2x].Phi();
	double dphi = GoodTrackFourVectorT_ALam[ik1x].DeltaPhi(map_GenALam_matched[ik2x]);//phiR - phiG;
    double deta = etaR - etaG;
    double dR = sqrt(dphi*dphi + deta*deta);
    double dpt = (pTG - pTR)/(pTG);
    if(dR < drdpt && fabs(dpt) < drdpt){
    if(GoodTrackFourVectorT_ALam[ik1x].M() > massLAL - nsigmapeak*sigmaLAL && GoodTrackFourVectorT_ALam[ik1x].M() < massLAL + nsigmapeak*sigmaLAL){map_GenALam_matchedT.push_back(map_GenALam_matched[ik2x]);}
    break;
}else{}}}

}


if(map_GenK0s_matchedT.size()>=2){

hbt(map_GenK0s_matchedT, (Double_t) aux_N_tk_offline, hGen_K0s_K0s_MatchedT);

ev_z_vtxK0sGenT.push_back(aux_vtxz);
ev_ntrkoff_vecGenK0sT.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecGenK0sT.push_back(map_GenK0s_matchedT); 

}

if(map_GenLam_matchedT.size()>=2){

hbt(map_GenLam_matchedT, (Double_t) aux_N_tk_offline, hGen_Lam_Lam_MatchedT);

ev_z_vtxLamGenT.push_back(aux_vtxz);
ev_ntrkoff_vecGenLamT.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecGenLamT.push_back(map_GenLam_matchedT); 

}

if(map_GenALam_matchedT.size()>=2){

hbt(map_GenALam_matchedT, (Double_t) aux_N_tk_offline, hGen_ALam_ALam_MatchedT);

ev_z_vtxALamGenT.push_back(aux_vtxz);
ev_ntrkoff_vecGenALamT.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecGenALamT.push_back(map_GenALam_matchedT); 

}

if(map_GenLam_matchedT.size()>=1 && map_GenALam_matchedT.size()>=1){

hbtOS(map_GenLam_matchedT, map_GenALam_matchedT, (Double_t) aux_N_tk_offline, hGen_Lam_ALam_MatchedT);

ev_z_vtxLamALamGenT.push_back(aux_vtxz);
ev_ntrkoff_vecGenLamALamT.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecGenLamALamT.push_back(map_GenLam_matchedT); 
ev_GoodTrackFourVector_vecGenALamLamT.push_back(map_GenALam_matchedT); 

}

if(map_GenK0s_matchedT.size()>=1 && map_GenLam_matchedT.size()>=1){

hbtOS(map_GenK0s_matchedT, map_GenLam_matchedT, (Double_t) aux_N_tk_offline, hGen_K0s_Lam_MatchedT);

ev_z_vtxK0sLamGenT.push_back(aux_vtxz);
ev_ntrkoff_vecGenK0sLamT.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecGenK0sLamT.push_back(map_GenK0s_matchedT); 
ev_GoodTrackFourVector_vecGenLamK0sT.push_back(map_GenLam_matchedT); 

}

if(map_GenK0s_matchedT.size()>=1 && map_GenALam_matchedT.size()>=1){

hbtOS(map_GenK0s_matchedT, map_GenALam_matchedT, (Double_t) aux_N_tk_offline, hGen_K0s_ALam_MatchedT);

ev_z_vtxK0sALamGenT.push_back(aux_vtxz);
ev_ntrkoff_vecGenK0sALamT.push_back(aux_N_tk_offline);
ev_GoodTrackFourVector_vecGenK0sALamT.push_back(map_GenK0s_matchedT); 
ev_GoodTrackFourVector_vecGenALamK0sT.push_back(map_GenALam_matchedT); 

}

}

map_GenK0s_matched.clear();
map_GenLam_matched.clear();
map_GenALam_matched.clear();
map_GenK0s_matchedT.clear();
map_GenLam_matchedT.clear();
map_GenALam_matchedT.clear();




/*========================= DONE =========================*/
}

if(i==n_ev-1)cout << "RECO is DONE!"  << endl;  

}//end the loop over events

cout << "" << endl;

cout << "------------Time for mixing------------" << endl;

cout << "" << endl;


cout << "===========================K0sK0s===========================" << endl;

call_mix(mult, Evt_mix, ev_ntrkoff_vecK0s, ev_z_vtxK0s, vtxcut,  ev_GoodTrackFourVector_vecK0s, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_T_mix, hside_K0s_K0s_T_mix, hpeakside_K0s_K0s_T_mix, hsideL_K0s_K0s_T_mix, hsideR_K0s_K0s_T_mix, hpeaksideL_K0s_K0s_T_mix, hpeaksideR_K0s_K0s_T_mix, p, effhisto_K0s, nev_K0s_ssT_mix, nev_K0s_bbT_mix, nev_K0s_sbT_mix,weight_K0sK0s);

call_mix_eta(mult, ev_GoodTrackFourVector_etaMixWeight_vecK0s, ev_ntrkoff_etaMixWeight_vecK0s, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_T_etamix, hside_K0s_K0s_T_etamix, hpeakside_K0s_K0s_T_etamix, hsideL_K0s_K0s_T_etamix, hsideR_K0s_K0s_T_etamix, hpeaksideL_K0s_K0s_T_etamix, hpeaksideR_K0s_K0s_T_etamix, p, effhisto_K0s, nev_K0s_ssT_etamix, nev_K0s_bbT_etamix, nev_K0s_sbT_etamix,weight_K0sK0s);

if(!isMC)call_mix_hdibaryon(mult, Evt_mix, ev_ntrkoff_vecK0s, ev_z_vtxK0s, vtxcut,  ev_GoodTrackFourVector_vecK0s, massK0s, sigmaK0s, nsigmapeak, nsigmaside, h_Hf2_K0s_K0s_mix, h_Hf2_K0s_K0s_mix_side, h_Hf2_K0s_K0s_mix_peakside, weight_K0sK0s);

cout << "===========================LamLam===========================" << endl;

call_mix(mult, Evt_mix, ev_ntrkoff_vecLam, ev_z_vtxLam, vtxcut,  ev_GoodTrackFourVector_vecLam, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_T_mix, hside_Lam_Lam_T_mix, hpeakside_Lam_Lam_T_mix, hsideL_Lam_Lam_T_mix, hsideR_Lam_Lam_T_mix, hpeaksideL_Lam_Lam_T_mix, hpeaksideR_Lam_Lam_T_mix, p, effhisto_Lam, nev_Lam_ssT_mix, nev_Lam_bbT_mix, nev_Lam_sbT_mix,weight_LL);

call_mix_eta(mult, ev_GoodTrackFourVector_etaMixWeight_vecLam, ev_ntrkoff_etaMixWeight_vecLam, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_T_etamix, hside_Lam_Lam_T_etamix, hpeakside_Lam_Lam_T_etamix, hsideL_Lam_Lam_T_etamix, hsideR_Lam_Lam_T_etamix, hpeaksideL_Lam_Lam_T_etamix, hpeaksideR_Lam_Lam_T_etamix, p, effhisto_Lam, nev_Lam_ssT_etamix, nev_Lam_bbT_etamix, nev_Lam_sbT_etamix,weight_LL);

if(!isMC)call_mix_hdibaryon(mult, Evt_mix, ev_ntrkoff_vecLam, ev_z_vtxLam, vtxcut,  ev_GoodTrackFourVector_vecLam, massLAL, sigmaLAL, nsigmapeak, nsigmaside, h_HDibaryon_Lam_Lam_mix, h_HDibaryon_Lam_Lam_mix_side, h_HDibaryon_Lam_Lam_mix_peakside,weight_LL);

cout << "===========================ALamALam===========================" << endl;

call_mix(mult, Evt_mix, ev_ntrkoff_vecALam, ev_z_vtxALam, vtxcut,  ev_GoodTrackFourVector_vecALam, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_T_mix, hside_ALam_ALam_T_mix, hpeakside_ALam_ALam_T_mix, hsideL_ALam_ALam_T_mix, hsideR_ALam_ALam_T_mix, hpeaksideL_ALam_ALam_T_mix, hpeaksideR_ALam_ALam_T_mix, p, effhisto_ALam, nev_ALam_ssT_mix, nev_ALam_bbT_mix, nev_ALam_sbT_mix,weight_ALAL);

call_mix_eta(mult, ev_GoodTrackFourVector_etaMixWeight_vecALam, ev_ntrkoff_etaMixWeight_vecALam, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_T_etamix, hside_ALam_ALam_T_etamix, hpeakside_ALam_ALam_T_etamix, hsideL_ALam_ALam_T_etamix, hsideR_ALam_ALam_T_etamix, hpeaksideL_ALam_ALam_T_etamix, hpeaksideR_ALam_ALam_T_etamix, p, effhisto_ALam, nev_ALam_ssT_etamix, nev_ALam_bbT_etamix, nev_ALam_sbT_etamix,weight_ALAL);

if(!isMC)call_mix_hdibaryon(mult, Evt_mix, ev_ntrkoff_vecALam, ev_z_vtxALam, vtxcut,  ev_GoodTrackFourVector_vecALam, massLAL, sigmaLAL, nsigmapeak, nsigmaside, h_HDibaryon_ALam_ALam_mix, h_HDibaryon_ALam_ALam_mix_side, h_HDibaryon_ALam_ALam_mix_peakside,weight_ALAL);

cout << "===========================LamALam===========================" << endl;

call_mix_cross(mult, Evt_mix, ev_ntrkoff_vecLamALam, ev_z_vtxLamALam, vtxcut, ev_GoodTrackFourVector_vecLamALam, ev_GoodTrackFourVector_vecALamLam, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  hpeak_Lam_ALam_T_mix, hside_Lam_ALam_T_mix, hpeakside_Lam_ALam_T_mix, hsideL_Lam_ALam_T_mix, hsideR_Lam_ALam_T_mix, hpeaksideL_Lam_ALam_T_mix, hpeaksideR_Lam_ALam_T_mix, p, effhisto_Lam, effhisto_ALam, nev_LAL_ssT_mix, nev_LAL_bbT_mix, nev_LAL_sbT_mix,weight_LAL);

call_mix_eta_cross(mult, ev_GoodTrackFourVector_etaMixWeight_vecLamALam, ev_GoodTrackFourVector_etaMixWeight_vecALamLam, ev_ntrkoff_etaMixWeight_vecLamALam, massLAL, sigmaLAL,massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_T_etamix, hside_Lam_ALam_T_etamix, hpeakside_Lam_ALam_T_etamix, hsideL_Lam_ALam_T_etamix, hsideR_Lam_ALam_T_etamix, hpeaksideL_Lam_ALam_T_etamix, hpeaksideR_Lam_ALam_T_etamix, p,  effhisto_Lam, effhisto_ALam, nev_LAL_ssT_etamix, nev_LAL_bbT_etamix, nev_LAL_sbT_etamix,weight_LAL);

cout << "===========================KLam===========================" << endl;

call_mix_cross(mult, Evt_mix, ev_ntrkoff_vecK0sLam, ev_z_vtxK0sLam, vtxcut, ev_GoodTrackFourVector_vecK0sLam, ev_GoodTrackFourVector_vecLamK0s,  massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  hpeak_K0s_Lam_T_mix, hside_K0s_Lam_T_mix, hpeakside_K0s_Lam_T_mix, hsideL_K0s_Lam_T_mix, hsideR_K0s_Lam_T_mix, hpeaksideL_K0s_Lam_T_mix, hpeaksideR_K0s_Lam_T_mix, p, effhisto_K0s, effhisto_Lam, nev_KL_ssT_mix, nev_KL_bbT_mix, nev_KL_sbT_mix,weight_K0sL);

call_mix_eta_cross(mult, ev_GoodTrackFourVector_etaMixWeight_vecK0sLam, ev_GoodTrackFourVector_etaMixWeight_vecLamK0s, ev_ntrkoff_etaMixWeight_vecK0sLam, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_T_etamix, hside_K0s_Lam_T_etamix, hpeakside_K0s_Lam_T_etamix, hsideL_K0s_Lam_T_etamix, hsideR_K0s_Lam_T_etamix, hpeaksideL_K0s_Lam_T_etamix, hpeaksideR_K0s_Lam_T_etamix, p, effhisto_K0s, effhisto_Lam, nev_KL_ssT_etamix, nev_KL_bbT_etamix, nev_KL_sbT_etamix,weight_K0sL);

if(!isMC)call_mix_hdibaryon_OS(mult, Evt_mix, ev_ntrkoff_vecK0sLam, ev_z_vtxK0sLam, vtxcut, ev_GoodTrackFourVector_vecK0sLam, ev_GoodTrackFourVector_vecLamK0s,  massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  h_Hcasc1820_K0s_Lam_mix, h_Hcasc1820_K0s_Lam_mix_side, h_Hcasc1820_K0s_Lam_mix_peakside,weight_K0sL);

cout << "===========================KALam===========================" << endl;

call_mix_cross(mult, Evt_mix, ev_ntrkoff_vecK0sALam, ev_z_vtxK0sALam, vtxcut, ev_GoodTrackFourVector_vecK0sALam, ev_GoodTrackFourVector_vecALamK0s,  massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  hpeak_K0s_ALam_T_mix, hside_K0s_ALam_T_mix, hpeakside_K0s_ALam_T_mix, hsideL_K0s_ALam_T_mix, hsideR_K0s_ALam_T_mix, hpeaksideL_K0s_ALam_T_mix, hpeaksideR_K0s_ALam_T_mix, p, effhisto_K0s, effhisto_ALam, nev_KAL_ssT_mix, nev_KAL_bbT_mix, nev_KAL_sbT_mix,weight_K0sAL);

call_mix_eta_cross(mult, ev_GoodTrackFourVector_etaMixWeight_vecK0sALam, ev_GoodTrackFourVector_etaMixWeight_vecALamK0s, ev_ntrkoff_etaMixWeight_vecK0sALam, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_T_etamix, hside_K0s_ALam_T_etamix, hpeakside_K0s_ALam_T_etamix, hsideL_K0s_ALam_T_etamix, hsideR_K0s_ALam_T_etamix, hpeaksideL_K0s_ALam_T_etamix, hpeaksideR_K0s_ALam_T_etamix, p, effhisto_K0s, effhisto_ALam, nev_KAL_ssT_etamix, nev_KAL_bbT_etamix, nev_KAL_sbT_etamix,weight_K0sAL);

if(!isMC)call_mix_hdibaryon_OS(mult, Evt_mix, ev_ntrkoff_vecK0sALam, ev_z_vtxK0sALam, vtxcut, ev_GoodTrackFourVector_vecK0sALam, ev_GoodTrackFourVector_vecALamK0s,  massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside,  h_Hcasc1820_K0s_ALam_mix, h_Hcasc1820_K0s_ALam_mix_side, h_Hcasc1820_K0s_ALam_mix_peakside,weight_K0sAL);


if(DoFeedDown){

cout << "------------------feed down------------------" << endl;

//LL
cout << endl;
cout << "Lam(feed)+Lam(feed)" << endl;
cout << endl;

cout << endl;
call_mix_test(mult, Evt_mix, ev_ntrkoff_vec_LamFD_LamFD, ev_z_vtx_LamFD_LamFD, vtxcut, ev_GoodTrackFourVector_vec_LamFD_LamFD, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_LamFD_LamFD_mix,weight_LFDLFD);
cout << endl;

cout << endl;
cout << "Lam+Lam(feed)" << endl;
cout << endl;

cout << endl;
call_mix_cross_test(mult, Evt_mix, ev_ntrkoff_vec_Lam_LamFD, ev_z_vtx_Lam_LamFD, vtxcut, ev_GoodTrackFourVector_vec_Lam_LamFD, ev_GoodTrackFourVector_vec_LamFD_Lam, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_LamFD_mix,weight_LLFD);
cout << endl;

//ALAL
cout << endl;
cout << "ALam(feed)+ALam(feed)" << endl;
cout << endl;

cout << endl;
call_mix_test(mult, Evt_mix, ev_ntrkoff_vec_ALamFD_ALamFD, ev_z_vtx_ALamFD_ALamFD, vtxcut, ev_GoodTrackFourVector_vec_ALamFD_ALamFD, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALamFD_ALamFD_mix,weight_ALFDALFD);
cout << endl;

cout << endl;
cout << "ALam+ALam(feed)" << endl;
cout << endl;

cout << endl;
call_mix_cross_test(mult, Evt_mix, ev_ntrkoff_vec_ALam_ALamFD, ev_z_vtx_ALam_ALamFD, vtxcut, ev_GoodTrackFourVector_vec_ALam_ALamFD, ev_GoodTrackFourVector_vec_ALamFD_ALam, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALamFD_mix,weight_ALALFD);
cout << endl;

//LAL

cout << endl;
cout << "Lam(feed)+ALam(feed)" << endl;
cout << endl;

cout << endl;
call_mix_cross_test(mult, Evt_mix, ev_ntrkoff_vec_LFDALFD, ev_z_vtx_LFDALFD, vtxcut, ev_GoodTrackFourVector_vec_LFDALFD, ev_GoodTrackFourVector_vec_ALFDLFD, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_LFDALFD_mix,weight_LFDALFD);
cout << endl;

cout << endl;
cout << "Lam+ALam(feed)" << endl;
cout << endl;

cout << endl;
call_mix_cross_test(mult, Evt_mix, ev_ntrkoff_vec_LALFD, ev_z_vtx_LALFD, vtxcut, ev_GoodTrackFourVector_vec_LALFD, ev_GoodTrackFourVector_vec_LFDAL, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_LALFD_mix,weight_LALFD);
cout << endl;

cout << endl;
cout << "Lam(feed)+ALam" << endl;
cout << endl;

cout << endl;
call_mix_cross_test(mult, Evt_mix, ev_ntrkoff_vec_ALLFD, ev_z_vtx_ALLFD, vtxcut, ev_GoodTrackFourVector_vec_ALLFD, ev_GoodTrackFourVector_vec_ALFDL, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_LALFD_mix,weight_ALLFD);
cout << endl;

cout << endl;
cout << "K0s+Lam(feed)" << endl;
cout << endl;

cout << endl;
call_mix_cross_test(mult, Evt_mix, ev_ntrkoff_vec_KLFD, ev_z_vtx_KLFD, vtxcut, ev_GoodTrackFourVector_vec_KLFD, ev_GoodTrackFourVector_vec_LFDK, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_KLFD_mix,weight_KLFD);
cout << endl;

cout << endl;
cout << "K0s+ALam(feed)" << endl;
cout << endl;

cout << endl;
call_mix_cross_test(mult, Evt_mix, ev_ntrkoff_vec_KALFD, ev_z_vtx_KALFD, vtxcut, ev_GoodTrackFourVector_vec_KALFD, ev_GoodTrackFourVector_vec_ALFDK, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_KALFD_mix,weight_KALFD);
cout << endl;

}


if(isMC){

cout << "------------------shape check------------------" << endl;

//K0sK0s
cout << endl;
cout << "K0sK0s" << endl;
cout << endl;

call_mix_testxxx(mult, Evt_mix, ev_ntrkoff_vecK0s_TM, ev_z_vtxK0s_TM, vtxcut, ev_GoodTrackFourVector_vecK0s_TM, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_TM_mix);
cout << endl;
call_mix_testxxx(mult, Evt_mix, ev_ntrkoff_vecK0s_TU, ev_z_vtxK0s_TU, vtxcut, ev_GoodTrackFourVector_vecK0s_TU, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_TU_mix);
cout << endl;
call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecK0s_TM_TU, ev_z_vtxK0s_TM_TU, vtxcut, ev_GoodTrackFourVector_vecK0s_TM_TU, ev_GoodTrackFourVector_vecK0s_TU_TM, massK0s, sigmaK0s, massK0s, sigmaK0s, nsigmapeak, nsigmaside, hpeak_K0s_K0s_TM_TU_mix);

//LamLam

cout << endl;
cout << "LamLam" << endl;
cout << endl;


call_mix_testxxx(mult, Evt_mix, ev_ntrkoff_vecLam_TM, ev_z_vtxLam_TM, vtxcut, ev_GoodTrackFourVector_vecLam_TM, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_TM_mix);
cout << endl;
call_mix_testxxx(mult, Evt_mix, ev_ntrkoff_vecLam_TU, ev_z_vtxLam_TU, vtxcut, ev_GoodTrackFourVector_vecLam_TU, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_TU_mix);
cout << endl;
call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecLam_TM_TU, ev_z_vtxLam_TM_TU, vtxcut, ev_GoodTrackFourVector_vecLam_TM_TU, ev_GoodTrackFourVector_vecLam_TU_TM, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_Lam_TM_TU_mix);

//ALamALam

cout << endl;
cout << "ALamALam" << endl;
cout << endl;


call_mix_testxxx(mult, Evt_mix, ev_ntrkoff_vecALam_TM, ev_z_vtxALam_TM, vtxcut, ev_GoodTrackFourVector_vecALam_TM, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_TM_mix);
cout << endl;
call_mix_testxxx(mult, Evt_mix, ev_ntrkoff_vecALam_TU, ev_z_vtxALam_TU, vtxcut, ev_GoodTrackFourVector_vecALam_TU, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_TU_mix);
cout << endl;
call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecALam_TM_TU, ev_z_vtxALam_TM_TU, vtxcut, ev_GoodTrackFourVector_vecALam_TM_TU, ev_GoodTrackFourVector_vecALam_TU_TM, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_ALam_ALam_TM_TU_mix);

//LAL

cout << endl;
cout << "LamALam" << endl;
cout << endl;

call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecLAL_TM, ev_z_vtxLAL_TM, vtxcut, ev_GoodTrackFourVector_vecLAL_TM, ev_GoodTrackFourVector_vecALL_TM, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TM_mix);
cout << endl;
call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecLAL_TU, ev_z_vtxLAL_TU, vtxcut, ev_GoodTrackFourVector_vecLAL_TU, ev_GoodTrackFourVector_vecALL_TU, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TU_mix);
cout << endl;
call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecLAL_TM_TU, ev_z_vtxLAL_TM_TU, vtxcut, ev_GoodTrackFourVector_vecLAL_TM_TU, ev_GoodTrackFourVector_vecALL_TM_TU, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TM_TU_mix);
cout << endl;
call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecLAL_TU_TM, ev_z_vtxLAL_TU_TM, vtxcut, ev_GoodTrackFourVector_vecLAL_TU_TM, ev_GoodTrackFourVector_vecALL_TU_TM, massLAL, sigmaLAL, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_Lam_ALam_TM_TU_mix);

//KL

cout << endl;
cout << "K0sLam" << endl;
cout << endl;


call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecKL_TM, ev_z_vtxKL_TM, vtxcut, ev_GoodTrackFourVector_vecKL_TM, ev_GoodTrackFourVector_vecLK_TM, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TM_mix);
cout << endl;
call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecKL_TU, ev_z_vtxKL_TU, vtxcut, ev_GoodTrackFourVector_vecKL_TU, ev_GoodTrackFourVector_vecLK_TU, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TU_mix);
cout << endl;
call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecKL_TM_TU, ev_z_vtxKL_TM_TU, vtxcut, ev_GoodTrackFourVector_vecKL_TM_TU, ev_GoodTrackFourVector_vecLK_TM_TU, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TM_TU_mix);
cout << endl;
call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecKL_TU_TM, ev_z_vtxKL_TU_TM, vtxcut, ev_GoodTrackFourVector_vecKL_TU_TM, ev_GoodTrackFourVector_vecLK_TU_TM, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_Lam_TM_TU_mix);

//KAL

cout << endl;
cout << "K0sALam" << endl;
cout << endl;


call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecKAL_TM, ev_z_vtxKAL_TM, vtxcut, ev_GoodTrackFourVector_vecKAL_TM, ev_GoodTrackFourVector_vecALK_TM, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TM_mix);
cout << endl;
call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecKAL_TU, ev_z_vtxKAL_TU, vtxcut, ev_GoodTrackFourVector_vecKAL_TU, ev_GoodTrackFourVector_vecALK_TU, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TU_mix);
cout << endl;
call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecKAL_TM_TU, ev_z_vtxKAL_TM_TU, vtxcut, ev_GoodTrackFourVector_vecKAL_TM_TU, ev_GoodTrackFourVector_vecALK_TM_TU, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TM_TU_mix);
cout << endl;
call_mix_cross_testxxx(mult, Evt_mix, ev_ntrkoff_vecKAL_TU_TM, ev_z_vtxKAL_TU_TM, vtxcut, ev_GoodTrackFourVector_vecKAL_TU_TM, ev_GoodTrackFourVector_vecALK_TU_TM, massK0s, sigmaK0s, massLAL, sigmaLAL, nsigmapeak, nsigmaside, hpeak_K0s_ALam_TM_TU_mix);



}

//jets

if(isMC){


cout << "===========================Jet-Jet===========================" << endl;

cout << endl;
call_mix_Jet(ev_ntrkoff_vecJet,ev_GoodTrackFourVector_vecJet,Evt_mix,hMix_jet_jet,ev_z_vtxJet, vtxcut);
cout << endl;
cout << "======================================================" << endl;


cout << "===========================V0+Jets===========================" << endl;
cout << endl;

cout << "-> K0s+Jet" << endl;
call_mix_gen_OS_basics(ev_ntrkoff_vecK0sJet,ev_GoodTrackFourVector_vecK0sJet,ev_GoodTrackFourVector_vecJetK0s,Evt_mix,hMix_K0s_Jet,ev_z_vtxK0sJet,vtxcut,massK0s,sigmaK0s,nsigmapeak);
cout << endl;

cout << "-> K0s_NotaJet+Jet" << endl;
call_mix_gen_OS_basics(ev_ntrkoff_vecK0snojetJet,ev_GoodTrackFourVector_vecK0snojetJet,ev_GoodTrackFourVector_vecJetK0snojet,Evt_mix,hMix_K0snojet_Jet,ev_z_vtxK0snojetJet,vtxcut,massK0s,sigmaK0s,nsigmapeak);
cout << endl;

cout << "-> K0s_Jet+Jet" << endl;
call_mix_gen_OS_basics(ev_ntrkoff_vecK0sfromjetJet,ev_GoodTrackFourVector_vecK0sfromjetJet,ev_GoodTrackFourVector_vecJetK0sfromjet,Evt_mix,hMix_K0sfromjet_Jet,ev_z_vtxK0sfromjetJet,vtxcut,massK0s,sigmaK0s,nsigmapeak);
cout << endl;

cout << "-> Lam+Jet" << endl;
call_mix_gen_OS_basics(ev_ntrkoff_vecLamJet,ev_GoodTrackFourVector_vecLamJet,ev_GoodTrackFourVector_vecJetLam,Evt_mix,hMix_Lam_Jet,ev_z_vtxLamJet,vtxcut,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> Lam_NotaJet+Jet" << endl;
call_mix_gen_OS_basics(ev_ntrkoff_vecLamnojetJet,ev_GoodTrackFourVector_vecLamnojetJet,ev_GoodTrackFourVector_vecJetLamnojet,Evt_mix,hMix_Lamnojet_Jet,ev_z_vtxLamnojetJet,vtxcut,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> Lam_Jet+Jet" << endl;
call_mix_gen_OS_basics(ev_ntrkoff_vecLamfromjetJet,ev_GoodTrackFourVector_vecLamfromjetJet,ev_GoodTrackFourVector_vecJetLamfromjet,Evt_mix,hMix_Lamfromjet_Jet,ev_z_vtxLamfromjetJet,vtxcut,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> ALam+Jet" << endl;
call_mix_gen_OS_basics(ev_ntrkoff_vecALamJet,ev_GoodTrackFourVector_vecALamJet,ev_GoodTrackFourVector_vecJetALam,Evt_mix,hMix_ALam_Jet,ev_z_vtxALamJet,vtxcut,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> ALam_NotaJet+Jet" << endl;
call_mix_gen_OS_basics(ev_ntrkoff_vecALamnojetJet,ev_GoodTrackFourVector_vecALamnojetJet,ev_GoodTrackFourVector_vecJetALamnojet,Evt_mix,hMix_ALamnojet_Jet,ev_z_vtxALamnojetJet,vtxcut,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> ALam_Jet+Jet" << endl;
call_mix_gen_OS_basics(ev_ntrkoff_vecALamfromjetJet,ev_GoodTrackFourVector_vecALamfromjetJet,ev_GoodTrackFourVector_vecJetALamfromjet,Evt_mix,hMix_ALamfromjet_Jet,ev_z_vtxALamfromjetJet,vtxcut,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "======================================================" << endl;


cout << "===========================Not-a-Jet-Study===========================" << endl;
cout << endl;

//V0NotaJet+V0NotaJet


cout << "-> K0sK0s" << endl;
call_mix_gen_V0jet(ev_ntrkoff_vecK0s_nojet,ev_GoodTrackFourVector_vecK0s_nojet,Evt_mix,hMix_K0s_K0s_nojet,ev_z_vtxK0s_nojet,vtxcut,massK0s,sigmaK0s,nsigmapeak);
cout << endl;

cout << "-> LamLam" << endl;
call_mix_gen_V0jet(ev_ntrkoff_vecLam_nojet,ev_GoodTrackFourVector_vecLam_nojet,Evt_mix,hMix_Lam_Lam_nojet,ev_z_vtxLam_nojet,vtxcut,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> ALamALam" << endl;
call_mix_gen_V0jet(ev_ntrkoff_vecALam_nojet,ev_GoodTrackFourVector_vecALam_nojet,Evt_mix,hMix_ALam_ALam_nojet,ev_z_vtxALam_nojet,vtxcut,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> LamALam" << endl;
call_mix_gen_V0jet_OS(ev_ntrkoff_vecLamALam_nojet,ev_GoodTrackFourVector_vecLamALam_nojet,ev_GoodTrackFourVector_vecALamLam_nojet,Evt_mix,hMix_Lam_ALam_nojet,ev_z_vtxLamALam_nojet,vtxcut,massLAL,sigmaLAL,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> KLam" << endl;
call_mix_gen_V0jet_OS(ev_ntrkoff_vecK0sLam_nojet,ev_GoodTrackFourVector_vecK0sLam_nojet,ev_GoodTrackFourVector_vecLamK0s_nojet,Evt_mix,hMix_K0s_Lam_nojet,ev_z_vtxK0sLam_nojet,vtxcut,massK0s,sigmaK0s,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> KALam" << endl;
call_mix_gen_V0jet_OS(ev_ntrkoff_vecK0sALam_nojet,ev_GoodTrackFourVector_vecK0sALam_nojet,ev_GoodTrackFourVector_vecALamK0s_nojet,Evt_mix,hMix_K0s_ALam_nojet,ev_z_vtxK0sALam_nojet,vtxcut,massK0s,sigmaK0s,massLAL,sigmaLAL,nsigmapeak);
cout << endl;


//V0NotaJet+V0Jet

cout << "===========================Not-a-Jet+Jet===========================" << endl;

cout << "-> K0sK0s" << endl;
call_mix_gen_V0jet_OS(ev_ntrkoff_vecK0s_only1jet,ev_GoodTrackFourVector_vecK0s_only1jet,ev_GoodTrackFourVector_vecK0s_only1jetX,Evt_mix,hMix_K0s_K0s_only1jet,ev_z_vtxK0s_only1jet,vtxcut,massK0s,sigmaK0s,massK0s,sigmaK0s,nsigmapeak);
cout << endl;

cout << "-> LamLam" << endl;
call_mix_gen_V0jet_OS(ev_ntrkoff_vecLam_only1jet,ev_GoodTrackFourVector_vecLam_only1jet,ev_GoodTrackFourVector_vecLam_only1jetX,Evt_mix,hMix_Lam_Lam_only1jet,ev_z_vtxLam_only1jet,vtxcut,massLAL,sigmaLAL,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> ALamALam" << endl;
call_mix_gen_V0jet_OS(ev_ntrkoff_vecALam_only1jet,ev_GoodTrackFourVector_vecALam_only1jet,ev_GoodTrackFourVector_vecALam_only1jetX,Evt_mix,hMix_ALam_ALam_only1jet,ev_z_vtxALam_only1jet,vtxcut,massLAL,sigmaLAL,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> LamALam" << endl;
call_mix_gen_V0jet_OS(ev_ntrkoff_vecLamALam_only1jet,ev_GoodTrackFourVector_vecLamALam_only1jet,ev_GoodTrackFourVector_vecALamLam_only1jet,Evt_mix,hMix_Lam_ALam_only1jet,ev_z_vtxLamALam_only1jet,vtxcut,massLAL,sigmaLAL,massLAL,sigmaLAL,nsigmapeak);
call_mix_gen_V0jet_OS(ev_ntrkoff_vecLamALam_only1jetX,ev_GoodTrackFourVector_vecLamALam_only1jetX,ev_GoodTrackFourVector_vecALamLam_only1jetX,Evt_mix,hMix_Lam_ALam_only1jet,ev_z_vtxLamALam_only1jetX,vtxcut,massLAL,sigmaLAL,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> KLam" << endl;
call_mix_gen_V0jet_OS(ev_ntrkoff_vecK0sLam_only1jet,ev_GoodTrackFourVector_vecK0sLam_only1jet,ev_GoodTrackFourVector_vecLamK0s_only1jet,Evt_mix,hMix_K0s_Lam_only1jet,ev_z_vtxK0sLam_only1jet,vtxcut,massK0s,sigmaK0s,massLAL,sigmaLAL,nsigmapeak);
call_mix_gen_V0jet_OS(ev_ntrkoff_vecK0sLam_only1jetX,ev_GoodTrackFourVector_vecK0sLam_only1jetX,ev_GoodTrackFourVector_vecLamK0s_only1jetX,Evt_mix,hMix_K0s_Lam_only1jet,ev_z_vtxK0sLam_only1jetX,vtxcut,massLAL,sigmaLAL,massK0s,sigmaK0s,nsigmapeak);
cout << endl;

cout << "-> KALam" << endl;
call_mix_gen_V0jet_OS(ev_ntrkoff_vecK0sALam_only1jet,ev_GoodTrackFourVector_vecK0sALam_only1jet,ev_GoodTrackFourVector_vecALamK0s_only1jet,Evt_mix,hMix_K0s_ALam_only1jet,ev_z_vtxK0sALam_only1jet,vtxcut,massK0s,sigmaK0s,massLAL,sigmaLAL,nsigmapeak);
call_mix_gen_V0jet_OS(ev_ntrkoff_vecK0sALam_only1jetX,ev_GoodTrackFourVector_vecK0sALam_only1jetX,ev_GoodTrackFourVector_vecALamK0s_only1jetX,Evt_mix,hMix_K0s_ALam_only1jet,ev_z_vtxK0sALam_only1jetX,vtxcut,massLAL,sigmaLAL,massK0s,sigmaK0s,nsigmapeak);
cout << endl;

cout << "===========================Single and Double Jets===========================" << endl;

cout << "-> K0sK0s" << endl;
call_mix_gen_V0_jet(ev_ntrkoff_vecK0s_samediffjet,ev_GoodTrackFourVector_vecK0s_samediffjet,ev_njet_vecK0s_samediffjet,Evt_mix,hMix_K0s_K0s_samejet,ev_z_vtxK0s_samediffjet,vtxcut,true,massK0s,sigmaK0s,nsigmapeak);
call_mix_gen_V0_jet(ev_ntrkoff_vecK0s_samediffjet,ev_GoodTrackFourVector_vecK0s_samediffjet,ev_njet_vecK0s_samediffjet,Evt_mix,hMix_K0s_K0s_diffjet,ev_z_vtxK0s_samediffjet,vtxcut,false,massK0s,sigmaK0s,nsigmapeak);
cout << endl;


cout << "-> LamLam" << endl;
call_mix_gen_V0_jet(ev_ntrkoff_vecLam_samediffjet,ev_GoodTrackFourVector_vecLam_samediffjet,ev_njet_vecLam_samediffjet,Evt_mix,hMix_Lam_Lam_samejet,ev_z_vtxLam_samediffjet,vtxcut,true,massLAL,sigmaLAL,nsigmapeak);
call_mix_gen_V0_jet(ev_ntrkoff_vecLam_samediffjet,ev_GoodTrackFourVector_vecLam_samediffjet,ev_njet_vecLam_samediffjet,Evt_mix,hMix_Lam_Lam_diffjet,ev_z_vtxLam_samediffjet,vtxcut,false,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> ALamALam" << endl;
call_mix_gen_V0_jet(ev_ntrkoff_vecALam_samediffjet,ev_GoodTrackFourVector_vecALam_samediffjet,ev_njet_vecALam_samediffjet,Evt_mix,hMix_ALam_ALam_samejet,ev_z_vtxALam_samediffjet,vtxcut,true,massLAL,sigmaLAL,nsigmapeak);
call_mix_gen_V0_jet(ev_ntrkoff_vecALam_samediffjet,ev_GoodTrackFourVector_vecALam_samediffjet,ev_njet_vecALam_samediffjet,Evt_mix,hMix_ALam_ALam_diffjet,ev_z_vtxALam_samediffjet,vtxcut,false,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> LamALam" << endl;
call_mix_gen_V0_jet_OS(ev_ntrkoff_vecLamALam_samediffjet,ev_GoodTrackFourVector_vecLamALam_samediffjet,ev_njet_vecLamALam_samediffjet,ev_GoodTrackFourVector_vecALamLam_samediffjet,ev_njet_vecALamLam_samediffjet,Evt_mix,hMix_Lam_ALam_samejet,ev_z_vtxLamALam_samediffjet,vtxcut,true,massLAL,sigmaLAL,massLAL,sigmaLAL,nsigmapeak);
call_mix_gen_V0_jet_OS(ev_ntrkoff_vecLamALam_samediffjet,ev_GoodTrackFourVector_vecLamALam_samediffjet,ev_njet_vecLamALam_samediffjet,ev_GoodTrackFourVector_vecALamLam_samediffjet,ev_njet_vecALamLam_samediffjet,Evt_mix,hMix_Lam_ALam_diffjet,ev_z_vtxLamALam_samediffjet,vtxcut,false,massLAL,sigmaLAL,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> KLam" << endl;
call_mix_gen_V0_jet_OS(ev_ntrkoff_vecK0sLam_samediffjet,ev_GoodTrackFourVector_vecK0sLam_samediffjet,ev_njet_vecK0sLam_samediffjet,ev_GoodTrackFourVector_vecLamK0s_samediffjet,ev_njet_vecLamK0s_samediffjet,Evt_mix,hMix_K0s_Lam_samejet,ev_z_vtxK0sLam_samediffjet,vtxcut,true,massK0s,sigmaK0s,massLAL,sigmaLAL,nsigmapeak);
call_mix_gen_V0_jet_OS(ev_ntrkoff_vecK0sLam_samediffjet,ev_GoodTrackFourVector_vecK0sLam_samediffjet,ev_njet_vecK0sLam_samediffjet,ev_GoodTrackFourVector_vecLamK0s_samediffjet,ev_njet_vecLamK0s_samediffjet,Evt_mix,hMix_K0s_Lam_diffjet,ev_z_vtxK0sLam_samediffjet,vtxcut,false,massK0s,sigmaK0s,massLAL,sigmaLAL,nsigmapeak);
cout << endl;

cout << "-> KALam" << endl;
call_mix_gen_V0_jet_OS(ev_ntrkoff_vecK0sALam_samediffjet,ev_GoodTrackFourVector_vecK0sALam_samediffjet,ev_njet_vecK0sALam_samediffjet,ev_GoodTrackFourVector_vecALamK0s_samediffjet,ev_njet_vecALamK0s_samediffjet,Evt_mix,hMix_K0s_ALam_samejet,ev_z_vtxK0sALam_samediffjet,vtxcut,true,massK0s,sigmaK0s,massLAL,sigmaLAL,nsigmapeak);
call_mix_gen_V0_jet_OS(ev_ntrkoff_vecK0sALam_samediffjet,ev_GoodTrackFourVector_vecK0sALam_samediffjet,ev_njet_vecK0sALam_samediffjet,ev_GoodTrackFourVector_vecALamK0s_samediffjet,ev_njet_vecALamK0s_samediffjet,Evt_mix,hMix_K0s_ALam_diffjet,ev_z_vtxK0sALam_samediffjet,vtxcut,false,massK0s,sigmaK0s,massLAL,sigmaLAL,nsigmapeak);
cout << endl;


}

if(isMC){

cout << "===========================Gen-tests===========================" << endl;

cout << "-> K0s+K0s" << endl;
call_mix_gen_SS(ev_ntrkoff_vecGenK0s,ev_GoodTrackFourVector_vecGenK0s,Evt_mix,hGen_K0s_K0s_Matched_mix,ev_z_vtxK0sGen,vtxcut);
cout << endl;

cout << "-> Lam+Lam" << endl;
call_mix_gen_SS(ev_ntrkoff_vecGenLam,ev_GoodTrackFourVector_vecGenLam,Evt_mix,hGen_Lam_Lam_Matched_mix,ev_z_vtxLamGen,vtxcut);
cout << endl;

cout << "-> ALam+ALam" << endl;
call_mix_gen_SS(ev_ntrkoff_vecGenALam,ev_GoodTrackFourVector_vecGenALam,Evt_mix,hGen_ALam_ALam_Matched_mix,ev_z_vtxALamGen,vtxcut);
cout << endl;

cout << "-> Lam+ALam" << endl;
call_mix_gen_OS(ev_ntrkoff_vecGenLamALam,ev_GoodTrackFourVector_vecGenLamALam,ev_GoodTrackFourVector_vecGenALamLam,Evt_mix,hGen_Lam_ALam_Matched_mix,ev_z_vtxLamALamGen,vtxcut);
cout << endl;

cout << "-> K0s+Lam" << endl;
call_mix_gen_OS(ev_ntrkoff_vecGenK0sLam,ev_GoodTrackFourVector_vecGenK0sLam,ev_GoodTrackFourVector_vecGenLamK0s,Evt_mix,hGen_K0s_Lam_Matched_mix,ev_z_vtxK0sLamGen,vtxcut);
cout << endl;

cout << "-> K0s+ALam" << endl;
call_mix_gen_OS(ev_ntrkoff_vecGenK0sALam,ev_GoodTrackFourVector_vecGenK0sALam,ev_GoodTrackFourVector_vecGenALamK0s,Evt_mix,hGen_K0s_ALam_Matched_mix,ev_z_vtxK0sALamGen,vtxcut);
cout << endl;


cout << "===========================Gen-tests-true===========================" << endl;

cout << "-> K0s+K0s" << endl;
call_mix_gen_SS(ev_ntrkoff_vecGenK0sT,ev_GoodTrackFourVector_vecGenK0sT,Evt_mix,hGen_K0s_K0s_MatchedT_mix,ev_z_vtxK0sGenT,vtxcut);
cout << endl;

cout << "-> Lam+Lam" << endl;
call_mix_gen_SS(ev_ntrkoff_vecGenLamT,ev_GoodTrackFourVector_vecGenLamT,Evt_mix,hGen_Lam_Lam_MatchedT_mix,ev_z_vtxLamGenT,vtxcut);
cout << endl;

cout << "-> ALam+ALam" << endl;
call_mix_gen_SS(ev_ntrkoff_vecGenALamT,ev_GoodTrackFourVector_vecGenALamT,Evt_mix,hGen_ALam_ALam_MatchedT_mix,ev_z_vtxALamGenT,vtxcut);
cout << endl;

cout << "-> Lam+ALam" << endl;
call_mix_gen_OS(ev_ntrkoff_vecGenLamALamT,ev_GoodTrackFourVector_vecGenLamALamT,ev_GoodTrackFourVector_vecGenALamLamT,Evt_mix,hGen_Lam_ALam_MatchedT_mix,ev_z_vtxLamALamGenT,vtxcut);
cout << endl;

cout << "-> K0s+Lam" << endl;
call_mix_gen_OS(ev_ntrkoff_vecGenK0sLamT,ev_GoodTrackFourVector_vecGenK0sLamT,ev_GoodTrackFourVector_vecGenLamK0sT,Evt_mix,hGen_K0s_Lam_MatchedT_mix,ev_z_vtxK0sLamGenT,vtxcut);
cout << endl;

cout << "-> K0s+ALam" << endl;
call_mix_gen_OS(ev_ntrkoff_vecGenK0sALamT,ev_GoodTrackFourVector_vecGenK0sALamT,ev_GoodTrackFourVector_vecGenALamK0sT,Evt_mix,hGen_K0s_ALam_MatchedT_mix,ev_z_vtxK0sALamGenT,vtxcut);
cout << endl;

}

MyFile->mkdir("Info");
MyFile->cd("Info"); 

write_info_histos();

MyFile->mkdir("HBT_K0s");
MyFile->cd("HBT_K0s"); 

write_K0s_histos();

MyFile->mkdir("HBT_Lam");
MyFile->cd("HBT_Lam"); 

write_Lam_histos();

MyFile->mkdir("HBT_ALam");
MyFile->cd("HBT_ALam"); 

write_ALam_histos();

MyFile->mkdir("HBT_LAL");
MyFile->cd("HBT_LAL"); 

write_LAL_histos();

MyFile->mkdir("HBT_K0sLam");
MyFile->cd("HBT_K0sLam"); 

write_KL_histos();

MyFile->mkdir("HBT_K0sALam");
MyFile->cd("HBT_K0sALam"); 

write_KAL_histos();


if(DoFeedDown){
MyFile->mkdir("HBT_FeedDown");
MyFile->cd("HBT_FeedDown"); 
write_feeddown();
}

if(!isMC){

MyFile->mkdir("HBT_HDibaryon");
MyFile->cd("HBT_HDibaryon"); 
write_hdibaryon();

MyFile->mkdir("HBT_f2");
MyFile->cd("HBT_f2"); 
write_f2();

MyFile->mkdir("HBT_cas1820");
MyFile->cd("HBT_cas1820"); 
write_cas1820();

}

if(isMC){

MyFile->mkdir("HBT_shape");
MyFile->cd("HBT_shape"); 
write_shapes();

MyFile->mkdir("HBT_jets");
MyFile->cd("HBT_jets"); 
write_jets();

}

MyFile->mkdir("HBT_cplots");
MyFile->cd("HBT_cplots"); 
write_cp();


MyFile->Close();

sec1 = clock();

cout << "Time: " << (double)(sec1 - sec0)/CLOCKS_PER_SEC << " [s]"<< endl;

cout << "-------------------------------------DONE--------------------------------------" << endl;

}//macro final


int main(int argc, char** argv){
        TString firstArgument(argv[1]);
        int number  = atoi(argv[2]);
        float number2 = atof(argv[3]);
        int number3 = atoi(argv[4]);
        int number4 = atoi(argv[5]);
        int number5 = atoi(argv[6]);
        int number6 = atoi(argv[7]);
        int number7 = atoi(argv[8]);
        int number8 = atoi(argv[9]);
        hbtV0(firstArgument,number,number2,number3,number4,number5,number6,number7,number8);
}
