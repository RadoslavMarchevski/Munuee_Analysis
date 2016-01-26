/* COmPACT user routine: user_init()                    */

#include "cmpio.h"
#include "user.h"
#include "reader.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TLorentzVector.h"

#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>

#include "Hist_dir.h"

using namespace std;

superMcEvent* mcevent;

char rootfilename[200];
char rootFileName[200];

char histoname[200];
char histoname2[200];

TFile *file1;

TH1F* HistnoBursts;
TH1F* HistnoGoodBursts;

TH1F* h_cutflow;

TH1I* h_evttypebytrks;
TH1I* h_ntrack;
TH1I* h_nclust;
TH1I* h_nelectrack;

TH1I* h_mutrack;
TH1I* h_mutrackq;
TH1I* h_e1trackq;
TH1I* h_e2trackq;

TH2I* h_neltrk_mutrk;
TH1I* h_Q1;
TH1I* h_1vtx_mcut;
TH1I* h_2vtx;
TH1I* h_CPRE;
TH1I* h_full_trig;
TH1I* h_k3pi_1vtx_mcut;
TH1I* h_k3pi_2vtx;
TH1I* h_k3pi_CPRE;
TH1I* h_k3pi_full_trig;
TH1I* h_mu1;
TH1I* h_nrun;
TH1I* h_Q2;
TH1I* h_Q1_NAKL;
//TH2I* h_neltrk_nphclus;

TH1F* h_dtrkcl;
TH1F* h_muchtrktime;
TH1F* h_muchtrkhodtime;
TH1F* h_muchcltime;
TH1F* h_muchtrkcltimediff;
TH1F* h_muchtrktimediff;

TH1F* h_el1trkhodtime;
TH1F* h_el1trktime;
TH1F* h_el1trktimediff;
TH1F* h_el2trkhodtime;

TH1F* h_el2trktime;
TH1F* h_el2trktimediff;
TH1F* h_e1trkcltimediff;
TH1F* h_e2trkcltimediff;
TH1F* h_ecltime1;
TH1F* h_ecltime2;
TH1F* h_gcltime;
TH1F* h_e1e2timediff;
// new
TH1F* h_e1e2hodtimediff;
TH1F* h_e1e2cltimediff;
TH1F* h_e1muhodtimediff;
TH1F* h_e2muhodtimediff;
TH1F* h_e1gtimediff;
TH1F* h_e2gtimediff;

TH1F* h_eldistDCH1plane;
//TH2F* h_eldistDCH1_ee_vz;
//TH2F* h_eldistDCH1_pch_vz;
TH1F* h_eldistLKrplane;
TH1F* h_gel1distLKrplane;
TH1F* h_gel2distLKrplane;
TH1F* h_muche1distDHC1plane;

TH1F* h_muche2distDHC1plane;
TH1F* h_muche1distLKrplane;
TH1F* h_muche2distLKrplane;
TH1F* h_gmuchdistLKrplane;

TH1F* h_eE        ;
TH1F* h_eP        ;
TH1F* h_ePt     ;
TH1F* h_elsE    ;
TH1F* h_elsPt    ;

TH1F* h_gE         ;
TH1F* h_gP         ;
TH1F* h_gPt     ;
TH1F* h_gtrkRadDHC1;

TH1F* h_muchlvE    ;
TH1F* h_muchlvPt    ;

TH1F* h_pi0E    ;
TH1F* h_mueeE    ;
TH1F* h_gmueeE    ;
TH1F* h_KE        ;

TH1F* h_nu_P;
TH1F* h_nu_E;
TH1F* h_nu_Pt;
TH1F* h_nu_m;
TH1F* h_missing_mass;
TH1F* h_missing_mass_kp;
TH1F* h_missing_mass_km;
//K3pi charged
TH1I* h_k3pi_evttypebytrks;
TH1I* h_k3pi_pitrack;
TH1I* h_k3pi_pi1trackq;
//TH1F* h_k3pi_pi2track;
TH1I* h_k3pi_pi2trackq;
//TH1F* h_k3pi_pi3track;
TH1I* h_k3pi_pi3trackq;
TH1F* h_eop_pi_track;
TH1F* h_eop_pi1_track;
TH1F* h_eop_pi2_track;
TH1F* h_eop_pi3_track;
TH1F* h_piodd_eovp    ;
TH1F* h_pi1trktime;
TH1F* h_pi1trkhodtime;
TH1F* h_pi1cltime;
TH1F* h_pi1trkcltimediff;
TH1F* h_pi1trktimediff;
TH1F* h_ecltimepi1;

TH1F* h_pi2trkhodtime;
TH1F* h_pi2trktime;
TH1F* h_pi2trktimediff;
TH1F* h_pi2trkcltimediff;
TH1F* h_ecltimepi2;

TH1F* h_pi3trkhodtime;
TH1F* h_pi3trktime;
TH1F* h_pi3trktimediff;
TH1F* h_pi3trkcltimediff;
TH1F* h_ecltimepi3;

TH1F* h_pi12hodtimediff;
TH1F* h_pi23hodtimediff;
TH1F* h_pi13hodtimediff;

TH1F* h_k3pi_3track_P;
TH1F* h_k3pi_3track_m;
TH1F* h_k3pi_3track_E;
TH1F* h_k3pi_3track_Pt;
TH1F* h_k3pi_pi1_P;
TH1F* h_k3pi_pi1_m;
TH1F* h_k3pi_pi1_E;
TH1F* h_k3pi_pi1_Pt;
TH1F* h_k3pi_pi2_P;
TH1F* h_k3pi_pi2_m;
TH1F* h_k3pi_pi2_E;
TH1F* h_k3pi_pi2_Pt;
TH1F* h_k3pi_pi3_P;
TH1F* h_k3pi_pi3_m;
TH1F* h_k3pi_pi3_E;
TH1F* h_k3pi_pi3_Pt;

TH1F* h_k3pi_missing_mass;
TH1F* h_k3pi_nu_P;
TH1F* h_k3pi_nu_m;
TH1F* h_k3pi_nu_E;
TH1F* h_k3pi_nu_Pt;

TH1F* h_k3pi_vx_12;
TH1F* h_k3pi_vy_12;
TH1F* h_k3pi_vz_12;
TH1F* h_cda_k3pic12;
TH1F* h_k3pi_vx_13;
TH1F* h_k3pi_vy_13;
TH1F* h_k3pi_vz_13;
TH1F* h_cda_k3pic13;
TH1F* h_k3pi_vx_23;
TH1F* h_k3pi_vy_23;
TH1F* h_k3pi_vz_23;
TH1F* h_cda_k3pic23;
TH1F* h_vz_pi12_13_diff;
TH1F* h_vz_pi12_23_diff;
TH1F* h_vz_pi13_23_diff;

//end K3pi charged
TH1F* h_3trk_invm;
TH1F* h_eop;
TH1F* h_eop_el;
TH1F* h_eop_mu;
TH2F* h_eop_vs_p;
TH2F* h_eop_el_vs_p;
TH2F* h_eop_mu_vs_p;
TH2F* h_eop_vs_Pt;
TH2F* h_eop_vs_piodd;
TH1F* h_e1p;
TH1F* h_e2p;
TH1F* h_e1m;
TH1F* h_e2m;
TH1F* h_elsP;
TH1F* h_mofpi0els;
TH1F* h_mofpi0;
TH2F* h_pi0_dalitz;
TH1F* h_pi0P;
TH1F* h_pi0Pt;
TH1F* h_pi0Pcm;
TH1F* h_mofpi0diff;
TH1F* h_muchlvP;
TH1F* h_KP;
TH1F* h_KPt;
TH1F* h_KPcm;
TH1F* h_mofK;
TH1F* h_mueeP;
TH1F* h_mueePt;
TH1F* h_mueeM;
TH1F* h_muee_piee_M;
TH1F* h_muee_piee_P;
TH1F* h_muee_piee_E;
TH1F* h_muee_piee_Pt;
TH1F* h_gmueeP;
TH1F* h_gmueePt;
TH1F* h_gmueeM;

TH1F* h_bx;
TH1F* h_by;
TH1F* h_bdxdz;
TH1F* h_bdydz;
TH2F* h_bx_vs_by;
TH2F* h_bdxdz_vs_bdydz;

TH1F* h_bx_much;
TH1F* h_by_much;
TH1F* h_bdxdz_much;
TH1F* h_bdydz_much;
TH2F* h_bx_vs_by_much;
TH2F* h_bdxdz_vs_bdydz_much;
TH2F* h_vx_vs_vy;
TH1F* h_vx;
TH1F* h_vy;
TH1F* h_vz;
TH1F* h_vz_initial;
TH1F* h_cda_ee;

TH2F* h_vx_vs_vy_e1much;
TH2F* h_vx_vs_vy_e2much;

TH1F* h_vx_e1much;
TH1F* h_vy_e1much;
TH1F* h_vz_e1much;

TH1F* h_vx_e2much;
TH1F* h_vy_e2much;
TH1F* h_vz_e2much;
TH1F* h_cda_e1much;
TH1F* h_cda_e2much;
TH1F* h_vz_e1much_diff;
TH1F* h_vz_e2much_diff;
TH1F* h_vz_mu1mu2_diff;
TH1F* h_vz_e12_sevt_diff;
TH1F* h_vz_e1mu_sevt_diff;
TH1F* h_vz_e2mu_sevt_diff;


TH1F* h_angle_e1e2;
TH1F* h_angle_e1e2_mu;
TH1F* h_angle_ph_e1e2;
TH1F* h_angle_much_pi0;
TH2F* h_a_e1e2_phe1e2    ;
TH2F* h_a_e1e2_muchpi0    ;
TH2F* h_a_phe1e2_muchpi0    ;

TH1F* h_angle_much_pi0_Kcm;

//Monte Carlo true
TH1F* h_true_M_0;
TH1F* h_true_P_0;
TH1F* h_true_E_0;
TH1F* h_true_Pt_0;
TH1F* h_true_M_1;
TH1F* h_true_P_1;
TH1F* h_true_E_1;
TH1F* h_true_Pt_1;
TH1F* h_true_M_2;
TH1F* h_true_P_2;
TH1F* h_true_E_2;
TH1F* h_true_Pt_2;
TH1F* h_true_M_3;
TH1F* h_true_P_3;
TH1F* h_true_E_3;
TH1F* h_true_Pt_3;
TH1F* h_true_M_4;
TH1F* h_true_P_4;
TH1F* h_true_E_4;
TH1F* h_true_Pt_4;
TH1F* h_true_M_5;
TH1F* h_true_P_5;
TH1F* h_true_E_5;
TH1F* h_true_Pt_5;
TH1F* h_true_M_6;
TH1F* h_true_P_6;
TH1F* h_true_E_6;
TH1F* h_true_Pt_6;
TH1F* h_true_M_7;
TH1F* h_true_P_7;
TH1F* h_true_E_7;
TH1F* h_true_Pt_7;
TH1F* h_mc_elsP;
TH1F* h_mc_elsM;
TH1F* h_mc_muvee_M;
TH1F* h_mc_muvee_P;
TH1F* h_mc_muvee_Pt;

TH1F* h_mc_nu_M;
TH1F* h_mc_nu_M2;
TH1F* h_mc_nu_P;
TH1F* h_mc_nu_Pt;

TH1I* h_mc_Npart_tr;
TH1I* h_mc_pt4;
TH1I* h_mc_pt5;
TH1I* h_mc_pt6;
TH1I* h_mc_pt7;

TH1F* h_mc_pvtx_pi;
TH1F* h_mc_dvtx_pi;

TH1F* h_mc_pvtx_pt0;
TH1F* h_mc_pvtx_pt1;
TH1F* h_mc_pvtx_pt2;
TH1F* h_mc_pvtx_pt3;
TH1F* h_mc_pvtx_pt4;
TH1F* h_mc_pvtx_pt5;
TH1F* h_mc_pvtx_pt6;
TH1F* h_mc_pvtx_pt7;

TH1F* h_mc_dvtx_pt0;
TH1F* h_mc_dvtx_pt1;
TH1F* h_mc_dvtx_pt2;
TH1F* h_mc_dvtx_pt3;
TH1F* h_mc_dvtx_pt4;
TH1F* h_mc_dvtx_pt5;
TH1F* h_mc_dvtx_pt6;
TH1F* h_mc_dvtx_pt7;

TH1F* h_mc_vtxdiff_pi_mu0;
TH1F* h_mc_vtxdiff_pi_mu1;
TH1F* h_mc_vtxdiff_pi_mu2;

//TH1F* h_mc_vtxdiff_5;
//TH1F* h_mc_vtxdiff_6;
//TH1F* h_mc_vtxdiff_7;
//TH1F* h_mc_vtxdiff_47;
//TH1F* h_mc_vtxdiff_57;
//TH1F* h_mc_vtxdiff_67;




//vector<int> pType;
//vector<float> pvtx;
//vector<float> dvtx;
//vector<float> mu_pvtx;
int pType[20];
float pvtx[20];
float dvtx[20];
float k_dvtx[3];
vector<string> cuts;
vector<TDirectory*> dirs;



float CPDpos_leftDownCorner[256][2];
float CELLpos_leftDownCorner[256][64][2];
float CPDineff[256][50];     // 0-9: total, 10-19: Ptrack 15-20GeV, 20-29: Ptrack 20-25GeV, 30-39: Ptrack 25-35GeV, 40-49: Ptrack 35-65GeV
float CELLineff[256][64][3]; // 0: E/p>0.6, 1: E/p 0.6-0.95, 2: E/p>1.1
int   CPDindex, CELLindex;
float CELLlength = 1.975;
float CPDlength = 8 * CELLlength;


int runNo;
int magnetCurrent;
int noBursts;
int burstCounter;
int IS_DATA;          // data type == 1
int IS_MC;            // data type == 2

int Npart;
int p0Type;
int p1Type;
int p2Type;
int p3Type;
int p4Type;
int p5Type;
int p6Type;
int p7Type;
int first_Nevt=0;
int second_Nevt=0;
int nEvt_counter=0;
int rep_counter=0;

float  zVertexTrue;

int f_Ntrack;
int s_Ntrack;
int f_NCl;
int s_NCl;
int counter_ofmuons;
float p0E ;
float p0Px;
float p0Py;
float p0Pz;
float pkgen ;

float p1E ;
float p1Px;
float p1Py;
float p1Pz;

float p2E ;
float p2Px;
float p2Py;
float p2Pz;

float p3E ;
float p3Px;
float p3Py;
float p3Pz;

float p4E ;
float p4Px;
float p4Py;
float p4Pz;

float p5E ;
float p5Px;
float p5Py;
float p5Pz;

float p6E ;
float p6Px;
float p6Py;
float p6Pz;

float p7E ;
float p7Px;
float p7Py;
float p7Pz;

float p0m ;
float p1m ;
float p2m ;
float p3m ;
float p4m ;
float p5m ;
float p6m ;
float p7m ;



TLorentzVector K_true;
TLorentzVector Pi0_true;
TLorentzVector El_true;
TLorentzVector Nu_true;
TLorentzVector true_0;
TLorentzVector true_1;
TLorentzVector true_2;
TLorentzVector true_3;
TLorentzVector true_4;
TLorentzVector true_5;
TLorentzVector true_6;
TLorentzVector true_7;
TLorentzVector Muee_three_track;
TLorentzVector mc_ee;

TLorentzVector *true_V[20];
//(*true_V) = new TLorentzVector[8];

// E/p corrections for each cell
FILE *EopCorrfile;
char  EopCorrfilename[200];
float EopCorr[100][256][64]; // ab v65: corrections for different (sub-) periods (up to 100)
int   periodFlag;            // ab v65: period to be defined in superBurst.c


void replicateTH1I(TH1I*h, int replicas)    ;
void replicateTH2I(TH2I*h, int replicas)    ;
void replicateTH1F(TH1F*h, int replicas)    ;
void replicateTH2F(TH2F*h, int replicas)    ;
//New Method
Hist_dir *dir1,*dir2;

int user_init() {
  /* WARNING: do not alter things before this line */
  /*---------- Add user C code here ----------*/
  dir1 =  new  Hist_dir("Test1");
  dir2 =  new  Hist_dir("Test2");

  int   i, j, k;
  int   l, m, n;


  fprt=fopen("compact.txt","w");

  if (strcmp (rootfilename, "") == 0)
    {
      sprintf (rootFileName, "run2007.root");
    }
  else
    {
      sprintf (rootFileName, "%s.root", rootfilename);
    }
  printf ("### saving histograms to file %s\n", rootFileName);

  file1 = new TFile(rootFileName, "RECREATE");
  //  file1 = new TFile("run2007.root", "RECREATE");


//###############################################################################################
  //
  //      Cuts
  //
  //###############################################################################################

  //object selection

  //cuts.push_back(string("no_cuts"));
  //cuts.push_back(string("3tracks"));
  //cuts.push_back(string("3goodtracks"));
  //cuts.push_back(string("cluster_isolation"));
  //cuts.push_back(string("Normalization channel selection"));
  cuts.push_back(string("Initial cuts"));
  //cuts.push_back(string("Geometry cuts"));
  //cuts.push_back(string("Vertex matching"));

  //cuts.push_back(string("Momentum cuts"));
  //cuts.push_back(string("K3pi selection"));
  //cuts.push_back(string("Muee Pt cut"));
  //cuts.push_back(string("Final k3pi selection"));
  //cuts.push_back(string("Particle id"));
  //cuts.push_back(string("Inv mass - 140"));
  ////Tests
  ////cuts.push_back(string("Energy of mu momentum  > 60."));
  //cuts.push_back(string("After k3pi 3trm cut"));
  //cuts.push_back(string("Side band substraction"));
  //cuts.push_back(string("Time matching"));
  //cuts.push_back(string("electron E/p between 0.9 and 0.95 (both)"));

  //cuts.push_back(string("electron_positron"));
  //cuts.push_back(string("e1e2_hod_Tmatch"));

  //cuts.push_back(string("mu_track"));
  //cuts.push_back(string("mu_hod_Tmatch"));

  //cuts.push_back(string("e1e2_cl_Tmatch"));
  //cuts.push_back(string("e1e2_DCHinnerR"));
  //cuts.push_back(string("e1e2_DCHouterR"));
  //cuts.push_back(string("e1e2_vtx_found"));
  //cuts.push_back(string("ee_mass_140_pcut"));
  //cuts.push_back(string("e1e2_vtx_good"));
  //cuts.push_back(string("e1e2_dist_DCH1"));
  //cuts.push_back(string("e1e2_dist_LKr"));
  //

  //cuts.push_back(string("mu_ee_Tmatch"));
  //cuts.push_back(string("mu_vtx_good"));
  //


  //cuts.push_back(string("mu_ee_dist_LKr"));



  // cuts.push_back(string("ph_cluster"));
  //cuts.push_back(string("ph_energy"));
  //cuts.push_back(string("ph_RatDHC1"));


  //electron - positron matching

  //cuts.push_back(string("e1e2_DCHinnerR"));
  //cuts.push_back(string("e1e2_DCHouterR"));

  //cuts.push_back(string("e1e2_dist_DCH1"));
  //cuts.push_back(string("e1e2_dist_LKr"));

  //photon - electron - positron matching
  //cuts.push_back(string("ph_ee_Tmatch"));
  //cuts.push_back(string("ph_ee_dist_LKr"));

  //electron - positron - charged pion matching

  //cuts.push_back(string("mu_ee_vtx_match"));

  //

  //photon - charged pion matching
  //cuts.push_back(string("ph_pi_dist_LKr"));

  //higher level cuts
  //cuts.push_back(string("pi0_mass_wndw"));
  //cuts.push_back(string("K_mom_wndw"));
  //cuts.push_back(string("K_tmom_wndw"));
  //cuts.push_back(string("K_mass_wndw"));

  //MC vs. data tuning cuts
  //cuts.push_back(string("pi_P"));
  //cuts.push_back(string("pi_eop"));
  //cuts.push_back(string("pi0pi_angle"));


  int cutsize = cuts.size();

  file1->cd();
  for (int i = 0; i < cutsize; ++i) {
    TDirectory* d = file1->mkdir( cuts[i].c_str(),  cuts[i].c_str() );
    dirs.push_back(d);
  }
  file1->cd();



  //###############################################################################################
  //
  //      Histogramme
  //
  //###############################################################################################


  // Histo mit Anzahl der bearbeiteten bursts
  HistnoBursts                 = new TH1F("noBursts", "number of bursts read", 1,0,1);
  HistnoGoodBursts             = new TH1F("noGoodBursts", "number of good bursts read", 1,0,1);

  ///////////////////////////////////////////////
  // Allgemeine cut Histogramme

  //general histos
  h_cutflow                = new TH1F("h_cutflow",         "cut flow histogram",                           50,1,51);

  //replicated histos
  dirs[0]->cd();

  h_ntrack                = new TH1I("h_ntrack",             "total number of tracks",                         20,0,20);
  h_ntrack              ->GetXaxis()->SetTitle(Form("Number of tracks"));


  h_nclust                = new TH1I("h_nclust",             "total number of clusters",                     20,0,20);
  h_nclust              ->GetXaxis()->SetTitle(Form("Number of clusters"));

  h_k3pi_evttypebytrks            = new TH1I("h_k3pi_evttypebytrks",    "0-(+++) ,1-(---),2-(2+),3-(2-) ", 7,-1,6);
  h_k3pi_evttypebytrks       ->GetXaxis()->SetTitle(Form("Type of 3track events k3pi selection"));

  h_evttypebytrks            = new TH1I("h_evttypebytrks",    "0=#mu+e+e- , 3=#mu-e+e- ", 7,-1,6);
  h_evttypebytrks       ->GetXaxis()->SetTitle(Form("Type of 3track events"));


  h_nelectrack             = new TH1I("h_nelectrack",         "number of electron candidates",                 10,0,10);
  h_nelectrack          ->GetXaxis()->SetTitle(Form("Number of electron tracks"));

  h_mutrack                = new TH1I("h_mutrack",          "number of mu^{#pm} candidates",                10,0,10);
  h_mutrack             ->GetXaxis() ->SetTitle(Form("Number of muons"));

  h_k3pi_pitrack                 = new TH1I("h_k3pi_pitrack",          "number of #pi^{#pm} candidates",                10,0,10);
  h_k3pi_pitrack          ->GetXaxis() ->SetTitle(Form("Number of pions"));


  h_mutrackq                = new TH1I("h_mutrackq",          "charge of mu^{#pm} candidates",                10,-5,5);
  h_mutrackq            ->GetXaxis()->SetTitle(Form("Charge of selected muons"));
  h_e1trackq                = new TH1I("h_e1trackq",          "charge of first electron candidate",                10,-5,5);
  h_e1trackq            ->GetXaxis()->SetTitle(Form("Charge of first selected electron"));
  h_e2trackq                = new TH1I("h_e2trackq",          "charge of second electron candidate",                10,-5,5);
  h_e2trackq            ->GetXaxis()->SetTitle(Form("Charge of second selected electron"));

  h_k3pi_pi1trackq                = new TH1I("h_k3pi_pi1trackq",          "charge of first pion candidate",                10,-5,5);
  h_k3pi_pi1trackq            ->GetXaxis()->SetTitle(Form("Charge of first selected pion"));
  h_k3pi_pi2trackq                = new TH1I("h_k3pi_pi2trackq",          "charge of second pion candidate",                10,-5,5);
  h_k3pi_pi2trackq            ->GetXaxis()->SetTitle(Form("Charge of second selected pion"));
  h_k3pi_pi3trackq                = new TH1I("h_k3pi_pi3trackq",          "charge of third pion candidate",                10,-5,5);
  h_k3pi_pi3trackq            ->GetXaxis()->SetTitle(Form("Charge of third selected pion"));

  h_neltrk_mutrk            = new TH2I("h_neltrk_mutrk",    "number of electron vs. mu^{#pm} candidates",   4,0,4, 4,0,4);
  h_neltrk_mutrk        ->GetXaxis()->SetTitle(Form("Number of electrons"));
  h_neltrk_mutrk        ->GetYaxis()->SetTitle(Form("Number of muons"));

  h_Q1               = new TH1I("h_Q1",          "Q2*!AKL + Q1/100 refference trigger PU1 ch5 bit 4",                9,0,9);
  h_Q2               = new TH1I("h_Q2",          "Q2 trigger PU1 ch7 bit 6",                9,0,9);
  h_mu1              = new TH1I("h_mu1",        "2-vtx or 1-vtx_mass cut trigger",                9,0,9);
  h_Q1_NAKL          = new TH1I("h_Q1_NAKL",    "MB-1TRK-P PU1 channel 5 bit 13",            9,0,9);
  h_CPRE               = new TH1I("h_CPRE",          "Q2*!AKL + Q1/100 refference trigger  bit 3 of the trigger word",                9,0,9);
  h_k3pi_CPRE               = new TH1I("h_k3pi_CPRE",          "Q2*!AKL + Q1/100 refference trigger  bit 3 of the trigger word",                9,0,9);
  h_2vtx               = new TH1I("h_2vtx",          "MB-2VTX trigger  bit 1 of the trigger word",                9,0,9);
  h_k3pi_2vtx               = new TH1I("h_k3pi_2vtx",          "MB-2VTX trigger  bit 1 of the trigger word",                9,0,9);
  h_1vtx_mcut               = new TH1I("h_1vtx_mcut",          "MB-1VTX trigger  bit 2 of the trigger word",                9,0,9);
  h_k3pi_1vtx_mcut               = new TH1I("h_k3pi_1vtx_mcut",          "MB-1VTX trigger  bit 2 of the trigger word",                9,0,9);
  h_k3pi_full_trig               = new TH1I("h_k3pi_full_trig",          "MB-1VTX or  MB-2VTX or MB-1TRK-P and CPRE",                9,0,9);
  h_full_trig               = new TH1I("h_full_trig",          "MB-1VTX or  MB-2VTX  and CPRE",                9,0,9);
  h_nrun               = new TH1I("h_nrun",          "number of run",                1000,15000,16000);


  h_dtrkcl                = new TH1F("h_dtrkcl",             "distance between matching track and cluster at cluster's z",    500, 0,250);
  h_dtrkcl              ->GetXaxis()->SetTitle(Form("Distance between matching track and cluster"));

  h_el1trktime            = new TH1F("h_el1trktime",         "first electron DCH time ",                    200,0,200);

  h_el1trkhodtime            = new TH1F("h_el1trkhodtime",     "first electron track time in hodoscope",        300,0,300);

  h_el1trktimediff        = new TH1F("h_el1trktimediff",     "first electron track time difference DCH - HoD",            400,-200,200);

  h_el2trktime            = new TH1F("h_el2trktime",         "second track DCH time ",                            200,0,200);

  h_el2trkhodtime            = new TH1F("h_el2trkhodtime",     "second electron track time in hodoscope",        300,0,300);

  h_el2trktimediff        = new TH1F("h_el2trktimediff",     "second electron track time difference DCH - HoD",        400,-200,200);


  h_ecltime1                = new TH1F("h_ecltime1",         "first electron cluster time ",                    400,-200,200);

  h_ecltime2                = new TH1F("h_ecltime2",         "second electron cluster time ",                400,-200,200);

  h_gcltime                = new TH1F("h_gcltime",         "photon cluster time ",                            2000,-1000,1000);

  h_e1trkcltimediff        = new TH1F("h_e1trkcltimediff", "electron track DCH & cluster time difference",        160,-40,40);
  h_e2trkcltimediff        = new TH1F("h_e2trkcltimediff", "electron track DCH & cluster time difference",        160,-40,40);

  h_e1e2timediff            = new TH1F("h_e1e2timediff",     "e_{1} & e_{2} track DCH time difference ",            80,-20,20);
  h_e1e2cltimediff        = new TH1F("h_e1e2cltimediff",     "e_{1} & e_{2} cluster time difference ",        80,-20,20);

  h_e1e2hodtimediff            = new TH1F("h_e1e2hodtimediff",     "e_{1} & e_{2} track hodoscope time difference ",            80,-20,20);

  h_e1muhodtimediff            = new TH1F("h_e1muhodtimediff",     "e_{1} & #mu track hodoscope time difference ",            80,-20,20);

  h_e2muhodtimediff            = new TH1F("h_e2muhodtimediff",     "e_{2} & #mu track hodoscope time difference ",            80,-20,20);

  h_pi12hodtimediff            = new TH1F("h_pi12hodtimediff",     "#pi_{1} & #pi_{2} track hodoscope time difference ",            80,-20,20);

  h_pi13hodtimediff            = new TH1F("h_pi13hodtimediff",     "#pi_{1} &  #pi_{3} track hodoscope time difference ",            80,-20,20);

  h_pi23hodtimediff            = new TH1F("h_pi23hodtimediff",     "#pi_{2} & #pi_{3} track hodoscope time difference ",            80,-20,20);


  h_pi1trktime          = new TH1F("h_pi1trktime",         "first track DCH time ",                            200,0,200);

  h_pi1trkhodtime            = new TH1F("h_pi1trkhodtime",     "first track time in hodoscope",        300,0,300);
  h_pi1trktimediff        = new TH1F("h_pi1trktimediff",     "first track time difference DCH - HoD",        400,-200,200);
  h_ecltimepi1                = new TH1F("h_ecltimepi1",         "first track cluster time ",       400,-200,200);


  h_pi2trktime          = new TH1F("h_pi2trktime",         "second track DCH time ",                            200,0,200);

  h_pi2trkhodtime            = new TH1F("h_pi2trkhodtime",     "second track time in hodoscope",        300,0,300);
  h_pi2trktimediff        = new TH1F("h_pi2trktimediff",     "second track time difference DCH - HoD",        400,-200,200);
  h_ecltimepi2                = new TH1F("h_ecltimepi2",         "second track cluster time ",       400,-200,200);

  h_pi3trktime          = new TH1F("h_pi3trktime",         "third track DCH time ",                            200,0,200);

  h_pi3trkhodtime            = new TH1F("h_pi3trkhodtime",     "third track time in hodoscope",        300,0,300);
  h_pi3trktimediff        = new TH1F("h_pi3trktimediff",     "third track time difference DCH - HoD",        400,-200,200);
  h_ecltimepi3                = new TH1F("h_ecltimepi3",         "third track cluster time ",       400,-200,200);


  h_e1gtimediff            = new TH1F("h_e1gtimediff",     "e_{1} & photon cluster time difference ",        2000,-500,500);

  h_e2gtimediff            = new TH1F("h_e2gtimediff",     "e_{2} & photon cluster time difference ",        2000,-500,500);

  h_muchtrktime            = new TH1F("h_muchtrktime",         "#mu^{#pm} track DCH time ",                        200,0,200);

  h_muchtrkhodtime            = new TH1F("h_muchtrkhodtime",     "#mu^{#pm} track time in hodoscope",            500,-250,250);

  h_muchcltime                = new TH1F("h_muchcltime",         "#mu^{#pm} cluster time ",                        2000,-1000,1000);
  h_muchtrktimediff        = new TH1F("h_muchtrktimediff",    "#mu^{#pm} track time difference",                400,-200,200);

  h_muchtrkcltimediff        = new TH1F("h_muchtrkcltimediff","Track HoD time difference",    200,-100,100);

  h_eop                 = new TH1F("h_eop",         "track E/p",                            120,0.,1.2);


  h_eop_el                 = new TH1F("h_eop_el",         "track E/p for electrons",                            120,0.,1.2);

  h_eop_mu                 = new TH1F("h_eop_mu",         "track E/p for muons",                            120,0.,1.2);

  h_eop_pi_track                 = new TH1F("h_eop_pi_track",         "track E/p for pions",                            120,0.,1.2);

  h_eop_pi1_track                 = new TH1F("h_eop_pi1_track",         "track E/p for first pion",                            120,0.,1.2);
  h_eop_pi2_track                 = new TH1F("h_eop_pi2_track",         "track E/p for second pion",                            120,0.,1.2);
  h_eop_pi3_track                 = new TH1F("h_eop_pi3_track",         "track E/p for third pion",                            120,0.,1.2);

  h_eop_vs_p             = new TH2F("h_eop_vs_p",     "track E/p vs track momentum",            120,0.,1.2,     100,0.,100.);

  h_eop_el_vs_p             = new TH2F("h_eop_el_vs_p",     "electron E/p vs track momentum",            120,0.,1.2,     100,0.,100.);


  h_eop_mu_vs_p             = new TH2F("h_eop_mu_vs_p",     "muon E/p vs track momentum",            120,0.,1.2,     100,0.,100.);

  h_eop_vs_Pt             = new TH2F("h_eop_vs_Pt",     " E/p vs track transverse momentum for the low E/p Track",      120,0.,1.2 ,  200,0,1);
  h_eop_vs_piodd             = new TH2F("h_eop_vs_piodd",     "E/p vs track momentum for odd pion track",            120,0.,1.2,     100,0.,100.);
  h_e1p                 = new TH1F("h_e1p",         "first electron track momentum",        100,0.,100.);
  h_e1p         ->GetXaxis()->SetTitle(Form("First electron momentum [GeV]"));

  h_e2p                 = new TH1F("h_e2p",         "second electron track momentum",        100,0.,100.);
  h_e2p             ->GetXaxis()->SetTitle(Form("Second electron momentum [GeV]"));

  h_k3pi_pi1_P                 = new TH1F("h_k3pi_pi1_P",         "first pion track momentum",        100,0.,100.);
  h_k3pi_pi1_P             ->GetXaxis()->SetTitle(Form("first pion momentum [GeV]"));

  h_k3pi_pi2_P                 = new TH1F("h_k3pi_pi2_P",         "second pion track momentum",        100,0.,100.);
  h_k3pi_pi2_P             ->GetXaxis()->SetTitle(Form("second pion momentum [GeV]"));

  h_k3pi_pi3_P                 = new TH1F("h_k3pi_pi3_P",         "third pion track momentum",        100,0.,100.);
  h_k3pi_pi3_P             ->GetXaxis()->SetTitle(Form("third pion momentum [GeV]"));

  h_k3pi_pi1_E                 = new TH1F("h_k3pi_pi1_E",         "first pion track energy",        100,0.,100.);
  h_k3pi_pi1_E             ->GetXaxis()->SetTitle(Form("first pion energy [GeV]"));

  h_k3pi_pi2_E                 = new TH1F("h_k3pi_pi2_E",         "second pion track energy",        100,0.,100.);
  h_k3pi_pi2_E             ->GetXaxis()->SetTitle(Form("second pion energy [GeV]"));

  h_k3pi_pi3_E                 = new TH1F("h_k3pi_pi3_E",         "third pion track energy",        100,0.,100.);
  h_k3pi_pi3_E             ->GetXaxis()->SetTitle(Form("third pion energy [GeV]"));

  h_k3pi_pi1_Pt                 = new TH1F("h_k3pi_pi1_Pt",         "first pion track transverse momentum",        200,0.,2.);
  h_k3pi_pi1_Pt             ->GetXaxis()->SetTitle(Form("first pion transverse momentum [GeV]"));

  h_k3pi_pi2_Pt                 = new TH1F("h_k3pi_pi2_Pt",         "second pion track transverse momentum",        200,0.,2.);
  h_k3pi_pi2_Pt             ->GetXaxis()->SetTitle(Form("second pion transverse momentum [GeV]"));

  h_k3pi_pi3_Pt                 = new TH1F("h_k3pi_pi3_Pt",         "third pion track transverse momentum",        200,0.,2.);
  h_k3pi_pi3_Pt             ->GetXaxis()->SetTitle(Form("third pion transverse momentum [GeV]"));

  h_k3pi_pi1_m                 = new TH1F("h_k3pi_pi1_m",         "first pion track mass",        100,0.,0.3);
  h_k3pi_pi1_m             ->GetXaxis()->SetTitle(Form("first pion mass [GeV]"));

  h_k3pi_pi2_m                 = new TH1F("h_k3pi_pi2_m",         "second pion track mass",        100,0.,0.3);
  h_k3pi_pi2_m             ->GetXaxis()->SetTitle(Form("second pion mass [GeV]"));

  h_k3pi_pi3_m                 = new TH1F("h_k3pi_pi3_m",         "third pion track mass",        100,0.,0.3);
  h_k3pi_pi3_m             ->GetXaxis()->SetTitle(Form("third pion mass [GeV]"));

  h_k3pi_3track_P                 = new TH1F("h_k3pi_3track_P",         "three track momentum",        100,0.,100.);
  h_k3pi_3track_P             ->GetXaxis()->SetTitle(Form("three track momentum [GeV]"));
  h_k3pi_3track_E                 = new TH1F("h_k3pi_3track_E",         "three track energy",        100,0.,100.);
  h_k3pi_3track_E             ->GetXaxis()->SetTitle(Form("three track energy [GeV]"));
  h_k3pi_3track_Pt                 = new TH1F("h_k3pi_3track_Pt",         "three track transverse momentum",        200,0.,2.);
  h_k3pi_3track_Pt             ->GetXaxis()->SetTitle(Form("three track transverse momentum [GeV]"));
  h_k3pi_3track_m                 = new TH1F("h_k3pi_3track_m",         "three track mass",        400,0.,1.2);
  h_k3pi_3track_m             ->GetXaxis()->SetTitle(Form("three track mass [GeV]"));



  h_e1m                 = new TH1F("h_e1m",         "first electron mass",        100,0.,0.3);
  h_e1m         ->GetXaxis()->SetTitle(Form("First electron mass [GeV]"));

  h_e2m                 = new TH1F("h_e2m",         "second electron mass",        100,0.,0.3);
  h_e2m             ->GetXaxis()->SetTitle(Form("Second electron mass [GeV]"));

  h_eE                 = new TH1F("h_eE",             "electron energy",        100,0.,100.);
  h_eE        ->GetXaxis()->SetTitle(Form("Energy of electron [GeV]"));


  h_eP                 = new TH1F("h_eP",             "electron track momentum",        100,0.,100.);
  h_eP        ->GetXaxis()->SetTitle(Form("Momentum of electron [GeV]"));


  h_ePt                 = new TH1F("h_ePt",         "electron track transverse momentum",        200,0.,2.);
  h_ePt       ->GetXaxis()->SetTitle(Form("Second electron momentum [GeV]"));

  h_gE                 = new TH1F("h_gE",             "photon energy",                100, 0.,100.);
  h_gP                 = new TH1F("h_gP",             "photon momentum",                100, 0.,100.);
  h_gPt                = new TH1F("h_gPt",         "photon transverse momentum",                1000, 0.,1.);
  h_gtrkRadDHC1         = new TH1F("h_gtrkRadDHC1",    "photon track radius on DHC1 plane",        500, 0.,500.);

  h_muchlvE            = new TH1F("h_muchlvE",      "charged muon energy",                        100,0.,100.);
  h_muchlvE       ->GetXaxis()->SetTitle(Form("Muon energy [GeV]"));

  h_muchlvP            = new TH1F("h_muchlvP",          "charged muon momentum",                        100,0.,100.);
  h_muchlvP       ->GetXaxis()->SetTitle(Form("Muon momentum [GeV]"));

  h_muchlvPt            = new TH1F("h_muchlvPt",     "charged muon transverse momentum",                        200,0.,2.);
  h_muchlvPt       ->GetXaxis()->SetTitle(Form("Muon transverse momentum [GeV]"));

  h_elsE                = new TH1F("h_elsE",         "e^{+}e^{-} system energy",                    100,0.,100.);
  h_elsE       ->GetXaxis()->SetTitle(Form("Energy of the e^{+}e^{-} system [GeV]"));

  h_elsP                = new TH1F("h_elsP",         "e^{+}e^{-} system momentum",            100,0.,100.);
  h_elsP      ->GetXaxis()->SetTitle(Form("Momentum of the e^{+}e^{-} system [GeV]"));

  h_piodd_eovp            = new TH1F("h_piodd_eovp",          "Odd pion momentum",                        120,0.,1.2);
  h_piodd_eovp       ->GetXaxis()->SetTitle(Form("Odd pion momentum [GeV]"));


  h_elsPt                = new TH1F("h_elsPt",         "e^{+}e^{-} system transverse momentum",            500,0.,1.);
  h_elsPt      ->GetXaxis()->SetTitle(Form(" e^{+}e^{-} system transverse momentum [GeV]"));

  h_pi0E                = new TH1F("h_pi0E",            "#pi^{0} energy",                              100,0.,100.);
  h_pi0P                = new TH1F("h_pi0P",            "#pi^{0} momentum",                              100,0.,100.);
  h_pi0Pt                = new TH1F("h_pi0Pt",           "#pi^{0} transvers momentum",                      1000,0.,1.);
  h_pi0Pcm            = new TH1F("h_pi0Pcm",          "#pi^{0} momentum in cm frame",                     100,0.,1.);
  h_mofpi0            = new TH1F("h_mofpi0",         "#pi^{0} mass reconstructed from e^{+}e^{-}#gamma",2000, 0.,2.);

  h_mofpi0diff        = new TH1F("h_mofpi0diff",      "difference of #pi^{0} mass, lab vs. cm",        2000,-2.,2.);

  h_mofpi0els            = new TH1F("h_mofpi0els",     "mass of e^{+}e^{-} couple",            200, 0.,2.);
  h_mofpi0els        ->GetXaxis()->SetTitle(Form("Mass of e^{+}e^{-} couple [GeV]"));

  h_pi0_dalitz          = new TH2F("h_pi0_dalitz",  "#pi^{0} decay Dalitz distribution",                200, 0., 0.2, 200, 0.,0.2);

  h_mueeE                = new TH1F("h_mueeE",           "Energy of #mu^{#pm} e^{+}e^{-} system",             100,0.,100.);
  h_mueeE         ->GetXaxis()->SetTitle(Form("Energy of #mu^{#pm} e^{+}e^{-} system [GeV]"));

  h_mueeP                = new TH1F("h_mueeP",           "Momentum of #mu^{#pm} e^{+}e^{-} system",             100,0.,100.);
  h_mueeP         ->GetXaxis()->SetTitle(Form("Momentum of #mu^{#pm} e^{+}e^{-} system [GeV]"));

  h_mueePt            = new TH1F("h_mueePt",          "Transvers momentum of #mu^{#pm}e^{+}e^{-} system",     200,0.,2.);
  h_mueePt        ->GetXaxis()->SetTitle(Form("Transverse momentum of #mu^{#pm} e^{+}e^{-} system [GeV]"));

  h_mueeM                = new TH1F("h_mueeM",              "Mass of #mu^{#pm}e^{+}e^{-} system",             200,0.,2.);
  h_mueeM          ->GetXaxis()->SetTitle(Form("Mass of #mu^{#pm}e^{+}e^{-} system [GeV]"));

  h_muee_piee_E               = new TH1F("h_muee_piee_E","three track system  energy  #pi e e hypothesis",100,0,100);
  h_muee_piee_M               = new TH1F("h_muee_piee_M","three track system mass #pi e e hypothesis",200,0,2.);
  h_muee_piee_P               = new TH1F("h_muee_piee_P","three track system  momentum #pi e e hypothesis",100,0,100);
  h_muee_piee_Pt              = new TH1F("h_muee_piee_Pt","three track system transverse momentum #pi e e hypothesis",200,0,2);

  h_gmueeE            = new TH1F("h_gmueeE",          "Energy of #mu^{#pm} e^{+}e^{-}#gamma system",           200,0.,200.);
  h_gmueeP            = new TH1F("h_gmueeP",          "Momentum of #mu^{#pm} e^{+}e^{-}#gamma system",           200,0.,200.);
  h_gmueePt            = new TH1F("h_gmueePt",         "Transvers momentum of #mu^{#pm}e^{+}e^{-}#gamma system",   1000,0.,1.);
  h_gmueeM            = new TH1F("h_gmueeM",            "Mass of  #mu^{#pm}e^{+}e^{-}#gamma system",             200,0.,2.);

  h_KE                = new TH1F("h_KE",               "Kaon energy",                             100,0.,100.);
  h_KE         ->GetXaxis()->SetTitle(Form("Kaon energy  [GeV]"));

  h_KP                = new TH1F("h_KP",               "Kaon momentum",                             100,0.,100.);
  h_KP         ->GetXaxis()->SetTitle(Form("Kaon momentum  [GeV]"));

  h_KPt                = new TH1F("h_KPt",               "Kaon transvers momentum",                     300,0.,1.5);
  h_KPt         ->GetXaxis()->SetTitle(Form("Kaon transverse momentum  [GeV]"));

  h_KPcm                = new TH1F("h_KPcm",              "Kaon momentum in cm frame",                    100, 0.,2.);

  h_mofK                = new TH1F("h_mofK",              "Kaon mass reconstructed from #pi^{#pm} e^{+} e^{-}#gamma",2000,0.,1.);

  h_nu_P                = new TH1F("h_nu_P",               "#nu reconstructed momentum",                     100,0.,100);
  h_nu_P              ->GetXaxis()->SetTitle(Form("#nu_P [GeV]"));

  h_nu_E                = new TH1F("h_nu_E",               "#nu ",                     100,-50,50.);
  h_nu_E              ->GetXaxis()->SetTitle(Form("#nu_E [GeV]"));

  h_nu_Pt                = new TH1F("h_nu_Pt",               "#nu transvers momentum",                     200,0.,2.);
  h_nu_Pt              ->GetXaxis()->SetTitle(Form("#nu_Pt [GeV]"));

  h_nu_m                = new TH1F("h_nu_m",               "#nu missing mass",                     200,-1.,1);
  h_nu_m              ->GetXaxis()->SetTitle(Form("#nu_m [GeV]"));

  h_missing_mass        = new TH1F("h_missing_mass",    "Missing mass squared of the K -> #mu #nu e^{+} e^{-} selection",100,-0.05,0.05);

  h_missing_mass         ->GetXaxis()->SetTitle(Form("M_{K - ( mu + e^{+} + e^{-} )}^{2}  [GeV]"));

  h_missing_mass_kp        = new TH1F("h_missing_mass_kp",    "Missing mass squared of the K^+ -> #mu^+ #nu e^{+} e^{-} selection",100,-0.05,0.05);

  h_missing_mass_kp         ->GetXaxis()->SetTitle(Form("M_{K^+ - ( mu^+ + e^{+} + e^{-} )}^{2}  [GeV]"));

  h_missing_mass_km        = new TH1F("h_missing_mass_km",    "Missing mass squared of the K^- -> #mu^- #nu e^{+} e^{-} selection",100,-0.05,0.05);

  h_missing_mass_km         ->GetXaxis()->SetTitle(Form("M_{K^- - ( mu^- + e^{+} + e^{-} )}^{2}  [GeV]"));


  h_k3pi_nu_P                = new TH1F("h_k3pi_nu_P",               "#nu reconstructed momentum",                     100,0.,100);
  h_k3pi_nu_P              ->GetXaxis()->SetTitle(Form("#nu_P [GeV]"));

  h_k3pi_nu_E                = new TH1F("h_k3pi_nu_E",               "#nu ",                     100,0.,100.);
  h_k3pi_nu_E              ->GetXaxis()->SetTitle(Form("#nu_E [GeV]"));

  h_k3pi_nu_Pt                = new TH1F("h_k3pi_nu_Pt",               "#nu transvers momentum",                     200,0.,2.);
  h_k3pi_nu_Pt              ->GetXaxis()->SetTitle(Form("#nu_Pt [GeV]"));

  h_k3pi_nu_m                = new TH1F("h_k3pi_nu_m",               "#nu missing mass",                     200,-1.,1);
  h_k3pi_nu_m              ->GetXaxis()->SetTitle(Form("#nu_m [GeV]"));

  h_k3pi_missing_mass        = new TH1F("h_k3pi_missing_mass",    "Missing mass squared of the K -> #mu #nu e^{+} e^{-} selection",100,-0.05,0.05);

  h_k3pi_missing_mass         ->GetXaxis()->SetTitle(Form("M_{K - ( mu + e^{+} + e^{-} )}^{2}  [GeV]"));


  h_3trk_invm        = new TH1F("h_3trk_invm",    "Invariant mass squared of the three track system",1000,-1.,1.);

  h_3trk_invm         ->GetXaxis()->SetTitle(Form("M_{mu + e^{+} + e^{-}}^{2}  [GeV]"));



  h_eldistDCH1plane    = new TH1F("h_eldistDCH1plane",        "distance between e^{+} and e^{-} in the DCH1 plane",200,0,200);
  h_eldistDCH1plane      ->GetXaxis()->SetTitle(Form("Distance between e^{+} and e^{-} in the DCH1 plane"));

  //    h_eldistDCH1_ee_vz    = new TH2F("h_eldistDCH1plane_ee_vz",     "distance btw. e^{+} and e^{-} at DCH1 vs. z of e^{+}e^{-} vtx",        250, 0,250, 500,-2000.,8000.);//500, 0,500, 1500,-5000.,10000.);
  //    h_eldistDCH1_pch_vz    = new TH2F("h_eldistDCH1plane_pch_vz",    "distance btw. e^{+} and e^{-} at DCH1 vs. z of #pi^{#pm}e_{1} vtx",    250, 0,250, 500,-2000.,8000.);//500, 0,500, 1500,-5000.,10000.);
  h_eldistLKrplane    = new TH1F("h_eldistLKrplane",         "distance between e^{+} and e^{-} in the LKr plane",250,0,250);
  h_eldistLKrplane        ->GetXaxis()->SetTitle(Form("Distance between e^{+} and e^{-} in the LKr plane"));


  h_gel1distLKrplane    = new TH1F("h_gel1distLKrplane",    "distance between e_{1} and photon in the LKr plane",500,0,500);

  h_gel2distLKrplane    = new TH1F("h_gel2distLKrplane",    "distance between e_{2} and photon in the LKr plane",500,0,500);

  h_muche1distDHC1plane= new TH1F("h_muche1distDHC1plane",    "distance between e_{1} and #mu^{+} in the DHC1 plane",500,0,500);


  h_muche2distDHC1plane= new TH1F("h_muche2distDHC1plane",    "distance between e_{2} and #mu^{+} in the DHC1 plane",500,0,500);


  h_muche1distLKrplane    = new TH1F("h_muche1distLKrplane",    "distance between e_{1} and #mu^{+} in the LKr plane",500,0,500);


  h_muche2distLKrplane    = new TH1F("h_muche2distLKrplane",    "distance between e_{2} and #mu^{+} in the LKr plane",500,0,500);


  h_gmuchdistLKrplane    = new TH1F("h_gmuchdistLKrplane",    "distance between #mu^{#pm} and photon in the LKr plane",500,0,500);

  h_bx                 = new TH1F("h_bx",                 "electron x coordinate before magnet",        300,-150.,150.);
  h_by                 = new TH1F("h_by",                 "electron y coordinate before magnet",        300,-150.,150.);
  h_bdxdz                = new TH1F("h_bdxdz",             "electron slope on x-z plane before magnet",    200,-0.5,0.5);
  h_bdydz                = new TH1F("h_bdydz",             "electron slope on y-z plane before magnet",    200,-0.5,0.5);
  h_bx_vs_by             = new TH2F("h_bx_vs_by",         "electron x vs y coordinates before magnet",            300,-150.,150., 300,-150.,150.);
  h_bdxdz_vs_bdydz    = new TH2F("h_bdxdz_vs_bdydz",     "electron slope on x-z vs y-z plane before magnet",        200,-0.1,0.1, 200,-0.1,0.1);

  h_bx_much            = new TH1F("h_bx_much",             "#mu^{#pm} x coordinate before magnet",        300,-150.,150.);
  h_by_much            = new TH1F("h_by_much",             "#mu^{#pm} y coordinate before magnet",        300,-150.,150.);
  h_bdxdz_much            = new TH1F("h_bdxdz_much",         "#mu^{#pm} slope on x-z plane before magnet",    200,-0.5,0.5);
  h_bdydz_much            = new TH1F("h_bdydz_much",         "#mu^{#pm} slope on y-z plane before magnet",    200,-0.5,0.5);
  h_bx_vs_by_much        = new TH2F("h_bx_vs_by_much",     "#mu^{#pm} x vs y coordinates before magnet",    300,-150.,150., 300,-150.,150.);
  h_bdxdz_vs_bdydz_much= new TH2F("h_bdxdz_vs_bdydz_much","#mu^{#pm} slope on x-z vs y-z plane",        200,-0.1,0.1, 200,-0.1,0.1);

  h_cda_ee             = new TH1F("h_cda_ee",             "closest distance between e^{+} and e^{-}",        200, 0.,20.);
  h_cda_ee       ->GetXaxis()->SetTitle(Form("Closest distance between e^{+} and e^{-} [cm]"));

  h_vx                 = new TH1F("h_vx",                 "e^{+}e^{-} reconstructed vertex x",        100,-25.,25.);


  h_vy                 = new TH1F("h_vy",                 "e^{+}e^{-} reconstructed vertex y",        40,-10.,10.);


  h_vz                 = new TH1F("h_vz",                 "e^{+}e^{-} reconstructed vertex z",        1500,-5000.,10000.);
  h_vz_initial                 = new TH1F("h_vz_initial",                 " z vertex before (MC)",        2100,-5000.,16000.);

  h_k3pi_vx_12                 = new TH1F("h_k3pi_vx_12",                 "#pi_1 , #pi_2 reconstructed vertex x",        100,-25.,25.);
  h_k3pi_vy_12                 = new TH1F("h_k3pi_vy_12",                 "#pi_1 , #pi_2 reconstructed vertex y",        40,-10.,10.);
  h_k3pi_vz_12                 = new TH1F("h_k3pi_vz_12",                 "#pi_1 , #pi_2 reconstructed vertex z",        1500,-5000.,10000.);

  h_k3pi_vx_13                 = new TH1F("h_k3pi_vx_13",                 "#pi_{1} , #pi_{3} reconstructed vertex x",        100,-25.,25.);
  h_k3pi_vy_13                 = new TH1F("h_k3pi_vy_13",                 "#pi_{1} , #pi_{3} reconstructed vertex y",        40,-10.,10.);
  h_k3pi_vz_13                 = new TH1F("h_k3pi_vz_13",                 "#pi_{1} , #pi_{3} reconstructed vertex z",        1500,-5000.,10000.);

  h_k3pi_vx_23                 = new TH1F("h_k3pi_vx_23",                 "#pi_{2} , #pi_{3} reconstructed vertex x",        100,-25.,25.);
  h_k3pi_vy_23                 = new TH1F("h_k3pi_vy_23",                 "#pi_{2} , #pi_{3} reconstructed vertex y",        40,-10.,10.);
  h_k3pi_vz_23                 = new TH1F("h_k3pi_vz_23",                 "#pi_{2} , #pi_{3} reconstructed vertex z",        1500,-5000.,10000.);
  h_cda_k3pic12         = new TH1F("h_cda_k3pic12",         "closest distance between #pi_{1} and pi_{2} ",        200, 0.,50.);
  h_cda_k3pic13         = new TH1F("h_cda_k3pic13",         "closest distance between #pi_{1} and pi_{3} ",        200, 0.,50.);
  h_cda_k3pic23         = new TH1F("h_cda_k3pic23",         "closest distance between #pi_{2} and pi_{3} ",        200, 0.,50.);
  h_vz_pi12_13_diff     = new TH1F("h_vz_pi12_13_diff",   "#pi_{1},#pi_{2} and #pi_{1},#pi_{3} reconstructed vertex z difference",        100,-5000.,5000.);
  h_vz_pi12_23_diff     = new TH1F("h_vz_pi12_23_diff",   "#pi_{1},#pi_{2} and #pi_{2},#pi_{3} reconstructed vertex z difference",        100,-5000.,5000.);
  h_vz_pi13_23_diff     = new TH1F("h_vz_pi13_23_diff",   "#pi_{1},#pi_{3} and #pi_{2},#pi_{3} reconstructed vertex z difference",        100,-5000.,5000.);


  h_vz_e12_sevt_diff = new TH1F("h_vz_e12_sevt_diff",   "e_{1},e_{2} and sevt->vtx[0].z  z difference",        100,-5000.,5000.);

  h_vz_e1mu_sevt_diff = new TH1F("h_vz_e1mu_sevt_diff",   "e_{1},#mu and sevt->vtx[0].z  z difference",        100,-5000.,5000.);

  h_vz_e2mu_sevt_diff = new TH1F("h_vz_e2mu_sevt_diff",   "e_{2},#mu and sevt->vtx[0].z  z difference",        100,-5000.,5000.);

  h_vx_vs_vy             = new TH2F("h_vx_vs_vy",         "e^{+}e^{-} reconstructed vertex x vs y",        100,-25.,25., 40,-10.,10.);
  h_vx_vs_vy      ->GetXaxis()->SetTitle(Form("e^{+}e^{-} vertex x"));
  h_vx_vs_vy      ->GetYaxis()->SetTitle(Form("e^{+}e^{-} vertex y"));


  h_cda_e1much         = new TH1F("h_cda_e1much",         "closest distance between #mu^{#pm} and e_{1} ",        50, 0.,25.);
  h_cda_e1much         ->GetXaxis()->SetTitle(Form("Closest distance between #mu^{#pm} and e_{1} [cm]"));

  h_vx_e1much             = new TH1F("h_vx_e1much",         "#mu^{#pm} + e_{1}  reconstructed vertex x",        200,-50.,50.);


  h_vy_e1much             = new TH1F("h_vy_e1much",         "#mu^{#pm} + e_{1}  reconstructed vertex y",        200,-50.,50.);
  h_vz_e1much             = new TH1F("h_vz_e1much",         "#mu^{#pm} + e_{1}  reconstructed vertex z",        1500,-5000.,10000.);
  h_vx_vs_vy_e1much    = new TH2F("h_vx_vs_vy_e1much",     "#mu^{#pm} + e_{1}  reconstructed vertex x vs y",        200,-50.,50., 200,-50.,50.);

  h_vz_e1much_diff     = new TH1F("h_vz_e1much_diff",   "e^{+}e^{-} and #mu^{#pm}e_{1} reconstructed vertex z difference",        100,-5000.,5000.);


  h_cda_e2much         = new TH1F("h_cda_e2much",         "closest distance between #mu^{#pm} and e_{2} ",        50, 0.,25.);
  h_cda_e2much         ->GetXaxis()->SetTitle(Form("Closest distance between #mu^{#pm} and e_{2} [cm]"));

  h_vx_e2much             = new TH1F("h_vx_e2much",         "#mu^{#pm} + e_{2}  reconstructed vertex x",        200,-50.,50.);
  h_vy_e2much             = new TH1F("h_vy_e2much",         "#mu^{#pm} + e_{2}  reconstructed vertex y",        200,-50.,50.);
  h_vz_e2much             = new TH1F("h_vz_e2much",         "#mu^{#pm} + e_{2}  reconstructed vertex z",        1500,-5000.,10000.);

  h_vz_e2much_diff     = new TH1F("h_vz_e2much_diff",   "e^{+}e^{-} and #mu^{#pm}e_{2} reconstructed vertex z difference",        100,-5000.,5000.);

  h_vz_mu1mu2_diff     = new TH1F("h_vz_mu1mu2_diff",   "#mu^{#pm}e_{1}  and mu^{#pm}e_{2} reconstructed vertex z difference",        100,-5000.,5000.);


  h_angle_e1e2        = new TH1F("h_angle_e1e2",         "angle between e^{-} and e^{+} on their momentum plane",                         (int)(200*M_PI), 0., M_PI/30.);
  h_angle_e1e2_mu        = new TH1F("h_angle_e1e2_mu",         "angle between e^{-} and e^{+} momentum plane and the muon",                         (int)(200*M_PI), 0., M_PI/30.);

  h_angle_ph_e1e2        = new TH1F("h_angle_ph_e1e2",     "angle between #gamma and e^{-}e^{+} system on their momentum plane",             (int)(200*M_PI), 0., M_PI/30.);
  h_angle_much_pi0        = new TH1F("h_angle_much_pi0",     "angle between #mu^{#pm} and #pi^{0} on their momentum plane",                     (int)(200*M_PI), 0., M_PI/30.);
  h_a_e1e2_phe1e2        = new TH2F("h_a_e1e2_phe1e2",    "angle between e^{-} and e^{+} vs. angle #gamma and e^{-}e^{+} system",            (int)(200*M_PI), 0., M_PI/30., (int)(200*M_PI), 0., M_PI/30.);
  h_a_e1e2_muchpi0        = new TH2F("h_a_e1e2_muchpi0",    "angle between e^{-} and e^{+} vs. angle #mu^{#pm} and #pi^{0}",                (int)(200*M_PI), 0., M_PI/30., (int)(200*M_PI), 0., M_PI/30.);
  h_a_phe1e2_muchpi0    = new TH2F("h_a_phe1e2_muchpi0",    "angle between #gamma and e^{-}e^{+} system vs. angle #mu^{#pm} and #pi^{0}",    (int)(200*M_PI), 0., M_PI/30., (int)(200*M_PI), 0., M_PI/30.);

  //    h_angle_pch_pi0_Kcm    = new TH1F("h_angle_pch_pi0_Kcm", "angle between #pi^{#pm} and #pi^{0} in the K cm frame",                         (int)(400*M_PI), -2*M_PI, 2*M_PI);


  h_true_M_0               = new TH1F("true_M_0","mass of the particle type kaon",3000,0,0.5);
  h_true_P_0               = new TH1F("true_P_0","momentum of the particle type kaon",100,0,100);
  h_true_E_0               = new TH1F("true_E_0","energy of the particle type kaon",120,-20,100);
  h_true_Pt_0               = new TH1F("true_Pt_0","transverse momentum of the particle type kaon",200,0,2);

  h_true_M_1               = new TH1F("true_M_1","mass of the particle type     Pi0",3000,0,0.5);
  h_true_P_1               = new TH1F("true_P_1","momentum of the particle type Pi0",100,0,100);
  h_true_E_1               = new TH1F("true_E_1","energy of the particle type   Pi0",120,-20,100);
  h_true_Pt_1               = new TH1F("true_Pt_1","transverse momentum of the particle type Pi0",200,0,2);

  h_true_M_2               = new TH1F("true_M_2","mass of the particle     type muon",3000,0,0.5);
  h_true_P_2               = new TH1F("true_P_2","momentum of the particle type muon",100,0,100);
  h_true_E_2               = new TH1F("true_E_2","energy of the particle   type muon",120,-20,100);
  h_true_Pt_2              = new TH1F("true_Pt_2","transverse momentum of the particle type muon",200,0,2);

  h_true_M_3               = new TH1F("true_M_3","mass of the particle     type neutrino",3000,0,0.5);
  h_true_P_3               = new TH1F("true_P_3","momentum of the particle type neutrino",100,0,100);
  h_true_E_3               = new TH1F("true_E_3","energy of the particle   type neutrino",120,-20,100);
  h_true_Pt_3               = new TH1F("true_Pt_3","transverse momentum of the particle type neutrino",200,0,2);

  h_true_M_4               = new TH1F("true_M_4","mass of the particle     type positron",3000,0,0.5);
  h_true_P_4               = new TH1F("true_P_4","momentum of the particle type positron",100,0,100);
  h_true_E_4               = new TH1F("true_E_4","energy of the particle   type positron",120,-20,100);
  h_true_Pt_4               = new TH1F("true_Pt_4","transverse momentum of the particle type positron",200,0,2);

  h_true_M_5               = new TH1F("true_M_5","mass of the particle     type electron",3000,0,0.5);
  h_true_P_5               = new TH1F("true_P_5","momentum of the particle type electron",100,0,100);
  h_true_E_5               = new TH1F("true_E_5","energy of the particle   type electron",120,-20,100);
  h_true_Pt_5               = new TH1F("true_Pt_5","transverse momentum of the particle type electron",200,0,2);

  h_true_M_6               = new TH1F("true_M_6","mass of the particle     type 6",3000,0,0.5);
  h_true_P_6               = new TH1F("true_P_6","momentum of the particle type 6",100,0,100);
  h_true_E_6               = new TH1F("true_E_6","energy of the particle   type 6",120,-20,100);
  h_true_Pt_6               = new TH1F("true_Pt_6","transverse momentum of the particle type 6",200,0,2);

  h_true_M_7               = new TH1F("true_M_7","mass of the particle     type 7",3000,0,0.5);
  h_true_P_7               = new TH1F("true_P_7","momentum of the particle type 7",100,0,100);
  h_true_E_7               = new TH1F("true_E_7","energy of the particle   type 7",120,-20,100);
  h_true_Pt_7               = new TH1F("true_Pt_7","transverse momentum of the particle type 7",200,0,2);

  h_mc_muvee_M               = new TH1F("h_mc_muvee_M","three track system  true mass",200,0,2.);
  h_mc_muvee_P               = new TH1F("h_mc_muvee_P","three track system  true momentum",100,0,100);
  h_mc_muvee_Pt              = new TH1F("h_mc_muvee_Pt","three track system true transverse momentum ",200,0,2);

  h_mc_elsP               = new TH1F("h_mc_elsP"," e^{+} e^{-} system momentum",100,0,100);
  h_mc_elsM               = new TH1F("h_mc_elsM"," e^{+} e^{-} system mass",200,0,2);

  h_mc_nu_M               = new TH1F("h_mc_nu_M","neutrino mass using true values for the three tracks and the kaon",200,-0.05,0.05);
  h_mc_nu_M2               = new TH1F("h_mc_nu_M2","neutrino invariant mass squared using true values for the three tracks and the kaon",100,-0.05,0.05);
  h_mc_nu_P               = new TH1F("h_mc_nu_P","neutrino momentum using true values for the three tracks and the kaon",100,0,100);
  h_mc_nu_Pt              = new TH1F("h_mc_nu_Pt","neutrino momentum using true values for the three tracks and the kaon ",200,0,2);

  h_mc_pvtx_pi = new TH1F("h_mc_pvtx_pi",   "Production vertex for true #pi(if #piee)",         120,-3000.,9000.);
  h_mc_dvtx_pi = new TH1F("h_mc_dvtx_pi",   "Decay vertex for true #pi(if #piee)",        120,-3000.,9000.);

  h_mc_pvtx_pt0 = new TH1F("h_mc_pvtx_pt0",   "Production vertex for kaon ",        300,-14000.,16000.);
  h_mc_pvtx_pt1 = new TH1F("h_mc_pvtx_pt1",   "Production vertex for first particle (#pi if #piee),(#mu if #mu#nuee)",        300,-14000.,16000.);
  h_mc_pvtx_pt2 = new TH1F("h_mc_pvtx_pt2",   "Production vertex for e+(both #piee and #mu#nu e e)",        300,-14000.,16000.);
  h_mc_pvtx_pt3 = new TH1F("h_mc_pvtx_pt3",   "Production vertex for  e-(both #piee and #mu#nu e e)",        300,-14000.,16000.);
  h_mc_pvtx_pt4 = new TH1F("h_mc_pvtx_pt4",   "Production vertex ",        300,-14000.,16000.);
  h_mc_pvtx_pt5 = new TH1F("h_mc_pvtx_pt5",   "Production vertex ",        300,-14000.,16000.);
  h_mc_pvtx_pt6 = new TH1F("h_mc_pvtx_pt6",   "Production vertex ",        300,-14000.,16000.);
  h_mc_pvtx_pt7 = new TH1F("h_mc_pvtx_pt7",   "Production vertex ",        300,-14000.,16000.);

  h_mc_dvtx_pt0 = new TH1F("h_mc_dvtx_pt0",   "Decay vertex for kaon",        300,-14000.,16000.);
  h_mc_dvtx_pt1 = new TH1F("h_mc_dvtx_pt1",   "Decay vertex for first particle (#pi if #piee),(#mu if #mu#nuee)",        300,-14000.,16000.);
  h_mc_dvtx_pt2 = new TH1F("h_mc_dvtx_pt2",   "Decay vertex for for e+(both #piee and #mu#nu e e)",        300,-14000.,16000.);
  h_mc_dvtx_pt3 = new TH1F("h_mc_dvtx_pt3",   "Decay vertex for for e-(both #piee and #mu#nu e e)",        300,-14000.,16000.);
  h_mc_dvtx_pt4 = new TH1F("h_mc_dvtx_pt4",   "Decay vertex ",        300,-14000.,16000.);
  h_mc_dvtx_pt5 = new TH1F("h_mc_dvtx_pt5",   "Decay vertex ",        300,-14000.,16000.);
  h_mc_dvtx_pt6 = new TH1F("h_mc_dvtx_pt6",   "Decay vertex ",        300,-14000.,16000.);
  h_mc_dvtx_pt7 = new TH1F("h_mc_dvtx_pt7",   "Decay vertex ",        300,-14000.,16000.);

  h_mc_vtxdiff_pi_mu0 = new TH1F("h_mc_vtxdiff_pi_mu0",   "Vertex diff between #pi and #mu(if #piee)",         120,-3000.,9000.);
  h_mc_vtxdiff_pi_mu1 = new TH1F("h_mc_vtxdiff_pi_mu1",   "Vertex diff between #pi and #mu(if #piee)",         120,-3000.,9000.);
  h_mc_vtxdiff_pi_mu2 = new TH1F("h_mc_vtxdiff_pi_mu2",   "Vertex diff between #pi and #mu(if #piee)",         120,-3000.,9000.);
  h_mc_Npart_tr = new TH1I("h_mc_Npart_tr",          "Number of true particles",                10,0,10);
  h_mc_pt4 = new TH1I("h_mc_pt4",          "Index of particle 4 ",               550,0,550);
  h_mc_pt5 = new TH1I("h_mc_pt5",          "Index of particle 5 ",               550,0,550);
  h_mc_pt6 = new TH1I("h_mc_pt6",          "Index of particle 6 ",               550,0,550);
  h_mc_pt7 = new TH1I("h_mc_pt7",          "Index of particle 7 ",               550,0,550);


  h_cutflow->Sumw2();









  h_ntrack                    ->Sumw2();
  h_nclust                    ->Sumw2();
  h_evttypebytrks                ->Sumw2();
  h_k3pi_evttypebytrks                ->Sumw2();
  h_nelectrack                 ->Sumw2();
  h_mutrack                    ->Sumw2();
  h_k3pi_pitrack                    ->Sumw2();
  h_mutrackq                    ->Sumw2();
  h_e1trackq                    ->Sumw2();
  h_e2trackq                    ->Sumw2();

  h_k3pi_pi1trackq                    ->Sumw2();
  h_k3pi_pi2trackq                    ->Sumw2();
  h_k3pi_pi3trackq                    ->Sumw2();

  h_neltrk_mutrk                ->Sumw2();
  h_mu1                        ->Sumw2();
  h_Q1                        ->Sumw2();
  h_Q2                        ->Sumw2();
  h_Q1_NAKL                        ->Sumw2();
  h_nrun                      ->Sumw2();


  h_1vtx_mcut                         ->Sumw2();
  h_2vtx                             ->Sumw2();
  h_CPRE                             ->Sumw2();
  h_k3pi_1vtx_mcut                        ->Sumw2();
  h_k3pi_2vtx                        ->Sumw2();
  h_k3pi_CPRE                             ->Sumw2();
  h_k3pi_full_trig                             ->Sumw2();
  h_full_trig                             ->Sumw2();


  h_dtrkcl                    ->Sumw2();
  h_el1trkhodtime                ->Sumw2();
  h_el1trktime                ->Sumw2();
  h_el1trktimediff               ->Sumw2();
  h_el2trkhodtime                ->Sumw2();

  h_el2trktime                ->Sumw2();
  h_el2trktimediff               ->Sumw2();
  h_e1trkcltimediff              ->Sumw2();
  h_e2trkcltimediff              ->Sumw2();

  h_ecltime1                    ->Sumw2();
  h_ecltime2                    ->Sumw2();
  h_gcltime                    ->Sumw2();

  h_e1e2timediff                ->Sumw2();
  h_e1e2cltimediff            ->Sumw2();

  h_e1e2hodtimediff            ->Sumw2();
  h_e1muhodtimediff            ->Sumw2();
  h_e2muhodtimediff            ->Sumw2();

  h_pi12hodtimediff            ->Sumw2();
  h_pi23hodtimediff            ->Sumw2();
  h_pi13hodtimediff            ->Sumw2();

  h_pi1trktime                    ->Sumw2();
  h_pi1trkhodtime                    ->Sumw2();
  h_pi1trktimediff                    ->Sumw2();
  h_ecltimepi1                    ->Sumw2();
  h_pi2trktime                    ->Sumw2();
  h_pi2trkhodtime                    ->Sumw2();
  h_pi2trktimediff                    ->Sumw2();
  h_ecltimepi2                    ->Sumw2();
  h_pi3trktime                    ->Sumw2();
  h_pi3trkhodtime                    ->Sumw2();
  h_pi3trktimediff                    ->Sumw2();
  h_ecltimepi3                    ->Sumw2();


  h_e1gtimediff                ->Sumw2();
  h_e2gtimediff                ->Sumw2();

  h_muchtrktime                   ->Sumw2();
  h_muchcltime                    ->Sumw2();
  h_muchtrkhodtime                ->Sumw2();
  h_muchtrkcltimediff             ->Sumw2();
  h_muchtrktimediff            ->Sumw2();

  h_eop                      ->Sumw2();
  h_eop_el             ->Sumw2() ;
  h_eop_mu             ->Sumw2() ;
  h_eop_pi_track             ->Sumw2() ;
  h_eop_pi1_track             ->Sumw2() ;
  h_eop_pi2_track             ->Sumw2() ;
  h_eop_pi3_track             ->Sumw2() ;
  h_eop_vs_p                  ->Sumw2();
  h_eop_mu_vs_p        ->Sumw2() ;
  h_eop_el_vs_p        ->Sumw2() ;
  h_eop_vs_Pt          ->Sumw2() ;
  h_eop_vs_piodd          ->Sumw2() ;


  h_e1p                      ->Sumw2();
  h_e2p                      ->Sumw2();

  h_k3pi_pi1_P                ->Sumw2();
  h_k3pi_pi2_P                ->Sumw2();
  h_k3pi_pi3_P                ->Sumw2();

  h_k3pi_pi1_E                ->Sumw2();
  h_k3pi_pi2_E                ->Sumw2();
  h_k3pi_pi3_E                ->Sumw2();

  h_k3pi_pi1_Pt                ->Sumw2();
  h_k3pi_pi2_Pt                ->Sumw2();
  h_k3pi_pi3_Pt                ->Sumw2();

  h_k3pi_pi1_m                ->Sumw2();
  h_k3pi_pi2_m                ->Sumw2();
  h_k3pi_pi3_m                ->Sumw2();

  h_k3pi_3track_P                ->Sumw2();
  h_k3pi_3track_P                ->Sumw2();
  h_k3pi_3track_P                ->Sumw2();

  h_k3pi_3track_E                ->Sumw2();
  h_k3pi_3track_E                ->Sumw2();
  h_k3pi_3track_E                ->Sumw2();

  h_k3pi_3track_Pt                ->Sumw2();
  h_k3pi_3track_Pt                ->Sumw2();
  h_k3pi_3track_Pt                ->Sumw2();

  h_k3pi_3track_m                ->Sumw2();
  h_k3pi_3track_m                ->Sumw2();
  h_k3pi_3track_m                ->Sumw2();

  h_e1m                      ->Sumw2();
  h_e2m                      ->Sumw2();

  h_eE                     ->Sumw2();
  h_eP                     ->Sumw2();
  h_ePt                     ->Sumw2();

  h_gE                     ->Sumw2();
  h_gP                     ->Sumw2();
  h_gPt                     ->Sumw2();
  h_gtrkRadDHC1            ->Sumw2();

  h_muchlvE                ->Sumw2();
  h_muchlvP                 ->Sumw2();
  h_muchlvPt                ->Sumw2();
  h_piodd_eovp              ->Sumw2();
  h_elsE                    ->Sumw2();
  h_elsP                     ->Sumw2();
  h_elsPt                    ->Sumw2();

  h_pi0E                    ->Sumw2();
  h_pi0P                     ->Sumw2();
  h_pi0Pt                     ->Sumw2();
  h_pi0Pcm                 ->Sumw2();
  h_mofpi0                 ->Sumw2();
  h_mofpi0diff             ->Sumw2();

  h_mofpi0els                 ->Sumw2();
  h_pi0_dalitz               ->Sumw2();

  h_mueeE                    ->Sumw2();
  h_mueeP                     ->Sumw2();
  h_mueePt                 ->Sumw2();
  h_mueeM                     ->Sumw2();

  h_muee_piee_E                    ->Sumw2();
  h_muee_piee_P                     ->Sumw2();
  h_muee_piee_Pt                 ->Sumw2();
  h_muee_piee_M                     ->Sumw2();


  h_gmueeE                ->Sumw2();
  h_gmueeP                 ->Sumw2();
  h_gmueePt                 ->Sumw2();
  h_gmueeM                 ->Sumw2();

  h_KE                    ->Sumw2();
  h_KP                     ->Sumw2();
  h_KPt                     ->Sumw2();
  h_KPcm                     ->Sumw2();
  h_mofK                     ->Sumw2();

  h_nu_P                    ->Sumw2();
  h_nu_m                    ->Sumw2();
  h_nu_E                   ->Sumw2();
  h_nu_Pt                   ->Sumw2();
  h_missing_mass            ->Sumw2();
  h_missing_mass_kp            ->Sumw2();
  h_missing_mass_km            ->Sumw2();

  h_k3pi_nu_P                    ->Sumw2();
  h_k3pi_nu_m                    ->Sumw2();
  h_k3pi_nu_E                   ->Sumw2();
  h_k3pi_nu_Pt                   ->Sumw2();
  h_k3pi_missing_mass            ->Sumw2();

  h_3trk_invm                ->Sumw2();
  h_eldistDCH1plane              ->Sumw2();
  //    h_eldistDCH1_ee_vz            ->Sumw2();
  //    h_eldistDCH1_pch_vz            ->Sumw2();
  h_eldistLKrplane               ->Sumw2();
  h_gel1distLKrplane             ->Sumw2();
  h_gel2distLKrplane             ->Sumw2();
  h_muche1distDHC1plane        ->Sumw2();
  h_muche2distDHC1plane        ->Sumw2();
  h_muche1distLKrplane            ->Sumw2();
  h_muche2distLKrplane            ->Sumw2();
  h_gmuchdistLKrplane            ->Sumw2();

  h_bx                      ->Sumw2();
  h_by                      ->Sumw2();
  h_bdxdz                     ->Sumw2();
  h_bdydz                     ->Sumw2();
  h_bx_vs_by                  ->Sumw2();
  h_bdxdz_vs_bdydz        ->Sumw2();

  h_bx_much                ->Sumw2();
  h_by_much                ->Sumw2();
  h_bdxdz_much                ->Sumw2();
  h_bdydz_much                ->Sumw2();
  h_bx_vs_by_much            ->Sumw2();
  h_bdxdz_vs_bdydz_much    ->Sumw2();

  h_cda_ee                  ->Sumw2();
  h_vx                      ->Sumw2();
  h_vy                      ->Sumw2();
  h_vz                      ->Sumw2();
  h_vz_initial              ->Sumw2();
  h_vx_vs_vy                  ->Sumw2();

  h_k3pi_vx_12              ->Sumw2();
  h_k3pi_vy_12              ->Sumw2();
  h_k3pi_vz_12              ->Sumw2();
  h_cda_k3pic12             ->Sumw2();
  h_k3pi_vx_13              ->Sumw2();
  h_k3pi_vy_13              ->Sumw2();
  h_k3pi_vz_13              ->Sumw2();
  h_cda_k3pic13             ->Sumw2();
  h_k3pi_vx_23              ->Sumw2();
  h_k3pi_vy_23              ->Sumw2();
  h_k3pi_vz_23              ->Sumw2();
  h_cda_k3pic23             ->Sumw2();
  h_vz_pi12_13_diff         ->Sumw2();
  h_vz_pi12_23_diff	      ->Sumw2();
  h_vz_pi13_23_diff	      ->Sumw2();
  h_vz_e12_sevt_diff        ->Sumw2();
  h_vz_e1mu_sevt_diff       ->Sumw2();
  h_vz_e2mu_sevt_diff       ->Sumw2();
  h_cda_e1much              ->Sumw2();
  h_vx_e1much                  ->Sumw2();
  h_vy_e1much                  ->Sumw2();
  h_vz_e1much                  ->Sumw2();
  h_vx_vs_vy_e1much        ->Sumw2();

  h_cda_e2much              ->Sumw2();
  h_vx_e2much                  ->Sumw2();
  h_vy_e2much                  ->Sumw2();
  h_vz_e2much                  ->Sumw2();


  h_vz_e1much_diff         ->Sumw2();
  h_vz_e2much_diff         ->Sumw2();
  h_vz_mu1mu2_diff         ->Sumw2();

  h_angle_e1e2            ->Sumw2();
  h_angle_e1e2_mu            ->Sumw2();
  h_angle_ph_e1e2            ->Sumw2();
  h_angle_much_pi0            ->Sumw2();
  h_a_e1e2_phe1e2            ->Sumw2();
  h_a_e1e2_muchpi0            ->Sumw2();
  h_a_phe1e2_muchpi0        ->Sumw2();


  h_true_M_0    ->Sumw2();
  h_true_P_0    ->Sumw2();
  h_true_E_0    ->Sumw2();
  h_true_Pt_0    ->Sumw2();
  h_true_M_1    ->Sumw2();
  h_true_P_1    ->Sumw2();
  h_true_E_1    ->Sumw2();
  h_true_Pt_1    ->Sumw2();
  h_true_M_2    ->Sumw2();
  h_true_P_2    ->Sumw2();
  h_true_E_2    ->Sumw2();
  h_true_Pt_2    ->Sumw2();
  h_true_M_3    ->Sumw2();
  h_true_P_3    ->Sumw2();
  h_true_E_3    ->Sumw2();
  h_true_Pt_3    ->Sumw2();
  h_true_M_4    ->Sumw2();
  h_true_P_4    ->Sumw2();
  h_true_E_4    ->Sumw2();
  h_true_Pt_4    ->Sumw2();
  h_true_M_5    ->Sumw2();
  h_true_P_5    ->Sumw2();
  h_true_E_5    ->Sumw2();
  h_true_Pt_5    ->Sumw2();
  h_true_M_6    ->Sumw2();
  h_true_P_6    ->Sumw2();
  h_true_E_6    ->Sumw2();
  h_true_Pt_6    ->Sumw2();
  h_true_M_7    ->Sumw2();
  h_true_P_7    ->Sumw2();
  h_true_E_7    ->Sumw2();
  h_true_Pt_7    ->Sumw2();

  h_mc_muvee_M      ->Sumw2();
  h_mc_muvee_P      ->Sumw2();
  h_mc_muvee_Pt     ->Sumw2();
  h_mc_elsP   ->Sumw2();
  h_mc_elsM   ->Sumw2();

  h_mc_nu_M      ->Sumw2();
  h_mc_nu_M2      ->Sumw2();
  h_mc_nu_P      ->Sumw2();
  h_mc_nu_Pt     ->Sumw2();

  h_mc_pvtx_pi       ->Sumw2();
  h_mc_dvtx_pi       ->Sumw2();
  h_mc_vtxdiff_pi_mu0->Sumw2();
  h_mc_vtxdiff_pi_mu1->Sumw2();
  h_mc_vtxdiff_pi_mu2->Sumw2();

  h_mc_dvtx_pt0    ->Sumw2();
  h_mc_dvtx_pt1    ->Sumw2();
  h_mc_dvtx_pt2    ->Sumw2();
  h_mc_dvtx_pt3    ->Sumw2();
  h_mc_dvtx_pt4    ->Sumw2();
  h_mc_dvtx_pt5    ->Sumw2();
  h_mc_dvtx_pt6    ->Sumw2();
  h_mc_dvtx_pt7    ->Sumw2();

  h_mc_Npart_tr ->Sumw2();
  h_mc_pt4->Sumw2();
  h_mc_pt5->Sumw2();
  h_mc_pt6->Sumw2();
  h_mc_pt7->Sumw2();

  //    h_angle_pch_pi0_Kcm        ->Sumw2();


  replicateTH1I( h_mc_Npart_tr                        , cutsize -1 );
  replicateTH1I(h_mc_pt4                         , cutsize -1 );
  replicateTH1I(h_mc_pt5                         , cutsize -1 );
  replicateTH1I(h_mc_pt6                         , cutsize -1 );
  replicateTH1I(h_mc_pt7                         , cutsize -1 );
  replicateTH1I( h_ntrack                        , cutsize -1 );
  replicateTH1I( h_nclust                        , cutsize -1 );

  replicateTH1I( h_evttypebytrks                , cutsize -1 );
  replicateTH1I( h_k3pi_evttypebytrks                , cutsize -1 );
  replicateTH1I( h_nelectrack                    , cutsize -1 );
  replicateTH1I( h_mutrack                       , cutsize -1 );
  replicateTH1I( h_k3pi_pitrack                       , cutsize -1 );
  replicateTH1I( h_mutrackq                       , cutsize -1 );
  replicateTH1I( h_e1trackq                       , cutsize -1 );
  replicateTH1I( h_e2trackq                       , cutsize -1 );
  replicateTH2I( h_neltrk_mutrk                   , cutsize -1 );

  replicateTH1I( h_mu1                           , cutsize -1 );
  replicateTH1I( h_Q1                            , cutsize -1 );
  replicateTH1I( h_Q2                            , cutsize -1 );
  replicateTH1I( h_Q1_NAKL                       , cutsize -1 );
  replicateTH1I( h_nrun                       , cutsize -1 );

  replicateTH1I(h_1vtx_mcut                           , cutsize -1 );
  replicateTH1I(h_2vtx                                , cutsize -1 );
  replicateTH1I(h_CPRE                                , cutsize -1 );
  replicateTH1I(h_k3pi_1vtx_mcut                      , cutsize -1 );
  replicateTH1I(h_k3pi_2vtx                           , cutsize -1 );
  replicateTH1I(h_k3pi_CPRE                           , cutsize -1 );
  replicateTH1I(h_k3pi_full_trig                           , cutsize -1 );
  replicateTH1I(h_full_trig                      , cutsize -1 );


  replicateTH1I( h_k3pi_pi1trackq                       , cutsize -1 );
  replicateTH1I( h_k3pi_pi2trackq                       , cutsize -1 );
  replicateTH1I( h_k3pi_pi3trackq                       , cutsize -1 );

  replicateTH1F( h_dtrkcl                        , cutsize -1 );
  replicateTH1F( h_el1trkhodtime                , cutsize -1 );
  replicateTH1F( h_el1trktime                    , cutsize -1 );
  replicateTH1F( h_el1trktimediff               , cutsize -1 );
  replicateTH1F( h_el2trkhodtime                , cutsize -1 );
  replicateTH1F( h_el2trktime                    , cutsize -1 );
  replicateTH1F( h_el2trktimediff               , cutsize -1 );
  replicateTH1F( h_e1trkcltimediff              , cutsize -1 );
  replicateTH1F( h_e2trkcltimediff              , cutsize -1 );

  replicateTH1F( h_ecltime1                    , cutsize -1 );
  replicateTH1F( h_ecltime2                    , cutsize -1 );
  replicateTH1F( h_gcltime                       , cutsize -1 );

  replicateTH1F( h_e1e2timediff                   , cutsize -1 );
  replicateTH1F( h_e1e2cltimediff                , cutsize -1 );

  replicateTH1F( h_e1e2hodtimediff                   , cutsize -1 );
  replicateTH1F( h_e1muhodtimediff                   , cutsize -1 );
  replicateTH1F( h_e2muhodtimediff                   , cutsize -1 );

  replicateTH1F( h_pi12hodtimediff                   , cutsize -1 );
  replicateTH1F( h_pi13hodtimediff                   , cutsize -1 );
  replicateTH1F( h_pi23hodtimediff                   , cutsize -1 );

  replicateTH1F( h_pi1trktime                   , cutsize -1 );
  replicateTH1F( h_pi1trkhodtime                   , cutsize -1 );
  replicateTH1F( h_pi1trktimediff                   , cutsize -1 );
  replicateTH1F( h_ecltimepi1                   , cutsize -1 );

  replicateTH1F( h_pi2trktime                   , cutsize -1 );
  replicateTH1F( h_pi2trkhodtime                   , cutsize -1 );
  replicateTH1F( h_pi2trktimediff                   , cutsize -1 );
  replicateTH1F( h_ecltimepi2                   , cutsize -1 );

  replicateTH1F( h_pi3trktime                   , cutsize -1 );
  replicateTH1F( h_pi3trkhodtime                   , cutsize -1 );
  replicateTH1F( h_pi3trktimediff                   , cutsize -1 );
  replicateTH1F( h_ecltimepi3                   , cutsize -1 );

  replicateTH1F( h_e1gtimediff                   , cutsize -1 );
  replicateTH1F( h_e2gtimediff                  , cutsize -1 );

  replicateTH1F( h_muchtrktime                , cutsize -1 );
  replicateTH1F( h_muchcltime                , cutsize -1 );
  replicateTH1F( h_muchtrkhodtime            , cutsize -1 );
  replicateTH1F( h_muchtrkcltimediff         , cutsize -1 );
  replicateTH1F( h_muchtrktimediff         , cutsize -1 );

  replicateTH1F( h_eop                      , cutsize -1 );
  replicateTH1F( h_eop_el                      , cutsize -1 );
  replicateTH1F( h_eop_mu                      , cutsize -1 );

  replicateTH1F( h_eop_pi_track                      , cutsize -1 );
  replicateTH1F( h_eop_pi1_track                      , cutsize -1 );
  replicateTH1F( h_eop_pi2_track                      , cutsize -1 );
  replicateTH1F( h_eop_pi3_track                      , cutsize -1 );

  replicateTH2F( h_eop_vs_p                  , cutsize -1 );
  replicateTH2F( h_eop_el_vs_p                  , cutsize -1 );
  replicateTH2F( h_eop_mu_vs_p                  , cutsize -1 );
  replicateTH2F( h_eop_vs_Pt                  , cutsize -1 );
  replicateTH2F( h_eop_vs_piodd                  , cutsize -1 );
  replicateTH1F( h_e1p                      , cutsize -1 );
  replicateTH1F( h_e2p                      , cutsize -1 );

  replicateTH1F( h_k3pi_pi1_P                      , cutsize -1 );
  replicateTH1F( h_k3pi_pi2_P                      , cutsize -1 );
  replicateTH1F( h_k3pi_pi3_P                      , cutsize -1 );

  replicateTH1F( h_k3pi_pi1_E                      , cutsize -1 );
  replicateTH1F( h_k3pi_pi2_E                      , cutsize -1 );
  replicateTH1F( h_k3pi_pi3_E                      , cutsize -1 );

  replicateTH1F( h_k3pi_pi1_Pt                      , cutsize -1 );
  replicateTH1F( h_k3pi_pi2_Pt                      , cutsize -1 );
  replicateTH1F( h_k3pi_pi3_Pt                      , cutsize -1 );

  replicateTH1F( h_k3pi_pi1_m                      , cutsize -1 );
  replicateTH1F( h_k3pi_pi2_m                      , cutsize -1 );
  replicateTH1F( h_k3pi_pi3_m                      , cutsize -1 );

  replicateTH1F( h_k3pi_3track_P                      , cutsize -1 );
  replicateTH1F( h_k3pi_3track_P                      , cutsize -1 );
  replicateTH1F( h_k3pi_3track_P                      , cutsize -1 );

  replicateTH1F( h_k3pi_3track_E                      , cutsize -1 );
  replicateTH1F( h_k3pi_3track_E                      , cutsize -1 );
  replicateTH1F( h_k3pi_3track_E                      , cutsize -1 );

  replicateTH1F( h_k3pi_3track_Pt                      , cutsize -1 );
  replicateTH1F( h_k3pi_3track_Pt                      , cutsize -1 );
  replicateTH1F( h_k3pi_3track_Pt                      , cutsize -1 );

  replicateTH1F( h_k3pi_3track_m                      , cutsize -1 );
  replicateTH1F( h_k3pi_3track_m                      , cutsize -1 );
  replicateTH1F( h_k3pi_3track_m                      , cutsize -1 );

  replicateTH1F( h_e1m                      , cutsize -1 );
  replicateTH1F( h_e2m                      , cutsize -1 );
  replicateTH1F( h_eE                     , cutsize -1 );
  replicateTH1F( h_eP                     , cutsize -1 );
  replicateTH1F( h_ePt                     , cutsize -1 );

  replicateTH1F( h_gE                     , cutsize -1 );
  replicateTH1F( h_gP                     , cutsize -1 );
  replicateTH1F( h_gPt                     , cutsize -1 );
  replicateTH1F( h_gtrkRadDHC1            , cutsize -1 );

  replicateTH1F( h_muchlvE                    , cutsize -1 );
  replicateTH1F( h_muchlvP                     , cutsize -1 );
  replicateTH1F( h_muchlvPt                , cutsize -1 );
  replicateTH1F( h_piodd_eovp                , cutsize -1 );
  replicateTH1F( h_elsE                    , cutsize -1 );
  replicateTH1F( h_elsP                     , cutsize -1 );
  replicateTH1F( h_elsPt                    , cutsize -1 );

  replicateTH1F( h_pi0E                    , cutsize -1 );
  replicateTH1F( h_pi0P                     , cutsize -1 );
  replicateTH1F( h_pi0Pt                     , cutsize -1 );
  replicateTH1F( h_pi0Pcm                     , cutsize -1 );
  replicateTH1F( h_mofpi0                     , cutsize -1 );
  replicateTH1F( h_mofpi0diff                 , cutsize -1 );

  replicateTH1F( h_mofpi0els                 , cutsize -1 );
  replicateTH2F( h_pi0_dalitz              , cutsize -1 );

  replicateTH1F( h_mueeE                    , cutsize -1 );
  replicateTH1F( h_mueeP                     , cutsize -1 );
  replicateTH1F( h_mueePt                     , cutsize -1 );
  replicateTH1F( h_mueeM                     , cutsize -1 );

  replicateTH1F( h_muee_piee_E                    , cutsize -1 );
  replicateTH1F( h_muee_piee_P                     , cutsize -1 );
  replicateTH1F( h_muee_piee_Pt                     , cutsize -1 );
  replicateTH1F( h_muee_piee_M                     , cutsize -1 );


  replicateTH1F( h_gmueeE                    , cutsize -1 );
  replicateTH1F( h_gmueeP                     , cutsize -1 );
  replicateTH1F( h_gmueePt                 , cutsize -1 );
  replicateTH1F( h_gmueeM                     , cutsize -1 );

  replicateTH1F( h_KE                        , cutsize -1 );
  replicateTH1F( h_KP                         , cutsize -1 );
  replicateTH1F( h_KPt                     , cutsize -1 );
  replicateTH1F( h_KPcm                     , cutsize -1 );
  replicateTH1F( h_mofK                     , cutsize -1 );

  replicateTH1F( h_nu_P                    , cutsize -1 );
  replicateTH1F( h_nu_m                    , cutsize -1 );
  replicateTH1F( h_nu_E                     , cutsize -1 );
  replicateTH1F( h_nu_Pt                     , cutsize -1 );

  replicateTH1F( h_missing_mass             , cutsize - 1);
  replicateTH1F( h_missing_mass_kp             , cutsize - 1);
  replicateTH1F( h_missing_mass_km             , cutsize - 1);

  replicateTH1F( h_k3pi_nu_P                    , cutsize -1 );
  replicateTH1F( h_k3pi_nu_m                    , cutsize -1 );
  replicateTH1F( h_k3pi_nu_E                     , cutsize -1 );
  replicateTH1F( h_k3pi_nu_Pt                     , cutsize -1 );

  replicateTH1F( h_k3pi_missing_mass             , cutsize - 1);

  replicateTH1F( h_3trk_invm             , cutsize - 1);
  replicateTH1F( h_eldistDCH1plane          , cutsize -1 );
  //    replicateTH2F( h_eldistDCH1_ee_vz        , cutsize -1 );
  //    replicateTH2F( h_eldistDCH1_pch_vz        , cutsize -1 );
  replicateTH1F( h_eldistLKrplane           , cutsize -1 );
  replicateTH1F( h_gel1distLKrplane         , cutsize -1 );
  replicateTH1F( h_gel2distLKrplane         , cutsize -1 );
  replicateTH1F( h_muche1distDHC1plane        , cutsize -1 );
  replicateTH1F( h_muche2distDHC1plane        , cutsize -1 );
  replicateTH1F( h_muche1distLKrplane         , cutsize -1 );
  replicateTH1F( h_muche2distLKrplane         , cutsize -1 );
  replicateTH1F( h_gmuchdistLKrplane        , cutsize -1 );

  replicateTH1F( h_bx                      , cutsize -1 );
  replicateTH1F( h_by                      , cutsize -1 );
  replicateTH1F( h_bdxdz                     , cutsize -1 );
  replicateTH1F( h_bdydz                     , cutsize -1 );
  replicateTH2F( h_bx_vs_by                  , cutsize -1 );
  replicateTH2F( h_bdxdz_vs_bdydz            , cutsize -1 );

  replicateTH1F( h_bx_much                    , cutsize -1 );
  replicateTH1F( h_by_much                    , cutsize -1 );
  replicateTH1F( h_bdxdz_much                , cutsize -1 );
  replicateTH1F( h_bdydz_much                , cutsize -1 );
  replicateTH2F( h_bx_vs_by_much            , cutsize -1 );
  replicateTH2F( h_bdxdz_vs_bdydz_much        , cutsize -1 );

  replicateTH1F( h_cda_ee                  , cutsize -1 );
  replicateTH1F( h_vx                      , cutsize -1 );
  replicateTH1F( h_vy                      , cutsize -1 );
  replicateTH1F( h_vz                      , cutsize -1 );
  replicateTH1F( h_vz_initial              , cutsize -1 );
  replicateTH2F( h_vx_vs_vy                  , cutsize -1 );


  replicateTH1F( h_k3pi_vx_12                   , cutsize -1 );
  replicateTH1F( h_k3pi_vy_12                   , cutsize -1 );
  replicateTH1F( h_k3pi_vz_12                   , cutsize -1 );
  replicateTH1F( h_cda_k3pic12                   , cutsize -1 );

  replicateTH1F( h_k3pi_vx_13                   , cutsize -1 );
  replicateTH1F( h_k3pi_vy_13                   , cutsize -1 );
  replicateTH1F( h_k3pi_vz_13                   , cutsize -1 );
  replicateTH1F( h_cda_k3pic13                   , cutsize -1 );

  replicateTH1F( h_k3pi_vx_23                   , cutsize -1 );
  replicateTH1F( h_k3pi_vy_23                   , cutsize -1 );
  replicateTH1F( h_k3pi_vz_23                   , cutsize -1 );
  replicateTH1F( h_cda_k3pic23                   , cutsize -1 );

  replicateTH1F( h_vz_pi12_13_diff                 , cutsize -1 );
  replicateTH1F( h_vz_pi12_23_diff                 , cutsize -1 );
  replicateTH1F( h_vz_pi13_23_diff                 , cutsize -1 );

  replicateTH1F( h_vz_e12_sevt_diff              , cutsize -1 );
  replicateTH1F( h_vz_e1mu_sevt_diff             , cutsize -1 );
  replicateTH1F( h_vz_e2mu_sevt_diff              , cutsize -1 );


  replicateTH1F( h_cda_e1much                  , cutsize -1 );
  replicateTH1F( h_vx_e1much                  , cutsize -1 );
  replicateTH1F( h_vy_e1much                  , cutsize -1 );
  replicateTH1F( h_vz_e1much                  , cutsize -1 );
  replicateTH2F( h_vx_vs_vy_e1much         , cutsize -1 );

  replicateTH1F( h_cda_e2much                  , cutsize -1 );
  replicateTH1F( h_vx_e2much                  , cutsize -1 );
  replicateTH1F( h_vy_e2much                  , cutsize -1 );
  replicateTH1F( h_vz_e2much                  , cutsize -1 );

  replicateTH1F( h_vz_e1much_diff             , cutsize -1 );
  replicateTH1F( h_vz_e2much_diff             , cutsize -1 );
  replicateTH1F( h_vz_mu1mu2_diff             , cutsize -1 );

  replicateTH1F( h_angle_e1e2             , cutsize -1 );
  replicateTH1F( h_angle_e1e2_mu             , cutsize -1 );
  replicateTH1F( h_angle_ph_e1e2            , cutsize -1 );
  replicateTH1F( h_angle_much_pi0            , cutsize -1 );
  replicateTH2F( h_a_e1e2_phe1e2            , cutsize -1 );
  replicateTH2F( h_a_e1e2_muchpi0            , cutsize -1 );
  replicateTH2F( h_a_phe1e2_muchpi0        , cutsize -1 );


  replicateTH1F( h_true_M_0            , cutsize -1 );
  replicateTH1F( h_true_P_0            , cutsize -1 );
  replicateTH1F( h_true_E_0            , cutsize -1 );
  replicateTH1F( h_true_Pt_0           , cutsize -1 );
  replicateTH1F( h_true_M_1            , cutsize -1 );
  replicateTH1F( h_true_P_1            , cutsize -1 );
  replicateTH1F( h_true_E_1            , cutsize -1 );
  replicateTH1F( h_true_Pt_1           , cutsize -1 );
  replicateTH1F( h_true_M_2            , cutsize -1 );
  replicateTH1F( h_true_P_2            , cutsize -1 );
  replicateTH1F( h_true_E_2            , cutsize -1 );
  replicateTH1F( h_true_Pt_2           , cutsize -1 );
  replicateTH1F( h_true_M_3            , cutsize -1 );
  replicateTH1F( h_true_P_3            , cutsize -1 );
  replicateTH1F( h_true_E_3            , cutsize -1 );
  replicateTH1F( h_true_Pt_3           , cutsize -1 );
  replicateTH1F( h_true_M_4            , cutsize -1 );
  replicateTH1F( h_true_P_4            , cutsize -1 );
  replicateTH1F( h_true_E_4            , cutsize -1 );
  replicateTH1F( h_true_Pt_4           , cutsize -1 );
  replicateTH1F( h_true_M_5            , cutsize -1 );
  replicateTH1F( h_true_P_5            , cutsize -1 );
  replicateTH1F( h_true_E_5            , cutsize -1 );
  replicateTH1F( h_true_Pt_5           , cutsize -1 );
  replicateTH1F( h_true_M_6            , cutsize -1 );
  replicateTH1F( h_true_P_6            , cutsize -1 );
  replicateTH1F( h_true_E_6            , cutsize -1 );
  replicateTH1F( h_true_Pt_6           , cutsize -1 );
  replicateTH1F( h_true_M_7            , cutsize -1 );
  replicateTH1F( h_true_P_7            , cutsize -1 );
  replicateTH1F( h_true_E_7            , cutsize -1 );
  replicateTH1F( h_true_Pt_7           , cutsize -1 );

  replicateTH1F(  h_mc_muvee_M                , cutsize -1 );
  replicateTH1F(  h_mc_muvee_P                , cutsize -1 );
  replicateTH1F(  h_mc_muvee_Pt               , cutsize -1 );
  replicateTH1F(  h_mc_elsP             , cutsize -1 );
  replicateTH1F(  h_mc_elsM             , cutsize -1 );

  replicateTH1F(  h_mc_nu_M                , cutsize -1 );
  replicateTH1F(  h_mc_nu_M2                , cutsize -1 );
  replicateTH1F(  h_mc_nu_P                , cutsize -1 );
  replicateTH1F(  h_mc_nu_Pt               , cutsize -1 );

  replicateTH1F(   h_mc_pvtx_pi                  , cutsize -1 );
  replicateTH1F(   h_mc_dvtx_pi                   , cutsize -1 );

  replicateTH1F(   h_mc_vtxdiff_pi_mu0            , cutsize -1 );
  replicateTH1F(   h_mc_vtxdiff_pi_mu1           , cutsize -1 );
  replicateTH1F(   h_mc_vtxdiff_pi_mu2           , cutsize -1 );

  replicateTH1F(  h_mc_dvtx_pt0                  , cutsize -1 );
  replicateTH1F(  h_mc_dvtx_pt1                  , cutsize -1 );
  replicateTH1F(  h_mc_dvtx_pt2                  , cutsize -1 );
  replicateTH1F(  h_mc_dvtx_pt3                  , cutsize -1 );
  replicateTH1F(  h_mc_dvtx_pt4                  , cutsize -1 );
  replicateTH1F(  h_mc_dvtx_pt5                  , cutsize -1 );
  replicateTH1F(  h_mc_dvtx_pt6                  , cutsize -1 );
  replicateTH1F(  h_mc_dvtx_pt7                 , cutsize -1 );

  replicateTH1F(  h_mc_pvtx_pt0                  , cutsize -1 );
  replicateTH1F(  h_mc_pvtx_pt1                  , cutsize -1 );
  replicateTH1F(  h_mc_pvtx_pt2                  , cutsize -1 );
  replicateTH1F(  h_mc_pvtx_pt3                  , cutsize -1 );
  replicateTH1F(  h_mc_pvtx_pt4                  , cutsize -1 );
  replicateTH1F(  h_mc_pvtx_pt5                  , cutsize -1 );
  replicateTH1F(  h_mc_pvtx_pt6                  , cutsize -1 );
  replicateTH1F(  h_mc_pvtx_pt7                 , cutsize -1 );

  //    replicateTH1F( h_angle_pch_pi0_Kcm        , cutsize -1 );


  noBursts = 0;
  burstCounter = 0;




  /*----------- End of user C code -----------*/
  return 0;
}



void replicateTH1I(TH1I*h, int replicas)    {

    file1->cd();

    replicas = (replicas>dirs.size()) ? 0 : dirs.size() - replicas;

    for (unsigned int i = dirs.size(); i>replicas; i--)     {

        dirs[i-1]->cd();

        TH1I* h_clone = new TH1I(*h);

        file1->cd();
    }
}


void replicateTH2I(TH2I*h, int replicas)    {

    file1->cd();

    replicas = (replicas>dirs.size()) ? 0 : dirs.size() - replicas;

    for (unsigned int i = dirs.size(); i>replicas; i--)     {

        dirs[i-1]->cd();

        TH2I* h_clone = new TH2I(*h);

        file1->cd();
    }
}


void replicateTH1F(TH1F*h, int replicas)    {

    file1->cd();

    replicas = (replicas>dirs.size()) ? 0 : dirs.size() - replicas;

    for (unsigned int i = dirs.size(); i>replicas; i--)     {

        dirs[i-1]->cd();

        TH1F* h_clone = new TH1F(*h);

        file1->cd();
    }
}


void replicateTH2F(TH2F*h, int replicas)    {

    file1->cd();

    replicas = (replicas>dirs.size()) ? 0 : dirs.size() - replicas;

    for (unsigned int i = dirs.size(); i>replicas; i--)     {

        dirs[i-1]->cd();

        TH2F* h_clone = new TH2F(*h);

        file1->cd();
    }
}
