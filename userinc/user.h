extern FILE *fprt;
extern char gString[50];

extern char rootfilename[200];
extern char rootFileName[200];


extern float CPDpos_leftDownCorner[256][2];
extern float CELLpos_leftDownCorner[256][64][2];
extern int   CPDindex, CELLindex;
extern float CPDlength, CELLlength;

extern int runNo;
extern int magnetCurrent;
extern int noBursts;
extern int burstCounter;
extern int IS_DATA, IS_MC;



// missing compact declarations
extern "C" int accep_(int* run, float* x, float* y);
extern "C" int closap_double_(double p1[3], double p2[3], double v1[3], double v2[3], double *dmin, double vtx[3]);
extern "C" void blue_tack_(int *nchar, float *tmom, float Vxyz[3], float vpnt[2], float vdir[2]);

//extern int LKr_acc(int nrun, double pos_x, double pos_y, float par);
extern int LKr_acc(int nrun, float pos_x, float pos_y, float par);        // offizielle Routine in lkraccep_2007.c benutzt ab v26

extern void GetCpdCellIndex(float pos_x, float pos_y, int *cpd_index, int *cell_index);

extern char histoname[200];
extern char histoname2[200];

// E/p corrections for each cell
extern FILE *EopCorrfile;
extern char  EopCorrfilename[200];
extern float EopCorr[100][256][64]; // ab v65: corrections for different (sub-) periods (up to 100)
extern int   periodFlag;            // ab v65: period to be defined in superBurst.c


#define SQR(x) ((x)*(x))

/*************************************************************************/
/*                        root include files                             */
/*************************************************************************/
#include "TROOT.h"
#include "TMath.h"
#include "TH3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TList.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TDirectory.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
using namespace std;
extern superMcEvent * mcevent;

extern TFile *file1;
extern TH1F* HistnoBursts;
extern TH1F* HistnoGoodBursts;
extern TH1F* h_cutflow;
extern TH1I* h_evttypebytrks;
extern TH1I* h_ntrack;
extern TH1I* h_nclust;
extern TH1I* h_nelectrack;
extern TH1I* h_mutrack;
extern TH1I* h_mutrackq;
extern TH1I* h_e1trackq;
extern TH1I* h_e2trackq;
extern TH2I* h_neltrk_mutrk;
//extern TH1I* h_nphcluster;
//extern TH2I* h_neltrk_nphclus;
extern TH1I* h_Q1;
extern TH1I* h_1vtx_mcut;
extern TH1I* h_2vtx;
extern TH1I* h_CPRE;
extern TH1I* h_full_trig;
extern TH1I* h_k3pi_1vtx_mcut;
extern TH1I* h_k3pi_2vtx;
extern TH1I* h_k3pi_CPRE;  
extern TH1I* h_k3pi_full_trig;
extern TH1I* h_mu1;
extern TH1I* h_nrun;
extern TH1I* h_Q2;
extern TH1I* h_Q1_NAKL;
extern TH1F* h_dtrkcl;
extern TH1F* h_muchtrktime;
extern TH1F* h_muchtrkhodtime;
extern TH1F* h_muchcltime;
extern TH1F* h_muchtrkcltimediff;
extern TH1F* h_muchtrktimediff;
extern TH1F* h_el1trkhodtime;
extern TH1F* h_el1trktime;
extern TH1F* h_el1trktimediff;
extern TH1F* h_el2trkhodtime;
extern TH1F* h_el2trktime;
extern TH1F* h_el2trktimediff;
extern TH1F* h_e1trkcltimediff;
extern TH1F* h_e2trkcltimediff;
extern TH1F* h_ecltime1;
extern TH1F* h_ecltime2;
extern TH1F* h_gcltime;
extern TH1F* h_e1e2timediff;
extern TH1F* h_e1e2cltimediff;
extern TH1F* h_e1e2hodtimediff;
extern TH1F* h_e1muhodtimediff;
extern TH1F* h_e2muhodtimediff;
extern TH1F* h_e1gtimediff;
extern TH1F* h_e2gtimediff;

extern TH1F* h_eldistDCH1plane;
//extern //TH2F* h_eldistDCH1_ee_vz;
//extern //TH2F* h_eldistDCH1_pch_vz;
extern TH1F* h_eldistLKrplane;
extern TH1F* h_gel1distLKrplane;
extern TH1F* h_gel2distLKrplane;
extern TH1F* h_muche1distDHC1plane;
extern TH1F* h_muche2distDHC1plane;
extern TH1F* h_muche1distLKrplane;
extern TH1F* h_muche2distLKrplane;
extern TH1F* h_gmuchdistLKrplane;

extern TH1F* h_eE        ;
extern TH1F* h_eP        ;
extern TH1F* h_ePt     ;
extern TH1F* h_elsE    ;
extern TH1F* h_elsPt    ;

extern TH1F* h_gE         ;
extern TH1F* h_gP         ;
extern TH1F* h_gPt     ;
extern TH1F* h_gtrkRadDHC1;

extern TH1F* h_muchlvE    ;
extern TH1F* h_muchlvPt    ;
extern TH1F* h_pi0E    ;
extern TH1F* h_mueeE    ;
extern TH1F* h_gmueeE    ;
extern TH1F* h_KE        ;


extern TH1F* h_nu_P; 
extern TH1F* h_nu_m; 
extern TH1F* h_nu_E; 
extern TH1F* h_nu_Pt;
extern TH1F* h_missing_mass;
extern TH1F* h_missing_mass_kp;
extern TH1F* h_missing_mass_km;
//K3pi charged
extern TH1F* h_k3pi_missing_mass;
extern TH1F* h_k3pi_nu_P; 
extern TH1F* h_k3pi_nu_m; 
extern TH1F* h_k3pi_nu_E; 
extern TH1F* h_k3pi_nu_Pt;

extern TH1I* h_k3pi_evttypebytrks;
extern TH1I* h_k3pi_pitrack;
extern TH1I* h_k3pi_pi1trackq;
//extern TH1F* h_k3pi_pi2track;
extern TH1I* h_k3pi_pi2trackq;
//extern TH1F* h_k3pi_pi3track;
extern TH1I* h_k3pi_pi3trackq;
extern TH1F* h_eop_pi_track;
extern TH1F* h_eop_pi1_track;
extern TH1F* h_eop_pi2_track;
extern TH1F* h_eop_pi3_track;

extern TH1F* h_piodd_eovp    ;
extern TH1F* h_pi1trktime;
extern TH1F* h_pi1trkhodtime;
extern TH1F* h_pi1cltime;
extern TH1F* h_pi1trkcltimediff;
extern TH1F* h_pi1trktimediff;
extern TH1F* h_ecltimepi1;
extern TH1F* h_pi12hodtimediff;
extern TH1F* h_pi23hodtimediff;
extern TH1F* h_pi13hodtimediff;

extern TH1F* h_pi2trkhodtime;
extern TH1F* h_pi2trktime;
extern TH1F* h_pi2trktimediff;
extern TH1F* h_pi2trkcltimediff;
extern TH1F* h_ecltimepi2;

extern TH1F* h_pi3trkhodtime;
extern TH1F* h_pi3trktime;
extern TH1F* h_pi3trktimediff;
extern TH1F* h_pi3trkcltimediff;
extern TH1F* h_ecltimepi3;

extern TH1F* h_k3pi_3track_P; 
extern TH1F* h_k3pi_3track_m; 
extern TH1F* h_k3pi_3track_E; 
extern TH1F* h_k3pi_3track_Pt;
extern TH1F* h_k3pi_pi1_P; 
extern TH1F* h_k3pi_pi1_m; 
extern TH1F* h_k3pi_pi1_E; 
extern TH1F* h_k3pi_pi1_Pt;
extern TH1F* h_k3pi_pi2_P; 
extern TH1F* h_k3pi_pi2_m; 
extern TH1F* h_k3pi_pi2_E; 
extern TH1F* h_k3pi_pi2_Pt;
extern TH1F* h_k3pi_pi3_P; 
extern TH1F* h_k3pi_pi3_m; 
extern TH1F* h_k3pi_pi3_E; 
extern TH1F* h_k3pi_pi3_Pt;

extern TH1F* h_k3pi_vx_12;
extern TH1F* h_k3pi_vy_12;
extern TH1F* h_k3pi_vz_12;
extern TH1F* h_cda_k3pic12;
extern TH1F* h_k3pi_vx_13;
extern TH1F* h_k3pi_vy_13;
extern TH1F* h_k3pi_vz_13;
extern TH1F* h_cda_k3pic13;
extern TH1F* h_k3pi_vx_23;
extern TH1F* h_k3pi_vy_23;
extern TH1F* h_k3pi_vz_23;
extern TH1F* h_cda_k3pic23;
extern TH1F* h_vz_pi12_13_diff;
extern TH1F* h_vz_pi12_23_diff;
extern TH1F* h_vz_pi13_23_diff;

//end K3pi charged
extern TH1F* h_3trk_invm;
extern TH1F* h_eop;
extern TH1F* h_eop_el;
extern TH1F* h_eop_mu;
extern TH2F* h_eop_vs_p;
extern TH2F* h_eop_el_vs_p;
extern TH2F* h_eop_mu_vs_p;
extern TH2F* h_eop_vs_Pt;
extern TH2F* h_eop_vs_piodd;
extern TH1F* h_e1p;
extern TH1F* h_e2p;
extern TH1F* h_e1m;
extern TH1F* h_e2m;
extern TH1F* h_elsP;
extern TH1F* h_mofpi0els;
extern TH1F* h_mofpi0;
extern TH2F* h_pi0_dalitz;
extern TH1F* h_pi0P;
extern TH1F* h_pi0Pt;
extern TH1F* h_pi0Pcm;
extern TH1F* h_mofpi0diff;
extern TH1F* h_muchlvP;
extern TH1F* h_KP;
extern TH1F* h_KPt;
extern TH1F* h_KPcm;
extern TH1F* h_mofK;
extern TH1F* h_mueeP;
extern TH1F* h_mueePt;
extern TH1F* h_mueeM;
extern TH1F* h_muee_piee_E;
extern TH1F* h_muee_piee_M;
extern TH1F* h_muee_piee_P;
extern TH1F* h_muee_piee_Pt;
extern TH1F* h_gmueeP;
extern TH1F* h_gmueePt;
extern TH1F* h_gmueeM;
extern TH1F* h_piodd_eovp;
extern TH1F* h_bx;
extern TH1F* h_by;
extern TH1F* h_bdxdz;
extern TH1F* h_bdydz;
extern TH2F* h_bx_vs_by;
extern TH2F* h_bdxdz_vs_bdydz;

extern TH1F* h_bx_much;
extern TH1F* h_by_much;
extern TH1F* h_bdxdz_much;
extern TH1F* h_bdydz_much;
extern TH2F* h_bx_vs_by_much;
extern TH2F* h_bdxdz_vs_bdydz_much;
extern TH2F* h_vx_vs_vy;
extern TH1F* h_vx;
extern TH1F* h_vy;
extern TH1F* h_vz;
extern TH1F* h_vz_initial;
extern TH1F* h_cda_ee;

extern TH2F* h_vx_vs_vy_e1much;
extern TH1F* h_vx_e1much;
extern TH1F* h_vy_e1much;
extern TH1F* h_vz_e1much;

extern TH1F* h_vx_e2much;
extern TH1F* h_vy_e2much;
extern TH1F* h_vz_e2much;

extern TH1F* h_cda_e1much;
extern TH1F* h_cda_e2much;
extern TH1F* h_vz_e1much_diff;
extern TH1F* h_vz_e2much_diff;
extern TH1F* h_vz_mu1mu2_diff;
extern TH1F* h_vz_e12_sevt_diff;
extern TH1F* h_vz_e1mu_sevt_diff;
extern TH1F* h_vz_e2mu_sevt_diff;

extern TH1F* h_angle_e1e2;
extern TH1F* h_angle_e1e2_mu;
extern TH1F* h_angle_ph_e1e2;
extern TH1F* h_angle_much_pi0;
extern TH2F* h_a_e1e2_phe1e2    ;
extern TH2F* h_a_e1e2_muchpi0    ;
extern TH2F* h_a_phe1e2_muchpi0    ;

extern TH1F* h_angle_much_pi0_Kcm;

extern TH1F* h_true_M_0;
extern TH1F* h_true_P_0;
extern TH1F* h_true_E_0;
extern TH1F* h_true_Pt_0;
extern TH1F* h_true_M_1;
extern TH1F* h_true_P_1;
extern TH1F* h_true_E_1;
extern TH1F* h_true_Pt_1;
extern TH1F* h_true_M_2;
extern TH1F* h_true_P_2;
extern TH1F* h_true_E_2;
extern TH1F* h_true_Pt_2;
extern TH1F* h_true_M_3;
extern TH1F* h_true_P_3;
extern TH1F* h_true_E_3;
extern TH1F* h_true_Pt_3;
extern TH1F* h_true_M_4;
extern TH1F* h_true_P_4;
extern TH1F* h_true_E_4;
extern TH1F* h_true_Pt_4;
extern TH1F* h_true_M_5;
extern TH1F* h_true_P_5;
extern TH1F* h_true_E_5;
extern TH1F* h_true_Pt_5;
extern TH1F* h_true_M_6;
extern TH1F* h_true_P_6;
extern TH1F* h_true_E_6;
extern TH1F* h_true_Pt_6;
extern TH1F* h_true_M_7;
extern TH1F* h_true_P_7;
extern TH1F* h_true_E_7;
extern TH1F* h_true_Pt_7;
extern TH1F* h_mc_elsP;
extern TH1F* h_mc_elsM;
extern TH1F* h_mc_muvee_M;
extern TH1F* h_mc_muvee_P;
extern TH1F* h_mc_muvee_Pt;

extern TH1F* h_mc_nu_M;
extern TH1F* h_mc_nu_P;
extern TH1F* h_mc_nu_Pt;
extern TH1F* h_mc_nu_M2;

extern TH1I* h_mc_Npart_tr;
extern TH1I* h_mc_pt4;
extern TH1I* h_mc_pt5;
extern TH1I* h_mc_pt6;
extern TH1I* h_mc_pt7;

extern TH1F* h_mc_pvtx_pi;
extern TH1F* h_mc_dvtx_pi;
extern TH1F* h_mc_pvtx_pt0;
extern TH1F* h_mc_pvtx_pt1;
extern TH1F* h_mc_pvtx_pt2;
extern TH1F* h_mc_pvtx_pt3;
extern TH1F* h_mc_pvtx_pt4;
extern TH1F* h_mc_pvtx_pt5;
extern TH1F* h_mc_pvtx_pt6;
extern TH1F* h_mc_pvtx_pt7;

extern TH1F* h_mc_dvtx_pt0;
extern TH1F* h_mc_dvtx_pt1;
extern TH1F* h_mc_dvtx_pt2;
extern TH1F* h_mc_dvtx_pt3;
extern TH1F* h_mc_dvtx_pt4;
extern TH1F* h_mc_dvtx_pt5;
extern TH1F* h_mc_dvtx_pt6;
extern TH1F* h_mc_dvtx_pt7;

extern TH1F* h_mc_vtxdiff_pi_mu0;
extern TH1F* h_mc_vtxdiff_pi_mu1;
extern TH1F* h_mc_vtxdiff_pi_mu2;

extern int Npart;
extern int p0Type;
extern int p1Type;
extern int p2Type;
extern int p3Type;
extern int p4Type;
extern int p5Type;
extern int p6Type;
extern int p7Type;
extern int first_Nevt;
extern int second_Nevt;
extern int nEvt_counter;
extern int rep_counter;

extern float  zVertexTrue;

extern int f_Ntrack;
extern int s_Ntrack;
extern int f_NCl;
extern int s_NCl;
extern int counter_ofmuons;
extern float p0E ;
extern float p0Px;
extern float p0Py;
extern float p0Pz;
extern float pkgen ;
   
extern float p1E ;
extern float p1Px;
extern float p1Py;
extern float p1Pz;
    
extern float p2E ;
extern float p2Px;
extern float p2Py;
extern float p2Pz;

extern float p3E ;
extern float p3Px;
extern float p3Py;
extern float p3Pz;

extern float p4E ;
extern float p4Px;
extern float p4Py;
extern float p4Pz;

extern float p5E ;
extern float p5Px;
extern float p5Py;
extern float p5Pz;

extern float p6E ;
extern float p6Px;
extern float p6Py;
extern float p6Pz;

extern float p7E ;
extern float p7Px;
extern float p7Py;
extern float p7Pz;

extern float p0m ;
extern float p1m ;
extern float p2m ;
extern float p3m ;
extern float p4m ;
extern float p5m ;
extern float p6m ;
extern float p7m ;


extern TLorentzVector K_true;
extern TLorentzVector Pi0_true;
extern TLorentzVector El_true;
extern TLorentzVector Nu_true;
extern TLorentzVector true_0;
extern TLorentzVector true_1;
extern TLorentzVector true_2;
extern TLorentzVector true_3;
extern TLorentzVector true_4;
extern TLorentzVector true_5;
extern TLorentzVector true_6;
extern TLorentzVector true_7;
extern TLorentzVector Muee_three_track;
extern TLorentzVector mc_ee;
extern TLorentzVector *true_V[20];
//extern TLorentzVector v_first_e1;
//extern TLorentzVector v_second_e1;
//extern TLorentzVector v_first_e2;
//extern TLorentzVector v_second_e2;
//extern TLorentzVector v_first_mu;
//extern TLorentzVector v_second_mu;
extern int pType[20];
extern float pvtx[20];
extern float dvtx[20];
extern float k_dvtx[3];
//extern vector<int> pType;
//extern vector<float> pvtx;
//extern vector<float> dvtx;
//extern vector<float> mu_pvtx;
extern vector<string> cuts;     
extern vector<TDirectory*> dirs;


/*
 *
 *
 * $Log: user.h,v $
 * Revision 1.2  2003/10/31 12:32:30  andrew
 * Added the -string option to main
 *
 * An arbitrary string can be passed to compact which is
 * saved in a global variable gString (C), COMMON/GSTRING/GSTRING (FORTRAN)
 *
 *
 * made -ndb the default. For people needin the compact database the -db option
 * was created
 *
 *
 *
 */
