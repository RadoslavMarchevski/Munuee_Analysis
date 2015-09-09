#ifndef __Hist_dir_h_
#define  __Hist_dir_h_
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include <iostream>
#include <string>
#include "Charged_Particle.h"
#include "Neutral_Particle.h"
#include "MC_Charged_Particle.h"
#include "Cuts_pi0d.h"

class Hist_dir {
 public:
    //Three type of directories one for monitoring of some main
    //analysis quantities and two for the two different analyses:
    //K->mu nu e e and the normalization channel K3pi
    Hist_dir(const std::string& dir_Name,int type);
    ~Hist_dir();
    void AddToFile(TFile* file);
    void FillCommonHist(superCmpEvent* sevt);
    void FillHist(Charged_Particle& part, std::string particle);
    void FillHist(Charged_Particle& p1,Charged_Particle& p2, std::string particles);
    void FillVertexHist(double* mu_e1, double cda_mu_e1,double* mu_e2, double cda_mu_e2, double* e1_e2, double cda_e1_e2,std::string decay_type);
    TLorentzVector& GetThreeTrackMomentum(){return Three_Track_Momentum; };
    TLorentzVector& GetTwoTrackMomentum(){return Two_Track_Momentum; };
    TLorentzVector& GetKaonMomentum(){return Kaon_Momentum; };
    TLorentzVector& GetMuNuMomentum(){return MuNu_Momentum; };
    TLorentzVector& GetNuMomentum(){return Nu_Momentum; };
    TLorentzVector Kaon_Momentum;
    TLorentzVector Two_Track_Momentum;
    TLorentzVector Three_Track_Momentum;
    TLorentzVector Nu_Momentum;
    TLorentzVector MuNu_Momentum;

    TH1I* fh_Ntracks           ;
    TH1I* fh_Nclusters         ;
    TH1I* fh_Nvtx              ;
    TH1I* fh_Kaon_Charge       ;
    TH1I* fh_el1_Charge        ;
    TH1I* fh_el2_Charge        ;
    TH1I* fh_pion_Charge       ;
    TH1F* fh_gamma_momentum    ;
    TH1F* fh_pi_momentum       ;
    TH1F* fh_el1_momentum      ;
    TH1F* fh_el2_momentum      ;
    TH1F* fh_k2pi0d_P          ;
    TH1F* fh_k2pi0d_Pt         ;
    TH1F* fh_k2pi0d_M          ;
    TH1F* fh_pi0d_P            ;
    TH1F* fh_pi0d_Pt           ;
    TH1F* fh_pi0d_M            ;
    TH1F* fh_EoP_el1 ;
    TH1F* fh_EoP_el2 ;
    TH1F* fh_EoP_pion;
    TH1F* fh_lda3_e1_plus      ;
    TH1F* fh_lda3_e1_minus     ;
    TH1F* fh_lda3_e1           ;
    TH1F* fh_lda3_e2           ;
    TH1I* fh_el1_plus_Charge   ;
    TH1F* fh_EoP_el1_plus      ;
    TH2F* fh_EoP_vs_p_el1_plus ;
    TH2F* fh_lda3_vs_p_e1_plus ;
    TH1I* fh_el1_minus_Charge  ;
    TH1F* fh_EoP_el1_minus     ;
    TH2F* fh_EoP_vs_p_el1_minus;
    TH2F* fh_EoP_vs_p_el1      ;
    TH2F* fh_lda3_vs_p_e1      ;
    TH2F* fh_lda3_vs_p_e1_minus;
    TH1F* fh_DCH_timediff_pi_e1;
    TH1F* fh_DCH_timediff_pi_e2;
    TH1F* fh_DCH_timediff_e1_e2;
    TH1F* fh_Hod_timediff_pi_e1;
    TH1F* fh_Hod_timediff_pi_e2;
    TH1F* fh_Hod_timediff_e1_e2;
    TH1F* fh_el1_el2_cltimediff;
    TH1F* fh_el1_g_cltimediff  ;
    TH1F* fh_el2_g_cltimediff  ;
    TH1F* fh_DCH1_distance_e1e2;
    TH1F* fh_Lkr_distance_e1e2 ;
    TH2F* fh_Lkr_extrap_x_vs_y ;
    TH1F* fh_xvtxdiff_pie1_pie2;
    TH1F* fh_xvtxdiff_pie1_e1e2;
    TH1F* fh_xvtxdiff_pie2_e1e2;
    TH1F* fh_yvtxdiff_pie1_pie2;
    TH1F* fh_yvtxdiff_pie1_e1e2;
    TH1F* fh_yvtxdiff_pie2_e1e2;
    TH1F* fh_zvtxdiff_pie1_pie2;
    TH1F* fh_zvtxdiff_pie1_e1e2;
    TH1F* fh_zvtxdiff_pie2_e1e2;
    TH1F* fh_cda_pi_e1         ;
    TH1F* fh_cda_pi_e2         ;
    TH1F* fh_cda_e1_e2         ;

 private:
    std::string fDir;
    int ftype;
    double mKaon;
};
extern Hist_dir* Initial_dir;
extern Hist_dir* K3pi_selection;
extern Hist_dir* dir1;
extern Hist_dir* dir2;
extern Hist_dir* dir3;
extern Hist_dir* dir4;
extern Hist_dir* dir5;
extern Hist_dir* dir6;
extern Hist_dir* dir7;
extern Hist_dir* dir8;
extern Hist_dir* dir9;
extern Hist_dir* dir10;
extern Hist_dir* dir11;
#endif
