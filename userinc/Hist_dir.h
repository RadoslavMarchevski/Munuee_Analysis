#ifndef __Hist_dir_h_
#define  __Hist_dir_h_
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include <iostream>
#include <string>
#include "Charged_Particle.h"
<<<<<<< HEAD
#include "MC_Charged_Particle.h"
#include "Cuts.h"

=======
#include "Cuts.h"
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
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
    void FillHist(TLorentzVector Three_Track_Momentum, TLorentzVector Nu_Momentum);
<<<<<<< HEAD
    void FillAngle(TLorentzVector muon, TLorentzVector Two_Track_Momentum);
    //void FillAngle();
=======
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
    void FillVertexHist(double* mu_e1, double cda_mu_e1,double* mu_e2, double cda_mu_e2, double* e1_e2, double cda_e1_e2,std::string decay_type);
    void ComputeThreeTrack(Charged_Particle& p1,Charged_Particle& p2, Charged_Particle& p3);
    TLorentzVector& GetThreeTrackMomentum(){return Three_Track_Momentum; };
    TLorentzVector& GetTwoTrackMomentum(){return Two_Track_Momentum; };
    TLorentzVector& GetKaonMomentum(){return Kaon_Momentum; };
<<<<<<< HEAD
    TLorentzVector& GetMuNuMomentum(){return MuNu_Momentum; };
    TLorentzVector& GetNuMomentum(){return Nu_Momentum; };
=======
    TLorentzVector& GetNuMomentum(){return Nu_Momentum; };

>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
    TLorentzVector Kaon_Momentum;
    TLorentzVector Two_Track_Momentum;
    TLorentzVector Three_Track_Momentum;
    TLorentzVector Nu_Momentum;
<<<<<<< HEAD
    TLorentzVector MuNu_Momentum;
=======
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb

    TH1I* fh_Ntracks           ;
    TH1I* fh_Nclusters         ;
    TH1I* fh_Nvtx              ;
    TH1I* fh_Kaon_Charge       ;
    TH1I* fh_Mu_charge         ;
    TH1I* fh_Event_Type        ;
    TH1F* fh_Track_Momentum    ;
    TH1F* fh_Pion_Momentum     ;
    TH1F* fh_Mu_momentum       ;
    TH1F* fh_Electron_Momentum ;
    TH1F* fh_eop               ;
    TH1F* fh_Mu_eop            ;
    TH1F* fh_Electron_eop      ;
    TH1F* fh_odd_eop           ;
    TH1F* fh_COmPaCt_Z_Vertex  ;
    TH1F* fh_Mu_Zvtx_min_COmPaCt_Zvtx;
    TH1F* fh_xvtxdiff_mue1_mue2;
    TH1F* fh_xvtxdiff_mue1_e1e2;
    TH1F* fh_xvtxdiff_mue2_e1e2;
    TH1F* fh_yvtxdiff_mue1_mue2;
    TH1F* fh_yvtxdiff_mue1_e1e2;
    TH1F* fh_yvtxdiff_mue2_e1e2;
    TH1F* fh_zvtxdiff_mue1_mue2;
    TH1F* fh_zvtxdiff_mue1_e1e2;
    TH1F* fh_zvtxdiff_mue2_e1e2;
    TH1F* fh_xvtxdiff_pi1pi2_pi2pi3;
    TH1F* fh_xvtxdiff_pi1pi2_pi1pi3;
    TH1F* fh_xvtxdiff_pi1pi3_pi2pi3;
    TH1F* fh_yvtxdiff_pi1pi2_pi2pi3;
    TH1F* fh_yvtxdiff_pi1pi2_pi1pi3;
    TH1F* fh_yvtxdiff_pi1pi3_pi2pi3;
    TH1F* fh_zvtxdiff_pi1pi2_pi2pi3;
    TH1F* fh_zvtxdiff_pi1pi2_pi1pi3;
    TH1F* fh_zvtxdiff_pi1pi3_pi2pi3;
    TH1F* fh_cda_pi1_pi2;
    TH1F* fh_cda_pi1_pi3;
    TH1F* fh_cda_pi2_pi3;
    TH1F* fh_DCHtime_mu;
    TH1F* fh_DCHtime_e1;
    TH1F* fh_DCHtime_e2;
    TH1F* fh_HodTime_mu;
    TH1F* fh_HodTime_e1;
    TH1F* fh_HodTime_e2;
    TH1F* fh_DCH_timediff_mu_e1;
    TH1F* fh_DCH_timediff_mu_e2;
    TH1F* fh_DCH_timediff_e1_e2;
    TH1F* fh_Hod_timediff_mu_e2;
    TH1F* fh_Hod_timediff_e1_e2;
    TH1F* fh_Hod_timediff_mu_e1;
    TH1F* fh_mee;
    TH1F* fh_muee_M;
    TH1F* fh_missing_mass;
    TH1F* fh_muee_Pt;
    TH1F* fh_muee_P;
    TH1F* fh_cda_mu_e2;
    TH1F* fh_cda_e1_e2;
    TH1F* fh_cda_mu_e1;
    TH1F* fh_muon_bx      ;
    TH1F* fh_muon_by      ;
    TH1F* fh_electron1_bx ;
    TH1F* fh_electron1_by ;
    TH1F* fh_electron2_bx ;
    TH1F* fh_electron2_by ;
    TH1F* fh_pion1_bx     ;
    TH1F* fh_pion1_by     ;
    TH1F* fh_pion2_bx     ;
    TH1F* fh_pion2_by     ;
    TH1F* fh_pion3_bx     ;
    TH1F* fh_pion3_by     ;
    TH2F* fh_bx_vs_by_muon;
    TH2F* fh_bx_vs_by_el1;
    TH2F* fh_bx_vs_by_el2;
    TH2F* fh_bx_vs_by_pi1;
    TH2F* fh_bx_vs_by_pi2;
    TH2F* fh_bx_vs_by_pi3;
<<<<<<< HEAD
    TH1F* fh_el1_cluster_x;
    TH1F* fh_el1_cluster_y;
    TH1F* fh_el2_cluster_x;
    TH1F* fh_el2_cluster_y;
    TH1F* fh_mu_cluster_x;
    TH1F* fh_mu_cluster_y;
    TH1F* fh_el1_cluster_time;
    TH1F* fh_el2_cluster_time;
    TH1F* fh_el1_el2_cluster_timediff;
    TH1F* fh_el1_dtrk_cl;
    TH1F* fh_el2_dtrk_cl;
    TH1F* fh_muon_dtrk_cl;
    TH2F* fh_el1_cluster_x_y;
    TH2F* fh_el2_cluster_x_y;
    TH2F* fh_mu_cluster_x_y;
    TH1F* fh_deadcell_distance;
    TH1F* fh_muv2_trk_cl_diff;
    TH2F* fh_muv_x_y_position;
    TH1F* fh_DCH1_distance_e1e2;
    TH1F* fh_Lkr_distance_e1e2;
    TH1F* fh_angle_el_mu;
    TH1F* fh_mc_KDzvtx;
    TH1F* fh_mc_P1_Pzvtx;
    TH1F* fh_mc_P1_mass;
    TH1F* fh_mc_P1_momentum;
    TH1F* fh_mc_three_track_123_momentum;
    TH1F* fh_mc_P2_mass;
    TH1F* fh_mc_P2_Pzvtx;
    TH1F* fh_mc_P2_momentum;
    TH1F* fh_mc_P3_Pzvtx;
    TH1F* fh_mc_P3_mass;
    TH1F* fh_mc_P3_momentum;
    TH1F* fh_mc_three_track_123_mass;
    TH1F* fh_mc_two_track_23_momentum;
    TH1F* fh_mc_two_track_23_mass;
    TH1F* fh_mc_P4_Pzvtx;
    TH1F* fh_mc_P4_mass;
    TH1F* fh_mc_P4_momentum;
    TH1F* fh_mc_four_track_1234_momentum;
    TH1F* fh_mc_four_track_1234_mass;
    TH1F* fh_missing_mass_z_variable;

=======
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb


private:
    std::string fDir;
    int ftype;
<<<<<<< HEAD
    double mKaon;
=======
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
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
#endif
