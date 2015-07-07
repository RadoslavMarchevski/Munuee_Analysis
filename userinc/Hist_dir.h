#ifndef __Hist_dir_h_
#define  __Hist_dir_h_
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include <iostream>
#include <string>
#include "Charged_Particle.h"
#include "Cuts.h"
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
    void FillHist(Charged_Particle& p1,Charged_Particle& p2, Charged_Particle& p3);
    void FillVertexHist(double* mu_e1, double cda_mu_e1,double* mu_e2, double cda_mu_e2, double* e1_e2, double cda_e1_e2);
    TLorentzVector& GetThreeTrackMomentum(){return Three_Track_Momentum; };
    TLorentzVector& GetTwoTrackMomentum(){return Two_Track_Momentum; };
    TLorentzVector& GetKaonMomentum(){return Kaon_Momentum; };
    TLorentzVector& GetNuMomentum(){return Nu_Momentum; };

    TLorentzVector Kaon_Momentum;
    TLorentzVector Two_Track_Momentum;
    TLorentzVector Three_Track_Momentum;
    TLorentzVector Nu_Momentum;

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
    TH1F* fh_Z_Vertex          ;
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


private:
    std::string fDir;
    int ftype;
};
extern Hist_dir* Initial_dir;
extern Hist_dir* K3pi_selection;
extern Hist_dir* dir1;
extern Hist_dir* dir2;
#endif
