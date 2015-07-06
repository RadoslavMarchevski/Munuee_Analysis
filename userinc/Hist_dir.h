#ifndef __Hist_dir_h_
#define  __Hist_dir_h_
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include <iostream>
#include <string>
#include "Charged_Particle.h"

class Hist_dir {
public:
    //Three type of directories one for monitoring of some main
    //analysis quantities and two for the two different analyses:
    //K->mu nu e e and the normalization channel K3pi
    Hist_dir(const std::string& dir_Name,int type);
    ~Hist_dir();
    void AddToFile(TFile* file);
    void FillHist(Charged_Particle& part);
    TH1I* fh_Ntracks           ;
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


private:
    std::string fDir;
    int ftype;
};
extern Hist_dir* Initial_dir;
extern Hist_dir* K3pi_selection;
extern Hist_dir* dir1;
extern Hist_dir* dir2;
#endif
