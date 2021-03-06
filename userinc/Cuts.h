#ifndef __Cuts_h_
#define __Cuts_h_
#include "Charged_Particle.h"
#include "MC_Charged_Particle.h"
#include "Hist_dir.h"
#include "user_NEW.h"
#include "Cuts.h"


class Cuts{
public:
    /* Cuts(superCmpEvent* sevt, superBurst *sburst); */
    Cuts();
    ~Cuts();
    double DCH_e1e2;
    double DCH_mue1;
    double DCH_mue2;
    double DCH_Radius_mu;
    double DCH_Radius_el1;
    double DCH_Radius_el2;
    double DCH_Radius_pi1;
    double DCH_Radius_pi2;
    double DCH_Radius_pi3;

    double Lkr_x_el1 ;
    double Lkr_x_el2 ;
    double Lkr_x_mu  ;
    double Lkr_y_el1 ;
    double Lkr_y_el2 ;
    double Lkr_y_mu  ;
    double Lkr_cut_el1;
    double Lkr_cut_el2;
    double Lkr_cut_mu;

    double MUV_x_el1 ;
    double MUV_x_el2 ;
    double MUV_x_mu  ;
    double MUV_y_el1 ;
    double MUV_y_el2 ;
    double MUV_y_mu  ;


    double MUV_x_pi1 ;
    double MUV_x_pi2 ;
    double MUV_x_pi3 ;
    double MUV_y_pi1 ;
    double MUV_y_pi2 ;
    double MUV_y_pi3 ;

    double Lkr_x_pi1 ;
    double Lkr_x_pi2 ;
    double Lkr_x_pi3 ;
    double Lkr_y_pi1 ;
    double Lkr_y_pi2 ;
    double Lkr_y_pi3  ;
    double Lkr_cut_pi1;
    double Lkr_cut_pi2;
    double Lkr_cut_pi3;


    double Hod_e1e2;
    double cluster_e1e2;
    double cluster_mue1;
    double cluster_mue2;
    double Hod_mue1;
    double Hod_mue2;
    double Mu_P;
    double E1_P;
    double E2_P;
    double mee;
    double muee_P;
    double muee_M;
    double muee_Pt;
    double zvtx_mue1_mue2;
    double zvtx_mue1_e1e2;
    double zvtx_mue2_e1e2;
    double yvtx_mue1_mue2;
    double yvtx_mue1_e1e2;
    double yvtx_mue2_e1e2;
    double xvtx_mue1_mue2;
    double xvtx_mue1_e1e2;
    double xvtx_mue2_e1e2;

    double zvtx_pi1pi2_pi2pi3;
    double zvtx_pi1pi2_pi1pi3;
    double zvtx_pi1pi3_pi2pi3;
    double yvtx_pi1pi2_pi2pi3;
    double yvtx_pi1pi2_pi1pi3;
    double yvtx_pi1pi3_pi2pi3;
    double xvtx_pi1pi2_pi2pi3;
    double xvtx_pi1pi2_pi1pi3;
    double xvtx_pi1pi3_pi2pi3;
};
#endif
