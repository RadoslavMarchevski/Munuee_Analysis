#include "Hist_dir_pi0d.h"
#include "Charged_Particle.h"
#include "Neutral_Particle.h"
#include "Cuts_pi0d.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH1I.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <string>
Hist_dir::Hist_dir(const std::string& dir_Name, int type){

    ftype = type;
    fDir  = dir_Name;
    mKaon = 0.493677;
    //Creating directories with different histograms depending
    //on the directory type:
    //0 - Initial
    //1 - K3pi selection
    //2 - Kmunuee selection


    fh_Ntracks           = new TH1I("Ntracks","Number of charged tracks;Ntrk;Nevents",10,0,10);
    fh_Nclusters         = new TH1I("Nclusters","Number of clusters in the Lkr;Nclusters;Nevents",10,0,10);
    fh_Nvtx              = new TH1I("Nvtx","Number of vertexes;Nvtx;Nevents",10,0,10);
    fh_Kaon_Charge       = new TH1I("Kaon_Charge","Charge of the Kaon;Q;Nevents",4,-2,2);
    fh_el1_Charge        = new TH1I("el1_Charge","Charge of the electron 1;Q;Nevents",4,-2,2);
    fh_el2_Charge        = new TH1I("el2_Charge","Charge of the electron 2;Q;Nevents",4,-2,2);
    fh_pion_Charge       = new TH1I("pion_Charge","Charge of the #pi;Q;Nevents",4,-2,2);
    fh_gamma_momentum    = new TH1F("gamma_momentum","#gamma momentum ;#gamma_P[GeV];Nevents ",100.,0.,100.);
    fh_pi_momentum       = new TH1F("pi_momentum","#gamma momentum ;#pi_P[GeV];Nevents ",100.,0.,100.);
    fh_el1_momentum      = new TH1F("el1_momentum","el1 momentum ;el1_P[GeV];Nevents ",100.,0.,100.);
    fh_el2_momentum      = new TH1F("el2_momentum","el2 momentum ;el2_P[GeV];Nevents ",100.,0.,100.);
    fh_k2pi0d_P          = new TH1F ("k2pi0d_P", "k2pi0d momentum ;P_{#pi e e #gamma}[GeV];Nevents",100,0.,100);
    fh_k2pi0d_Pt         = new TH1F("k2pi0d_Pt","k2pi0d system momentum ;Pt_{#pi e e #gamma}[GeV];Nevents",200,0,2.);
    fh_k2pi0d_M          = new TH1F ("k2pi0d_M", "k2pi0d invariant mass ;M_{#pi e e #gamma}[GeV];Nevents",200,0.4,0.6);

    fh_pi0d_P            = new TH1F ("pi0d_P", "pi0d momentum ;P_{#pi e e #gamma}[GeV];Nevents",100,0.,100);
    fh_pi0d_Pt           = new TH1F("pi0d_Pt","pi0d system momentum ;Pt_{#pi e e #gamma}[GeV];Nevents",200,0,2.);
    fh_pi0d_M            = new TH1F ("pi0d_M", "pi0d invariant mass ;M_{e e #gamma}[GeV];Nevents",100,0.1,0.2);
    fh_EoP_el1           = new TH1F("EoP_el1","Electron 1 E/p ;E/p;Nevents",120.,0.,1.2);
    fh_EoP_el2           = new TH1F("EoP_el2","Electron 2 E/p ;E/p;Nevents",120.,0.,1.2);
    fh_EoP_pion          = new TH1F("EoP_pion","#pion E/p ;E/p;Nevents",120.,0.,1.2);

    fh_lda3_e1_plus      = new TH1F("lda3_e1_plus","lda3 variable for electron 1 +", 150,0.,1.5);
    fh_lda3_e1_minus     = new TH1F("lda3_e2_minus","lda3 variable for electron 1 -", 150,0.,1.5);
    fh_lda3_e2           = new TH1F("lda3_e2","lda3 variable for electron 2", 200,0.,2);
    fh_el1_plus_Charge   = new TH1I("el1_plus_Charge","Charge of the electron 1 +;Q;Nevents",4,-2,2);
    fh_EoP_el1_plus      = new TH1F("EoP_el1_plus","Electron 1 + E/p ;E/p;Nevents",120.,0.,1.2);
    fh_EoP_vs_p_el1_plus = new TH2F("EoP_vs_p_el1_plus","E/p vs Momentum for the electron 1 +",120,0.,1.2, 20,0., 100.);
    fh_lda3_vs_p_e1_plus = new TH2F("lda3_vs_p_e1_plus","lda3 variable vs P for electron 1 +", 150,0.,1.5,20,0.,100.);
    fh_el1_minus_Charge  = new TH1I("el1_minus_Charge","Charge of the electron 1 -;Q;Nevents",4,-2,2);
    fh_EoP_el1_minus     = new TH1F("EoP_el1_minus","Electron 1 - E/p ;E/p;Nevents",120.,0.,1.2);
    fh_EoP_vs_p_el1_minus= new TH2F("EoP_vs_p_el1_minus","E/p vs Momentum for the electron 1 -",120,0.,1.2, 20,0., 100.);
    fh_lda3_vs_p_e1_minus= new TH2F("lda3_vs_p_e1_minus","lda3 variable vs P for electron 1 -", 150,0.,1.5,20,0.,100.);
    fh_DCH_timediff_pi_e1= new TH1F("DCH_timediff_pi_e1","DCHtime difference between #pi and e1;DCHtimediff[ns];Nevents",100,-50.,50.);
    fh_DCH_timediff_pi_e2= new TH1F("DCH_timediff_pi_e2","DCHtime difference between #pi and e2;DCHtimediff[ns];Nevents",100,-50.,50.);
    fh_DCH_timediff_e1_e2= new TH1F("DCH_timediff_e1_e2","DCHtime difference between e1 and e2;DCHtimediff[ns];Nevents",100,-50.,50.);
    fh_Hod_timediff_pi_e1= new TH1F("Hod_timediff_pi_e1","Hodoscope time difference between #pi and e1;Hodtimediff[ns];Nevents",100,-50.,50.);
    fh_Hod_timediff_pi_e2= new TH1F("Hod_timediff_pi_e2","Hodoscope time difference between #pi and e2;Hodtimediff[ns];Nevents",100,-50.,50.);
    fh_Hod_timediff_e1_e2= new TH1F("Hod_timediff_e1_e2","Hodoscope time difference between e1 and e2 ;Hodtimediff[ns];Nevents",100,-50.,50.);
    fh_el1_el2_cltimediff= new TH1F("el1_el2_cltimediff","Time difference between the clustes of the two electrons;Tdiff_{e1e2}[ns];Nevents",200,-50.,50.);
    fh_el1_g_cltimediff  = new TH1F("el1_g_cltimediff","Time difference between the clustes of the electron 1 and the #gamma;Tdiff_{e1 #gamma}[ns];Nevents",200,-50.,50.);
    fh_el2_g_cltimediff  = new TH1F("el2_g_cltimediff","Time difference between the clustes of the electron 2 and the #gamma;Tdiff_{#gamma e2}[ns];Nevents",200,-50.,50.);
    fh_DCH1_distance_e1e2= new TH1F("DCH1_distance_e1e2","Distance between the two electrons in the DCH1",200,0,200.);
    fh_Lkr_distance_e1e2 = new TH1F("Lkr_distance_e1e2","Distance between the two electrons in the Lkr",200,0,200.);
    fh_Lkr_extrap_x_vs_y = new TH2F("Lkr_extrap_x_vs_y","Extrapolated tracks to the Lkr ;Lkr_x[cm];Lkr_y[cm]",400,-200., 200., 400, -200., 200.);

    //Reconstructed Vtx difference
    fh_xvtxdiff_pie1_pie2= new TH1F("xvtxdiff_pie1_pie2","Difference between #pi e1 and #pi e2 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
    fh_xvtxdiff_pie1_e1e2= new TH1F("xvtxdiff_pie1_e1e2","Difference between #pi e1 and e1 e2 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
    fh_xvtxdiff_pie2_e1e2= new TH1F("xvtxdiff_pie2_e1e2","Difference between #pi e2 and e1 e2 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
    fh_yvtxdiff_pie1_pie2= new TH1F("yvtxdiff_pie1_pie2","Difference between #pi e1 and #pi e2 reconstructed Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
    fh_yvtxdiff_pie1_e1e2= new TH1F("yvtxdiff_pie1_e1e2","Difference between #pi e1 and e1 e2 reconstructed  Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
    fh_yvtxdiff_pie2_e1e2= new TH1F("yvtxdiff_pie2_e1e2","Difference between #pi e2 and e1 e2 reconstructed  Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
    fh_zvtxdiff_pie1_pie2= new TH1F("zvtxdiff_pie1_pie2","Difference between #pi e1 and #pi e2 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
    fh_zvtxdiff_pie1_e1e2= new TH1F("zvtxdiff_pie1_e1e2","Difference between #pi e1 and e1 e2 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
    fh_zvtxdiff_pie2_e1e2= new TH1F("zvtxdiff_pie2_e1e2","Difference between #pi e2 and e1 e2 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
    fh_cda_pi_e1         = new TH1F("cda_pi_e1","Closest distance approach between the #pi and e1;cda[cm];Nevents",50, 0., 50.);
    fh_cda_pi_e2         = new TH1F("cda_pi_e2","Closest distance approach between the #pi and e2;cda[cm];Nevents",50, 0., 50.);
    fh_cda_e1_e2         = new TH1F("cda_e1_e2","Closest distance approach between the e1 and e2;cda[cm];Nevents" ,50, 0., 50.);
}

void Hist_dir::AddToFile(TFile* file){
    fh_Ntracks           ->Write();
    fh_Nclusters         ->Write();
    fh_Nvtx              ->Write();
    fh_Kaon_Charge       ->Write();
    fh_el1_Charge        ->Write();
    fh_el2_Charge        ->Write();
    fh_pion_Charge       ->Write();
    fh_EoP_el1           ->Write();
    fh_EoP_el2           ->Write();
    fh_EoP_pion          ->Write();
    fh_gamma_momentum    ->Write();
    fh_pi_momentum       ->Write();
    fh_el1_momentum      ->Write();
    fh_el2_momentum      ->Write();
    fh_k2pi0d_P          ->Write();
    fh_k2pi0d_Pt         ->Write();
    fh_k2pi0d_M          ->Write();
    fh_pi0d_P            ->Write();
    fh_pi0d_Pt           ->Write();
    fh_pi0d_M            ->Write();
    fh_lda3_e1_plus      ->Write();
    fh_lda3_e1_minus     ->Write();
    fh_lda3_e2           ->Write();
    fh_el1_plus_Charge   ->Write();
    fh_EoP_el1_plus      ->Write();
    fh_EoP_vs_p_el1_plus ->Write();
    fh_lda3_vs_p_e1_plus ->Write();
    fh_el1_minus_Charge  ->Write();
    fh_EoP_el1_minus     ->Write();
    fh_EoP_vs_p_el1_minus->Write();
    fh_lda3_vs_p_e1_minus->Write();
    fh_DCH_timediff_pi_e1->Write();
    fh_DCH_timediff_pi_e2->Write();
    fh_DCH_timediff_e1_e2->Write();
    fh_Hod_timediff_pi_e1->Write();
    fh_Hod_timediff_pi_e2->Write();
    fh_Hod_timediff_e1_e2->Write();
    fh_el1_el2_cltimediff->Write();
    fh_el1_g_cltimediff  ->Write();
    fh_el2_g_cltimediff  ->Write();
    fh_DCH1_distance_e1e2->Write();
    fh_Lkr_distance_e1e2 ->Write();
    fh_Lkr_extrap_x_vs_y ->Write();
    fh_xvtxdiff_pie1_pie2->Write();
    fh_xvtxdiff_pie1_e1e2->Write();
    fh_xvtxdiff_pie2_e1e2->Write();
    fh_yvtxdiff_pie1_pie2->Write();
    fh_yvtxdiff_pie1_e1e2->Write();
    fh_yvtxdiff_pie2_e1e2->Write();
    fh_zvtxdiff_pie1_pie2->Write();
    fh_zvtxdiff_pie1_e1e2->Write();
    fh_zvtxdiff_pie2_e1e2->Write();
    fh_cda_pi_e1         ->Write();
    fh_cda_pi_e2         ->Write();
    fh_cda_e1_e2         ->Write();

}

Hist_dir::~Hist_dir(){
    delete fh_Ntracks          ;
    delete fh_Nclusters        ;
    delete fh_Nvtx             ;
    delete fh_Kaon_Charge      ;
    delete fh_el1_Charge       ;
    delete fh_el2_Charge       ;
    delete fh_pion_Charge      ;
    delete fh_gamma_momentum   ;
    delete fh_pi_momentum      ;
    delete fh_el1_momentum     ;
    delete fh_el2_momentum     ;
    delete fh_k2pi0d_P         ;
    delete fh_k2pi0d_Pt        ;
    delete fh_k2pi0d_M         ;
    delete fh_pi0d_P           ;
    delete fh_pi0d_Pt          ;
    delete fh_pi0d_M           ;
    delete fh_EoP_el1 ;
    delete fh_EoP_el2 ;
    delete fh_EoP_pion;
    delete fh_lda3_e1_plus      ;
    delete fh_lda3_e1_minus     ;
    delete fh_lda3_e2           ;
    delete fh_el1_plus_Charge   ;
    delete fh_EoP_el1_plus      ;
    delete fh_EoP_vs_p_el1_plus ;
    delete fh_lda3_vs_p_e1_plus ;
    delete fh_el1_minus_Charge  ;
    delete fh_EoP_el1_minus     ;
    delete fh_EoP_vs_p_el1_minus;
    delete fh_lda3_vs_p_e1_minus;
    delete fh_DCH_timediff_pi_e1;
    delete fh_DCH_timediff_pi_e2;
    delete fh_DCH_timediff_e1_e2;
    delete fh_Hod_timediff_pi_e1;
    delete fh_Hod_timediff_pi_e2;
    delete fh_Hod_timediff_e1_e2;
    delete fh_el1_el2_cltimediff;
    delete fh_el1_g_cltimediff  ;
    delete fh_el2_g_cltimediff  ;
    delete fh_DCH1_distance_e1e2;
    delete fh_Lkr_distance_e1e2 ;
    delete fh_Lkr_extrap_x_vs_y ;
    delete fh_xvtxdiff_pie1_pie2;
    delete fh_xvtxdiff_pie1_e1e2;
    delete fh_xvtxdiff_pie2_e1e2;
    delete fh_yvtxdiff_pie1_pie2;
    delete fh_yvtxdiff_pie1_e1e2;
    delete fh_yvtxdiff_pie2_e1e2;
    delete fh_zvtxdiff_pie1_pie2;
    delete fh_zvtxdiff_pie1_e1e2;
    delete fh_zvtxdiff_pie2_e1e2;
    delete fh_cda_pi_e1         ;
    delete fh_cda_pi_e2         ;
    delete fh_cda_e1_e2         ;







}
