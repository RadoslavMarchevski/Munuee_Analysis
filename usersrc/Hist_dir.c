#include "Hist_dir.h"
#include "Charged_Particle.h"
#include "Neutral_Particle.h"
#include "Cuts.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
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
    //3 - MC final selection


    if(ftype==0){



        fh_Nrun              = new TH1I("Nrun","Run Number;Run Number",2000,15000,17000);
        fh_Ntracks           = new TH1I("Ntracks","Number of charged tracks;Ntrk;Nevents",10,0,10);
        fh_Nclusters         = new TH1I("Nclusters","Number of clusters in the Lkr;Nclusters;Nevents",10,0,10);
        fh_Nvtx              = new TH1I("Nvtx","Number of vertexes;Nvtx;Nevents",10,0,10);
        fh_Track_Momentum    = new TH1F("Track_Momentum","Track momentum after Nvtx,Ntrack and ZVtx cut;Track_P[GeV];Nevents ",100.,0.,100.);
        fh_eop               = new TH1F("Track_E/p","Track E/p after Nvtx,Ntrack and ZVtx cut;E/p;Nevents",120.,0.,1.2);
        fh_COmPaCt_Z_Vertex  = new TH1F("COmPaCt_Z_Vertex","Three track vertex from COmPaCt",150,-4000.,11000.);



//MC true histograms
        fh_mc_KDzvtx = new TH1F("mc_KDzvtx","MC decayed kaon Z vertex",150,-4000.,11000.);
        fh_mc_P1_Pzvtx = new TH1F("mc_P1_Pzvtx","MC Particle 1 production Z vertex",150,-4000.,11000.);
        fh_mc_P1_mass = new TH1F("mc_P1_mass","MC Particle 1 Mass",550,0.,0.55);
        fh_mc_P1_momentum = new TH1F("mc_P1_momentum","MC Particle 1 Momentum",100,0.,100.);
        fh_mc_P1_dist_prod_dec = new TH1F("mc_P1_dist_prod_dec","MC Particle1 distance between production and decay vertex; P1Zvtx_{Production - Decay};Nevents",-150,-4000,11000);
        fh_mc_P2_Pzvtx = new TH1F("mc_P2_Pzvtx","MC Particle 2 production Z vertex",150,-4000.,11000.);
        fh_mc_P2_mass = new TH1F("mc_P2_mass","MC Particle 2 Mass",1100,0.,0.55);
        fh_mc_P2_momentum = new TH1F("mc_P2_momentum","MC Particle 2 Momentum",100,0.,100.);
        fh_mc_P2_dist_prod_dec = new TH1F("mc_P2_dist_prod_dec","MC Particle2 distance between production and decay vertex; P2Zvtx_{Production - Decay};Nevents",-150,-4000,11000);
        fh_mc_P3_Pzvtx = new TH1F("mc_P3_Pzvtx","MC Particle 3 production Z vertex",150,-4000.,11000.);
        fh_mc_P3_mass = new TH1F("mc_P3_mass","MC Particle 3 Mass",1100,0.,0.55);
        fh_mc_P3_momentum = new TH1F("mc_P3_momentum","MC Particle 3 Momentum",100,0.,100.);
        fh_mc_P3_dist_prod_dec = new TH1F("mc_P3_dist_prod_dec","MC Particle3 distance between production and decay vertex; P3Zvtx_{Production - Decay};Nevents",-150,-4000,11000);
        fh_mc_P4_Pzvtx = new TH1F("mc_P4_Pzvtx","MC Particle 4 production Z vertex",150,-4000.,11000.);
        fh_mc_three_track_123_momentum = new TH1F("mc_three_track_123_momentum","MC Momentum of the track 1+2+3",100,0.,100.);
        fh_mc_three_track_123_mass = new TH1F("mc_three_track_123_mass","MC Mass of the track 1+2+3 ",50,0,0.5);
        fh_mc_two_track_23_momentum = new TH1F("mc_two_track_23_momentum","MC Momentum of the track 2+3",100,0.,100.);
        fh_mc_two_track_23_mass = new TH1F("mc_two_track_23_mass","MC Mass of the track 2+3 ",50,0,0.5);
        fh_mc_two_track_23_mass_z_variable = new TH1F("mc_two_track_23_mass_z_variable","MC Mass of the track 2+3 in terms of the z variable;z;Nevents",25,0,0.5);
        fh_mc_P4_mass = new TH1F("mc_P4_mass","MC Particle 4 Mass",550,0.,0.55);
        fh_mc_P4_momentum = new TH1F("mc_P4_momentum","MC Particle 4 Momentum",100,0.,100.);
        fh_mc_four_track_1234_momentum = new TH1F("mc_four_track_1234_momentum","MC Momentum of the track 1+2+3+4",100,0.,100.);
        fh_mc_four_track_1234_mass = new TH1F("mc_four_track_1234_mass","MC Mass of the track 1+2+3+4 ",700,0,0.7);

    }
    if(ftype==1){

//Trigger efficiency plots for the different supersamples
        //SS0
        fh_Trigger_bits          = new TH1I("Trigger_bits"    ,"Trigger word ;Trigger word bits",0 ,0, 20);
        fh_SS0_CPRE              = new TH1I("SS0_h_CPRE"      ," SS0 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS0_MB_1TRK_P         = new TH1I("SS0_MB_1TRK_P"   ," SS0 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS0_MB_1VTX           = new TH1I("SS0_MB_1VTX"     ," SS0 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS0_full_trig         = new TH1I("SS0_MB_full_trig"," SS0 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS1_CPRE              = new TH1I("SS1_h_CPRE"      ," SS1 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS1_MB_1TRK_P         = new TH1I("SS1_MB_1TRK_P"   ," SS1 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS1_MB_1VTX           = new TH1I("SS1_MB_1VTX"     ," SS1 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS1_full_trig         = new TH1I("SS1_MB_full_trig"," SS1 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS2_CPRE              = new TH1I("SS2_h_CPRE"      ," SS2 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS2_MB_1TRK_P         = new TH1I("SS2_MB_1TRK_P"   ," SS2 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS2_MB_1VTX           = new TH1I("SS2_MB_1VTX"     ," SS2 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS2_full_trig         = new TH1I("SS2_MB_full_trig"," SS2 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS3_CPRE              = new TH1I("SS3_h_CPRE"      ," SS3 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS3_MB_1TRK_P         = new TH1I("SS3_MB_1TRK_P"   ," SS3 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS3_MB_1VTX           = new TH1I("SS3_MB_1VTX"     ," SS3 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS3_full_trig         = new TH1I("SS3_MB_full_trig"," SS3 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS4_CPRE              = new TH1I("SS4_h_CPRE"      ," SS4 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS4_MB_1TRK_P         = new TH1I("SS4_MB_1TRK_P"   ," SS4 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS4_MB_1VTX           = new TH1I("SS4_MB_1VTX"     ," SS4 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS4_full_trig         = new TH1I("SS4_MB_full_trig"," SS4 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS5_CPRE              = new TH1I("SS5_h_CPRE"      ," SS5 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS5_MB_1TRK_P         = new TH1I("SS5_MB_1TRK_P"   ," SS5 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS5_MB_1VTX           = new TH1I("SS5_MB_1VTX"     ," SS5 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS5_full_trig         = new TH1I("SS5_MB_full_trig"," SS5 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS6_CPRE              = new TH1I("SS6_h_CPRE"      ," SS6 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS6_MB_1TRK_P         = new TH1I("SS6_MB_1TRK_P"   ," SS6 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS6_MB_1VTX           = new TH1I("SS6_MB_1VTX"     ," SS6 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS6_full_trig         = new TH1I("SS6_MB_full_trig"," SS6 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS7_CPRE              = new TH1I("SS7_h_CPRE"      ," SS7 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS7_MB_1TRK_P         = new TH1I("SS7_MB_1TRK_P"   ," SS7 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS7_MB_1VTX           = new TH1I("SS7_MB_1VTX"     ," SS7 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS7_full_trig         = new TH1I("SS7_MB_full_trig"," SS7 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS8_CPRE              = new TH1I("SS8_h_CPRE"      ," SS8 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS8_MB_1TRK_P         = new TH1I("SS8_MB_1TRK_P"   ," SS8 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS8_MB_1VTX           = new TH1I("SS8_MB_1VTX"     ," SS8 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS8_full_trig         = new TH1I("SS8_MB_full_trig"," SS8 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_Ntracks           = new TH1I("Ntracks","Number of charged tracks;Ntrk;Nevents",10,0,10);
        fh_Nclusters         = new TH1I("Nclusters","Number of clusters in the Lkr;Nclusters;Nevents",10,0,10);
        fh_Nvtx              = new TH1I("Nvtx","Number of vertexes;Nvtx;Nevents",10,0,10);
        fh_Kaon_Charge       = new TH1I("Kaon_Charge","Charge of the Kaon;Q;Nevents",4,-2,2);
        fh_Event_Type        = new TH1I("Event_Type","Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_Pion_Momentum     = new TH1F("Pion_Momentum","Pion Momentum for the K3pi selection;#pi_P[GeV];Nevents",100,0.,100.);
        fh_eop               = new TH1F("Track_E/p","Track E/p ;E/p;Nevents",120.,0.,1.2);
        fh_odd_eop           = new TH1F("odd_E/p","E/p for the odd track in the k3pi selection",120,0.,1.2);
        fh_EoP_vs_p_odd_tr = new TH2F("EoP_vs_p_odd_tr","E/p vs Momentum for the odd track",120,0.,1.2, 20,0., 100.);
        fh_lda3_p1 = new TH1F("lda3_p1","lda3 variable for pion 1", 200,0.,2);
        fh_lda3_p2 = new TH1F("lda3_p2","lda3 variable for pion 2", 200,0.,2);
        fh_lda3_p3 = new TH1F("lda3_p3","lda3 variable for pion 3", 200,0.,2);

        fh_COmPaCt_Z_Vertex  = new TH1F("COmPaCt_Z_Vertex","Three track vertex from COmPaCt",150,-4000.,11000.);
        //Reconstructed Vtx difference
        fh_xvtxdiff_pi1pi2_pi2pi3 = new TH1F("xvtxdiff_pi1pi2_pi2pi3","Difference between #pi_1 #pi_2 and #pi_2 #pi_3 reconstructed X vtx; X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_xvtxdiff_pi1pi2_pi1pi3 = new TH1F("xvtxdiff_pi1pi2_pi1pi3","Difference between  #pi_1 #pi_2 and #pi_1 #pi_3 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_xvtxdiff_pi1pi3_pi2pi3 = new TH1F("xvtxdiff_pi1pi3_pi2pi3","Difference between  #pi_1 #pi_3 and #pi_2 #pi_3 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_pi1pi2_pi2pi3 = new TH1F("yvtxdiff_pi1pi2_pi2pi3","Difference between #pi_1 #pi_2 and #pi_2 #pi_3 reconstructed  Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_pi1pi2_pi1pi3 = new TH1F("yvtxdiff_pi1pi2_pi1pi3","Difference between  #pi_1 #pi_2 and #pi_1 #pi_3 reconstructed Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_pi1pi3_pi2pi3 = new TH1F("yvtxdiff_pi1pi3_pi2pi3","Difference between  #pi_1 #pi_3 and #pi_2 #pi_3 reconstructed Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_zvtxdiff_pi1pi2_pi2pi3 = new TH1F("zvtxdiff_pi1pi2_pi2pi3","Difference between #pi_1 #pi_2 and #pi_2 #pi_3 reconstructed  Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_zvtxdiff_pi1pi2_pi1pi3 = new TH1F("zvtxdiff_pi1pi2_pi1pi3","Difference between  #pi_1 #pi_2 and #pi_1 #pi_3 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_zvtxdiff_pi1pi3_pi2pi3 = new TH1F("zvtxdiff_pi1pi3_pi2pi3","Difference between  #pi_1 #pi_3 and #pi_2 #pi_3 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_cda_pi1_pi2 = new TH1F("cda_pi1_pi2","Closest distance approach between the #pi_1 and #pi_3;cda[cm];Nevents",50, 0., 50.);
        fh_cda_pi1_pi3 = new TH1F("cda_pi1_pi3","Closest distance approach between the #pi_1 and #pi_3;cda[cm];Nevents",50, 0., 50.);
        fh_cda_pi2_pi3 = new TH1F("cda_pi2_pi3","Closest distance approach between the #pi_2 and #pi_3;cda[cm];Nevents" ,50, 0., 50.);
        //ENDOF Vertex Reconstruction

        fh_muee_P = new TH1F("muee_P","Three track momentum ;P_{3#pi}[GeV];Nevents",100,0.,100);
        fh_muee_Pt = new TH1F("muee_Pt","Transverse momentum of 3#pi ;Pt_{#mu e e}[GeV];Nevents",200,0,2.);
        fh_muee_M = new TH1F("muee_M","Three track invariant mass ;M_{3#pi}[GeV];Nevents",1000,-1.,1.);

        fh_missing_mass = new TH1F("missing_mass","Missing mass squared;M^{2}_{miss};Nevents",100,-0.05,0.05);

        //Pions
        fh_pion1_bx = new TH1F("pion1_bx","#pi_1 x coordinate on DCH1 before magnet;DCH1_x[cm];Nevents",300,-150.,150.);
        fh_pion1_by = new TH1F("pion1_by","#pi_1 y coordinate on DCH1 before magnet;DCH1_y[cm];Nevents",300,-150.,150.);
        fh_pion2_bx = new TH1F("pion2_bx","#pi_2 x coordinate on DCH1 before magnet;DCH1_x[cm];Nevents",300,-150.,150.);
        fh_pion2_by = new TH1F("pion2_by","#pi_2 y coordinate on DCH1 before magnet;DCH1_y[cm];Nevents",300,-150.,150.);
        fh_pion3_bx = new TH1F("pion3_bx","#pi_3 x coordinate on DCH1 before magnet;DCH1_x[cm];Nevents",300,-150.,150.);
        fh_pion3_by = new TH1F("pion3_by","#pi_3 y coordinate on DCH1 before magnet;DCH1_y[cm];Nevents",300,-150.,150.);
        fh_bx_vs_by_pi1 = new TH2F("bx_vs_by_pi1","#pi_1 x vs y coordinate on DCH1 before the magnet;DCH_x[cm];DCH_y[cm]",300,-150.,150.,300,-150,150);
        fh_bx_vs_by_pi2 = new TH2F("bx_vs_by_pi2","#pi_2 x vs y coordinate on DCH1 before the magnet;DCH_x[cm];DCH_y[cm]",300,-150.,150.,300,-150,150);
        fh_bx_vs_by_pi3 = new TH2F("bx_vs_by_pi3","#pi_3 x vs y coordinate on DCH1 before the magnet;DCH_x[cm];DCH_y[cm]",300,-150.,150.,300,-150,150);
        fh_Lkr_extrap_tracks_x_vs_y = new TH2F("Lkr_extrap_tracks_x_vs_y","Extrapolated tracks to the Lkr ;Lkr_x[cm];Lkr_y[cm]",400,-200., 200., 400, -200., 200.);
        fh_muv_x_y_position = new TH2F("muv_x_y_position","MUV2 extrapolated position x vs y;Track_x[cm];Track_y[cm]",400,-200.,200.,400,-200.,200.);
        fh_MM2_plus = new TH1F("MM2_plus", "Missing mass squared for K+",100, -0.05, 0.05);
        fh_MM2_minus = new TH1F("MM2_plus", "Missing mass squared for K-",100, -0.05, 0.05);
    }
    if(ftype==2){
        fh_Ntracks           = new TH1I("Ntracks","Number of charged tracks;Ntrk;Nevents",10,0,10);
        fh_Nclusters         = new TH1I("Nclusters","Number of clusters in the Lkr;Nclusters;Nevents",10,0,10);
        fh_Nvtx              = new TH1I("Nvtx","Number of vertexes;Nvtx;Nevents",10,0,10);
        fh_Kaon_Charge       = new TH1I("Kaon_Charge","Charge of the Kaon;Q;Nevents",4,-2,2);
        fh_Mu_charge         = new TH1I("Mu_charge","Muon charge ;Mu_Q;Nevents",4.,-2.,2.);
        fh_Event_Type        = new TH1I("Event_Type","Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_Track_Momentum    = new TH1F("Track_Momentum","Track momentum ;Track_P[GeV];Nevents ",100.,0.,100.);
        fh_Mu_momentum       = new TH1F("Mu_momentum","Muon momentum ;Mu_P[GeV];Nevents ",100.,0.,100.);
        fh_Electron_Momentum = new TH1F("Electron_Momentum", "Electron momentum ;El_P[GeV];Nevents",100,0.,100.);
        fh_eop               = new TH1F("Track_E/p","Track E/p ;E/p;Nevents",120.,0.,1.2);
        fh_Mu_eop            = new TH1F("Mu_E/p","Muon E/p ;E/p;Nevents",120.,0.,1.2);
        fh_Electron_eop      = new TH1F("Electron_E/p", "Electron E/p ;E/p;Nevents",120,0.,1.2);
        fh_odd_eop           = new TH1F("odd_E/p","E/p for the odd track in the k3pi selection",120,0.,1.2);
        fh_COmPaCt_Z_Vertex  = new TH1F("COmPaCt_Z_Vertex","Three track vertex from COmPaCt",150,-4000.,11000.);
        fh_Mu_Zvtx_min_COmPaCt_Zvtx = new TH1F("Mu_Zvtx_min_COmPaCt_Zvtx","Calculated #mu Zvtx - COmPact three track Zvtx;(Z_#mu - Z_vtx)[cm];Nevents ",100,-10000.,10000.);
        //Reconstructed Vtx difference
        fh_xvtxdiff_mue1_mue2 = new TH1F("xvtxdiff_mue1_mue2","Difference between #mu e1 and #mu e2 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_xvtxdiff_mue1_e1e2 = new TH1F("xvtxdiff_mue1_e1e2","Difference between #mu e1 and e1 e2 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_xvtxdiff_mue2_e1e2 = new TH1F("xvtxdiff_mue2_e1e2","Difference between #mu e2 and e1 e2 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_mue1_mue2 = new TH1F("yvtxdiff_mue1_mue2","Difference between #mu e1 and #mu e2 reconstructed Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_mue1_e1e2 = new TH1F("yvtxdiff_mue1_e1e2","Difference between #mu e1 and e1 e2 reconstructed  Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_mue2_e1e2 = new TH1F("yvtxdiff_mue2_e1e2","Difference between #mu e2 and e1 e2 reconstructed  Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_zvtxdiff_mue1_mue2 = new TH1F("zvtxdiff_mue1_mue2","Difference between #mu e1 and #mu e2 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_zvtxdiff_mue1_e1e2 = new TH1F("zvtxdiff_mue1_e1e2","Difference between #mu e1 and e1 e2 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_zvtxdiff_mue2_e1e2 = new TH1F("zvtxdiff_mue2_e1e2","Difference between #mu e2 and e1 e2 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);

        //ENDOF Reconstructed Vtx difference

        //Time matching
        fh_DCHtime_mu        = new TH1F("DCHtime_mu","Drift Chamber time for the #mu ;DCHtime[ns];Nevents",200,50,250);
        fh_DCHtime_e1        = new TH1F("DCHtime_e1","Drift Chamber time for the e1 ;DCHtime[ns];Nevents",200,50,250);
        fh_DCHtime_e2        = new TH1F("DCHtime_e2","Drift Chamber time for the e2 ;DCHtime[ns];Nevents",200,50,250);
        fh_HodTime_mu        = new TH1F("HodTime_mu","Hodoscope time for the #mu ;HodTime[ns];Nevents",200,50,250);
        fh_HodTime_e1        = new TH1F("HodTime_e1","Hodoscope time for the e1 ;HodTime[ns];Nevents",200,50,250);
        fh_HodTime_e2        = new TH1F("HodTime_e2","Hodoscope time for the e2 ;HodTime[ns];Nevents",200,50,250);
        fh_DCH_timediff_mu_e1 = new TH1F("DCH_timediff_mu_e1","DCHtime difference between #mu and e1;DCHtimediff[ns];Nevents",100,-50.,50.);
        fh_DCH_timediff_mu_e2 = new TH1F("DCH_timediff_mu_e2","DCHtime difference between #mu and e2;DCHtimediff[ns];Nevents",100,-50.,50.);
        fh_DCH_timediff_e1_e2 = new TH1F("DCH_timediff_e1_e2","DCHtime difference between e1 and e2;DCHtimediff[ns];Nevents",100,-50.,50.);
        fh_Hod_timediff_mu_e1 = new TH1F("Hod_timediff_mu_e1","Hodoscope time difference between #mu and e1;Hodtimediff[ns];Nevents",100,-50.,50.);
        fh_Hod_timediff_mu_e2 = new TH1F("Hod_timediff_mu_e2","Hodoscope time difference between #mu and e2;Hodtimediff[ns];Nevents",100,-50.,50.);
        fh_Hod_timediff_e1_e2 = new TH1F("Hod_timediff_e1_e2","Hodoscope time difference between e1 and e2 ;Hodtimediff[ns];Nevents",100,-50.,50.);
        //ENDOF Time matching
        fh_mee = new TH1F("mee","Invariant mass of the electron pair ",50,0.,0.5);
        fh_mee_z_variable = new TH1F("mee_z_variable","Invariant mass of the electron pair in terms of the z variable;z;Nevents",25,0,0.5);
        fh_muee_P = new TH1F("muee_P","Three track momentum ;P_{#mu e e}[GeV];Nevents",100,0.,100);
        fh_muee_Pt = new TH1F("muee_Pt","Transverse momentum of #mu^{#pm} e^{+}e^{-} system;Pt_{#mu e e}[GeV];Nevents",200,0,2.);
        fh_muee_M = new TH1F("muee_M","Three track invariant mass ;M_{#mu e e}[GeV];Nevents",100,0.,1.);
        fh_Muee_M_3pi_assumption = new TH1F("Muee_M_3pi_assumption","Three track invariant mass with 3 #pi assumtion ;M_{3#pi};Nevents",100,0.,1.);

        fh_missing_mass = new TH1F("missing_mass","Missing mass squared;M^{2}_{miss}[GeV^2];Nevents",100,-0.05,0.05);
        fh_cda_mu_e1 = new TH1F("cda_mu_e1","Closest distance approach between the #mu and e1;cda[cm];Nevents",50, 0., 50.);
        fh_cda_mu_e2 = new TH1F("cda_mu_e2","Closest distance approach between the #mu and e2;cda[cm];Nevents",50, 0., 50.);
        fh_cda_e1_e2 = new TH1F("cda_e1_e2","Closest distance approach between the e1 and e2;cda[cm];Nevents" ,50, 0., 50.);
        //Drift chambers geometry
        //Muons
        fh_muon_bx = new TH1F("muon_bx","#mu x coordinate on DCH1 before magnet;DCH1_x[cm];Nevents",300,-150,150);
        fh_muon_by = new TH1F("muon_by","#mu y coordinate on DCH1 before magnet;DCH1_y[cm];Nevents",300,-150,150);
        //Electrons
        fh_electron1_bx = new TH1F("electron1_bx","el1 x coordinate on DCH1 before magnet;DCH1_x[cm];Nevents",300,-150.,150.);
        fh_electron1_by = new TH1F("electron1_by","el1 y coordinate on DCH1 before magnet;DCH1_y[cm];Nevents",300,-150.,150.);
        fh_electron2_bx = new TH1F("electron2_bx","el2 x coordinate on DCH1 before magnet;DCH1_x[cm];Nevents",300,-150.,150.);
        fh_electron2_by = new TH1F("electron2_by","el2 y coordinate on DCH1 before magnet;DCH1_y[cm];Nevents",300,-150.,150.);
        fh_bx_vs_by_muon = new TH2F("bx_vs_by_muon","#mu x vs y coordinate on DCH1 before the magnet;DCH_x[cm];DCH_y[cm]",300,-150.,150.,300,-150,150);
        fh_bx_vs_by_el1 = new TH2F("bx_vs_by_el1","el1 x vs y coordinate on DCH1 before the magnet;DCH_x[cm];DCH_y[cm]",300,-150.,150.,300,-150,150);
        fh_bx_vs_by_el2 = new TH2F("bx_vs_by_el2","el2 x vs y coordinate on DCH1 before the magnet;DCH_x[cm];DCH_y[cm]",300,-150.,150.,300,-150,150);
        fh_el1_cluster_x = new TH1F("el1_cluster_x","x position of el1 cluster;Cl_x[cm];Nevents",400,-200.,200.);
        fh_el1_cluster_y = new TH1F("el1_cluster_y","y position of el1 cluster;Cl_y[cm];Nevents",400,-200.,200.);
        fh_el2_cluster_x = new TH1F("el2_cluster_x","x position of el2 cluster;Cl_x[cm];Nevents",400,-200.,200.);
        fh_el2_cluster_y = new TH1F("el2_cluster_y","y position of el2 cluster;Cl_y[cm];Nevents",400,-200.,200.);
        fh_mu_cluster_x = new TH1F("mu_cluster_x","x position of #mu cluster;Cl_x[cm];Nevents",400,-200.,200.);
        fh_mu_cluster_y = new TH1F("mu_cluster_y","y position of #mu cluster;Cl_y[cm];Nevents",400,-200.,200.);
        fh_el1_cluster_time = new TH1F("el1_cluster_time","Time of the el1 cluster;Cl_time[ns];Nevents",200,0.,200.);
        fh_el2_cluster_time = new TH1F("el2_cluster_time","Time of the el2 cluster;Cl_time[ns];Nevents",200,0.,200.);
        fh_el1_el2_cluster_timediff = new TH1F("el1_el2_cluster_timediff","Time difference between the clustes of the two electrons;Tdiff_{e1e2}[ns];Nevents",200,-50.,50.);
        fh_el1_dtrk_cl = new TH1F("el1_dtrk_cl","Difference between extrapolated track from DCH and to the Lkr and cluster for el1   ;dtrk_cl[cm];Nevents",200,0.,200.);
        fh_el2_dtrk_cl = new TH1F("el2_dtrk_cl","Difference between extrapolated track from DCH and to the Lkr and cluster for el2   ;dtrk_cl[cm];Nevents",200,0.,200.);
        fh_muon_dtrk_cl = new TH1F("muon_dtrk_cl","Difference between extrapolated track from DCH and to the Lkr and cluster for #mu ;dtrk_cl[cm];Nevents",200,0.,200.);
        //fh_el1_cluster_x_y = new TH2F("el1_cluster_x_y","Cluster position x vs y for el1;cl_x[cm];cl_y[cm]",400,-200.,200.,400,-200.,-200.);
        //fh_el2_cluster_x_y = new TH2F("el2_cluster_x_y","Cluster position x vs y for el2;cl_x[cm];cl_y[cm]",400,-200.,200.,400,-200.,-200.);
        //fh_mu_cluster_x_y = new TH2F("mu_cluster_x_y","Cluster position x vs y;cl_x[cm];cl_y[cm]",400,-200.,200.,400,-200.,-200.);
        fh_deadcell_distance = new TH1F("deadcell_distance","Distance to deadcell",150,0.,150.);
        fh_muv2_trk_cl_diff = new TH1F("muv2_trk_cl_diff","Difference between extrapolated track and MUV2 cluster position;Trk_cl_diff[cm];Nevents",200,0.,200.);
        fh_muv_xpos = new TH1F("muv_xpos","Muon X position at MUV2",400,-200.,200.);
        fh_muv_ypos = new TH1F("muv_ypos","Muon Y position at MUV1",400,-200.,200.);
        fh_muv_x_y_position = new TH2F("muv_x_y_position","MUV2 extrapolated position x vs y;Track_x[cm];Track_y[cm]",400,-200.,200.,400,-200.,200.);

        fh_DCH1_distance_e1e2 = new TH1F("DCH1_distance_e1e2","Distance between the two electrons in the DCH1",200,0,200.);
        fh_Lkr_distance_e1e2 = new TH1F("Lkr_distance_e1e2","Distance between the two electrons in the Lkr",200,0,200.);
        fh_Lkr_extrap_tracks_x_vs_y = new TH2F("Lkr_extrap_tracks_x_vs_y","Extrapolated tracks to the Lkr ;Lkr_x[cm];Lkr_y[cm]",400,-200., 200., 400, -200., 200.);
        fh_angle_el_mu = new TH1F("angle_el_mu","Angle between e+e- plane and the #mu",200,0.0,0.2);

        //MC true histograms
        fh_mc_KDzvtx = new TH1F("mc_KDzvtx","MC decayed kaon Z vertex",150,-4000.,11000.);
        fh_mc_P1_Pzvtx = new TH1F("mc_P1_Pzvtx","MC Particle 1 production Z vertex",150,-4000.,11000.);
        fh_mc_P1_mass = new TH1F("mc_P1_mass","MC Particle 1 Mass",550,0.,0.55);
        fh_mc_P1_momentum = new TH1F("mc_P1_momentum","MC Particle 1 Momentum",100,0.,100.);
        fh_mc_P2_Pzvtx = new TH1F("mc_P2_Pzvtx","MC Particle 2 production Z vertex",150,-4000.,11000.);
        fh_mc_P2_mass = new TH1F("mc_P2_mass","MC Particle 2 Mass",1100,0.,0.55);
        fh_mc_P2_momentum = new TH1F("mc_P2_momentum","MC Particle 2 Momentum",100,0.,100.);
        fh_mc_P1_dist_prod_dec = new TH1F("mc_P1_dist_prod_dec","MC Particle1 distance between production and decay vertex; P1Zvtx_{Production - Decay};Nevents",-150,-4000,11000);
        fh_mc_P2_dist_prod_dec = new TH1F("mc_P2_dist_prod_dec","MC Particle2 distance between production and decay vertex; P2Zvtx_{Production - Decay};Nevents",-150,-4000,11000);
        fh_mc_P3_dist_prod_dec = new TH1F("mc_P3_dist_prod_dec","MC Particle3 distance between production and decay vertex; P3Zvtx_{Production - Decay};Nevents",-150,-4000,11000);
        fh_mc_P3_Pzvtx = new TH1F("mc_P3_Pzvtx","MC Particle 3 production Z vertex",150,-4000.,11000.);
        fh_mc_P3_mass = new TH1F("mc_P3_mass","MC Particle 3 Mass",1100,0.,0.55);
        fh_mc_P3_momentum = new TH1F("mc_P3_momentum","MC Particle 3 Momentum",100,0.,100.);
        fh_mc_P4_Pzvtx = new TH1F("mc_P4_Pzvtx","MC Particle 4 production Z vertex",150,-4000.,11000.);
        fh_mc_three_track_123_momentum = new TH1F("mc_three_track_123_momentum","MC Momentum of the track 1+2+3",100,0.,100.);
        fh_mc_three_track_123_mass = new TH1F("mc_three_track_123_mass","MC Mass of the track 1+2+3 ",70,0,0.7);
        fh_mc_two_track_23_momentum = new TH1F("mc_two_track_23_momentum","MC Momentum of the track 2+3",100,0.,100.);
        fh_mc_two_track_23_mass = new TH1F("mc_two_track_23_mass","MC Mass of the track 2+3;M_{ee}[GeV]; ",50,0,0.5);
        fh_mc_two_track_23_mass_z_variable = new TH1F("mc_two_track_23_mass_z_variable","MC Mass of the track 2+3 in terms of the z variable;z;Nevents",25,0,0.5);
        fh_mc_P4_mass = new TH1F("mc_P4_mass","MC Particle 4 Mass",550,0.,0.55);
        fh_mc_P4_momentum = new TH1F("mc_P4_momentum","MC Particle 4 Momentum",100,0.,100.);
        fh_mc_four_track_1234_momentum = new TH1F("mc_four_track_1234_momentum","MC Momentum of the track 1+2+3+4",100,0.,100.);
        fh_mc_four_track_1234_mass = new TH1F("mc_four_track_1234_mass","MC Mass of the track 1+2+3+4 ",70,0,0.7);
        fh_lda3_e1 = new TH1F("lda3_e1","lda3 variable for electron 1", 200,0.,2);
        fh_lda3_e2 = new TH1F("lda3_e2","lda3 variable for electron 2", 200,0.,2);
        fh_Mee_before = new TH1F("z_vtx_before","ComPaCt z vtx before z vtx cut",150,-4000.,11000.);
        fh_muon_status = new TH1I("muon_status","1 - MUV 1 & 2 & 3, 2 - 1 & 2 &! 3, 3 - !1 &2 &3, 4 - 1 & &!2 &3", 10, 0, 10);
        fh_MM2_plus = new TH1F("MM2_plus", "Missing mass squared for K+",100, -0.05, 0.05);
        fh_MM2_minus = new TH1F("MM2_minus", "Missing mass squared for K-",100, -0.05, 0.05);


//Trigger efficiency plots for the different supersamples
        //SS0
        fh_Trigger_bits          = new TH1I("Trigger_bits"    ,"Trigger word ;Trigger word bits",0 ,0, 20);
        fh_SS0_CPRE              = new TH1I("SS0_h_CPRE"      ," SS0 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS0_MB_1TRK_P         = new TH1I("SS0_MB_1TRK_P"   ," SS0 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS0_MB_1VTX           = new TH1I("SS0_MB_1VTX"     ," SS0 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS0_full_trig         = new TH1I("SS0_MB_full_trig"," SS0 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS1_CPRE              = new TH1I("SS1_h_CPRE"      ," SS1 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS1_MB_1TRK_P         = new TH1I("SS1_MB_1TRK_P"   ," SS1 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS1_MB_1VTX           = new TH1I("SS1_MB_1VTX"     ," SS1 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS1_full_trig         = new TH1I("SS1_MB_full_trig"," SS1 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS2_CPRE              = new TH1I("SS2_h_CPRE"      ," SS2 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS2_MB_1TRK_P         = new TH1I("SS2_MB_1TRK_P"   ," SS2 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS2_MB_1VTX           = new TH1I("SS2_MB_1VTX"     ," SS2 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS2_full_trig         = new TH1I("SS2_MB_full_trig"," SS2 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS3_CPRE              = new TH1I("SS3_h_CPRE"      ," SS3 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS3_MB_1TRK_P         = new TH1I("SS3_MB_1TRK_P"   ," SS3 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS3_MB_1VTX           = new TH1I("SS3_MB_1VTX"     ," SS3 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS3_full_trig         = new TH1I("SS3_MB_full_trig"," SS3 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS4_CPRE              = new TH1I("SS4_h_CPRE"      ," SS4 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS4_MB_1TRK_P         = new TH1I("SS4_MB_1TRK_P"   ," SS4 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS4_MB_1VTX           = new TH1I("SS4_MB_1VTX"     ," SS4 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS4_full_trig         = new TH1I("SS4_MB_full_trig"," SS4 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS5_CPRE              = new TH1I("SS5_h_CPRE"      ," SS5 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS5_MB_1TRK_P         = new TH1I("SS5_MB_1TRK_P"   ," SS5 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS5_MB_1VTX           = new TH1I("SS5_MB_1VTX"     ," SS5 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS5_full_trig         = new TH1I("SS5_MB_full_trig"," SS5 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS6_CPRE              = new TH1I("SS6_h_CPRE"      ," SS6 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS6_MB_1TRK_P         = new TH1I("SS6_MB_1TRK_P"   ," SS6 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS6_MB_1VTX           = new TH1I("SS6_MB_1VTX"     ," SS6 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS6_full_trig         = new TH1I("SS6_MB_full_trig"," SS6 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS7_CPRE              = new TH1I("SS7_h_CPRE"      ," SS7 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS7_MB_1TRK_P         = new TH1I("SS7_MB_1TRK_P"   ," SS7 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS7_MB_1VTX           = new TH1I("SS7_MB_1VTX"     ," SS7 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS7_full_trig         = new TH1I("SS7_MB_full_trig"," SS7 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS8_CPRE              = new TH1I("SS8_h_CPRE"      ," SS8 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS8_MB_1TRK_P         = new TH1I("SS8_MB_1TRK_P"   ," SS8 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS8_MB_1VTX           = new TH1I("SS8_MB_1VTX"     ," SS8 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS8_full_trig         = new TH1I("SS8_MB_full_trig"," SS8 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);


    }
    if(ftype==3){
        fh_missing_mass = new TH1F("missing_mass","Missing mass squared;M^{2}_{miss};Nevents",100,-0.05,0.05);
        fh_mee = new TH1F("mee","Invariant mass of the electron pair ",50,0.,0.5);
        fh_mee_z_variable = new TH1F("mee_z_variable","Invariant mass of the electron pair in terms of the z variable;z;Nevents",25,0,0.5);
        fh_MM2_plus = new TH1F("MM2_plus", "Missing mass squared for K+",100, -0.05, 0.05);
        fh_MM2_minus = new TH1F("MM2_minus", "Missing mass squared for K-",100, -0.05, 0.05);

        fh_mee_z_variable->Sumw2();
        fh_mee->Sumw2();
        fh_missing_mass->Sumw2();
    }
}

Hist_dir::Hist_dir(const std::string& dir_Name, int type, const std::string& tag){

    ftype = type;
    fDir  = dir_Name;
    mKaon = 0.493677;
    //Creating directories with different histograms depending
    //on the directory type:
    //0 - Initial
    //1 - K3pi selection
    //2 - Kmunuee selection
    //3 - MC final selection

    //std::string hname ;

    if(ftype==0){



        fh_Nrun              = new TH1I(Form("%s_%s_Nrun",tag.c_str()),"Run Number;Run Number",2000,15000,17000);
        fh_Ntracks           = new TH1I(Form("%s_Ntracks",tag.c_str()),"Number of charged tracks;Ntrk;Nevents",10,0,10);
        fh_Nclusters         = new TH1I(Form("%s_Nclusters",tag.c_str()),"Number of clusters in the Lkr;Nclusters;Nevents",10,0,10);
        fh_Nvtx              = new TH1I(Form("%s_Nvtx",tag.c_str()),"Number of vertexes;Nvtx;Nevents",10,0,10);
        fh_Track_Momentum    = new TH1F(Form("%s_Track_Momentum",tag.c_str()),"Track momentum after Nvtx,Ntrack and ZVtx cut;Track_P[GeV];Nevents ",100.,0.,100.);
        fh_eop               = new TH1F(Form("%s_Track_E/p",tag.c_str()),"Track E/p after Nvtx,Ntrack and ZVtx cut;E/p;Nevents",120.,0.,1.2);
        fh_COmPaCt_Z_Vertex  = new TH1F(Form("%s_COmPaCt_Z_Vertex",tag.c_str()),"Three track vertex from COmPaCt",150,-4000.,11000.);



//MC true histograms
        fh_mc_KDzvtx = new TH1F(Form("%s_mc_KDzvtx",tag.c_str()),"MC decayed kaon Z vertex",150,-4000.,11000.);
        fh_mc_P1_Pzvtx = new TH1F(Form("%s_mc_P1_Pzvtx",tag.c_str()),"MC Particle 1 production Z vertex",150,-4000.,11000.);
        fh_mc_P1_mass = new TH1F(Form("%s_mc_P1_mass",tag.c_str()),"MC Particle 1 Mass",550,0.,0.55);
        fh_mc_P1_momentum = new TH1F(Form("%s_mc_P1_momentum",tag.c_str()),"MC Particle 1 Momentum",100,0.,100.);
        fh_mc_P1_dist_prod_dec = new TH1F(Form("%s_mc_P1_dist_prod_dec",tag.c_str()),"MC Particle1 distance between production and decay vertex; P1Zvtx_{Production - Decay};Nevents",-150,-4000,11000);
        fh_mc_P2_Pzvtx = new TH1F(Form("%s_mc_P2_Pzvtx",tag.c_str()),"MC Particle 2 production Z vertex",150,-4000.,11000.);
        fh_mc_P2_mass = new TH1F(Form("%s_mc_P2_mass",tag.c_str()),"MC Particle 2 Mass",1100,0.,0.55);
        fh_mc_P2_momentum = new TH1F(Form("%s_mc_P2_momentum",tag.c_str()),"MC Particle 2 Momentum",100,0.,100.);
        fh_mc_P2_dist_prod_dec = new TH1F(Form("%s_mc_P2_dist_prod_dec",tag.c_str()),"MC Particle2 distance between production and decay vertex; P2Zvtx_{Production - Decay};Nevents",-150,-4000,11000);
        fh_mc_P3_Pzvtx = new TH1F(Form("%s_mc_P3_Pzvtx",tag.c_str()),"MC Particle 3 production Z vertex",150,-4000.,11000.);
        fh_mc_P3_mass = new TH1F(Form("%s_mc_P3_mass",tag.c_str()),"MC Particle 3 Mass",1100,0.,0.55);
        fh_mc_P3_momentum = new TH1F(Form("%s_mc_P3_momentum",tag.c_str()),"MC Particle 3 Momentum",100,0.,100.);
        fh_mc_P3_dist_prod_dec = new TH1F(Form("%s_mc_P3_dist_prod_dec",tag.c_str()),"MC Particle3 distance between production and decay vertex; P3Zvtx_{Production - Decay};Nevents",-150,-4000,11000);
        fh_mc_P4_Pzvtx = new TH1F(Form("%s_mc_P4_Pzvtx",tag.c_str()),"MC Particle 4 production Z vertex",150,-4000.,11000.);
        fh_mc_three_track_123_momentum = new TH1F(Form("%s_mc_three_track_123_momentum",tag.c_str()),"MC Momentum of the track 1+2+3",100,0.,100.);
        fh_mc_three_track_123_mass = new TH1F(Form("%s_mc_three_track_123_mass",tag.c_str()),"MC Mass of the track 1+2+3 ",50,0,0.5);
        fh_mc_two_track_23_momentum = new TH1F(Form("%s_mc_two_track_23_momentum",tag.c_str()),"MC Momentum of the track 2+3",100,0.,100.);
        fh_mc_two_track_23_mass = new TH1F(Form("%s_mc_two_track_23_mass",tag.c_str()),"MC Mass of the track 2+3 ",50,0,0.5);
        fh_mc_two_track_23_mass_z_variable = new TH1F(Form("%s_mc_two_track_23_mass_z_variable",tag.c_str()),"MC Mass of the track 2+3 in terms of the z variable;z;Nevents",25,0,0.5);
        fh_mc_P4_mass = new TH1F(Form("%s_mc_P4_mass",tag.c_str()),"MC Particle 4 Mass",550,0.,0.55);
        fh_mc_P4_momentum = new TH1F(Form("%s_mc_P4_momentum",tag.c_str()),"MC Particle 4 Momentum",100,0.,100.);
        fh_mc_four_track_1234_momentum = new TH1F(Form("%s_mc_four_track_1234_momentum",tag.c_str()),"MC Momentum of the track 1+2+3+4",100,0.,100.);
        fh_mc_four_track_1234_mass = new TH1F(Form("%s_mc_four_track_1234_mass",tag.c_str()),"MC Mass of the track 1+2+3+4 ",700,0,0.7);

    }
    if(ftype==1){
        //Trigger efficiency plots for the different supersamples
        //SS0
        fh_Trigger_bits          = new TH1I(Form("%s_Trigger_bits" ,tag.c_str())   ,"Trigger word ;Trigger word bits",0 ,0, 20);
        fh_SS0_CPRE              = new TH1I(Form("%s_SS0_h_CPRE"   ,tag.c_str())   ," SS0 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS0_MB_1TRK_P         = new TH1I(Form("%s_SS0_MB_1TRK_P",tag.c_str())   ," SS0 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS0_MB_1VTX           = new TH1I(Form("%s_SS0_MB_1VTX"  ,tag.c_str())   ," SS0 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS0_full_trig         = new TH1I(Form("%s_SS0_MB_full_trig",tag.c_str())," SS0 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS1_CPRE              = new TH1I(Form("%s_SS1_h_CPRE"     ,tag.c_str()) ," SS1 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS1_MB_1TRK_P         = new TH1I(Form("%s_SS1_MB_1TRK_P"  ,tag.c_str()) ," SS1 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS1_MB_1VTX           = new TH1I(Form("%s_SS1_MB_1VTX"    ,tag.c_str()) ," SS1 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS1_full_trig         = new TH1I(Form("%s_SS1_MB_full_trig",tag.c_str())," SS1 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS2_CPRE              = new TH1I(Form("%s_SS2_h_CPRE"     ,tag.c_str()) ," SS2 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS2_MB_1TRK_P         = new TH1I(Form("%s_SS2_MB_1TRK_P"  ,tag.c_str()) ," SS2 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS2_MB_1VTX           = new TH1I(Form("%s_SS2_MB_1VTX"    ,tag.c_str()) ," SS2 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS2_full_trig         = new TH1I(Form("%s_SS2_MB_full_trig",tag.c_str())," SS2 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS3_CPRE              = new TH1I(Form("%s_SS3_h_CPRE"     ,tag.c_str()) ," SS3 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS3_MB_1TRK_P         = new TH1I(Form("%s_SS3_MB_1TRK_P"  ,tag.c_str()) ," SS3 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS3_MB_1VTX           = new TH1I(Form("%s_SS3_MB_1VTX"    ,tag.c_str()) ," SS3 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS3_full_trig         = new TH1I(Form("%s_SS3_MB_full_trig",tag.c_str())," SS3 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS4_CPRE              = new TH1I(Form("%s_SS4_h_CPRE"     ,tag.c_str()) ," SS4 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS4_MB_1TRK_P         = new TH1I(Form("%s_SS4_MB_1TRK_P"  ,tag.c_str()) ," SS4 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS4_MB_1VTX           = new TH1I(Form("%s_SS4_MB_1VTX"    ,tag.c_str()) ," SS4 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS4_full_trig         = new TH1I(Form("%s_SS4_MB_full_trig",tag.c_str())," SS4 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS5_CPRE              = new TH1I(Form("%s_SS5_h_CPRE"     ,tag.c_str()) ," SS5 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS5_MB_1TRK_P         = new TH1I(Form("%s_SS5_MB_1TRK_P"  ,tag.c_str()) ," SS5 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS5_MB_1VTX           = new TH1I(Form("%s_SS5_MB_1VTX"    ,tag.c_str()) ," SS5 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS5_full_trig         = new TH1I(Form("%s_SS5_MB_full_trig",tag.c_str())," SS5 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS6_CPRE              = new TH1I(Form("%s_SS6_h_CPRE"     ,tag.c_str()) ," SS6 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS6_MB_1TRK_P         = new TH1I(Form("%s_SS6_MB_1TRK_P"  ,tag.c_str()) ," SS6 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS6_MB_1VTX           = new TH1I(Form("%s_SS6_MB_1VTX"    ,tag.c_str()) ," SS6 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS6_full_trig         = new TH1I(Form("%s_SS6_MB_full_trig",tag.c_str())," SS6 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS7_CPRE              = new TH1I(Form("%s_SS7_h_CPRE"    ,tag.c_str())  ," SS7 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS7_MB_1TRK_P         = new TH1I(Form("%s_SS7_MB_1TRK_P" ,tag.c_str())  ," SS7 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS7_MB_1VTX           = new TH1I(Form("%s_SS7_MB_1VTX"   ,tag.c_str())  ," SS7 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS7_full_trig         = new TH1I(Form("%s_SS7_MB_full_trig",tag.c_str())," SS7 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_SS8_CPRE              = new TH1I(Form("%s_SS8_h_CPRE"    ,tag.c_str())  ," SS8 CPRE trigger word - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS8_MB_1TRK_P         = new TH1I(Form("%s_SS8_MB_1TRK_P" ,tag.c_str())  ," SS8 MB_1TRK_P trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS8_MB_1VTX           = new TH1I(Form("%s_SS8_MB_1VTX"   ,tag.c_str())  ," SS8 MB_1VTX trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_SS8_full_trig         = new TH1I(Form("%s_SS8_MB_full_trig",tag.c_str())," SS8 MB_full_trigger - Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);

        fh_Ntracks           = new TH1I(Form("%s_Ntracks",tag.c_str()),"Number of charged tracks;Ntrk;Nevents",10,0,10);
        fh_Nclusters         = new TH1I(Form("%s_Nclusters",tag.c_str()),"Number of clusters in the Lkr;Nclusters;Nevents",10,0,10);
        fh_Nvtx              = new TH1I(Form("%s_Nvtx",tag.c_str()),"Number of vertexes;Nvtx;Nevents",10,0,10);
        fh_Kaon_Charge       = new TH1I(Form("%s_Kaon_Charge",tag.c_str()),"Charge of the Kaon;Q;Nevents",4,-2,2);
        fh_Event_Type        = new TH1I(Form("%s_Event_Type",tag.c_str()),"Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_Pion_Momentum     = new TH1F(Form("%s_Pion_Momentum",tag.c_str()),"Pion Momentum for the K3pi selection;#pi_P[GeV];Nevents",100,0.,100.);
        fh_eop               = new TH1F(Form("%s_Track_E/p",tag.c_str()),"Track E/p ;E/p;Nevents",120.,0.,1.2);
        fh_odd_eop           = new TH1F(Form("%s_odd_E/p",tag.c_str()),"E/p for the odd track in the k3pi selection",120,0.,1.2);
        fh_EoP_vs_p_odd_tr = new TH2F(Form("%s_EoP_vs_p_odd_tr",tag.c_str()),"E/p vs Momentum for the odd track",120,0.,1.2, 20,0., 100.);
        fh_lda3_p1 = new TH1F(Form("%s_lda3_p1",tag.c_str()),"lda3 variable for pion 1", 200,0.,2);
        fh_lda3_p2 = new TH1F(Form("%s_lda3_p2",tag.c_str()),"lda3 variable for pion 2", 200,0.,2);
        fh_lda3_p3 = new TH1F(Form("%s_lda3_p3",tag.c_str()),"lda3 variable for pion 3", 200,0.,2);

        fh_COmPaCt_Z_Vertex  = new TH1F(Form("%s_COmPaCt_Z_Vertex",tag.c_str()),"Three track vertex from COmPaCt",150,-4000.,11000.);
        //Reconstructed Vtx difference
        fh_xvtxdiff_pi1pi2_pi2pi3 = new TH1F(Form("%s_xvtxdiff_pi1pi2_pi2pi3",tag.c_str()),"Difference between #pi_1 #pi_2 and #pi_2 #pi_3 reconstructed X vtx; X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_xvtxdiff_pi1pi2_pi1pi3 = new TH1F(Form("%s_xvtxdiff_pi1pi2_pi1pi3",tag.c_str()),"Difference between  #pi_1 #pi_2 and #pi_1 #pi_3 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_xvtxdiff_pi1pi3_pi2pi3 = new TH1F(Form("%s_xvtxdiff_pi1pi3_pi2pi3",tag.c_str()),"Difference between  #pi_1 #pi_3 and #pi_2 #pi_3 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_pi1pi2_pi2pi3 = new TH1F(Form("%s_yvtxdiff_pi1pi2_pi2pi3",tag.c_str()),"Difference between #pi_1 #pi_2 and #pi_2 #pi_3 reconstructed  Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_pi1pi2_pi1pi3 = new TH1F(Form("%s_yvtxdiff_pi1pi2_pi1pi3",tag.c_str()),"Difference between  #pi_1 #pi_2 and #pi_1 #pi_3 reconstructed Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_pi1pi3_pi2pi3 = new TH1F(Form("%s_yvtxdiff_pi1pi3_pi2pi3",tag.c_str()),"Difference between  #pi_1 #pi_3 and #pi_2 #pi_3 reconstructed Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_zvtxdiff_pi1pi2_pi2pi3 = new TH1F(Form("%s_zvtxdiff_pi1pi2_pi2pi3",tag.c_str()),"Difference between #pi_1 #pi_2 and #pi_2 #pi_3 reconstructed  Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_zvtxdiff_pi1pi2_pi1pi3 = new TH1F(Form("%s_zvtxdiff_pi1pi2_pi1pi3",tag.c_str()),"Difference between  #pi_1 #pi_2 and #pi_1 #pi_3 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_zvtxdiff_pi1pi3_pi2pi3 = new TH1F(Form("%s_zvtxdiff_pi1pi3_pi2pi3",tag.c_str()),"Difference between  #pi_1 #pi_3 and #pi_2 #pi_3 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_cda_pi1_pi2 = new TH1F(Form("%s_cda_pi1_pi2",tag.c_str()),"Closest distance approach between the #pi_1 and #pi_3;cda[cm];Nevents",50, 0., 50.);
        fh_cda_pi1_pi3 = new TH1F(Form("%s_cda_pi1_pi3",tag.c_str()),"Closest distance approach between the #pi_1 and #pi_3;cda[cm];Nevents",50, 0., 50.);
        fh_cda_pi2_pi3 = new TH1F(Form("%s_cda_pi2_pi3",tag.c_str()),"Closest distance approach between the #pi_2 and #pi_3;cda[cm];Nevents" ,50, 0., 50.);
        //ENDOF Vertex Reconstruction

        fh_muee_P = new TH1F(Form("%s_muee_P",tag.c_str()),"Three track momentum ;P_{3#pi}[GeV];Nevents",100,0.,100);
        fh_muee_Pt = new TH1F(Form("%s_muee_Pt",tag.c_str()),"Transverse momentum of 3#pi ;Pt_{#mu e e}[GeV];Nevents",200,0,2.);
        fh_muee_M = new TH1F(Form("%s_muee_M",tag.c_str()),"Three track invariant mass ;M_{3#pi}[GeV];Nevents",1000,-1.,1.);

        fh_missing_mass = new TH1F(Form("%s_missing_mass",tag.c_str()),"Missing mass squared;M^{2}_{miss};Nevents",100,-0.05,0.05);

        //Pions
        fh_pion1_bx = new TH1F(Form("%s_pion1_bx",tag.c_str()),"#pi_1 x coordinate on DCH1 before magnet;DCH1_x[cm];Nevents",300,-150.,150.);
        fh_pion1_by = new TH1F(Form("%s_pion1_by",tag.c_str()),"#pi_1 y coordinate on DCH1 before magnet;DCH1_y[cm];Nevents",300,-150.,150.);
        fh_pion2_bx = new TH1F(Form("%s_pion2_bx",tag.c_str()),"#pi_2 x coordinate on DCH1 before magnet;DCH1_x[cm];Nevents",300,-150.,150.);
        fh_pion2_by = new TH1F(Form("%s_pion2_by",tag.c_str()),"#pi_2 y coordinate on DCH1 before magnet;DCH1_y[cm];Nevents",300,-150.,150.);
        fh_pion3_bx = new TH1F(Form("%s_pion3_bx",tag.c_str()),"#pi_3 x coordinate on DCH1 before magnet;DCH1_x[cm];Nevents",300,-150.,150.);
        fh_pion3_by = new TH1F(Form("%s_pion3_by",tag.c_str()),"#pi_3 y coordinate on DCH1 before magnet;DCH1_y[cm];Nevents",300,-150.,150.);
        fh_bx_vs_by_pi1 = new TH2F(Form("%s_bx_vs_by_pi1",tag.c_str()),"#pi_1 x vs y coordinate on DCH1 before the magnet;DCH_x[cm];DCH_y[cm]",300,-150.,150.,300,-150,150);
        fh_bx_vs_by_pi2 = new TH2F(Form("%s_bx_vs_by_pi2",tag.c_str()),"#pi_2 x vs y coordinate on DCH1 before the magnet;DCH_x[cm];DCH_y[cm]",300,-150.,150.,300,-150,150);
        fh_bx_vs_by_pi3 = new TH2F(Form("%s_bx_vs_by_pi3",tag.c_str()),"#pi_3 x vs y coordinate on DCH1 before the magnet;DCH_x[cm];DCH_y[cm]",300,-150.,150.,300,-150,150);
        fh_Lkr_extrap_tracks_x_vs_y = new TH2F(Form("%s_Lkr_extrap_tracks_x_vs_y",tag.c_str()),"Extrapolated tracks to the Lkr ;Lkr_x[cm];Lkr_y[cm]",400,-200., 200., 400, -200., 200.);
        fh_muv_x_y_position = new TH2F(Form("%s_muv_x_y_position",tag.c_str()),"MUV2 extrapolated position x vs y;Track_x[cm];Track_y[cm]",400,-200.,200.,400,-200.,200.);
        fh_MM2_plus = new TH1F(Form("%s_MM2_plus",tag.c_str()),"Missing mass squared for K+",100, -0.05, 0.05);
        fh_MM2_minus = new TH1F(Form("%s_MM2_minus",tag.c_str()),"Missing mass squared for K-",100, -0.05, 0.05);
    }
    if(ftype==2){
        fh_Ntracks           = new TH1I(Form("%s_Ntracks",tag.c_str()),"Number of charged tracks;Ntrk;Nevents",10,0,10);
        fh_Nclusters         = new TH1I(Form("%s_Nclusters",tag.c_str()),"Number of clusters in the Lkr;Nclusters;Nevents",10,0,10);
        fh_Nvtx              = new TH1I(Form("%s_Nvtx",tag.c_str()),"Number of vertexes;Nvtx;Nevents",10,0,10);
        fh_Kaon_Charge       = new TH1I(Form("%s_Kaon_Charge",tag.c_str()),"Charge of the Kaon;Q;Nevents",4,-2,2);
        fh_Mu_charge         = new TH1I(Form("%s_Mu_charge",tag.c_str()),"Muon charge ;Mu_Q;Nevents",4.,-2.,2.);
        fh_Event_Type        = new TH1I(Form("%s_Event_Type",tag.c_str()),"Mu e e : 0- ++-, 1- +++, 2- +--, 3- -+-, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_Track_Momentum    = new TH1F(Form("%s_Track_Momentum",tag.c_str()),"Track momentum ;Track_P[GeV];Nevents ",100.,0.,100.);
        fh_Mu_momentum       = new TH1F(Form("%s_Mu_momentum",tag.c_str()),"Muon momentum ;Mu_P[GeV];Nevents ",100.,0.,100.);
        fh_Electron_Momentum = new TH1F(Form("%s_Electron_Momentum",tag.c_str()),"Electron momentum ;El_P[GeV];Nevents",100,0.,100.);
        fh_eop               = new TH1F(Form("%s_Track_E/p",tag.c_str()),"Track E/p ;E/p;Nevents",120.,0.,1.2);
        fh_Mu_eop            = new TH1F(Form("%s_Mu_E/p",tag.c_str()),"Muon E/p ;E/p;Nevents",120.,0.,1.2);
        fh_Electron_eop      = new TH1F(Form("%s_Electron_E/p",tag.c_str()),"Electron E/p ;E/p;Nevents",120,0.,1.2);
        fh_odd_eop           = new TH1F(Form("%s_odd_E/p",tag.c_str()),"E/p for the odd track in the k3pi selection",120,0.,1.2);
        fh_COmPaCt_Z_Vertex  = new TH1F(Form("%s_COmPaCt_Z_Vertex",tag.c_str()),"Three track vertex from COmPaCt",150,-4000.,11000.);
        fh_Mu_Zvtx_min_COmPaCt_Zvtx = new TH1F(Form("%s_Mu_Zvtx_min_COmPaCt_Zvtx",tag.c_str()),"Calculated #mu Zvtx - COmPact three track Zvtx;(Z_#mu - Z_vtx)[cm];Nevents ",100,-10000.,10000.);
        //Reconstructed Vtx difference
        fh_xvtxdiff_mue1_mue2 = new TH1F(Form("%s_xvtxdiff_mue1_mue2",tag.c_str()),"Difference between #mu e1 and #mu e2 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_xvtxdiff_mue1_e1e2 = new TH1F(Form("%s_xvtxdiff_mue1_e1e2",tag.c_str()),"Difference between #mu e1 and e1 e2 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_xvtxdiff_mue2_e1e2 = new TH1F(Form("%s_xvtxdiff_mue2_e1e2",tag.c_str()),"Difference between #mu e2 and e1 e2 reconstructed X vtx;X_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_mue1_mue2 = new TH1F(Form("%s_yvtxdiff_mue1_mue2",tag.c_str()),"Difference between #mu e1 and #mu e2 reconstructed Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_mue1_e1e2 = new TH1F(Form("%s_yvtxdiff_mue1_e1e2",tag.c_str()),"Difference between #mu e1 and e1 e2 reconstructed  Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_mue2_e1e2 = new TH1F(Form("%s_yvtxdiff_mue2_e1e2",tag.c_str()),"Difference between #mu e2 and e1 e2 reconstructed  Y vtx;Y_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_zvtxdiff_mue1_mue2 = new TH1F(Form("%s_zvtxdiff_mue1_mue2",tag.c_str()),"Difference between #mu e1 and #mu e2 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_zvtxdiff_mue1_e1e2 = new TH1F(Form("%s_zvtxdiff_mue1_e1e2",tag.c_str()),"Difference between #mu e1 and e1 e2 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_zvtxdiff_mue2_e1e2 = new TH1F(Form("%s_zvtxdiff_mue2_e1e2",tag.c_str()),"Difference between #mu e2 and e1 e2 reconstructed Z vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);

        //ENDOF Reconstructed Vtx difference

        //Time matching
        fh_DCHtime_mu        = new TH1F(Form("%s_DCHtime_mu",tag.c_str()),"Drift Chamber time for the #mu ;DCHtime[ns];Nevents",200,50,250);
        fh_DCHtime_e1        = new TH1F(Form("%s_DCHtime_e1",tag.c_str()),"Drift Chamber time for the e1 ;DCHtime[ns];Nevents",200,50,250);
        fh_DCHtime_e2        = new TH1F(Form("%s_DCHtime_e2",tag.c_str()),"Drift Chamber time for the e2 ;DCHtime[ns];Nevents",200,50,250);
        fh_HodTime_mu        = new TH1F(Form("%s_HodTime_mu",tag.c_str()),"Hodoscope time for the #mu ;HodTime[ns];Nevents",200,50,250);
        fh_HodTime_e1        = new TH1F(Form("%s_HodTime_e1",tag.c_str()),"Hodoscope time for the e1 ;HodTime[ns];Nevents",200,50,250);
        fh_HodTime_e2        = new TH1F(Form("%s_HodTime_e2",tag.c_str()),"Hodoscope time for the e2 ;HodTime[ns];Nevents",200,50,250);
        fh_DCH_timediff_mu_e1 = new TH1F(Form("%s_DCH_timediff_mu_e1",tag.c_str()),"DCHtime difference between #mu and e1;DCHtimediff[ns];Nevents",100,-50.,50.);
        fh_DCH_timediff_mu_e2 = new TH1F(Form("%s_DCH_timediff_mu_e2",tag.c_str()),"DCHtime difference between #mu and e2;DCHtimediff[ns];Nevents",100,-50.,50.);
        fh_DCH_timediff_e1_e2 = new TH1F(Form("%s_DCH_timediff_e1_e2",tag.c_str()),"DCHtime difference between e1 and e2;DCHtimediff[ns];Nevents",100,-50.,50.);
        fh_Hod_timediff_mu_e1 = new TH1F(Form("%s_Hod_timediff_mu_e1",tag.c_str()),"Hodoscope time difference between #mu and e1;Hodtimediff[ns];Nevents",100,-50.,50.);
        fh_Hod_timediff_mu_e2 = new TH1F(Form("%s_Hod_timediff_mu_e2",tag.c_str()),"Hodoscope time difference between #mu and e2;Hodtimediff[ns];Nevents",100,-50.,50.);
        fh_Hod_timediff_e1_e2 = new TH1F(Form("%s_Hod_timediff_e1_e2",tag.c_str()),"Hodoscope time difference between e1 and e2 ;Hodtimediff[ns];Nevents",100,-50.,50.);
        //ENDOF Time matching
        fh_mee = new TH1F(Form("%s_mee",tag.c_str()),"Invariant mass of the electron pair ",50,0.,0.5);
        fh_mee_z_variable = new TH1F(Form("%s_mee_z_variable",tag.c_str()),"Invariant mass of the electron pair in terms of the z variable;z;Nevents",25,0,0.5);

        //2D Mee and MM2 plots
        fh_mee_vs_Pmu  = new TH2F(Form("%s_mee_vs_Pmu",tag.c_str()),"Invariant mass of the electron pair vs #mu momentum ",50,0.,0.5, 100, 0, 100);
        fh_mee_vs_Pel1 = new TH2F(Form("%s_mee_vs_Pel1",tag.c_str()),"Invariant mass of the electron pair vs electron1 momentum ",50,0.,0.5, 100, 0, 100);
        fh_mee_vs_Pel2 = new TH2F(Form("%s_mee_vs_Pel2",tag.c_str()),"Invariant mass of the electron pair vs electron1 momentum ",50,0.,0.5, 100, 0, 100);
        fh_mee_vs_Pt   = new TH2F(Form("%s_mee_vs_Pt",tag.c_str()),"Invariant mass of the electron pair vs Total transverse momentum ",50,0.,0.5, 100, 0, 1);
        fh_mee_vs_MM2  = new TH2F(Form("%s_mee_vs_MM2",tag.c_str()),"Invariant mass of the electron pair vs MM2 ",50,0.,0.5, 100, -0.05, 0.05);
        fh_mee_vs_zvtx = new TH2F(Form("%s_mee_vs_zvtx",tag.c_str()),"Invariant mass of the electron pair vs zvtx ",50,0.,0.5, 100, -4000, 11000);
        fh_P_vs_dtrkcl = new TH2F(Form("%s_P_vs_dtrkcl",tag.c_str()),"Momentum vs distance between track and cluster ",100,0.,100,200,0.,200.);
        fh_Ecl_vs_dtrkcl = new TH2F(Form("%s_mee_vs_zvtx",tag.c_str()),"Invariant mass of the electron pair vs zvtx ",100,0.,100,100,0.,50.);

        fh_HoDTotTime = new TH1F(Form("%s_HoDTotTime",tag.c_str()),"(T^{HoD}_{1} + T^{HoD}_{2} + T^{HoD}_{3} )/3;EventTime [ns];Nevents",200,0.,200.);

        fh_muee_P = new TH1F(Form("%s_muee_P",tag.c_str()),"Three track momentum ;P_{#mu e e}[GeV];Nevents",100,0.,100);
        fh_muee_Pt = new TH1F(Form("%s_muee_Pt",tag.c_str()),"Transverse momentum of #mu^{#pm} e^{+}e^{-} system;Pt_{#mu e e}[GeV];Nevents",200,0,2.);
        fh_muee_M = new TH1F(Form("%s_muee_M",tag.c_str()),"Three track invariant mass ;M_{#mu e e}[GeV];Nevents",100,0.,1.);
        fh_Muee_M_3pi_assumption = new TH1F(Form("%s_Muee_M_3pi_assumption",tag.c_str()),"Three track invariant mass with 3 #pi assumtion ;M_{3#pi};Nevents",100,0.,1.);

        fh_missing_mass = new TH1F(Form("%s_missing_mass",tag.c_str()),"Missing mass squared;M^{2}_{miss}[GeV^2];Nevents",100,-0.05,0.05);
        fh_cda_mu_e1 = new TH1F(Form("%s_cda_mu_e1",tag.c_str()),"Closest distance approach between the #mu and e1;cda[cm];Nevents",50, 0., 50.);
        fh_cda_mu_e2 = new TH1F(Form("%s_cda_mu_e2",tag.c_str()),"Closest distance approach between the #mu and e2;cda[cm];Nevents",50, 0., 50.);
        fh_cda_e1_e2 = new TH1F(Form("%s_cda_e1_e2",tag.c_str()),"Closest distance approach between the e1 and e2;cda[cm];Nevents" ,50, 0., 50.);
        //Drift chambers geometry
        //Muons
        fh_muon_bx = new TH1F(Form("%s_muon_bx",tag.c_str()),"#mu x coordinate on DCH1 before magnet;DCH1_x[cm];Nevents",300,-150,150);
        fh_muon_by = new TH1F(Form("%s_muon_by",tag.c_str()),"#mu y coordinate on DCH1 before magnet;DCH1_y[cm];Nevents",300,-150,150);
        //Electrons
        fh_electron1_bx = new TH1F(Form("%s_electron1_bx",tag.c_str()),"el1 x coordinate on DCH1 before magnet;DCH1_x[cm];Nevents",300,-150.,150.);
        fh_electron1_by = new TH1F(Form("%s_electron1_by",tag.c_str()),"el1 y coordinate on DCH1 before magnet;DCH1_y[cm];Nevents",300,-150.,150.);
        fh_electron2_bx = new TH1F(Form("%s_electron2_bx",tag.c_str()),"el2 x coordinate on DCH1 before magnet;DCH1_x[cm];Nevents",300,-150.,150.);
        fh_electron2_by = new TH1F(Form("%s_electron2_by",tag.c_str()),"el2 y coordinate on DCH1 before magnet;DCH1_y[cm];Nevents",300,-150.,150.);
        fh_bx_vs_by_muon = new TH2F(Form("%s_bx_vs_by_muon",tag.c_str()),"#mu x vs y coordinate on DCH1 before the magnet;DCH_x[cm];DCH_y[cm]",300,-150.,150.,300,-150,150);
        fh_bx_vs_by_el1 = new TH2F(Form("%s_bx_vs_by_el1",tag.c_str()),"el1 x vs y coordinate on DCH1 before the magnet;DCH_x[cm];DCH_y[cm]",300,-150.,150.,300,-150,150);
        fh_bx_vs_by_el2 = new TH2F(Form("%s_bx_vs_by_el2",tag.c_str()),"el2 x vs y coordinate on DCH1 before the magnet;DCH_x[cm];DCH_y[cm]",300,-150.,150.,300,-150,150);
        fh_el1_cluster_x = new TH1F(Form("%s_el1_cluster_x",tag.c_str()),"x position of el1 cluster;Cl_x[cm];Nevents",400,-200.,200.);
        fh_el1_cluster_y = new TH1F(Form("%s_el1_cluster_y",tag.c_str()),"y position of el1 cluster;Cl_y[cm];Nevents",400,-200.,200.);
        fh_el2_cluster_x = new TH1F(Form("%s_el2_cluster_x",tag.c_str()),"x position of el2 cluster;Cl_x[cm];Nevents",400,-200.,200.);
        fh_el2_cluster_y = new TH1F(Form("%s_el2_cluster_y",tag.c_str()),"y position of el2 cluster;Cl_y[cm];Nevents",400,-200.,200.);
        fh_mu_cluster_x = new TH1F(Form("%s_mu_cluster_x",tag.c_str()),"x position of #mu cluster;Cl_x[cm];Nevents",400,-200.,200.);
        fh_mu_cluster_y = new TH1F(Form("%s_mu_cluster_y",tag.c_str()),"y position of #mu cluster;Cl_y[cm];Nevents",400,-200.,200.);
        fh_el1_cluster_time = new TH1F(Form("%s_el1_cluster_time",tag.c_str()),"Time of the el1 cluster;Cl_time[ns];Nevents",200,0.,200.);
        fh_el2_cluster_time = new TH1F(Form("%s_el2_cluster_time",tag.c_str()),"Time of the el2 cluster;Cl_time[ns];Nevents",200,0.,200.);
        fh_el1_el2_cluster_timediff = new TH1F(Form("%s_el1_el2_cluster_timediff",tag.c_str()),"Time difference between the clustes of the two electrons;Tdiff_{e1e2}[ns];Nevents",200,-50.,50.);
        fh_el1_dtrk_cl = new TH1F(Form("%s_el1_dtrk_cl",tag.c_str()),"Difference between extrapolated track from DCH and to the Lkr and cluster for el1   ;dtrk_cl[cm];Nevents",200,0.,200.);
        fh_el2_dtrk_cl = new TH1F(Form("%s_el2_dtrk_cl",tag.c_str()),"Difference between extrapolated track from DCH and to the Lkr and cluster for el2   ;dtrk_cl[cm];Nevents",200,0.,200.);
        fh_muon_dtrk_cl = new TH1F(Form("%s_muon_dtrk_cl",tag.c_str()),"Difference between extrapolated track from DCH and to the Lkr and cluster for #mu ;dtrk_cl[cm];Nevents",200,0.,200.);
        //fh_el1_cluster_x_y = new TH2F(Form("%s_el1_cluster_x_y",tag.c_str()),"Cluster position x vs y for el1;cl_x[cm];cl_y[cm]",400,-200.,200.,400,-200.,-200.);
        //fh_el2_cluster_x_y = new TH2F(Form("%s_el2_cluster_x_y",tag.c_str()),"Cluster position x vs y for el2;cl_x[cm];cl_y[cm]",400,-200.,200.,400,-200.,-200.);
        //fh_mu_cluster_x_y = new TH2F(Form("%s_mu_cluster_x_y",tag.c_str()),"Cluster position x vs y;cl_x[cm];cl_y[cm]",400,-200.,200.,400,-200.,-200.);
        fh_deadcell_distance = new TH1F(Form("%s_deadcell_distance",tag.c_str()),"Distance to deadcell",150,0.,150.);
        fh_muv2_trk_cl_diff = new TH1F(Form("%s_muv2_trk_cl_diff",tag.c_str()),"Difference between extrapolated track and MUV2 cluster position;Trk_cl_diff[cm];Nevents",200,0.,200.);
        fh_muv_xpos = new TH1F(Form("%s_muv_xpos",tag.c_str()),"Muon X position at MUV2",400,-200.,200.);
        fh_muv_ypos = new TH1F(Form("%s_muv_ypos",tag.c_str()),"Muon Y position at MUV1",400,-200.,200.);
        fh_muv_x_y_position = new TH2F(Form("%s_muv_x_y_position",tag.c_str()),"MUV2 extrapolated position x vs y;Track_x[cm];Track_y[cm]",400,-200.,200.,400,-200.,200.);
        fh_DCH1_distance_e1e2 = new TH1F(Form("%s_DCH1_distance_e1e2",tag.c_str()),"Distance between the two electrons in the DCH1",200,0,200.);
        fh_Lkr_distance_e1e2 = new TH1F(Form("%s_Lkr_distance_e1e2",tag.c_str()),"Distance between the two electrons in the Lkr",200,0,200.);
        fh_Lkr_extrap_tracks_x_vs_y = new TH2F(Form("%s_Lkr_extrap_tracks_x_vs_y",tag.c_str()),"Extrapolated tracks to the Lkr ;Lkr_x[cm];Lkr_y[cm]",400,-200., 200., 400, -200., 200.);
        fh_angle_el_mu = new TH1F(Form("%s_angle_el_mu",tag.c_str()),"Angle between e+e- plane and the #mu",200,0.0,0.2);

        //MC true histograms
        fh_mc_KDzvtx = new TH1F(Form("%s_mc_KDzvtx",tag.c_str()),"MC decayed kaon Z vertex",150,-4000.,11000.);
        fh_mc_P1_Pzvtx = new TH1F(Form("%s_mc_P1_Pzvtx",tag.c_str()),"MC Particle 1 production Z vertex",150,-4000.,11000.);
        fh_mc_P1_mass = new TH1F(Form("%s_mc_P1_mass",tag.c_str()),"MC Particle 1 Mass",550,0.,0.55);
        fh_mc_P1_momentum = new TH1F(Form("%s_mc_P1_momentum",tag.c_str()),"MC Particle 1 Momentum",100,0.,100.);
        fh_mc_P2_Pzvtx = new TH1F(Form("%s_mc_P2_Pzvtx",tag.c_str()),"MC Particle 2 production Z vertex",150,-4000.,11000.);
        fh_mc_P2_mass = new TH1F(Form("%s_mc_P2_mass",tag.c_str()),"MC Particle 2 Mass",1100,0.,0.55);
        fh_mc_P2_momentum = new TH1F(Form("%s_mc_P2_momentum",tag.c_str()),"MC Particle 2 Momentum",100,0.,100.);
        fh_mc_P1_dist_prod_dec = new TH1F(Form("%s_mc_P1_dist_prod_dec",tag.c_str()),"MC Particle1 distance between production and decay vertex; P1Zvtx_{Production - Decay};Nevents",-150,-4000,11000);
        fh_mc_P2_dist_prod_dec = new TH1F(Form("%s_mc_P2_dist_prod_dec",tag.c_str()),"MC Particle2 distance between production and decay vertex; P2Zvtx_{Production - Decay};Nevents",-150,-4000,11000);
        fh_mc_P3_dist_prod_dec = new TH1F(Form("%s_mc_P3_dist_prod_dec",tag.c_str()),"MC Particle3 distance between production and decay vertex; P3Zvtx_{Production - Decay};Nevents",-150,-4000,11000);
        fh_mc_P3_Pzvtx = new TH1F(Form("%s_mc_P3_Pzvtx",tag.c_str()),"MC Particle 3 production Z vertex",150,-4000.,11000.);
        fh_mc_P3_mass = new TH1F(Form("%s_mc_P3_mass",tag.c_str()),"MC Particle 3 Mass",1100,0.,0.55);
        fh_mc_P3_momentum = new TH1F(Form("%s_mc_P3_momentum",tag.c_str()),"MC Particle 3 Momentum",100,0.,100.);
        fh_mc_P4_Pzvtx = new TH1F(Form("%s_mc_P4_Pzvtx",tag.c_str()),"MC Particle 4 production Z vertex",150,-4000.,11000.);
        fh_mc_three_track_123_momentum = new TH1F(Form("%s_mc_three_track_123_momentum",tag.c_str()),"MC Momentum of the track 1+2+3",100,0.,100.);
        fh_mc_three_track_123_mass = new TH1F(Form("%s_mc_three_track_123_mass",tag.c_str()),"MC Mass of the track 1+2+3 ",70,0,0.7);
        fh_mc_two_track_23_momentum = new TH1F(Form("%s_mc_two_track_23_momentum",tag.c_str()),"MC Momentum of the track 2+3",100,0.,100.);
        fh_mc_two_track_23_mass = new TH1F(Form("%s_mc_two_track_23_mass",tag.c_str()),"MC Mass of the track 2+3;M_{ee}[GeV]; ",50,0,0.5);
        fh_mc_two_track_23_mass_z_variable = new TH1F(Form("%s_mc_two_track_23_mass_z_variable",tag.c_str()),"MC Mass of the track 2+3 in terms of the z variable;z;Nevents",25,0,0.5);
        fh_mc_P4_mass = new TH1F(Form("%s_mc_P4_mass",tag.c_str()),"MC Particle 4 Mass",550,0.,0.55);
        fh_mc_P4_momentum = new TH1F(Form("%s_mc_P4_momentum",tag.c_str()),"MC Particle 4 Momentum",100,0.,100.);
        fh_mc_four_track_1234_momentum = new TH1F(Form("%s_mc_four_track_1234_momentum",tag.c_str()),"MC Momentum of the track 1+2+3+4",100,0.,100.);
        fh_mc_four_track_1234_mass = new TH1F(Form("%s_mc_four_track_1234_mass",tag.c_str()),"MC Mass of the track 1+2+3+4 ",70,0,0.7);
        fh_lda3_e1 = new TH1F(Form("%s_lda3_e1",tag.c_str()),"lda3 variable for electron 1", 200,0.,2);
        fh_lda3_e2 = new TH1F(Form("%s_lda3_e2",tag.c_str()),"lda3 variable for electron 2", 200,0.,2);
        fh_Mee_before = new TH1F(Form("%s_z_vtx_before",tag.c_str()),"ComPaCt z vtx before z vtx cut",150,-4000.,11000.);
        fh_muon_status = new TH1I(Form("%s_muon_status",tag.c_str()),"1 - MUV 1 & 2 & 3, 2 - 1 & 2 &! 3, 3 - !1 &2 &3, 4 - 1 & &!2 &3", 10, 0, 10);
        fh_MM2_plus = new TH1F(Form("%s_MM2_plus",tag.c_str()),"Missing mass squared for K+",100, -0.05, 0.05);
        fh_MM2_minus = new TH1F(Form("%s_MM2_minus",tag.c_str()),"Missing mass squared for K-",100, -0.05, 0.05);

    }
    if(ftype==3){
        fh_missing_mass = new TH1F(Form("%s_missing_mass",tag.c_str()),"Missing mass squared;M^{2}_{miss};Nevents",100,-0.05,0.05);
        fh_mee = new TH1F(Form("%s_mee",tag.c_str()),"Invariant mass of the electron pair ",50,0.,0.5);
        fh_mee_z_variable = new TH1F(Form("%s_mee_z_variable",tag.c_str()),"Invariant mass of the electron pair in terms of the z variable;z;Nevents",25,0,0.5);
        fh_MM2_plus = new TH1F(Form("%s_MM2_plus",tag.c_str()),"Missing mass squared for K+",100, -0.05, 0.05);
        fh_MM2_minus = new TH1F(Form("%s_MM2_minus",tag.c_str()),"Missing mass squared for K-",100, -0.05, 0.05);
    }
}


void Hist_dir::AddToFile(TFile* file){
    TDirectory* dir = file->mkdir(fDir.c_str());
    dir->cd();
    if(ftype==0){
        fh_Nrun              ->Write();
        fh_Ntracks           ->Write();
        fh_Nclusters         ->Write();
        fh_Nvtx              ->Write();
        fh_Track_Momentum    ->Write();
        fh_eop               ->Write();
        fh_COmPaCt_Z_Vertex  ->Write();
        fh_mc_KDzvtx->Write();
        fh_mc_P1_Pzvtx->Write();
        fh_mc_P1_mass->Write();
        fh_mc_P1_momentum->Write();
        fh_mc_three_track_123_momentum->Write();
        //fh_mc_P2_mass->Write();
        fh_mc_P2_Pzvtx->Write();
        fh_mc_P2_momentum->Write();
        fh_mc_P3_Pzvtx->Write();
        fh_mc_P3_mass->Write();
        fh_mc_P3_momentum->Write();
        fh_mc_three_track_123_mass->Write();
        fh_mc_two_track_23_momentum->Write();
        fh_mc_two_track_23_mass->Write();
        fh_mc_two_track_23_mass_z_variable->Write();
        fh_mc_P4_Pzvtx->Write();
        fh_mc_P4_mass->Write();
        fh_mc_P4_momentum->Write();
        fh_mc_four_track_1234_momentum->Write();
        fh_mc_four_track_1234_mass->Write();

    }
    if(ftype==1){
        fh_Trigger_bits      ->Write();
        fh_SS0_CPRE         ->Write();
        fh_SS0_MB_1TRK_P    ->Write();
        fh_SS0_MB_1VTX      ->Write();
        fh_SS0_full_trig    ->Write();
        fh_SS1_CPRE         ->Write();
        fh_SS1_MB_1TRK_P    ->Write();
        fh_SS1_MB_1VTX      ->Write();
        fh_SS1_full_trig    ->Write();
        fh_SS2_CPRE         ->Write();
        fh_SS2_MB_1TRK_P    ->Write();
        fh_SS2_MB_1VTX      ->Write();
        fh_SS2_full_trig    ->Write();
        fh_SS3_CPRE         ->Write();
        fh_SS3_MB_1TRK_P    ->Write();
        fh_SS3_MB_1VTX      ->Write();
        fh_SS3_full_trig    ->Write();
        fh_SS4_CPRE         ->Write();
        fh_SS4_MB_1TRK_P    ->Write();
        fh_SS4_MB_1VTX      ->Write();
        fh_SS4_full_trig    ->Write();
        fh_SS5_CPRE         ->Write();
        fh_SS5_MB_1TRK_P    ->Write();
        fh_SS5_MB_1VTX      ->Write();
        fh_SS5_full_trig    ->Write();
        fh_SS6_CPRE         ->Write();
        fh_SS6_MB_1TRK_P    ->Write();
        fh_SS6_MB_1VTX      ->Write();
        fh_SS6_full_trig    ->Write();
        fh_SS7_CPRE         ->Write();
        fh_SS7_MB_1TRK_P    ->Write();
        fh_SS7_MB_1VTX      ->Write();
        fh_SS7_full_trig    ->Write();
        fh_SS8_CPRE         ->Write();
        fh_SS8_MB_1TRK_P    ->Write();
        fh_SS8_MB_1VTX      ->Write();
        fh_SS8_full_trig    ->Write();
        fh_Ntracks           ->Write();
        fh_Nvtx              ->Write();
        fh_Nclusters         ->Write();
        fh_Kaon_Charge       ->Write();
        fh_Event_Type        ->Write();
        fh_Pion_Momentum     ->Write();
        fh_eop               ->Write();
        fh_odd_eop           ->Write();
        fh_EoP_vs_p_odd_tr   ->Write();
        fh_lda3_p1->Write();
        fh_lda3_p2->Write();
        fh_lda3_p3->Write() ;
        fh_COmPaCt_Z_Vertex  ->Write();
        fh_zvtxdiff_pi1pi2_pi2pi3->Write();
        fh_zvtxdiff_pi1pi2_pi1pi3->Write();
        fh_zvtxdiff_pi1pi3_pi2pi3->Write();
        fh_xvtxdiff_pi1pi2_pi2pi3->Write();
        fh_xvtxdiff_pi1pi2_pi1pi3->Write();
        fh_xvtxdiff_pi1pi3_pi2pi3->Write();
        fh_yvtxdiff_pi1pi2_pi2pi3->Write();
        fh_yvtxdiff_pi1pi2_pi1pi3->Write();
        fh_yvtxdiff_pi1pi3_pi2pi3->Write();
        fh_cda_pi1_pi2           ->Write();
        fh_cda_pi1_pi3           ->Write();
        fh_cda_pi2_pi3           ->Write();
        fh_pion1_bx              ->Write();
        fh_pion1_by              ->Write();
        fh_pion2_bx              ->Write();
        fh_pion2_by              ->Write();
        fh_pion3_bx              ->Write();
        fh_pion3_by              ->Write();
        fh_muee_P                ->Write();
        fh_muee_Pt               ->Write();
        fh_muee_M                ->Write();
        fh_missing_mass          ->Write();
        fh_bx_vs_by_pi1->Write();
        fh_bx_vs_by_pi2->Write();
        fh_bx_vs_by_pi3->Write();
        fh_Lkr_extrap_tracks_x_vs_y->Write();
        fh_muv_x_y_position->Write();
        fh_MM2_plus->Write();
        fh_MM2_minus->Write();
    }
    if(ftype==2){
        fh_Ntracks           ->Write();
        fh_Nvtx              ->Write();
        fh_Nclusters         ->Write();
        fh_Kaon_Charge       ->Write();
        fh_Mu_charge         ->Write();
        fh_Event_Type        ->Write();
        fh_Track_Momentum    ->Write();
        fh_Mu_momentum       ->Write();
        fh_Electron_Momentum ->Write();
        fh_eop               ->Write();
        fh_Mu_eop            ->Write();
        fh_Electron_eop      ->Write();
        fh_odd_eop           ->Write();
        fh_COmPaCt_Z_Vertex  ->Write();
        fh_Mu_Zvtx_min_COmPaCt_Zvtx->Write();
        fh_zvtxdiff_mue1_mue2->Write();
        fh_zvtxdiff_mue1_e1e2->Write();
        fh_zvtxdiff_mue2_e1e2->Write();
        fh_xvtxdiff_mue1_mue2->Write();
        fh_xvtxdiff_mue1_e1e2->Write();
        fh_xvtxdiff_mue2_e1e2->Write();
        fh_yvtxdiff_mue1_mue2->Write();
        fh_yvtxdiff_mue1_e1e2->Write();
        fh_yvtxdiff_mue2_e1e2->Write();
        fh_cda_mu_e1->Write();
        fh_cda_mu_e2->Write();
        fh_cda_e1_e2->Write();
        fh_muon_bx      ->Write();
        fh_muon_by      ->Write();
        fh_electron1_bx ->Write();
        fh_electron1_by ->Write();
        fh_electron2_bx ->Write();
        fh_electron2_by ->Write();
        fh_DCHtime_mu->Write();
        fh_DCHtime_e1->Write();
        fh_DCHtime_e2->Write();
        fh_HodTime_mu->Write();
        fh_HodTime_e1->Write();
        fh_HodTime_e2->Write();
        fh_DCH_timediff_mu_e1->Write();
        fh_DCH_timediff_mu_e2->Write();
        fh_DCH_timediff_e1_e2->Write();
        fh_Hod_timediff_mu_e1->Write();
        fh_Hod_timediff_mu_e2->Write();
        fh_Hod_timediff_e1_e2->Write();
        fh_HoDTotTime->Write();
        fh_mee         ->Write();
        fh_mee_z_variable->Write();
        fh_muee_P      ->Write();
        fh_muee_M      ->Write();
        fh_Muee_M_3pi_assumption->Write();
        fh_missing_mass->Write();
        fh_mee_vs_Pmu  ->Write();
        fh_mee_vs_Pel1 ->Write();
        fh_mee_vs_Pel2 ->Write();
        fh_mee_vs_Pt   ->Write();
        fh_mee_vs_MM2  ->Write();
        fh_mee_vs_zvtx ->Write();
        fh_P_vs_dtrkcl ->Write();
        fh_Ecl_vs_dtrkcl ->Write();
        fh_muee_Pt      ->Write();
        fh_bx_vs_by_muon->Write();
        fh_bx_vs_by_el1 ->Write();
        fh_bx_vs_by_el2 ->Write();
        fh_el1_cluster_x->Write();
        fh_el1_cluster_y->Write();
        fh_el2_cluster_x->Write();
        fh_el2_cluster_y->Write();
        fh_mu_cluster_x->Write();
        fh_mu_cluster_y->Write();
        fh_el1_cluster_time->Write();
        fh_el2_cluster_time->Write();
        fh_el1_el2_cluster_timediff->Write();
        fh_el1_dtrk_cl->Write();
        fh_el2_dtrk_cl->Write();
        fh_muon_dtrk_cl->Write();
        fh_muon_status->Write();
        //fh_el1_cluster_x_y->Write();
        //fh_el2_cluster_x_y->Write();
        //fh_mu_cluster_x_y->Write();
        fh_deadcell_distance->Write();
        fh_muv2_trk_cl_diff->Write();
        fh_muv_xpos->Write();
        fh_muv_ypos->Write();
        fh_muv_x_y_position->Write();
        fh_DCH1_distance_e1e2->Write();
        fh_Lkr_distance_e1e2->Write();
        fh_Lkr_extrap_tracks_x_vs_y->Write();
        fh_angle_el_mu->Write();
        fh_mc_KDzvtx->Write();
        fh_mc_P1_Pzvtx->Write();
        fh_mc_P1_mass->Write();
        fh_mc_P1_momentum->Write();
        fh_mc_three_track_123_momentum->Write();
        fh_mc_P2_mass->Write();
        fh_mc_P2_Pzvtx->Write();
        fh_mc_P2_momentum->Write();
        fh_mc_P3_Pzvtx->Write();
        fh_mc_P3_mass->Write();
        fh_mc_P3_momentum->Write();
        fh_mc_three_track_123_mass->Write();
        fh_mc_two_track_23_momentum->Write();
        fh_mc_two_track_23_mass->Write();
        fh_mc_two_track_23_mass_z_variable->Write();
        fh_mc_P4_Pzvtx->Write();
        fh_mc_P4_mass->Write();
        fh_mc_P4_momentum->Write();
        fh_mc_four_track_1234_momentum->Write();
        fh_mc_four_track_1234_mass->Write();
        fh_mc_P1_dist_prod_dec->Write();
        fh_mc_P2_dist_prod_dec->Write();
        fh_mc_P3_dist_prod_dec->Write();
        fh_lda3_e1->Write();
        fh_lda3_e2->Write();
        fh_Mee_before->Write();
        fh_MM2_plus->Write();
        fh_MM2_minus->Write();
        fh_SS0_CPRE         ->Write();
        fh_SS0_MB_1TRK_P    ->Write();
        fh_SS0_MB_1VTX      ->Write();
        fh_SS0_full_trig    ->Write();
        fh_SS1_CPRE         ->Write();
        fh_SS1_MB_1TRK_P    ->Write();
        fh_SS1_MB_1VTX      ->Write();
        fh_SS1_full_trig    ->Write();
        fh_SS2_CPRE         ->Write();
        fh_SS2_MB_1TRK_P    ->Write();
        fh_SS2_MB_1VTX      ->Write();
        fh_SS2_full_trig    ->Write();
        fh_SS3_CPRE         ->Write();
        fh_SS3_MB_1TRK_P    ->Write();
        fh_SS3_MB_1VTX      ->Write();
        fh_SS3_full_trig    ->Write();
        fh_SS4_CPRE         ->Write();
        fh_SS4_MB_1TRK_P    ->Write();
        fh_SS4_MB_1VTX      ->Write();
        fh_SS4_full_trig    ->Write();
        fh_SS5_CPRE         ->Write();
        fh_SS5_MB_1TRK_P    ->Write();
        fh_SS5_MB_1VTX      ->Write();
        fh_SS5_full_trig    ->Write();
        fh_SS6_CPRE         ->Write();
        fh_SS6_MB_1TRK_P    ->Write();
        fh_SS6_MB_1VTX      ->Write();
        fh_SS6_full_trig    ->Write();
        fh_SS7_CPRE         ->Write();
        fh_SS7_MB_1TRK_P    ->Write();
        fh_SS7_MB_1VTX      ->Write();
        fh_SS7_full_trig    ->Write();
        fh_SS8_CPRE         ->Write();
        fh_SS8_MB_1TRK_P    ->Write();
        fh_SS8_MB_1VTX      ->Write();
        fh_SS8_full_trig    ->Write();
    }
    if(ftype==3){
        fh_missing_mass          ->Write();
        fh_mee         ->Write();
        fh_mee_z_variable->Write();
        fh_MM2_plus->Write();
        fh_MM2_minus->Write();
    }
}

Hist_dir::~Hist_dir(){
    if(ftype==0){
        delete fh_Nrun;
        delete fh_Ntracks           ;
        delete fh_Nclusters         ;
        delete fh_Nvtx              ;
        delete fh_Track_Momentum    ;
        delete fh_eop               ;
        delete fh_COmPaCt_Z_Vertex  ;
        delete fh_mc_KDzvtx;
        delete fh_mc_P1_Pzvtx;
        delete fh_mc_P1_mass;
        delete fh_mc_P1_momentum;
        delete fh_mc_three_track_123_momentum;
        delete fh_mc_P2_mass;
        delete fh_mc_P2_Pzvtx;
        delete fh_mc_P2_momentum;
        delete fh_mc_P3_Pzvtx;
        delete fh_mc_P3_mass;
        delete fh_mc_P3_momentum;
        delete fh_mc_three_track_123_mass;
        delete fh_mc_two_track_23_momentum;
        delete fh_mc_two_track_23_mass;
        delete fh_mc_P4_Pzvtx;
        delete fh_mc_P4_mass;
        delete fh_mc_P4_momentum;
        delete fh_mc_four_track_1234_momentum;
        delete fh_mc_four_track_1234_mass;


    }
    if(ftype==1){
        delete fh_Ntracks           ;
        delete fh_Nclusters         ;
        delete fh_Nvtx              ;
        delete fh_Kaon_Charge       ;
        delete fh_Event_Type        ;
        delete fh_Pion_Momentum     ;
        delete fh_eop               ;
        delete fh_odd_eop           ;
        delete fh_COmPaCt_Z_Vertex  ;
        delete fh_zvtxdiff_pi1pi2_pi2pi3   ;
        delete fh_zvtxdiff_pi1pi2_pi1pi3   ;
        delete fh_zvtxdiff_pi1pi3_pi2pi3   ;
        delete fh_xvtxdiff_pi1pi2_pi2pi3   ;
        delete fh_xvtxdiff_pi1pi2_pi1pi3   ;
        delete fh_xvtxdiff_pi1pi3_pi2pi3   ;
        delete fh_yvtxdiff_pi1pi2_pi2pi3   ;
        delete fh_yvtxdiff_pi1pi2_pi1pi3   ;
        delete fh_yvtxdiff_pi1pi3_pi2pi3   ;
        delete fh_cda_pi1_pi2              ;
        delete fh_cda_pi1_pi3              ;
        delete fh_cda_pi2_pi3              ;
        delete fh_muee_P                ;
        delete fh_muee_Pt               ;
        delete fh_muee_M                ;
        delete fh_missing_mass          ;
        delete fh_pion1_bx;
        delete fh_pion1_by;
        delete fh_pion2_bx;
        delete fh_pion2_by;
        delete fh_pion3_bx;
        delete fh_pion3_by;
        delete fh_bx_vs_by_pi1;
        delete fh_bx_vs_by_pi2;
        delete fh_bx_vs_by_pi3;
        delete fh_Trigger_bits;
        delete fh_SS0_CPRE     ;
        delete fh_SS0_MB_1TRK_P;
        delete fh_SS0_MB_1VTX  ;
        delete fh_SS0_full_trig;
        delete fh_SS1_CPRE     ;
        delete fh_SS1_MB_1TRK_P;
        delete fh_SS1_MB_1VTX  ;
        delete fh_SS1_full_trig;
        delete fh_SS2_CPRE     ;
        delete fh_SS2_MB_1TRK_P;
        delete fh_SS2_MB_1VTX  ;
        delete fh_SS2_full_trig;
        delete fh_SS3_CPRE     ;
        delete fh_SS3_MB_1TRK_P;
        delete fh_SS3_MB_1VTX  ;
        delete fh_SS3_full_trig;
        delete fh_SS4_CPRE     ;
        delete fh_SS4_MB_1TRK_P;
        delete fh_SS4_MB_1VTX  ;
        delete fh_SS4_full_trig;
        delete fh_SS5_CPRE     ;
        delete fh_SS5_MB_1TRK_P;
        delete fh_SS5_MB_1VTX  ;
        delete fh_SS5_full_trig;
        delete fh_SS6_CPRE     ;
        delete fh_SS6_MB_1TRK_P;
        delete fh_SS6_MB_1VTX  ;
        delete fh_SS6_full_trig;
        delete fh_SS7_CPRE     ;
        delete fh_SS7_MB_1TRK_P;
        delete fh_SS7_MB_1VTX  ;
        delete fh_SS7_full_trig;
        delete fh_SS8_CPRE     ;
        delete fh_SS8_MB_1TRK_P;
        delete fh_SS8_MB_1VTX  ;
        delete fh_SS8_full_trig;
        delete fh_MM2_minus;
        delete fh_MM2_plus;
    }
    if(ftype==2){
        delete fh_Ntracks           ;
        delete fh_Nclusters         ;
        delete fh_Nvtx              ;
        delete fh_Kaon_Charge       ;
        delete fh_Mu_charge         ;
        delete fh_Event_Type        ;
        delete fh_Track_Momentum    ;
        delete fh_Mu_momentum       ;
        delete fh_Electron_Momentum ;
        delete fh_eop               ;
        delete fh_Mu_eop            ;
        delete fh_Electron_eop      ;
        delete fh_odd_eop           ;
        delete fh_COmPaCt_Z_Vertex  ;
        delete fh_Mu_Zvtx_min_COmPaCt_Zvtx;
        delete fh_zvtxdiff_mue1_mue2;
        delete fh_zvtxdiff_mue1_e1e2;
        delete fh_zvtxdiff_mue2_e1e2;
        delete fh_xvtxdiff_mue1_mue2;
        delete fh_xvtxdiff_mue1_e1e2;
        delete fh_xvtxdiff_mue2_e1e2;
        delete fh_yvtxdiff_mue1_mue2;
        delete fh_yvtxdiff_mue1_e1e2;
        delete fh_yvtxdiff_mue2_e1e2;
        delete fh_DCHtime_mu;
        delete fh_HodTime_mu;
        delete fh_DCHtime_e1;
        delete fh_DCHtime_e2;
        delete fh_HodTime_e1;
        delete fh_HodTime_e2;
        delete fh_DCH_timediff_mu_e1;
        delete fh_DCH_timediff_mu_e2;
        delete fh_Hod_timediff_mu_e1;
        delete fh_DCH_timediff_e1_e2;
        delete fh_Hod_timediff_mu_e2;
        delete fh_Hod_timediff_e1_e2;
        delete fh_muee_P;
        delete fh_muee_Pt;
        delete fh_mee;
        delete fh_muee_M;
        delete fh_missing_mass;
        delete fh_mee_vs_Pmu  ;
        delete fh_mee_vs_Pel1 ;
        delete fh_mee_vs_Pel2 ;
        delete fh_mee_vs_Pt   ;
        delete fh_mee_vs_MM2  ;
        delete fh_mee_vs_zvtx ;
        delete fh_P_vs_dtrkcl ;
        delete fh_Ecl_vs_dtrkcl ;
        delete fh_HoDTotTime;
        delete fh_cda_mu_e1;
        delete fh_cda_mu_e2;
        delete fh_cda_e1_e2;
        delete fh_muon_bx     ;
        delete fh_muon_by     ;
        delete fh_electron1_bx;
        delete fh_electron1_by;
        delete fh_electron2_bx;
        delete fh_electron2_by;
        delete fh_bx_vs_by_muon;
        delete fh_bx_vs_by_el1;
        delete fh_bx_vs_by_el2;
        delete fh_el1_cluster_x;
        delete fh_el1_cluster_y;
        delete fh_el2_cluster_x;
        delete fh_el2_cluster_y;
        delete fh_mu_cluster_x;
        delete fh_mu_cluster_y;
        delete fh_el1_cluster_time;
        delete fh_el2_cluster_time;
        delete fh_el1_el2_cluster_timediff;
        delete fh_el1_dtrk_cl;
        delete fh_el2_dtrk_cl;
        delete fh_muon_dtrk_cl;
        //delete fh_el1_cluster_x_y;
        //delete fh_el2_cluster_x_y;
        //delete fh_mu_cluster_x_y;
        delete fh_deadcell_distance;
        delete fh_muv2_trk_cl_diff;
        delete fh_muv_x_y_position;
        delete fh_DCH1_distance_e1e2;
        delete fh_Lkr_distance_e1e2;
        delete fh_mc_KDzvtx;
        delete fh_mc_P1_Pzvtx;
        delete fh_mc_P1_mass;
        delete fh_mc_P1_momentum;
        delete fh_mc_three_track_123_momentum;
        delete fh_mc_P2_mass;
        delete fh_mc_P2_Pzvtx;
        delete fh_mc_P2_momentum;
        delete fh_mc_P3_Pzvtx;
        delete fh_mc_P3_mass;
        delete fh_mc_P3_momentum;
        delete fh_mc_three_track_123_mass;
        delete fh_mc_two_track_23_momentum;
        delete fh_mc_two_track_23_mass;
        delete fh_mc_P4_Pzvtx;
        delete fh_mc_P4_mass;
        delete fh_mc_P4_momentum;
        delete fh_mc_four_track_1234_momentum;
        delete fh_mc_four_track_1234_mass;
        delete fh_angle_el_mu;
        delete fh_mee_z_variable;
        delete fh_muv_xpos;
        delete fh_muv_ypos;
        delete fh_mc_P1_dist_prod_dec;
        delete fh_mc_P2_dist_prod_dec;
        delete fh_mc_P3_dist_prod_dec;
        delete fh_mc_two_track_23_mass_z_variable;
        delete fh_Muee_M_3pi_assumption;
        delete fh_Lkr_extrap_tracks_x_vs_y;
        delete fh_EoP_vs_p_odd_tr;
        delete fh_lda3_e1;
        delete fh_lda3_e2;
        delete fh_lda3_p1;
        delete fh_Mee_before;
        delete fh_lda3_p2;
        delete fh_lda3_p3;
        delete fh_muon_status;
        delete fh_MM2_minus;
        delete fh_MM2_plus;
    }
    if(ftype==3){
        delete fh_missing_mass          ;
        delete fh_mee;
        delete fh_mee_z_variable;
        delete fh_MM2_plus;
        delete fh_MM2_minus;
    }
}
