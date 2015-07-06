#include "Hist_dir.h"
#include "Charged_Particle.h"
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
    //Creating directories with different histograms depending
    //on the directory type:
    //0 - Initial
    //1 - K3pi selection
    //2 - Kmunuee selection


    if(ftype==0){
        fh_Ntracks           = new TH1I("Ntracks","Number of charged tracks;Ntrk;Nevents",10,0,10);
        fh_Nvtx              = new TH1I("Nvtx","Number of vertexes;Nvtx;Nevents",10,0,10);
        fh_Track_Momentum    = new TH1F("Track_Momentum","Track momentum after Nvtx,Ntrack and ZVtx cut;Track_P[GeV];Nevents ",100.,0.,100.);
        fh_eop               = new TH1F("Track_E/p","Track E/p after Nvtx,Ntrack and ZVtx cut;E/p;Nevents",120.,0.,1.2);
        fh_COmPaCt_Z_Vertex  = new TH1F("COmPaCt_Z_Vertex","Three track vertex from COmPaCt",150,-4000.,11000.);
    }
    if(ftype==1){
        fh_Ntracks           = new TH1I("Ntracks","Number of charged tracks;Ntrk;Nevents",10,0,10);
        fh_Nvtx              = new TH1I("Nvtx","Number of vertexes;Nvtx;Nevents",10,0,10);
        fh_Kaon_Charge       = new TH1I("Kaon_Charge","Charge of the Kaon;Q;Nevents",4,-2,2);
        fh_Event_Type        = new TH1I("Event_Type","Mu e e : 0- ++-, 1- +++, 2- +--, 3- -++, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_Pion_Momentum     = new TH1F("Pion_Momentum","Pion Momentum for the K3pi selection;#pi_P[GeV];Nevents",100,0.,100.);
        fh_eop               = new TH1F("Track_E/p","Track E/p ;E/p;Nevents",120.,0.,1.2);
        fh_odd_eop           = new TH1F("odd_E/p","E/p for the odd track in the k3pi selection",120,0.,1.2);
        fh_COmPaCt_Z_Vertex  = new TH1F("COmPaCt_Z_Vertex","Three track vertex from COmPaCt",150,-4000.,11000.);
        fh_Z_Vertex          = new TH1F("Z_vertex","Reconstructed Z Vtx position ;Z_Vtx[cm];Nevents ",1500.,-5000.,10000.);

    }
    if(ftype==2){
        fh_Ntracks           = new TH1I("Ntracks","Number of charged tracks;Ntrk;Nevents",10,0,10);
        fh_Nvtx              = new TH1I("Nvtx","Number of vertexes;Nvtx;Nevents",10,0,10);
        fh_Kaon_Charge       = new TH1I("Kaon_Charge","Charge of the Kaon;Q;Nevents",4,-2,2);
        fh_Mu_charge         = new TH1I("Mu_charge","Muon charge ;Mu_Q;Nevents",4.,-2.,2.);
        fh_Event_Type        = new TH1I("Event_Type","Mu e e : 0- ++-, 1- +++, 2- +--, 3- -++, 4- -++, 5- ---;Event Type;Nevents",7,-1,6);
        fh_Track_Momentum    = new TH1F("Track_Momentum","Track momentum ;Track_P[GeV];Nevents ",100.,0.,100.);
        fh_Mu_momentum       = new TH1F("Mu_momentum","Muon momentum ;Mu_P[GeV];Nevents ",100.,0.,100.);
        fh_Electron_Momentum = new TH1F("Electron_Momentum", "Electron momentum ;El_P[GeV];Nevents",100,0.,100.);
        fh_eop               = new TH1F("Track_E/p","Track E/p ;E/p;Nevents",120.,0.,1.2);
        fh_Mu_eop            = new TH1F("Mu_E/p","Muon E/p ;E/p;Nevents",120.,0.,1.2);
        fh_Electron_eop      = new TH1F("Electron_E/p", "Electron E/p ;E/p;Nevents",120,0.,1.2);
        fh_odd_eop           = new TH1F("odd_E/p","E/p for the odd track in the k3pi selection",120,0.,1.2);
        fh_COmPaCt_Z_Vertex  = new TH1F("COmPaCt_Z_Vertex","Three track vertex from COmPaCt",150,-4000.,11000.);
        fh_Z_Vertex          = new TH1F("Z_vertex","Reconstructed Z Vtx position ;Z_Vtx[cm];Nevents ",150.,-5000.,10000.);
        fh_Mu_Zvtx_min_COmPaCt_Zvtx = new TH1F("Mu_Zvtx_min_COmPaCt_Zvtx","Calculated #mu Zvtx - COmPact three track Zvtx;(Z_#mu - Z_vtx)[cm];Nevents ",100,-10000.,10000.);
        //Reconstructed Vtx difference
        fh_xvtxdiff_mue1_mue2 = new TH1F("xvtxdiff_mue1_mue2","Difference between #mu e1 and #mu e2 reconstructed X vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_xvtxdiff_mue1_e1e2 = new TH1F("xvtxdiff_mue1_e1e2","Difference between #mu e1 and e1 e2 reconstructed X vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_xvtxdiff_mue2_e1e2 = new TH1F("xvtxdiff_mue2_e1e2","Difference between #mu e2 and e1 e2 reconstructed X vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_mue1_mue2 = new TH1F("yvtxdiff_mue1_mue2","Difference between #mu e1 and #mu e2 reconstructed Y vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_mue1_e1e2 = new TH1F("yvtxdiff_mue1_e1e2","Difference between #mu e1 and e1 e2 reconstructed  Y vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
        fh_yvtxdiff_mue2_e1e2 = new TH1F("yvtxdiff_mue2_e1e2","Difference between #mu e2 and e1 e2 reconstructed  Y vtx;Z_Vtx_diff[cm];Nevents ",400.,-2000.,2000.);
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
        fh_mee = new TH1F("mee","Invariant mass of the electron pair ",500,0.,0.5);
        fh_muee_P = new TH1F("muee_P","Three track momentum ;P_{#mu e e}[GeV];Nevents",100,0.,100);
        fh_muee_Pt = new TH1F("muee_Pt","Transverse momentum of #mu^{#pm} e^{+}e^{-} system;Pt_{#mu e e}[GeV];Nevents",200,0,2.);
        fh_muee_M = new TH1F("muee_M","Three track invariant mass ;M_{#mu e e}[GeV];Nevents",1000,-1.,1.);

        fh_missing_mass = new TH1F("missing_mass","Missing mass squared;M^{2}_{miss};Nevents",100,-0.05,0.05);

    }
}


void Hist_dir::AddToFile(TFile* file){
    TDirectory* dir = file->mkdir(fDir.c_str());
    dir->cd();
    if(ftype==0){
        fh_Ntracks           ->Write();
        fh_Nvtx              ->Write();
        fh_Track_Momentum    ->Write();
        fh_eop               ->Write();
        fh_COmPaCt_Z_Vertex  ->Write();
    }
    if(ftype==1){
        fh_Ntracks           ->Write();
        fh_Nvtx              ->Write();
        fh_Kaon_Charge       ->Write();
        fh_Event_Type        ->Write();
        fh_Pion_Momentum     ->Write();
        fh_eop               ->Write();
        fh_odd_eop           ->Write();
        fh_COmPaCt_Z_Vertex  ->Write();
        fh_Z_Vertex          ->Write();
    }
    if(ftype==2){
        fh_Ntracks           ->Write();
        fh_Nvtx              ->Write();
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
        fh_Z_Vertex          ->Write();
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
        fh_DCHtime_mu->Write();
        fh_DCHtime_mu->Write();
        fh_DCHtime_mu->Write();
        fh_HodTime_mu->Write();
        fh_DCHtime_e1->Write();
        fh_DCHtime_e2->Write();
        fh_HodTime_e1->Write();
        fh_HodTime_e2->Write();
        fh_DCH_timediff_mu_e1->Write();
        fh_DCH_timediff_mu_e2->Write();
        fh_Hod_timediff_mu_e1->Write();
        fh_DCH_timediff_e1_e2->Write();
        fh_Hod_timediff_mu_e2->Write();
        fh_Hod_timediff_e1_e2->Write();
        fh_mee->Write();
        fh_muee_P->Write();
        fh_muee_M->Write();
        fh_missing_mass->Write();
        fh_muee_Pt->Write();
    }
}

Hist_dir::~Hist_dir(){
    if(ftype==0){
        delete fh_Ntracks           ;
        delete fh_Nvtx              ;
        delete fh_Track_Momentum    ;
        delete fh_eop               ;
        delete fh_COmPaCt_Z_Vertex  ;
    }
    if(ftype==1){
        delete fh_Ntracks           ;
        delete fh_Nvtx              ;
        delete fh_Kaon_Charge       ;
        delete fh_Event_Type        ;
        delete fh_Pion_Momentum     ;
        delete fh_eop               ;
        delete fh_odd_eop           ;
        delete fh_COmPaCt_Z_Vertex  ;
        delete fh_Z_Vertex          ;
    }
    if(ftype==2){
        delete fh_Ntracks           ;
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
        delete fh_Z_Vertex          ;
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
    }
}


void Hist_dir::FillHist(Charged_Particle& part){
    fh_muee_M->Fill(part.GetMass());
    return;

}
