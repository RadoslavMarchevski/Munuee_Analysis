#include "Hist_dir.h"
#include "Charged_Particle.h"
#include "Cuts.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH1I.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <string>


void Hist_dir::FillHist(Charged_Particle& part, std::string particle){
    if(particle.compare("muon") == 0 ){
        //Particle specific histograms
        fh_muon_bx->Fill(part.Position[0]);
        fh_muon_by->Fill(part.Position[1]);
        fh_bx_vs_by_muon->Fill(part.Position[0],part.Position[1]);
        fh_DCHtime_mu->Fill( part.GetDCHtime());
        fh_Mu_momentum->Fill(part.GetMomentum());
        fh_Mu_charge->Fill(part.GetCharge() );
        fh_Mu_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
        fh_Track_Momentum->Fill(part.GetMomentum());
        fh_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );

    }
    if(particle.compare("electron1") == 0 ){
        fh_electron1_bx->Fill(part.Position[0]);
        fh_electron1_by->Fill(part.Position[1]);
        fh_bx_vs_by_el1->Fill(part.Position[0],part.Position[1]);
        fh_Track_Momentum->Fill(part.GetMomentum());
        fh_Electron_Momentum->Fill(part.GetMomentum());
        fh_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
        fh_Electron_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
        fh_DCHtime_e1->Fill( part.GetDCHtime());
    }
    if(particle.compare("electron2") == 0 ){
        fh_electron2_bx->Fill(part.Position[0]);
        fh_electron2_by->Fill(part.Position[1]);
        fh_bx_vs_by_el2->Fill(part.Position[0],part.Position[1]);
        fh_Track_Momentum->Fill(part.GetMomentum());
        fh_Electron_Momentum->Fill(part.GetMomentum());
        fh_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
        fh_Electron_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
        fh_DCHtime_e2->Fill( part.GetDCHtime());
    }
    if(particle.compare("pion1") == 0 ){
        //Particle specific histograms
        fh_pion1_bx->Fill(part.Position[0]);
        fh_pion1_by->Fill(part.Position[1]);
        fh_bx_vs_by_pi1->Fill(part.Position[0],part.Position[1]);
        fh_Pion_Momentum->Fill(part.GetMomentum());
        fh_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
    }

    if(particle.compare("pion2") == 0 ){
        //Particle specific histograms
        fh_pion2_bx->Fill(part.Position[0]);
        fh_pion2_by->Fill(part.Position[1]);
        fh_bx_vs_by_pi2->Fill(part.Position[0],part.Position[1]);
        fh_Pion_Momentum->Fill(part.GetMomentum());
        fh_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
    }

    if(particle.compare("pion3") == 0 ){
        //Particle specific histograms
        fh_pion3_bx->Fill(part.Position[0]);
        fh_pion3_by->Fill(part.Position[1]);
        fh_bx_vs_by_pi3->Fill(part.Position[0],part.Position[1]);
        fh_Pion_Momentum->Fill(part.GetMomentum());
        fh_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
    }
}

void Hist_dir::FillHist(Charged_Particle& p1,Charged_Particle& p2, std::string particles){


    if(particles.compare("mue1") == 0 ){
        fh_DCH_timediff_mu_e1->Fill( p1.GetDCHtime() - p2.GetDCHtime());
        fh_Hod_timediff_mu_e1->Fill( p1.GetHodTime() - p2.GetHodTime());
    }
    if(particles.compare("mue2") == 0 ){
        fh_DCH_timediff_mu_e2->Fill( p1.GetDCHtime() - p2.GetDCHtime());
        fh_Hod_timediff_mu_e2->Fill( p1.GetHodTime() - p2.GetHodTime());
    }
    if(particles.compare("e1e2") == 0 ){
        //TLorentzVector Two_Track_Momentum;
        Two_Track_Momentum = p1.Momentum + p2.Momentum;
        fh_mee->Fill( Two_Track_Momentum.M());
        fh_DCH_timediff_e1_e2->Fill(p1.GetDCHtime() -  p2.GetDCHtime());
        fh_Hod_timediff_e1_e2->Fill(p1.GetHodTime() - p2.GetHodTime());
    }


}
void Hist_dir::FillHist(TLorentzVector Three_Track_Momentum, TLorentzVector Nu_Momentum){

    //fh_Kaon_Charge->Fill(p1.GetCharge() + p2.GetCharge() + p3.GetCharge());
    fh_missing_mass->Fill(Nu_Momentum.M2());
    fh_muee_P->Fill(Three_Track_Momentum.P());
    fh_muee_Pt->Fill(Three_Track_Momentum.Pt());
    fh_muee_M->Fill(Three_Track_Momentum.M());

    //Make_Cuts();
}

void Hist_dir::FillCommonHist(superCmpEvent* sevt){
    int    ntrack = sevt->Ntrack; // number of tracks
    int    nclust = sevt->Ncluster; // number of clusters
    int    nvtx   = sevt->Nvtx;   // number of vtx
    double COmPaCt_Z_Vertex = sevt->vtx[0].z;

    fh_Ntracks->Fill(ntrack);
    fh_Nclusters->Fill(nclust);
    fh_Nvtx->Fill(nvtx);
    fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex);

}
void Hist_dir::FillVertexHist(double* mu_e1, double cda_mu_e1,double* mu_e2, double cda_mu_e2, double* e1_e2, double cda_e1_e2,std::string decay_type){
    if(decay_type.compare("munuee") ==0){
        fh_zvtxdiff_mue1_mue2->Fill( mu_e1[2] - mu_e2[2]);
        fh_zvtxdiff_mue1_e1e2->Fill( mu_e1[2] - e1_e2[2]);
        fh_zvtxdiff_mue2_e1e2->Fill( mu_e2[2] - e1_e2[2]);
        fh_yvtxdiff_mue1_mue2->Fill( mu_e1[1] - mu_e2[1]);
        fh_yvtxdiff_mue1_e1e2->Fill( mu_e1[1] - e1_e2[1]);
        fh_yvtxdiff_mue2_e1e2->Fill( mu_e2[1] - e1_e2[1]);
        fh_xvtxdiff_mue1_mue2->Fill( mu_e1[0] - mu_e2[0]);
        fh_xvtxdiff_mue1_e1e2->Fill( mu_e1[0] - e1_e2[0]);
        fh_xvtxdiff_mue2_e1e2->Fill( mu_e2[0] - e1_e2[0]);
        fh_cda_mu_e1->Fill(cda_mu_e1);
        fh_cda_mu_e2->Fill(cda_mu_e2);
        fh_cda_e1_e2->Fill(cda_e1_e2);
    }
    if(decay_type.compare("K3pi") == 0 ){
        fh_zvtxdiff_pi1pi2_pi2pi3->Fill(mu_e1[2] - mu_e2[2]);
        fh_zvtxdiff_pi1pi2_pi1pi3->Fill(mu_e1[2] - e1_e2[2]);
        fh_zvtxdiff_pi1pi3_pi2pi3->Fill(mu_e2[2] - e1_e2[2]);
        fh_xvtxdiff_pi1pi2_pi2pi3->Fill(mu_e1[1] - mu_e2[1]);
        fh_xvtxdiff_pi1pi2_pi1pi3->Fill(mu_e1[1] - e1_e2[1]);
        fh_xvtxdiff_pi1pi3_pi2pi3->Fill(mu_e2[1] - e1_e2[1]);
        fh_yvtxdiff_pi1pi2_pi2pi3->Fill(mu_e1[0] - mu_e2[0]);
        fh_yvtxdiff_pi1pi2_pi1pi3->Fill(mu_e1[0] - e1_e2[0]);
        fh_yvtxdiff_pi1pi3_pi2pi3->Fill(mu_e2[0] - e1_e2[0]);
        fh_cda_pi1_pi2           ->Fill(cda_mu_e1);
        fh_cda_pi1_pi3           ->Fill(cda_mu_e2);
        fh_cda_pi2_pi3           ->Fill(cda_e1_e2);

    }
}

void Hist_dir::ComputeThreeTrack(Charged_Particle& p1,Charged_Particle& p2, Charged_Particle& p3){
    static double massKaonC = 0.493677;
    //TLorentzVector Three_Track_Momentum;
    //TLorentzVector Kaon_Momentum;
    //TLorentzVector Nu_Momentum;
    Three_Track_Momentum = p1.Momentum + p2.Momentum + p3.Momentum;
    Kaon_Momentum.SetPxPyPzE(0,0,60.,TMath::Sqrt(60*60 + massKaonC*massKaonC));
    Nu_Momentum = Kaon_Momentum - Three_Track_Momentum;

}
