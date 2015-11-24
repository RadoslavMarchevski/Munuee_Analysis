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
        if(part.cluster_exists){
            fh_mu_cluster_x->Fill(part.cluster_position[0]);
            fh_mu_cluster_y->Fill(part.cluster_position[1]);
            fh_muon_dtrk_cl->Fill(part.GetDistanceTrackCluster());
            //fh_mu_cluster_x_y->Fill(part.cluster_position[0],part.cluster_position[1]);
            fh_deadcell_distance->Fill(part.GetDistanceDeadcell());
        }
        fh_muv2_trk_cl_diff->Fill(part.MUV2_distance_trk_cl);
        fh_muv_xpos->Fill(part.extrapolated_track_MUV2[0]);
        fh_muv_ypos->Fill(part.extrapolated_track_MUV2[1]);
        fh_muv_x_y_position->Fill(part.extrapolated_track_MUV2[0],part.extrapolated_track_MUV2[1]);
        fh_Lkr_extrap_tracks_x_vs_y->Fill(part.extrapolated_track_Lkr[0],part.extrapolated_track_Lkr[1]);

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
        fh_el1_cluster_x->Fill(part.cluster_position[0]);
        fh_el1_cluster_y->Fill(part.cluster_position[1]);
        //fh_el1_cluster_x_y->Fill(part.cluster_position[0],part.cluster_position[1]);
        fh_el1_cluster_time->Fill(part.GetClusterTime());
        fh_el1_dtrk_cl->Fill(part.GetDistanceTrackCluster());
        fh_Lkr_extrap_tracks_x_vs_y->Fill(part.extrapolated_track_Lkr[0],part.extrapolated_track_Lkr[1]);

        fh_deadcell_distance->Fill(part.GetDistanceDeadcell());
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
        //fh_el2_cluster_x_y->Fill(part.cluster_position[0],part.cluster_position[1]);
        fh_el2_cluster_x->Fill(part.cluster_position[0]);
        fh_el2_cluster_y->Fill(part.cluster_position[1]);
        fh_el2_cluster_time->Fill(part.GetClusterTime());
        fh_el2_dtrk_cl->Fill(part.GetDistanceTrackCluster());
        fh_deadcell_distance->Fill(part.GetDistanceDeadcell());
        fh_Track_Momentum->Fill(part.GetMomentum());
        fh_Electron_Momentum->Fill(part.GetMomentum());
        fh_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
        fh_Electron_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
        fh_DCHtime_e2->Fill( part.GetDCHtime());
        fh_Lkr_extrap_tracks_x_vs_y->Fill(part.extrapolated_track_Lkr[0],part.extrapolated_track_Lkr[1]);
    }
    if(particle.compare("pion1") == 0 ){
        //Particle specific histograms
        fh_pion1_bx->Fill(part.Position[0]);
        fh_pion1_by->Fill(part.Position[1]);
        fh_bx_vs_by_pi1->Fill(part.Position[0],part.Position[1]);
        fh_Pion_Momentum->Fill(part.GetMomentum());
        fh_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
        fh_muv_x_y_position->Fill(part.extrapolated_track_MUV2[0],part.extrapolated_track_MUV2[1]);
        fh_Lkr_extrap_tracks_x_vs_y->Fill(part.extrapolated_track_Lkr[0],part.extrapolated_track_Lkr[1]);
    }

    if(particle.compare("pion2") == 0 ){
        //Particle specific histograms
        fh_pion2_bx->Fill(part.Position[0]);
        fh_pion2_by->Fill(part.Position[1]);
        fh_bx_vs_by_pi2->Fill(part.Position[0],part.Position[1]);
        fh_Pion_Momentum->Fill(part.GetMomentum());
        fh_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
        fh_muv_x_y_position->Fill(part.extrapolated_track_MUV2[0],part.extrapolated_track_MUV2[1]);
        fh_Lkr_extrap_tracks_x_vs_y->Fill(part.extrapolated_track_Lkr[0],part.extrapolated_track_Lkr[1]);

    }

    if(particle.compare("pion3") == 0 ){
        //Particle specific histograms
        fh_pion3_bx->Fill(part.Position[0]);
        fh_pion3_by->Fill(part.Position[1]);
        fh_bx_vs_by_pi3->Fill(part.Position[0],part.Position[1]);
        fh_Pion_Momentum->Fill(part.GetMomentum());
        fh_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
        fh_muv_x_y_position->Fill(part.extrapolated_track_MUV2[0],part.extrapolated_track_MUV2[1]);
        fh_Lkr_extrap_tracks_x_vs_y->Fill(part.extrapolated_track_Lkr[0],part.extrapolated_track_Lkr[1]);
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
        fh_el1_el2_cluster_timediff->Fill(p1.GetClusterTime() - p2.GetClusterTime());
        fh_DCH1_distance_e1e2->Fill( sqrt( pow(p1.Position[0] - p2.Position[0], 2 ) + pow(p1.Position[1] - p2.Position[1],2)));
        fh_Lkr_distance_e1e2->Fill( sqrt( pow(p1.extrapolated_track_Lkr[0] - p2.extrapolated_track_Lkr[0], 2 ) + pow(p1.extrapolated_track_Lkr[1] - p2.extrapolated_track_Lkr[1],2)));
    }


}
void Hist_dir::FillHist(TLorentzVector Three_Track_Momentum, TLorentzVector Nu_Momentum, int Kcharge){

    //fh_Kaon_Charge->Fill(p1.GetCharge() + p2.GetCharge() + p3.GetCharge());

    fh_missing_mass->Fill(Nu_Momentum.M2());

    fh_muee_P->Fill(Three_Track_Momentum.P());
    fh_muee_Pt->Fill(Three_Track_Momentum.Pt());
    fh_muee_M->Fill(Three_Track_Momentum.M());
    if(Kcharge == 1)
        fh_MM2_plus->Fill(Nu_Momentum.M2());
    if(Kcharge == -1)
        fh_MM2_minus->Fill(Nu_Momentum.M2());


    //Make_Cuts();
}

void Hist_dir::FillAngle(TLorentzVector muon, TLorentzVector Two_Track_Momentum){
    /* std::cout << muon.Vect().Angle(Two_Track_Momentum.Vect()) << "AND "; */
    //std::cout << "With E= " << Two_Track_Momentum.E() << "With p= " << Two_Track_Momentum.P() ;

    //TLorentzVector pKaon;
    //pKaon.SetPxPyPzE(0,0,60, 60*60 + 0.493677*0.493677);
    //muon.Boost(pKaon.BoostVector());
    //Two_Track_Momentum.Boost(pKaon.BoostVector());
    //std::cout << "AND With E= " << Two_Track_Momentum.E() << "With p= " << Two_Track_Momentum.P() << endl;
    //std::cout << muon.Vect().Angle(Two_Track_Momentum.Vect()) << endl;
    //muon.Boost(0,0,pKaon.Beta());
    //Two_Track_Momentum.Boost(0.,0.,pKaon.Beta());
    /* fh_angle_el_mu->Fill(muon.Vect().Angle(Two_Track_Momentum.Vect())); */
    //double cos_theta;
    // cos_theta = ( pow(MuNu_Momentum.M(),2) + pow(Two_Track_Momentum.M(),2) - pow(mKaon,2) + 2*MuNu_Momentum.E()*Two_Track_Momentum.M())/ (2*MuNu_Momentum.P()*Two_Track_Momentum.P()) ;
    fh_angle_el_mu->Fill(MuNu_Momentum.Angle(Two_Track_Momentum.Vect()));
    //fh_angle_el_mu->Fill(cos_theta);
    fh_mee_z_variable->Fill( Two_Track_Momentum.M()*Two_Track_Momentum.M() / (mKaon*mKaon));
//std::cout  << "Just muon and ee = " << muon.Angle(Two_Track_Momentum.Vect()) << std::endl;
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
    Nu_Momentum   = Kaon_Momentum - Three_Track_Momentum;
    MuNu_Momentum = Kaon_Momentum - Two_Track_Momentum;
    //std::cout  << "Using combined vectors" << MuNu_Momentum.Angle(Two_Track_Momentum.Vect()) << std::endl;
}

void Hist_dir::Fill3pi(TLorentzVector Three_Track_Momentum ){
    fh_Muee_M_3pi_assumption->Fill(Three_Track_Momentum.M());
}
