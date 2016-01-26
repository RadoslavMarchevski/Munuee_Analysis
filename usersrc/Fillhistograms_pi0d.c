#include "Hist_dir_pi0d.h"
#include "Charged_Particle.h"
#include "Cuts_pi0d.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH1I.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <string>


void Hist_dir::FillHist(Charged_Particle& part, std::string particle){
    //if(particle.compare("pion") == 0 )
    if(particle.compare("electron1") == 0 ){
        //fh_electron1_bx->Fill(part.Position[0]);
        //fh_electron1_by->Fill(part.Position[1]);
        //fh_el1_cluster_x->Fill(part.cluster_position[0]);
        //fh_el1_cluster_y->Fill(part.cluster_position[1]);
        //fh_el1_cluster_x_y->Fill(part.cluster_position[0],part.cluster_position[1]);
        //fh_el1_cluster_time->Fill(part.GetClusterTime());
        //fh_el1_dtrk_cl->Fill(part.GetDistanceTrackCluster());
        fh_Lkr_extrap_x_vs_y->Fill(part.extrapolated_track_Lkr[0],part.extrapolated_track_Lkr[1]);
        //
        //fh_deadcell_distance->Fill(part.GetDistanceDeadcell());
        //fh_bx_vs_by_el1->Fill(part.Position[0],part.Position[1]);
        //fh_Track_Momentum->Fill(part.GetMomentum());
        //fh_Electron_Momentum->Fill(part.GetMomentum());
        //fh_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
        //fh_Electron_eop->Fill(part.GetEnergyLeftInEcal()/ part.GetMomentum() );
        //fh_DCHtime_e1->Fill( part.GetDCHtime());
    }
    if(particle.compare("electron2") == 0 ){
        //fh_electron2_bx->Fill(part.Position[0]);
        //fh_electron2_by->Fill(part.Position[1]);
        //fh_bx_vs_by_el2->Fill(part.Position[0],part.Position[1]);
        //fh_el2_cluster_x_y->Fill(part.cluster_position[0],part.cluster_position[1]);
        //fh_el2_cluster_x->Fill(part.cluster_position[0]);
        //fh_el2_cluster_y->Fill(part.cluster_position[1]);
        //fh_el2_cluster_time->Fill(part.GetClusterTime());
        //fh_el2_dtrk_cl->Fill(part.GetDistanceTrackCluster());
        //fh_deadcell_distance->Fill(part.GetDistanceDeadcell());
        fh_Lkr_extrap_x_vs_y->Fill(part.extrapolated_track_Lkr[0],part.extrapolated_track_Lkr[1]);
    }

}

void Hist_dir::FillHist(Charged_Particle& p1,Charged_Particle& p2, std::string particles){


    if(particles.compare("pie1") == 0 ){
        fh_DCH_timediff_pi_e1->Fill( p1.GetDCHtime() - p2.GetDCHtime());
        fh_Hod_timediff_pi_e1->Fill( p1.GetHodTime() - p2.GetHodTime());
    }
    if(particles.compare("pie2") == 0 ){
        fh_DCH_timediff_pi_e2->Fill( p1.GetDCHtime() - p2.GetDCHtime());
        fh_Hod_timediff_pi_e2->Fill( p1.GetHodTime() - p2.GetHodTime());
    }
    if(particles.compare("e1e2") == 0 ){
        fh_DCH_timediff_e1_e2->Fill(p1.GetDCHtime() -  p2.GetDCHtime());
        fh_Hod_timediff_e1_e2->Fill(p1.GetHodTime() - p2.GetHodTime());
        fh_el1_el2_cltimediff->Fill(p1.GetClusterTime() - p2.GetClusterTime());
        fh_DCH1_distance_e1e2->Fill( sqrt( pow(p1.Position[0] - p2.Position[0], 2 ) + pow(p1.Position[1] - p2.Position[1],2)));
        fh_Lkr_distance_e1e2->Fill( sqrt( pow(p1.extrapolated_track_Lkr[0] - p2.extrapolated_track_Lkr[0], 2 ) + pow(p1.extrapolated_track_Lkr[1] - p2.extrapolated_track_Lkr[1],2)));
    }


}

void Hist_dir::FillVertexHist(double* pi_e1, double cda_pi_e1,double* pi_e2, double cda_pi_e2, double* e1_e2, double cda_e1_e2,std::string decay_type){
    if(decay_type.compare("k2pid") ==0){
        fh_zvtxdiff_pie1_pie2->Fill( pi_e1[2] - pi_e2[2]);
        fh_zvtxdiff_pie1_e1e2->Fill( pi_e1[2] - e1_e2[2]);
        fh_zvtxdiff_pie2_e1e2->Fill( pi_e2[2] - e1_e2[2]);
        fh_yvtxdiff_pie1_pie2->Fill( pi_e1[1] - pi_e2[1]);
        fh_yvtxdiff_pie1_e1e2->Fill( pi_e1[1] - e1_e2[1]);
        fh_yvtxdiff_pie2_e1e2->Fill( pi_e2[1] - e1_e2[1]);
        fh_xvtxdiff_pie1_pie2->Fill( pi_e1[0] - pi_e2[0]);
        fh_xvtxdiff_pie1_e1e2->Fill( pi_e1[0] - e1_e2[0]);
        fh_xvtxdiff_pie2_e1e2->Fill( pi_e2[0] - e1_e2[0]);
        fh_cda_pi_e1->Fill(cda_pi_e1);
        fh_cda_pi_e2->Fill(cda_pi_e2);
        fh_cda_e1_e2->Fill(cda_e1_e2);
    }
}
