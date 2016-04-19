/* Charged_Particle.c
 *
 * Created on: 16.06.2015
 * User defined function for the class Charged_Particle
 * defined in usersrc/Charged_Particle.c
 *
 * Created by: Radoslav Marchevski
 */

#include "Charged_Particle.h"
#include <iostream>
#include <math.h>

Charged_Particle::Charged_Particle(superCmpEvent* sevt, superBurst *sburst, int PDGCode, int pindex){
  fSevt = sevt;
  fsburst = sburst;
  fpindex = pindex;
  fPID = PDGCode;

  distance_trk_cl = 0;
  CompactMomentum = (double)sevt->track[pindex].p;
  if(fPID == -13)/*Muon+ choosen*/{
    fMass = 0.105658369; // GeV
  }
  else if(fPID == -11)/*Positron*/{
    fMass = 0.000510998918;
  }
  else if(fPID == 211)/*Pi+*/{
    fMass = 0.13957018;
  }
  else if(fPID == 321)/*K+*/{
    fMass = 0.493677;
  }
  else {
    std::cout << "+++ ERROR +++ ERROR +++ NO OR WRONG CHARGED PARTICLE TYPE CHOOSEN +++ ERROR +++ ERROR" << std::endl;
  }

  icluster = -1;
  imuon    = -1;
  cluster_exists   = false;

  icluster= fSevt->track[fpindex].iClus;
  imuon   = fSevt->track[fpindex].iMuon;
  EcalEnergy = fSevt->cluster[icluster].energy;

  CalculateMomentum();
  MakeTime();

  if(icluster != -1){
    cluster_exists = true;
    MakeLkrCluster(cluster_exists);
  } else {
    cluster_exists = false;
    MakeLkrCluster(cluster_exists);

  }
  if(imuon != -1){
    muvcl_exits = true;
    MakeMUVExtrapolation(muvcl_exits);
  } else {
    muvcl_exits = false;
    MakeMUVExtrapolation(muvcl_exits);

  }
}



void Charged_Particle::CalculateMomentum(){
  //Call routine for vertex
  //Alpha and beta corrections
  fMomentum_AB_Corr = p_corr_ab (CompactMomentum, fSevt->track[fpindex].q);

  //Getting the correct positions and slopes for the three track reconstruction
  Position[0]= fSevt->track[fpindex].bx;
  Position[1]= fSevt->track[fpindex].by;
  Position[2]= Geom->DCH.bz;

  Slopes[0]= fSevt->vtx[0].vtxtrack[fpindex].bdxdz;
  Slopes[1]= fSevt->vtx[0].vtxtrack[fpindex].bdydz;
  Slopes[2]= 1;

  bxdz = Slopes[0];
  bydz = Slopes[1];
  //ChargedVertexPosVectorNonBFCorr.SetXYZ(Vertex[0],Vertex[1],Vertex[2]);


  //Calculating the correct momentum for the track
  Momentum3.SetX( (1.0 / TMath::Sqrt( pow(bxdz,2) + pow(bydz,2) +1 )) * bxdz * fMomentum_AB_Corr);
  Momentum3.SetY( (1.0 / TMath::Sqrt( pow(bxdz,2) + pow(bydz,2) +1 )) * bydz * fMomentum_AB_Corr);
  Momentum3.SetZ( (1.0 / TMath::Sqrt( pow(bxdz,2) + pow(bydz,2) +1 )) * fMomentum_AB_Corr);

  Momentum.SetPx(Momentum3.X());
  Momentum.SetPy(Momentum3.Y());
  Momentum.SetPz(Momentum3.Z());
  Momentum.SetE( TMath::Sqrt( fMass*fMass + Momentum.Vect()*Momentum.Vect()));

}

void Charged_Particle::MakeTime(){
  DCHtime = fSevt->track[fpindex].time;
  HodTime = fSevt->track[fpindex].hodTime;
  ClusterTime = fSevt->cluster[icluster].time;
}

void Charged_Particle::MakeLkrCluster(bool cluster){
  double DCHz = Geom->DCH.z;
  double Lkrz = Geom->Lkr.z;
  dxdz = fSevt->track[fpindex].dxdz;
  dydz = fSevt->track[fpindex].dydz;

  //cluster_position[0] = fSevt->cluster[fpindex].x;
  //cluster_position[1] = fSevt->cluster[fpindex].y;
  cluster_position[0] = fSevt->cluster[icluster].x;
  cluster_position[1] = fSevt->cluster[icluster].y;
  if(cluster) {cluster_position[2] = Lkrz  + 16.5 + 4.3*EcalEnergy;}
  else {cluster_position[2] = Lkrz;}



  extrapolated_track_Lkr[0] = fSevt->track[fpindex].x + (cluster_position[2] - DCHz) * dxdz;
  extrapolated_track_Lkr[1] = fSevt->track[fpindex].y + (cluster_position[2] - DCHz) * dydz;

  distance_trk_cl = sqrt( pow(cluster_position[0] - extrapolated_track_Lkr[0],2) + pow(cluster_position[1] - extrapolated_track_Lkr[1],2) ) ;
  deadcell_distance = fSevt->track[fpindex].dDeadCell;

}

void Charged_Particle::MakeMUVExtrapolation(bool cluster){
  double DCHz = Geom->DCH.z;
  double MUV2z = Geom->Muv2.z;

  dxdz = fSevt->track[fpindex].dxdz;
  dydz = fSevt->track[fpindex].dydz;

  if(cluster){
    MUV2_position[0] = fSevt->muon[imuon].x;
    MUV2_position[1] = fSevt->muon[imuon].y;
    MUV2_position[2] = MUV2z;
  }else {

    MUV2_position[0] = 0;
    MUV2_position[1] = 0;
    MUV2_position[2] = MUV2z;
  }

  extrapolated_track_MUV2[0] = fSevt->track[fpindex].x + (MUV2_position[2] - DCHz) * dxdz;
  extrapolated_track_MUV2[1] = fSevt->track[fpindex].y + (MUV2_position[2] - DCHz) * dydz;

  MUV2_distance_trk_cl = sqrt( pow(MUV2_position[0] - extrapolated_track_MUV2[0],2) + pow(MUV2_position[1] - extrapolated_track_MUV2[1],2) ) ;
}
