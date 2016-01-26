/* Charged_Particle.c
 *
 * Created on: 28.08.2015
 * User defined function for the class Neutral_Particle
 * defined in usersrc/Neutral_Particle.c
 *
 * Created by: Radoslav Marchevski
 */

#include "Neutral_Particle.h"
#include <iostream>
#include <math.h>


Neutral_Particle::Neutral_Particle(superCmpEvent* sevt, superBurst *sburst, int icluster){

    fSevt     = sevt;
    fsburst   = sburst;
    ficluster = icluster;
    CompactEnergy    = fSevt->cluster[ficluster].energy;
    CompactMomentum  = CompactEnergy;

    CalculateMomentum();

}

void Neutral_Particle::CalculateMomentum(){
    double Lkrz = Geom->Lkr.z;
    double Norm;

    Cluster_Position[0] = fSevt->cluster[ficluster].x;
    Cluster_Position[1] = fSevt->cluster[ficluster].y;
    Cluster_Position[2] = Lkrz + 16.5 + 4.3*CompactEnergy;

    Vtx_Position[0] = fSevt->vtx[0].x;
    Vtx_Position[1] = fSevt->vtx[0].y;
    Vtx_Position[2] = fSevt->vtx[0].z;

    Slopes[0] = (Vtx_Position[0] - Cluster_Position[0])/(Vtx_Position[2] - Cluster_Position[2]);
    Slopes[1] = (Vtx_Position[1] - Cluster_Position[1])/(Vtx_Position[2] - Cluster_Position[2]);
    Slopes[2] = 1.;

    Norm = 1./sqrt(1. + pow(Slopes[0],2) + pow(Slopes[1],2));

    Momentum.SetPx(Norm*CompactMomentum*Slopes[0]);
    Momentum.SetPy(Norm*CompactMomentum*Slopes[1]);
    Momentum.SetPz(Norm*CompactMomentum*Slopes[2]);
    Momentum.SetE(CompactMomentum);

}
