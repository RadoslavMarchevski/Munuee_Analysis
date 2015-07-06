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
    else {
        std::cout << "+++ ERROR +++ ERROR +++ NO OR WRONG CHARGED PARTICLE TYPE CHOOSEN +++ ERROR +++ ERROR" << std::endl;
    }

    EcalEnergy = fSevt->cluster[sevt->track[pindex].iClus].energy;

    CalculateMomentum();
    MakeTime();
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
    ClusterTime = fSevt->cluster[fSevt->track[fpindex].iClus].time;
}
