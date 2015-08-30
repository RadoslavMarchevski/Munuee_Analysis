#include <iostream>
#include <math.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TVector.h>
#include <TVector3.h>
#include "Charged_track.h"
using namespace std;


Charged_track::Charged_track(superCmpEvent *sevt, superBurst *sburst, int PDGCode){
  evt = sevt;
  burst = sburst;
  PID = PDGCode;
  CompactMomentum = (double) evt->track[0].p;
  if(PID == -13)/*Mu */{
    Mass = 0.105658369; // GeV
  }
  else if(PID == -11)/*Electron*/{
    Mass = 0.000510998918; //GeV
  }
  else if(PID == 211)/*Pi+*/{
    Mass = 0.13957018;
  }
  else {
    std::cout << "+++ ERROR +++ ERROR +++ NO OR WRONG CHARGED PARTICLE TYPE CHOOSEN +++ ERROR +++ ERROR" << std::endl;
  }
  
  CorrectedMomentumAlphaBeta();
  CalculateMomentum3();
  CalculateMomentum();
}

void CorrectedMomentumAlphaBeta(){

  double MomentumAlphaBetaCorr = (double) p_corr_ab(CompactMomentum,evt->track[0].q);

  Bdxdz = (double) evt->track[0].bdxdz;
  Bdydz = (double) evt->track[0].bdydz;

  Momentum3.SetX( (1./TMath::Sqrt( pow(Bdxdz,2) + pow(Bdydz,2) + 1 ) ) * Bdxdz * MomentumAlphaBetaCorr);
  Momentum3.SetY( (1./TMath::Sqrt( pow(Bdxdz,2) + pow(Bdydz,2) + 1 ) ) * Bdydz * MomentumAlphaBetaCorr);
  Momentum3.SetY( (1./TMath::Sqrt( pow(Bdxdz,2) + pow(Bdydz,2) + 1 ) ) * MomentumAlphaBetaCorr);
}

  
void Charged_track::CalculateMomentum(){
    
  Momentum.SetPx(Momentum3.X());
  Momentum.SetPy(Momentum3.Y());
  Momentum.SetPz(Momentum3.Z());
  Momentum.SetE(TMath::Sqrt( Mass*Mass + Momentum.Vect()*Momentum.Vect()));
}
