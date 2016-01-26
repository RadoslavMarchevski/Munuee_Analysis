/*
 *  Charged_track.h
 *       Created on 19.03.2015
 *       Author:Radoslav Marchevski
 */


#include <iostream>
#include <math.h>
#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include "compact.h"
#include <constants.h>

class Charged_track {
 public:
  Charged_track(superCmpEvent *sevt, superBurst *sburst, int PDGCode);
  void CorrectedMomentumAlphaBeta();
  void CalculateMomentum();


  double CompactMomentum;
  double Bdxdz;
  double Bdydz;
  double Mass;
  TVector3 Momentum3;
  TLorentzVector Momentum;
  
 private:
  superCmpEvent *evt;
  superBurst *burst;
  int PID;
};
