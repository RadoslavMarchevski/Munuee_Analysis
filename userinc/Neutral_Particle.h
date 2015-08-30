#ifndef __Neutral_Particle_h_
#define  __Neutral_Particle_h_
#include <math.h>
#include <TAxis.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TDirectory.h>
#include <TVector3.h>
#include "cmpio.h"
#include "user_NEW.h"
#include "reader.h"
#include <iostream>
#include <vector>
#include <math.h>


class Neutral_Particle {
public:
    Neutral_Particle(superCmpEvent* sevt, superBurst *sburst, int icluster);
    void CalculateMomentum();

    double CompactEnergy;
    double CompactMomentum;
    double Slopes[3];
    double Vtx_Position[3];
    double Cluster_Position[3];
    TLorentzVector Momentum;

private:
    superCmpEvent* fSevt;
    superBurst* fsburst;
    int ficluster;
};

#endif
