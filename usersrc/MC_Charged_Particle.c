#include "Charged_Particle.h"
#include "MC_Charged_Particle.h"
#include "Hist_dir.h"

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


MC_Charged_Particle::MC_Charged_Particle(){}

MC_Charged_Particle::MC_Charged_Particle(superMcEvent* evt, superBurst* sburst, int pindex ){

    fEvt    = evt;
    fBurst  = sburst;
    fpindex = pindex;


}

void MC_Charged_Particle::Make4Momentum(){

    Get4Momentum.SetPx(fEvt->part[fpindex].p[1]);
    Get4Momentum.SetPy(fEvt->part[fpindex].p[2]);
    Get4Momentum.SetPz(fEvt->part[fpindex].p[3]);
    Get4Momentum.SetE(fEvt->part[fpindex].p[0]);

}
