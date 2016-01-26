#ifndef __MC_Charged_Particle_h_
#define  __MC_Charged_Particle_h_

//#include "Charged_Particle.h"
//#include "Hist_dir.h"

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

class MC_Charged_Particle{
public:
    MC_Charged_Particle();
    MC_Charged_Particle(superMcEvent* evt, superBurst* sburst, int pindex);
    const int GetType() const{ return fEvt->part[fpindex].type; };
    int GetI() { return fpindex;};
    double GetProductionZVertex(){return fEvt->part[fpindex].pvertex[2]; };
    double GetDecayZVertex(){return fEvt->part[fpindex].dvertex[2]; };
    void Make4Momentum();
    TLorentzVector Get4Momentum;
private:
    superMcEvent* fEvt;
    superBurst* fBurst;
    //int fNpart;
    int fpindex;
};
#endif
