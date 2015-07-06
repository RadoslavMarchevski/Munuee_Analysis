#ifndef __Charged_Particle_h_
#define  __Charged_Particle_h_

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

class Charged_Particle {
public:
    Charged_Particle(superCmpEvent* sevt, superBurst *sburst, int PDGCode, int pindex);
    void CalculateMomentum();
    void MakeTime();
    //Methods for getting objects properties
    const int    GetCharge()      const {return fSevt->track[fpindex].q; };
    const double GetMass()        const {return fMass;}
    const double GetMomentum()    const {return Momentum.P();}
    const double GetEnergyLeftInEcal() const{return EcalEnergy;}
    const double GetDCHtime() const{return DCHtime;}
    const double GetHodTime() const{return HodTime;}
    const double GetClusterTime() const{return ClusterTime;}

    //Class Variables
    double CompactMomentum;
    TVector3 Momentum3NoAlphaNoBetaCorrection;
    TVector3 Momentum3NoNoBetaCorrection;
    TVector3 Momentum3;
    TLorentzVector Momentum;
    double EcalEnergy;
    double Slopes[3];
    double Position[3];
    double bxdz;
    double bydz;
    double DCHtime;
    double HodTime;
    double ClusterTime;
private:
    superCmpEvent* fSevt;
    superBurst* fsburst;
    double fMass;
    /* std::vector<double> fslopes; */
    /* std::vector<double> fposition; */

    double fMomentum_AB_Corr;
    int fPID;
    int fpindex;
};
#endif
