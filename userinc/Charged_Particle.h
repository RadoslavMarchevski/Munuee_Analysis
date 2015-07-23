#ifndef __Charged_Particle_h_
#define  __Charged_Particle_h_
<<<<<<< HEAD
#include "MC_Charged_Particle.h"
=======

>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
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
<<<<<<< HEAD
    void MakeLkrCluster();
    void MakeMUVExtrapolation();
=======
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
    //Methods for getting objects properties
    const int    GetCharge()      const {return fSevt->track[fpindex].q; };
    const double GetMass()        const {return fMass;}
    const double GetMomentum()    const {return Momentum.P();}
    const double GetEnergyLeftInEcal() const{return EcalEnergy;}
    const double GetDCHradius() const{return TMath::Sqrt(pow(Position[0],2) + pow(Position[1],2));}
    const double GetDCHtime() const{return DCHtime;}
    const double GetHodTime() const{return HodTime;}
<<<<<<< HEAD
    double* GetClusterPosition() {return cluster_position;}
    const double GetClusterIndex() const{return icluster;}
    const double GetClusterTime() const{return ClusterTime;}
    const double GetDistanceTrackCluster() const{return distance_trk_cl;}
    const double GetDistanceDeadcell() const{return deadcell_distance;}

    //Class Variables
    int icluster;
    int imuon;
    bool cluster_exists;
=======
    const double GetClusterTime() const{return ClusterTime;}

    //Class Variables
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
    double CompactMomentum;
    TVector3 Momentum3NoAlphaNoBetaCorrection;
    TVector3 Momentum3NoNoBetaCorrection;
    TVector3 Momentum3;
    TLorentzVector Momentum;
    double EcalEnergy;
<<<<<<< HEAD
    double Slopes[3];
    double Position[3];
    double cluster_position[3];
    double MUV2_position[3];
    double extrapolated_track_Lkr[2];
    double extrapolated_track_MUV2[2];
    double deadcell_distance;
    double bxdz;
    double bydz;
    double dxdz;
    double dydz;
    double distance_trk_cl;
    double MUV2_distance_trk_cl;
=======
    //double R_DCH;
    double Slopes[3];
    double Position[3];
    double bxdz;
    double bydz;
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
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
