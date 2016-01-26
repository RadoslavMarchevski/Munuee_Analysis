#include "Hist_dir.h"
#include "Charged_Particle.h"
#include "user_NEW.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH1I.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <string>


class Charged_Vertex{
public:
    Charged_Vertex();
    ~Charged_Vertex();
    double momentum;
    //TVector3& GetChargedVertex() {};

private:
}
