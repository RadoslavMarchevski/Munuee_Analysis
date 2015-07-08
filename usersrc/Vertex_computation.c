#include "Vertex_computation.h"
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

Charged_Vertex::Charged_Vertex(){
    //std::cout << muon.GetMomentum() << std::endl;
    //cda = 0;
    //vertex[3]= {0.};
    //fp1 = particle1;
    //fp2 = particle2;
    //closap_double_(fp1.Position,fp2.Position,fp1.Slopes,fp2.Slopes,&cda,vertex);
}

~Charged_Vertex();
