/********************************************************/
/* COmPACT user routine: user_exit()                    */
/*                                                      */
/* User routine called once all data has been read. This*/
/* allows output files to be closed and any other user  */
/* resources to be tidied up.                           */
/*                                          RWM 20/6/97 */
/********************************************************/

#include <math.h>
#include "cmpio.h"
#include "user_NEW.h"
//#include "user.h"
#include "reader.h"
#include "vscompact.h"
#include "Hist_dir.h"

int user_exit() {

/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/

    TFile* file = new TFile("output.root", "RECREATE");
    Initial_dir->AddToFile(file);
    K3pi_selection->AddToFile(file);
    dir1->AddToFile(file);
    dir3->AddToFile(file);
    dir4->AddToFile(file);
    dir7->AddToFile(file);
    dir10->AddToFile(file);
    dir11->AddToFile(file);
    dir2->AddToFile(file);
    dir5->AddToFile(file);
    dir6->AddToFile(file);
    dir9->AddToFile(file);
    dir20->AddToFile(file);
    dir21->AddToFile(file);
    dir22->AddToFile(file);
    dir23->AddToFile(file);
    dir19->AddToFile(file);
    dir12->AddToFile(file);
    dir13->AddToFile(file);
    //dir14->AddToFile(file);
    dir15->AddToFile(file);
    dir16->AddToFile(file);
    dir17->AddToFile(file);
    dir24->AddToFile(file);
    dir25->AddToFile(file);
    dir26->AddToFile(file);
    dir27->AddToFile(file);

    dir18->AddToFile(file);
    //MC_reweight->AddToFile(file);
    printf ("### saving histograms to file output.root \n");
    //delete dir1;
    //delete dir2;
    HistnoBursts->Write();
    HistnoGoodBursts->Write();
    file->Write();
    file->Close();

/*----------- End of user C code -----------*/
    return 0;
}
