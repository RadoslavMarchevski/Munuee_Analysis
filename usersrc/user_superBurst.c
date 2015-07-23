/***********************************************************/
/* COmPACT user routine: user_superBurst(superBurst *sbur) */
/*                                                         */
/* User routine called everytime a SuperCOmPACT burst      */
/* `*bur' is loaded.                                       */
/* A return value of greater than zero denotes an error    */
/* condition has occured.                                  */
/*                                             BH 2/3/98   */
/***********************************************************/

#include "cmpio.h"
//#include "user.h"
#include "user_NEW.h"
#include "reader.h"

int user_superBurst(superBurst *sbur) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
    sbur->BadB.Skip = 0; /* see user_superBurst.example.c to learn to use it */

<<<<<<< HEAD
    noBursts++;            // total number of bursts
    burstCounter++;        // counter set to 0 in superEob if certain number of bursts reached to write the root file
=======
    //noBursts++;            // total number of bursts
    //burstCounter++;        // counter set to 0 in superEob if certain number of bursts reached to write the root file
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
    magnetCurrent =  sbur->MNP33current;
    runNo = sbur->nrun;

    if (runNo==20486)          // exclude run 20486 in P5
        sbur->BadB.Skip = 1;
    if (sbur->BadB.Phys!=0)    // exclude bad physics bursts
        sbur->BadB.Skip = 1;
    if (sbur->BadB.Dch!=0)    // exclude bad DCH bursts
        sbur->BadB.Skip = 1;


    //if (sbur->time != 1190177600)
    //sbur->BadB.Skip = 1;

    // Ab v65: Determine in which period we are and set flag
    //
    if (runNo >= 20114 && runNo <= 20203)       // we are in period 1
        periodFlag = 1;
    else if (runNo >= 20204 && runNo <= 20238)  // we are in period 2, runs 20209-20238
        periodFlag = 21;
    else if (runNo >= 20239 && runNo <= 20285)  // we are in period 2, runs 20254-20285
        periodFlag = 22;
    else if (runNo >= 20286 && runNo <= 20291)  // we are in period 3, runs 20286-20291
        periodFlag = 31;
    else if (runNo >= 20292 && runNo <= 20324)  // we are in period 3, runs 20296-20324
        periodFlag = 32;
    else if (runNo >= 20325 && runNo <= 20404)  // we are in period 4
        periodFlag = 4;
    else if (runNo >= 20410 && runNo <= 20415)  // we are in period 5, runs 20410-20415
        periodFlag = 51;
    else if (runNo >= 20416 && runNo <= 20485)  // we are in period 5, runs 20416-20485
        periodFlag = 52;
    else if (runNo >= 20487 && runNo <= 20531)  // we are in period 6
        periodFlag = 6;
    else if (runNo >= 20613 && runNo <= 20695)  // we are in period 7
        periodFlag = 7;
    else if (runNo >= 21082 && runNo <= 21103)  // we are in period 8, runs 21082-21103
        periodFlag = 81;
    else if (runNo >= 21106 && runNo <= 21120)  // we are in period 8, runs 21106-21120
        periodFlag = 82;

    // Determine Data type
    IS_DATA = false;
    IS_MC = false;
    if (sbur->brtype == 1)
        IS_DATA = true;
    else if (sbur->brtype == 2)
        IS_MC = true;
    //printf ("IS_DATA = %i \t IS_MC = %i\n", IS_DATA, IS_MC);

    //  printf ("burst number = %d  -  run %d  -  BadB.Phys = %i\n", noBursts, runNo, sbur->BadB.Phys);
    printf ("burst number = %d  -  run %d  -  period %d  -  burst time %d  -  BadB.Skip = %i  -  magnet current = %.1f  -  dataType = %i\n",
            noBursts, runNo, periodFlag, sbur->time, sbur->BadB.Skip, (float)(magnetCurrent/1000.), sbur->brtype);
    //printf ("magnet current = %d\n", magnetCurrent);

    // HistnoBursts->Fill(0);
    //if (sbur->BadB.Skip == 0)
    //HistnoGoodBursts->Fill(0);

    //  printf ("n2track_total = %i \t n2track_total_cut05 = %i \t nKe3_total = %i \t nKe3_total_psum60_120 = %i\n",
    //  n2track_total, n2track_total_cut05, nKe3_total, nKe3_total_psum60_120);

/*----------- End of user C code -----------*/
    return 0;
}
