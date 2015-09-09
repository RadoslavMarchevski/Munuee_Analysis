/*****************************************************************/
/* COmPACT user routine: user_superCmpEvent(superCmpEvent *sevt) */
/*                                                               */
/* User routine called everytime an event `*sevt' is             */
/* loaded. A return value of greater than zero denotes           */
/* an error condition has occured.                               */
/*                                     BH 13/2/98    RWM 20/6/97 */
/*****************************************************************/

#include <math.h>
#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include <constants.h>
// Datei iostream einbinden
#include <iostream>

// using Anweisung
using std::cout;
using std::endl;

int user_superCmpEvent(superBurst *sbur, superCmpEvent *sevt) {
	/* WARNING: do not alter things before this line */
	/*---------- Add user C code here ----------*/
	static int nuserevt = 0;
	int a, i, j, k;
	int l, m, n;

	if (nuserevt < 20) {
		printSuperCmpEvent(sevt, fprt);
	}

	nuserevt++;


	/*----------- End of user C code -----------*/
	return 0;
}