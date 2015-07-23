#include "stdio.h"
#include <math.h>
#include "user.h"

void GetCpdCellIndex(float pos_x, float pos_y, int *cpd_index, int *cell_index)
{
  int   i, j, k;
  int   l, m, n;
  float CELLlength = 1.975;
  float CPDlength = 8 * CELLlength;
  
  *cpd_index  = -1;
  *cell_index = -1;
  
  // Define CPD index	  
  for (i=0; i<16; i++)
    for (j=0; j<16; j++)
      {
	k = i*16 + j;
	if (pos_x <= CPDpos_leftDownCorner[k][0] && pos_x > (CPDpos_leftDownCorner[k][0]-CPDlength) &&
	    pos_y >= CPDpos_leftDownCorner[k][1] && pos_y < (CPDpos_leftDownCorner[k][1]+CPDlength) )
	  {		  
	    //printf ("CPD %d: position left down corner = %.2f, \t%.2f \t track position = %.2f, \t %.2f\n", 
	    //  k, CPDpos_leftDownCorner[k][0], CPDpos_leftDownCorner[k][1], pos_x, pos_y);
	    *cpd_index = k;		  
	    break;
	  }	      
      }
  // Define Cell index if CPD has been found
  if (*cpd_index >= 0)
    for (m=0; m<8; m++)
      for (n=0; n<8; n++)
	{
	  l = m*8 + n;
	  //printf ("CELL %d in CPD %d: position left down corner = %.2f, \t%.2f\n", 
	  //      l, *cpd_index, CELLpos_leftDownCorner[*cpd_index][l][0], CELLpos_leftDownCorner[*cpd_index][l][1]);
	  
	  if (pos_x <= CELLpos_leftDownCorner[*cpd_index][l][0] && pos_x > (CELLpos_leftDownCorner[*cpd_index][l][0]-CELLlength) &&
	      pos_y >= CELLpos_leftDownCorner[*cpd_index][l][1] && pos_y < (CELLpos_leftDownCorner[*cpd_index][l][1]+CELLlength) )
	    {		  
	      //printf ("Cell %d in CPD %d: position left down corner = %.2f, \t%.2f \t track position = %.2f, \t %.2f\n", 
	      //  l, *cpd_index, CELLpos_leftDownCorner[*cpd_index][l][0], CELLpos_leftDownCorner[*cpd_index][l][1], pos_x, pos_y);
	      *cell_index = l;		  
	    }
	}
  //printf ("cpd_index = %3d \t cell_index = %3d\n", *cpd_index, *cell_index);
  
}
