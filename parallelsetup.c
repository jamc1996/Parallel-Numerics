#include "parallelsetup.h"


// This function was originally written for 5611, and slightly editted for this assignment.

int decomp1d(int offset, int n, int nprocs, int myid, int *sp, int *ep)
/* Function for each processor to find it's block of the matrix to qr factorize. */
{
  // partsize is the number ALL processors must take.
  // extras is the number of processors that take an extra row.
  int partsize = (n-offset)/nprocs;
  int extras = (n-offset)%nprocs;
   
  //These processors take an 1 extra row.
  if (myid<extras)
  {
    partsize += 1;
    *sp = myid*partsize + offset; 
    *ep = myid*partsize + offset + partsize - 1;
    return 0;
  }

  //These processes do not need to take extra row.
  else
  {
    *sp = myid*partsize + offset + extras;
    *ep = myid*partsize + offset + extras + partsize - 1;
    return 0;
  }
}
