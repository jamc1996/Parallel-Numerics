#ifndef STRUCTS_H
#define STRUCTS_H

/* 		structs.h -- declaration of structs for implementation of the conjugate gradient method in serial and parallel.
 *
 *		Author: John Cormican
 */

// structure for serial implementation
typedef struct
{
	int num_sections;
	int len[4];
	int start[4];
	int size;
	int total_size;
	double * data1d;
	double ** data2d;

} Grid;


// structure for parallel implementation
typedef struct
{
	double * data1d;
	double ** data2d;

	double* halos1d;
	double** halos2d;
	int halo_size;	

	int id;
	int nbrup, nbrdown, nbrleft, nbrright;

	int width;
	int height;
	int size;

} Block;



#endif
