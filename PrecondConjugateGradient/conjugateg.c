#include "conjugateg.h"

void parallel_precond_cg(Block *x, Block*b, MPI_Comm comm, MPI_Datatype vec)
{
	int i;
	Block r,p,Ap,z;
	double local_rp, local_zrp, local_App, localnew_zrp;
	double r_prod, zr_prod, new_zrprod, App;
	double alpha, beta;
	int it_count=0;

	alloc_block(&r,x->size,x->id);
	alloc_block(&p,x->size,x->id);
	alloc_block(&Ap,x->size, x->id);
	alloc_block(&z,x->size, x->id);
	
	par_calc_error_vector(&r, x, b);

	local_rp = par_inner_prod(&r,&r);
	MPI_Allreduce(&local_rp, &r_prod, 1, MPI_DOUBLE, MPI_SUM, comm); 


	par_rb_gauss_seidel(&z, &r, comm, vec);

	copy_block(&p, &r);

	local_zrp = par_inner_prod(&z,&r);
	MPI_Allreduce(&local_zrp, &zr_prod, 1, MPI_DOUBLE, MPI_SUM, comm);	
	while(sqrt(r_prod) >0.0001 && it_count<63)
	{
		par_A_mul(&Ap, &p);

		local_App = par_inner_prod(&Ap, &p);	

		MPI_Allreduce(&local_App, &App, 1, MPI_DOUBLE, MPI_SUM, comm); 
		alpha = zr_prod/App;
		
		par_linear_transform(1.0, x, alpha, &p);
		par_linear_transform(1.0, &r, -alpha, &Ap);

		par_rb_gauss_seidel(&z, &r, comm, vec);

		localnew_zrp = par_inner_prod(&z, &r);
		MPI_Allreduce(&localnew_zrp, &new_zrprod, 1, MPI_DOUBLE, MPI_SUM, comm);

		beta = new_zrprod/zr_prod;
		par_linear_transform(beta, &p, 1.0, &z);

		zr_prod = new_zrprod;
		it_count ++;
		local_rp = par_inner_prod(&r,&r);
		MPI_Allreduce(&local_rp, &r_prod, 1, MPI_DOUBLE, MPI_SUM, comm);
		
		halo_swapping(&p, comm, vec);
		halo_swapping(x, comm, vec);
		halo_swapping(&r, comm, vec);
	}
	printf("\n\nrprod is %lf, iterations = %d\n\n", r_prod, it_count);
	free_block(&z);
	free_block(&r);
	free_block(&p);
	free_block(&Ap);
}


void precond_conjugate_gradient(Grid *x, Grid *b)
{
	Grid r,p, Ap, z;
	alloc_serial_grid(&r, x->size);
	alloc_serial_grid(&p, x->size);
	alloc_serial_grid(&Ap, x->size);
	alloc_serial_grid(&z, x->size);
	calc_error_vector(&r, x, b);
	rb_gauss_seidel(&z, &r);
	copy_grid(&p, &r);
	double zr_prod, r_prod, new_zrprod;
	double alpha, beta;
	r_prod = inner_prod(&r, &r);
	zr_prod = inner_prod(&z,&r);
	printf("r_prod is %lf\n",r_prod);
	int it_count = 0;
	while(sqrt(r_prod) > 0.0001 && it_count<1000)
	{
		A_mul(&Ap,&p);	

		alpha = zr_prod/inner_prod(&Ap, &p);

		linear_transform(1.0, x, alpha, &p);

		linear_transform(1.0, &r, -alpha, &Ap);

		rb_gauss_seidel(&z, &r);

		new_zrprod = inner_prod(&z, &r);

		beta = new_zrprod/zr_prod;

		linear_transform(beta, &p, 1.0, &z);

		zr_prod = new_zrprod;
		it_count ++;
		r_prod = inner_prod(&r,&r);
	}
	printf("\n\nrprod is %lf, iterations = %d\n\n", r_prod, it_count);
	free_grid(&r);
	free_grid(&p);
	free_grid(&Ap);
}


void conjugate_gradient(Grid *x, Grid *b)
{
	Grid r,p, Ap;
	alloc_serial_grid(&r, x->size);
	alloc_serial_grid(&p, x->size);
	alloc_serial_grid(&Ap, x->size);

	calc_error_vector(&r, x, b);
	copy_grid(&p, &r);


	double r_prod, new_rprod, pr_prod, new_prprod;
	double alpha, beta;

	

	r_prod = inner_prod(&r,&r);
	printf("r_prod is %lf\n",r_prod);
	int it_count = 0;


	while(sqrt(r_prod) > 0.0001)
	{
		A_mul(&Ap,&p);	

		alpha = r_prod/inner_prod(&Ap, &p);

		linear_transform(1.0, x, alpha, &p);

		linear_transform(1.0, &r, -alpha, &Ap);

		new_prprod = inner_prod(&r, &r);

		beta = new_rprod/r_prod;

		linear_transform(beta, &p, 1.0, &r);

		pr_prod = new_prprod;
		it_count ++;
		r_prod = inner_prod(&r,&r);
	}
	printf("\n\nrprod is %lf, iterations = %d\n\n", r_prod, it_count);
	free_grid(&r);
	free_grid(&p);
	free_grid(&Ap);
}
