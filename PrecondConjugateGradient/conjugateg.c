#include "conjugateg.h"

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


	while(r_prod > 0.001 && it_count<1000)
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


	while(r_prod > 0.0001)
	{
		A_mul(&Ap,&p);	

		alpha = r_prod/inner_prod(&Ap, &p);

		linear_transform(1.0, x, alpha, &p);

		linear_transform(1.0, &r, -alpha, &Ap);

		new_rprod = inner_prod(&r, &r);

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
