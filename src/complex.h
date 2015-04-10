/* This file is part of gTMMa.
 * Copyright (c) 2004, 2013, Luc Jaouen 
 * under the terms of the new BSD License. 
 */

struct aComplex {
  double re;
  double im;
};

struct aComplexMatrix {
  struct aComplex** data;
  int nb_rows;
  int nb_columns;
};

struct aComplex cpx( double, double );
struct aComplex cpx_add( struct aComplex, struct aComplex );
struct aComplex cpx_sub( struct aComplex, struct aComplex );
struct aComplex cpx_real_mult( struct aComplex, double ); 
struct aComplex cpx_mult( struct aComplex, struct aComplex );
struct aComplex cpx_conj( struct aComplex );
struct aComplex cpx_div( struct aComplex, struct aComplex );
struct aComplex cpx_inv( struct aComplex );
double  cpx_abs( struct aComplex );
struct aComplex cpx_sqrt( struct aComplex );

struct aComplex cpx_cos( struct aComplex );
struct aComplex cpx_sin( struct aComplex );
struct aComplex cpx_tan( struct aComplex );
struct aComplex cpx_tanh( struct aComplex );

/*
struct aComplex RCmul ( double, struct aComplex );
struct aComplex Cexp ( struct aComplex );
struct aComplex Clog ( struct aComplex );
struct aComplex **  pscmatrix (int dim);   // mallocs a square struct aComplex matrix
void        free_pscmatrix (struct aComplex **m); // frees the square struct aComplex matrix
void        dump_pscmatrix (struct aComplex **m, int dim); // prints a square struct aComplex matrix
void        copy_pscmatrix (struct aComplex **from, struct aComplex **to, int dim); // copies
void        dump_struct aComplexVector (struct aComplex *vec, int dim);
*/
int invert_matrix( struct aComplex**, int, double*, int*, struct aComplex* );
int LU_decompose( struct aComplex**, int, double*, int*, double* );
void LU_back_substitution( struct aComplex**, int, int*, struct aComplex* );

void cpx_matrices_mult( struct aComplexMatrix*, 
		        struct aComplexMatrix*, 
		        struct aComplexMatrix* );

void cpx_copy_matrix( struct aComplexMatrix*, 
		      struct aComplexMatrix*  );

