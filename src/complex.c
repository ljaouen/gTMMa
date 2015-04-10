/* This file is part of gTMMa.
 * Copyright (c) 2004, 2013, Luc Jaouen 
 * under the terms of the new BSD License. 
 *
 * The complex/math.c file of the Gnu Scientific Library (GSL)
 * was a great source of inspiration as well as the references its
 * authors provide.
 *
 * GSL/complex/math.c is Copyright (C) 1996, 1997, 1998, 1999, 2000
 * Jorma Olavi Tähtinen, Brian Gough.
 *  
 * References found in GSL/complex/math.c:
 *
 *   T. E. Hull and Thomas F. Fairgrieve and Ping Tak Peter Tang,
 *   "Implementing Complex Elementary Functions Using Exception
 *   Handling", ACM Transactions on Mathematical Software, Volume 20
 *   (1994), pp 215-244, Corrigenda, p553
 *
 *   Hull et al, "Implementing the complex arcsin and arccosine
 *   functions using exception handling", ACM Transactions on
 *   Mathematical Software, Volume 23 (1997) pp 299-335
 *
 *   Abramowitz and Stegun, Handbook of Mathematical Functions, "Inverse
 *   Circular Functions in Terms of Real and Imaginary Parts", Formulas
 *   4.4.37, 4.4.38, 4.4.39
 */

#include <stdio.h>
#include <math.h>

#include "gtmma.h"
#include "complex.h"

/* definition: c = a + i x b with a and b reals */
struct aComplex 
cpx ( double a, 
      double b  ) {

  struct aComplex c;

  c.re = a;
  c.im = b;
  return c;
}

/* c = a + b */
struct aComplex
cpx_add ( struct aComplex a, 
	  struct aComplex b  ) {

  struct aComplex c;

  c.re = a.re + b.re;
  c.im = a.im + b.im;
  return c;
}

/* c = a - b */
struct aComplex
cpx_sub ( struct aComplex a, 
	  struct aComplex b  ) {

  struct aComplex c;

  c.re = a.re - b.re;
  c.im = a.im - b.im;
  return c;
}

/* c = b x a with b real */
struct aComplex
cpx_real_mult ( struct aComplex a,
		double b ) {

  struct aComplex c;

  c.re = a.re * b;
  c.im = a.im * b;
  return c;
}

/* c = a x b */
struct aComplex
cpx_mult ( struct aComplex a, 
	   struct aComplex b  ) {

  struct aComplex c;

  c.re = a.re * b.re - a.im * b.im;
  c.im = a.im * b.re + a.re * b.im;
  return c;
}

/* c = c* */
struct aComplex
cpx_conj ( struct aComplex a ) {

  struct aComplex c;

  c.re = a.re;
  c.im = -a.im;
  return c;
}

/* c = a/b */
struct aComplex
cpx_div ( struct aComplex a, 
	  struct aComplex b  ) {

  struct aComplex c;
  double r, tmp;

  if( fabs(b.re) >= fabs(b.im) ) {
    r = b.im / b.re;
    tmp = b.re + r * b.im;
    c.re = (a.re + r * a.im) / tmp;
    c.im = (a.im - r * a.re) / tmp;
  } else {
    r = b.re / b.im;
    tmp = b.im + r * b.re;
    c.re = (a.re * r + a.im) / tmp;
    c.im = (a.im * r - a.re) / tmp;
  }
  return c;
}

/* c = 1/a */
struct aComplex
cpx_inv ( struct aComplex a) {  

  double b = 1.0 / cpx_abs( a );

  struct aComplex c;
  c.re = (a.re * b) * b;
  c.im = -(a.im * b) * b;
  return c;
}

/* c = | a | */
double 
cpx_abs ( struct aComplex a ) {

  double x, y, ans, tmp;

  x = fabs(a.re);
  y = fabs(a.im);
  if (x == 0.0) 
    ans = y;
  else if (y == 0.0)
    ans = x;
  else if (x > y) {
    tmp = y / x;
    ans = x * sqrt(1.0 + tmp * tmp);
  } else {
    tmp = x / y;
    ans = y * sqrt(1.0 + tmp * tmp);
  }
  return ans;
}

/* c = sqrt(a) */
struct aComplex
cpx_sqrt ( struct aComplex a ) {

  struct aComplex c;
  double x, y, w, r;

  if( (a.re == 0.0) && (a.im == 0.0) ) {
    c.re = 0.0;
    c.im = 0.0;
    return c;
  } else {
    x = fabs(a.re);
    y = fabs(a.im);
    if(x >= y) {
      r = y / x;
      w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r * r)));
    } else {
      r = x / y;
      w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r * r)));
    }
    if(a.re >= 0.0) {
      c.re = w;
      c.im = a.im / (2.0 * w);
    } else {
      c.im = (a.im >= 0.0) ?  w : -w;
      c.re = a.im / (2.0 * c.im);
    }
    return c;
  }
}

/*
 * TRIGONOMETRY 
 */
/* c = cos(a) */
struct aComplex
cpx_cos( struct aComplex a ) {

  struct aComplex c;
  double R = a.re, I = a.im;

  if (I == 0.0) {
    c.re = cos(R);
    c.im = 0.0;
  } 
  else {
    c.re = cos (R) * cosh (I);
    c.im = sin (R) * sinh (-I);
  }

  return c;
}

/* c = sin(a) */
struct aComplex
cpx_sin( struct aComplex a ) {
 
  struct aComplex c;
  double R = a.re, I = a.im;

  if (I == 0.0) {
    c.re = sin(R);
    c.im = 0.0;
  } 
  else {
    c.re = sin (R) * cosh (I);
    c.im = cos (R) * sinh (I);
  }

  return c;
}

/* c = tan(a) */
struct aComplex
cpx_tan( struct aComplex a ) {

  struct aComplex c;
  double R = a.re, I = a.im;

  if ( fabs(I) < 1 ) {
      double D = pow( cos(R), 2.0 ) + pow( sinh(I), 2.0 );

      c.re = 0.5*sin( 2*R)/D;
      c.im = 0.5*sinh( 2*I )/D;
  }
  else {
      double u = exp (-I);
      double C = 2 * u / ( 1 - pow( u, 2.0 ) );
      double D = 1 + pow( cos(R), 2.0 ) * pow( C, 2.0 );

      double S = pow( C, 2.0 );
      double T = 1.0/tanh(I);

      c.re = 0.5*sin( 2*R )*S/D;
      c.im = T/D;
  }
  return c;
}

/* c = tanh(a) */
struct aComplex
cpx_tanh( struct aComplex a ) {

  struct aComplex c;
  double R = a.re, I = a.im;

  if (fabs(R) < 1.0) {
    double D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);
      
    c.re = sinh (R) * cosh (R) / D;
    c.im = 0.5 * sin (2 * I) / D;
  }
  else {
    double D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);
    double F = 1 + pow (cos (I) / sinh (R), 2.0);

    c.re = 1.0 / (tanh (R) * F);
    c.im = 0.5 * sin (2 * I) / D;
  }

  return c;
}



/* 

 cpx_exp(complex a)
{
   complex c;
   c.re = exp(a.re);
   if (fabs(a.im)==0) c.im = 0; 
   else  { c.im = c.re*sin(a.im); c.re*=cos(a.im); }
   return (c);
}



complex cpx_log(complex a)
{
    complex c;
    c.re = log(cpx_abs(a));
    if(a.re == 0.0) {
        c.im = PIOVER2;
    } else {
        c.im = atan2(a.im, a.re);
    }
    return c;
}

*/

/*
 * Invert matrix a using the LU decomposition technique. 
 * Inverse matrix is a_inv, matrix a is destroyed.
 * Returns 1 if matrix is singular, 0 otherwise.
 *
 * This function is adapted from paragraph 2.3 of the
 * Numerical Recipes in C (LU Decomposition and its applications)
 * cf. http://www.nr.com/
 *
 * complex **a:     matrix represented as vector of row pointers 
 * int     n:       order of matrix 
 * double  *dwork:  work vector of size n
 * int     *indx:   work vector of size n 
 * complex **a_inv: inverse of input matrix a (matrix a is destroyed) 
 * complex *col:    for the second part, involving back subst 
 */

int 
invert_matrix ( struct aComplex** a, 
		int n, 
		double* dwork, 
		int* indx, 
		struct aComplex* col     ) {

  int pb;
    
  pb = LU_decompose( a, n, dwork, indx, (double*) NULL );

  if ( pb == 0 ) {
    LU_back_substitution( a, n, indx, col );
  }
  return pb;
}

/*
 * LU decomposition of a matrix a.
 * Returns 1 if matrix is singular, 0 otherwise.
 *
 * This function is adapted from paragraph 2.3 of the
 * Numerical Recipes in C (LU Decomposition and its applications)
 * cf. http://www.nr.com/
 *
 * complex **a:   the matrix whose LU-decomposition is wanted 
 * int     n:     order of a
 * double  *vv:   work vector of size n (stores implicit
 *                scaling of each row)
 * int     *indx: => row permutation according to partial
 *                   pivoting sequence
 * double  *pd:   => 1 if number of row interchanges was even,
 *                   -1 if odd (NULL OK)
 */

int 
LU_decompose ( struct aComplex** a, 
	       int n, 
	       double* vv, 
	       int* indx, 
	       double* pd   ) {

  int i, imax, j, k;
  double big, dum, temp, d;
  struct aComplex sum, cdum;

  d = 1.0;
  imax = 0;

  for ( i = 0; i < n; i++ ) {
    big = 0.0;
    for ( j = 0; j < n; j++ ) {
      if ( (temp = cpx_abs(a[i][j])) > big ) {
	big = temp;
      }
    }
    if ( big == 0.0 ) {
      printf("Singular matrix in routine LU_decompose\n");
      return 1;
    }
    vv[i] = 1.0 / big;
  }
    
  for ( j = 0; j < n; j++ ) {
    for ( i = 0; i < j; i++ ) {
      sum = a[i][j];
      for ( k = 0; k < i; k++ ) {
	sum = cpx_sub( sum, cpx_mult( a[i][k], a[k][j] ) );
      }
      a[i][j] = sum;
    }
    big = 0.0;
    for ( i = j; i < n; i++ ) {
      sum = a[i][j];
      for ( k = 0; k < j; k++ ) {
	sum = cpx_sub( sum, cpx_mult( a[i][k], a[k][j] ) );
      }
      a[i][j] = sum;
      dum = vv[i] * cpx_abs( sum );
      if ( dum >= big ) {
	big = dum;
	imax = i;
      }
    }
    if ( j != imax ) {
      for (k = 0; k < n; k++) {
	cdum = a[imax][k];
	a[imax][k] = a[j][k];
	a[j][k] = cdum;
      }       
      d = -d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if ( a[j][j].re == 0.0 && a[j][j].im == 0.0 ) {
      a[j][j] = cpx( 1.0e-20, 1.0e-20 );
    }
    if (j != n - 1) {
      cdum = cpx_div( cpx(1.0, 0.0), a[j][j] );
      for (i = j + 1; i < n; i++) {
	a[i][j] = cpx_mult( a[i][j], cdum );
      }
    }
  }

  if ( pd != NULL ) {
    *pd = d;
  }
  return 0;
}

/*
 * Back-substitution into LU-decomposed matrix
 *
 * This function is adapted from paragraph 2.3 of the
 * Numerical Recipes in C (LU Decomposition and its applications)
 * cf. http://www.nr.com/
 */

void 
LU_back_substitution ( struct aComplex **a, 
		       int n, 
		       int *indx, 
		       struct aComplex *b   ) {

  int i, ip, j;
  int ii = -1;

  struct aComplex sum;

  for ( i = 0; i < n; i++ ) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if ( ii >= 0 ) {
      for ( j = ii; j <= i - 1; j++ ) {
	sum = cpx_sub( sum, cpx_mult( a[i][j], b[j] ) );
      }
    } 
    else if ( (sum.re != 0.0) || (sum.im != 0.0) ) {
      ii = i;
    }
    b[i] = sum;
  }
    
  for ( i = n - 1; i >= 0; i-- ) {
    sum = b[i];
    for (j = i + 1; j < n; j++) {
      sum = cpx_sub( sum, cpx_mult( a[i][j], b[j] ) );
    }
    b[i] = cpx_div( sum, a[i][i] );
  }
}


/*
 * Complex vectors and matrices
 */

/* Compute (a x b) and store the result in c */
void
cpx_matrices_mult(struct aComplexMatrix* a, 
		  struct aComplexMatrix* b, 
		  struct aComplexMatrix* c) {

  int i, j, k; 
  int nb_rows_a = a->nb_rows;
  int nb_columns_a = a->nb_columns; 
  int nb_columns_b = b->nb_columns; 
  struct aComplex cpx_tmp;

  for ( i = 0; i < nb_rows_a; i++ ) {
    for ( j = 0; j < nb_columns_b; j++ ) {
      cpx_tmp = cpx(0.0e00,0.0e00);
      for ( k = 0; k < nb_columns_a; k++ ) {
	cpx_tmp = cpx_add( cpx_tmp, cpx_mult(a->data[i][k], b->data[k][j]) );
      }
      c->data[i][j] = cpx_tmp;
    }
  }
}

/* Copy matrix b in matrix a (a = b) */
void
cpx_copy_matrix( struct aComplexMatrix* a, 
		 struct aComplexMatrix* b  ) {

  int i, j;
  int nb_rows = b->nb_rows;
  int nb_columns = b->nb_columns;

  a->nb_rows = nb_rows;
  a->nb_columns = nb_columns;

  for ( i = 0; i < nb_rows; i++ ) {
    for ( j = 0; j < nb_columns; j++ ) {
      a->data[i][j] = b->data[i][j];
    }
  }

}
