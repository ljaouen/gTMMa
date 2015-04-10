/* This file is part of gTMMa.
 * Copyright (c) 2004, 2013, Luc Jaouen 
 * under the terms of the new BSD License. 
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gtmma.h"

/* Compute the number of frequency points in defined spectrum */
int 
get_nb_frequency_points ( struct aSpectrum* spectrum ) {

  int i = 0;
  double tmp = spectrum->min;

  while ( tmp <= spectrum->max ) {
    i++;
    tmp += spectrum->step; 
  }
  
  return i;

}

/* Initialize problem with default values */
void
problem_initialize( struct aProblem* problem ) {

  problem->spectrum->min  = 10.0;
  problem->spectrum->max  = 4000.0;
  problem->spectrum->step = 10.0;

  problem->fluid->T     = 20.0;
  problem->fluid->P_0   = 101310;
  problem->fluid->H     = 40;
  problem->fluid->c_0   = 0.0;
  problem->fluid->rho_0 = 0.0;
  problem->fluid->Z_0   = 0.0;
  problem->fluid->eta   = 0.0;
  problem->fluid->gamma = 0.0;
  problem->fluid->Pr    = 0.0;


  problem->conditions->diffuse_field = 0;
  problem->conditions->angle         = 0.0;

}

/* Initialize material with default values */
void 
material_initialize( struct aMaterial* material ) {

  material->type = AIR_GAP;
  strncpy( material->name, "\0", 1 );
  material->L           = 0.0;
  material->sigma       = 0.0;
  material->phi         = 0.0;
  material->alpha_infty = 0.0;
  material->lambda      = 0.0;
  material->lambda_p    = 0.0;
  material->k_p_0       = 0.0;
}

/* MOVED TO GUI */
/*
struct aMaterial*
material_duplicate( struct aMaterial* material ) {

  struct aMaterial* new_material;

  new_material = (struct aMaterial*) malloc( sizeof(struct aMaterial) );
  
  new_material->type  = material->type;
  strncpy(new_material->name, material->name, 40);
  new_material->L = material->L;
  new_material->phi   = material->phi;
  new_material->sigma = material->sigma;
  new_material->alpha_infty = material->alpha_infty;
  new_material->lambda = material->lambda;
  new_material->lambda_p = material->lambda_p;
  new_material->k_p_0 = material->k_p_0;

  return( new_material );

}
*/


/*
 * function [X W] = gauss_rule(N)
 * Gauss-Legendre integration points and weights
 * 
 * X: coordinates of the N integration points
 * W: corresponding weights
 *
 * TODO: comment and clean up this awful code
 */

void
gauss_legendre_rule( int N, 
		     double* X, 
		     double* W  ) {

  int i,k;

  int M, E1;
  double T, X0, PKM1, T1, PKP1, PK;

  double DEN, D1, DPN, D2PN, D3PN, D4PN, U, V, H, P, DP, FX;

  M  = (int)floor((N+1)/2);
  E1 = N*(N+1);

  for ( i = 1; i <= M; i++ ) {
    T    = (4*i-1) * PI / ( 4*N+2 );
    X0   = ( 1.0- (1.0-1.0/N) / (8.0*N*N) ) * cos(T);
    PKM1 = 1.0;
    PK   = X0;
    for ( k = 2; k <= N; k++ ) {
      T1   = X0*PK;
      PKP1 = T1 - PKM1 - (T1-PKM1)/k + T1;
      PKM1 = PK;
      PK   = PKP1;
    }

    DEN  = 1.0-X0*X0;
    D1   = N * ( PKM1 - X0*PK );
    DPN  = D1/DEN;
    D2PN = ( 2.0*X0*DPN  - E1*PK) / DEN;
    D3PN = ( 4.0*X0*D2PN + (2.0-E1)*DPN) / DEN;
    D4PN = ( 6.0*X0*D3PN + (6.0-E1)*D2PN) / DEN;
    U    = PK / DPN;
    V    = D2PN / DPN;
    H    = -U*(1.0+0.50*U*(V+U*(V*V-D3PN/(3.0*DPN))));
    P    = PK+H*(DPN+0.50*H*(D2PN+H/3.0*(D3PN+0.250*H*D4PN)));
    DP   = DPN+H*(D2PN+0.50*H*(D3PN+H*D4PN/3.0));
    H    = H - P / DP;
    X[i-1] = X0+H;
    FX   = D1-H*E1 * ( PK + 0.50 * H * ( DPN+ H/3.0 * ( D2PN + 0.250 * H * (D3PN+0.20*H*D4PN) ) ) );
    W[i-1] = 2.0 * ( 1.0-X[i-1]*X[i-1] ) / ( FX*FX );
  }

  if ( (M+M) == N ) { /* even number of points: N */
    for ( i = 1; i <= M; i++ ) {
      X[N-i] =  X[i-1];
      X[i-1] = -X[i-1];
      W[N-i] =  W[i-1];
    }
  }
  else {            /* odd number of points: N */
    X[M-1] = 0.0;     /* 0 is a first point */
    for ( i = 1; i <= M-1; i++ ) {
      X[N-i] =  X[i-1];
      X[i-1] = -X[i-1];
      W[N-i] =  W[i-1];  
    }
  }

}

/* Find frequency points in the range [f_min, f_max].
   The result is a vector, idx, which stores index of these frequency points.
 */
void
find_freq_points ( double f_min, 
		   double f_max,
		   double* frequencies,
		   int nb_frequency_points,
		   int* idx,
		   int* length_idx ) {
  int a, length;
  int* tmp_idx;

  tmp_idx = (int*) malloc ( sizeof(int) );
  length = 0;

  for ( a = 0; a < nb_frequency_points; a++ ) {

    if ( frequencies[a] >= f_min && frequencies[a] <= f_max ) { 
      length++;

      if ((tmp_idx = realloc( (int*)idx, length * sizeof(int) )) == NULL) {
	printf("Re-allocation error of idx in function find_freq_points\n");
      }
      tmp_idx[length] = a;
      idx = tmp_idx;
    }
  }
  length_idx = &length;
  free( tmp_idx );
}
