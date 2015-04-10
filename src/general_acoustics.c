/* This file is part of gTMMa.
 * Copyright (c) 2004, 2013, Luc Jaouen 
 * under the terms of the new BSD License. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "complex.h"
#include "gtmma.h"
#include "misc_operations.h"

/* Properties of the air at temperature and pressure 
 * specified by the user in parameters file 
 */
void 
air_properties( struct aFluid* fluid ) {

  double R;

  fluid->T += 273.16; /* Temperature in K */ 
  R = 287.031;        /* J.Kg-1.K-1: constant (or 8.314 J.mol-1.K-1) */
  /* specific heat per unit mass at constant pressure in  J.kg-1.K-1 (260 K < T < 600 K): */
  fluid->Cp = 4168.8 * ( 0.249679 - 7.55179e-05 * fluid->T 
		         + 1.69194e-07 * pow(fluid->T,2.) 
		         - 6.46128e-11 * pow(fluid->T,3.) ); 
  fluid->gamma = fluid->Cp / ( fluid->Cp - R ); /* Ratio of specific heats */
  fluid->c_0 = sqrt(fluid->gamma*R*fluid->T);   /* m.s-1 : speed of acoustic waves */
  fluid->rho_0 = fluid->P_0/(R*fluid->T);       /* kg.m-3: mass per unit volume */
  /* Dynamic viscosity in  N.s.m-2 (100 K < T < 600 K): */
  fluid->eta = 7.72488e-08 * fluid->T
               - 5.95238e-11 * pow(fluid->T,2.) 
               + 2.71368e-14 * pow(fluid->T,3.);
  /* Thermal conductivity - cf Pierce p 513 - (W.m-1.K-1) */
  fluid->kappa =  2.624e-02 * 
    ( pow(fluid->T/300,1.5) * (300+245.4*exp(-27.6/300))/(fluid->T+245.4*exp(-27.6/fluid->T)) );

  fluid->Pr = ( fluid->eta * fluid->Cp ) / fluid->kappa; /* Prandtl number */

}


/* 
 * Compute the characteristic impedance Zc and the wave number k of a
 * layer from its equivalent volumic mass rho and equivalent bulk
 * modulus K 
 */
void
Zc_k_from_rho_K( struct aComplex* Zc, 
		 struct aComplex* k, 
		 struct aComplex* rho, 
		 struct aComplex* K,
		 double omega          ) {

  struct aComplex tmp1;

  tmp1 = cpx_mult( *rho, *K );
  *Zc = cpx_sqrt( tmp1 );

  tmp1 = cpx_div( *rho, *K );
  tmp1 = cpx_sqrt( tmp1 );
  *k = cpx_real_mult( tmp1, omega );

}


/* 
 * Compute the impedance Z of a multilayer under the
 * hypothesis: 
 *  o the first layer is backed by a rigid wall and
 *  o the last layer is a fluid or equivalent fluid 
 */
void
backed_layer_impedance( struct aComplexMatrix* a,
			struct aComplex* Z,
			struct aComplex* k         ) {

  *Z = cpx_div( a->data[0][0], a->data[1][0] );

}

/* 
 * Compute the sound absorption coefficient alpha
 * from the impedance Z, the incident angle and the
 * characteristic impedance of air Z0 = rho0*c0 
 *
 *
 *
 *             | Z - Z0/cos(angle) | 2
 * alpha = 1 - | ----------------- |
 *             | Z + Z0/cos(angle) |
 *
 */
double
compute_alpha( struct aComplex* Z, 
	       double Z0,
               double angle        ) {

  struct aComplex ctmp1, ctmp2;
  double dtmp, alpha;

  ctmp1.re = Z->re - Z0/cos(angle);
  ctmp1.im = Z->im;
  ctmp2.re = Z->re + Z0/cos(angle);
  ctmp2.im = Z->im;
  ctmp1 = cpx_div( ctmp1, ctmp2 );

  dtmp = cpx_abs( ctmp1 ) * cpx_abs( ctmp1 );

  alpha = 1 - dtmp;

  return alpha;
}

/* Compute the sound absorption coefficient alpha in diffuse field
 * from the knowledge of the impedance at normal incidence by using
 * London's formula (locally reacting material)
 *
 * Bibliography:
 * A. London, J. Acoust. Soc. Am. 22(2), 1950.
 *
 */
double
compute_alpha_diffuse_london( struct aComplex* Z,
			      double Z0           ) {

  double r, x, alpha_diffuse;
  r = Z->re / Z0; 
  x = Z->im / Z0;

  alpha_diffuse = 8*r / ( pow(x,2) + pow(r,2) ) * 
      ( 1+(pow(r,2)-pow(x,2))/(x*(pow(x,2)+pow(r,2))) * atan(x/(1+r)) 
	- r/(pow(x,2)+pow(r,2))*log(1+2*r+pow(r,2)+pow(x,2))          );

  return alpha_diffuse;
}


/*
 * Narrow bands to one-third octave bands representation.
 *
 */

void
one_third_octave( double* frequencies,
		  double* alpha_diffuse,
		  int nb_frequency_points,
		  /*		  double one_third_freq[12], */
		  double bands[12]) {

  int a, b;
  double one_third_bands[2][12];

  int* idx;
  int length_idx = 0;
  idx = (int*) malloc ( sizeof(int) );
  

  double one_third_freq[12] = { 200., 250., 315., 400., 500., 630., 800., 1000., 1250., 1600., 2000., 2500};

  /* Determine lower and upper limits of each 1/3 octave band */
  for ( a = 0; a < 12; a++) {
    one_third_bands[0][a] = one_third_freq[a]/pow(2,1/6); 
    one_third_bands[1][a] = one_third_freq[a]*pow(2,1/6); 
  }

  /* Compute the acoustic absorption coefficient per 1/3 octave band */
  for ( a = 0; a < 12; a++) {
    bands[a] = 0;
    find_freq_points ( one_third_bands[0][a], 
		       one_third_bands[1][a], 
		       frequencies, 
		       nb_frequency_points,
		       idx,
		       &length_idx );

    /* If we have no 'measurement' point in this band: */
    if ( length_idx == 0 ) {
      printf("Warning: no point found in band centered at %4.0f\n",one_third_freq[a] );
    }
    /* If we have only 1 'measurement' point in this band: */
    else if ( length_idx == 1 ) {
      printf("Warning: only one point found in band centered at %4.0f\n",one_third_freq[a] );
      bands[a] = alpha_diffuse[idx[0]];
    }
    /* If we have more than 1 'measurement' point in this band: */
    else if ( length_idx > 1 ) {
      for ( b = 0; b < length_idx-2; b++ ) {
	bands[a] = bands[a] +
		  ( frequencies[idx[0]+b+1] - frequencies[idx[0]+b] ) * 
		  fabs( alpha_diffuse[idx[0]+b+1] + alpha_diffuse[idx[0]+b] ) / 2;
      }
      bands[a] = bands[a] / ( frequencies[idx[length_idx-1]]-frequencies[idx[0]] );
    }
  }

  free( idx );
}



/* Build the 2x2 transfert matrix for a fluid from the knowledge of
 * Zc, k and L.
 *
 *                                              k
 *           cos( k3*L )              j * Zc * ---- * sin( k3*L )
 *                                              k3 
 *
 *        1     k3 
 *   j * --- * ---- * sin( k3*L )            cos( k3*L )
 *        Zc    k
 *
 */
void 
build_fluid_matrix( struct aComplex Zc,
		    struct aComplex k,
		    struct aFluid* fluid,
		    double angle,
		    double L,
		    double omega,
		    struct aComplexMatrix* c ) {

  struct aComplex cpx_tmp1, cpx_tmp2, cpx_tmp3;
  struct aComplex k3;
  double dtmp;

  double c_0 = fluid->c_0;

  k3     = cpx_mult( k, k );
  dtmp   = omega/c_0*sin(angle);
  k3.re -= pow(dtmp,2);
  k3     = cpx_sqrt( k3 );

  /* k3 = cpx_real_mult( k, cos(angle) ); */

  cpx_tmp1 = cpx_real_mult( k3, L );

  c->data[0][0] = cpx_cos( cpx_tmp1 );
  c->data[1][1] = c->data[0][0];

  cpx_tmp1 = cpx_sin( cpx_tmp1 );
  cpx_tmp2 = cpx_mult( Zc, cpx_tmp1 );
  cpx_tmp3 = cpx_div( k, k3 );
  cpx_tmp2 = cpx_mult( cpx_tmp3, cpx_tmp2 );

  c->data[0][1].re = -cpx_tmp2.im;
  c->data[0][1].im = cpx_tmp2.re;

  cpx_tmp2 = cpx_inv( Zc );
  cpx_tmp2 = cpx_mult( cpx_tmp2, cpx_tmp1 );
  cpx_tmp3 = cpx_inv( cpx_tmp3 );
  cpx_tmp2 = cpx_mult( cpx_tmp2, cpx_tmp3 );

  c->data[1][0].re = -cpx_tmp2.im;
  c->data[1][0].im = cpx_tmp2.re;

  /*
  cpx_tmp1 = cpx_real_mult( k, L );

  c->data[0][0] = cpx_cos( cpx_tmp1 );
  c->data[1][1] = c->data[0][0];

  cpx_tmp1 = cpx_sin( cpx_tmp1 );
  cpx_tmp2 = cpx_mult( Zc, cpx_tmp1 );

  c->data[0][1].re = -cpx_tmp2.im;
  c->data[0][1].im = cpx_tmp2.re;

  cpx_tmp2 = cpx_inv( Zc );
  cpx_tmp2 = cpx_mult( cpx_tmp2, cpx_tmp1 );

  c->data[1][0].re = -cpx_tmp2.im;
  c->data[1][0].im = cpx_tmp2.re;
  */

}

