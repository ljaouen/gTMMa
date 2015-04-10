/* This file is part of gTMMa.
 * Copyright (c) 2004, 2013, Luc Jaouen 
 * under the terms of the new BSD License. 
 *
 * Garai-Pompoli model added to previous version from a request of and thanks to Davide Borelli
 */

#include <math.h>
#include <stdio.h>

#include "complex.h"
#include "gtmma.h"
#include "general_acoustics.h"


/* Delany and Bazley model
 * M. E. Delany and E. N. Bazley, 
 * Acoustical properties of fibrous absorbent materials,
 * Applied Acoustics (3), 1970, pp. 105-116
 */
void
DB_Zc_k ( struct aComplex* rho,
	  struct aComplex* K,
	  struct aComplex* Zc,
	  struct aComplex* k,
	  double omega,
	  struct aFluid* fluid,
	  double sigma          ) {

  struct aComplex cpx_tmp;
  /* TODO: Warn about boundaries: X = rho_0*f/sigma, 0.01 < X < 1.0 */

  double rho_0 = fluid->rho_0;
  double c_0 = fluid->c_0;
  double f_over_sigma = omega/(sigma*2*PI);
  
  Zc->re =  rho_0*c_0*( 1 + 9.08 * pow(f_over_sigma*1000,-0.75) );
  Zc->im = -rho_0*c_0*11.9 * pow(f_over_sigma*1000,-0.73); 

  k->re  =  omega/c_0 * ( 1 + 10.8*pow(f_over_sigma*1000,-0.70) );
  k->im  = -omega/c_0 * 10.3*pow(f_over_sigma*1000,-0.59);

  cpx_tmp = cpx_mult( *Zc, *k );
  *rho = cpx_real_mult ( cpx_tmp, 1/omega );
  cpx_tmp = cpx_div( *Zc, *k );
  *K = cpx_real_mult( cpx_tmp, omega );

}

/* Garai-Pompoli model
 * M. Garai and F. Pompoli
 * A simple empirical model of polyester fibre materials for acoustical applications
 * Applied Acoustics 66, 2005, pp. 1383–1398
 */

void
GP_Zc_k ( struct aComplex* rho,
	  struct aComplex* K,
	  struct aComplex* Zc,
	  struct aComplex* k,
	  double omega,
	  struct aFluid* fluid,
	  double sigma          ) {

  struct aComplex cpx_tmp;
  /* TODO: Warn about boundaries */

  double rho_0 = fluid->rho_0;
  double c_0 = fluid->c_0;
  double X = omega*rho_0/(sigma*2*PI);
  
  Zc->re =  rho_0*c_0*( 1 + 0.078 * pow(X,-0.623) );
  Zc->im = -rho_0*c_0*0.074 * pow(X,-0.660); 

  k->re  =  omega/c_0 * ( 1 + 0.121*pow(X,-0.530) );
  k->im  = -omega/c_0 * 0.159*pow(X,-0.571);

  cpx_tmp = cpx_mult( *Zc, *k );
  *rho = cpx_real_mult ( cpx_tmp, 1/omega );
  cpx_tmp = cpx_div( *Zc, *k );
  *K = cpx_real_mult( cpx_tmp, omega );

}

struct aComplex 
function_Gj ( double rho_0,
	      double eta,
	      double phi,
	      double alpha_infty,
	      double sigma,
	      double lambda,
	      double omega          ) {

  struct aComplex tmp;
  struct aComplex Gj;

  tmp.re = 1.0;
  tmp.im = ( 4*pow(alpha_infty,2)*eta*rho_0*omega ) /
    ( pow(sigma,2.0)*pow(lambda,2.0)*pow(phi,2.0) );
  
  Gj = cpx_sqrt ( tmp );

  return Gj;
  
}

/* Johnson Champoux Allard model
 * 
 * Johnson et al 
 * Champoux & Allard
 *
 */
void
JCA_Zc_k ( struct aComplex* rho,
	   struct aComplex* K,
	   struct aComplex* Zc,
	   struct aComplex* k,
	   double omega,
	   struct aFluid* fluid,
	   struct aMaterial* material ) {

  struct aComplex tmp1, tmp2;
  struct aComplex Gj;

  double P_0 = fluid->P_0;
  double rho_0 = fluid->rho_0;
  double eta = fluid->eta;
  double gamma = fluid->gamma;
  double Pr = fluid->Pr;

  double phi = material->phi;
  double sigma = material->sigma;
  double alpha_infty = material->alpha_infty;
  double lambda = material->lambda;
  double lambda_p = material->lambda_p;

  /* Frequency dependence function of viscous effects */
  Gj = function_Gj ( rho_0, eta, phi, alpha_infty, sigma, lambda, omega );

  /* Equivalent mass per unit volume rho */
  tmp1.re = 0.0;
  tmp1.im = -sigma*phi / ( omega*rho_0*alpha_infty );
  tmp1 = cpx_mult ( Gj, tmp1 );
  tmp1.re += 1.0; 
  *rho = cpx_real_mult ( tmp1, alpha_infty*rho_0/phi );

  /* Equivalent bulk modulus K */
  tmp1.re = 1.0;
  tmp1.im = rho_0*omega*Pr*pow(lambda_p,2.0) / ( 16*eta) ;
  tmp1 = cpx_sqrt ( tmp1 );
  tmp2.re = 0.0;
  tmp2.im = -8*eta / ( pow(lambda_p,2.0)*Pr*omega*rho_0 ); 
  tmp1 = cpx_mult ( tmp1, tmp2 );
  tmp1.re += 1.0;
  tmp1 = cpx_inv ( tmp1 );
  tmp1 = cpx_real_mult ( tmp1, gamma-1 );
  tmp1.re = gamma - tmp1.re;
  tmp1.im = -tmp1.im;
  tmp2.re = gamma*P_0/phi;
  tmp2.im = 0.0;
  *K = cpx_div ( tmp2, tmp1 );

  /*
  printf("K = %f +i %f\n",K->re,K->im);
  printf("rho = %f +i %f\n",rho->re,rho->im);

  K->re = 1.624224486183634e+05;
  K->im = 3.959935572555331e+03;*/

  Zc_k_from_rho_K( Zc, k, rho, K, omega );

}

/* Johnson Champoux Allard Lafarge model
 * 
 * Johnson et al
 * Champoux & Allard
 * Lafarge et al
 */

void
JL_Zc_k ( struct aComplex* rho,
	    struct aComplex* K,
	    struct aComplex* Zc,
	    struct aComplex* k,
	    double omega,
	    struct aFluid* fluid,
	    struct aMaterial* material ) {

  struct aComplex tmp1, tmp2;
  struct aComplex Gj;

  double P_0 = fluid->P_0;
  double rho_0 = fluid->rho_0;
  double eta = fluid->eta;
  double gamma = fluid->gamma;
  double Cp = fluid->Cp;
  double kappa = fluid->kappa;

  double phi = material->phi;
  double sigma = material->sigma;
  double alpha_infty = material->alpha_infty;
  double lambda = material->lambda;
  double lambda_p = material->lambda_p;
  double k_p_0 = material->k_p_0;

  /* Frequency dependence function of viscous effects */
  Gj = function_Gj ( rho_0, eta, phi, alpha_infty, sigma, lambda, omega );

  /* Equivalent mass per unit volume rho */
  tmp1.re = 0.0;
  tmp1.im = -sigma*phi / ( omega*rho_0*alpha_infty );
  tmp1 = cpx_mult ( Gj, tmp1 );
  tmp1.re += 1.0; 
  *rho = cpx_real_mult ( tmp1, alpha_infty*rho_0/phi );

  /* Equivalent bulk modulus K */
  tmp1.re = 1.0;
  tmp1.im = 4*rho_0*omega*Cp*pow(k_p_0,2) / ( kappa*pow(lambda_p,2)*pow(phi,2) );
  tmp1 = cpx_sqrt ( tmp1 );
  tmp2.re = 0.0;
  tmp2.im = - phi*kappa / ( k_p_0*Cp*omega*rho_0 );
  tmp1 = cpx_mult ( tmp1, tmp2 );
  tmp1.re += 1.0;
  tmp2.re = 1.0;
  tmp2.im = 0.0;
  tmp1 = cpx_div ( tmp2, tmp1 );
  tmp1 = cpx_real_mult ( tmp1, gamma-1 );
  tmp1.re = gamma - tmp1.re;
  tmp1.im = -tmp1.im;
  tmp2.re = gamma*P_0/phi;
  tmp2.im = 0.0;
  *K = cpx_div ( tmp2, tmp1 );

  Zc_k_from_rho_K( Zc, k, rho, K, omega );

}

