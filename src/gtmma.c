/* 

This code, distributed under the terms of the
new BSD License, is an example of use of
the transfert matrix method in the field of acoustics
 
Copyright (c) 2004, 2013, Luc Jaouen.
All rights reserved.
 
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met: 
    * Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.  
    * Redistributions in binary form must reproduce the above
  copyright notice, this list of conditions and the following
  disclaimer in the documentation and/or other materials provided
  with the distribution.
    * Neither the name of the <ORGANIZATION> nor the names of its
  contributors may be used to endorse or promote products derived
  from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gtmma.h"
#include "parser.h"
#include "complex.h"
#include "misc_operations.h"
#include "general_acoustics.h"
#include "motionless_skel_models.h"

/*
 * MAIN FUNCTION
 *
 */
int
main ( int argc, 
       char** argv ) {

  struct aProblem* problem;
  struct aComplex* rho;
  struct aComplex* K;
  struct aComplex* Zc;
  struct aComplex* k;

  struct aComplex* Z;

  struct aComplexMatrix* fluid_matrix_1;
  struct aComplexMatrix* fluid_matrix_2;
  struct aComplexMatrix* merge_fluid_matrices;

  double Z0;
  double omega;
  double alpha;

  int nb_frequency_points;


  /* File stuff */
  FILE* pfile;
  char* parameters_file;
  char* file_name;
  char* file_name_out;


  /* Counters and temporary variables */
  int i, j;
  int tmp;





  /* Parse command line */
  if (argc == 1) {
    printf("Usage: gtmma parameters_file.in\n");
    return 1;
  } else {
    tmp = strlen(argv[1]);
    parameters_file = malloc ( sizeof(char) * (tmp+1) );
    strncpy(parameters_file, argv[1], sizeof(char) * (tmp+1) );
    parameters_file[tmp] = '\0';
  }

  tmp = strlen(parameters_file)-3;
  file_name = malloc ( sizeof(char) * (tmp+1) );
  strncpy ( file_name, parameters_file, tmp );
  file_name[tmp] = '\0';

  file_name_out = malloc ( sizeof(char) * (tmp+4+1) );
  strncpy ( file_name_out, file_name, tmp+1 );
  strncat ( file_name_out, ".out", 4);
  file_name_out[tmp+4] = '\0';


  /* Memory allocation */
  rho = (struct aComplex*) malloc ( sizeof (struct aComplex) );
  K   = (struct aComplex*) malloc ( sizeof (struct aComplex) );
  Zc  = (struct aComplex*) malloc ( sizeof (struct aComplex) );
  k   = (struct aComplex*) malloc ( sizeof (struct aComplex) );
  Z   = (struct aComplex*) malloc ( sizeof (struct aComplex) );

  problem                  = (struct aProblem*) malloc ( sizeof(struct aProblem) );
  problem->fluid           = (struct aFluid*) malloc( sizeof(struct aFluid) );
  problem->spectrum        = (struct aSpectrum*) malloc( sizeof(struct aSpectrum) );
  problem->mat             = (struct aMaterial**) malloc( 1 * sizeof(struct aMaterial*) );
  problem->nb_of_materials = 0;
  problem->conditions      = (struct aSetOfConditions*) malloc( sizeof(struct aSetOfConditions) );

  problem_initialize( problem );

  fluid_matrix_1 = (struct aComplexMatrix*) malloc( sizeof(struct aComplexMatrix) );
  fluid_matrix_1->data =  (struct aComplex**) malloc( 2 * sizeof(struct aComplex*) );
  fluid_matrix_2 = (struct aComplexMatrix*) malloc( sizeof(struct aComplexMatrix) );
  fluid_matrix_2->data =  (struct aComplex**) malloc( 2 * sizeof(struct aComplex*) );
  merge_fluid_matrices = (struct aComplexMatrix*) malloc( sizeof(struct aComplexMatrix) );
  merge_fluid_matrices->data =  (struct aComplex**) malloc( 2 * sizeof(struct aComplex*) );
  for ( i = 0; i < 2; i++ ) {
    fluid_matrix_1->data[i] = (struct aComplex*) malloc( 2 * sizeof(struct aComplex) );
    fluid_matrix_2->data[i] = (struct aComplex*) malloc( 2 * sizeof(struct aComplex) );
    merge_fluid_matrices->data[i] = (struct aComplex*) malloc( 2 * sizeof(struct aComplex) );
  }
  fluid_matrix_1->nb_rows = 2;
  fluid_matrix_1->nb_columns = 2;
  fluid_matrix_2->nb_rows = 2;
  fluid_matrix_2->nb_columns = 2;
  merge_fluid_matrices->nb_rows = 2;
  merge_fluid_matrices->nb_columns = 2;


  /* Pre-processing:
   * Read parameters file and compute some properties */

  read_parameters ( parameters_file, problem );

  nb_frequency_points = get_nb_frequency_points ( problem->spectrum );
  air_properties( problem->fluid );
  Z0 = problem->fluid->rho_0*problem->fluid->c_0;

  /* FREQUENCY LOOP */
  for ( i = 0; i < nb_frequency_points; i++ ) {
    omega = 2*PI* ( problem->spectrum->min + i*problem->spectrum->step );

    if ( problem->nb_of_materials > 0 ) {
      for ( j = 0; j < problem->nb_of_materials; j++ ) {
	switch( problem->mat[j]->type ) {
	case AIR_GAP:
	  Zc->re = Z0;
	  Zc->im = 0.0;
	  k->re  = omega/problem->fluid->c_0;
	  k->im  = 0.0;
	  break;
	case DELANY_BAZLEY:
	  DB_Zc_k( rho, K, Zc, k, omega, problem->fluid, problem->mat[j]->sigma );
	  break;
	case GARAI_POMPOLI:
	  GP_Zc_k( rho, K, Zc, k, omega, problem->fluid, problem->mat[j]->sigma );
	  break;
	case JOHNSON_CHAMPOUX_ALLARD:
	  JCA_Zc_k( rho, K, Zc, k, omega, problem->fluid, problem->mat[j] );
	  break;
	case JOHNSON_LAFARGE:
	  JL_Zc_k( rho, K, Zc, k, omega, problem->fluid, problem->mat[j] );
	  break;
	default:
	  printf("mat nb: %d, total nb mat: %d, type: %d\n",j+1, problem->nb_of_materials,problem->mat[j]->type );
	  printf("Unknown material type !\n");
	  exit(0);
	  break;
	}
	build_fluid_matrix( *Zc, *k, problem->fluid, problem->conditions->angle, problem->mat[j]->L, omega, fluid_matrix_2 );
	if ( j == 0 ) { /*problem->nb_of_materials == 1 ) {*/
	  /* Build a 2x2 identity matrix */
	  fluid_matrix_1->data[0][0].re = 1.;
	  fluid_matrix_1->data[0][0].im = 0.;
	  fluid_matrix_1->data[0][1].re = 0.;
	  fluid_matrix_1->data[0][1].im = 0.;
	  fluid_matrix_1->data[1][0].re = 0.;
	  fluid_matrix_1->data[1][0].im = 0.;
	  fluid_matrix_1->data[1][1].re = 1.;
	  fluid_matrix_1->data[1][1].im = 0.;
	}
	cpx_matrices_mult( fluid_matrix_2, fluid_matrix_1, merge_fluid_matrices );  
	cpx_copy_matrix( fluid_matrix_1, merge_fluid_matrices );
      }
    }


    backed_layer_impedance( fluid_matrix_1, Z, k );


    /*TODO: switch( probem->conditions->diffuse_field ) {
     * case STATISTICAL:
     * case LOCALLY_REACTING: 
     */
    if ( problem->conditions->diffuse_field != 0 ) {
      alpha = compute_alpha_diffuse_london( Z, Z0 );
    }
    else {
      alpha = compute_alpha( Z, Z0, problem->conditions->angle );
    }

    if ( i == 0 ) {
      pfile = fopen( file_name_out, "w");
      fprintf( pfile, "%% Frequency Re(Rho) Im(Rho) Re(K) Im(K) alpha\n" );
      fclose( pfile );
    }

    pfile = fopen( file_name_out, "a" );
    fprintf( pfile, "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", 
	     omega/(2*PI), 
	     rho->re/problem->fluid->rho_0, 
	     rho->im/problem->fluid->rho_0, 
	     K->re/problem->fluid->P_0, 
	     K->im/problem->fluid->P_0, 
	     alpha );
    fclose( pfile );

  }
  /* END OF FREQUENCY LOOP */


  /* Before leaving, free all dynamic memories */
  free( parameters_file );
  free( file_name );
  free( file_name_out );
  
  for ( i = 0; i < 2; i++ ) {
    free( fluid_matrix_1->data[i] );
    free( fluid_matrix_2->data[i] );
    free( merge_fluid_matrices->data[i] );
  }
  free( fluid_matrix_1 );
  free( fluid_matrix_2 );
  free( merge_fluid_matrices );

  free( problem->fluid );
  for ( i = 0; i < problem->nb_of_materials; i++ ) {
    free( problem->mat[i] );
  } 
  free( problem->mat );
  free( problem->spectrum );
  free( problem->conditions );
  free( problem );

  free( rho );
  free( K );
  free( Z );
  free( k );
  free( Zc );

  return 0;

}
