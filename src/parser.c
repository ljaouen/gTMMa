/* This file is part of gTMMa.
 * Copyright (c) 2004, 2013, Luc Jaouen 
 * under the terms of the new BSD License. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gtmma.h"
#include "misc_operations.h"

/*
static void 
error_msg ( char* errmsg ) {
  fprintf(stderr,"Oops : %s\n",errmsg);
  return;
}
*/

/* Parser for parameters_file ('.in')  */
void
read_parameters ( char* parameters_file,
		  struct aProblem* problem )
{
  
  struct aMaterial** tmp_problem_mat = NULL;
  
  double tmp_double;

  char tmp_char[50] = "\0";

  FILE* pfile;
  char line[255];
  char header[255];
  int end_tag = 0;

  int idx;

  pfile = fopen( parameters_file, "r");

  do {
    fgets( line, 254, pfile ); sscanf(line, "%s", header);

    /* fluid tag */
    if ( strncmp("<fluid>",header,7) == 0 ) {
      printf("Find tag fluid\n");
      do {
	fgets( line, 254, pfile ); sscanf(line, "%s", header);
	if ( strncmp("</fluid>",header,8) == 0 ) end_tag = 1;
	else if ( sscanf(line, "%s %lf", tmp_char, &tmp_double) == 2 ) {
	  if ( strncmp("T",tmp_char,1) == 0 ) problem->fluid->T = tmp_double;
	  else if ( strncmp("P",tmp_char,1) == 0 ) problem->fluid->P_0 = tmp_double;
	  else if ( strncmp("H",tmp_char,1) == 0 ) problem->fluid->H = tmp_double;
	}
      }
      while ( !feof(pfile) && !end_tag );
      end_tag = 0;
    }

    /* DB material tag */
    else if ( strncmp("<DB_material>",header,13) == 0 ) {
      printf("Find tag DB_material\n");
      
      problem->nb_of_materials += 1;
      idx = problem->nb_of_materials;
      if ((tmp_problem_mat = realloc( (struct aMaterial**)problem->mat, idx * sizeof(struct aMaterial*) )) == NULL) {
	printf("Re-allocation error of problem->mat in parser.c (DB mat tag)\n");
      }
      idx--;
      tmp_problem_mat[idx] = (struct aMaterial*) malloc( sizeof(struct aMaterial) );
      problem->mat = tmp_problem_mat;

      material_initialize( problem->mat[idx] );

      printf("total_nb mat = %d\n",problem->nb_of_materials);

      problem->mat[idx]->type = DELANY_BAZLEY;

      do {
	fgets( line, 254, pfile ); sscanf(line, "%s", header);
	if ( strncmp("</DB_material>",header,14) == 0 ) end_tag = 1;
	else if ( sscanf(line, "%s %lf", tmp_char, &tmp_double) == 2 ) {
	  if ( strncmp("sigma",tmp_char,5) == 0 ) problem->mat[idx]->sigma = tmp_double;
	  else if ( strcmp("L",tmp_char) == 0 ) problem->mat[idx]->L = tmp_double;
	}
      }
      while ( !feof(pfile) && !end_tag );
      end_tag = 0;
    }


    /* DB material tag */
    else if ( strncmp("<GP_material>",header,13) == 0 ) {
      printf("Find tag GP_material\n");
      
      problem->nb_of_materials += 1;
      idx = problem->nb_of_materials;
      if ((tmp_problem_mat = realloc( (struct aMaterial**)problem->mat, idx * sizeof(struct aMaterial*) )) == NULL) {
	printf("Re-allocation error of problem->mat in parser.c (GP mat tag)\n");
      }
      idx--;
      tmp_problem_mat[idx] = (struct aMaterial*) malloc( sizeof(struct aMaterial) );
      problem->mat = tmp_problem_mat;

      material_initialize( problem->mat[idx] );

      printf("total_nb mat = %d\n",problem->nb_of_materials);

      problem->mat[idx]->type = GARAI_POMPOLI;

      do {
	fgets( line, 254, pfile ); sscanf(line, "%s", header);
	if ( strncmp("</GP_material>",header,14) == 0 ) end_tag = 1;
	else if ( sscanf(line, "%s %lf", tmp_char, &tmp_double) == 2 ) {
	  if ( strncmp("sigma",tmp_char,5) == 0 ) problem->mat[idx]->sigma = tmp_double;
	  else if ( strcmp("L",tmp_char) == 0 ) problem->mat[idx]->L = tmp_double;
	}
      }
      while ( !feof(pfile) && !end_tag );
      end_tag = 0;
    }

    /* JCA material tag */
    else if ( strncmp("<JCA_material>",header,14) == 0 ) {
      printf("Find tag JCA_material\n");
      
      problem->nb_of_materials = problem->nb_of_materials + 1;
      idx = problem->nb_of_materials;
      if ((tmp_problem_mat = realloc( (struct aMaterial**)problem->mat, idx * sizeof(struct aMaterial*) )) == NULL) {
	printf("Re-allocation error of problem->mat in parser.c (JCA mat tag)\n");
      }
      idx--;
      tmp_problem_mat[idx] = (struct aMaterial*) malloc( sizeof(struct aMaterial) );
      problem->mat = tmp_problem_mat;
      /*
      free( tmp_problem_mat[idx] );
      free( tmp_problem_mat );
      */
      /*
      for ( i = 0; i < problem->nb_of_materials-1; i++ ) {
	free( tmp_problem_mat[i] );
      } 
      free( tmp_problem_mat );
      */
      material_initialize( problem->mat[idx] );

      printf("total_nb mat = %d\n",problem->nb_of_materials);

      problem->mat[idx]->type = JOHNSON_CHAMPOUX_ALLARD;

      do {
	fgets( line, 254, pfile ); sscanf(line, "%s", header);
	if ( strncmp("</JCA_material>",header,15) == 0 ) end_tag = 1;
	else if ( sscanf(line, "%s %lf", tmp_char, &tmp_double) == 2 ) {
	  if ( strncmp("phi",tmp_char,3) == 0 ) problem->mat[idx]->phi = tmp_double;
	  else if ( strncmp("sigma",tmp_char,5) == 0 ) problem->mat[idx]->sigma = tmp_double;
	  else if ( strncmp("alpha",tmp_char,5) == 0 ) problem->mat[idx]->alpha_infty = tmp_double;
	  else if ( strcmp("lambda",tmp_char) == 0 )   problem->mat[idx]->lambda = tmp_double;
	  else if ( strncmp("lambda_p",tmp_char,8) == 0 ) problem->mat[idx]->lambda_p = tmp_double;
	  else if ( strcmp("L",tmp_char) == 0 ) problem->mat[idx]->L = tmp_double;
	}
      }
      while ( !feof(pfile) && !end_tag );
      end_tag = 0;
    }

    /* JCAL material tag */
    else if ( strncmp("<JL_material>",header,13) == 0 ) {
      printf("Find tag JL_material\n");
      
      problem->nb_of_materials = problem->nb_of_materials + 1;
      idx = problem->nb_of_materials;
      if ((tmp_problem_mat = realloc( (struct aMaterial**)problem->mat, idx * sizeof(struct aMaterial*) )) == NULL) {
	printf("Re-allocation error of problem->mat in parser.c (JCA mat tag)\n");
      }
      idx--;
      tmp_problem_mat[idx] = (struct aMaterial*) malloc( sizeof(struct aMaterial) );
      problem->mat = tmp_problem_mat;

      material_initialize( problem->mat[idx] );

      printf("total_nb mat = %d\n",problem->nb_of_materials);

      problem->mat[idx]->type = JOHNSON_LAFARGE;

      do {
	fgets( line, 254, pfile ); sscanf(line, "%s", header);
	if ( strncmp("</JL_material>",header,13) == 0 ) end_tag = 1;
	else if ( sscanf(line, "%s %lf", tmp_char, &tmp_double) == 2 ) {
	  if ( strncmp("phi",tmp_char,3) == 0 ) problem->mat[idx]->phi = tmp_double;
	  else if ( strncmp("sigma",tmp_char,5) == 0 ) problem->mat[idx]->sigma = tmp_double;
	  else if ( strncmp("alpha",tmp_char,5) == 0 ) problem->mat[idx]->alpha_infty = tmp_double;
	  else if ( strcmp("lambda",tmp_char) == 0 )   problem->mat[idx]->lambda = tmp_double;
	  else if ( strncmp("lambda_p",tmp_char,8) == 0 ) problem->mat[idx]->lambda_p = tmp_double;
	  else if ( strncmp("k_p_0",tmp_char,5) == 0 ) problem->mat[idx]->k_p_0 = tmp_double;
	  else if ( strcmp("L",tmp_char) == 0 ) problem->mat[idx]->L = tmp_double;
	}
      }
      while ( !feof(pfile) && !end_tag );
      end_tag = 0;
    }

    /* Air gap tag */
    else if ( strncmp("<air_gap>",header,14) == 0 ) {
      printf("Find tag air gap\n");
      
      problem->nb_of_materials = problem->nb_of_materials + 1;
      idx = problem->nb_of_materials;
      if ((tmp_problem_mat = realloc( (struct aMaterial**)problem->mat, idx * sizeof(struct aMaterial*) )) == NULL) {
	printf("Re-allocation error of problem->mat in parser.c (air gap tag)\n");
      }
      idx--;
      tmp_problem_mat[idx] = (struct aMaterial*) malloc( sizeof(struct aMaterial) );
      problem->mat = tmp_problem_mat;
      material_initialize( problem->mat[idx] );

      /*
      problem->nb_of_materials = problem->nb_of_materials + 1;
      realloc( (struct aMaterial**)problem->mat, problem->nb_of_materials * sizeof(struct aMaterial*) );
      idx = problem->nb_of_materials - 1;
      problem->mat[idx] = (struct aMaterial*) malloc( sizeof(struct aMaterial) );
      */

      printf("total_nb mat = %d\n",problem->nb_of_materials);

      problem->mat[idx]->type = AIR_GAP;

      do {
	fgets( line, 254, pfile ); sscanf(line, "%s", header);
	if ( strncmp("</air_gap>",header,15) == 0 ) end_tag = 1;
	else if ( sscanf(line, "%s %lf", tmp_char, &tmp_double) == 2 ) {
	  if ( strcmp("L",tmp_char) == 0 ) problem->mat[idx]->L = tmp_double;
	}
      }
      while ( !feof(pfile) && !end_tag );
      end_tag = 0;
    }

    /* Conditions tag */
    else if ( strncmp("<conditions>",header,12) == 0 ) {
      printf("Find tag conditions\n");
      do {
	fgets( line, 254, pfile ); sscanf(line, "%s", header);
	if ( strncmp("</conditions>",header,13) == 0 ) end_tag = 1;
	else if ( sscanf(line, "%s %lf", tmp_char, &tmp_double) == 2 ) {
	  if ( strcmp("incidence",tmp_char) == 0 ) {
	    problem->conditions->diffuse_field = 0;
	    problem->conditions->angle = tmp_double*PI/180;
	  }
	}
	else if ( sscanf(line, "%s %lf", tmp_char, &tmp_double) == 1 ) {
	  if ( strcmp("diffuse_field",tmp_char) == 0 ) {
	    problem->conditions->diffuse_field = 1;
	    problem->conditions->angle = 0.0;
	    printf("Diffuse field condition\n");
	  }
	}
      }
      while ( !feof(pfile) && !end_tag );
      end_tag = 0;
    }

    /* Spectrum tag */
    else if ( strncmp("<spectrum>",header,7) == 0 ) {
      printf("Find tag spectrum\n");
      do {
	fgets( line, 254, pfile ); sscanf(line, "%s", header);
	if ( strncmp("</spectrum>",header,8) == 0 ) end_tag = 1;
	else if ( sscanf(line, "%s %lf", tmp_char, &tmp_double) == 2 ) {
	  if ( strncmp("min",tmp_char,3) == 0 ) problem->spectrum->min = tmp_double;
	  else if ( strncmp("max",tmp_char,3) == 0 ) problem->spectrum->max = tmp_double;
	  else if ( strncmp("step",tmp_char,3) == 0 ) problem->spectrum->step = tmp_double;
	}
      }
      while ( !feof(pfile) && !end_tag );
      end_tag = 0;
    }
  } while ( !feof(pfile) );

  fclose(pfile);

}
