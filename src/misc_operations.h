/* This file is part of gTMMa.
 * Copyright (c) 2004, 2013, Luc Jaouen 
 * under the terms of the new BSD License. 
 */

int 
get_nb_frequency_points ( struct aSpectrum* );

void
problem_initialize( struct aProblem* );

void 
material_initialize( struct aMaterial* );

void
gauss_legendre_rule( int, 
		     double*, 
		     double* );

void
find_freq_points ( double,
		   double,
		   double*,
		   int,
		   int*,
		   int* );
