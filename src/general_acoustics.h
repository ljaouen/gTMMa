/* This file is part of gTMMa.
 * Copyright (c) 2004, 2013, Luc Jaouen 
 * under the terms of the new BSD License. 
 */

void 
air_properties( struct aFluid* );

void
Zc_k_from_rho_K( struct aComplex*,
		 struct aComplex*,
		 struct aComplex*,
		 struct aComplex*,
		 double            );

void
backed_layer_impedance( struct aComplexMatrix*,
			struct aComplex*, 
			struct aComplex*       );

double
compute_alpha( struct aComplex*, 
	       double,
	       double  );

double
compute_alpha_diffuse_london( struct aComplex*,
			      double );

void
one_third_octave( double*,
		  double*,
		  int,
		  /*  double one_third_freq[12], */
		  double[12] );

void 
build_fluid_matrix( struct aComplex,
		    struct aComplex,
		    struct aFluid*,
		    double,
		    double,
		    double,
		    struct aComplexMatrix* );

