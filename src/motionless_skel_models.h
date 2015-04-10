/* This file is part of gTMMa.
 * Copyright (c) 2004, 2013, Luc Jaouen 
 * under the terms of the new BSD License. 
 */

void
DB_Zc_k ( struct aComplex*,
	  struct aComplex*,
	  struct aComplex*,
	  struct aComplex*,
	  double,
	  struct aFluid*,
	  double           );

void
GP_Zc_k ( struct aComplex*,
	  struct aComplex*,
	  struct aComplex*,
	  struct aComplex*,
	  double,
	  struct aFluid*,
	  double           );

void
JCA_Zc_k ( struct aComplex*,
	   struct aComplex*,
	   struct aComplex*,
	   struct aComplex*,
	   double,
	   struct aFluid*,
	   struct aMaterial* );

void
JL_Zc_k ( struct aComplex*,
	    struct aComplex*,
	    struct aComplex*,
	    struct aComplex*,
	    double,
	    struct aFluid*,
	    struct aMaterial* );

struct aComplex 
function_Gj ( double,
	      double,
	      double,
	      double,
	      double,
	      double,
	      double );
