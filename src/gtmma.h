/* This file is part of gTMMa.
 * Copyright (c) 2004, 2013, Luc Jaouen 
 * under the terms of the new BSD License. 
 */

#define PI 3.1415926535897932385

/* Definition of a frequency spectrum */
struct aSpectrum{
  double min;
  double max;
  double step;
};

/* Parameters of a fluid */
struct aFluid {
  double T;
  double P_0;
  double H;
  double c_0;
  double rho_0;
  double Z_0;
  double eta;
  double gamma;
  double Cp;
  double kappa;
  double Pr;
};

typedef enum {
  AIR_GAP,
  DELANY_BAZLEY,
  GARAI_POMPOLI,
  JOHNSON_CHAMPOUX_ALLARD,
  JOHNSON_LAFARGE
} enumMaterialType;

/* Parameters for all kinds of materials */
struct aMaterial {
  enumMaterialType type;
  char   name[40];
  double L;
  double phi;
  double sigma;
  double alpha_infty;
  double lambda;
  double lambda_p;
  double k_p_0;
};

/*
struct allMaterials {
  int total_nb;
  struct aMaterial** mat;
};
*/
struct aSetOfConditions {
  int diffuse_field;
  double angle;
};

struct aProblem {
  struct aSpectrum* spectrum;
  struct aFluid* fluid;
  struct aSetOfConditions* conditions;
  int nb_of_materials;
  struct aMaterial** mat;
};
