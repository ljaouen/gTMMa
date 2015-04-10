gTMMa allows to predict the sound absorption properties of porous materials such as polymer foams or mineral wools.
gTMMa is written in ANSI C and is designed to run under multiple OS.

# License

This code is distributed under the terms of the new BSD License.
Copyleft Luc Jaouen.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    The names of the authors may not be used to endorse or promote products derived from this software without specific prior written permission.

# Installing gTMMa

Unzip gtmma.zip.
4 sub-directories should be created:
- doc contains the documentation of gTMMa,
- examples contains some input files (.in),
- bin contains pre-compiled versions of gTMMa (a Windows XP one and a Linux kernel 2.4/libc6 one).
- src contains the ANSI C sources of gTMMa.


# Using gTMMa

At this time, only a console version of gTMMa is available: you need to open a console to lunch gTMMa (button start/execute and type "cmd" for Windows).
Change directory to the directory where the binary "gtmma" (lower case) is or add this directory to your path.
Simply type "gtmma path_to_an_input_file" to compute the response for an input file. (example: gtmma JCA_gap_JCA.in).
Input file example

Each input file (identified by a ".in" extension) must contain 4 tags: a fluid, a layers, a spectrum and a conditions one.

% Datafile for gTMMa Headers/comments
%
% 2004.10.09

<fluid>
    T  20 Room temperature T in Celsius degrees
    P  101300 Atmospheric static pressure P in Pascal
</fluid>

<layers>
  <ridig_impervious_wall> Required in "absorption mode"
  <JCA_material> Beginning of the 1st layer tag (Johnson Champoux Allard mat.)
    alpha  2.1       High frequency limit of the tortuosity
    phi  0.80        Open porosity
    sigma  40000     Static air flow resistivity in N.s.m-4
    lambda  125e-06  Viscoous characteristic length in meters
    lambda_p 350e-06 Thermal characteristic length in meters
    L  0.01          Layer thickness in meters
  </JCA_material>    End of the layer tag

  <air_gap> Beginning of the 2nd layer tag (air gap/plenum)
  L 0.02    Layer thickness in meters
  </air_gap>
 
  <JCA_material> 3rd layer: Johnson Champoux Allard material
    alpha  1.6
    phi 0.82
    sigma  20000        (N.s.m-4)
    lambda  80e-06      (m)
    lambda_p  160e-06   (m)
    L  0.01             (m)
  </JCA_material>
</layers>

<spectrum>
    min    10.0 Lower frequency in Hertz
    max  4000.0 Upper frequency in Hertz
    step   10.0 Frequency step in Hertz
</spectrum>

<conditions>
  incidence 0.0   Incidence angle in radian (here: normal incidence)
</conditions>

# Material tags allowed

<air_gap>: an air gap/plenum. The only parameter required is the layer thickness L.
<DB_material>: Delany-Bazley model. Parameters are the static air flow resistivity sigma (in N.s.m-4) and the layer thickness L.
<GP_material>: Garai-Pompoli model. Parameters are the static air flow resistivity sigma (in N.s.m-4) and the layer thickness L.
<JCA_material>: Johnson-Champoux-Allard model. Parameters are the static air flow resistivity sigma (in N.s.m-4), the open porosity phi, the high frequency limit of the tortuosity alpha, the viscous and thermal characteristic lengths (in meters) lambda and lambda_p respectively and the layer thickness L.
<JL_material>: Johnson-Champoux-Allard-Lafarge model. Parameters are the static air flow resistivity sigma (in N.s.m-4), the open porosity phi, the high frequency limit of the tortuosity alpha, the viscous and thermal characteristic lengths (in meters) lambda and lambda_p respectively, the thermal permeability k_p_0 (in m^2) and finally the layer thickness L.

