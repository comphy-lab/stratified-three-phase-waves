/** Title: stratified.c
 * waves on a stratified interface
# Authors: Ilies Haouche & Vatsal Sanjay
# Version 0.1
# Last Updated: Nov 5 2024

# Changelog: 
## v0.1 (Nov 5, 2024)
- introduce a stratified layer such that:
  - fluid 1 (property index): (f1*f2). This is the water (ocean).
  - fluid 2 (property index):(f1*(1-f2)) forms a precursor film. This is the oil layer.
  - fluid 3 (property index):(1-f1). this is air. 
- **Note:** Precursor layer means that even in the worst case, the minimum thickness of fluid 2 (f1*(1-f2)) will be one grid cell. In general it will form a macroscopic thick layer between fluids 1 and 3.
*/

// 1 is coating, 2 is bulk and 3 is air
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "three-phase-stratified.h"
#include "tension.h"
#include "reduced-three-phase-stratified.h"

#define MINlevel 3                                              // maximum level

#define tsnap (1e-2)

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-2)                            // error tolerances in velocity

// boundary conditions

double OhBulk, mu21, mu31, rhoBulk, rho21, rho31, Hf, sigma12, sigma23, Bo, tmax, Ldomain;
int MAXlevel;

int main(int argc, char const *argv[]) {
  
  /**
   * The parameters below can be filled in here or passed in from the command line after necessary changes. 
  */

  /*
  Be careful with zero viscosity in one-fluid approximation. See: [[2406.05416] Focusing of concentric free-surface waves](https://arxiv.org/abs/2406.05416)
  */
  // bulk is fluid 1
  OhBulk = 0.01; // atof(argv[1]) // this is the dimesionless viscosity of bulk fluid 1 (water): \mu/sqrt(\rho \gamma L)
  mu21 = 1.0; // atof(argv[2]) // \mu_2/\mu_1
  mu31 = 1e-2; // atof(argv[3]) // \mu_3/\mu_1
  
  rhoBulk = 1e0; // this is the density of fluid 1 and is the repeating variable in the normalization. 
  rho21 = 1e0; // atof(argv[4]) // density ratio: \rho_2/\rho_1
  rho31 = 1e-3; // atof(argv[5]) // density ratio: \rho_3/\rho_1

  sigma12 = 1.0; // this is the surface tension coefficient between fluid 1 (bulk) and 2 (oil layer).. this is dimensionless
  sigma23 = 1.0; // this is the surface tension coefficient between fluid 2 (oil layer) and 3 (air).. this is dimensionless

  Bo = 0.1; // this is the dimensionless gravity: \rho g L^2/\sigma (where L is the relevant length scale)

  tmax = 1e0; // atof(argv[6]) // total simulation time
  Ldomain = 4e0; // atof(argv[7]) // domain size
  MAXlevel = 8; // atof(argv[8]) // maximum level of refinement

  Hf = 0.1; // atof(argv[9]) // height of fluid 2

  L0=Ldomain;
  X0=-L0/2.; Y0=0.; // 0,0 is the mean height of fluid 1

  init_grid (1 << (8));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = rhoBulk; mu1 = OhBulk;
  rho2 = rho21*rhoBulk; mu2 = mu21*OhBulk;
  rho3 = rho31*rhoBulk; mu3 = mu31*OhBulk; 

  Bf1.x = -Bo; Bf2.x = -Bo;

  /*
  The two lines below are the main considerations in deciding the configuration. 
  - f1 = 1 for fluids 1 (water) and 2 (oil), and 0 for fluid 3. so the surface tension coefficient associated with f1 is that of oil - air
  - And f2 = 1 for fluid 1 and 0 for fluids 2 and 3. so the interface coefficient associated with f2 is that of oil - water. 
  */

  f1.sigma = sigma23; // surface tension oil-air!!!
  f2.sigma = sigma12; // interfacial tension water-oil!!!

  fprintf(ferr, "Level %d tmax %g. OhBulk %g, mu21 %g, mu31 %g, rhoBulk %g, rho21 %g, rho31 %g, sigma12 %g, sigma23 %g\n", MAXlevel, tmax, OhBulk, mu21, mu31, rhoBulk, rho21, rho31, sigma12, sigma23);
  run();

}

event init(t = 0){
  if(!restore (file = "dump")){
    /**
    Initialization for f1 and f2
    */
    fraction (f2, -x + 0.05*cos(2*pi*y/L0)); // #TODO: Ilies, this needs to be changed
    fraction (f1, -(x-Hf) + 0.05*cos(2*pi*y/L0)); // #TODO: Ilies, this needs to be changed
  }
  // dump (file = "dump");
  // return 1;
}

scalar KAPPA1[], KAPPA2[];
event adapt(i++) {
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);

  adapt_wavelet ((scalar *){f1, f2, u.x, u.y, KAPPA1, KAPPA2},
    (double[]){fErr, fErr, VelErr, VelErr, KErr, KErr},
    MAXlevel, MINlevel);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += sq(Delta)*(sq(u.x[]) + sq(u.y[]))*rho(f1[],f2[]);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf(fp, "Level %d tmax %g. OhBulk %g, mu21 %g, mu31 %g, rhoBulk %g, rho21 %g, rho31 %g, sigma12 %g, sigma23 %g\n", MAXlevel, tmax, OhBulk, mu21, mu31, rhoBulk, rho21, rho31, sigma12, sigma23);
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
}