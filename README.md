# Stratified Three-Phase Waves

A numerical framework for simulating three-phase stratified interfacial flows using the Basilisk flow solver.

## Overview

This code simulates three immiscible fluid layers with the following configuration:
- Fluid 1 (bottom): Water (ocean)
- Fluid 2 (middle): Oil layer that forms a precursor film
- Fluid 3 (top): Air

The key feature is that Fluid 2 (oil) always maintains a minimum thickness of one grid cell (when we call it a precursor film, it is the worst case scenario) between Fluids 1 and 3.

## Key Files

### Core Simulation Files
- `stratified.c`: Main simulation file that sets up the three-phase flow configuration
- `three-phase-stratified.h`: Header file implementing the three-phase flow solver with Volume-Of-Fluid method
- `reduced-three-phase-stratified.h`: Implements reduced gravity formulation for the three-phase system

### Post-processing Tools
- `getData-threePhase.c`: Extracts field data from simulation snapshots
- `getFacet-threePhase.c`: Extracts interface positions
- `Video.py`: Python script for visualization and video generation

## Key Parameters

The main simulation parameters can be found in `stratified.c`:
- `OhBulk`: Ohnesorge number for bulk fluid (water)
- `mu21`, `mu31`: Viscosity ratios
- `rho21`, `rho31`: Density ratios
- `sigma12`, `sigma23`: Surface tension coefficients
- `Bo`: Bond number
- `MAXlevel`: Maximum grid refinement level

## Implementation Details

The three-phase implementation uses two Volume-Of-Fluid fields:
- `f1`: Distinguishes between fluids 1+2 (value 1) and fluid 3 (value 0)
- `f2`: Distinguishes between fluid 1 (value 1) and fluids 2+3 (value 0)

The property indices for each fluid are:
- Fluid 1 (water): `f1*f2`
- Fluid 2 (oil): `f1*(1-f2)`
- Fluid 3 (air): `1-f1`



## Running the Code

1. Compile the main simulation:

- Serial
```bash
qcc -O2 -Wall stratified.c -o stratified -lm
./stratified
```

- Using OpenMP
```bash
qcc -O2 -Wall stratified.c -o stratified -lm -fopenmp
export OMP_NUM_THREADS=8
./stratified
```

Here, 8 is the number of threads to use. 

- Using OpenMPI
```bash
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions stratified.c -o stratified -lm 
mpirun -np 8 stratified
```

Here, 8 is the number of MPI processes to use.

Post processing: 

```bash
python Video.py --CPUs 4
```

Here, 4 is the number of processors to use.
