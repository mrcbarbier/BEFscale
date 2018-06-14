Simulation experiments for multi-scale ecological interactions.
=======

## Requirements:
 Python 2.7

 Libraries: numpy, scipy, pandas, matplotlib

## Examples:

 *Basic run*

    python landscapesimu.py

 will create folder ./RESULTS/beflandscape and store simulation results in it, then display various plots.

 *Save results in folder ./RESULTS/xxx*

    python landscapesimu.py path=xxx

 *Rerun simulations instead of opening saved results*

    python landscapesimu.py rerun

 *Run with/out multiscale interactions (without is faster, useful for tests)* (NB: brackets = optional/alternative)

    python landscapesimu.py [multiscale=0]  # (None)
    python landscapesimu.py multiscale[=all]  # (All interactions)
    python landscapesimu.py multiscale=[trophic,competition]  # (Only one interaction)

 *Change other parameters (simulation runtime, number of snapshots taken)*

    python landscapesimu.py tmax=500 nsample=50

 *Plot abundances over time in plots/xxx/movie folder (to make an animation)*

    python landscapesimu.py movie


## Where to edit:

   Default parameter values:  default_parameters.dat

   Parameter values that are iterated over:  loop_parameters.dat

   Plots: plots.py