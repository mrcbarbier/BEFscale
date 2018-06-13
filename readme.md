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

 *Change other parameters*

    python landscapesimu.py tmax=500 tsample=10


## Where to edit:

   Parameters and values to iterate over in axes: landscapesimu.py
   Plots: plots.py