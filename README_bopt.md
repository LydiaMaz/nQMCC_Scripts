# BoundOpt :: Bound State Parameter Optimization for nQMCC

This script automates the optimization of variational parameters for bound state calculations in nuclear Quantum Monte Carlo with Continuum (nQMCC).


## Usage

To run bound state optimization:

```bash
python bopt.py --utility /absolute/path/to/file.util
```

> ✅ Make sure your `utility` file includes the following fields:
- `ESEP_SCALE`: Fraction of Current ESEP in Decks to set bounds for optimization (generally smaller when using alpha-seeded deck)
- `OPT_SCALE`: scaling for 2b and 3b correlations
    - flat opt_scale value used for 3b correlations (delta, epsilon, theta, upsilon, rscal, uscal) since these are typically 0 in initial decks
    - scale opt_scale value used for 2b correlations (fscal, alpha, beta) as fraction of current non-zero values in deck 

> ⚠️ `OPTIONAL PARAMETERS`:
- `ALPHA_DIR`: Path to directory containing pre-optimized alpha-core parameters
  - If provided, enables alpha-core seeding for improved convergence
  - If not found or inaccessible, the script continues with unseeded optimization

Outputs include optimized nQMCC decks (`.dk`), and a summary of results with energy improvements and optimization details; saved in json format to the working directory.

