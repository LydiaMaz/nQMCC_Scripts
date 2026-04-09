# AutoOpt :: A Parameter Optimization Framework for nQMCC

This project automates the optimizations of variational parameters for nuclear Quantum Monte Carlo with Continuum (nQMCC).

## Maintainer
Lydia Mazeeva (WASHU)

## AutoOpt Contributors
Lydia Mazeeva (WASHU)\
Abraham R. Flores (WASHU)\
Juan Jose Silva (WASHU)\
Maria Piarulli (WASHU)

## Todo
1. Generate zero deck main program
2. Improve bounds checking in scan
3. Update logic for +1 node
4. Update Documentation
5. Automate tests

## Usage

To run a full optimization:

```bash
python /absolute/path/to/AutoOpt.py --utility /absolute/path/to/file.util
```
> ✅ Make sure your `utility` file includes the following fields:
- `NAME`: label for outputs
- `SYSTEM_TYPE`: bound, sc_scattering
- `RUN_CMD`: prefix to run binary (mpirun -np X)
- `NQMCC_DIR`: /path/to/nQMCC
- `BIN_DIR`: /path/to/nQMCC/bin
- `WORKING_DIR`:/path/to/your_logs
- `CTRL_FILE`:/path/to/sys.ctrl
- `NUM_POTENTIALS`:number of potentials you want to loop over
- `CONSTANT_FILES`:constant file for each 2/3 body potential listed
- `TWO_BODY_POT_FILES`: list of 2 body files assumed in nQMCC/pots/
- `THREE_BODY_POT_FILES`: list of 3 body files assumed in nQMCC/pots/
- `NUM_BLOCKS BLOCK_SIZE WALKERS_PER_NODE`: Monte Carlo Sample information and topology (0 walkers per node unless Abe says otherwise)
- `OPT_SCALE NUM_OPT_EVALUATIONS`: Fraction of Current Value in Decks to set bounds and number of total samples in a single COBYLA optimization
> ⚠️ `IF SYSTEM_TYPE == bound`:
- `ESEP_START ESEP_STOP ESEP_NUM`: ESEP scale range and number of grid points
- `ETA_FLAT`: flat step size for ETA parameters
- `THREE_BODY_FLAT`: flat step size for 3b correlations (delta, epsilon, theta, upsilon, rscal, uscal)
- `SS_FLAT`: flat step size for spatial symmetry (SS) block parameters
- `ALPHA_DIR` *(optional)*: path to directory containing pre-optimized He-4 alpha-core decks; enables alpha-core seeding for improved convergence. If not found the run continues with the ctrl-file default deck.

Optimization runs four phases per ESEP scale, each building on the previous:
1. **Phase 0** — main correlations only (ESEP + global corr) at 50% scale
2. **Phase 1** — first SS block only at 100% scale
3. **Phase 2** — remaining SS blocks + main at 50% scale *(skipped if only one SS block)*
4. **Phase 3** — full optimization (main + all SS) at 100% scale

The best result across all ESEP scales is kept as the final deck. If no scale improves the energy, the original deck is preserved.

Outputs include optimized nQMCC decks (`.dk`), per-pair results in `{pair}_results.json`, a combined summary in `combined_results.json`, and a full run log in `optimization_log.out`; all saved to the working directory.

> ⚠️ `IF SYSTEM_TYPE == sc_scattering`:
- `NUM_CHANNELS`: total number of channels to loop over
- `SCATTERING_CTRL_FILES`: list of /path/to/scattering.ctrl
- `SPATIAL_SYMMETRIES_INDEXS`: index of spatial symmetry
- `OPTIMIZE_TARGET`: Optimize the A-1 body (0: no  1: yes)
- `ENERGY_LOWER_BOUND ENERGY_UPPER_BOUND DELTA_ENERGY`: scan will terminate if out of bounds, delta_e is approximated by foward difference scheme
- `INITIAL_BSCAT INITIAL_DELTA_BSCAT MAX_BSCAT_SLOPE`: first bscat +- delta, used to calculate db/de
- `MAX_NUMBER_OF_DECKS_GENERATED`: scan will terminate if number of evals is greater or equal to

Outputs include optimized nQMCC decks (`.dk`), logs of nQMCC calculations, and a summary of results with `BSCAT`, `EREL`, `VREL`, `WSE`, and others; saved in json format to the working directory.
