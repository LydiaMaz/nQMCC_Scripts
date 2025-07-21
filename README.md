# nQMCC Parameter Optimization Framework

This project automates variational Monte Carlo (VMC) optimizations for He4+n scattering using the nQMCC framework.

## Usage

To run a full optimization:

```bash
python /absolute/path/to/AutoOpt.py --utility /absolute/path/to/file.util --system /absolute/path/to/file.sys --run_cmd "<executable command>"
```
> ⚠️ Make sure your `.util` file includes the following fields:

- `working_dir`: Absolute path where output files will be saved  
- `input_dir`: Directory containing nuclear input files (e.g., decks, params)  
- `run_cmd`: Executable command for running nQMCC (depends on your system, e.g. srun -n 256 -c 2 --cpu_bind=cores)

## Modules

### `util.py`

Parses `.util` files that specify global run configuration for QMC automation.

- Uses the `Utility` class to:
  - Read `.util` files and expose fields like `run_cmd`, `working_dir`, etc.
  - Copy and update control files for target and scattering systems.
- Integrates with the `Control` class to adjust file paths.

---

### `nuc_system.py`

Parses and stores nuclear system configuration from `.system` input files.

- Uses the `NuclearSystem` class to:
  - Load system-wide parameters such as interaction potentials, energy bounds, and symmetry channels.
  - Access and modify parameters via dictionary-like interface.
  - Write updated system configuration back to disk.

---

### `control.py`

Parses `.ctrl` input files used to control QMC simulations.

- Uses the `Control` class to:
  - Read/write simulation control files.
  - Store all simulation parameters in a Python dictionary for easy lookup.
  - Write updated control files with new paths, optimization inputs, etc.
- Internally uses:
  - `WavefunctionInputs` – a helper class to manage wavefunction-specific parameters like deck/param paths.

---

### `wavefunction_input.py`

Parses wavefunction blocks present in QMC deck files. Supports different wavefunction types via subclasses.

- **`WavefunctionInput`** *(abstract base class)*  
  Base framework for handling groups of wavefunction parameters, including reading, setting, and writing parameter blocks. Params stored as Python attributes.

- **`VariationalWavefunctionInput`**  
  Parses standard variational wavefunction input sections, including param and deck files plus parameter groups.

- **`ProductWavefunctionInput`**  
  Supports product-type wavefunctions with multiple param and deck files and parameter groups.

- **`DeuteronWavefunctionInput`**  
  Minimal subclass for deuteron wavefunction input where no parameter blocks are required.

- Used internally by `control.py`.

---

### `parameters.py`

Parses `.params` files which define structural layout for decks.

- Provides `Parameters` class used throughout `deck.py` and `zero_params.py`.
- Params stored as Python attributes

---

### `deck.py`

Parses and writes QMC `.dk` (deck) files containing variational parameters.

- Uses the `Deck` class to:
  - Read scalar, array, and symmetry-based parameters (stored in Python dictionary).
  - Track grouped parameter lines for correct output formatting.
  - Automatically detect and parse spatial symmetries.

- Internally uses:
  - `Parameters` (from `parameters.py`) for deck layout info.
  - `read_spatial_symmetries()` from `spatial_symmetry.py`.

- Provides:
  - `read_params_and_deck()` – loads both `.params` and `.dk` files into a `Deck`.

---

### `spatial_symmetry.py`

Parses spatial symmetry blocks present in QMC deck files.

- Spatial symmetry related parameters stored as Python attributes.
- Provides `read_spatial_symmetries()` to extract group-structured float parameters.
- Used internally by `deck.py`.

---

### `zero_params.py`

Provides tools to set variational parameters before optimization.

- `zero_var_params()` sets specific correlation groups to specific value/scale factor as instructed by mode
- `save_opt_file()` saves modified decks (.opt) for use in optimization.

---

### `bscat_optimizer.py`

Performs a forward and backward scan over `bscat`:

- **Step 1:** Computes initial energy and slope (`db/de`) near a given `bscat` using central difference method.
- **Step 2:** Walks forward in `bscat` to target higher relative energies.
- **Step 3:** Walks backward to reach lower relative energies.
- **Step 4:** Collects results (energy, variance, WSE, deck path) at each point.
- **Step 5:** Outputs sorted results to `optimized_decks.json` and optionally saves `.opt` deck files.

Intermediate `.dk` files are cleaned from the `scratch/` directory after the run.

---

### `AutoOpt.py`

Automates optimization of a scattering system (He⁴+n):

- **Step 0:** Initializes system; copies necessary files to working directory.   
- **Step 1:** Optimizes the bound target (He⁴).
- **Step 2:** Copies target parameters into the scattering deck.
- **Step 3:** Searches `bscat` values to find initial relative energy near 3 MeV.
- **Step 4:** Runs `bscat_optimizer` to dynamically search range of `bscat` values.
- **Step 5:** Saves optimized decks and results to JSON file.

Outputs include optimized `.dk` files, logs of run outputs, and a json summary of results with `bscat`, `E_rel`, `wse`, and variance. All files saved to working directory.