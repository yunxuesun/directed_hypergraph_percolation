# Directed Hypergraph Percolation

This repository contains the source code and analytical tools for the paper:
**"Directionality and node heterogeneity reshape criticality in hypergraph percolation"**.

The codebase implements a message-passing framework and Monte Carlo simulations to evaluate the emergence of the Hypergraph Giant Weakly Connected Component (HGWCC), Hypergraph Giant In-Component (HGIN), Giant Out-Component (HGOUT), and Giant Strongly Connected Component (HGSCC) in directed hypergraphs with anchor nodes.

## Repository Structure

To ensure computational efficiency and ease of use, the core simulations are written in modularized **C**, while data processing, pipeline automation, and visualization are handled in **Python** and **MATLAB**.

```text
.
|-- src/
|   |-- globals.h        # Global variable declarations, macros, and configuration parameters
|   |-- globals.c        # Global variable definitions and memory allocations
|   |-- hypergraph.c     # Hypergraph topology generation, file I/O, and reverse mapping
|   |-- percolation.c    # Core algorithms: Anchor assignment, Message-Passing, and Monte Carlo DFS
|   `-- main.c           # Main execution flow and memory cleanup
|-- data/                # Directory for empirical network data processing
|   |-- iJO1366.mat               # Raw data for the E. coli metabolic network
|   `-- hypergraph_construction.m # MATLAB script to parse stoichiometric matrix
|-- seq_generator.c      # Auxiliary C generator for degree/cardinality sequences
|-- Makefile             # Automated compilation script
|-- result/              # Directory for output data (created automatically during runtime)
|-- fig/                 # Directory for generated plots (created automatically during runtime)
|-- README.md            # This documentation
`-- run.ipynb            # End-to-end Python pipeline and additional analysis plots
```

## Prerequisites

To compile and run the code, you will need the following installed on your system:

1. **C Compiler:** GCC (GNU Compiler Collection) or any standard C99 compatible compiler.
2. **Python 3.x:** With the following scientific packages installed:
   ```bash
   pip install jupyter pandas matplotlib numpy scipy
   ```
3. **MATLAB (Optional):** Only required if you intend to process empirical network data (e.g., `.mat` files).

## Installation and Compilation

The repository includes a `Makefile` for automated compilation. Open your terminal, navigate to the repository root, and run:

```bash
make clean
make
```

If successful, this will generate `percolation_sim` and the auxiliary `seq_generator` executable in the root directory. The notebook also falls back to direct `gcc` compilation when `make` is not available.

## Usage

This simulator supports two independent workflows: **Synthetic Network Modeling** and **Real Empirical Network Analysis**. You can switch between them by toggling the `USE_REAL_DATA` macro in `src/globals.h`.

### Workflow A: Synthetic Network Modeling

In this mode, the program generates a synthetic directed hypergraph from scratch based on predefined structural parameters.

1. **Configuration:** Open `src/globals.h` and ensure the real data mode is turned off:
   ```c
   #define USE_REAL_DATA 0
   ```
   You can also adjust the synthetic graph scale `N_SYNTHETIC`, `Q_SYNTHETIC`, and average cardinalities `M_IN_SYNTHETIC`, `M_OUT_SYNTHETIC` here.
2. **Recompile:** Run `make` or use the compile cell in `run.ipynb`.
3. **Run and analyze:**
   * **Recommended:** Open `run.ipynb` and execute the cells sequentially. It compiles the C programs, runs the simulation, streams progress, and plots theoretical versus numerical percolation curves into the `fig/` folder.
   * **CLI:** Alternatively, run `./percolation_sim` in your terminal. By default, the simulation saves data to `result/ph_data.txt`. You can change this by modifying `output_filename` in `src/globals.c`.

### Workflow B: Real Empirical Network Analysis

In this mode, the program loads the genome-scale metabolic network of *E. coli* (iJO1366). Metabolites act as nodes, and biochemical reactions act as directed hyperedges.

1. **Data preprocessing (MATLAB):**
   * Ensure `iJO1366.mat` and `hypergraph_construction.m` are located in the `data/` directory.
   * Run the `hypergraph_construction.m` script in MATLAB.
   * The script parses the stoichiometric matrix, filters self-loops, separates reactants/products, and generates the topology file required by the C code.
2. **Configuration:** Open `src/globals.h` and turn on real data mode:
   ```c
   #define USE_REAL_DATA 1
   ```
   The empirical `N` and `Q` values are read directly from the parsed text file, overriding the synthetic size macros.
3. **Recompile:** Run `make` or use the compile cell in `run.ipynb`.
4. **Run and analyze:** Use `run.ipynb` for the automated end-to-end pipeline, or execute `./percolation_sim` directly in the terminal.

### Additional Notebook Analysis

The final two notebook cells provide optional figure-generation workflows:

* `Figure_giant_component.pdf`: compares Poisson and scale-free topologies under correlated and uncorrelated settings when the four expected result files are available in `result/`.
* `Figure_delta_pc.pdf`: calls `seq_generator` to sample realized degree/cardinality sequences and estimate the directed-versus-undirected critical threshold gap across anchor fractions. The notebook uses `N_SEQUENCE_SAMPLE = 2000` for interactive runtime; increase this value for larger samples.

`seq_generator` can also be called directly:

```bash
./seq_generator <topology_id> <lambda> <gamma> [N]
```

where `topology_id` is `0` for uncorrelated Poisson, `1` for correlated Poisson, `2` for uncorrelated scale-free, and `3` for correlated scale-free. The optional `N` argument defaults to `10000`.

## Advanced Configuration

You can configure physical parameters and percolation sweep settings by modifying the macros defined at the top of `src/globals.h`:

* `ANCHOR_SELECTION_MODE`: Defines the candidate pool for anchor nodes. Options include `ANCHOR_MODE_INPUT` (1), `ANCHOR_MODE_OUTPUT` (2), or `ANCHOR_MODE_GLOBAL` (3).
* `THETA`: The probability of a candidate node acting as an anchor.
* `SWEEP_MODE`: Determines how percolation parameters are varied during the simulation:
  * `1`: Sweep node retention probability (`p_N`), keep hyperedge probability fixed (`FIXED_P_H`).
  * `2`: Sweep hyperedge retention probability (`p_H`), keep node probability fixed (`FIXED_P_N`).
  * `3`: Sweep both parameters simultaneously (`p_N = p_H`).
* `MONTE_CARLO_RUNS`: The number of independent iterations for numerical averaging.
