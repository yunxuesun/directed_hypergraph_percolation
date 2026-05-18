# Directed Hypergraph Percolation

This repository contains the source code and analytical tools for the paper: 
**"Directionality and node heterogeneity reshape criticality in hypergraph percolation"**. 

The codebase implements a message-passing framework and Monte Carlo simulations to evaluate the emergence of the Hypergraph Giant Weakly Connected Component (HGWCC), Hypergraph Giant In-Component (HGIN), Giant Out-Component (HGOUT), and Giant Strongly Connected Component (HGSCC) in directed hypergraphs with anchor nodes.

## Repository Structure

To ensure computational efficiency and ease of use, the core simulations are written in modularized **C**, while data processing, pipeline automation, and visualization are handled in **Python** and **MATLAB**.

```text
.
├── src/
│   ├── globals.h        # Global variable declarations, macros, and configuration parameters
│   ├── globals.c        # Global variable definitions and memory allocations
│   ├── hypergraph.c     # Hypergraph topology generation, file I/O, and reverse mapping
│   ├── percolation.c    # Core algorithms: Anchor assignment, Message-Passing, and Monte Carlo DFS
│   └── main.c           # Main execution flow and memory cleanup
├── data_processing/     # Directory for empirical network data processing
│   ├── iJO1366.mat               # Raw data for the E. coli metabolic network
│   └── hypergraph_construction.m # MATLAB script to parse stoichiometric matrix
├── Makefile             # Automated compilation script
├── result/              # Directory for output data (created automatically during runtime)
├── fig/                 # Directory for generated plots (created automatically during runtime)
├── README.md            # This documentation
└── run.ipynb            # End-to-end Python pipeline (Compile -> Run -> Plot)
```

## Prerequisites

To compile and run the code, you will need the following installed on your system:

1. **C Compiler:** GCC (GNU Compiler Collection) or any standard C99 compatible compiler.
2. **Python 3.x:** With the following scientific packages installed:
   ```bash
   pip install jupyter pandas matplotlib numpy
   ```
3. **MATLAB (Optional):** Only required if you intend to process empirical network data (e.g., `.mat` files).

## Installation and Compilation

The repository includes a `Makefile` for automated, one-click compilation. Open your terminal, navigate to the repository root, and run:

```bash
# Clean up any previous builds and compile the C code
make clean
make
```

If successful, this will generate an executable file named `percolation_sim` in the root directory.

## Usage

This simulator supports two independent workflows: **Synthetic Network Modeling** and **Real Empirical Network Analysis**. You can easily switch between them by toggling the `USE_REAL_DATA` macro in `src/globals.h`.

### Workflow A: Synthetic Network Modeling

In this mode, the program generates a synthetic directed hypergraph from scratch based on your predefined structural parameters.

1. **Configuration:** Open `src/globals.h` and ensure the real data mode is turned off:
   ```c
   #define USE_REAL_DATA 0
   ```
   *(You can also adjust the synthetic graph scale `N_SYNTHETIC`, `Q_SYNTHETIC` and average cardinalities `M_IN_SYNTHETIC`, `M_OUT_SYNTHETIC` here).*
2. **Recompile:** Run `make` in your terminal.
3. **Run & Analyze:**
   * **(Recommended)** Open `run.ipynb` and execute the cells sequentially. It will automatically run the C program, stream the progress, and plot the theoretical vs. numerical percolation curves into the `fig/` folder.
   * **(CLI)** Alternatively, run `./percolation_sim` in your terminal. By default, the simulation saves the resulting data to `result/results.txt`. You can easily change this destination by modifying the `output_filename` variable located at the top of `src/globals.c`:

### Workflow B: Real Empirical Network Analysis

In this mode, the program loads the genome-scale metabolic network of *E. coli* (iJO1366). Metabolites act as nodes, and biochemical reactions act as directed hyperedges.

1. **Data Preprocessing (MATLAB):**
   * Ensure `iJO1366.mat` and `hypergraph_construction.m` are located in the `data_processing/` directory.
   * Run the `hypergraph_construction.m` script in MATLAB.
   * The script parses the stoichiometric matrix, filters out self-loops, separates reactants/products, and generates the topology file `iJO1366_hypergraph.txt` required by the C code.
2. **Configuration:** Open `src/globals.h` and turn on the real data mode:
   ```c
   #define USE_REAL_DATA 1
   ```
   *(Note: The empirical $N$ and $Q$ will be read directly from the parsed text file, overriding any synthetic size macros).*
3. **Recompile:** Run `make` in your terminal.
4. **Run & Analyze:**
   * Use `run.ipynb` for the automated end-to-end pipeline, or execute `./percolation_sim` directly in the terminal.

## Advanced Configuration

You can configure various physical parameters and percolation sweep settings by modifying the macros defined at the top of `src/globals.h`:

* `ANCHOR_SELECTION_MODE`: Defines the candidate pool for anchor nodes. Options include `ANCHOR_MODE_INPUT` (1), `ANCHOR_MODE_OUTPUT` (2), or `ANCHOR_MODE_GLOBAL` (3).
* `THETA`: The probability of a candidate node acting as an anchor.
* `SWEEP_MODE`: Determines how percolation parameters are varied during the simulation:
    * `1`: Sweep node retention probability (p_N), keep hyperedge probability fixed (`FIXED_P_H`).
    * `2`: Sweep hyperedge retention probability (p_H), keep node probability fixed (`FIXED_P_N`).
    * `3`: Sweep both parameters simultaneously (p_N = p_H).
* `MONTE_CARLO_RUNS`: The number of independent iterations for numerical averaging to eliminate finite-size random fluctuations.
