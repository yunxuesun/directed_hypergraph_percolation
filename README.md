# Directed Hypergraph Percolation

This repository contains the source code and analytical tools for the paper: 
**"Directionality and node heterogeneity reshape criticality in hypergraph percolation"**. 

The codebase implements a message-passing framework and Monte Carlo simulations to evaluate the emergence of the Hypergraph Giant In-Component (HGIN), Giant Out-Component (HGOUT), and Giant Strongly Connected Component (HGSCC) in directed hypergraphs with anchor nodes.

## Repository Structure

To ensure computational efficiency and ease of use, the core simulations are written in modularized **C**, while data processing, pipeline automation, and visualization are handled in **Python**.

```text
.
├── src/
│   ├── globals.h                # Global variable declarations, macros, and function prototypes
│   ├── globals.c                # Global variable definitions and memory allocations
│   ├── hypergraph.c             # Hypergraph topology generation, file I/O, and reverse mapping
│   ├── percolation.c            # Core algorithms: Anchor assignment, Message-Passing, and Monte Carlo SCC/DFS
│   └── main.c                   # Main execution flow and memory cleanup
├── Makefile                     # Automated compilation script
├── result/                      # Directory for output data (created automatically during runtime)
├── README.md                    # This documentation
└── Analysis_and_Plotting.ipynb  # End-to-end Python pipeline (Compile -> Run -> Plot)
```

## Prerequisites

To compile and run the code, you will need the following installed on your system:

1. **C Compiler:** GCC (GNU Compiler Collection) or any standard C99 compatible compiler.
2. **Python 3.x:** With the following packages installed:
   ```bash
   pip install jupyter pandas matplotlib numpy
   ```

## Installation and Compilation

The repository includes a `Makefile` for automated, one-click compilation. Open your terminal, navigate to the repository root, and run:

```bash
# Clean up any previous builds and compile the C code
make clean
make
```

If successful, this will generate an executable file named `percolation_sim` in the root directory.

## Usage

### Option 1: End-to-End Analysis Pipeline (Recommended)

For the easiest experience and to perfectly reproduce the analysis workflow shown in the paper, we provide an interactive Jupyter Notebook. 

Open **`Analysis_and_Plotting.ipynb`**. This notebook acts as an automated pipeline:
1. It automatically triggers the `Makefile` to compile the C source code.
2. It executes the simulation (`./percolation_sim`) and streams the progress.
3. It loads the generated simulation data from `result/ph_data.txt`.
4. It plots the theoretical (Message-Passing) vs. numerical (Monte Carlo) percolation curves side-by-side.

### Option 2: Manual Command-Line Execution

If you prefer to run the simulation manually or integrate it into bash scripts:

1. Ensure the code is compiled using `make`.
2. Create the result directory (if it doesn't exist):
   ```bash
   mkdir -p result
   ```
3. Run the executable:
   ```bash
   ./percolation_sim
   ```
4. The simulation progress will be printed to the terminal, and the final data will be saved to `result/ph_data.txt`. You can then process this CSV-formatted text file with your preferred tools.

## Configuration

You can adjust the simulation parameters (e.g., synthetic network size, anchor probabilities, Monte Carlo runs) by modifying the macros defined at the top of `src/globals.h`. After making any changes, simply re-run the Jupyter Notebook or use `make` to recompile.

- `USE_REAL_DATA`: Toggle between synthetic hypergraph generation (`0`) and loading empirical networks (`1`).
- `SET_ALL_INPUTS_AS_ANCHORS`: Toggle the biologically-motivated asymmetric anchor mode (`1`) vs. purely random anchor assignment (`0`).
- `THETA`: The probability of a node acting as an anchor (used when random assignment is active).
