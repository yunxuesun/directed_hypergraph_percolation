#!/usr/bin/env python3
"""Construct the directed iJO1366 hypergraph from the supplied MATLAB model.

This is an open-source Python replacement for ``hypergraph_construction.m``.
For each reaction (a column of ``Model.S``), metabolites with negative
stoichiometric coefficients are treated as input nodes and metabolites with
positive coefficients as output nodes.  Reactions without both types of node,
or with a node on both sides, are excluded.  Node identifiers in the output
are zero-based, matching the format used by the accompanying C code.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterator

import numpy as np
from scipy import io, sparse


SCRIPT_DIRECTORY = Path(__file__).resolve().parent


def load_stoichiometric_matrix(input_path: Path):
    """Load ``Model.S`` from a MATLAB v5 MAT file."""
    contents = io.loadmat(input_path, struct_as_record=False, squeeze_me=True)
    if "Model" not in contents:
        raise KeyError("The input MAT file does not contain a 'Model' variable.")

    model = contents["Model"]
    if not hasattr(model, "S"):
        raise KeyError("The 'Model' variable does not contain an 'S' matrix.")
    return model.S


def directed_hyperedges(stoichiometric_matrix) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield zero-based input/output node arrays for each valid reaction."""
    for reaction_index in range(stoichiometric_matrix.shape[1]):
        if sparse.issparse(stoichiometric_matrix):
            reaction = stoichiometric_matrix.getcol(reaction_index).toarray().ravel()
        else:
            reaction = np.asarray(stoichiometric_matrix[:, reaction_index]).ravel()

        input_nodes = np.flatnonzero(reaction < 0)
        output_nodes = np.flatnonzero(reaction > 0)
        has_self_loop = np.intersect1d(input_nodes, output_nodes).size > 0

        if input_nodes.size and output_nodes.size and not has_self_loop:
            yield input_nodes, output_nodes


def write_hypergraph(output_path: Path, num_nodes: int, hyperedges: list[tuple[np.ndarray, np.ndarray]]) -> None:
    """Write the text format consumed by the accompanying C implementation."""
    with output_path.open("w", encoding="utf-8", newline="\n") as output_file:
        output_file.write(f"{num_nodes} {len(hyperedges)}\n")
        for input_nodes, output_nodes in hyperedges:
            output_file.write(f"{input_nodes.size} {output_nodes.size}\n")
            output_file.write(" ".join(map(str, input_nodes)) + "\n")
            output_file.write(" ".join(map(str, output_nodes)) + "\n")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=SCRIPT_DIRECTORY / "iJO1366.mat",
        help="Input MATLAB model file (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=SCRIPT_DIRECTORY / "iJO1366_hypergraph.txt",
        help="Output directed-hypergraph file (default: %(default)s)",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    matrix = load_stoichiometric_matrix(arguments.input)
    hyperedges = list(directed_hyperedges(matrix))

    write_hypergraph(arguments.output, matrix.shape[0], hyperedges)
    print(f"Wrote {len(hyperedges)} directed hyperedges for {matrix.shape[0]} nodes to {arguments.output}")


if __name__ == "__main__":
    main()
