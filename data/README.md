# iJO1366 directed-hypergraph construction

`hypergraph_construction.py` is an open-source Python replacement for
`hypergraph_construction.m`.  It reads the stoichiometric matrix `Model.S`
from `iJO1366.mat` and writes the directed-hypergraph format used by the
accompanying C implementation.

## Requirements

Python 3.9 or later, plus the open-source packages listed in
`requirements.txt`:

```bash
python -m pip install -r requirements.txt
```

## Reproduce the supplied hypergraph

From this directory, run:

```bash
python hypergraph_construction.py
```

The command writes `iJO1366_hypergraph.txt`.  The output contains 1805 nodes
and 2253 directed hyperedges and reproduces the supplied file.  To keep the
existing file unchanged while testing, choose a different output path:

```bash
python hypergraph_construction.py --output reproduced_iJO1366_hypergraph.txt
```

For each reaction, metabolites with negative stoichiometric coefficients are
input nodes and those with positive coefficients are output nodes.  Reactions
without at least one input and one output are excluded.  Node identifiers are
zero-based in the exported text file.
