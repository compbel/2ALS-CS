# 2ALS-CS

`2ALS-CS` infers rooted binary phylogenies from binary character matrices by
exact 2-adic local search for Camin–Sokal parsimony.

The program is a heuristic local-search method: every active neighborhood is
enumerated completely and every accepted move is a strict improvement, but a
reported local optimum is not a certificate of the global parsimony optimum.
There are no restarts, plateau moves, ratchets, annealing steps, or hybrid
tree-search moves.

## Requirements

- MATLAB R2024b or newer
- Symbolic Math Toolbox
- Parallel Computing Toolbox when `ParallelFamilies=true`

The serial mode does not require Parallel Computing Toolbox.

## Installation

Clone the repository and add its `matlab` directory to the MATLAB path:

```bash
git clone https://github.com/compbel/2ALS-CS.git
cd 2ALS-CS
```

```matlab
addpath(fullfile(pwd, 'matlab'));
```

No compilation or package installation is required.

## Input

The input is an `n`-by-`m` numeric or logical matrix. Rows are taxa, columns
are binary irreversible characters, and every entry must be exactly `0` or
`1`. Duplicate columns are supported and retain their objective
multiplicities.

```matlab
M = logical([ ...
    1 0 0 1; ...
    1 1 0 1; ...
    0 1 1 0; ...
    0 0 1 0]);
```

Taxon `i` is represented by row `i` and by leaf label `i` in returned trees.

## Quick start

```matlab
addpath('matlab');
rng(1, 'twister');                 % reproducible generated initial tree

result = twoALSCS(M, ...
    'PairMode', 'ALWAYS', ...
    'ParallelFamilies', true, ...
    'ParallelWorkers', 8, ...
    'ParallelPoolType', 'processes');

result.initialScore
result.finalScore
result.finalTree
```

For a machine without Parallel Computing Toolbox:

```matlab
result = twoALSCS(M, 'ParallelFamilies', false);
```

The bundled example runs in serial mode:

```matlab
run('examples/run_example.m');
```

## Main options

| Option | Default | Meaning |
|---|---:|---|
| `InitialTree` | empty `digraph` | User-supplied rooted binary tree. |
| `InitialTreeModel` | `YULE` | Generated start: `YULE` or `UNIFORM`. |
| `PairMode` | `ALWAYS` | Pair-family schedule: `ALWAYS`, `LOCAL`, or `NEVER`. |
| `MaxIterations` | `Inf` | Maximum number of complete neighborhood scans. |
| `ParallelFamilies` | `true` | Process transformation families in parallel. |
| `ParallelWorkers` | `4` | Exact requested worker count. |
| `ParallelPoolType` | `processes` | `processes` or `threads`. |
| `ValidateNeighborhood` | `false` | Enable expensive internal consistency checks. |

`PairMode` controls the search neighborhood:

- `ALWAYS`: scan single-column and allelic-state pair families jointly at
  every iteration. This is the accuracy-first default.
- `LOCAL`: scan single-column families first; scan pair families only when
  the single-column search is locally optimal. After a successful pair move,
  return to single-column families.
- `NEVER`: use only single-column families.

## Supplying an initial tree

Trees are MATLAB `digraph` objects. Edges point from parent to child. The
node table must contain a `Taxon` variable: leaves have labels `1,...,n` and
internal vertices have label `0`.

```matlab
result = twoALSCS(M, ...
    'InitialTree', myTree, ...
    'PairMode', 'ALWAYS', ...
    'ParallelFamilies', false);
```

The tree is validated before the search starts. If `InitialTree` is empty, a
random tree is generated independently of the alignment.

## Result structure

Important fields are:

| Field | Meaning |
|---|---|
| `initialTree`, `finalTree` | Rooted binary `digraph` objects. |
| `initialScore`, `finalScore` | Camin–Sokal parsimony scores. |
| `converged` | True only after a complete nonimproving scan. |
| `terminationReason` | `local_optimum` or `maximum_iterations`. |
| `trajectory` | Complete iteration summaries and selected trees. |
| `acceptedMoves` | Replayable strict-improvement witnesses. |
| `pairMode` | Effective pair-family schedule. |
| `singleColumnFamilyCount`, `pairFamilyCount` | Active family counts. |

Each accepted move records its source and target scores and trees, family
identifier, exact representative `theta`, evaluated `gValues` and `hValues`,
and transformed coordinates. A move can be replayed as

```matlab
coordinates = move.gValues - move.theta .* move.hValues;
```

## Algorithm outline

1. Canonically embed the current rooted binary tree into distinct integers.
2. Construct one transformation for each unique alignment column and the
   applicable transformations for each unordered column pair.
3. Compute all exact collision parameters shared by these transformations.
4. Enumerate the complete 2-adic parameter-region neighborhood of every
   active family.
5. Score candidate embeddings directly from pairwise 2-adic valuations,
   without constructing a tree for each candidate.
6. Accept the globally best strict improvement, reconstruct that one tree,
   and repeat.

Parallel workers receive a read-only collision table. Each worker returns
only its best family candidate and replay witness. Completion order does not
affect deterministic tie breaking.

## Reproducibility

When the initial tree is generated, set MATLAB's random generator before the
call:

```matlab
rng(20261110, 'twister');
result = twoALSCS(M, 'ParallelFamilies', false);
```

For exact replication across machines, save and reuse `result.initialTree`.

## Tests

Run the release tests from the repository root:

```matlab
addpath('tests');
report = run_release_tests();
```

The tests check the public API, input rejection, user-supplied initial trees,
strict score descent, deterministic serial replay, duplicate columns,
allelic-state pair families, final-tree rescoring, and accepted-move
witnesses.

## Limitations

- The method returns a strict local optimum, not a global-optimality proof.
- Complete pair-family exploration can be expensive as the number of unique
  columns grows.
- The exact integer fast path enforces conservative magnitude bounds; an
  unusually deep tree can exceed them.
- Performance depends on the number of pair families and parameter regions,
  not only on the numbers of taxa and characters.
