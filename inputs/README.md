# Input Files

## Parton distribution grid

The CT18NLO parton distribution grid, member 0, in the `lhagrid1` format that
`src/pepper-lhapdfLite.c` reads. `CT18NLO.info` carries the set metadata and
the `α_s` knots, and `CT18NLO_0000.dat` carries the grid itself.

The demo resolves this directory relative to its working directory. Running
from the repository root it reads `inputs/CT18NLO/`, and on the Signaloid Cloud
Compute Engine, where the working directory is `inputs/`, it reads `CT18NLO/`.
No environment variable is needed in either case.

The grid is distributed with LHAPDF by the CTEQ-TEA collaboration. See
[`docs/README.md`](../docs/README.md) for the citation.
