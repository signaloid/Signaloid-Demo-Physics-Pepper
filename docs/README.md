# Reference Documents for the Pepper Cross-Section Demo

This folder indexes the works cited in the demo source, and holds any assets
the READMEs use. The index points at each work's canonical source rather than
redistributing it, so every entry below links to a DOI or an arXiv record.

The top-level [README](../README.md) maps each reference to the file that
uses it.

## Matrix elements and phase space

| Citation | Source |
|---|---|
| Berends, F.A. and Giele, W.T. (1988). "Recursive calculations for processes with n gluons." *Nuclear Physics B*, 306, 759-808. | [doi:10.1016/0550-3213(88)90442-7](https://doi.org/10.1016/0550-3213(88)90442-7) |
| Plätzer, S. and Sjödahl, M. (2023). "Sampling the phase space of hadron collider events." | [arXiv:2302.10449](https://arxiv.org/abs/2302.10449) |

## Parton distribution functions

| Citation | Source |
|---|---|
| Hou, T.-J. et al. (2021). "New CTEQ global analysis of quantum chromodynamics with high-precision data from the LHC." *Physical Review D*, 103, 014013. | [arXiv:1912.10053](https://arxiv.org/abs/1912.10053) |
| Buckley, A. et al. (2015). "LHAPDF6: parton density access in the LHC precision era." *European Physical Journal C*, 75, 132. | [arXiv:1412.7420](https://arxiv.org/abs/1412.7420) |

## Reference by source file

| File | References |
|---|---|
| `pepper-berendsGiele.cpp` | Berends and Giele (1988) |
| `pepper-phaseSpace.cpp` | Plätzer and Sjödahl (2023) |
| `pepper-lhapdfLite.c` | Hou et al. (2021), Buckley et al. (2015) |

## Upstream tutorial

The demo is a C++ port of the parton-level example from the MCnet Schools
Stirling 2025 tutorial, tracked as the `submodules/pepper_python` submodule and
available at
[gitlab.com/hepcedar/mcnet-schools/stirling-2025](https://gitlab.com/hepcedar/mcnet-schools/stirling-2025).

