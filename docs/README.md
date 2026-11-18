# Reference Documents for the Pepper Cross-Section Demo

This folder indexes the works cited in the demo source, and holds any assets
the README files use. The index points at each work's canonical source rather than
redistributing it, so every entry below links to a DOI or an arXiv record.

The top-level [README](../README.md) maps each reference to the file that
uses it.

## Parton distribution functions

| Citation | Source |
|---|---|
| Hou, T.-J. et al. (2021). "New CTEQ global analysis of quantum chromodynamics with high-precision data from the LHC." *Physical Review D*, 103, 014013. | [arXiv:1912.10053](https://arxiv.org/abs/1912.10053) |
| Buckley, A. et al. (2015). "LHAPDF6: parton density access in the LHC precision era." *European Physical Journal C*, 75, 132. | [arXiv:1412.7420](https://arxiv.org/abs/1412.7420) |

## Reference by source file

| File | References |
|---|---|
| `pepper-lhapdfLite.c` | Hou et al. (2021), Buckley et al. (2015) |

## Upstream tutorial

The demo is a C++ port of the parton-level example from the MCnet Schools
Stirling 2025 tutorial, tracked as the `submodules/pepper_python` submodule and
available at
[gitlab.com/hepcedar/mcnet-schools/stirling-2025](https://gitlab.com/hepcedar/mcnet-schools/stirling-2025).

