# SmoQyHankelCorrCleaner.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://SmoQySuite.github.io/SmoQyHankelCorrCleaner.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://SmoQySuite.github.io/SmoQyHankelCorrCleaner.jl/dev/)
[![Build Status](https://github.com/SmoQySuite/SmoQyHankelCorrCleaner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SmoQySuite/SmoQyHankelCorrCleaner.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/SmoQySuite/SmoQyHankelCorrCleaner.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/SmoQySuite/SmoQyHankelCorrCleaner.jl)

[SmoQyHankelCorrCleaner.jl](https://github.com/SmoQySuite/SmoQyHankelCorrCleaner.jl) implements the method introduced in the paper
["Denoising of Imaginary Time Response Functions with Hankel Projections"](https://arxiv.org/abs/2403.12349)
for denoising imaginary time correlation data, the citation for which is given below:

```bibtex
@article{Yu2024Denoising,
  title={Denoising of Imaginary Time Response Functions with Hankel projections},
  author={Yu, Yang and Kemper, Alexander F and Yang, Chao and Gull, Emanuel},
  journal={arXiv preprint arXiv:2403.12349},
  year={2024}
}
```

## Installation

To install [`SmoQyHankelCorrCleaner.jl`](https://github.com/SmoQySuite/SmoQyHankelCorrCleaner.jl.git),
simply open the Julia REPL and run the commands
```julia
julia> ]
pkg> add SmoQyHankelCorrCleaner
```
or equivalently via `Pkg` do
```julia
julia> using Pkg; Pkg.add("SmoQyHankelCorrCleaner")
```

## Documentation

- [STABLE](https://SmoQySuite.github.io/SmoQyHankelCorrCleaner.jl/stable/): Documentation for the latest version of the code published to the Julia [`General`](https://github.com/JuliaRegistries/General.git) registry.
- [DEV](https://SmoQySuite.github.io/SmoQyHankelCorrCleaner.jl/dev/): Documentation for the latest commit to the `main` branch.