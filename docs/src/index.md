```@meta
CurrentModule = SmoQyHankelCorrCleaner
```

# SmoQyHankelCorrCleaner.jl

[SmoQyHankelCorrCleaner.jl](https://github.com/SmoQySuite/SmoQyHankelCorrCleaner.jl) the method introduced in the paper
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