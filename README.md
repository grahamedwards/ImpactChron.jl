# ImpactChron.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://grahamedwards.github.io/ImpactChron.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://grahamedwards.github.io/ImpactChron.jl/dev/)
[![Build Status](https://github.com/grahamedwards/ImpactChron.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/grahamedwards/ImpactChron.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/grahamedwards/ImpactChron.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/grahamedwards/ImpactChron.jl)

A Bayesian framework to deconvolve the dual contributions of (1.) cooling in primary planetesimals and (2.) impact heating to thermochronologic ages in chondritic meteorites.

Developed from the statistical framework of [Chron.jl](https://github.com/brenhinkeller/Chron.jl) ([doi:10.17605/osf.io/TQX3F](https://doi.org/10.17605/osf.io/TQX3F)).

📄 Associated pre-print coming soon! 📄

---
---
## Installation
`ImpactChron.jl` is written in the Julia programming language. To install, open an instance of Julia, enter the package manager (type `]` in the REPL), and type:
```julia
pkg> add https://github.com/grahamedwards/ImpactChron.jl
```

## Usage
The primary function in this package is `thermochron_metropolis`. This executes a Markov chain Monte Carlo algorithm (an adaptation of the Metropolis algorithm) that uses a hierarchical set of priors (see DATA.md for additional information) to constrain the thermal histories of chondrite meteorite parent bodies. 

Support functions like `planetesimal_temperature`, `planetesimal_cooling_dates`, and `asteroid_agedist!` may be used to explore specific asteroid thermal histories. 

The package documentation outlines how to implement this and other functions. To access documentation from Julia's REPL, type `?`, followed by the function name, e.g.
```julia
help?> thermochron_metropolis
```

## Example
An example notebook to show how to use the major features of the package will be available soon! 