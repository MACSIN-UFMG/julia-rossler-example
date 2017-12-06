# julia-rossler-example
Rossler simulation and bifurcation diagram in Julia.

## Requirements

This example requires the packages "PlotlyJS", "PyPlot", "Plots", "DifferentialEquations", "LaTeXStrings" and "ProgressMeter". The instalation of those can be done from julia REPL:
```julia
julia> Pkg.add("PlotlyJS")
julia> Pkg.add("PyPlot")
julia> Pkg.add("Plots")
julia> Pkg.add("DifferentialEquations")
julia> Pkg.add("LaTeXStrings")
julia> Pkg.add("ProgressMeter")
```

## Files Desciption

- ``rossler_v1.jl``: Simulation of Rossler using inefficient implementation of 4-th order Runge-Kutta.
- ``rossler_v2.jl``: Simulation of Rossler using more efficient implementation of 4-th order Runge-Kutta. It uses inplace operations and the BLAS library for matricial operations.
- ``rossler_v3.jl``: Simulation of Rossler using "DifferentialEquations" library.
- ``rossler_v4.jl``: Example of how to plot bifurcation diagram for Rossler atractor.

