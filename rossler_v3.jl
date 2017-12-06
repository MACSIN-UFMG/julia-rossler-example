using Plots
using DifferentialEquations
plotlyjs()

# Rossler space equation
function rossler(t::Float64, x::Vector{Float64}, dx::Vector{Float64})
    dx[1] = -x[2] - x[3]
    dx[2] = x[1] + 0.2*x[2]
    dx[3] = 0.2 + x[3]*(x[1]-5.7)
    return dx
end

# Simulate Rossler
time_span = (0.0, 1000.0)
x0 = [1.0, 1.0, 2.0]
prob = ODEProblem(rossler, x0, time_span)
sol = solve(prob, RK4(), reltol=1e-8)

plot(sol, vars=(1, 2, 3))
