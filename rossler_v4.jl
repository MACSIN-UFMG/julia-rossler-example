using Plots
using DifferentialEquations
using LaTeXStrings
using ProgressMeter
plotlyjs()

# Rossler space equation
g = @ode_def Rossler begin
    dx = -y - z
    dy = x + a*y
    dz = b + z*(x-c)
end a=>0.2 b=>0.2 c=>5.7

# Simulate
time_span = (0.00, 3000.0)
x0 = [1.0, 1.0, 2.0]
prob = ODEProblem(g, x0, time_span)
sol = solve(prob, RK4(), dt=0.1)

plot(sol, vars=(1, 2))

# Plot Poincare section
xconst = 0
plot!(zeros(100), linspace(0.0, 10.0, 100), lw=4, color="black")

# Get Poincare section
init = 600
poincare_y = []
poincare_z = []
p = zeros(3)
for i = init:length(sol)-1
    if sol[i][1] > 0 && sol[i+1][1] < 0
        p = (sol[i+1] - sol[i])/(sol[i+1][1] - sol[i][1])*(xconst - sol[i][1]) + sol[i]
        push!(poincare_y, p[2])
        push!(poincare_z, p[3])
    end
end

scatter(poincare_y, poincare_z)

# Bifurcation diagram
list_y = []
list_params = []
@showprogress 1 "Computing..." for param = 0.01:0.01:2.00
    # Solve Rossler
    h = Rossler(b=param)
    prob = ODEProblem(h, x0, time_span)
    sol = solve(prob, RK4(), dt=0.1)
    # Find variables
    for i = init:length(sol)-1
        if sol[i][1] > 0 && sol[i+1][1] < 0
            p = (sol[i+1] - sol[i])/(sol[i+1][1] - sol[i][1])*(xconst - sol[i][1]) + sol[i]
            push!(list_y, p[2])
            push!(list_params, param)
        end
    end
end

pyplot()
scatter(list_params, list_y, color="black", markersize=0.01, markerstrokewidth=0.5)
