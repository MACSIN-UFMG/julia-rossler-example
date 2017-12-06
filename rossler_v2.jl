using Plots
plotlyjs()

# Rossler space equation
function rossler(t::Float64, x::Vector{Float64}, dx::Vector{Float64},
                 a::Float64, b::Float64, c::Float64)
    dx[1] = -x[2] - x[3]
    dx[2] = x[1] + a*x[2]
    dx[3] = b + x[3]*(x[1]-c)
    return dx
end

"""
    rk4(f, x0, time_span[, h, args])

4-th order Runge-Kutta. ``f`` is a handler to the
space state funtion, ``x0`` is a vector of initial
conditions, ``time_span`` is a tuple containing
the start and end simulation times and ``h`` is the
integration interval.

The function ``f`` is called as:

    dx = f(t, x, args...)
"""
function rk4(f::Function, x0::Vector{Float64}, time_span; h=(time_span[2]-time_span[1])/1e4, args=())
    t = time_span[1]:h:time_span[2]
    N = length(t)
    n = length(x0)
    # Allocate buffers
    x = Vector{Vector{Float64}}(N)
    for i = 1:N
        x[i] = Vector{Float64}(n)
    end
    xd = Vector{Float64}(n)
    x_aux = Vector{Float64}(n)
    # Simulate
    copy!(x[1], x0)
    @inbounds for i = 1:N-1
        copy!(x[i+1], x[i])
        # 1-st evaluation
        f(t[i], x[i], xd, args...)
        Base.LinAlg.axpy!(h/6, xd, x[i+1])
        # 2-nd evaluation
        copy!(x_aux, x[i])
        Base.LinAlg.axpy!(h/2, xd, x_aux)
        f(t[i]+h/2, x_aux, xd, args...)
        Base.LinAlg.axpy!(h/6, xd, x[i+1])
        # 3-rd evaluation
        copy!(x_aux, x[i])
        Base.LinAlg.axpy!(h/2, xd, x_aux)
        f(t[i]+h/2, x_aux, xd, args...)
        Base.LinAlg.axpy!(h/6, xd, x[i+1])
        # 4-th evaluation
        copy!(x_aux, x[i])
        Base.LinAlg.axpy!(h, xd, x_aux)
        f(t[i]+h, x_aux, xd, args...)
        Base.LinAlg.axpy!(h/6, xd, x[i+1])
    end
    return x
end

# Simulate Rossler
time_span = (0, 1000)
x0 = [1.0, 1.0, 2.0]
args = (0.2, 0.2, 5.7)

@time x = rk4(rossler, x0, time_span, h=0.01, args=args)

x = hcat(x...)'

plot(x)

plot(x[:, 1], x[:,2], x[:, 3])
