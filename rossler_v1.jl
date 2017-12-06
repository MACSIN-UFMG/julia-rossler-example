using Plots
plotlyjs()

# Rossler space equation
function rossler(t, state, a, b, c)
    x, y, z = state
    dx = -y - z
    dy = x + a*y
    dz = b + z*(x-c)
    return [dx, dy, dz]
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
function rk4_old(f, x0, time_span; h=(time_span[2]-time_span[1])/1e4, args=())
    t = time_span[1]:h:time_span[2]
    N = length(t)
    n = length(x0)
    # Allocate buffers
    x = zeros(N, n)
    xd = zeros(4, n)
    # Simulate
    x[1, :] = x0
    for i = 1:N-1
        # 1-st evaluation
        xd[1, :] = f(t[i], x[i, :], args...)
        # 2-nd evaluation
        xd[2, :] = f(t[i]+h/2, x[i, :]+h/2*xd[1, :], args...)
        # 3-rd evaluation
        xd[3, :] = f(t[i]+h/2, x[i, :]+h/2*xd[2, :], args...)
        # 4-th evaluation
        xd[4, :] = f(t[i]+h, x[i, :]+h*xd[3, :], args...)
        # Next point
        x[i+1, :] = x[i, :]
        x[i+1, :] += (xd[1, :] + 2*xd[2, :] + 2*xd[3, :] + xd[4, :]) * h/6
    end
    return x
end

# Simulate Rossler
time_span = (0, 1000)
x0 = [1.0, 1.0, 2.0]
args = (0.2, 0.2, 5.7)

@time x = rk4_old(rossler, x0, time_span, h=0.001, args=args)

plot(x)

plot(x[:, 1], x[:,2], x[:, 3])
