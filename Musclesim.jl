module MuscleSim
function RK4(u::Float64, du, t::Float64, dt::Float64)
    k1 = du(t, u)
    k2 = du(t + dt/2, u + dt/2*k1)
    k3 = du(t + dt/2, u + dt/2*k2)
    k4 = du(t + dt, u + dt*k3)
    dU = k1 + 2*k2 + 2*k2 + k4
    return U = u + dU*dt/6.0
end

function Int4(u::Vector{Float64}, du, t, dt)
    A = 18.0/11.0
    B = -9.0/11.0
    C = 2.0/11.0
    D = 6.0/11.0
    
    print(length(u))
    # U = A*u[n-1] + B*[n-2] + C*[n-3] + D*du*dt
    return U = A*u[3] + B*u[2] + C*u[1] + D*du(t, u[4])*dt
end

#audot :: t -> (t -> b) -> b -> b
function audot(t, u, a)
    beta = 0.2
    tau = 0.02
    numerator = u(t) - (beta + (1 - beta) * u(t))*a
    return numerator/tau
end

function excite(t::Float64)
    T0 = 0.5
    T1 = 1
    LVL = 0.5
    ret = 0
    if t > T1
        return 0
    end
    if t > T0
        return LVL
    end
    return 0
end

adot = (t, a) -> audot(t, excite, a)

#norm_length_tension :: norm_length -> norm_force
function norm_length_tension(norm_length)
    A = -3.0508;
    B = 5.9758;
    C = -1.9597;
    if 0.42 <= norm_length <= 1.54
        return A.*(norm_length).^2 + B*norm_length + C;
    end

    return 0
end

#norm_inv_fv :: y -> norm_vel
function norm_inv_fv(y)
    if y > 0
        V_m = 0.995242 * exp(13.8817 * (y - 1.39214)) - 0.996815 * exp(-3.91442 * y)
    else
        V_m = -0.9968
    end
end
end
ms = MuscleSim
