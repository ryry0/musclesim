module MuscleSim

type HillMuscleModel
    # Activation constants
    dt::Float64
    tau::Float64
    beta::Float64

    # Maximum muscle velocity and force
    V_max::Float64
    F_max::Float64
    L_max::Float64

    # Tendon constants
    K_t::Float64
    L_st::Float64

    # Length and Velocity of the Muscle and Tendon
    L_mt::Vector{Float64}
    V_mt::Vector{Float64}

    # Force, normalized length, and Velocity of the muscle
    F_m::Vector{Float64}
    norm_L_m::Vector{Float64}
    V_m::Vector{Float64}
    activation::Vector{Float64}
    excitation::Vector{Float64}
    time::Vector{Float64}

    # function of time that describes the excitation
    excitation_func
end

function CreateModel(;dt = 0.0, tau = 0.0, beta = 0.0, V_max = 0.0, F_max = 0.0,
                     L_max = 0.0, K_t = 0.0, L_st = 0.0, L_mt::Vector{Float64} =
                     [], V_mt::Vector{Float64} = [], time::Vector{Float64} = [],
                     excitation_func = excite)
return model = HillMuscleModel(
        dt,
        tau,
        beta,
        V_max,
        F_max,
        L_max,
        K_t,
        L_st,
        L_mt,
        V_mt,
        [], # F_m
        [], # norm_L_m
        [], # V_m
        [], # a
        [],
        time, # time
        excitation_func)
end

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

    # U = A*u[n-1] + B*[n-2] + C*[n-3] + D*du*dt
    return U = A*u[end-1] + B*u[end-2] + C*u[end-3] + D*du(t, u[end])*dt
end

#audot :: t -> (t -> b) -> b -> b
function audot(t, u, a; beta=BETA, tau=TAU)
    numerator = u(t) - (beta + (1 - beta) * u(t))*a
    return numerator/tau
end

function excite(t::Float64)
    T0 = 0.5
    T1 = 2
    HIGH = 0.5
    LOW = 0.01
    ret = 0
    if t > T1
        return LOW
    end
    if t > T0
        return HIGH
    end
    return LOW
end

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

function calcV_m(F_0, V_max, F_m, a, norm_L_m)
    y = F_m/(F_0*a*norm_length_tension(norm_L_m))
    return V_max*norm_inv_fv(y)
end

function calcL_m(L_st, K_t, F_m)
    return L_st + F_m/K_t
end

function interp(x_basis::Vector{Float64}, y_basis::Vector{Float64}, t::Float64)
    #print("time $t\n")
    next_index = findfirst(x -> x > t, x_basis)
    if t >= x_basis[end]
        next_index = length(x_basis)
    end
    prev_index = next_index -1
    nearest_time = x_basis[prev_index]
    #=
    print("nearest time $nearest_time\n")
    =#
    avg_constant = t - x_basis[prev_index]

    #=
    print("avg_constant ")
    print(avg_constant)
    print("\n")

    print("1 -avg_constant ")
    print(1- avg_constant)
    print("\n")
    print("1 -avg_constant ")
    print(1- avg_constant)
    print("\n")
    =#
    out = (1 - avg_constant)*y_basis[prev_index] + avg_constant*y_basis[next_index]

    #=
    print("out ")
    print(out)
    print("\n\n")
    =#

    return out
end


function interp_activation(model::HillMuscleModel, t::Float64)
    return interp(model.time, model.activation, t);
end

function interp_length(model::HillMuscleModel, t::Float64)
    return interp(model.time, model.L_mt, t);
end

function interp_velocity(model::HillMuscleModel, t::Float64)
    return interp(model.time, model.V_mt, t);
end

function calcF_mdot(model::HillMuscleModel, t, F_m)
    norm_L_m = (interp_length(model, t) - F_m/model.K_t - model.L_st)/model.L_max
    y =
        F_m/(model.F_max * interp_activation(model, t) * norm_length_tension(norm_L_m))

    V_m = model.V_max*norm_inv_fv(y)

    V_t = interp_velocity(model, t) - V_m

    print("F_m $F_m\n")
    print("numerator ")
    print(interp_length(model, t) - F_m/model.K_t - model.L_st)
    print("\n")
    print("norm Lm $norm_L_m\n")
    print("y $y\n")
    print("norm_inf_fvy $(norm_inv_fv(y))\n")
    print("V_m $V_m\n")
    print("V_t $V_t\n")
    print("\n")

    return F_mdot = model.K_t*V_t
end

function gen_activation!(model::HillMuscleModel)
    model.excitation = map(model.excitation_func, model.time)
    a_dot =
        (t, a) ->
            audot(t, model.excitation_func, a, tau=model.tau, beta=model.beta)

    model.activation =
        foldl((acc, t) ->
            vcat(acc, RK4(acc[end], a_dot, t, model.dt)), 
                model.excitation_func(model.time[1]), model.time) # integrate

    model.activation =
        collect(Iterators.take(model.activation, length(model.time)))
end

function simulate(model::HillMuscleModel)
    gen_activation!(model)
    F_m_dot = (t, F_m) -> calcF_mdot(model, t, F_m)
    model.F_m =
        foldl((acc, t) ->
            vcat(acc, RK4(acc[end], F_m_dot, t, model.dt)), model.time)
end

end
ms = MuscleSim
