module MuscleSim
include("jules/Jules.jl")
js = Jules

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
    K_sp::Float64
    L_st::Float64

    # Length and Velocity of the Muscle and Tendon
    L_mt_init::Float64
    L_total::Float64
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
                     L_max = 0.0, K_t = 0.0, K_sp = 0.0, L_st = 0.0, L_mt_init = 0.0, 
                     L_total = 0.0, time::Vector{Float64} = [], excitation_func = excite)
return model = HillMuscleModel(
        dt,
        tau,
        beta,
        V_max,
        F_max,
        L_max,
        K_t,
        K_sp,
        L_st,
        #L_load,
        L_mt_init,
        L_total,
        zeros(length(time)),
        zeros(length(time)),
        zeros(length(time)), # F_m
        [], # norm_L_m
        [], # V_m
        [], # a
        [],
        time, # time
        excitation_func)
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
        if V_m > 1.0
            V_m = 1.0
        end
    else
        V_m = -0.9968
    end
    return V_m
end

function calcV_m(F_0, V_max, F_m, a, norm_L_m)
    y = F_m/(F_0*a*norm_length_tension(norm_L_m))
    return V_max*norm_inv_fv(y)
end

function calcL_m(L_st, K_t, F_m)
    return L_st + F_m/K_t
end

function interp_activation(model::HillMuscleModel, t::Float64)
    return js.interp(model.time, model.activation, t);
end

function interp_length(model::HillMuscleModel, t::Float64)
    return js.interp(model.time, model.L_mt, t);
end

function interp_velocity(model::HillMuscleModel, t::Float64)
    return js.interp(model.time, model.V_mt, t);
end

function calcF_mdot(model::HillMuscleModel, t, F_m, V_mt)
    L_t = model.L_st + F_m/model.K_t
    L_m = (model.L_total - L_t - F_m/model.K_sp)
    L_mt = L_m + L_t

    norm_L_m = L_m/model.L_max

    y =
        F_m/(model.F_max * interp_activation(model, t) * norm_length_tension(norm_L_m))

    V_m = model.V_max*norm_inv_fv(y)

    V_t = V_mt - V_m

    #=
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
    =#

    return F_mdot = model.K_t*V_t
end

function calcFmdot_Vmt(model::HillMuscleModel, t, F_m, V_mt)
    F_mdot = calcF_mdot(model, t, F_m, V_mt)
    V_mt = -F_mdot/K_sp

    return (F_mdot, V_mt)
end

function gen_activation!(model::HillMuscleModel)
    model.excitation = map(model.excitation_func, model.time)
    a_dot = # Create anonymous function for activation.
        (t, a) ->
            audot(t, model.excitation_func, a, tau=model.tau, beta=model.beta)

    model.activation = # integrate to form activation
        foldl((acc, t) ->
            vcat(acc, js.RK4(acc[end], a_dot, t, model.dt)), 
                model.excitation_func(model.time[1]), model.time) # integrate

    model.activation =
        collect(Iterators.take(model.activation, length(model.time)))
end

function simulate(model::HillMuscleModel)
    gen_activation!(model)
    #F_m_dot = (t, F_m) -> calcF_mdot(model, t, F_m)

    #=
    model.F_m =
        foldl((acc, t) ->
            vcat(acc, js.RK4(acc[end], F_m_dot, t, model.dt)), model.time)
        =#

    model.L_mt[1] = model.L_mt_init

    # take initial stab at calculating F_m
    temp_force = collect(0:0.01:2)
    temp_vel = map(norm_inv_fv, temp_force)
    FV = js.interp(temp_vel, temp_force, 0.0)
    activation = interp_activation(model, model.time[1])

    model.F_m[1] = 
        model.F_max * activation * norm_length_tension((model.L_mt[1] - model.L_st)/model.L_max)*FV

    model.V_mt[1] = 0.0

    for i in 1:length(model.time)-1
    model.F_m[i+1] = js.RK4(model.F_m[i], (t, F_m) -> calcF_mdot(model, t, F_m, model.V_mt[i]), model.time[i], model.dt)

    model.V_mt[i+1] = calcF_mdot(model, model.time[i], model.F_m[i], model.V_mt[i])/model.K_sp

    print(model.F_m[i])
    print("\n")

    L_t = model.L_st + model.F_m[i]/model.K_t
    L_m = (model.L_total - L_t - model.F_m[i]/model.K_sp)
    model.L_mt[i+1] = L_m + L_t
    end
end

end
