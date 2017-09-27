module MuscleSim
include("jules/Jules.jl")
js = Jules

type HillMuscleModel
    start_time::Float64
    end_time::Float64
    dt::Float64

    # Activation constants
    tau::Float64
    beta::Float64

    # Maximum muscle velocity and force
    V_max::Float64
    F_max::Float64
    L_optimal::Float64

    # Tendon constants
    K_t::Float64
    L_st::Float64

    # function of time that describes the excitation
    excitation_func
end


# for attaching to compliant spring
type HillExternalModel
    L_load::Float64
    K_load::Float64
    L_total::Float64
end

type HillModelOutputs
    # Length and Velocity of the Muscle and Tendon
    L_mt::Vector{Float64} #
    V_mt::Vector{Float64} #

    # Force, normalized length, and Velocity of the muscle
    F_m::Vector{Float64} #
    F_mdot::Vector{Float64} #
    L_m::Vector{Float64} #
    V_m::Vector{Float64} #

    V_t::Vector{Float64} #
    L_t::Vector{Float64} #

    activation::Vector{Float64} #
    excitation::Vector{Float64} #

    time::Vector{Float64}
end

function CreateModel(;start_time = 0.0, end_time = 0.0, dt = 0.0, tau = 0.0, beta = 0.0, V_max = 0.0, F_max = 0.0,
                     L_optimal = 0.0, K_t = 0.0, K_sp = 0.0, L_st = 0.0, L_load = 0.0,
                     L_total = 0.0, excitation_func = excite)

    model = HillMuscleModel(
        start_time,
        end_time,
        dt,

        tau,
        beta,

        V_max,
        F_max,
        L_optimal,

        K_t,
        L_st,
        excitation_func)

    external_model = HillExternalModel(
        L_load,
        K_sp,
        L_total)

    return (model, external_model)
end

#audot :: t -> (t -> b) -> b -> b
function audot(t, u, a; beta=BETA, tau=TAU)
    numerator = u(t) - (beta + (1 - beta) * u(t))*a
    return numerator/tau
end

function zero_excite(t::Float64)
    return 0
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

function calcL_t(model::HillMuscleModel, F_m)
    return model.L_st + F_m/model.K_t
end

function interp_activation(model::HillMuscleModel, t::Float64)
    return js.interp(model.time[1:length(model.activation)], model.activation, t);
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


# Model params -> activation -> muscle force -> mucle length -> muscle velocity
function calcMuscleVelocity(model::HillMuscleModel, activation, F_m, L_m)
    y = F_m/(model.F_max * activation * norm_length_tension(L_m/model.L_optimal))
    V_m = model.V_max*norm_inv_fv(y)
    return V_m
end

function calcFmdotnew(model::HillMuscleModel, V_m, V_mt)
    V_t = V_mt - V_m
    F_mdot = V_t*model.K_t
    return F_mdot
end

function calcV_mtExternal(external_model::HillExternalModel, model::HillMuscleModel, F_mdot)
    V_mt = -F_mdot/external_model.K_load
    return V_mt
end

function calcL_mtExternal(external_model::HillExternalModel, model::HillMuscleModel, F_m)
    L_mt = external_model.L_total - (F_m/external_model.K_load + external_model.L_load)
    return L_mt
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

function simulateExternal(model::HillMuscleModel, external_model::HillExternalModel)

length_time = length(model.start_time:model.dt:model.end_time)

initial_value_outputs = HillModelOutputs(
    zeros(length_time), #L_mt
    zeros(length_time), #V_mt

    zeros(length_time), #F_m
    zeros(length_time), #F_mdot
    zeros(length_time), #L_m
    zeros(length_time), #V_m

    zeros(length_time), #V_t
    zeros(length_time), #L_t

    zeros(length_time), #activation
    zeros(length_time), #excitation

    zeros(length_time) #time
    )

    temp_excitation = model.excitation_func
    model.excitation_func = zero_excite

# with low activation, simulate for a while to get initial F_m
    loopSimulation(model, external_model, initial_value_outputs)

#=
outputs = HillModelOutputs(
    zeros(length_time), #L_mt
    zeros(length_time), #V_mt

    zeros(length_time), #F_m
    zeros(length_time), #F_mdot
    zeros(length_time), #L_m
    zeros(length_time), #V_m

    zeros(length_time), #V_t
    zeros(length_time), #L_t

    zeros(length_time), #activation
    zeros(length_time), #excitation

    zeros(length_time) #time
    )


# loop over simulation, then allow to settle
    loopSimulation(model, external_model, outputs)

    =#
return initial_value_outputs
end

function loopSimulation(model::HillMuscleModel, external_model::HillExternalModel, outputs::HillModelOutputs)
    time = model.start_time
    iteration = 1

    while time < model.end_time - model.dt

# calculate iteration + 1 values
        outputs.time[iteration] = time
        simulateStep(model, external_model, outputs, iteration, time)

        time += model.dt
        iteration += 1
    end
end

function calcActivation(model::HillMuscleModel, previous_value, t)
#### based on excitation, calculate the next activation
# calculate current excitation

    a_dot = # Create anonymous function for activation.
        (t, a) ->
            audot(t, model.excitation_func, a, tau=model.tau, beta=model.beta)

    activation = # integrate to form activation
            js.RK4(previous_value, a_dot, t, model.dt)

    activation = js.constrain(activation, 0.1, 1.0)

    return activation
end

function simulateStep(model::HillMuscleModel, external_model::HillExternalModel, 
    outputs::HillModelOutputs, iteration, time)

#### based on excitation, calculate the next activation
# calculate current excitation

    outputs.excitation[iteration] = model.excitation_func(time)

    if iteration == 1
        outputs.activation[iteration + 1] = calcActivation(model, model.excitation_func(time), time)
    else
        outputs.activation[iteration + 1] = calcActivation(model, outputs.activation[iteration], time)
    end

#=
 use the activation to calculate the next muscle model state
 this includes:
 F_m = F_t -
 F_mdot = F_tdot -
 V_m -
 V_t - 
 L_t -
 L_m -
 =#

    outputs.V_m[iteration + 1] =
        calcMuscleVelocity(model,
            outputs.activation[iteration], outputs.F_m[iteration], outputs.L_m[iteration])

    outputs.F_mdot[iteration + 1] = calcFmdotnew(model, outputs.V_m[iteration], outputs.V_mt[iteration])

    outputs.V_t[iteration + 1] = outputs.V_mt[iteration] - outputs.V_m[iteration]

    slope = (outputs.activation[iteration + 1] - outputs.activation[iteration])/model.dt
    interpolate_activation = t -> slope*t + outputs.activation[iteration] - slope*time

    outputs.F_m[iteration + 1] = js.RK4(
        outputs.F_m[iteration], 
        (t, u) -> calcFmdotnew(
                    model, 
                    calcMuscleVelocity(
                        model, 
                        interpolate_activation(t), 
                        u, 
                        outputs.L_m[iteration]),
                    outputs.V_mt[iteration]),
        time, 
        model.dt)

    outputs.L_t[iteration + 1] = calcL_t(model, outputs.F_m[iteration])
    outputs.L_m[iteration + 1] = outputs.L_mt[iteration] - outputs.L_t[iteration]

#=
 interact with the outside world, calculating:
 L_mt
 V_mt
=#

    outputs.V_mt[iteration + 1] = calcV_mtExternal(external_model, model, outputs.F_mdot[iteration])
    outputs.L_mt[iteration + 1] = calcL_mtExternal(external_model, model, outputs.F_m[iteration])
end
end
