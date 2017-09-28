module MuscleSim
include("jules/Jules.jl")
js = Jules

allow_debug = false

function print_debug(message)
    if allow_debug
        print(message)
    end
end
type HillMuscleModel
    start_time::Float64
    end_time::Float64
    dt::Float64

    # Activation constants
    tau::Float64
    beta::Float64
    activ_lower_bound::Float64
    activ_upper_bound::Float64

    # Maximum muscle velocity and force
    V_max::Float64
    F_max::Float64
    L_optimal::Float64
    L_mt_initial::Float64

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
    activation_dot::Vector{Float64}
    excitation::Vector{Float64} #

    time::Vector{Float64}
end

function CreateModel(;start_time = 0.0, end_time = 0.0, dt = 0.0, tau = 0.0, beta = 0.0, V_max = 0.0, F_max = 0.0,
                     L_optimal = 0.0, K_t = 0.0, K_sp = 0.0, L_st = 0.0, L_load = 0.0, L_mt_initial = 0.0,
                     L_total = 0.0, excitation_func = excite, activ_lower_bound = 0.1, activ_upper_bound = 1.0)

    model = HillMuscleModel(
        start_time,
        end_time,
        dt,

        tau,
        beta,
        activ_lower_bound,
        activ_upper_bound,

        V_max,
        F_max,
        L_optimal,
        L_mt_initial,

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
    return 0.0
end

function excite(t::Float64)
    T0 = 1
    T1 = 6
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
    print_debug("F_m $F_m\n")
    print_debug("numerator ")
    print_debug(interp_length(model, t) - F_m/model.K_t - model.L_st)
    print_debug("\n")
    print_debug("norm Lm $norm_L_m\n")
    print_debug("y $y\n")
    print_debug("norm_inf_fvy $(norm_inv_fv(y))\n")
    print_debug("V_m $V_m\n")
    print_debug("V_t $V_t\n")
    print_debug("\n")
    =#

    return F_mdot = model.K_t*V_t
end


# Model params -> activation -> muscle force -> mucle length -> muscle velocity
function calcMuscleVelocity(model::HillMuscleModel, activation, F_m, L_m)
    y = F_m/(model.F_max * activation * norm_length_tension(L_m/model.L_optimal))
    print_debug("calcMuscleVelocity ")
    print_debug("y $y ")
    print_debug("activation $activation ")
    print_debug("F_m $F_m ")
    print_debug("L_m $L_m\n")
    V_m = model.V_max*norm_inv_fv(y)
    return V_m
end

function calcFmdotnew(model::HillMuscleModel, V_m, V_mt)
    V_t = V_mt - V_m
    F_mdot = V_t*model.K_t

    print_debug("calcfmdot ")
    print_debug("V_m $V_m ")
    print_debug("V_mt $V_mt ")
    print_debug("V_t $V_t ")
    print_debug("F_mdot $F_mdot\n")
    return F_mdot
end

function calcV_mtExternal(external_model::HillExternalModel, model::HillMuscleModel, F_mdot)
    V_mt = -F_mdot/external_model.K_load # TODO: SHOULD THIS BE NEGATIVE?
    return V_mt
end

function calcL_mtExternal(external_model::HillExternalModel, model::HillMuscleModel, F_m)
    L_mt = external_model.L_total - (F_m/external_model.K_load + external_model.L_load)
    if (F_m/external_model.K_load + external_model.L_load) > external_model.L_total
        print("$F_m\n")
    end
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

function simulate(model::HillMuscleModel, external_model::HillExternalModel)

    #F_m_dot = (t, F_m) -> calcF_mdot(model, t, F_m)

    #=
    model.F_m =
        foldl((acc, t) ->
            vcat(acc, js.RK4(acc[end], F_m_dot, t, model.dt)), model.time)
        =#


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
        zeros(length_time), #activation_dot
        zeros(length_time), #excitation

        zeros(length_time) #time
        )

    initial_value_outputs.L_mt[1] = model.L_mt_init

    # take initial stab at calculating F_m
    temp_force = collect(0:0.01:2)
    temp_vel = map(norm_inv_fv, temp_force)
    FV = js.interp(temp_vel, temp_force, 0.0)
    initial_value_outputs.activation[1] = model.activ_lower_bound

    initial_value_outputsmodel.F_m[1] =
        model.F_max * activation * norm_length_tension((initial_value_outputsmodel.L_mt[1] - model.L_st)/model.L_max)*FV

    model.V_mt[1] = 0.0

    for i in 1:length(model.time)-1
    initial_value_outputs.F_m[i+1] =
        js.RK4(initial_value_outputs.F_m[i], (t, F_m) -> calcF_mdot(model, t, F_m, model.V_mt[i]), model.time[i], model.dt)

    model.V_mt[i+1] = calcF_mdot(model, model.time[i], model.F_m[i], model.V_mt[i])/model.K_sp

    print_debug(model.F_m[i])
    print_debug("\n")

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
    zeros(length_time), #activation_dot
    zeros(length_time), #excitation

    zeros(length_time) #time
    )

    temp_excitation = model.excitation_func
    model.excitation_func = temp_excitation # t -> 0.1

    initial_value_outputs.L_mt[1] = model.L_mt_initial
    initial_value_outputs.L_t[1] = model.L_st
    initial_value_outputs.L_m[1] = model.L_mt_initial - model.L_st
    initial_value_outputs.activation[1] = js.constrain(model.excitation_func(model.start_time), model.activ_lower_bound, model.activ_upper_bound)


    # take initial stab at calculating F_m
    temp_force = collect(0:0.01:2)
    temp_vel = map(norm_inv_fv, temp_force)
    FV = js.interp(temp_vel, temp_force, 0.0)

    initial_value_outputs.F_m[1] =
        model.F_max * initial_value_outputs.activation[1] *
            norm_length_tension((initial_value_outputs.L_m[1])/model.L_optimal)*FV

    print_debug("F_m $(initial_value_outputs.F_m[1])\n")
    print_debug("L_t $(initial_value_outputs.L_t[1])\n")
    print_debug("L_m $(initial_value_outputs.L_m[1])\n")
    print_debug("L_mt $(initial_value_outputs.L_mt[1])\n")


# with low activation, simulate for a while to get initial F_m
    loopSimulationAB(model, external_model, initial_value_outputs)

#=
# the main simulation
    model.excitation_func = temp_excitation

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

    outputs.activation[1] = model.activ_lower_bound

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


function loopSimulationAB(model::HillMuscleModel, external_model::HillExternalModel, outputs::HillModelOutputs)
    time = model.start_time
    iteration = 1

    for i in 1:4
        outputs.time[iteration] = time
        simulateStep(model, external_model, outputs, iteration, time)

        time += model.dt
        iteration += 1
    end

    while time < model.end_time - model.dt

# calculate iteration + 1 values
        outputs.time[iteration] = time
        simulateStepAB(model, external_model, outputs, iteration, time)

        time += model.dt
        iteration += 1
    end
end

function calcActivationAB(model::HillMuscleModel, previous_value, previous_adot, t)
#### based on excitation, calculate the next activation
# calculate current excitation

    a_dot = # Create anonymous function for activation.
        (t, a) ->
            audot(t, model.excitation_func, a, tau=model.tau, beta=model.beta)

    activation = # integrate to form activation
            js.adams_bashforth_moulton(previous_value, previous_adot, a_dot, t, model.dt)

    activation = js.constrain(activation, model.activ_lower_bound, model.activ_upper_bound)

    return (activation, a_dot(t, previous_value))
end

function simulateStepAB(model::HillMuscleModel, external_model::HillExternalModel,
    outputs::HillModelOutputs, iteration, time)

#### based on excitation, calculate the next activation
# calculate current excitation

    outputs.excitation[iteration] = model.excitation_func(time)

    if iteration == 1
        (outputs.activation[iteration + 1], outputs.activation_dot[iteration + 1]) =
        calcActivation(model,
            js.constrain(model.excitation_func(time), model.activ_lower_bound, model.activ_upper_bound), time)
    else
        (outputs.activation[iteration + 1], outputs.activation_dot[iteration + 1]) =
        calcActivation(model, outputs.activation[iteration], time)
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

    print_debug("\nactual calc ")
    outputs.V_m[iteration + 1] =
        calcMuscleVelocity(model,
            outputs.activation[iteration], outputs.F_m[iteration], outputs.L_m[iteration])

    outputs.F_mdot[iteration + 1] = calcFmdotnew(model, outputs.V_m[iteration], outputs.V_mt[iteration])

    outputs.V_t[iteration + 1] = outputs.V_mt[iteration] - outputs.V_m[iteration]

    slope = (outputs.activation[iteration + 1] - outputs.activation[iteration])/model.dt
    interpolate_activation = t -> slope*t + outputs.activation[iteration] - slope*time

    print_debug("\nintegral calc ")
    outputs.F_m[iteration + 1] = js.adams_bashforth_moulton(
        outputs.F_m[iteration],
        outputs.F_mdot,
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

    outputs.F_m[iteration + 1] = js.constrain(outputs.F_m[iteration + 1], 0.0, model.F_max) #outputs.F_m[iteration + 1])

    F_m = outputs.F_m[iteration + 1]
    print_debug("\nFinal Fm $F_m \n")

    calcLengths(model, external_model, outputs, iteration, time)

#=
 interact with the outside world, calculating:
 L_mt
 V_mt
=#
    calcExternal(model, external_model, outputs, iteration, time)
end

function calcActivation(model::HillMuscleModel, previous_value, t)
#### based on excitation, calculate the next activation
# calculate current excitation

    a_dot = # Create anonymous function for activation.
        (t, a) ->
            audot(t, model.excitation_func, a, tau=model.tau, beta=model.beta)

    activation = # integrate to form activation
            js.RK4(previous_value, a_dot, t, model.dt)

    activation = js.constrain(activation, model.activ_lower_bound, model.activ_upper_bound)

    return (activation, a_dot(t, previous_value))
end

function simulateStep(model::HillMuscleModel, external_model::HillExternalModel,
    outputs::HillModelOutputs, iteration, time)

#### based on excitation, calculate the next activation
# calculate current excitation

    outputs.excitation[iteration] = model.excitation_func(time)

    if iteration == 1
        (outputs.activation[iteration + 1], outputs.activation_dot[iteration + 1]) =
            calcActivation(model,
                js.constrain(model.excitation_func(time), model.activ_lower_bound, model.activ_upper_bound), time)
    else
        (outputs.activation[iteration + 1], outputs.activation_dot[iteration + 1]) =
            calcActivation(model, outputs.activation[iteration], time)
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

    print_debug("\nactual calc ")
    outputs.V_m[iteration + 1] =
        calcMuscleVelocity(model,
            outputs.activation[iteration], outputs.F_m[iteration], outputs.L_m[iteration])

    outputs.F_mdot[iteration + 1] = calcFmdotnew(model, outputs.V_m[iteration], outputs.V_mt[iteration])

    outputs.V_t[iteration + 1] = outputs.V_mt[iteration] - outputs.V_m[iteration]

    slope = (outputs.activation[iteration + 1] - outputs.activation[iteration])/model.dt
    interpolate_activation = t -> slope*t + outputs.activation[iteration] - slope*time

    print_debug("\nintegral calc ")
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

    outputs.F_m[iteration + 1] = js.constrain(outputs.F_m[iteration + 1], 0.0, outputs.F_m[iteration + 1])
    F_m = outputs.F_m[iteration + 1]
    print_debug("\nFinal Fm $F_m \n")

    calcLengths(model, external_model, outputs, iteration, time)
#=
 interact with the outside world, calculating:
 L_mt
 V_mt
=#
    calcExternal(model, external_model, outputs, iteration, time)
end

function calcLengths(model, external_model, outputs, iteration, time)

    outputs.L_t[iteration + 1] = calcL_t(model, outputs.F_m[iteration])
    Keq = external_model.K_load * model.K_t/(external_model.K_load + model.K_t)
    outputs.L_m[iteration + 1] = 
        outputs.L_mt[iteration] - outputs.L_t[iteration]
        #external_model.L_total - (Keq * outputs.F_m[iteration] + model.L_st + external_model.L_load)

end

function calcExternal(model, external_model, outputs, iteration, time)

    outputs.V_mt[iteration + 1] =
        #0
        calcV_mtExternal(external_model, model, outputs.F_mdot[iteration])
        #calcvmt(model, time)

    outputs.L_mt[iteration + 1] =
        #model.L_mt_initial
        calcL_mtExternal(external_model, model, outputs.F_m[iteration])
        #calclmt(model, time)

end

delay = 2
end_ramp = 4
slope = 0.5*10.0^-2

function calclmt(model, time)
    if time < delay
        return model.L_mt_initial
    end

    if time < end_ramp
        return model.L_mt_initial +  (time - delay)*slope
    end

    return model.L_mt_initial + (end_ramp - delay) * slope
end

function calcvmt(model, time)
    if (time > delay) && (time < end_ramp)
        return slope
    end
    return 0.0
end


end
