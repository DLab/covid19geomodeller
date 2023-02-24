include("./modules.jl")
import .SIRAbstractSpace as SIRa, .SIRGraphSpace as SIRg
import .SEIRAbstractSpace as SEIRa, .SEIRGraphSpace as SEIRg
import .SEIRHVDAbstractSpace as SEIRHVDa, .SEIRHVDGraphSpace as SEIRHVDg
import Pandas
using PyCall, Distributions, Graphs


function run_SEIR(
    S::Int,
    E::Int, 
    I::Int, 
    R::Int, 
    chanceMeetPy::PyObject, 
    chanceInfectPy::PyObject,
    tGetWell::Float64,
    tGetSick::Float64,
    stepsPerDay::Int
    ;
    days::Union{Int,Nothing} = nothing,
    startDay::Int = 0,
    isGraphSpace::Bool = false,
    )

    if isGraphSpace
        SEIR = SEIRg
        N = S + E + I + R
        graph = watts_strogatz(N, 4, 0.51)
    else
        SEIR = SEIRa
        graph = nothing
    end

    chanceMeet(day) = pycall(chanceMeetPy, Float64, day)
    chanceInfect(day) = pycall(chanceInfectPy, Float64, day)
    distGetWell() = rand(Normal(tGetWell, 2.0))
    distGetInfected() = rand(Normal(tGetSick, 2.0)) 
    population, exposed, infected = SEIR.build_population(S, E, I, R, ()-> distGetWell()/2, ()-> distGetInfected()/2)
    model = SEIR.Model(population, exposed, infected, chanceInfect, chanceMeet, distGetWell, distGetInfected, stepsPerDay ; graph = graph, startDay = startDay + 1)
    data = SEIR.run_sim!(model, days)

    totals = Pandas.DataFrame(data.totals)
    daily = Pandas.DataFrame(data.daily)

    return Dict(
        "totals" => totals,
        "daily" => daily
    )
end

function run_SIR(
    S::Int,
    I::Int, 
    R::Int, 
    chanceMeetPy::PyObject, 
    chanceInfectPy::PyObject,
    tInfected::Float64,
    stepsPerDay::Int
    ;
    days::Union{Int,Nothing} = nothing,
    startDay::Int = 0,
    isGraphSpace::Bool = false
    )

    if isGraphSpace
        SIR = SIRg
        N = S + I + R
        graph = watts_strogatz(N, 4, 0.51)
    else
        SIR = SIRa
        graph = nothing
    end

    chanceMeet(day) = pycall(chanceMeetPy, Float64, day)
    chanceInfect(day) = pycall(chanceInfectPy, Float64, day)
    distGetWell() = rand(Normal(tInfected, 2))

    population, infected = SIR.build_population(S, I, R, distGetWell)
    model = SIR.Model(
        population, infected, 
        chanceInfect, chanceMeet, distGetWell, 
        stepsPerDay ; 
        graph = graph,
        startDay = startDay + 1)
    data = SIR.run_sim!(model, days)

    totals = Pandas.DataFrame(data.totals)
    daily = Pandas.DataFrame(data.daily)
    
    return Dict(
        "totals" => totals,
        "daily" => daily
    )
end



function run_SEIRHVD(
    nStates::Dict{Any,Any},

    chanceMeetPy::PyObject,
    chanceInfectPy::Dict{Any, Any},
    chanceHospitalizePy::PyObject,
    probDieHospPy::PyObject,
    vaccinesPerDayPy::PyObject,

    tGetSick::Float64,
    tsRecover::Dict{Any, Any},
    tsCriticalDie::Dict{Any, Any},
    tHospitalDie::Float64,
    tHospitalRecover::Float64,
    tLooseInmunity::Float64,

    stepsPerDay::Int,
    days::Int
    ;
    startDay::Int = 0,
    isGraphSpace::Bool = false
    )
    
    if isGraphSpace
        SEIRHVD = SEIRHVDg
        N = sum(values(nStates))
        graph = watts_strogatz(N, 4, 0.51)
    else
        SEIRHVD = SEIRHVDa
        graph = nothing
    end

    chanceMeet(day) = pycall(chanceMeetPy, Float64, day)
    chanceInfect = Dict(
        (true, true) => (day::Int) -> pycall(chanceInfectPy[(true, true)], Float64, day),
        (true, false) => (day::Int) -> pycall(chanceInfectPy[(true, false)], Float64, day),
        (false, true) => (day::Int) -> pycall(chanceInfectPy[(false, true)], Float64, day),
        (false, false) => (day::Int) -> pycall(chanceInfectPy[(false, false)], Float64, day)
    )

    chanceHospitalize(day) = pycall(chanceHospitalizePy, Float64, day)
    probDieHosp(day) = pycall(probDieHospPy, Float64, day)
    vaccinesPerDay(day) = pycall(vaccinesPerDayPy, Int, day)

    #Distributions

    distGetSick() = rand(Normal(tGetSick, 2.0))

    distsRecover = Dict(
        true => () -> rand(Normal(tsRecover[true], 2.0)),
        false => () -> rand(Normal(tsRecover[false], 2.0))
    )

    distsCriticalDie = Dict(
        true => () ->rand(Normal(tsCriticalDie[true], 2.0)),
        false => () -> rand(Normal(tsCriticalDie[false], 2.0))
    )

    distHospitalDie() = rand(Normal(tHospitalDie, 2.0))
    distHospitalRecover() = rand(Normal(tHospitalRecover, 2.0))
    distLooseInmunity() = rand(Normal(tLooseInmunity, 2.0))

    #change nStates to Dict{Symbol,Int}
    nStatesSym = Dict{Symbol,Int}()
    for state in keys(nStates)
        nStatesSym[Symbol(state)] = nStates[state]
    end

    population, compartments = SEIRHVD.build_population(
        nStatesSym, 
        distGetSick, 
        distsRecover[false],
        distsCriticalDie[false], 
        distLooseInmunity)

    model = SEIRHVD.Model( 
        population,
        compartments,

        chanceMeet,
        chanceInfect,
        chanceHospitalize,
        probDieHosp,
        vaccinesPerDay,
        
        distGetSick,
        distsRecover,
        distsCriticalDie,
        distHospitalDie,
        distHospitalRecover,
        distLooseInmunity
        ;
        graph = graph,
        stepsPerDay = stepsPerDay,
        startDay = startDay + 1)

    data = SEIRHVD.run_sim!(model, days)


    return Dict(
        "totals" => data.totals,
        "daily" => data.daily
    )
end
















