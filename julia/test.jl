using Distributions, Plots

include("./modules.jl")

import .SEIRHVDAbstractSpace as SEIRHVDa
using Distributions, Plots

function unit_test_SEIRHVDa(;N = 100_000, ndays = 500)
    #parameters
    stepsPerDay = 5
    chanceInfect = Dict{Tuple{Bool,Bool}, Function}(
        (true, true) => (day) -> 0.3, 
        (true, false) => (day) -> 0.3, 
        (false, true) => (day) -> 0.3, 
        (false, false) => (day) -> 0.3)
    chanceMeet(day) = 0.2
    chanceHospitalize(day) = 0.3
    probDieHosp(day) = 0.3 #prob of recover, prob of die
    vaccinesPerDay(day) = div(N, 100)

    #distributions
    distGetSick() = rand(Normal(5.0, 2.0))
    distsRecover = Dict{Bool,Function}(
        true => () -> rand(Normal(10.0, 2.0)), 
        false => () -> rand(Normal(10.0, 2.0)))
    distsCriticalDie = Dict{Bool,Function}(
        true => () -> rand(Normal(5.0, 2.0)), 
        false => () -> rand(Normal(3.0, 2.0)))
    distHospitalDie() = rand(Normal(10.0, 2.0))
    distHospitalRecover() = rand(Normal(10.0, 2.0))
    distLooseInmunity() = rand(Normal(1_000_000.0, 2.0))

    #create population
    nStates = Dict{Symbol,Int64}(
        :S => N - 80,
        :Sv => 0,
        :E => 0,
        :Ev => 0,
        :Im => 10,
        :Ivm => 0,
        :Ivcr => 0,
        :Icr => 10,
        :R => 10,
        :H => 10,
        :D => 0)
    
    
    population, compartments = SEIRHVDa.build_population(
        nStates, 
        distGetSick, 
        distsRecover[false],
        distsCriticalDie[false], 
        distLooseInmunity)

    model = SEIRHVDa.Model( 
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
        distLooseInmunity;

        stepsPerDay = stepsPerDay)
    
    return SEIRHVDa.run_sim!(model, ndays)
end

unit_test_SEIRHVDa()