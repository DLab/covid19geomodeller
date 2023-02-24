using Graphs

mutable struct Agent
    id::Int64
    status::Symbol
    daysToChangeState::Union{Float64,Nothing}
end

mutable struct Model
    population::Dict{Int64,Agent}
    graph::Union{SimpleGraph, Nothing}
    exposed::Set{Agent}
    infected::Set{Agent}
    chanceInfect::Function # β , per contact
    chanceMeet::Function # ⍺ , per day 
    distGetWell::Function #function that returns a random number of days to pass from I to R
    distGetInfected::Function #function that returns a random number of days to pass from E to I
    stepsPerDay::Int

    currentDay::Int
    currentChanceInfect::Float64
    currentChanceMeet::Float64 # this one is per step = chanceMeet / stepsPerDay
end

function Model(
    population::Dict{Int64,Agent}, 
    exposed::Set{Agent}, 
    infected::Set{Agent}, 
    chanceInfect::Function, 
    chanceMeet::Function, 
    distGetWell::Function, 
    distGetInfected::Function,  
    stepsPerDay::Int
    ; 
    graph::Union{SimpleGraph, Nothing} = nothing,
    startDay::Int = 1)

    return Model(
        population, 
        graph, exposed, 
        infected, 
        chanceInfect, 
        chanceMeet, 
        distGetWell, 
        distGetInfected, 
        stepsPerDay, 
        startDay, 
        chanceInfect(startDay), 
        chanceMeet(startDay))
end