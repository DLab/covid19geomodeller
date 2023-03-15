mutable struct Agent
    id::Int64
    status::Symbol
    daysToGetWell::Union{Float64,Nothing}
end

mutable struct Model
    population::Dict{Int64, Agent} #id => agent
    infected::Set{Agent} 
    graph::Union{SimpleGraph, Nothing}
    chanceInfect::Function # β , per contact
    chanceMeet::Function # ⍺, per day
    distGetWell::Function
    stepsPerDay::Int

    currentDay::Int
    currentChanceInfect::Float64
    currentChanceMeet::Float64 # this one is per step = chanceMeet / stepsPerDay
end

function Model(
    population::Dict{Int64, Agent}, 
    infected::Set{Agent},
    chanceInfect::Function, 
    chanceMeet::Function, 
    distGetWell::Function, 
    stepsPerDay::Int
    ; 
    graph::Union{SimpleGraph, Nothing} = nothing,
    startDay::Int = 1)

    return Model(
        population, 
        infected, 
        graph, 
        chanceInfect, 
        chanceMeet, 
        distGetWell, 
        stepsPerDay, 
        startDay,
        chanceInfect(startDay),
        chanceMeet(startDay))
end