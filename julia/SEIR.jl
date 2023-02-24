
include( "./collect_SEIR.jl")
using .Collect_SEIR

include("./SEIRModel.jl")

function build_population(numSuceptible::Int, numInfected::Int, numExposed::Int, numRecovered::Int, distGetWell::Function, distGetInfected::Function)
    population = Dict{Int64,Agent}()
    exposed  = Set{Agent}()
    infected  = Set{Agent}()
    for i in 1:numSuceptible
        s = Agent(i,:S,nothing)
        population[i] = s
    end
    for i in numSuceptible+1:numSuceptible+numExposed
        s = Agent(i,:E,distGetInfected())
        push!(exposed,s)
        population[i] = s
    end
    for i in numSuceptible+numExposed+1:numSuceptible+numExposed+numInfected
        s = Agent(i,:I,distGetWell())
        push!(infected,s)
        population[i] = s
    end
    for i in numSuceptible+numExposed+numInfected+1:numSuceptible+numExposed+numInfected+numRecovered
        s = Agent(i,:R,nothing)
        population[i] = s
    end
    return (population, exposed, infected)
end

function run_sim!(model::Model, ndays::Union{Int,Nothing} = nothing)
    #structures to collect data
    collector = Collector(model.population)
    data = Data()
    
    #run the simulation until there are no more infected or exposed
    if ndays === nothing
        while collector.totals[:I] > 0 || collector.totals[:E] > 0
            run_day!(model, collector)
            collect!(data, collector)
            update!(model)
            restartDay!(collector)
        end
    #run the simulation for a fixed number of days
    else
        for t in 1:ndays
            run_day!(model, collector)
            collect!(data, collector)
            update!(model)
            restartDay!(collector)
        end
    end
    
    return data
end

function update!(model::Model)
    model.currentDay += 1
    model.currentChanceInfect = model.chanceInfect(model.currentDay)
    model.currentChanceMeet = model.chanceMeet(model.currentDay) / model.stepsPerDay
end

function run_day!(model::Model, collector::Collector)
    for step in 1:model.stepsPerDay
        for agent in model.infected
            update!(model, agent)
            spread!(model, agent, collector)
            maybeGetWell!(model, agent, collector)
        end
        for agent in model.exposed
            update!(model, agent)
            maybeGetSick!(model, agent, collector)
        end
    end
end

function update!(model::Model, agent::Agent)
    agent.daysToChangeState -= 1/model.stepsPerDay
end

function maybeInfect!(model::Model, agent::Agent, other::Agent, collector::Collector)
    if rand() < model.chanceInfect(model.currentDay)
        infect!(model, other, collector)
    end
end

function infect!(model::Model, agent::Agent, collector::Collector)
    agent.status = :E
    agent.daysToChangeState = model.distGetInfected()
    push!(model.exposed, agent)
    
    collector.totals[:S] -= 1
    collector.totals[:E] += 1
    collector.daily[:E] += 1
end

function maybeGetSick!(model::Model, agent::Agent, collector::Collector)
    if agent.daysToChangeState <= 0
        getSick!(model, agent, collector)
    end
end

function getSick!(model::Model, agent::Agent, collector::Collector)
    agent.status = :I
    agent.daysToChangeState = model.distGetWell()
    delete!(model.exposed, agent)
    push!(model.infected, agent)

    collector.totals[:E] -= 1
    collector.totals[:I] += 1
    collector.daily[:I] += 1
end

function maybeGetWell!(model::Model, agent::Agent, collector::Collector)
    if agent.daysToChangeState <= 0
        getWell!(model, agent, collector)
    end
end

function getWell!(model::Model, agent::Agent, collector::Collector)
    agent.status = :R
    delete!(model.infected, agent)

    collector.totals[:I] -= 1
    collector.totals[:R] += 1
    collector.daily[:R] += 1
end







