using Graphs, Distributions, DataFrames

include( "./collect_SIR.jl")
using .Collect_SIR

include("./SIRModel.jl")


function build_population(numAgents::Int ,numInfected::Int, distGetWell::Function)
    population = Dict{Int64,Agent}()
    infected  = Set{Agent}()
    for i in 1:numAgents
        if i <= numInfected
            s = Agent(i, :I, distGetWell())
            push!(infected,s)
        else
            s = Agent(i, :S, nothing)
        end
        population[i] = s
    end
    return population, infected
end

function build_population(nS::Int, nI::Int, nR::Int, distGetWell::Function)
    population = Dict{Int64,Agent}()
    infected  = Set{Agent}()
    for i in 1:nS
        s = Agent(i, :S, nothing)
        population[i] = s
    end
    for i in nS+1:nS+nI
        s = Agent(i, :I, distGetWell())
        push!(infected,s)
        population[i] = s
    end
    for i in nS+nI+1:nS+nI+nR
        s = Agent(i, :R, nothing)
        population[i] = s
    end
    return population, infected
end


function run_sim!(model::Model, ndays::Union{Int,Nothing} = nothing)
    collector = Collector(model.population)
    data = Data()

    #run the simulation until there are no more infected
    if ndays === nothing
        while collector.totals[:I] > 0
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

function run_day!(model::Model, collector::Collector) #SIR
    for step in 1:model.stepsPerDay
        for agent in model.infected
            update!(model, agent)
            spread!(model, agent, collector)
            maybeGetWell!(model, agent, collector)
        end
    end
end

function update!(model::Model, agent::Agent)
    agent.daysToGetWell -= 1/model.stepsPerDay
end

function maybeInfect!(model::Model, agent::Agent, other::Agent, collector::Collector) 
    if rand() < model.currentChanceInfect
        infect!(model, other, collector)
    end
end

function infect!(model::Model, agent::Agent, collector::Collector) 
    agent.status = :I
    agent.daysToGetWell = model.distGetWell()
    push!(model.infected, agent)

    collector.totals[:S] -= 1
    collector.totals[:I] += 1
    collector.daily[:I] += 1
end

function maybeGetWell!(model::Model, agent::Agent, collector::Collector ) 
    if agent.daysToGetWell <= 0
        getWell!(model, agent, collector)
    end
end

function getWell!(model::Model, agent::Agent, collector::Collector ) 
    agent.status = :R
    delete!(model.infected, agent)

    collector.totals[:I] -= 1
    collector.totals[:R] += 1
    collector.daily[:R] += 1
    
end





