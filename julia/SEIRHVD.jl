using Graphs, StatsBase, Distributions, Pipe

include( "./collect_SEIRHVD.jl")
using .Collect_SEIRHVD


include("./SEIRHVDModel.jl")

function build_population( 
        numStates::Dict{Symbol,Int}, # :S, :Sv, :E, :Ev, :Im, :Icr, :Ivm, :Ivcr, :R, :H, :D
        distGetSick::Function, 
        distGetWell::Function, 
        distDie::Function, 
        distLooseInmunity::Function)

    population = Dict{Int,Agent}()
    compartments = Dict{Symbol,Set{Agent}}()

    index = 1
    for i in 1:numStates[:S]
        population[i] = Agent(i, :S, false, nothing, nothing)
    end

    index += numStates[:S]

    for i in index : index + numStates[:Sv]
        population[i] = Agent(i, :S, true, nothing, nothing)
    end

    index += numStates[:Sv]

    for i in index : index + numStates[:E]
        population[i] = Agent(i, :E, false, sample([:Im, :Icr]), distGetSick())
    end

    index += numStates[:E]

    for i in index : index + numStates[:Ev]
        population[i] = Agent(i, :E, true, sample([:Im, :Icr]), distGetSick())
    end

    index += numStates[:Ev]

    for i in index : index + numStates[:Im]
        population[i] = Agent(i, :Im, false, :R, distGetWell())
    end

    index += numStates[:Im]

    for i in index : index + numStates[:Icr]
        population[i] = Agent(i, :Icr, false, :D, distDie())
    end

    index += numStates[:Icr]

    for i in index : index + numStates[:Ivm]
        population[i] = Agent(i, :Im, true, :R, distGetWell())
    end

    index += numStates[:Ivm]

    for i in index : index + numStates[:Ivcr]
        population[i] = Agent(i, :Icr, true, :D, distDie() )
    end

    index += numStates[:Ivcr]

    for i in index : index + numStates[:R]
        population[i] = Agent(i, :R, true, :S, distLooseInmunity())
    end

    index += numStates[:R]

    for i in index : index + numStates[:H]
        population[i] = Agent(i, :H, sample([true, false]), sample([:R, :D]), distDie())
    end

    index += numStates[:H]

    for i in index : index + numStates[:D]
        population[i] = Agent(i, :D, sample([true, false]), nothing, nothing)
    end

    compartments[:S] = Set{Agent}()
    compartments[:E] = Set{Agent}()
    compartments[:Im] = Set{Agent}()
    compartments[:Icr] = Set{Agent}()
    compartments[:R] = Set{Agent}()
    compartments[:H] = Set{Agent}()
    compartments[:D] = Set{Agent}()

    for agent in values(population)
        push!(compartments[agent.status], agent)
    end

    return population, compartments

end

function run_sim!(model::Model, ndays::Int)
    #structures to collect data
    collector = Collector(model.population)
    data = Data()
    
    #since SEIRHVD could run forever, it only can be run for a fixed number of days
    for i in 1:ndays
        run_day!(model, collector)
        collect!(data, collector)
        update!(model)
        restartDay!(collector)
    end

    return data
end

function run_day!(model::Model, collector::Collector)
    for step in model.stepsPerDay
        run_step!(model, collector)
        randVaccinate!(model, collector)
    end
end

function run_step!(model::Model, collector::Collector)
    for agent in values(model.population)
        update!(model, agent)
    end

    for agent in model.compartments[:E]
        maybeGetSick!(model, agent, collector)
    end

    for agent in model.compartments[:Im]
        spread!(model, agent, collector)
        maybeGetWell!(model, agent, collector)
    end

    for agent in model.compartments[:Icr]
        spread!(model, agent, collector)
        if maybeDie!(model, agent, collector)
            continue
        else
            randHospitalize!(model, agent, collector)
        end
    end

      for agent in model.compartments[:R]
        maybeLooseInmunity!(model, agent, collector)
    end

    for agent in model.compartments[:H]
        maybeGetWellOrDie!(model, agent, collector)
    end
end


function maybeInfect!(model::Model, agent::Agent, other::Agent, collector::Collector)
    if rand() < model.currentChanceInfect[agent.isVaccinated, other.isVaccinated]
        infect!(model, other, collector)
    end
end

function infect!(model::Model, agent::Agent, collector::Collector)
    collector.totals[agent.status, agent.isVaccinated] -= 1
    collector.totals[:E, agent.isVaccinated] += 1
    collector.daily[:E, agent.isVaccinated] += 1

    delete!(model.compartments[:S], agent)
    push!(model.compartments[:E], agent)

    agent.status = :E
    agent.daysToChangeState = model.distGetSick()
    agent.nextState = sample([:Im, :Icr])
end

function maybeGetWellOrDie!(model::Model, agent::Agent, collector::Collector)
    if agent.daysToChangeState <= 0
        if agent.nextState == :R
            getWell!(model, agent, collector)
        else
            die!(model, agent, collector)
        end
    end
end

function maybeLooseInmunity!(model::Model, agent::Agent, collector::Collector)
    if agent.daysToChangeState <= 0
        looseInmunity!(model, agent, collector)
    end
end

function looseInmunity!(model::Model, agent::Agent, collector::Collector)
    exchange!(collector, (:R, agent.isVaccinated), (:S, false))

    delete!(model.compartments[:R], agent)
    push!(model.compartments[:S], agent)

    agent.status = :S
    agent.isVaccinated = false
    agent.daysToChangeState = nothing
    agent.nextState = nothing
end

function maybeDie!(model::Model, agent::Agent, collector::Collector) # true if agent dies
    if agent.daysToChangeState <= 0
        die!(model, agent, collector)
        return true
    end
    return false
end

function die!(model::Model, agent::Agent, collector::Collector)
    exchange!(collector, (agent.status, agent.isVaccinated), (:D, agent.isVaccinated))

    delete!(model.compartments[agent.status], agent)
    push!(model.compartments[:D], agent)

    agent.status = :D
    agent.daysToChangeState = nothing
    agent.nextState = nothing
end

function randHospitalize!(model::Model, agent::Agent, collector::Collector)
    if rand() < model.currentChanceHospitalize
        hospitalize!(model, agent, collector)
    end
end

function hospitalize!(model::Model, agent::Agent, collector::Collector)
    exchange!(collector, (agent.status, agent.isVaccinated), (:H, agent.isVaccinated))

    delete!(model.compartments[agent.status], agent)
    push!(model.compartments[:H], agent)

    agent.status = :H
    agent.nextState = rand() < model.currentProbDieHosp ? :D : :R
    agent.daysToChangeState = agent.nextState == :D ? model.distHospitalDie() : model.distHospitalRecover()
end

function maybeGetWell!(model::Model, agent::Agent, collector::Collector)
    if agent.daysToChangeState <= 0
        getWell!(model, agent, collector)
    end
end

function getWell!(model::Model, agent::Agent, collector::Collector)
    exchange!(collector, (agent.status, agent.isVaccinated), (:R, agent.isVaccinated))
    
    delete!(model.compartments[agent.status], agent)
    push!(model.compartments[:R], agent)

    agent.status = :R
    agent.daysToChangeState = model.distLooseInmunity()
end

function maybeGetSick!(model::Model, agent::Agent, collector::Collector)
    if agent.daysToChangeState <= 0
        getSick!(model, agent, collector)
    end
end

function getSick!(model::Model, agent::Agent, collector::Collector)
    exchange!(collector, (agent.status, agent.isVaccinated), (agent.nextState, agent.isVaccinated))

    agent.status = agent.nextState #:Im or :Icr
    if agent.status == :Im
        agent.nextState = :R
        agent.daysToChangeState = model.distsRecover[agent.isVaccinated]()
        delete!(model.compartments[:E], agent)
        push!(model.compartments[:Im], agent)
    else  #:Icr
        agent.nextState = :D
        agent.daysToChangeState = model.distsCriticalDie[agent.isVaccinated]()
        delete!(model.compartments[:E], agent)
        push!(model.compartments[:Icr], agent)
    end
end

# vaccinates model.vaccinesPerDay randomly chosen suceptible agents
function randVaccinate!(model::Model, collector::Collector)
    notVaccinated = filter(a -> !a.isVaccinated, model.compartments[:S])
    sampleSize = minimum([model.vaccinesPerDay(model.currentDay), length(notVaccinated)])
    @pipe model.compartments[:S] |>
        collect |> # INEFICIENT, HOPEFULLY EVENTUALY CHANGE
        sample(_, sampleSize, replace=false) |>
        foreach(a -> vaccinate!(model, a, collector), _)
end

function vaccinate!(model::Model, agent::Agent, collector::Collector)
    exchange!(collector, (agent.status, agent.isVaccinated), (agent.status, true))
    agent.isVaccinated = true
end

function update!(model::Model, agent::Agent)
    if agent.daysToChangeState !== nothing
        agent.daysToChangeState -= 1/model.stepsPerDay
    end
end

function update!(model::Model)
    model.currentDay += 1
    model.currentChanceMeet = model.chanceMeet(model.currentDay) / model.stepsPerDay
    model.currentChanceHospitalize = model.chanceHospitalize(model.currentDay)
    model.currentProbDieHosp = model.probDieHosp(model.currentDay)

    for a1 in [true, false], a2 in [true, false]
        model.currentChanceInfect[a1, a2] = model.chanceInfect[a1, a2](model.currentDay)
    end
end


