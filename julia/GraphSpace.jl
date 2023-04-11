

function spread!(model::Model, agent::Agent, collector) #graph
    for other in neighbours(model, agent)
        if other.status == :S && rand() < model.currentChanceMeet
            maybeInfect!(model, agent, other, collector)
        end
    end
end

function neighbours(model::Model, agent::Agent) #graph
    return [model.population[i] for i in neighbors(model.graph, agent.id)]
end



