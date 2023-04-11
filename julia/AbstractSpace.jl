function random_neighbour(model, agent)
    while (neigh = rand(model.population)[2]) == agent  end
    return neigh
end

function spread!(model, agent, collector)
    other = random_neighbour(model, agent)
    if rand() < model.currentChanceMeet && other.status == :S
        maybeInfect!(model, agent, other, collector)
    end 
end


