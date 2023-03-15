
module Collect_SIR
using DataFrames

export Collector, Data, collect!, restartDay!, exchange!

mutable struct Collector
    totals::Dict{Symbol, Int}
    daily::Dict{Symbol, Int} #new each day, regardless of loses
end

# COLLECTOR

function Collector()
    totals = Dict{Symbol, Int}()
    daily = Dict{Symbol, Int}()
    for status in [:S, :I, :R]
        totals[status] = 0
        daily[status] = 0
    end
    Collector(totals, daily)
end

function Collector(population::Set)
    collector = Collector()
    for agent in population
        collector.totals[agent.status] += 1
    end
    return collector
end

function Collector(population::Dict)
    return Collector(Set(values(population)))
end

# DATA

struct Data
    totals::Dict{Symbol, Vector{Int}}
    daily::Dict{Symbol, Vector{Int}}
end

function Data()
    totals = Dict{Symbol, Vector{Int}}()
    daily = Dict{Symbol, Vector{Int}}()
    for status in [:S, :I, :R]
        totals[status] = Int[]
        daily[status] = Int[]
    end
    return Data(totals, daily)
end

# COLLECT!

function collect!(data::Data, collector::Collector)
    for status in [ :S, :I, :R]
        push!(data.totals[status], collector.totals[status])
        push!(data.daily[status], collector.daily[status])
    end
    return data
end

# RESTART DAY

function restartDay!(collector::Collector)
    for status in [:S, :I, :R]
        collector.daily[status] = 0
    end
end

# EXCHANGE!

function exchange!(collector::Collector, from::Symbol, to::Symbol)
    collector.totals[to] += 1
    collector.totals[from] -= 1
    collector.daily[to] += 1
end

end # module




