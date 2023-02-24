module Collect_SEIRHVD
using DataFrames

export Collector, Data, collect!, restartDay!, exchange!

## COLLECTOR

mutable struct Collector
    totals::Dict{Tuple{Symbol,Bool}, Int}
    daily::Dict{Tuple{Symbol,Bool}, Int} #new each day, regardless of loses
end

function Collector()
    totals = Dict{Tuple{Symbol,Bool}, Int}()
    daily = Dict{Tuple{Symbol,Bool}, Int}()
    for status in [:S, :E, :Im, :Icr, :H, :R, :D]
        for vaccinated in [false, true]
            totals[status, vaccinated] = 0
            daily[status, vaccinated] = 0
        end
    end
    Collector(totals, daily)
end

function Collector(population::Set)
    collector = Collector()
    for agent in population
        collector.totals[agent.status, agent.isVaccinated] += 1
    end
    return collector
end

function Collector(population::Dict)
    return Collector(Set(values(population)))
end


## DATA

struct Data
    totals::Dict{Tuple{Symbol, Bool}, Vector{Int}}
    daily::Dict{Tuple{Symbol, Bool}, Vector{Int}}
end

function Data()
    totals = Dict{Tuple{Symbol, Bool}, Vector{Int}}()
    daily = Dict{Tuple{Symbol, Bool}, Vector{Int}}()
    for status in [:S, :E, :Im, :Icr, :H, :R, :D]
        for vaccinated in [false, true]
            totals[status, vaccinated] = Int[]
            daily[status, vaccinated] = Int[]
        end
    end
    Data(totals, daily)
end

## COLLECT!

function collect!(data::Data, collector::Collector)
    for status in [:S, :E, :Im, :Icr, :H, :R, :D]
        for vaccinated in [false, true]
            push!(data.totals[status, vaccinated], collector.totals[status, vaccinated])
            push!(data.daily[status, vaccinated], collector.daily[status, vaccinated])
        end
    end
    return data
end

## DAILY

function restartDay!(collector::Collector)
    for status in [:S, :E, :Im, :Icr, :H, :R, :D]
        for vaccinated in [false, true]
            collector.daily[status, vaccinated] = 0
        end
    end
end

## EXCHANGE

function exchange!(collector::Collector, from::Tuple{Symbol,Bool}, to::Tuple{Symbol,Bool})
    collector.totals[from] -= 1
    collector.totals[to] += 1
    collector.daily[to] += 1
end

end #module


