using Graphs

mutable struct Agent
    id::Int
    status::Symbol
    isVaccinated::Bool
    nextState:: Union{Symbol,Nothing}
    daysToChangeState::Union{Float64,Nothing}
end

mutable struct Model
    population::Dict{Int,Agent}
    graph::Union{SimpleGraph, Nothing}
    compartments::Dict{Symbol,Set{Agent}}
    stepsPerDay::Int

    chanceMeet::Function # ⍺ , per day
    chanceInfect::Dict{Tuple{Bool,Bool}, Function} # β, per contact
    chanceHospitalize::Function #probability of getting hospitalized each day while in the state :Icr
    probDieHosp::Function #for :R or :D
    vaccinesPerDay::Function #probability of getting vaccinated each day while in the state :S

    distGetSick::Function
    distsRecover::Dict{Bool,Function} #distribution for the days it takes for an agent with state :Im to pass to :R. The boolean is for vaccinated or not
    distsCriticalDie::Dict{Bool,Function} #boolean for vaccinated or not
    distHospitalDie::Function #distribution for the days it takes for an agent in the hospital to die
    distHospitalRecover::Function #distribution for the days it takes for an agent in the hospital to recover
    distLooseInmunity::Function #distribution for the days it takes for an agent to pass from :R to :S

    currentDay::Int
    currentChanceMeet::Float64 ## this one is per step = chanceMeet / stepsPerDay
    currentChanceInfect::Dict{Tuple{Bool,Bool}, Float64}
    currentChanceHospitalize::Float64
    currentProbDieHosp::Float64

end

function Model(
    population::Dict{Int,Agent},
    compartments::Dict{Symbol,Set{Agent}},

    chanceMeet::Function,
    chanceInfect::Dict{Tuple{Bool,Bool}, Function},
    chanceHospitalize::Function,
    probDieHosp::Function,
    vaccinesPerDay::Function,

    distGetSick::Function,
    distsRecover::Dict{Bool,Function},
    distsCriticalDie::Dict{Bool,Function},
    distHospitalDie::Function,
    distHospitalRecover::Function,
    distLooseInmunity::Function;

    graph::Union{SimpleGraph,Nothing} = nothing,
    stepsPerDay::Int = 5,
    startDay::Int = 1
)

return Model(
    population,
    graph,
    compartments,
    stepsPerDay,

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
    distLooseInmunity,

    startDay, #currentDay
    chanceMeet(startDay), #currentChanceMeet
    Dict{Tuple{Bool,Bool}, Float64}( #currentChanceInfect
        (false,false) => chanceInfect[(false,false)](startDay),
        (false,true) => chanceInfect[(false,true)](startDay),
        (true,false) => chanceInfect[(true,false)](startDay),
        (true,true) => chanceInfect[(true,true)](startDay)
    ),
    chanceHospitalize(startDay),
    probDieHosp(startDay)
)
end