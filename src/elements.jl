using PeriodicTable

abstract type AbstractElement end

struct Element{T} <: AbstractElement 
    name :: String
    symbol :: Symbol
    atomic_number :: Int64
    mass :: Float64 
end

function Element(element::PeriodicTable.Element)
    Element{Atoms}(element.name,Symbol(element.symbol), element.number, (element.atomic_mass |> u"kg").val)
end

function Element(n, s ,a ,m; type= Atoms)
    Element{type}(n,s,a,m)
end
    