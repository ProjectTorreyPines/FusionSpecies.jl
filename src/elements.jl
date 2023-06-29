

struct Element{T<:Particles} <: AbstractElement
    name::String
    symbol::Symbol
    atomic_number::Int64
    mass::Float64
    density::Float64
end
dummy_element = Element{DummyParticles}("dummy", :dummy, 0, 0.0, 0.0)
Element(element::PeriodicTable.Element) = Element{Atoms}(element.name, Symbol(element.symbol), element.number, (element.atomic_mass |> u"kg").val, (element.density |> u"kg/m^3").val)
Element(n, s, a, m, d; type=Atoms) = Element{type}(n, s, a, m, d)

