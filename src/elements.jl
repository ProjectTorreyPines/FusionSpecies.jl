

struct Element{T<:Particles} <: AbstractElement
    name::String
    symbol::Symbol
    atomic_number::ElementAtomicNumber
    mass::ElementMass
    density::ElementDensity
end
dummy_element = Element{DummyParticles}("dummy", :dummy, ElementAtomicNumber(0), ElementMass(0.0), ElementDensity(0.0))
Element(element::PeriodicTable.Element) = Element{Atoms}(element.name, Symbol(element.symbol), ElementAtomicNumber(element.number), ElementMass((element.atomic_mass |> u"kg").val), ElementDensity((element.density |> u"kg/m^3").val))
Element(n, s, a, m, d; kw...) = Element(n, s, ElementAtomicNumber(a), ElementMass(m), ElementDensity(d); kw...)
Element(n, s, a::ElementAtomicNumber, m::ElementMass, d::ElementDensity; type=Atoms) = Element{type}(n, s, a, m, d)


