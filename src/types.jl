abstract type BaseSpecies end
abstract type BaseActiveSpecies end
abstract type BaseElement end
abstract type Atom end
abstract type Dummy end
abstract type Electron end
abstract type Molecule end
abstract type Neutral <: Atom end
abstract type Ion <: Atom end

struct Element{T} <: BaseElement 
    name :: String
    symbol :: Symbol
    atomic_number :: Int64
    mass :: Float64 

end

    function Element(element::PeriodicTable.Element)
        Element{Atom}(element.name,Symbol(element.symbol), element.number, (element.atomic_mass |> u"kg").val)
    end
    function Element(n, s ,a ,m; type= Atom)
        Element{type}(n,s,a,m)
    end
    
struct Species{T} <:BaseSpecies
    charge_state :: Float64
    name :: String
    symbol :: Symbol
    element :: BaseElement
    mass :: Float64
    atomic_number:: Int64
    
end

function Species(e::BaseElement,z::Int64)
        T = get_type_species(z) 
        s = get_species_symbol(e.symbol,z)
        Species{T}(z,e.name,s,e,e.mass,e.atomic_number)
end

struct ActiveSpecies{T} <: BaseActiveSpecies
    charge_state :: Float64
    name :: String
    symbol :: Symbol
    element :: BaseElement
    mass :: Float64
    atomic_number:: Int64
    index:: Int64
end

function ActiveSpecies(s::Species, idx::Int64)
        T = get_type_species(s.charge_state) 
        ActiveSpecies{T}(s.charge_state,s.name,s.symbol,s.element,s.mass,s.atomic_number,idx)
end

struct ActiveSpeciesSet
    list_species :: Vector{BaseActiveSpecies}
    dic_species :: Dict{Int64,BaseActiveSpecies}
end

function ActiveSpeciesSet(active_species_registry)
    # check id species
    check_active_species_index(active_species_registry)
    
    # make dict with species indexes
    dic_species = Dict{Int64,BaseActiveSpecies}()
    for s in [v for v in values(active_species_registry) if typeof(v) <: BaseActiveSpecies]
        dic_species[s.index] = s 
    end
    list_species = [v for v in values(dic_species) if typeof(v) <: BaseActiveSpecies]
    for s in values(dic_species)
        list_species[s.index] = s
    end
    ActiveSpeciesSet(list_species ,dic_species)
end