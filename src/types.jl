abstract type AbstractSpecies end
abstract type AbstractActiveSpecies end

abstract type Dummy end
abstract type Particles end
abstract type ChargedParticles <: Particles end
abstract type Electron <: ChargedParticles end
abstract type Ions <: ChargedParticles end
abstract type Neutrals <: Particles  end
abstract type Atoms <: Neutrals end
abstract type Molecules <: Neutrals end
abstract type ImpurityMolecule <: Molecules end
abstract type MainMolecule <: Molecules end
abstract type MainAtom <: Atoms end
abstract type ImpurityAtom <: Atoms end
abstract type ImpurityIon <: Ions end
abstract type MainIon <: Ions end




struct Species{T} <:AbstractSpecies
    charge_state :: Float64
    name :: String
    symbol :: Symbol
    element :: AbstractElement
    mass :: Float64
    atomic_number:: Int64
    
end

function Species(e::AbstractElement,z::Int64)
        T = get_type_species(z) 
        s = get_species_symbol(e.symbol,z)
        Species{T}(z,e.name,s,e,e.mass,e.atomic_number)
end

struct ActiveSpecies{T} <: AbstractActiveSpecies
    charge_state :: Float64
    name :: String
    symbol :: Symbol
    element :: AbstractElement
    mass :: Float64
    atomic_number:: Int64
    index:: Int64
end

function ActiveSpecies(s::Species, idx::Int64)
        T = get_species_type(s) 
       
        ActiveSpecies{T}(s.charge_state,s.name,s.symbol,s.element,s.mass,s.atomic_number,idx)
end

type(s::ActiveSpecies{T}) where T = T 
stype(s::ActiveSpecies)  = split(string(type(s)),".")[end]

struct ActiveSpeciesSet
    list_species :: Vector{AbstractActiveSpecies}
    dic_species :: Dict{Int64,AbstractActiveSpecies}
end

function ActiveSpeciesSet(active_species_registry)
    # check id species
    check_active_species_index(active_species_registry)
    
    # make dict with species indexes
    dic_species = Dict{Int64,AbstractActiveSpecies}()
    for s in [v for v in values(active_species_registry) if typeof(v) <: AbstractActiveSpecies]
        dic_species[s.index] = s 
    end
    list_species = [v for v in values(dic_species) if typeof(v) <: AbstractActiveSpecies]
    for s in values(dic_species)
        list_species[s.index] = s
    end
    ActiveSpeciesSet(list_species ,dic_species)
end

