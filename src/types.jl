abstract type AbstractSpecies end
abstract type AbstractLoadedSpecies end

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




struct BaseSpecies{T} <:AbstractSpecies
    charge_state :: Float64
    name :: String
    symbol :: Symbol
    element :: AbstractElement
    mass :: Float64
    atomic_number:: Int64
    
end

function BaseSpecies(e::AbstractElement,z::Int64)
        T = get_type_species(z) 
        s = get_species_symbol(e.symbol,z)
        BaseSpecies{T}(z,e.name,s,e,e.mass,e.atomic_number)
end

struct LoadedSpecies{T} <: AbstractLoadedSpecies
    charge_state :: Float64
    name :: String
    symbol :: Symbol
    element :: AbstractElement
    mass :: Float64
    atomic_number:: Int64
    index:: Int64
end

function LoadedSpecies(s::BaseSpecies, idx::Int64)
        T = get_species_type(s) 
       
        LoadedSpecies{T}(s.charge_state,s.name,s.symbol,s.element,s.mass,s.atomic_number,idx)
end

type(s::LoadedSpecies{T}) where T = T 
stype(s::LoadedSpecies)  = split(string(type(s)),".")[end]

struct LoadedSpeciesSet
    list_species :: Vector{AbstractLoadedSpecies}
    dic_species :: Dict{Int64,AbstractLoadedSpecies}
end

function LoadedSpeciesSet(species_registry)
    # check id species
    check_species_index(species_registry)
    
    # make dict with species indexes
    dic_species = Dict{Int64,AbstractLoadedSpecies}()
    for s in [v for v in values(species_registry) if typeof(v) <: AbstractLoadedSpecies]
        dic_species[s.index] = s 
    end
    list_species = [v for v in values(dic_species) if typeof(v) <: AbstractLoadedSpecies]
    for s in values(dic_species)
        list_species[s.index] = s
    end
    LoadedSpeciesSet(list_species ,dic_species)
end

