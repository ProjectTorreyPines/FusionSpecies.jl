
abstract type AbstractSpecies{T<:ParticleType} end

abstract type AbstractLoadedSpecies{T<:ParticleType} <: AbstractSpecies{T} end
abstract type AbstractBaseSpecies{T} <: AbstractSpecies{T} end


abstract type DummyParticleType <: ParticleType end
abstract type AbstractElement end


abstract type ChargedParticleType <: ParticleType end
abstract type Electron <: ChargedParticleType end
abstract type Ions <: ChargedParticleType end
abstract type Neutrals <: ParticleType end
abstract type Atoms <: Neutrals end
abstract type Molecules <: Neutrals end

abstract type ImpurityMolecule <: Molecules end
abstract type MainMolecule <: Molecules end
abstract type MainAtom <: Atoms end
abstract type ImpurityAtom <: Atoms end
abstract type ImpurityIon <: Ions end
abstract type MainIon <: Ions end

const ChargedParticleSpecies = AbstractSpecies{<:ChargedParticleType}
const ElectronSpecies = AbstractSpecies{<:Electron}
const IonSpecies = AbstractSpecies{<:Ions}
const NeutralSpecies = AbstractSpecies{<:Neutrals}
const AtomSpecies = AbstractSpecies{<:Atoms}
const MoleculeSpecies = AbstractSpecies{<:Molecules}


const IonsAtoms = Union{Ions,Atoms}


struct BaseSpecies{T} <: AbstractBaseSpecies{T}
    charge_state::SpeciesChargeState
    name::String
    symbol::Symbol
    element::AbstractElement
    mass::SpeciesMass
    atomic_number::SpeciesAtomicNumber
end

BaseSpecies(e::AbstractElement, z::ElementAtomicNumber) = BaseSpecies(e, z.value)
function BaseSpecies(e::AbstractElement, z::Int64)
    T = get_type_species(z)
    s = get_species_symbol(e.symbol, z)
    BaseSpecies{T}(SpeciesChargeState(z), e.name, s, e, SpeciesMass(e.mass), SpeciesAtomicNumber(e.atomic_number))
end

struct LoadedSpecies{T<:ParticleType} <: AbstractLoadedSpecies{T}
    charge_state::SpeciesChargeState
    name::String
    symbol::Symbol
    element::AbstractElement
    mass::SpeciesMass
    atomic_number::SpeciesAtomicNumber
    index::SpeciesIndex
    is_active::Vector{Bool}
end

function LoadedSpecies(s::BaseSpecies, idx::Int64)
    T = get_species_particle_type(s)
    LoadedSpecies{T}(s.charge_state, s.name, s.symbol, s.element, s.mass, s.atomic_number, SpeciesIndex(idx), [false])
end

type(s::AbstractLoadedSpecies{T}) where {T} = T
stype(s::LoadedSpecies) = split(string(type(s)), ".")[end]

struct SpeciesSet{T<:AbstractLoadedSpecies}
    list_species::Vector{T}
    dic_species::Dict{Int64,T}
    lock::Vector{Bool}
end
get_elements(species_set::SpeciesSet) = unique([s.element for s in species_set.list_species])

Base.isempty(s::SpeciesSet) = isempty(s.list_species)

function SpeciesSet(v::Vector{T}) where {T<:AbstractLoadedSpecies}
    dic_species = Dict{Int64,T}()
    for s in v
        dic_species[s.index.value] = s
    end
    SpeciesSet(v, dic_species, [true])
end
Base.getindex(s::SpeciesSet, i::Int64) = s.list_species[i]
SpeciesSet() = SpeciesSet{LoadedSpecies}(Vector{LoadedSpecies}(), Dict{Int64,LoadedSpecies}(), [false])

" species parameters in vector format "
struct SpeciesParameters
    mass::SpeciesMasses
    μ::SpeciesReducedMasses
    Z::SpeciesChargeStates
    idx_e⁻::ElectronIndex
    all::SpeciesIndexes
    ions::SpeciesIndexes
    neutrals::SpeciesIndexes
    atoms::SpeciesIndexes
    molecules::SpeciesIndexes
    imp_ions::SpeciesIndexes
    imp_atoms::SpeciesIndexes
    idx_main_ion::MainIonIndex
    idx_main_atom::MainAtomIndex
    species_set::SpeciesSet
end

sdoc(s::SpeciesParameters) = "nspc=$(length(s.Z))"


get_electron(sp::SpeciesParameters; kw...) = get_electron(sp.species_set; kw...)

get_main_ion_index(sp::SpeciesParameters; kw...) = get_main_ion_index(sp.species_set; kw...)

get_ions_index(sp::SpeciesParameters) = sp.ions
get_species_index(sp::SpeciesParameters) = get_species_index(sp.species_set)
Species = Union{Vector{<:AbstractSpecies},Vector{Symbol},Symbol,FusionSpecies.AbstractSpecies,Int64,Vector{Int64}}

#name(species::Species) = name(get_species(species))
get_species_indexes(sp::SpeciesParameters, s::Union{Symbol,Species,Vector,AbstractSpeciesIndexes}) = get_species_indexes(sp.species_set, s)
get_species_indexes(sp::SpeciesParameters, e::AbstractElement) = get_species_indexes(sp.species_set, get_species(sp.species_set, e))


get_species_charge_states(species_parameters::SpeciesParameters) = species_parameters.Z
get_species_masses(species_parameters::SpeciesParameters) = SpeciesMasses(species_parameters.mass)
get_species_reduced_masses(species_parameters::SpeciesParameters) = SpeciesReducedMasses(species_parameters.mass)
export get_species_masses, get_species_reduced_masses, get_species_charge_states