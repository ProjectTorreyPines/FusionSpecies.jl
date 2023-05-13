abstract type AbstractSpecies end
abstract type AbstractLoadedSpecies <: AbstractSpecies end
abstract type AbstractBaseSpecies <: AbstractSpecies end
abstract type BaseSpeciesParameters end 
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


struct BaseSpecies{T} <:AbstractBaseSpecies
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

struct SpeciesSet{V<:Vector{<:LoadedSpecies},U<:Dict{Int64,<:LoadedSpecies}}
    list_species :: V
    dic_species :: U 
    lock ::Vector{Bool}
end

is_set(s::SpeciesSet) = s.lock[1]

LoadedSpeciesSet() = SpeciesSet(Vector{LoadedSpecies}(),Dict{Int64,LoadedSpecies}(), [false])

function LoadedSpeciesSet(species_registry)
    # check id species
    check_species_index(species_registry)
    
    # make dict with species indexes
    dic_species = Dict{Int64,LoadedSpecies}()
    for s in [v for v in values(species_registry) if typeof(v) <: AbstractLoadedSpecies]
        dic_species[s.index] = s 
    end
    list_species = [v for v in values(dic_species) if typeof(v) <: AbstractLoadedSpecies]
    for s in values(dic_species)
        list_species[s.index] = s
    end
    SpeciesSet(list_species ,dic_species, [true])
end

" species parameters in vector format "
struct SpeciesParameters <:BaseSpeciesParameters
    mass :: Vector{Float64}
    Z    :: Vector{Float64}
    idx_e⁻ :: Int64
    all :: Vector{Int64} 
    ions :: Vector{Int64}
    neutrals :: Vector{Int64}
    atoms :: Vector{Int64}
    molecules :: Vector{Int64}
    imp_ions :: Vector{Int64}
    imp_atoms :: Vector{Int64}
    idx_main_ion ::Int64 
    idx_main_atom :: Int64
    species_set :: SpeciesSet
end 



function SpeciesParameters(species_set:: SpeciesSet; kw...)
    @assert is_set(species_set)
    d = Dict()

    d[:mass] = get_species_mass(species_set)
    d[:Z] = get_species_Z(species_set)
    d[:idx_e⁻] = get_electron_index(species_set;enforce=false)
    d[:all] = get_species_index(species_set)
    d[:ions] = get_ions_index(species_set)
    d[:neutrals] = get_neutrals_index(species_set)
    d[:atoms] = get_atoms_index(species_set)
    d[:molecules] = get_molecules_index(species_set)
    d[:imp_ions] = get_imp_ions_index(species_set)
    d[:imp_atoms] = get_imp_atoms_index(species_set)
    d[:idx_e⁻] = get_electron_index(species_set)
    d[:idx_main_ion] = get_main_ion_index(species_set;enforce = false)
    d[:idx_main_atom] = get_main_atom_index(species_set;enforce = false)
    d[:species_set] = species_set
    return SpeciesParameters([d[f] for f in fieldnames(SpeciesParameters)]...)
end
# function SpeciesParameters(; kw...)
#     check_status_species_registry(;kw...)
#     d = Dict()
#     d[:mass] = get_species_mass()
#     d[:Z] = get_species_Z()
#     d[:idx_e⁻] = get_electron_index(;enforce=false)
#     d[:all] = get_species_index()
#     d[:ions] = get_ions_index()
#     d[:neutrals] = get_neutrals_index()
#     d[:atoms] = get_atoms_index()
#     d[:molecules] = get_molecules_index()
#     d[:imp_ions] = get_imp_ions_index()
#     d[:imp_atoms] = get_imp_atoms_index()
#     d[:idx_e⁻] = get_electron_index()
#     d[:idx_main_ion] = get_main_ion_index(;enforce = false)
#     d[:idx_main_atom] = get_main_atom_index(;enforce = false)
#     return SpeciesParameters([d[f] for f in fieldnames(SpeciesParameters)]...)
# end
doc(s::SpeciesParameters) = "nspc=$(length(s.Z))"

get_electron_index(sp::SpeciesParameters) = sp.idx_e⁻ 


get_electron(sp::SpeciesParameters;kw...) = get_electron(sp.species_set; kw...)

get_main_ion_index(sp::SpeciesParameters; kw...) = get_main_ion_index(sp.species_set; kw...)

get_ions_index(sp::SpeciesParameters) = sp.ions
get_species_index(sp::SpeciesParameters) = get_species_index(sp.species_set)
Species = Union{Vector{<:FusionSpecies.AbstractSpecies},Vector{Symbol}, Symbol, FusionSpecies.AbstractSpecies, Int64, Vector{Int64}}
get_species_index(sp::SpeciesParameters, s::Species) = get_species_index(sp.species_set,s)

