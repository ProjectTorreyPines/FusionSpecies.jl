
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
Base.length(s::SpeciesSet) = length(s.list_species)
Base.iterate(s::SpeciesSet, args...) = iterate(s.list_species, args...)

SpeciesSet() = SpeciesSet{LoadedSpecies}()
SpeciesSet{T}() where T = SpeciesSet{T}(Vector{LoadedSpecies}(), Dict{Int64,LoadedSpecies}(), [false])
check_status(species_set::SpeciesSet; lock=false, message="") = @assert species_set.lock[1] == lock message * " | species_registry : $(species_set.lock[1])"


function add_species(obj::BaseSpecies, species_set::SpeciesSet)
    check_status(species_set)
    tmp = LoadedSpecies(obj, get_next_species_index(species_set))
    @assert tmp ∉ species_set.list_species "Species $obj already added.... \n List of current species in species set: $(species_set.list_species)"
    push!(species_set.list_species, tmp)
end

function add_species(obj::Symbol, species_set::SpeciesSet) 
    @assert obj ∈ keys(element_registry) ||  obj ∈ keys(species_registry) " cannot find the species/element: $obj ...\n Available elements : $(keys(element_registry)) \n Available species: $(keys(species_registry))"
    add_species(getfield(@__MODULE__, obj), species_set)
end

function get_species(obj::Symbol)
    species_set = SpeciesSet()
    add_species(getfield(@__MODULE__, obj), species_set)
    setup_species!(species_set)
    return species_set
end

function get_species(objs::Vector{Symbol})
    species_set = SpeciesSet()
    for obj in objs
        add_species(getfield(@__MODULE__, obj), species_set)
    end
    setup_species!(species_set)
    return species_set
end

function add_species(el::Element, species_set::SpeciesSet)
    check_status(species_set)
    for s in el.species
        add_species(s, species_set)
    end
end

macro species_set(objs...)
    species_set = SpeciesSet()
    for obj in objs
        add_species(obj, species_set)
    end
    setup_species!(species_set)
    return species_set
end


function setup_species!(species_set::SpeciesSet)
    list_species = species_set.list_species
    @assert length(unique(list_species)) == length(list_species)

    list_idx = [v.index for v in species_set.list_species]
    @assert length(unique(list_idx)) == length(list_idx)
    species_set.lock[1] = true
end


function check_species_index(species_set::SpeciesSet)
    for (k, v) in species_set.dic_species
        @assert k == v.index
    end
    indexes = sort([v.index for (k, v) in species_set.dic_species])
    # check that indexes start at 1
    @assert minimum(indexes) == 1
    # check that indexes are incremental by 1
    @assert minimum(diff(indexes)) == maximum(diff(indexes)) == 1
end

function Base.show(io::IO, ::MIME"text/plain", species_set::SpeciesSet)
    print(io, stringstyled("⟦", color=20) * (isempty(species_set) ? "∅" : join(name_.(species_set.list_species), "; ")) * stringstyled("⟧", color=20))
end

function Base.show(io::IO, species_set::SpeciesSet{T}) where {T}
    print(io, stringstyled("⟦", color=20) * (isempty(species_set) ? "∅" : join(name_.(species_set.list_species), "; ")) * stringstyled("⟧", color=20))
end

name_(species_set::SpeciesSet) = join(name_.(species_set.list_species), " ")
inline_summary(species_set::SpeciesSet) = join([name_(species) * "[$(species.index.value)]" for species in species_set.list_species], " ")

import_species!(species_set::SpeciesSet, mod::Module; force_import::Bool=false, verbose=true) = import_species!(species_set, mod, force_import; verbose)
function import_species!(species_set::SpeciesSet, mod::Module, force_import::Bool; verbose=true)
    list = []
    expr = quote end
    for (i, s) in enumerate(species_set.list_species)
        if !hasproperty(mod, s.symbol) || force_import
            push!(list, s)
            setproperty!(mod, s.symbol, species_set.list_species[i])
        else
            println("warning: cannot import species $(name_(s)) into module $mod ...")
        end
    end
    for e in get_elements(species_set)
        if !hasproperty(mod, e.symbol) || force_import
            push!(list, e)
            setproperty!(mod, e.symbol, get_elements(species_set))
        else
            println("warning: cannot import element $(name_(e)) into module $mod ...")
        end
    end
    msg = "importing species into module `$mod`: $(join(name_.(list)," ")) "
    verbose ? println(msg) : return msg
end

macro import_species(args...)
    aargs, kwargs = convert_macro_kwargs(args)
    :force_import ∉ keys(kwargs) ? kwargs[:force_import] = false : nothing
    @assert length(aargs) == 1 && length(kwargs) == 1
    force_import = kwargs[:force_import]
    :(import_species!($(aargs[1]), @__MODULE__, $force_import))
end

function set_main_species!(species_set::SpeciesSet, s) 
    for s in get_species(species_set, s)
    s.is_main_species.bool = true 
    end
end

# ---------------------------------------------------------------------------------------------- # 
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

Base.ones(sp::FusionSpecies.SpeciesParameters) = ones(get_nspecies(sp.species_set))