get_type_species(z::Int64) = get_type_species(float(z))

function get_species_particle_type(s::BaseSpecies)::Type{<:Particles}
    if s.symbol == :D¹⁺
        T = MainIon
    elseif s.symbol == :D⁰
        T = MainAtom
    elseif s.charge_state.value > 0
        return ImpurityIon
    elseif s.charge_state.value == 0
        return ImpurityAtom
    elseif s.charge_state.value < 0
        return Electron
    end
end

get_species_abstract_type(s::AbstractSpecies{T}) where T = AbstractSpecies{<:T}


function get_type_species(z::Float64)
    if z > 0
        return Ions
    elseif z == 0
        return Atoms
    else
        return Electron
    end
end
get_str_charge_state(z::ElementAtomicNumber) = get_str_charge_state(z.value)
function get_str_charge_state(z::Int64)
    if z == 0
        return dic_expo[0]
    elseif z == -1
        return dic_expo["-"]
        # elseif z==1
        #     return dic_expo["+"]
    elseif z > 0
        return join(vcat([dic_expo[d] for d in reverse(digits(z))], [dic_expo["+"]]), "")
    end
end

function get_species_symbol(e_name::Symbol, z)
    return Symbol(string(e_name) * get_str_charge_state(z))
end

get_element(s::Union{AbstractSpecies,AbstractLoadedSpecies}) = return s.element


get_element(s::AbstractElement) = return s


function get_element(e::Symbol)
    @assert e ∈ keys(element_registry) "Cannot find element $e in element_registry:$(keys(element_registry))"
    return element_registry[e]
end

function get_element_type(e::AbstractElement)
    return typeof(e).parameters[1]
end

function get_list_elements()
end

function get_next_species_index()
    if length([k for (k, v) in species_registry if typeof(v) <: AbstractLoadedSpecies]) > 0
        return maximum([k for (k, v) in species_registry if typeof(v) <: AbstractLoadedSpecies]) + 1
    else
        return 1
    end
end

get_next_species_index(species_set::SpeciesSet) = length(species_set.list_species) + 1

get_nspecies(species_set::SpeciesSet) = length(species_set.list_species)


function show_plasma_species()
    list_elements = [getfield(@__MODULE__, n).symbol for n in names(@__MODULE__, all=true) if getfield(@__MODULE__, n) isa AbstractElement]
    list_species = [getfield(@__MODULE__, n).symbol for n in names(@__MODULE__, all=true) if getfield(@__MODULE__, n) isa AbstractSpecies]
    println(string.(list_elements))
    println(string.(list_species))
end

function show_species()
    check_status_species_registry(lock=true)
    show(species_registry["species_set"])
end

function get_species_set(; enforce=true)
    if enforce
        check_status_species_registry(lock=true)
        return species_registry["species_set"]
    elseif "species_set" ∈ keys(species_registry)
        return species_registry["species_set"]
    else
        return nothing
    end
end

function get_nspecies()
    @assert "species_set" ∈ keys(species_registry)
    check_status_species_registry(lock=true)
    return get_nspecies(species_registry["species_set"])
end

" $TYPEDSIGNATURES return the electron species in active species"
function get_electron(species_set::SpeciesSet; enforce=true)
    ss = [s for s in get_species(species_set) if type(s) <: Electron]
    if enforce isa String || enforce
        @assert length(ss) == 1 enforce
    end

@assert length(ss) < 2 "none or more than one electron species found."
    if length(ss) == 1
        return ss[1]
    elseif length(ss) == 0
        return nothing
    elseif length(ss) > 1
        return ss
    end
end

" $TYPEDSIGNATURES return the main ion species in active species"
function get_main_ion(species_set::SpeciesSet; enforce=true)
    ss = [s for s in get_species(species_set) if type(s) <: MainIon]
    if enforce
        @assert length(ss) == 1 "none or more than one main ion species found among: \n $ss"
    end

    if length(ss) == 1
        return ss[1]
    elseif length(ss) == 0
        return nothing
    elseif length(ss) > 1
        return ss
    end
end

" $TYPEDSIGNATURES return the main ion species in active species"
function get_main_atom(species_set::SpeciesSet; enforce=true)
    ss = [s for s in get_species(species_set) if type(s) <: MainAtom]
    if enforce
        @assert length(ss) == 1 "none or more than one main atom species found among: \n $ss"
    end

    if length(ss) == 1
        return ss[1]
    elseif length(ss) == 0
        return 0
    elseif length(ss) > 1
        return ss
    end
end


(Vector{<:LoadedSpecies})(::Vector{Any}) = (Vector{LoadedSpecies})()

get_species(species_set::SpeciesSet{T}, type_species::Type{<:Particles}) where {T} = Vector{T}([s for s in species_set.list_species if type(s) <: type_species])
"$TYPEDSIGNATURES get a list of active species"
get_species(species_set::SpeciesSet{T}) where {T} = species_set.list_species::Vector{T}

"$TYPEDSIGNATURES get a list of active species"
get_all(species_set::SpeciesSet{T}) where {T} = species_set.list_species::Vector{T}

"$TYPEDSIGNATURES get a list of the mass of active species"
get_species_masses(species_set::SpeciesSet) :: SpeciesMasses = SpeciesMasses([s.mass for s in species_set.list_species])

"$TYPEDSIGNATURES get a list of the charge state of active species"

get_species_charge_states(s::AbstractSpecies)::SpeciesChargeStates = SpeciesChargeStates(s.charge_state)

"$TYPEDSIGNATURES get a list of the charge state of active species"
get_species_charge_states(species_set::SpeciesSet)::SpeciesChargeStates = SpeciesChargeStates([s.charge_state for s in species_set.list_species])

"$TYPEDSIGNATURES get a list of the charge state of active species"
get_active(species_set::SpeciesSet)::Vector{Float64} = filter(x -> x.is_active[1], species_set.list_species)

"$TYPEDSIGNATURES get a list of ions among active species"
get_ions(species_set::SpeciesSet) = get_species(species_set, Ions)

"$TYPEDSIGNATURES get a list of ions among active species"
get_imp_ions(species_set::SpeciesSet) = get_species(species_set, ImpurityIon)

"$TYPEDSIGNATURES get a list of ions among active species"
get_imp_atoms(species_set::SpeciesSet) = get_species(species_set, ImpurityAtom)

"$TYPEDSIGNATURES get a list of atoms among active species"
get_atoms(species_set::SpeciesSet) = get_species(species_set, Atoms)

"$TYPEDSIGNATURES get a list of atoms among active species"
get_ions_atoms(species_set::SpeciesSet) = get_species(species_set, IonsAtoms)

"$TYPEDSIGNATURES get a list of molecules among active species"
get_molecules(species_set::SpeciesSet) = get_species(species_set, Molecules)

"$TYPEDSIGNATURES get a list of atoms among active species"
get_neutrals(species_set::SpeciesSet) = get_species(species_set, Neutrals)

"$TYPEDSIGNATURES get a list of index of neutrals among active species"
get_neutrals_indexes(species_set::SpeciesSet)::SpeciesIndexes = get_species_indexes(get_neutrals(species_set))

"$TYPEDSIGNATURES get a list of index of ions among active species"
get_ions_indexes(species_set::SpeciesSet)::SpeciesIndexes = get_species_indexes(get_ions(species_set))

"$TYPEDSIGNATURES get a list of index of ions among active species"
get_main_ion_index(species_set::SpeciesSet; kw...)::MainIonIndex = MainIonIndex(get_species_index(get_main_ion(species_set; kw...)))

"$TYPEDSIGNATURES get the index of main atom among active species"
get_main_atom_index(species_set::SpeciesSet; kw...)::MainAtomIndex = MainAtomIndex(get_species_index(get_main_atom(species_set; kw...)))

"$TYPEDSIGNATURES get a list of index of impurity ions among active species"
get_imp_ions_indexes(species_set::SpeciesSet)::SpeciesIndexes = get_species_indexes(get_imp_ions(species_set))

"$TYPEDSIGNATURES get a list of molecules among active species"
get_molecules_indexes(species_set::SpeciesSet)::SpeciesIndexes = get_species_indexes(get_molecules(species_set))

"$TYPEDSIGNATURES get a list of index of impurity atoms among active species"
get_imp_atoms_indexes(species_set::SpeciesSet)::SpeciesIndexes = get_species_indexes(get_imp_atoms(species_set))

"$TYPEDSIGNATURES get a list of index of atoms among active species"
get_atoms_indexes(species_set::SpeciesSet)::SpeciesIndexes = get_species_indexes(get_atoms(species_set))

"$TYPEDSIGNATURES get a list of index of atoms among active species"
get_ions_atoms_indexes(species_set::SpeciesSet)::SpeciesIndexes = get_species_indexes(get_ions_atoms(species_set))


" $TYPEDSIGNATURES return the index of the electron species"
function get_electron_index(species_set::SpeciesSet; enforce=false)::ElectronIndex
    e = get_electron(species_set)
    if e isa LoadedSpecies{<:Electron}
        idx = ElectronIndex(e.index)
    else
        idx = ElectronIndex(0)
    end
    if enforce isa String || enforce
        @assert !is_missing(idx) (enforce isa String ? enforce : "cannot find electron among loaded species")
    end
    return idx
end



get_ions_Z(species_set::SpeciesSet)::Vector{Float64} = [s.charge_state for s in get_ions(species_set)]


"""$TYPEDSIGNATURES get a list of the index of active species
   $(METHODLIST)
"""
get_species_index(s::Nothing)::SpeciesIndex = SpeciesIndex(missing)
get_species_index(s::Int64)::SpeciesIndex = SpeciesIndex(s)
get_species_index(s::AbstractLoadedSpecies)::SpeciesIndex = SpeciesIndex(s.index)
get_species_index(species_set::SpeciesSet, s::Int64)::SpeciesIndex = get_species_index(get_species(species_set, s))
get_species_indexes(species::Vector{Bool})::SpeciesIndexes = SpeciesIndexes([i for (i, as) in enumerate(species) if as])
get_species_indexes(species::Vector{Int64})::SpeciesIndexes = SpeciesIndexes(species)
get_species_indexes(species::Vector{<:AbstractSpecies})::SpeciesIndexes = length(species) > 0 ? SpeciesIndexes([s.index for s in species]) : SpeciesIndexes()
get_species_indexes(species_set::SpeciesSet, s::Int64)::SpeciesIndexes = get_species_indexes(species_set, get_species(species_set, s))
get_species_indexes(species_set::SpeciesSet)::SpeciesIndexes= get_species_indexes(species_set.list_species)
get_species_indexes(species_set::SpeciesSet, s::Symbol)::SpeciesIndexes = get_species_indexes(get_species(species_set, s))

get_species_indexes(species_set::SpeciesSet, s::AbstractSpecies) = get_species_indexes(species_set, [s])
get_species_indexes(species_set::SpeciesSet, s::AbstractSpeciesIndexes) = get_species_indexes(species_set, s.value)
function get_species_indexes(species_set::SpeciesSet, s::Union{Vector{<:Any},Vector{<:AbstractSpecies},Vector{Symbol},Vector{Int64}})::SpeciesIndexes 
     length(get_species(species_set, s)) > 0 ? SpeciesIndexes([ss.index for ss in get_species(species_set, s)]) : SpeciesIndexes()
end





Base.length(s::SpeciesSet) = length(s.list_species)

Base.iterate(s::SpeciesSet, args...) = iterate(s.list_species, args...)
# " $(TYPEDSIGNATURES) return the species that are of type `element`"
# function get_species(symbol::Symbol; lock=true)
#     @assert "species_set" ∈ keys(species_registry)
#     check_status_species_registry(lock=lock)
#     species =  filter(x->Symbol(x.symbol) == symbol  ,species_registry["species_set"].list_species)
#     @assert length(species) <2 
#     if length(species) == 1
#         return species[1]
#     elseif length(species) == 0
#         return nothing
#     end

# end
# " 
# return the species that are of type `element`
# $(TYPEDSIGNATURES) 
# $(METHODLIST)
# "
# function get_species(v :: Vector{Any}) :: Vector{AbstractLoadedSpecies}
#     out = Vector{AbstractLoadedSpecies}()
#     for s in v 
#          s_ = get_species(s)
#          if s_ isa Array
#             append!(out, s_)
#          else
#             push!(out,s_)
#          end
#     end
#     return out 
# end

" $(TYPEDSIGNATURES) return the species that are of type `element`"
function get_species(species_set::SpeciesSet, s::Symbol; lock=true)
    check_status(species_set, lock=lock)
    if s in keys(species_category)
        return species_category[s](species_set)
    else
        @assert s ∈ [ss.symbol for ss in species_set.list_species] "species $s not loaded. Available species $(name_(species_set))"
        species = filter(x -> Symbol(x.symbol) == s, species_set.list_species)[1]
        return species
    end
end

function get_species(species_set::SpeciesSet, s::Int64; kw...)
    @assert s ∈ [ss.index.value for ss in species_set.list_species] "species $s not loaded. Available species $(name_(species_set))"
    species = filter(x -> x.index.value == s, species_set.list_species)[1]
    return species

end

get_species(species_set::SpeciesSet, element::Element) = filter(x -> x.element == element, species_set.list_species)
get_species(@nospecialize(species_set::SpeciesSet), @nospecialize(s::Vector))::Vector{AbstractSpecies} = vcat([get_species(species_set, sp) for sp in s]...)
get_species(species_set::SpeciesSet, s::SpeciesIndexes)::Vector{AbstractSpecies} =  get_species(species_set, s.value)
get_species(species_set::SpeciesSet, elements::Vector{<:Element}) = filter(x -> x.element ∈ elements, species_set.list_species)
get_species_except(species_set::SpeciesSet, species::LoadedSpecies) = filter(x -> x != species, species_set.list_species)
get_species(species_set::SpeciesSet, element::Element, Z::Int64) = get_species(species_set, get_species_symbol(element.symbol, Z))
get_species(species_set::SpeciesSet, species::LoadedSpecies) = species
get_species(species_set::SpeciesSet, s::Vector{Int64}) = [ss for sss in s for ss in filter(x -> x.index.value == sss, species_set.list_species)]
function get_species(species_set::SpeciesSet, s::Vector{T}) where {T<:Union{AbstractLoadedSpecies,Symbol}}
    out = Vector{LoadedSpecies}()
    for ss in s
        sp = get_species(species_set, ss)
        if sp isa Vector
            append!(out, sp)
        else
            push!(out, sp)
        end
    end
    return out
end
# get_species(s::LoadedSpecies) = [s]
get_species_name(s) = name_.(get_species(s))
get_species_name(s, idx) = name_.(get_species(s, idx))
get_electron_species(species_set::SpeciesSet) = get_electron(species_set)

get_species_parameters(; kw...) = SpeciesParameters(; kw...)

get_nspc(species_set::SpeciesSet) = length(species_set.list_species)

const species_category = Dict(:all => get_all,
    :ions => get_ions,
    :imp_ions => get_imp_ions,
    :main_ion => get_main_ion,
    :main_ions => get_main_ion,
    :imp_atoms => get_imp_atoms,
    :atoms => get_atoms,
    :molecules => get_molecules,
    :neutrals => get_neutrals,
    :electron => get_electron,
    :active => get_active
)

function SpeciesParameters(species_set::SpeciesSet; kw...)
    @assert is_set(species_set)
    d = Dict()

    d[:mass] = get_species_masses(species_set)
    d[:Z] = get_species_charge_states(species_set)
    d[:idx_e⁻] = get_electron_index(species_set; enforce=false)
    d[:all] = get_species_indexes(species_set)
    d[:ions] = get_ions_indexes(species_set)
    d[:neutrals] = get_neutrals_indexes(species_set)
    d[:atoms] = get_atoms_indexes(species_set)
    d[:molecules] = get_molecules_indexes(species_set)
    d[:imp_ions] = get_imp_ions_indexes(species_set)
    d[:imp_atoms] = get_imp_atoms_indexes(species_set)
    d[:idx_main_ion] = get_main_ion_index(species_set; enforce=false)
    d[:idx_main_atom] = get_main_atom_index(species_set; enforce=false)
    d[:species_set] = species_set
    return SpeciesParameters([d[f] for f in fieldnames(SpeciesParameters)]...)
end