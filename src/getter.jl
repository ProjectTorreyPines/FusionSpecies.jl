get_type_species(z::Int64) = get_type_species(float(z))

function get_species_type(s::BaseSpecies)
    if s.symbol == :D¹⁺
        T = MainIon
    elseif s.symbol == :D⁰
            T = MainAtom
    elseif s.charge_state > 0
        return ImpurityIon
    elseif s.charge_state == 0
        return ImpurityAtom
    elseif s.charge_state <1
        return Electron
    end
end

function get_type_species(z::Float64)
    if z>0
        return Ions
    elseif z == 0
        return Atoms
    else
        return Electron
    end
end

function get_str_charge_state(z::Int64)
    if z == 0
        return dic_expo[0]
    elseif z == -1
        return dic_expo["-"]
    # elseif z==1
    #     return dic_expo["+"]
    elseif z>0
        return join(vcat([dic_expo[d] for d in reverse(digits(z))],[dic_expo["+"]]),"")
    end
end

function get_species_symbol(e_name::Symbol,z)
    return Symbol(string(e_name) * get_str_charge_state(z)) 
end
function get_element(s::AbstractSpecies)
    return s.element
end

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
    if length([k for (k,v) in species_registry if typeof(v) <: AbstractLoadedSpecies])>0
        return maximum([k for (k,v) in species_registry if typeof(v) <: AbstractLoadedSpecies])+1
    else
        return 1
    end
end



function get_nspecies(species_set::LoadedSpeciesSet)
    return length(species_set.dic_species)
end

function show_plasma_species()
    list_elements =  [getfield(@__MODULE__,n).symbol for n in names(@__MODULE__, all=true) if getfield(@__MODULE__,n) isa AbstractElement]
    list_species =  [getfield(@__MODULE__,n).symbol for n in names(@__MODULE__, all=true) if getfield(@__MODULE__,n) isa AbstractSpecies]
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
function get_electron(; enforce=true)
    ss = [s for s in get_species() if type(s) <:Electron]
    if enforce
        @assert length(ss) == 1 "none or more than one electron species found."
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
function get_main_ion(; enforce=true)
    ss = [s for s in get_species() if type(s) <: MainIon]
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
function get_main_atom(; enforce=true)
    ss = [s for s in get_species() if type(s) <:MainAtom]
    if enforce
        @assert length(ss) == 1 "none or more than one main atom species found among: \n $ss"
    end

    if length(ss) == 1
        return ss[1]
    elseif length(ss) == 0
        return nothing
    elseif length(ss) > 1
        return ss
    end
end

"$TYPEDSIGNATURES get a list of active species"
get_species() ::Vector{<:AbstractLoadedSpecies} = get_species_set().list_species

"$TYPEDSIGNATURES get a list of the mass of active species"
get_species_mass() ::Vector{Float64} = [s.mass for s in get_species_set().list_species]

"$TYPEDSIGNATURES get a list of the charge state of active species"
get_species_Z(s::AbstractSpecies) = s.charge_state

"$TYPEDSIGNATURES get a list of the charge state of active species"
get_species_Z() = [s.charge_state for s in get_species_set().list_species]

"$TYPEDSIGNATURES get a list of ions among active species"
get_ions() ::Vector{<:LoadedSpecies}  = [s for s in get_species_set().list_species if type(s) <: Ions]

"$TYPEDSIGNATURES get a list of ions among active species"
get_imp_ions() ::Vector{<:AbstractLoadedSpecies}  = [s for s in get_species_set().list_species if type(s) <: ImpurityIon]

"$TYPEDSIGNATURES get a list of ions among active species"
get_imp_atoms() ::Vector{<:AbstractLoadedSpecies}  = [s for s in get_species_set().list_species if type(s) <: ImpurityAtom]

"$TYPEDSIGNATURES get a list of atoms among active species"
get_atoms() ::Vector{<:AbstractLoadedSpecies} = [s for s in get_species_set().list_species if type(s) <: Atoms]

"$TYPEDSIGNATURES get a list of molecules among active species"
get_molecules() ::Vector{<:AbstractLoadedSpecies} = [s for s in get_species_set().list_species if type(s) <: Molecules]


"$TYPEDSIGNATURES get a list of atoms among active species"
get_neutrals() = [s for s in get_species_set().list_species if type(s) <: Neutrals]

"$TYPEDSIGNATURES get a list of index of neutrals among active species"
get_neutrals_index() :: Vector{Int64} = get_species_index(get_neutrals())

"$TYPEDSIGNATURES get a list of index of ions among active species"
get_ions_index() :: Vector{Int64} = get_species_index(get_ions())

"$TYPEDSIGNATURES get a list of index of ions among active species"
get_main_ion_index(;kw...) :: Int64 = get_species_index(get_main_ion(;kw...))

"$TYPEDSIGNATURES get the index of main atom among active species"
get_main_atom_index(;kw...) :: Int64 = get_species_index(get_main_atom(;kw...))

"$TYPEDSIGNATURES get a list of index of ions among active species"
get_imp_ions_index() :: Vector{Int64} = get_species_index(get_imp_ions())

"$TYPEDSIGNATURES get a list of molecules among active species"
get_molecules_index() ::Vector{Int64} = get_species_index(get_molecules())

"$TYPEDSIGNATURES get a list of index of impurity atoms among active species"
get_imp_atoms_index() :: Vector{Int64} = get_species_index(get_imp_atoms())

"$TYPEDSIGNATURES get a list of index of atoms among active species"
get_atoms_index() :: Vector{Int64} = get_species_index(get_atoms())


" $TYPEDSIGNATURES return the index of the electron species"
function get_electron_index(;enforce=true)
    e = get_electron(;enforce=enforce)
    if e isa LoadedSpecies{<:Electron}
        return e.index
    elseif e isa Nothing
        error("no electron found")
    end
end


get_ions_Z() = [s.charge_state for s in get_ions()]


"""$TYPEDSIGNATURES get a list of the index of active species
   $(METHODLIST)
"""
get_species_index(s::AbstractLoadedSpecies) = s.index
get_species_index(species::Vector{Bool}) = Vector{Int64}([i for (i,as) in enumerate(species) if as])
get_species_index(species::Vector{Int64}) = species
get_species_index(species::Vector{<:AbstractLoadedSpecies}) = Vector{Int64}([s.index for s in species])
get_species_index() ::Vector{Int64} = get_species_index(get_species())

" $(TYPEDSIGNATURES) return the species that are of type `element`"
function get_species(symbol::Symbol; lock=true)
    @assert "species_set" ∈ keys(species_registry)
    check_status_species_registry(lock=lock)
    species =  filter(x->Symbol(x.symbol) == symbol  ,species_registry["species_set"].list_species)
    @assert length(species) <2 
    if length(species) == 1
        return species[1]
    elseif length(species) == 0
        return nothing
    end

end
" 
return the species that are of type `element`
$(TYPEDSIGNATURES) 
$(METHODLIST)
"
function get_species(v :: Vector{Any}) :: Vector{AbstractLoadedSpecies}
    out = Vector{AbstractLoadedSpecies}()
    for s in v 
         s_ = get_species(s)
         if s_ isa Array
            append!(out, s_)
         else
            push!(out,s_)
         end
    end
    return out 
end
    

get_species(element::Element) = filter(x->x.element == element  ,get_species_set().list_species)
get_species_except(species::LoadedSpecies) = filter(x->x != species, get_species_set().list_species)

get_species(element::Element, Z::Int64) = get_species(get_species_symbol(element.symbol,Z))
function get_species(s::Vector{T}) where T<:Union{AbstractLoadedSpecies,Symbol}
    out = Vector{LoadedSpecies}()
    for ss in s 
        append!(out,get_species(ss))
    end
    return out
end
get_species(s::LoadedSpecies) = [s]
name(s::LoadedSpecies) = string(s.symbol)
get_species_name(s) = name.(get_species(s))

get_species(is_species_active::Vector{Bool}) =  [get_species(i) for (i,b) in enumerate(is_species_active) if b]
get_species(species_index::Vector{Int64}; kw...) = [get_species(s; kw...) for s in species_index]

"""  return the species with index `species_index`
    $TYPEDSIGNATURES
    $METHODLIST 
"""
function get_species(species_index::Int64; enforce=true)
    set = get_species_set(;enforce=enforce)
    if set === nothing 
        return species_index
    else
        @assert species_index ∈ keys(set.dic_species) "cannot find species with index $species_index in \n $(set.dic_species)"
        return set.dic_species[species_index]
    end
end
get_electron_species() = get_electron()




