get_type_species(z::Int64) = get_type_species(float(z))

function get_species_type(s::Species)
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
function get_element(s::ActiveSpecies)
    return s.element
end

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
    if length([k for (k,v) in active_species_registry if typeof(v) <: AbstractActiveSpecies])>0
        return maximum([k for (k,v) in active_species_registry if typeof(v) <: AbstractActiveSpecies])+1
    else
        return 1
    end
end



function get_nspecies(active_species_set::ActiveSpeciesSet)
    return length(active_species_set.dic_species)
end

function show_active_species()
    check_status_active_species_registry(lock=true)
    show(active_species_registry["active_species_set"])
end

function get_active_species_set()
    check_status_active_species_registry(lock=true)
    return active_species_registry["active_species_set"]
end

function get_nspecies()
    @assert "active_species_set" ∈ keys(active_species_registry)
    check_status_active_species_registry(lock=true)
    return get_nspecies(active_species_registry["active_species_set"])
end

" $TYPEDSIGNATURES return the electron species in active species"
function get_electron(; enforce=true)
    ss = [s for s in get_active_species() if type(s) <:Electron]
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
    ss = [s for s in get_active_species() if type(s) <: MainIon]
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
    ss = [s for s in get_active_species() if type(s) <:MainAtom]
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
get_active_species() ::Vector{<:AbstractActiveSpecies} = get_active_species_set().list_species

"$TYPEDSIGNATURES get a list of the mass of active species"
get_active_species_mass() ::Vector{Float64} = [s.mass for s in get_active_species_set().list_species]


"$TYPEDSIGNATURES get a list of the charge state of active species"
get_active_species_Z() = [s.charge_state for s in get_active_species_set().list_species]

"$TYPEDSIGNATURES get a list of ions among active species"
get_ions() ::Vector{<:ActiveSpecies}  = [s for s in get_active_species_set().list_species if type(s) <: Ions]

"$TYPEDSIGNATURES get a list of ions among active species"
get_imp_ions() ::Vector{<:AbstractActiveSpecies}  = [s for s in get_active_species_set().list_species if type(s) <: ImpurityIon]

"$TYPEDSIGNATURES get a list of ions among active species"
get_imp_atoms() ::Vector{<:AbstractActiveSpecies}  = [s for s in get_active_species_set().list_species if type(s) <: ImpurityAtom]

"$TYPEDSIGNATURES get a list of atoms among active species"
get_atoms() ::Vector{<:AbstractActiveSpecies} = [s for s in get_active_species_set().list_species if type(s) <: Atoms]

"$TYPEDSIGNATURES get a list of molecules among active species"
get_molecules() ::Vector{<:AbstractActiveSpecies} = [s for s in get_active_species_set().list_species if type(s) <: Molecules]


"$TYPEDSIGNATURES get a list of atoms among active species"
get_neutrals() = [s for s in get_active_species_set().list_species if type(s) <: Neutrals]

"$TYPEDSIGNATURES get a list of index of neutrals among active species"
get_neutrals_index() :: Vector{Int64} = get_active_species_index(get_neutrals())

"$TYPEDSIGNATURES get a list of index of ions among active species"
get_ions_index() :: Vector{Int64} = get_active_species_index(get_ions())

"$TYPEDSIGNATURES get a list of index of ions among active species"
get_main_ion_index(;kw...) :: Int64 = get_active_species_index(get_main_ion(;kw...))

"$TYPEDSIGNATURES get the index of main atom among active species"
get_main_atom_index(;kw...) :: Int64 = get_active_species_index(get_main_atom(;kw...))

"$TYPEDSIGNATURES get a list of index of ions among active species"
get_imp_ions_index() :: Vector{Int64} = get_active_species_index(get_imp_ions())

"$TYPEDSIGNATURES get a list of molecules among active species"
get_molecules_index() ::Vector{Int64} = get_active_species_index(get_molecules())

"$TYPEDSIGNATURES get a list of index of impurity atoms among active species"
get_imp_atoms_index() :: Vector{Int64} = get_active_species_index(get_imp_atoms())

"$TYPEDSIGNATURES get a list of index of atoms among active species"
get_atoms_index() :: Vector{Int64} = get_active_species_index(get_atoms())


" $TYPEDSIGNATURES return the index of the electron species"
function get_electron_index(;enforce=true)
    e = get_electron(;enforce=enforce)
    if e isa ActiveSpecies{<:Electron}
        return e.index
    elseif e isa Nothing
        error("no electron found")
    end
end


get_ions_Z() = [s.charge_state for s in get_ions()]
get_species_Z() = [s.charge_state for s in get_active_species_set().list_species]


"""$TYPEDSIGNATURES get a list of the index of active species
   $(METHODLIST)
"""
get_active_species_index(s::AbstractActiveSpecies) = s.index
get_active_species_index(active_species::Vector{Bool}) = Vector{Int64}([i for (i,as) in enumerate(active_species) if as])
get_active_species_index(active_species::Vector{Int64}) = active_species
get_active_species_index(active_species::Vector{<:AbstractActiveSpecies}) = Vector{Int64}([s.index for s in active_species])
get_active_species_index() ::Vector{Int64} = get_active_species_index(get_active_species())

" $(TYPEDSIGNATURES) return the species that are of type `element`"
function get_species(symbol::Symbol; lock=true)
    @assert "active_species_set" ∈ keys(active_species_registry)
    check_status_active_species_registry(lock=lock)
    species =  filter(x->Symbol(x.symbol) == symbol  ,active_species_registry["active_species_set"].list_species)
    @assert length(species) <2 
    if length(species) == 1
        return species[1]
    elseif length(species) == 0
        return nothing
    end

end
" $(TYPEDSIGNATURES) return the species that are of type `element`"
get_active_species(element::Element) = get_species(element::Element)
get_species_except(species::ActiveSpecies) = filter(x->x != species, get_active_species_set().list_species)
get_species(element::Element) = filter(x->x.element == element  ,get_active_species_set().list_species)
get_species(element::Element, Z::Int64) = get_species(get_species_symbol(element.symbol,Z))
function get_active_species(s::Vector{T}) where T<:Union{AbstractActiveSpecies,Symbol}
    out = Vector{ActiveSpecies}()
    for ss in s 
        append!(out,get_active_species(ss))
    end
    return out
end
get_active_species(s::ActiveSpecies) = [s]
get_species(species::ActiveSpecies) = [species]
name(s::ActiveSpecies) = string(s.symbol)

get_active_species_name(s) = name.(get_active_species(s))
get_species(is_species_active::Vector{Bool}) =  [get_species(i) for (i,b) in enumerate(is_species_active) if b]
get_active_species(species_index::Vector{Int64}) = get_species(species_index)
get_species(species_index::Vector{Int64}) = [get_species(s) for s in species_index]
get_active_species(species_index::Int64) = get_species(species_index)

" $TYPEDSIGNATURES return the species with index `species_index`"
function get_species(species_index::Int64)
    set = get_active_species_set()
    @assert species_index ∈ keys(set.dic_species) 
    return set.dic_species[species_index]
end
get_electron_species() = get_electron()




