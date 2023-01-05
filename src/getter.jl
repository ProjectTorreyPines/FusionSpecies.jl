get_type_species(z::Int64) = get_type_species(float(z))
function get_type_species(z::Float64)
    if z>0
        return Ion
    elseif z == 0
        return Neutral
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

function get_element_type(e::BaseElement)
    return typeof(e).parameters[1]
end

function get_list_elements()
end

function get_next_species_index()
    if length([k for (k,v) in active_species_registry if typeof(v) <: BaseActiveSpecies])>0
        return maximum([k for (k,v) in active_species_registry if typeof(v) <: BaseActiveSpecies])+1
    else
        return 1
    end
end

function get_species_index(s::BaseActiveSpecies)
    return s.index
end

function get_species_index(s::Vector{<:BaseActiveSpecies})
    return [ss.index for ss in s] 
end

function get_nspecies(active_species_set::ActiveSpeciesSet)
    return length(active_species_set.dic_species)
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

function get_active_species()
    return get_active_species_set().list_species
end

function get_mass_active_species()
    return [s.mass for s in get_active_species_set().list_species]
end

function get_Z_active_species()
    return [s.charge_state for s in get_active_species_set().list_species]
end
get_ions() = [s for s in get_active_species_set().list_species if s isa ActiveSpecies{Ion}]

get_ions_index() = [s.index for s in get_ions()]
get_ions_Z() = [s.charge_state for s in get_ions()]
get_species_Z() = [s.charge_state for s in get_active_species_set().list_species]



get_index_active_species(active_species::Vector{Bool}) = Vector{Int64}([i for (i,as) in enumerate(active_species) if as])
get_index_active_species(active_species::Vector{Int64}) = active_species
get_index_active_species(active_species::Vector{<:BaseActiveSpecies}) = Vector{Int64}([s.index for s in active_species])
get_index_active_species() = get_index_active_species(get_active_species())
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
get_active_species(element::Element) = get_species(element::Element)
get_species_except(species::ActiveSpecies) = filter(x->x != species, get_active_species_set().list_species)
get_species(element::Element) = filter(x->x.element == element  ,get_active_species_set().list_species)
get_species(element::Element, Z::Int64) = get_species(get_species_symbol(element.symbol,Z))
function get_active_species(s::Vector{T}) where T<:Union{BaseActiveSpecies,Symbol}
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
function get_species(is_species_active::Vector{Bool})
    return [get_species(i) for (i,b) in enumerate(is_species_active) if b]
end
get_active_species(species_index::Vector{Int64}) = get_species(species_index)
get_species(species_index::Vector{Int64}) = [get_species(s) for s in species_index]
get_active_species(species_index::Int64) = get_species(species_index)
function get_species(species_index::Int64)
    set = get_active_species_set()
    @assert species_index ∈ keys(set.dic_species) 
    return set.dic_species[species_index]
end
get_electron_species() = get_electron()

function get_electron_index(;enforce=true)
    e = get_electron(;enforce=enforce)
    if e isa ActiveSpecies{<:Electron}
        return e.index
    elseif e isa Vector{<:ActiveSpecies{<:Electron}}
        return e 
    elseif e isa Nothing
        return 0
    end
end

function get_electron(; enforce=true)
    ss = [s for s in get_active_species() if s isa ActiveSpecies{<:Electron}]
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


