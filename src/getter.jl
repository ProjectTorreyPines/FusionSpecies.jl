function get_type_species(z::Int64)
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
function get_species_symbol(e_name,z)
    return Symbol(e_name * get_str_charge_state(z)) 
end

function get_element_type(e::BaseElement)
    return typeof(e).parameters[1]
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




function get_species(symbol::Symbol)
    @assert "active_species_set" ∈ keys(active_species_registry)
    check_status_active_species_registry(lock=true)
    species =  filter(x->Symbol(x.symbol) == symbol  ,active_species_registry["active_species_set"].list_species)
    @assert length(species) <2 
    if length(species) == 1
        return species[1]
    elseif length(species) == 0
        return nothing
    end

end

function get_species(element::Element)
    @assert "active_species_set" ∈ keys(active_species_registry)
    check_status_active_species_registry(lock=true)
    species =  filter(x->x.element == element  ,active_species_registry["active_species_set"].list_species)
    return species
end

function get_species(species::ActiveSpecies)
    return species
end

function get_species(is_species_active::Vector{Bool})
    return [get_species(i) for (i,b) in enumerate(is_species_active) if b]
end

function get_species(species_index::Int64)
    set = get_active_species_set()
    @assert species_index ∈ keys(set.dic_species) 
    return set.dic_species[species_index]
end

function get_electron_index()
    check_status_active_species_registry(lock=true)
    for s in values(active_species_registry)
        if typeof(s) <: ActiveSpecies{<:Electron}
            return s.index
        end
    end
    return 0
end