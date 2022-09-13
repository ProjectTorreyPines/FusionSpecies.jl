const ElementRegistry = Dict{Int64,BaseElement}
const SpeciesRegistry = Dict{Int64,BaseSpecies}
const ActiveSpeciesRegistry = Dict{Union{Int64,String},Any}
const element_species_registry = Dict{BaseElement,Vector{BaseSpecies}}()

const active_species_registry = ActiveSpeciesRegistry("locked"=>false)
const element_registry = ElementRegistry()
const species_registry = SpeciesRegistry()

function add2registry(e::Element)
    @assert e ∉ values(element_registry)
    element_registry[length(element_registry)+1] = e
end

function add2registry(s::BaseSpecies)
    @assert s ∉ values(species_registry)
    species_registry[length(species_registry)+1] = s
end
#function add2registry(s::getfield(@__MODULE__,:BaseActiveSpecies); registry=nothing)
    function add2registry(s::BaseActiveSpecies; registry=nothing)
    if registry isa Nothing
        registry = active_species_registry
    end
    @assert s ∉ values(registry)
    registry[length([v for v in values(registry) if typeof(v) <: BaseActiveSpecies])+1] = s
end
# function add2registry(s::BaseActiveSpecies; registry=nothing)
#     if registry isa Nothing
#         registry = active_species_registry
#     end
#     @assert s ∉ values(active_species_registry)
#     active_species_registry[length(active_species_registry)+1] = s
# end

function add2registry(e::BaseElement,s::BaseSpecies; registry=nothing)
    if registry isa Nothing
        registry = element_species_registry
    end
    if e ∉ keys(registry)
        registry[e] = Vector{BaseSpecies}()
    end

    if  s ∉ registry[e]
        push!(registry[e],s) 
    end
end

function clean_registry(registry)
    for r in keys(registry)
    delete!(registry,r)
    end
end

clean_registry(element_registry)

function check_status_active_species_registry(;lock=false)
    @assert active_species_registry["locked"] == lock "active_species_registry : $active_species_registry"
end

macro reset_species()
    for k in keys(active_species_registry)
        delete!(active_species_registry,k)
    end
    active_species_registry["locked"] = false
    end