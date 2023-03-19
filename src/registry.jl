const ElementRegistry = Dict{Symbol,AbstractElement}
const SpeciesRegistry = Dict{Int64,AbstractSpecies}
const ActiveSpeciesRegistry = Dict{Union{Int64,String},Any}
const element_species_registry = Dict{AbstractElement,Vector{AbstractSpecies}}()

const active_species_registry = ActiveSpeciesRegistry("locked"=>false)
const element_registry = ElementRegistry()
const species_registry = SpeciesRegistry()

function add2registry(e::Element)
    @assert e ∉ values(element_registry)
    element_registry[Symbol(e.symbol)] = e
end

function add2registry(s::AbstractSpecies)
    @assert s ∉ values(species_registry)
    species_registry[length(species_registry)+1] = s
end
#function add2registry(s::getfield(@__MODULE__,:AbstractActiveSpecies); registry=nothing)
    function add2registry(s::AbstractActiveSpecies; registry=nothing)
    if registry isa Nothing
        registry = active_species_registry
    end
    @assert s ∉ values(registry)
    registry[length([v for v in values(registry) if typeof(v) <: AbstractActiveSpecies])+1] = s
end
# function add2registry(s::AbstractActiveSpecies; registry=nothing)
#     if registry isa Nothing
#         registry = active_species_registry
#     end
#     @assert s ∉ values(active_species_registry)
#     active_species_registry[length(active_species_registry)+1] = s
# end

function add2registry(e::AbstractElement,s::AbstractSpecies; registry=nothing)
    if registry isa Nothing
        registry = element_species_registry
    end
    if e ∉ keys(registry)
        registry[e] = Vector{AbstractSpecies}()
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

function check_status_active_species_registry(;lock=false, message ="" )
    @assert active_species_registry["locked"] == lock message * " | active_species_registry : $active_species_registry"
end

macro reset_species()
    for k in keys(active_species_registry)
        delete!(active_species_registry,k)
    end
    active_species_registry["locked"] = false
    end

    function get_active_species_registry(;lock=true) 
        check_status_active_species_registry(;lock=lock, message="run first @setup_variables")
        return active_species_registry 
    end