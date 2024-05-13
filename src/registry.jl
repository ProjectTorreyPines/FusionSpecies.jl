const ElementRegistry = Dict{Symbol,AbstractElement}
const BaseSpeciesRegistry = Dict{Int64,AbstractSpecies}
const LoadedSpeciesRegistry = Dict{Union{Int64,String},Any}
const element_species_registry = Dict{AbstractElement,Vector{AbstractSpecies}}()

const species_registry = LoadedSpeciesRegistry("locked" => false)
const element_registry = ElementRegistry()
const base_species_registry = BaseSpeciesRegistry()

function add2registry(e::Element)
    @assert e ∉ values(element_registry)
    element_registry[Symbol(e.symbol)] = e
end

function add2registry(s::AbstractSpecies)
    @assert s ∉ values(species_registry)
    species_registry[length(species_registry)+1] = s
end
#function add2registry(s::getfield(@__MODULE__,:AbstractLoadedSpecies); registry=nothing)
function add2registry(s::AbstractLoadedSpecies; registry=nothing)
    if registry isa Nothing
        registry = species_registry
    end
    @assert s ∉ values(registry)
    registry[length([v for v in values(registry) if typeof(v) <: AbstractLoadedSpecies])+1] = s
end
# function add2registry(s::AbstractLoadedSpecies; registry=nothing)
#     if registry isa Nothing
#         registry = species_registry
#     end
#     @assert s ∉ values(species_registry)
#     species_registry[length(species_registry)+1] = s
# end

function add2registry(e::AbstractElement, s::AbstractSpecies; registry=nothing)
    if registry isa Nothing
        registry = element_species_registry
    end
    if e ∉ keys(registry)
        registry[e] = Vector{AbstractSpecies}()
    end

    if s ∉ registry[e]
        push!(registry[e], s)
    end
end

function clean_registry(registry)
    for r in keys(registry)
        delete!(registry, r)
    end
end

clean_registry(element_registry)

function check_status_species_registry(; lock=false, message="")
    @assert species_registry["locked"] == lock message * " | species_registry : $species_registry"
end

check_status(species_set::SpeciesSet; lock=false, message="") = @assert species_set.lock[1] == lock message * " | species_registry : $(species_set.lock[1])"


# macro reset_species()
#     for k in keys(species_registry)
#         delete!(species_registry, k)
#     end
#     species_registry["locked"] = false
# end

# function get_species_registry(; lock=true)
#     check_status_species_registry(; lock=lock, message="run first @setup_model")
#     return species_registry
# end