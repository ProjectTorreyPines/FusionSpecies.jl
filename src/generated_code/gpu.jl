abstract type AbstractSpeciesIndexesGPU <: AbstractSpeciesIndexes end
export AbstractSpeciesIndexesGPU
struct SpeciesIndexesGPU{AAA} <: AbstractSpeciesIndexesGPU
    value::AAA
end
#= /Users/jeromeguterl/development/StataMater/FusionSpecies.jl/src/gpu_generator.jl:82 =# @nospecialize
function AbstractSpeciesIndexesGPU(#= /Users/jeromeguterl/development/StataMater/FusionSpecies.jl/src/gpu_generator.jl:92 =# @nospecialize(value); kw...)
    SpeciesIndexesGPU{typeof(value)}(value; kw...)
end
get_base_type_gpu(::SpeciesIndexes) = begin
        SpeciesIndexesGPU
    end
function Adapt.adapt_structure(to, v::AbstractSpeciesIndexes; kw...)
    AbstractSpeciesIndexesGPU((Adapt.adapt_structure(to, getproperty(v, f)) for f = propertynames(v))...)
end
function special_copy(#= /Users/jeromeguterl/development/StataMater/FusionSpecies.jl/src/gpu_generator.jl:113 =# @nospecialize(v::AbstractSpeciesIndexes), #= /Users/jeromeguterl/development/StataMater/FusionSpecies.jl/src/gpu_generator.jl:113 =# @nospecialize(backend); )
    AbstractSpeciesIndexesGPU((special_copy(getproperty(v, f), backend) for f = propertynames(v))...)
end
#= /Users/jeromeguterl/development/StataMater/FusionSpecies.jl/src/gpu_generator.jl:120 =# @specialize
