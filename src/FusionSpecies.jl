module FusionSpecies

import PeriodicTable
using Logging
using Crayons.Box
using DocStringExtensions
using Unitful
using PrecompileTools
using MacroTools
using Adapt


m_e = 9.1093837015e-31
dic_expo = Dict{}()
dic_expo[1] = "¹"
dic_expo[2] = "²"
dic_expo[3] = "³"
dic_expo[4] = "⁴"
dic_expo[5] = "⁵"
dic_expo[6] = "⁶"
dic_expo[7] = "⁷"
dic_expo[8] = "⁸"
dic_expo[9] = "⁹"
dic_expo[0] = "⁰"
dic_expo["+"] = "⁺"
dic_expo["-"] = "⁻"
include("gpu_generator.jl")
include("properties.jl")
include("types.jl")
include("elements.jl")
include("registry.jl")
include("getter.jl")
include("iterators.jl")

dummy_loaded_species = LoadedSpecies{DummyParticles}(SpeciesChargeState(0.0), "dummy", :dummy, dummy_element, SpeciesMass(0.0), SpeciesAtomicNumber(0), SpeciesIndex(0), [false])
function convert_macro_kwargs(args)
    aargs = []
    aakws = Pair{Symbol,Any}[]
    for el in args
        if Meta.isexpr(el, :(=))
            push!(aakws, Pair(el.args...))
        else
            push!(aargs, el)
        end
    end
    kwargs = Dict(aakws)
    return aargs, kwargs
end

macro create_elements(names...)
    names_ = [string(name) for name in names]
    blk = Expr(:block)
    for name in names_
        n = Symbol(name)
        expr = :($n = Element(PeriodicTable.elements[Symbol($name)]))
        push!(blk.args, expr)
        expr = :(add2registry($n))
        push!(blk.args, expr)
    end
    esc(blk)
end

macro create_element(symb, args...)
    aargs, kwargs = convert_macro_kwargs(args)
    if :type ∉ keys(kwargs)
        kwargs[:type] = Atoms
    end
    n = Symbol(string(symb))
    s = string(symb)
    blk = Expr(:block)
    expr = :($n = Element($(string(kwargs[:name])), Symbol($s), $(kwargs[:atomic_number]), $(kwargs[:mass]), $(kwargs[:density]); type=$(kwargs[:type])))
    push!(blk.args, expr)
    expr = :(add2registry($n))
    push!(blk.args, expr)
    esc(blk)
end

macro create_species()
    blk = Expr(:block)
    for el in values(element_registry)

        if get_element_type(el) == Atoms
            for z in 0:el.atomic_number.value
                n = get_species_symbol(el.symbol, z)
                #println("Adding $n ------")
                expr = :($n = BaseSpecies($el, $z))
                push!(blk.args, expr)
                expr = :(add2registry($n))
                expr = :(add2registry($el, $n))
                push!(blk.args, expr)
            end
        elseif get_element_type(el) == Electron
            z = el.atomic_number
            n = get_species_symbol(el.symbol, z)
            #println("Adding $n ------")
            expr = :($n = BaseSpecies($el, $z))
            push!(blk.args, expr)
            expr = :(add2registry($n))
            expr = :(add2registry($el, $n))
            push!(blk.args, expr)
        end
    end
    esc(blk)
end
macro add_plasma_species()
    create_species()
end
function create_species()
    blk = Expr(:block)
    for (k, el) in element_registry

        if get_element_type(el) == Atoms
            for z in 0:el.atomic_number
                n = get_species_symbol(el.symbol, z)
                es = el.symbol
                #println("Adding $n ------")
                expr = :(const $n = BaseSpecies($el, $z))

                push!(blk.args, expr)
                # expr = :(add2registry($n))
                expr = :(const $es = $el)
                # expr = :(add2registry($el,$n))
                push!(blk.args, expr)
            end
            # elseif get_element_type(el) == Electron
            #         z = el.atomic_number
            #         n = get_species_symbol(el.symbol,z)
            #         #println("Adding $n ------")
            #         expr = :(const $n = Species($el,$z))
            #         push!(blk.args,expr)
            #         # expr = :(add2registry($n))
            #         # expr = :(add2registry($el,$n))
            #         push!(blk.args,expr)
        end
    end
    esc(blk)
end
import Base: (==)

function Base.:(==)(s1::AbstractLoadedSpecies, s2::AbstractLoadedSpecies)
    return all([getfield(s1, f) == getfield(s2, f) for f in fieldnames(typeof(s1)) if f != :index])
end

function add_species(obj::BaseSpecies)
    check_status_species_registry()
    tmp = LoadedSpecies(obj, get_next_species_index())
    @assert tmp ∉ collect(values(species_registry)) "Species $obj already added.... \n List of current species: $(collect(values(species_registry))) \n Use @reset_species to clear the species registry"
    add2registry(tmp)
    return tmp
end

function add_species(obj::BaseSpecies, species_set::SpeciesSet)
    check_status(species_set)
    tmp = LoadedSpecies(obj, get_next_species_index(species_set))
    @assert tmp ∉ species_set.list_species "Species $obj already added.... \n List of current species: $(species_set.list_species) \n Use @reset_species to clear the species registry"
    push!(species_set.list_species, tmp)
end

add_species(obj::Symbol, species_set::SpeciesSet) = add_species(getfield(@__MODULE__, obj), species_set)

function add_species(obj::Element, species_set::SpeciesSet)
    check_status(species_set)
    for s in element_species_registry[obj]
        add_species(s, species_set)
    end
end

macro add_species(objs...)
    for obj in objs
        name = Symbol(obj)
        obj = getfield(@__MODULE__, Symbol(obj))
        add_species(obj)
    end
end


@create_elements H He Be C Si W Ne Ar Xe
@create_element ∅ mass = 0 name = dummy atomic_number = 0 density = 0.0
@create_element D mass = 2 * H.mass name = deuterium atomic_number = 1 density = 0.0
@create_element T mass = 3 * H.mass name = tritium atomic_number = 1 density = 0.0
@create_element e mass = m_e name = electron atomic_number = -1 type = Electron density = 0.0
@create_species




function import_species()
    for L in names(@__MODULE__, all=true)
        o = getfield(@__MODULE__, L)
        if o isa AbstractSpecies
            getfield(@__MODULE__, :species_registry)[length(species_registry)+1] = o
            add2registry(o.element, o; registry=getfield(@__MODULE__, :element_species_registry))
        elseif o isa AbstractElement
            getfield(@__MODULE__, :element_registry)[o.symbol] = o
        end
    end
end
import_species()




function setup_species()
    reorder_species_index(species_registry)
    list_species = [v for (k, v) in species_registry if typeof(v) <: AbstractLoadedSpecies]
    @assert length(unique(list_species)) == length(list_species)

    list_idx = [v.index for (k, v) in species_registry if typeof(v) <: AbstractLoadedSpecies]
    @assert length(unique(list_idx)) == length(list_idx)

    species_registry["species_set"] = LoadedSpeciesSet(species_registry)
    species_registry["locked"] = true
    show(species_registry["species_set"])
end

function setup_species(species_set::SpeciesSet)
    reorder_species_index(species_set)

    list_species = species_set.list_species
    @assert length(unique(list_species)) == length(list_species)

    list_idx = [v.index for v in species_set.list_species]
    @assert length(unique(list_idx)) == length(list_idx)
    species_set.lock[1] = true
end

function reorder_species_index(species_registry)
    # list_species = [v for (k,v) in species_registry if typeof(v) <: AbstractLoadedSpecies]
    # nspc = length(list_species)
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

function check_species_index(species_registry::LoadedSpeciesRegistry)
    indexes = sort([v.index for v in values(species_registry) if typeof(v) <: AbstractLoadedSpecies])
    @debug begin
        "Indexes of active species: $indexes"
    end
    # check that indexes start at 1
    @assert minimum(indexes) == 1
    # check that indexes are incremental by 1
    if length(indexes) > 1
        @assert minimum(diff(indexes)) == maximum(diff(indexes)) == 1
    end
end


Base.ones(sp::FusionSpecies.SpeciesParameters) = ones(get_nspecies(sp.species_set))

function Base.show(io::IO, ::MIME"text/plain", species::LoadedSpecies)
    print(io, MAGENTA_FG("$(string(species.symbol))"), " [$(stype(species))][", LIGHT_MAGENTA_FG("$(string(species.element.symbol))"), "] ", " - index: $(species.index)")
end

function Base.show(io::IO, species::LoadedSpecies)
    print(io, MAGENTA_FG("$(string(species.symbol))"), " [$(stype(species))][$(species.index)]")
end

function Base.show(io::IO, ::MIME"text/plain", element::AbstractElement)
    print(io, BLUE_FG("$(string(element.symbol))"), BLUE_FG("Z = $(string(element.atomic_number))"))
end

function Base.show(io::IO, element::AbstractElement)
    print(io, BLUE_FG("$(string(element.symbol))"), BLUE_FG(" $(element.name) - Z = $(string(element.atomic_number))"))
end


function Base.show(io::IO, ::MIME"text/plain", species_set::SpeciesSet)
    println("species set")
    for s in species_set
        println(s)
    end
end

function Base.show(io::IO, species_set::SpeciesSet{T}) where {T}
    println("species set")
    for s in species_set
        println(s)
    end
end

function Base.show(io::IO, element_registry::ElementRegistry)
    t = Tree(element_registry, title="Available elements",
        title_style="blue",
        guides_style="blue")
    print(io, t)
end

"$TYPEDSIGNATURES display available elements"
show_elements() = show(element_registry)
name_(species::Vector{<:AbstractElement}) = join([stringstyled(string(s.symbol); color=:orange, bold=true) for s in species], " ")
name_(species::Vector{<:AbstractSpecies}) = join([stringstyled(string(s.symbol); color=:magenta, bold=true) for s in species], " ")
name_(species::Vector{Symbol}) = join([stringstyled(s; color=:magenta, bold=true) for s in species], " ")
name_(species_set::SpeciesSet) = join(name_.(species_set.list_species), " ")
name_(species::AbstractSpecies) = stringstyled(string(species.symbol); color=:magenta, bold=true)
name_(species::AbstractElement) = stringstyled(string(species.symbol); color=:orange, bold=true)
inline_summary(species_set::SpeciesSet) = join([name_(species) * "[$(species.index.value)]" for species in species_set.list_species], " ")
Base.:!(s::LoadedSpecies) = get_species_except(s)

Base.to_index(s::LoadedSpecies) = Base.to_index(s.index)
import Base: text_colors, disable_text_style
function stringstyled(str::AbstractString; color::Union{Int,Symbol}=:normal,
    bold::Bool=false, underline::Bool=false, blink::Bool=false,
    reverse::Bool=false, hidden::Bool=false)
    bold && color === :bold && (color = :nothing)
    underline && color === :underline && (color = :nothing)
    blink && color === :blink && (color = :nothing)
    reverse && color === :reverse && (color = :nothing)
    hidden && color === :hidden && (color = :nothing)
    enable_ansi = get(text_colors, color, text_colors[:default]) *
                  (bold ? text_colors[:bold] : "") *
                  (underline ? text_colors[:underline] : "") *
                  (blink ? text_colors[:blink] : "") *
                  (reverse ? text_colors[:reverse] : "") *
                  (hidden ? text_colors[:hidden] : "")

    disable_ansi = (hidden ? disable_text_style[:hidden] : "") *
                   (reverse ? disable_text_style[:reverse] : "") *
                   (blink ? disable_text_style[:blink] : "") *
                   (underline ? disable_text_style[:underline] : "") *
                   (bold ? disable_text_style[:bold] : "") *
                   get(disable_text_style, color, text_colors[:default])

    return enable_ansi * str * disable_ansi
end

include(get_generated_gpu_struct_filename())
const Elements = Union{Element,Vector{<:Element}}
export show_elements, get_element, add_species, show_species, create_species
export @reset_species, @add_plasma_species
export import_species, is_set
export base_species_registry, species_registry
export get_species, get_species, get_species_set, get_electron_species, setup_species, get_species_index, get_electron_index
export name_, check_status_species_registry, species_registry, get_nspecies, get_species_Z, get_species_masses, get_electron_index, get_species_abstract_type
export AbstractSpecies, BaseSpecies, SpeciesSet, LoadedSpeciesSet, AbstractLoadedSpecies, Species, Elements, AbstractElement, LoadedSpecies, SpeciesParameters
export ElectronSpecies, IonSpecies
end