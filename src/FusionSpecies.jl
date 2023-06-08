module FusionSpecies
using PeriodicTable
using Logging
using Crayons.Box
using DocStringExtensions
using Unitful
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
include("elements.jl")
include("types.jl")
include("registry.jl")
include("getter.jl")
include("iterators.jl")


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
        expr = :($n = Element(elements[Symbol($name)]))
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
            for z in 0:el.atomic_number
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

function add_species(obj::Element)
    check_status_species_registry()
    for s in element_species_registry[obj]
        tmp = LoadedSpecies(s, get_next_species_index())
        @assert tmp ∉ collect(values(species_registry))
        add2registry(tmp)
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

macro setup_species()
    expr = quote
        setup_species()
    end

    for s in [v for v in values(species_registry) if v isa AbstractLoadedSpecies]
        ss = s.symbol
        sss = string(s.symbol)
        push!(expr.args, :($ss = get_species(Symbol($sss))))
    end
    list_elements = []
    for s in [v for v in values(species_registry) if v isa AbstractLoadedSpecies]
        ss = get_element(s).symbol
        sss = string(ss)
        push!(list_elements, (ss, sss))
    end
    unique!(list_elements)
    for (n, s) in list_elements
        push!(expr.args, :($n = get_element(Symbol($s))))
    end
    esc(expr)
end



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

function Base.show(io::IO, ::MIME"text/plain", species::AbstractLoadedSpecies)
    print(io, MAGENTA_FG("$(string(species.symbol))"), " [$(stype(species))][", LIGHT_MAGENTA_FG("$(string(species.element.symbol))"), "] ", " - index: $(species.index)")
end

function Base.show(io::IO, species::AbstractLoadedSpecies)
    print(io, MAGENTA_FG("$(string(species.symbol))"), " [$(stype(species))][$(species.index)]")
end

function Base.show(io::IO, ::MIME"text/plain", element::AbstractElement)
    print(io, BLUE_FG("$(string(element.symbol))"), BLUE_FG("Z = $(string(element.atomic_number))"))
end

function Base.show(io::IO, element::AbstractElement)
    print(io, BLUE_FG("$(string(element.symbol))"), BLUE_FG(" $(element.name) - Z = $(string(element.atomic_number))"))
end


function Base.show(io::IO, ::MIME"text/plain", species::SpeciesSet)
    t = Tree(sort(species.dic_species), title="set of species",
        title_style="magenta",
        guides_style="yellow")
    print(io, t)
end

function Base.show(io::IO, species::SpeciesSet)
    t = Tree(sort(species.dic_species), title="set of species",
        title_style="magenta",
        guides_style="yellow")
    print(io, t)
end

function Base.show(io::IO, element_registry::ElementRegistry)
    t = Tree(element_registry, title="Available elements",
        title_style="blue",
        guides_style="blue")
    print(io, t)
end

"$TYPEDSIGNATURES display available elements"
show_elements() = show(element_registry)


Base.:!(s::LoadedSpecies) = get_species_except(s)

Base.to_index(s::LoadedSpecies) = Base.to_index(s.index)

export show_elements, get_element, add_species, show_species, create_species
export @reset_species, @add_plasma_species
export import_species
export AbstractSpecies, BaseSpecies
export base_species_registry
export species_registry
export get_species, get_species, get_species_set, get_electron_species, setup_species
end