module FusionSpecies
using PeriodicTable
using Unitful
using Logging
import Term.Trees: Tree
using Crayons.Box

export import_species
export BaseSpecies
#export @add_species
export species_registry
#export ActiveSpecies
##export BaseActiveSpecies
export active_species_registry
export get_species, @setup_species, get_active_species, get_active_species_set
#import PhysicalConstants.CODATA2018: m_e
m_e = 9.1093837015e−31
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

include("types.jl")
include("registry.jl")
include("getter.jl")


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
    kwargs= Dict(aakws)
    return aargs, kwargs
end

macro create_elements(names...)
    names_ = [string(name) for name in names]
     blk = Expr(:block)
     for name in names_
         n = Symbol(name)
         expr = :($n = Element(elements[Symbol($name)]))
         push!(blk.args,expr)
         expr = :(add2registry($n))
         push!(blk.args,expr)
     end
     esc(blk)
end

macro create_element(symb, args...)
    aargs, kwargs = convert_macro_kwargs(args)
    if :type ∉ keys(kwargs)
        kwargs[:type] = Atom
    end
    n = Symbol(string(symb))
    blk = Expr(:block)
    expr = :($n = Element($(string(kwargs[:name])), $(string(n)),$(kwargs[:atomic_number]), $(kwargs[:mass]); type =  $(kwargs[:type])))
    push!(blk.args,expr)
    expr = :(add2registry($n))
    push!(blk.args,expr)
    esc(blk)
end

macro create_species()
    blk = Expr(:block)
    for el in values(element_registry)

        if get_element_type(el) == Atom
            for z in 0:el.atomic_number
                n = get_species_symbol(el.symbol,z)
                #println("Adding $n ------")
                expr = :($n = Species($el,$z))
                push!(blk.args,expr)
                expr = :(add2registry($n))
                expr = :(add2registry($el,$n))
                push!(blk.args,expr)
            end
        elseif get_element_type(el) == Electron
                z = el.atomic_number
                n = get_species_symbol(el.symbol,z)
                #println("Adding $n ------")
                expr = :($n = Species($el,$z))
                push!(blk.args,expr)
                expr = :(add2registry($n))
                expr = :(add2registry($el,$n))
                push!(blk.args,expr)
        end    
    end
esc(blk)
end

import Base:(==)

function Base.:(==)(s1::BaseActiveSpecies,s2::BaseActiveSpecies)
    #println("COMPARE", s1,s2)
       return all([getfield(s1,f) == getfield(s2,f) for f in fieldnames(typeof(s1)) if f != :index])
       end

macro add_species(objs...)
    blk = Expr(:block)
    for obj in objs
        name = Symbol(obj)
        obj = getfield(@__MODULE__,Symbol(obj))
        expr = :(check_status_active_species_registry())
        push!(blk.args,expr)

        if obj isa Species
            expr = quote
                tmp = ActiveSpecies($obj,get_next_species_index())
                @assert tmp ∉ collect(values(active_species_registry))
                $name = tmp
            end
            push!(blk.args,expr)
            expr = :(add2registry($name))
            push!(blk.args,expr)

        elseif obj isa Element
            for s in element_species_registry[obj]
                name = Symbol(s.symbol)
                expr = quote
                    tmp = ActiveSpecies($s,get_next_species_index())
                    @assert tmp ∉ collect(values(active_species_registry))
                    $name = tmp
                end
                push!(blk.args,expr)
                expr = :(add2registry($name))
                push!(blk.args,expr)
            end

        end
    end
    esc(blk)
end
# end


@create_elements H Be C Si W Ne Ar Xe 
@create_element ∅ mass=0 name=dummy atomic_number=0
@create_element D mass=2*H.mass name=deuterium atomic_number=1
@create_element T mass=3*H.mass name=tritium atomic_number=1
#@create_element e mass=m_e.val name=electron atomic_number=-1 type=Electron
@create_element e mass=m_e name=electron atomic_number=-1 type=Electron
@create_species 

macro setup_species()
    expr = quote
        setup_species()
    end

    for s in [v for v in values(active_species_registry) if v isa BaseActiveSpecies]
        ss = Symbol(s.symbol)
        

    push!(expr.args,:($ss = get_active_species($ss)))
    end

esc(expr)    
end

function import_species()
    for L in names(@__MODULE__,all=true)
        o = getfield(@__MODULE__,L)
        if o isa BaseSpecies
            getfield(@__MODULE__,:species_registry)[length(species_registry)+1] = o
            add2registry(o.element,o;registry=getfield(@__MODULE__,:element_species_registry))
        elseif o isa BaseElement
            getfield(@__MODULE__,:element_registry)[length(species_registry)+1] = o

        end
    end
end
import_species()



function setup_species()
    list_species = [v for (k,v) in active_species_registry if typeof(v) <: BaseActiveSpecies]
    @assert length(unique( list_species)) == length(list_species)
    list_idx = [v.index for (k,v) in active_species_registry if typeof(v) <: BaseActiveSpecies]
    @assert length(unique( list_idx)) == length(list_idx)
    active_species_registry["active_species_set"] = ActiveSpeciesSet(active_species_registry)
    active_species_registry["locked"] = true
    show(active_species_registry["active_species_set"])
end



function check_active_species_index(active_species_set::ActiveSpeciesSet)
    for (k,v) in active_species_set.dic_species
        @assert k == v.index
    end
    indexes = sort([v.index for (k,v) in active_species_set.dic_species])
    # check that indexes start at 1
    @assert minimum(indexes) == 1
    # check that indexes are incremental by 1
    @assert minimum(diff(indexes)) == maximum(diff(indexes)) == 1
end

function check_active_species_index(active_species_registry::ActiveSpeciesRegistry)
    indexes = sort([v.index for v in values(active_species_registry) if typeof(v) <: BaseActiveSpecies])
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




function Base.show(io::IO, ::MIME"text/plain", species::BaseActiveSpecies)
    print(io, MAGENTA_FG("$(string(species.symbol))")," [",LIGHT_MAGENTA_FG("$(string(species.element.symbol))"), "] ", " - index: $(species.index)")
end

function Base.show(io::IO, species::BaseActiveSpecies)
    print(io, MAGENTA_FG("$(string(species.symbol))")," [",LIGHT_MAGENTA_FG("$(string(species.element.symbol))"), "] ", " - index: $(species.index)")
end

function Base.show(io::IO, ::MIME"text/plain", species::ActiveSpeciesSet)
    t = Tree(species.dic_species ,title="Active Species",
    title_style="magenta",
    guides_style="yellow")
    print(io, t)
end

function Base.show(io::IO, species::ActiveSpeciesSet)
    t = Tree(species.dic_species ,title="Active Species",
    title_style="magenta",
    guides_style="yellow")
    print(io, t)
end


end