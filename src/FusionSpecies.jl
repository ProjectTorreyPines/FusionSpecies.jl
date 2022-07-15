module FusionSpecies
using RedefStructs
using PeriodicTable
using Unitful
import PhysicalConstants.CODATA2018: m_e

abstract type BaseSpecies end
abstract type BaseActiveSpecies end
abstract type BaseElement end
abstract type Atom end
abstract type Dummy end
abstract type Electron end
abstract type Molecule end
abstract type Neutral <: Atom end
abstract type Ion <: Atom end


const ElementRegistry = Dict{Int64,BaseElement}
const SpeciesRegistry = Dict{Int64,BaseSpecies}
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

function get_element_type(e::BaseElement)
    return typeof(e).parameters[1]
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



@redef struct Element{T} <: BaseElement 
    name :: String
    symbol :: String
    atomic_number :: Int64
    mass :: Float64 
    function Element(element::PeriodicTable.Element)
        new{Atom}(element.name,element.symbol , element.number, (element.atomic_mass |> u"kg").val)
    end
    function Element(n, s ,a ,m; type= Atom)
        new{type}(n,s,a,m)
    end
end
function get_type_species(z::Int64)
    if z>0
        return Ion
    elseif z == 0
        return Neutral
    else
        return Electron
    end
end

@redef struct Species{T} <:BaseSpecies
    charge_state :: Int64
    name :: String
    symbol :: String
    element :: BaseElement
    mass :: Float64
    atomic_number:: Int64
    function Species(e::BaseElement,z::Int64)
        T = get_type_species(z) 
        s = string(get_species_symbol(e.symbol,z))
        new{T}(z,e.name,s,e,e.mass,e.atomic_number)
    end
end

struct ActiveSpecies{T} <: BaseActiveSpecies
    charge_state :: Int64
    name :: String
    symbol :: String
    element :: BaseElement
    mass :: Float64
    atomic_number:: Int64
    index:: Int64

    function ActiveSpecies(s::Species, idx::Int64)
        T = get_type_species(s.charge_state) 
        new{T}(s.charge_state,s.name,s.symbol,s.element,s.mass,s.atomic_number,idx)
    end
    
    # function ActiveSpecies(s::Nothing)
    #     new{Dummy}(0,"nothing","nothing",Ø,0.0,0,0)
    # end
end

# import Base:(==)

# function Base.:(==)(s1::BaseActiveSpecies,s2::BaseActiveSpecies)
#        return all([getfield(s1,f) == getfield(s2,f) for f in fieldnames(typeof(s1)) if f != :index && f != :rates])
#        end

function get_species_index(s::BaseActiveSpecies)
    return s.index
end

@redef struct ActiveSpeciesSet
    list_species :: Vector{BaseActiveSpecies}
    dic_species :: Dict{Int64,BaseActiveSpecies}
    function ActiveSpeciesSet(active_species_registry)
        # check id species
        check_active_species_index(active_species_registry)
        
        # make dict with species indexes
        dic_species = Dict{Int64,BaseActiveSpecies}()
        for s in [v for v in values(active_species_registry) if typeof(v) <: BaseActiveSpecies]
            dic_species[s.index] = s 
        end

        new([v for v in values(active_species_registry) if typeof(v) <: BaseActiveSpecies],dic_species)
    end
end



const ActiveSpeciesRegistry = Dict{Union{Int64,String},Any}
const active_species_registry = ActiveSpeciesRegistry()
const element_registry = ElementRegistry()
const species_registry = SpeciesRegistry()
active_species_registry["locked"] = false
const element_species_registry = Dict{BaseElement,Vector{BaseSpecies}}()
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
    #print(blk)
esc(blk)
end

function get_next_species_index()
    if length([k for (k,v) in active_species_registry if typeof(v) <: BaseActiveSpecies])>0
        return maximum([k for (k,v) in active_species_registry if typeof(v) <: BaseActiveSpecies])+1
    else
        return 1
    end
end

function check_status_active_species_registry(;lock=false)
    @assert active_species_registry["locked"] == lock
end

import Base:(==)

function Base.:(==)(s1::BaseActiveSpecies,s2::BaseActiveSpecies)
    #println("COMPARE", s1,s2)
       return all([getfield(s1,f) == getfield(s2,f) for f in fieldnames(typeof(s1)) if f != :index])
       end

macro add_species(obj)
    name = Symbol(obj)
    obj = getfield(@__MODULE__,Symbol(obj))
    blk = Expr(:block)
    expr = :(check_status_active_species_registry())
    push!(blk.args,expr)
    if obj isa Species
        expr = quote
            tmp = ActiveSpecies($obj,get_next_species_index())
            @assert tmp ∉ collect(values(active_species_registry))
            $name = tmp
        end
        # expr = :(println(typeof(ActiveSpecies($obj))))
        # expr = :(println(typeof(ActiveSpecies($obj)) <: BaseActiveSpecies))
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
    esc(blk)
end
# end


@create_elements H Be C Si W Ne Ar Xe 
@create_element ∅ mass=0 name=dummy atomic_number=0
@create_element D mass=2*H.mass name=deuterium atomic_number=1
@create_element T mass=3*H.mass name=tritium atomic_number=1
@create_element e mass=m_e.val name=electron atomic_number=-1 type=Electron
@create_species 

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
end

macro reset_species()
for k in keys(active_species_registry)
    delete!(active_species_registry,k)
end
active_species_registry["locked"] = false
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
    # check that indexes start at 1
    @assert minimum(indexes) == 1
    # check that indexes are incremental by 1
    @assert minimum(diff(indexes)) == maximum(diff(indexes)) == 1
end



function get_nspecies(active_species_set::ActiveSpeciesSet)
    return length(active_species_set.dic_species)
end

function get_nspecies()
    @assert "active_species_set" ∈ keys(active_species_registry)
    check_status_active_species_registry(lock=true)
    return get_nspecies(active_species_registry["active_species_set"])
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

export import_species
export BaseSpecies
#export @add_species
export species_registry
#export ActiveSpecies
##export BaseActiveSpecies
export active_species_registry
#export add2registry
end