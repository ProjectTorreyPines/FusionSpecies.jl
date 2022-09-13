module FusionSpecies
using PeriodicTable
using Unitful
using Logging
import Term.Trees: Tree
using Crayons.Box
#import PhysicalConstants.CODATA2018: m_e
m_e = 9.1093837015e−31
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

struct Element{T} <: BaseElement 
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

struct Species{T} <:BaseSpecies
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

function get_species_index(s::Vector{<:BaseActiveSpecies})
    return [ss.index for ss in s] 
end

struct ActiveSpeciesSet
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
    @assert active_species_registry["locked"] == lock "active_species_registry : $active_species_registry"
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

macro setup_species()
    expr = quote
        setup_species()
        # list_species = [v for (k,v) in active_species_registry if typeof(v) <: BaseActiveSpecies]
        # @assert length(unique( list_species)) == length(list_species)
        # list_idx = [v.index for (k,v) in active_species_registry if typeof(v) <: BaseActiveSpecies]
        # @assert length(unique( list_idx)) == length(list_idx)
        # active_species_registry["active_species_set"] = ActiveSpeciesSet(active_species_registry)
        # active_species_registry["locked"] = true
    end
esc(expr)    
end
function setup_species()
    list_species = [v for (k,v) in active_species_registry if typeof(v) <: BaseActiveSpecies]
    @assert length(unique( list_species)) == length(list_species)
    list_idx = [v.index for (k,v) in active_species_registry if typeof(v) <: BaseActiveSpecies]
    @assert length(unique( list_idx)) == length(list_idx)
    active_species_registry["active_species_set"] = ActiveSpeciesSet(active_species_registry)
    active_species_registry["locked"] = true
    show(active_species_registry["active_species_set"])
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

export import_species
export BaseSpecies
#export @add_species
export species_registry
#export ActiveSpecies
##export BaseActiveSpecies
export active_species_registry
export get_species, @setup_species, get_active_species, get_active_species_set

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