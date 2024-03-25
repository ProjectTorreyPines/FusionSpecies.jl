get_generated_gpu_struct_filename() = "$(@__DIR__)/generated_code/gpu.jl"
generate_gpu_struct = false

function clean_file!(fn)
    open(fn, "w") do io
        write(io, "")
    end
end

generate_gpu_struct && clean_file!(get_generated_gpu_struct_filename())
macro add_gpu_struct(ex)
    fn = get_generated_gpu_struct_filename()
    return add_gpu_struct(ex, fn)
end


function add_gpu_struct(ex, fn)
    !generate_gpu_struct && return ex
    expr = make_gpu_struct(ex)
    println("writing gpu struct into  $fn ....")
    write_expr(fn, expr; mode="a")
    return ex
end

function write_expr(fn, expr; mode="a")
    open(fn, mode) do io
        write(io, join(string.(s for s in expr.args), "\n"))
        write(io, "\n")
    end
end

function make_gpu_struct(ex)
    blk = Expr(:block)
    P = Symbol.([p * p * p for p in collect('A':'Z')])
    eq_struct = MacroTools.splitstructdef(ex)
    #new_struct = deepcopy(old_struct)
    #config_struct = deepcopy(old_struct)
    if eq_struct[:supertype] == :Any
        eq_struct[:supertype] = Symbol("Abstract" * string(eq_struct[:name]))
        push!(blk.args, :(abstract type $(eq_struct[:supertype]) end))
    end
    abstract_name = eq_struct[:supertype] isa Symbol ? eq_struct[:supertype] : copy(eq_struct[:supertype])
    if abstract_name isa Expr
        list_params_abstract = Tuple([abstract_name.args[2:end]]...)
    else
        list_params_abstract = ()
    end
    types = [f[2] for f in eq_struct[:fields]]
    fieldtypes = [p isa Expr ? p.args[1] : (p isa Symbol ? p : p[1]) for p in types]
    fieldnames = [f[1] for f in eq_struct[:fields]]
    new_fields = [(f, p) for (f, p) in zip(fieldnames, P)]
    #new_fields_config = [(f, p) for (f, p) in zip(fieldnames, P) if f != :log]
    old_params = [p isa Expr ? p.args[1] : (p isa Symbol ? p : p[1]) for p in eq_struct[:params]]

    fieldtypes = [p isa Expr ? p.args[1] : (p isa Symbol ? p : p[1]) for p in types]
    new_params = vcat([p for (t, p) in zip(types, P)], [p for p in old_params if p ∈ list_params_abstract])

    eq_struct[:fields] = new_fields
    eq_struct[:params] = new_params
    name = eq_struct[:name]
    eq_struct[:name] = Symbol(string(eq_struct[:name]) * "GPU")
    name_gpu = eq_struct[:name]

    if eq_struct[:supertype] isa Expr
        abstract_name_gpu_ = Symbol(string(eq_struct[:supertype].args[1]) * "GPU")
        eq_struct[:supertype].args[1] = abstract_name_gpu_
        extra_p = eq_struct[:supertype].args[2:end]
    else
        abstract_name_gpu_ = Symbol(string(eq_struct[:supertype]) * "GPU")
        eq_struct[:supertype] = abstract_name_gpu_
        extra_p = []
    end
    abstract_name_gpu = eq_struct[:supertype]

    push!(blk.args, :(abstract type $abstract_name_gpu <: $abstract_name end))
    push!(blk.args, :(export $(abstract_name_gpu isa Symbol ? abstract_name_gpu : abstract_name_gpu.args[1])))
    eq_struct[:mutable] = false
    push!(blk.args, MacroTools.combinestructdef(eq_struct))
    list_old_params = [p isa Expr ? p.args[1] : (p isa Symbol ? p : p[1]) for p in old_params]
    list_params = [p isa Expr ? p.args[1] : (p isa Symbol ? p : p[1]) for p in eq_struct[:params]]
    push!(blk.args, :(@nospecialize))
    body = quote end
    list_p = vcat([:(typeof($fn)) for fn in fieldnames], extra_p)
    # println("list_p=", list_p) 
    # println("fieldnames=", fieldnames)
    push!(body.args, Expr(:call, Expr(:curly, name_gpu, list_p...), :($(Expr(:parameters, :(kw...)))), fieldnames...))
    push!(body.args)
    func = Dict(
        :params => [],
        :name => abstract_name_gpu_,
        :args => vcat([:(@nospecialize($p)) for p in extra_p], [:(@nospecialize($f)) for (f, p) in new_fields]),#[Expr(:(::), f, p) for (f, p) in new_fields],
        :kwargs => Any[:(kw...)],
        :body => body,
        :whereparams => ()
    )
    push!(blk.args, MacroTools.combinedef(func))
    push!(blk.args, :(get_base_type_gpu(::$name) = $name_gpu))

    func = Dict(
        :params => [],
        :name => :(Adapt.adapt_structure),
        :args => [:to, :(v::$abstract_name)],
        :kwargs => Any[:(kw...)],
        :body => :($abstract_name_gpu_($(list_params_abstract...), (Adapt.adapt_structure(to, getproperty(v, f)) for f in propertynames(v))...)),
        :whereparams => list_params_abstract
    )
    push!(blk.args, MacroTools.combinedef(func))

    func = Dict(
        :params => [],
        :name => :special_copy,
        :args => [:(@nospecialize(v::$abstract_name)), :(@nospecialize(backend))],
        :kwargs => Any[],
        :body => :($abstract_name_gpu_($(list_params_abstract...), (special_copy(getproperty(v, f), backend) for f in propertynames(v))...)),
        :whereparams => list_params_abstract
    )

    push!(blk.args, MacroTools.combinedef(func))
    push!(blk.args, :(@specialize))

    expr = Base.remove_linenums!(blk)
    return expr
end
