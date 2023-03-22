abstract type AbstractLoadedSpeciesIterators end
struct LoadedSpeciesIterators{S,I,N,M,A,T,U} <: AbstractLoadedSpeciesIterators
    all :: S
    ions :: I
    neutrals :: N
    atoms :: A
    molecules :: M
    imp_ions :: T
    imp_atoms :: U
    idx_e⁻ :: Int64
    idx_main_ion :: Int64
    idx_main_atom :: Int64
end
function convert2iterator(v::Vector{Int64})
    if isempty(v)
        return v
    elseif is_continuous(v)
        return UnitRange(minimum(v),maximum(v))
    else 
        return v
    end
end
is_continuous(v) = all(diff(v) .== 1)
function LoadedSpeciesIterators()
    all = convert2iterator(get_species_index())
    ions = convert2iterator(get_ions_index())
    neutrals = convert2iterator(get_neutrals_index())
    atoms = convert2iterator(get_atoms_index())
    molecules = convert2iterator(get_molecules_index())
    imp_ions = convert2iterator(get_imp_ions_index())
    imp_atoms = convert2iterator(get_imp_atoms_index())
    idx_e⁻ = get_electron_index()
    idx_main_ion = get_main_ion_index(;enforce = false)
    idx_main_atom = get_main_atom_index(;enforce = false)
    LoadedSpeciesIterators(all, ions, neutrals, atoms, molecules, imp_ions, imp_atoms, idx_e⁻, idx_main_ion, idx_main_atom)
end

get_species_iterators() = LoadedSpeciesIterators()

function Base.show(io::IO, spc_iter::LoadedSpeciesIterators)
    println("LoadedSpeciesIterators")
    for f in propertynames(spc_iter)
        println("└─ $f : $(getfield(spc_iter,f)) ")
    end
end
