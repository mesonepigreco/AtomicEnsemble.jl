module AtomicEnsemble

using FileIO
using PyCall

using Unitful, UnitfulAtomic
using PhysicalConstants.CODATA2018: k_B, ħ, m_u

abstract type AbstractEnsemble end
abstract type AbstractStructure end

struct Structure{T, U} <: AbstractStructure
    positions :: Matrix{T}
    masses :: Vector{U}
    cell :: Matrix{T}
    atoms :: Vector{String}
end

positions(s :: Structure) = s.positions
masses(s :: Structure) = s.masses
cell(s :: Structure) = s.cell
atoms(s :: Structure) = s.atoms
Base.length(s :: Structure) = size(s.positions, 2)

@doc raw"""
    get_atomic_types_int(s :: Structure) :: Vector{Int}

Return the atomic types as integers.
"""
function get_atomic_types_int(s :: Structure)
    atomic_types = atoms(s)
    unique_types = unique(atomic_types)
    atomic_types_int = zeros(Int, length(atomic_types))
    for i in 1:length(atomic_types)
        atomic_types_int[i] = findfirst(x -> x == atomic_types[i], unique_types)
    end
    return atomic_types_int
end

@doc raw"""
    copy_structure!(target :: Structure, origin :: Structure)

Copy the structure `origin` to the structure `target`.
"""
function copy_structure!(target :: Structure, origin :: Structure) 
    @views positions(target) .= positions(origin)
    @views masses(target) .= masses(origin)
    @views cell(target) .= cell(origin)
    @views atoms(target) .= atoms(origin)
end


@doc raw"""
    StandardEnsemble

Stores an ensemble of structures. 
If PyCall is available, you can interact with the sscha.Ensemble object from
python.
"""
mutable struct StandardEnsemble <: AbstractEnsemble
    structures :: Vector{Structure}
    energies :: Vector # Energy for each configuration
    forces :: Array # Forces for each configuratio
end


energies(ensemble :: StandardEnsemble) = ensemble.energies
forces(ensemble :: StandardEnsemble) = ensemble.forces
structures(ensemble :: StandardEnsemble) = ensemble.structures
Base.length(ensemble :: StandardEnsemble) = length(ensemble.structures)

export energies, forces, structures

include("standard_ensemble.jl")
include("unitful_ensemble.jl")
include("asr.jl")

include("pycall_structure.jl")
include("pycall_ensemble.jl")


export Structure, StandardEnsemble, save, load_ensemble, 
       atomic_units_strip,
       generate_standard_ensemble, copy_structure!, apply_asr!,
       cell, positions, energies, forces, structures, atoms, n_atoms

end # module AtomicEnsemble
