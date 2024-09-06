module AtomicEnsemble

using FileIO
using PyCall

abstract type AbstractEnsemble end
abstract type AbstractStructure end

struct Structure{T} <: AbstractStructure
    positions :: Matrix{T}
    masses :: Vector{T}
    cell :: Matrix{T}
    atoms :: Vector{String}
end

positions(s :: Structure) = s.positions
masses(s :: Structure) = s.masses
cell(s :: Structure) = s.cell
atoms(s :: Structure) = s.atoms
Base.length(s :: Structure) = size(s.positions, 2)

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
mutable struct StandardEnsemble{T} <: AbstractEnsemble
    structures :: Vector{Structure{T}}
    energies :: Vector{T} # Energy for each configuration
    forces :: Array{T} # Forces for each configuratio
end


energies(ensemble :: StandardEnsemble) = ensemble.energies
forces(ensemble :: StandardEnsemble) = ensemble.forces
structures(ensemble :: StandardEnsemble) = ensemble.structures
Base.length(ensemble :: StandardEnsemble) = length(ensemble.structures)

export energies, forces, structures

include("standard_ensemble.jl")
include("asr.jl")

include("pycall_structure.jl")
include("pycall_ensemble.jl")

export Structure, StandardEnsemble, save, load_ensemble, 
       generate_standard_ensemble, copy_structure!, apply_asr!

end # module AtomicEnsemble
