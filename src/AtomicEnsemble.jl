module AtomicEnsemble

using FileIO

abstract type AbstractEnsemble end
abstract type AbstractStructure end

struct Structure{T} <: AbstractStructure
    positions :: Matrix{T}
    masses :: Vector{T}
    cell :: Matrix{T}
    atoms :: Vector{String}
end

position(s :: Structure) = s.positions
masses(s :: Structure) = s.masses
cell(s :: Structure) = s.cell
atoms(s :: Structure) = s.atoms
length(s :: Structure) = size(s.positions, 2)

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
n_structures(ensemble :: StandardEnsemble) = length(ensemble.structures)

include("standard_ensemble.jl")

# Import the pycall interaction if PyCall is defined
@static if isdefined(Main, :PyCall) 
    import PyCall
    
    # Check whether pycall can import cellconstructor and sscha
    try
        @pyimport cellconstructor as CC
        @pyimport cellconstructor.Structure as ST
        @pyimport cellconstructor.Phonons as PH
        @pyimport sscha
        @pyimport sscha.Ensemble as ENS
    catch
        @warn "PyCall can't import cellconstructor or sscha. Structure will not be able to interact with ASE"
    end

    include("pycall_structure.jl")
    include("pycall_ensemble.jl")
else
    @warn "PyCall not defined. Structure will not be able to interact with ASE"
end


export Structure, StandardEnsemble, save, load_ensemble

end # module AtomicEnsemble
