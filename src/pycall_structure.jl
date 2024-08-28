


@doc raw"""
    Structure

    create from a CellConstructor structure object.
    A structure object that contains the positions, masses, unit cell, and atoms of a structure.
"""
function Structure(s :: PyCall.PyObject) :: Structure
    positions = copy(s.coords')
    masses = s.get_masses_array()
    cell = copy(s.unit_cell')
    atoms = s.atoms
    return Structure(positions, masses, cell, atoms)
end

function set_ase_positions!(ase_atoms :: PyCall.PyObject, structure :: Structure)
    ase_atoms.set_positions(positions(structure)')
end

@doc raw"""
    get_ase_atoms(structure :: Structure, ATM) :: PyObject

Create an ASE atoms object from a Structure object.
ATM is the Atoms module of the ASE package.
It can be loaded with PyCall using
```julia
using PyCall
ATM = pyimport("ase.atoms")
```
Or, alternatively,
```julia
using PyCall
@pyimport ase.atoms as ATM
```
"""
function get_ase_atoms(structure :: Structure, ATM) :: PyCall.PyObject
    return ATM.Atoms(atoms(structure), positions = positions(structure)', pbc = true, cell = cell(structure)')
end

@doc raw"""
    load_scf(scf_file :: String) :: Structure{Float64}

Load a structure from a Quantum Espresso scf file.
This file is defined in the python module cellconstructor.

ST must be the Structure module of the cellconstructor package.
It can be loaded with PyCall using
```julia
using PyCall
@pyimport cellconstructor.Structure as ST
```
"""
function load_scf(scf_file :: String, ST) :: Structure{Float64}
    structure = ST.Structure()
    structure.read_scf(scf_file)
    structure.build_masses()

    # Now convert the structure to a Structure object
    return Structure(structure)
end

function get_energy(ase_atoms :: PyCall.PyObject, calculator :: PyCall.PyObject)
    ase_atoms.set_calculator(calculator)
    return ase_atoms.get_total_energy()
end

@doc raw"""
    get_force!(forces :: Matrix{T}, structure :: Structure, calculator :: PyObject; ase_atoms = nothing)

Get the forces using the ASE calculator.
It is possible to provide an ase_atoms object commensurate with the provided structure, to avoid creating a new one.

The first argument ``forces`` are modified in-place storing a 3xN_atoms matrix for the structure.
"""
function get_force!(forces :: AbstractMatrix{T}, structure :: Structure, calculator :: PyCall.PyObject; ase_atoms = nothing) where {T} 
    if ase_atoms == nothing
        ase_atoms = get_ase_atoms(structure)
    end
    ase_atoms.set_calculator(calculator)
    set_ase_positions!(ase_atoms, structure)

    @views forces .= ase_atoms.get_forces()'
end

function get_forces!(forces :: Matrix{T}, structures :: Vector{Structure{T}}, calculator :: PyCall.PyObject) where {T}
    n_forces = size(forces, 3)
    ase_atoms = get_ase_atoms(structures[1])
    ase_atoms.set_calculator(calculator)

    for i in 1:n_forces
        set_ase_positions!(ase_atoms, structures[i])
        @views forces[:, :, i] .= ase_atoms.get_forces()'
    end
end


function n_atoms(structure :: Structure)
    return size(structure.positions, 1)
end
