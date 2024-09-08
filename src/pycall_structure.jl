


@doc raw"""
    Structure

create from a CellConstructor structure object.
A structure object that contains the positions, masses, unit cell, and atoms of a structure.

Using Unitful to specify the correct units
"""
function Structure(s :: PyCall.PyObject) :: Structure
    nat = s.N_atoms
    positions = zeros(1.0u"Å", 3, nat)
    for i in 1:nat
        positions[:, i] .= s.coords[i, :] * u"Å"
    end

    masses = s.get_masses_array() .* auconvert(m_u)
    cell = copy(s.unit_cell') .* u"Å"
    atoms = s.atoms
    return Structure(positions, masses, cell, atoms)
end

function set_ase_positions!(ase_atoms :: PyCall.PyObject, structure :: Structure)
    ase_atoms.set_positions(ustrip(uconvert.(u"Å", positions(structure)')))
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
    return ATM.Atoms(atoms(structure), positions = ustrip(uconvert.(u"Å", positions(structure)')), pbc = true, cell = ustrip(uconvert.(u"Å", cell(structure)')))
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

@doc raw"""
get_energy(ase_atoms :: PyCall.PyObject, calculator :: PyCall.PyObject) :: Quantity

Use unitful for returning a quantity.
"""
function get_energy(ase_atoms :: PyCall.PyObject, calculator :: PyCall.PyObject)
    ase_atoms.set_calculator(calculator)
    return ase_atoms.get_total_energy() * u"eV"
end

@doc raw"""
    get_force!(forces :: Matrix{Quantity}, structure :: Structure, calculator :: PyObject; ase_atoms = nothing)

Get the forces using the ASE calculator.
It is possible to provide an ase_atoms object commensurate with the provided structure, to avoid creating a new one.

The first argument ``forces`` are modified in-place storing a 3xN_atoms matrix for the structure.

Notably, the forces must be a unitful type
"""
function get_force!(forces :: AbstractMatrix{Quantity}, structure :: Structure, calculator :: PyCall.PyObject; ase_atoms = nothing)
    if ase_atoms == nothing
        ase_atoms = get_ase_atoms(structure)
    end
    ase_atoms.set_calculator(calculator)
    set_ase_positions!(ase_atoms, structure)

    forces .= ase_atoms.get_forces()' * u"eV/Å"
end

function get_forces!(forces :: Matrix, structures :: Vector{Structure}, calculator :: PyCall.PyObject)
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
