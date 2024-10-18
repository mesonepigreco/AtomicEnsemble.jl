@doc raw"""
    atomic_units_strip(ens :: Ensemble) :: Ensemble

Convert the ensemble into atomic units.
"""
function atomic_units_strip(ens :: StandardEnsemble) :: StandardEnsemble
    # Generate the new structures
    new_structures = []
    for s in structures(ens)
        new_positions = ustrip.(auconvert.(positions(s)))
        new_cell = ustrip.(auconvert.(cell(s)))
        new_masses = ustrip.(auconvert.(masses(s)))
        new_structure = Structure(new_positions, new_masses, new_cell, atoms(s))
        push!(new_structures, new_structure)
    end

    new_energies = ustrip.(auconvert.(energies(ens)))
    new_forces = ustrip.(auconvert.(forces(ens)))

    StandardEnsemble(new_structures, new_energies, new_forces)
end
