function StandardEnsemble(py_ensemble :: PyCall.PyObject)
    structures = zeros(Structure{Float64}, py_ensemble.N)
    for i in 1:py_ensemble.N
        structures[i] = Structure(py_ensemble.structures[i])
    end

    forces = zeros(Float64, 3, n_atoms(structures[1]), py_ensemble.N)
    energies = zeros(Float64, py_ensemble.N)

    for i in 1:py_ensemble.N
        forces[:, :, i] = py_ensemble.forces[i, :].reshape((n_atoms(structures[1]), 3))'
    end
    energies .= py_ensemble.energies
    return StandardEnsemble(structures, energies, forces)
end

    
function load_ase_trajectory(filename :: String, ase_io) :: StandardEnsemble
    ase_ensemble = ase_io.read(filename)
    n_structures = length(ase_ensemble)

    positions = zeros(Structure{Float64}, n_structures)
    forces = zeros(Float64, 3, length(ase_ensemble[1]), n_structures)
    energies = zeros(Float64, n_structures)
    masses = zeros(Float64, length(ase_ensemble[1]))

    conv_forc = ustrip(auconvert(1.0u"eV/Å"))
    conv_ener = ustrip(auconvert(1.0u"eV"))

    for i in 1:n_structures
        tmp_struct = get_from_ase_atoms(ase_ensemble[i])
        tmp_force = ase_ensemble[i].get_forces()
        tmp_energy = ase_ensemble[i].get_potential_energy()
        
        positions[:, i] = ustrip.(auconvert.(tmp_struct.positions))
        forces[:, :, i] = tmp_force' .* conv_forc
        energies[i] = tmp_energy * conv_ener
    end

    return StandardEnsemble(positions, energies, forces)
end



