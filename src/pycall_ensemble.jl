
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

    
