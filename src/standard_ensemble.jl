@doc """
    enerate_standard_ensemble(n_configs :: Int, n_atoms :: Int, copy_structure :: Structure{T}) :: Ensemble{T}

Generate a standard ensemble with the same structure for all configurations.
"""
function generate_standard_ensemble(n_configs :: Int, copy_structure :: Structure{T}) :: StandardEnsemble{T} where {T}
    n_atoms = length(copy_structure)
    structures = Vector{Structure{T}}(undef, n_configs)
    energies = zeros(T, n_configs)
    forces = zeros(T, 3, n_atoms, n_configs)

    for i in 1:n_configs
        structure = Structure{T}(positions(copy_structure), masses(copy_structure), cell(copy_structure), atoms(copy_structure))
        structures[i] = structure
    end

    return StandardEnsemble(structures, energies, forces)
end


function convert_to_dict(ensemble :: AbstractEnsemble)
    structs = [convert_to_dict(structure) for structure in structures(ensemble)]
    frc = [forces(ensemble)[:, :, i] for i in 1:length(ensemble)]
    ensemble_dict = Dict(:structures => structs,
                         :energies => energies(ensemble),
                         :forces => frc)
    return ensemble_dict
end

function convert_to_dict(structure :: AbstractStructure)
    struct_dict = Dict(:atoms => atoms(structure), :cell => cell(structure), :masses => masses(structure), :positions => positions(structure)) 
    return struct_dict
end

function convert_to_structure(structure_dict :: Dict) :: Structure
    return Structure(structure_dict[:positions], structure_dict[:masses], structure_dict[:cell], structure_dict[:atoms])
end
function convert_to_ensemble(ensemble_dict :: Dict) :: StandardEnsemble
    structures = [convert_to_structure(structure) for structure in ensemble_dict[:structures]]
    forces = zeros(Float64, 3, n_atoms(structures[1]), length(structures))
    for i in 1:length(structures)
        forces[:, :, i] = ensemble_dict[:forces][i]
    end
    return StandardEnsemble(structures, ensemble_dict[:energies], forces)
end


@doc raw"""
    save(filename :: String, ensemble :: AbstractEnsemble)

Save an ensemble to a file.
The extension is managed by FileIO, so any supported version can be used.
One of the recommended formats is JLD2.
This performs a function overload on the save function from FileIO.

Note that to load the ensemble, the function load_ensemble should be used.
(load only loads the dictionary, not the ensemble object)

## Example
```julia
using FileIO
ensemble = StandardEnsemble([Structure(rand(3, 3), rand(3), rand(3, 3), ["H", "H", "H"])], [0.0], rand(3, 3, 1))
save("ensemble.jld2", ensemble)
```

"""
function FileIO.save(filename :: String, ensemble :: AbstractEnsemble)
    ensemble_dict = convert_to_dict(ensemble)
    save(filename, ensemble_dict)
end

@doc raw"""
    load_ensemble(filename :: String) :: AbstractEnsemble

Load an ensemble from a file.
"""
function load_ensemble(filename :: String) :: AbstractEnsemble
    ensemble_dict = load(filename)
    return convert_to_ensemble(ensemble_dict)
end
