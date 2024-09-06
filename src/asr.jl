@doc raw"""
    apply_asr!(ensemble :: AbstractEnsemble; apply_on_forces :: Bool = false)
    apply_asr!(s :: Structure)

Apply the acoustic sum rule to the structures (even the ensemble )(remove the total translations of the system)
If `apply_on_forces` is true, the forces are also modified. (Only if the ensemble is passed)
"""
function apply_asr!(ensemble :: AbstractEnsemble; apply_on_forces :: Bool = false)
    for s in structures(ensemble)
        apply_asr!(s; apply_on_forces = apply_on_forces)
    end
end
function apply_asr!(s :: Structure; apply_on_forces :: Bool = false)
    apply_asr!(positions(s))
end
function apply_asr!(positions :: Matrix{T}) where T
    n_atoms = size(positions, 2)
    n_dim = size(positions, 1)
    n_translations = n_dim
    for i in 1:n_translations
        @views positions[i, :] .-= mean(positions[i, :])
    end
end
