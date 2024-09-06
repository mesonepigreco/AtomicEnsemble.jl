@doc raw"""
    apply_asr!(ensemble :: AbstractEnsemble; apply_on_forces :: Bool = false)
    apply_asr!(s :: Structure)
    apply_asr!(positions :: Matrix{T}) where T
    apply_asr!(vector :: AbstractVecctor{T}; ndim :: Int = 3) where {T}

Apply the acoustic sum rule to the structures (even the ensemble )(remove the total translations of the system)
If `apply_on_forces` is true, the forces are also modified. (Only if the ensemble is passed)

In the case of a vector, the vector is assumed to be a vector of positions, 
with the following structure: [x1, y1, z1, x2, y2, z2, ...]. 
The vector is modified in place.
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
        @views positions[i, :] .-= sum(positions[i, :]) ./ n_atoms
    end
end
function apply_asr!(vector :: AbstractVector{T}; ndim :: Int = 3) where {T}
    n_atoms = length(vector) ÷ ndim
    global_translation = zeros(T, ndim)
    for i in 1:ndim
        for j in 1:n_atoms
            index = (j - 1) * ndim + i
            global_translation[i] += vector[index]
        end
    end
    global_translation ./= n_atoms

    for i in 1:ndim
        for j in 1:n_atoms
            index = (j - 1) * ndim + i
            vector[index] -= global_translation[i]
        end
    end
end
