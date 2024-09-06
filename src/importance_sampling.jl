@doc """
   get_response_function_is(ensemble :: StandardEnsemble{T}) :: Matrix{T} where T

Get the ionic static response function of the ensemble using importance sampling.
"""
function get_response_function_is!(chi :: Matrix{T}, ensemble :: AbstractEnsemble{T}, β :: T) where {T}
    n_structs = length(ensemble)
    n_atms = length(structures(ensemble)[1])
    n_dims = 3

    # Apply the asr on the ensemble
    apply_asr!(ensemble)

    function Δx(λ :: Vector{T}) :: Vector{T}
        Δx = zeros(T, n_dims * n_atms)
        weights = zeros(T, n_structs)

        for j in 1:n_structs
            energy = 0
            for i in 1:n_dims * n_atms
                i_at = div(i - 1, n_dims) + 1
                i_dir = mod(i - 1, n_dims) + 1

                energy += (λ[i] * structures(ensemble[j])[i_dir, i_at])
            end
            weights[j] = exp(-β * energy)
        end
        weights /= sum(weights)

        # Now compute the Δx
        Δx .= 0 
        for j in 1:n_structs
            for i in 1:n_dims * n_atms
                i_at = div(i - 1, n_dims) + 1
                i_dir = mod(i - 1, n_dims) + 1

                Δx[i] += (weights[j] - 1.0 / n_structs) * structures(ensemble[j])[i_dir, i_at]
            end
        end
        
        return Δx

    end

end
