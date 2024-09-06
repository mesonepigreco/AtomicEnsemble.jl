using Test
using AtomicEnsemble


@testset "Load And Save" begin
    filename = "ensemble.jld2"
    new_filename = "ensemble_bis.jld2"

    # Load the ensemble
    ensemble = load_ensemble(filename)

    # Save the ensemble
    save(new_filename, ensemble)

    # Load the new ensemble
    ensemble_bis = load_ensemble(new_filename)

    # Now compare the two ensemble
    @test length(ensemble) == length(ensemble_bis)
    for i in 1:length(ensemble)
        @test AtomicEnsemble.energies(ensemble)[i] ≈ AtomicEnsemble.energies(ensemble_bis)[i]
        for k in 1:length(structures(ensemble)[i])
            for h in 1:3
                @test AtomicEnsemble.positions(structures(ensemble)[i])[h, k] ≈ AtomicEnsemble.positions(structures(ensemble_bis)[i])[h, k]
            end
            @test AtomicEnsemble.masses(structures(ensemble)[i])[k] ≈ AtomicEnsemble.masses(structures(ensemble_bis)[i])[k]
            @test AtomicEnsemble.atoms(structures(ensemble)[i])[k] == AtomicEnsemble.atoms(structures(ensemble_bis)[i])[k]
        end

        for h in 1:3
            for k in 1:3
                @test AtomicEnsemble.cell(structures(ensemble)[i])[h, k] ≈ AtomicEnsemble.cell(structures(ensemble_bis)[i])[h, k]
            end
        end
    end
end

@testset "Apply ASR" include("test_asr.jl")

