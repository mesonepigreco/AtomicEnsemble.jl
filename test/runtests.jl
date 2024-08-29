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
        @test energies(ensemble)[i] ≈ energies(ensemble_bis)[i]
        for k in 1:length(structures(ensemble)[i])
            for h in 1:3
                @test cell(structures(ensemble)[i])[h, k] ≈ cell(structures(ensemble_bis)[i])[h, k]
                @test positions(structures(ensemble)[i])[h, k] ≈ positions(structures(ensemble_bis)[i])[h, k]
            end
            @test masses(structures(ensemble)[i])[k] ≈ masses(structures(ensemble_bis)[i])[k]
            @test atoms(structures(ensemble)[i])[k] == atoms(structures(ensemble_bis)[i])[k]
        end
    end
end
