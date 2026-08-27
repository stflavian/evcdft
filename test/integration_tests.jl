"""
Integration tests for EVC_DFT package.

These tests verify that the complete DFT workflow works as expected,
testing the interaction between different modules.
"""

using Test
using LinearAlgebra

# Include the main module
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "EVC_DFT.jl"))

# Use the module
using .EVC_DFT
using EVC_DFT.Constants
using EVC_DFT.Types
using EVC_DFT.PlaneWave
using EVC_DFT.XCFunctionals
using EVC_DFT.SelfConsistent

# =============================================================================
# INTEGRATION TESTS FOR COMPLETE DFT WORKFLOW
# =============================================================================

@testset "Complete DFT Workflow" begin
    @testset "Uniform Electron Gas - Small System" begin
        # Create a small uniform electron gas system
        a = 5.0  # Bohr
        n_electrons = 2
        cutoff = 10.0  # Hartree
        fft_size = (8, 8, 8)
        
        # Create system
        system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
        
        # Verify system properties
        @test system.lattice.volume ≈ a^3 rtol=1e-10
        @test system.electrons == n_electrons
        @test system.basis.cutoff ≈ cutoff rtol=1e-10
        @test system.basis.fft_size == fft_size
        
        # Initialize density
        initialize_uniform_density(system)
        expected_density = n_electrons / (a^3)
        @test isapprox(system.density.data, fill(expected_density, size(system.density.data)); atol=1e-10)
        
        # Compute total energy
        total_energy = compute_total_energy(system)
        @test isfinite(total_energy)
        @test total_energy < 0.0
        
        # Run SCF
        params = SCFParameters(max_iter=10, energy_tolerance=1e-6, density_tolerance=1e-6)
        converged_system = self_consistent_field(system, params)
        
        # Verify convergence
        @test isfinite(converged_system.energies.total)
        @test converged_system.energies.total < 0.0
        
        # For uniform electron gas, density should remain uniform
        @test isapprox(converged_system.density.data, fill(expected_density, size(converged_system.density.data)); atol=1e-8)
    end
    
    @testset "Uniform Electron Gas - Medium System" begin
        # Create a medium-sized uniform electron gas system
        a = 10.0  # Bohr
        n_electrons = 4
        cutoff = 15.0  # Hartree
        fft_size = (16, 16, 16)
        
        system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
        
        # Initialize and run SCF
        initialize_uniform_density(system)
        params = SCFParameters(max_iter=10, energy_tolerance=1e-6, density_tolerance=1e-6)
        converged_system = run_scf!(system, params)
        
        # Verify results
        @test isfinite(converged_system.energies.total)
        @test converged_system.energies.total < 0.0
        
        # Check energy components
        @test isfinite(converged_system.energies.hartree)
        @test isfinite(converged_system.energies.exchange)
    end
    
    @testset "Uniform Electron Gas - Different Densities" begin
        # Test with different electron densities
        densities = [1, 2, 4]  # electrons per cell
        
        for n_electrons in densities
            a = 5.0
            cutoff = 10.0
            fft_size = (8, 8, 8)
            
            system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
            initialize_uniform_density(system)
            
            params = SCFParameters(max_iter=5, energy_tolerance=1e-4, density_tolerance=1e-4)
            converged_system = run_scf!(system, params)
            
            @test isfinite(converged_system.energies.total)
            
            # Higher density should have higher (more positive) energy per electron
            # due to increased kinetic energy
            expected_density = n_electrons / (a^3)
            @test isapprox(converged_system.density.data, fill(expected_density, size(converged_system.density.data)); atol=1e-8)
        end
    end
    
    @testset "Uniform Electron Gas - Different Cell Sizes" begin
        # Test with different cell sizes (same density)
        n_electrons = 2
        density = 0.01  # electrons/Bohr^3
        
        for a in [5.0, 10.0]
            cutoff = 10.0
            fft_size = (8, 8, 8)
            
            system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
            
            # Expected density
            expected_density = n_electrons / (a^3)
            
            # For uniform electron gas, the energy per electron should be the same
            # for the same rs value
            rs = (3.0 / (4 * pi * expected_density))^(1/3)
            expected_energy_per_electron = jellium_energy_per_electron(rs)
            
            # Initialize and compute energy
            initialize_uniform_density(system)
            total_energy = jellium_total_energy(system)
            
            @test total_energy ≈ n_electrons * expected_energy_per_electron rtol=1e-10
        end
    end
end

@testset "Density Functional Components" begin
    @testset "Hartree + XC Potential Calculation" begin
        a = 5.0
        n_electrons = 2
        cutoff = 10.0
        fft_size = (8, 8, 8)
        
        system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
        initialize_uniform_density(system)
        
        # Compute Hartree potential
        hartree_pot = compute_hartree_potential(system.density, system.basis)
        system.potential.hartree .= hartree_pot
        
        # Compute XC potential
        volume = system.lattice.volume
        xc_energy, xc_pot = compute_lda_xc(system.density, volume)
        system.potential.exchange .= xc_pot
        system.potential.correlation .= xc_pot
        
        # Check that potentials are computed
        @test !all(system.potential.hartree .== 0.0)
        @test !all(system.potential.exchange .== 0.0)
        
        # For uniform density, XC potential should be uniform
        @test isapprox(xc_pot, fill(xc_pot[1,1,1], size(xc_pot)); atol=1e-10)
    end
    
    @testset "Energy Components Breakdown" begin
        a = 5.0
        n_electrons = 2
        cutoff = 10.0
        fft_size = (8, 8, 8)
        
        system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
        initialize_uniform_density(system)
        
        # Manually compute energy components
        volume = system.lattice.volume
        
        # Hartree energy
        hartree_pot = compute_hartree_potential(system.density, system.basis)
        hartree_energy = compute_hartree_energy(system.density.data, hartree_pot, volume)
        
        # XC energy
        xc_energy, xc_pot = compute_lda_xc(system.density, volume)
        
        # Total energy (simplified - no kinetic energy)
        total_energy = hartree_energy + xc_energy
        
        # Compare with system's compute_total_energy
        system_total = compute_total_energy(system)
        
        # They should be close (within numerical precision)
        @test total_energy ≈ system_total rtol=1e-10
    end
end

@testset "Mixing Strategies" begin
    @testset "Linear Mixing Convergence" begin
        a = 5.0
        n_electrons = 2
        cutoff = 10.0
        fft_size = (8, 8, 8)
        
        system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
        initialize_uniform_density(system)
        
        # Use linear mixing with different parameters
        for alpha in [0.3, 0.5, 0.7]
            params = SCFParameters(
                max_iter=10,
                energy_tolerance=1e-6,
                density_tolerance=1e-6,
                mixing_parameter=alpha,
                mixing_type="linear"
            )
            
            converged_system = run_scf!(system, params)
            
            @test isfinite(converged_system.energies.total)
            @test converged_system.energies.total < 0.0
        end
    end
    
    @testset "Kerker Mixing" begin
        a = 5.0
        n_electrons = 2
        cutoff = 10.0
        fft_size = (8, 8, 8)
        
        system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
        initialize_uniform_density(system)
        
        params = SCFParameters(
            max_iter=10,
            energy_tolerance=1e-6,
            density_tolerance=1e-6,
            mixing_parameter=0.5,
            mixing_type="kerker"
        )
        
        converged_system = run_scf!(system, params)
        
        @test isfinite(converged_system.energies.total)
    end
end

@testset "System Creation and Manipulation" begin
    @testset "Create and Modify System" begin
        # Create initial system
        a = 5.0
        n_electrons = 2
        cutoff = 10.0
        fft_size = (8, 8, 8)
        
        system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
        
        # Modify density
        system.density.data .= 0.2
        
        # Compute Hartree potential with modified density
        hartree_pot = compute_hartree_potential(system.density, system.basis)
        
        @test size(hartree_pot) == fft_size
        @test !all(hartree_pot .== 0.0)
        
        # Reset to uniform density
        initialize_uniform_density(system)
        expected_density = n_electrons / (a^3)
        @test isapprox(system.density.data, fill(expected_density, size(system.density.data)); atol=1e-10)
    end
    
    @testset "Multiple Systems" begin
        # Create multiple systems with different parameters
        systems = []
        for n_electrons in [1, 2, 4]
            a = 5.0
            cutoff = 10.0
            fft_size = (8, 8, 8)
            
            system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
            push!(systems, system)
        end
        
        # Verify all systems are independent
        for (i, system) in enumerate(systems)
            initialize_uniform_density(system)
            expected_density = system.electrons / (system.lattice.volume)
            @test isapprox(system.density.data, fill(expected_density, size(system.density.data)); atol=1e-10)
        end
    end
end

@testset "Physical Properties" begin
    @testset "Jellium Model" begin
        # Test jellium model properties
        a = 5.0
        n_electrons = 2
        
        lattice = SimpleCubic(a)
        ueg = UniformElectronGas(lattice, n_electrons)
        
        # Check Wigner-Seitz radius
        expected_rs = (3.0 / (4 * pi * (n_electrons / a^3)))^(1/3)
        @test ueg.rs ≈ expected_rs rtol=1e-10
        
        # Check energy per electron
        energy_per_e = jellium_energy_per_electron(ueg.rs)
        @test isfinite(energy_per_e)
        
        # Total energy
        total_energy = n_electrons * energy_per_e
        @test isfinite(total_energy)
    end
    
    @testset "Energy Scaling" begin
        # Test that energy scales correctly with system size
        n_electrons = 2
        
        for a in [5.0, 10.0]
            cutoff = 10.0
            fft_size = (8, 8, 8)
            
            system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
            initialize_uniform_density(system)
            
            # Compute jellium energy
            jellium_energy = jellium_total_energy(system)
            
            # The energy per electron should be the same for the same rs
            rs = (3.0 / (4 * pi * (n_electrons / a^3)))^(1/3)
            expected_energy_per_e = jellium_energy_per_electron(rs)
            expected_total = n_electrons * expected_energy_per_e
            
            @test jellium_energy ≈ expected_total rtol=1e-10
        end
    end
end

@testset "Numerical Stability" begin
    @testset "Small FFT Grids" begin
        # Test with very small FFT grids
        for fft_size in [(4, 4, 4), (2, 2, 2)]
            a = 5.0
            n_electrons = 1
            cutoff = 5.0
            
            system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
            initialize_uniform_density(system)
            
            params = SCFParameters(max_iter=5, energy_tolerance=1e-3, density_tolerance=1e-3)
            converged_system = run_scf!(system, params)
            
            @test isfinite(converged_system.energies.total)
        end
    end
    
    @testset "Large FFT Grids" begin
        # Test with larger FFT grids (within memory limits)
        a = 5.0
        n_electrons = 2
        cutoff = 10.0
        fft_size = (32, 32, 32)
        
        system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
        initialize_uniform_density(system)
        
        params = SCFParameters(max_iter=3, energy_tolerance=1e-3, density_tolerance=1e-3)
        converged_system = run_scf!(system, params)
        
        @test isfinite(converged_system.energies.total)
    end
    
    @testset "Extreme Densities" begin
        # Test with very high and low densities
        a = 5.0
        cutoff = 10.0
        fft_size = (8, 8, 8)
        
        for n_electrons in [1, 2, 4]
            system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
            initialize_uniform_density(system)
            
            params = SCFParameters(max_iter=3, energy_tolerance=1e-3, density_tolerance=1e-3)
            converged_system = run_scf!(system, params)
            
            @test isfinite(converged_system.energies.total)
        end
    end
end

@testset "Error Handling" begin
    @testset "Invalid Parameters" begin
        # Test with invalid parameters (should still work due to defaults)
        a = 5.0
        n_electrons = 2
        cutoff = 10.0
        fft_size = (8, 8, 8)
        
        # Negative mixing parameter (should still work)
        params = SCFParameters(mixing_parameter=-0.5)
        system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
        initialize_uniform_density(system)
        
        # This might not converge well, but shouldn't crash
        try
            converged_system = run_scf!(system, params)
            @test true  # If it completes without error
        catch
            @test false  # Should not throw
        end
    end
    
    @testset "Zero Electrons" begin
        # Test with zero electrons
        a = 5.0
        n_electrons = 0
        cutoff = 10.0
        fft_size = (8, 8, 8)
        
        system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
        initialize_uniform_density(system)
        
        @test all(system.density.data .== 0.0)
        
        # XC functionals should return 0 for zero density
        @test lda_xc_energy(0.0) == 0.0
        @test lda_xc_potential(0.0) == 0.0
    end
end

# Print summary
println("\n=== Integration Tests Complete ===")
