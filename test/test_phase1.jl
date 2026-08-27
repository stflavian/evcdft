"""
Tests for Phase 1 of the DFT implementation.

This file contains unit tests for:
- Core types and constants
- Plane wave basis
- LDA exchange-correlation functionals
- SCF loop for uniform electron gas

For comprehensive unit tests, see unit_tests.jl
For integration tests, see integration_tests.jl
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

@testset "Core Constants" begin
    # Test that constants are defined
    @test hbar == 1.0
    @test m_e == 1.0
    @test e == 1.0
    
    # Test conversion factors
    @test hartree_to_ev ≈ 27.211386245988
    @test ev_to_hartree ≈ 1.0 / 27.211386245988
    
    # Test Bohr radius
    @test a0 ≈ 0.5291772109038427 * 1e-10 * 1.8897259886
end

@testset "Core Types" begin
    # Test Lattice construction
    a1 = [1.0, 0.0, 0.0]
    a2 = [0.0, 1.0, 0.0]
    a3 = [0.0, 0.0, 1.0]
    lattice = Lattice(a1, a2, a3)
    
    @test lattice.volume ≈ 1.0
    @test lattice.b1 ≈ [2 * pi, 0.0, 0.0]
    @test lattice.b2 ≈ [0.0, 2 * pi, 0.0]
    @test lattice.b3 ≈ [0.0, 0.0, 2 * pi]
    
    # Test SimpleCubic
    sc_lattice = SimpleCubic(2.0)
    @test sc_lattice.volume ≈ 8.0
    
    # Test UniformElectronGas
    ueg = UniformElectronGas(sc_lattice, 2)
    @test ueg.n_electrons == 2
    @test ueg.density ≈ 2.0 / 8.0
    @test ueg.rs ≈ (3.0 / (4 * pi * (2.0 / 8.0)))^(1/3)
    
    # Test ElectronDensity
    density = ElectronDensity((16, 16, 16))
    @test size(density.data) == (16, 16, 16)
    @test density.grid_size == (16, 16, 16)
    
    # Test KohnShamPotential
    potential = KohnShamPotential((16, 16, 16))
    @test size(potential.hartree) == (16, 16, 16)
    @test size(potential.exchange) == (16, 16, 16)
    
    # Test SCFParameters
    params = SCFParameters()
    @test params.max_iter == 100
    @test params.energy_tolerance ≈ 1e-6
    @test params.density_tolerance ≈ 1e-6
    @test params.mixing_parameter == 0.5
    @test params.mixing_type == "linear"
    
    # Test custom SCFParameters
    custom_params = SCFParameters(max_iter=50, energy_tolerance=1e-8, 
                                   density_tolerance=1e-8, mixing_parameter=0.7,
                                   mixing_type="kerker")
    @test custom_params.max_iter == 50
    @test custom_params.energy_tolerance ≈ 1e-8
    @test custom_params.mixing_parameter == 0.7
    @test custom_params.mixing_type == "kerker"
end

@testset "Plane Wave Basis" begin
    # Create a simple cubic lattice
    lattice = SimpleCubic(5.0)
    
    # Create plane wave basis with small cutoff
    cutoff = 5.0  # Hartree
    fft_size = (16, 16, 16)
    basis = PlaneWaveBasis(lattice, cutoff, fft_size)
    
    # Test that basis was created
    @test basis.lattice == lattice
    @test basis.cutoff ≈ cutoff
    @test basis.fft_size == fft_size
    
    # Test that G vectors were generated
    @test length(basis.g_vectors) > 0
    @test length(basis.g2) == length(basis.g_vectors)
    
    # Test kinetic energy computation
    kinetic = compute_kinetic_energy(basis)
    @test length(kinetic) == length(basis.g2)
    @test all(kinetic .>= 0.0)
    
    # Test FFT operations
    test_data = rand(Float64, 8, 8, 8)
    test_recip = fft_forward(test_data)
    @test size(test_recip) == size(test_data)
    
    test_back = fft_backward(test_recip)
    @test size(test_back) == size(test_data)
    
    # Test Hartree potential computation
    density = ElectronDensity((8, 8, 8))
    density.data .= 0.1  # Uniform density
    
    # Create a basis for this density
    small_lattice = SimpleCubic(5.0)
    small_basis = PlaneWaveBasis(small_lattice, 10.0, (8, 8, 8))
    
    hartree_pot = compute_hartree_potential(density, small_basis)
    @test size(hartree_pot) == size(density.data)
    
    # For uniform density, Hartree potential should be uniform (in periodic boundary conditions)
    # Actually, for a uniform density in a periodic box, the Hartree potential
    # should be zero because the positive background cancels the electron density
    # This is a property of jellium
end

@testset "LDA Exchange-Correlation" begin
    # Test exchange energy
    density = 0.1  # Bohr^-3
    ex_energy = lda_exchange_energy(density)
    @test ex_energy < 0.0  # Exchange energy should be negative
    
    # Test correlation energy
    cor_energy = lda_correlation_energy(density)
    @test cor_energy < 0.0  # Correlation energy should be negative
    
    # Test XC energy
    xc_energy = lda_xc_energy(density)
    @test xc_energy < 0.0
    @test xc_energy ≈ ex_energy + cor_energy
    
    # Test exchange potential
    ex_pot = lda_exchange_potential(density)
    @test ex_pot < 0.0
    
    # Test correlation potential
    cor_pot = lda_correlation_potential(density)
    @test cor_pot < 0.0
    
    # Test XC potential
    xc_pot = lda_xc_potential(density)
    @test xc_pot < 0.0
    @test xc_pot ≈ ex_pot + cor_pot
    
    # Test with zero density
    @test lda_exchange_energy(0.0) == 0.0
    @test lda_correlation_energy(0.0) == 0.0
    @test lda_exchange_potential(0.0) == 0.0
    @test lda_correlation_potential(0.0) == 0.0
    
    # Test with negative density (should return 0)
    @test lda_exchange_energy(-0.1) == 0.0
    @test lda_correlation_energy(-0.1) == 0.0
end

@testset "Density Mixing" begin
    # Test linear mixing
    old_dens = fill(0.1, 8, 8, 8)
    new_dens = fill(0.2, 8, 8, 8)
    
    mixed = linear_mixing(new_dens, old_dens, 0.5)
    @test all(mixed .≈ 0.15)
    
    mixed = linear_mixing(new_dens, old_dens, 0.7)
    @test all(mixed .≈ 0.17)
    
    # Test Kerker mixing (should fall back to linear for now)
    mixed = kerker_mixing(new_dens, old_dens, 0.5)
    @test all(mixed .≈ 0.15)
    
    # Test apply_mixing
    params = SCFParameters(mixing_type="linear", mixing_parameter=0.5)
    mixed = apply_mixing(new_dens, old_dens, params)
    @test all(mixed .≈ 0.15)
end

@testset "Jellium Energy" begin
    # Test jellium energy calculation
    rs = 2.0  # Wigner-Seitz radius
    energy = jellium_energy_per_electron(rs)
    
    # Known values for rs=2:
    # Kinetic: 2.8376 / 4 = 0.7094
    # Exchange: -0.9163 / 2 = -0.45815
    # Correlation: approximately -0.14 (from Perdew-Zunger)
    # Total: approximately 0.7094 - 0.45815 - 0.14 ≈ 0.111
    @test energy ≈ 0.111
    
    # Test with different rs values
    rs = 1.0
    energy = jellium_energy_per_electron(rs)
    # For rs=1, energy should be higher (more kinetic, less exchange/correlation)
    @test energy > 0.0
end

@testset "SCF Loop" begin
    # Create a simple uniform electron gas system
    a = 5.0  # Bohr
    n_electrons = 2
    cutoff = 10.0  # Hartree
    fft_size = (8, 8, 8)
    
    system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
    
    # Initialize density
    initialize_uniform_density(system)
    
    # Check that density is uniform
    expected_density = n_electrons / (a^3)
    @test isapprox(system.density.data, fill(expected_density, size(system.density.data));)
    
    # Run SCF with few iterations
    params = SCFParameters(max_iter=5, energy_tolerance=1e-4, density_tolerance=1e-4)
    
    # Run SCF
    converged_system = self_consistent_field(system, params)
    
    # Check that energy is computed
    @test converged_system.energies.total < 0.0  # Should be negative for bound system
    
    # For uniform electron gas, the density should remain uniform
    @test isapprox(converged_system.density.data, fill(expected_density, size(converged_system.density.data));)
end

@testset "System Creation" begin
    # Test create_uniform_electron_gas
    a = 5.0
    n_electrons = 4
    cutoff = 10.0
    fft_size = (16, 16, 16)
    
    system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)
    
    @test system.lattice.volume ≈ a^3
    @test system.electrons == n_electrons
    @test system.basis.cutoff ≈ cutoff
    @test system.basis.fft_size == fft_size
    
    # Test jellium_total_energy
    total_energy = jellium_total_energy(system)
    @test total_energy < 0.0  # Should be negative for bound system
end

# Print test summary
println("\nPhase 1 Tests Complete!")
println("All core functionality for uniform electron gas is working.")
