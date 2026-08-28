"""
Simple test script for EVC_DFT package.
This can be run directly to verify the package works.
"""

# Include the main module
using Test
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "EVC_DFT.jl"))

using .EVC_DFT
using .EVC_DFT.Constants
using .EVC_DFT.Types
using .EVC_DFT.PlaneWave
using .EVC_DFT.XCFunctionals
using .EVC_DFT.SelfConsistent

println("Testing EVC_DFT package...")

# Test 1: Basic module loading
println("✓ Module loaded successfully")

# Test 2: Constants
println("Testing constants...")
@assert hbar == 1.0
@assert m_e == 1.0
@assert e == 1.0
@assert pi ≈ 3.141592653589793
println("✓ Constants work")

# Test 3: Types
println("Testing types...")
lattice = SimpleCubic(5.0)
@assert lattice.volume ≈ 125.0

basis = PlaneWaveBasis(lattice, 10.0, (8, 8, 8))
@assert length(basis.g_vectors) > 0

system = create_uniform_electron_gas(5.0, 2, 10.0, (8, 8, 8))
@assert system.electrons == 2
println("✓ Types work")

# Test 4: XC Functionals
println("Testing XC functionals...")
density = 0.1
@assert lda_exchange_energy(density) < 0.0
@assert lda_correlation_energy(density) < 0.0
@assert lda_xc_energy(density) < 0.0
println("✓ XC functionals work")

# Test 5: SCF
println("Testing SCF...")
system = create_uniform_electron_gas(5.0, 2, 10.0, (8, 8, 8))
initialize_uniform_density(system)
params = SCFParameters(max_iter=5, energy_tolerance=1e-4, density_tolerance=1e-4)
converged_system = run_scf!(system, params)
@assert isfinite(converged_system.energies.total)
@assert converged_system.energies.total < 0.0
println("✓ SCF works")

# Test 6: Jellium energy
println("Testing jellium energy...")
rs = 2.0
energy = jellium_energy_per_electron(rs)
@assert isfinite(energy)
@assert energy ≈ 0.111
println("✓ Jellium energy works")

println("\n=== All tests passed! ===")
println("EVC_DFT package is working correctly.")
