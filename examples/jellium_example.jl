"""
Example: Uniform Electron Gas (Jellium) Calculation

This example demonstrates the basic usage of EVC_DFT for a uniform electron gas
(jellium model) calculation. The jellium model is a simple but important test case
for DFT implementations.

The uniform electron gas consists of electrons moving in a uniform positive
background charge, which makes the system overall neutral.
"""

using EVC_DFT
using Printf

println("="^60)
println("EVC_DFT - Uniform Electron Gas Example")
println("="^60)

# Parameters for the calculation
const a = 5.0       # Lattice constant (Bohr)
const n_electrons = 2  # Number of electrons
const cutoff = 10.0  # Energy cutoff for plane wave basis (Hartree)
const fft_size = (16, 16, 16)  # FFT grid size

println("\nCalculation Parameters:")
println("  Lattice constant (a): $a Bohr")
println("  Number of electrons: $n_electrons")
println("  Energy cutoff: $cutoff Hartree")
println("  FFT grid size: $fft_size")

# Create the uniform electron gas system
println("\nCreating uniform electron gas system...")
system = create_uniform_electron_gas(a, n_electrons, cutoff, fft_size)

# Calculate Wigner-Seitz radius
volume = system.lattice.volume
rs = (3.0 / (4π * (n_electrons / volume)))^(1/3)
println("  Wigner-Seitz radius (rs): $rs Bohr")

# Calculate expected energy per electron using analytical formula
energy_per_electron = jellium_energy_per_electron(rs)
println("  Expected energy per electron: $energy_per_electron Hartree")

# Expected total energy
expected_total_energy = n_electrons * energy_per_electron
println("  Expected total energy: $expected_total_energy Hartree")

# Initialize the electron density
println("\nInitializing electron density...")
initialize_uniform_density(system)

# Check initial density
expected_density = n_electrons / volume
@printf("  Initial density: %.6f electrons/Bohr³\n", expected_density)
@printf("  Max density: %.6f electrons/Bohr³\n", maximum(system.density.data))
@printf("  Min density: %.6f electrons/Bohr³\n", minimum(system.density.data))

# Set up SCF parameters
params = SCFParameters(
    max_iter=50,
    energy_tolerance=1e-8,
    density_tolerance=1e-8,
    mixing_parameter=0.5,
    mixing_type="linear"
)

println("\nRunning self-consistent field calculation...")
println("  Maximum iterations: $(params.max_iter)")
println("  Energy tolerance: $(params.energy_tolerance) Hartree")
println("  Density tolerance: $(params.density_tolerance) electrons/Bohr³")
println("  Mixing: $(params.mixing_type) with α=$(params.mixing_parameter)")

# Run the SCF calculation
converged_system = self_consistent_field(system, params)

# Print results
println("\n" * "-"^60)
println("Results:")
println("-"^60)

@printf("Total energy: %.8f Hartree\n", converged_system.energies.total)
@printf("Total energy: %.4f eV\n", converged_system.energies.total * hartree_to_ev)
@printf("Hartree energy: %.8f Hartree\n", converged_system.energies.hartree)
@printf("Exchange energy: %.8f Hartree\n", converged_system.energies.exchange)
@printf("Correlation energy: %.8f Hartree\n", converged_system.energies.correlation)

@printf("\nFinal density:\n")
@printf("  Max: %.6f electrons/Bohr³\n", maximum(converged_system.density.data))
@printf("  Min: %.6f electrons/Bohr³\n", minimum(converged_system.density.data))
@printf("  Mean: %.6f electrons/Bohr³\n", mean(converged_system.density.data))

# Compare with analytical result
error = abs(converged_system.energies.total - expected_total_energy)
@printf("\nComparison with analytical result:\n")
@printf("  Analytical total energy: %.8f Hartree\n", expected_total_energy)
@printf("  Numerical total energy: %.8f Hartree\n", converged_system.energies.total)
@printf("  Absolute error: %.8f Hartree\n", error)
@printf("  Relative error: %.4f%%\n", 100 * error / abs(expected_total_energy))

# Additional analysis for jellium
println("\n" * "-"^60)
println("Jellium Analysis:")
println("-"^60)

# Kinetic energy (free electron gas)
kinetic_per_electron = 2.8376 / (rs^2)
@printf("Kinetic energy per electron: %.6f Hartree\n", kinetic_per_electron)

# Exchange energy (Dirac)
exchange_per_electron = -0.9163 / rs
@printf("Exchange energy per electron: %.6f Hartree\n", exchange_per_electron)

# Correlation energy
correlation_per_electron = lda_correlation_energy(n_electrons / volume)
@printf("Correlation energy per electron: %.6f Hartree\n", correlation_per_electron)

# Sum
total_per_electron = kinetic_per_electron + exchange_per_electron + correlation_per_electron
@printf("Total energy per electron: %.6f Hartree\n", total_per_electron)

println("\n" * "="^60)
println("Example complete!")
println("="^60)
