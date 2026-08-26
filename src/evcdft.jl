"""
EVC_DFT - A Density Functional Theory calculator in Julia.

This module provides a basic implementation of DFT for uniform electron gas
and periodic systems using plane wave basis sets.
"""

module EVC_DFT

# Include core modules
include("Core/constants.jl")
include("Core/units.jl")
include("Core/types.jl")

# Include basis modules
include("Basis/plane_wave.jl")

# Include potential modules
include("Potential/xc_functionals.jl")

# Include SCF modules
include("SCF/self_consistent.jl")

# Export commonly used types and functions
using .Constants
using .Units
using .Types
using .PlaneWave
using .XCFunctionals
using .SelfConsistent

export Lattice, SimpleCubic, FCC

export UniformElectronGas

export PlaneWaveBasis

export ElectronDensity, ElectronDensityReciprocal

export KohnShamPotential

export EnergyComponents

export SCFParameters

export DFTSystem

export lda_exchange_energy, lda_correlation_energy, lda_xc_energy

export lda_exchange_potential, lda_correlation_potential, lda_xc_potential

export compute_lda_xc

export compute_hartree_potential, compute_hartree_energy

export run_scf!, self_consistent_field

export initialize_uniform_density

"""
    Create a simple uniform electron gas system for testing.
    
    Args:
    - a: Lattice constant (in Bohr)
    - n_electrons: Number of electrons
    - cutoff: Energy cutoff for plane wave basis (Hartree)
    - fft_size: FFT grid size (nx, ny, nz)
    
    Returns:
    - system: DFTSystem for uniform electron gas
"""
function create_uniform_electron_gas(a::Float64, n_electrons::Int, 
                                     cutoff::Float64, fft_size::Tuple{Int, Int, Int})
    # Create lattice
    lattice = SimpleCubic(a)
    
    # Create plane wave basis
    basis = PlaneWaveBasis(lattice, cutoff, fft_size)
    
    # Create DFT system
    system = DFTSystem(lattice, basis, n_electrons)
    
    return system
end

"""
    Calculate the total energy for a uniform electron gas.
    
    For a uniform electron gas (jellium), the total energy per electron is:
    E/N = (2.8376/rs²) - (0.9163/rs) + ε_c(rs)
    where rs is the Wigner-Seitz radius.
    
    Args:
    - rs: Wigner-Seitz radius (in Bohr)
    
    Returns:
    - energy_per_electron: Total energy per electron (Hartree)
"""
function jellium_energy_per_electron(rs::Float64)
    # Kinetic energy (free electron gas)
    kinetic = 2.8376 / (rs^2)
    
    # Exchange energy (Dirac)
    exchange = -0.9163 / rs
    
    # Correlation energy (Perdew-Zunger)
    correlation = lda_correlation_energy(3.0 / (4π * rs^3))
    
    return kinetic + exchange + correlation
end

"""
    Calculate the total energy for a uniform electron gas system.
    
    Args:
    - system: DFTSystem for uniform electron gas
    
    Returns:
    - total_energy: Total energy (Hartree)
"""
function jellium_total_energy(system::DFTSystem)
    # Get Wigner-Seitz radius from the uniform electron gas
    if system.basis isa PlaneWaveBasis
        # For uniform electron gas
        volume = system.lattice.volume
        rs = (3.0 / (4π * (system.electrons / volume)))^(1/3)
        return system.electrons * jellium_energy_per_electron(rs)
    else
        error("System is not a uniform electron gas")
    end
end

end # module EVC_DFT
