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

# Include IO modules
include("IO/hsd_parser.jl")

# Import types needed for helper functions
# These need to be defined before including input_parser.jl and cli.jl
using .Types: Lattice, SimpleCubic, PlaneWaveBasis, DFTSystem, ElectronDensity, 
               ElectronDensityReciprocal, KohnShamPotential, EnergyComponents, AtomicSystem
using .PlaneWave: PlaneWaveBasis

# Define helper functions that are used by submodules
# These need to be defined before including input_parser.jl and cli.jl

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
    Create an atomic system from lattice, species, and coordinates.
    
    Args:
    - lattice_constant_or_vectors: Either a lattice constant (Float64) for simple cubic,
      or a 3x3 matrix of lattice vectors
    - species_names: Vector of element symbols (e.g., ["Ga", "As"])
    - types_and_coords: N x 4 matrix where each row is [type_index, x, y, z]
      Coordinates are assumed to be in Angstrom and will be converted to Bohr
    - net_charge: Net charge of the system (default: 0)
    - periodic: Whether the system has periodic boundary conditions (default: true)
    - cutoff: Energy cutoff for plane wave basis (default: 10.0 Ha)
    - fft_size: FFT grid size (default: (16, 16, 16))
    
    Returns:
    - system: DFTSystem for atomic system
"""
function create_atomic_system(lattice_constant_or_vectors::Union{Float64, Matrix{Float64}},
                            species_names::Vector{String},
                            types_and_coords::Matrix{Float64};
                            net_charge::Int=0,
                            periodic::Bool=true,
                            cutoff::Float64=10.0,
                            fft_size::Tuple{Int,Int,Int}=(16,16,16))
    # If lattice_constant_or_vectors is a Float64, create simple cubic lattice
    if lattice_constant_or_vectors isa Float64
        lattice = SimpleCubic(lattice_constant_or_vectors)
    else
        # Assume it's a 3x3 matrix
        a1 = lattice_constant_or_vectors[:, 1]
        a2 = lattice_constant_or_vectors[:, 2]
        a3 = lattice_constant_or_vectors[:, 3]
        lattice = Lattice(a1, a2, a3)
    end
    
    # Convert types_and_coords to Bohr if needed (assume Angstrom input)
    types_and_coords_bohr = copy(types_and_coords)
    types_and_coords_bohr[:, 2:4] .*= angstrom_to_bohr
    
    # Create atomic system
    atomic_system = AtomicSystem(lattice, species_names, types_and_coords_bohr,
                               net_charge=net_charge, periodic=periodic)
    
    # Create basis
    basis = PlaneWaveBasis(lattice, cutoff, fft_size)
    
    # Create DFT system
    return DFTSystem(lattice, basis, atomic_system)
end

"""
    Calculate the total energy for a uniform electron gas.
    
    For a uniform electron gas (jellium), the total energy per electron is:
    E/N = (2.8376/rs^2) - (0.9163/rs) + \u03b5_c(rs)
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
    correlation = lda_correlation_energy(3.0 / (4\u03c0 * rs^3))
    
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
        rs = (3.0 / (4\u03c0 * (system.electrons / volume)))^(1/3)
        return system.electrons * jellium_energy_per_electron(rs)
    else
        error("System is not a uniform electron gas")
    end
end

# Now include modules that depend on the above functions
include("IO/input_parser.jl")
include("cli.jl")

# Export commonly used types and functions
using .Constants
using .Units
using .Types
using .Types: Lattice, SimpleCubic, FCC, UniformElectronGas, PlaneWaveBasis, ElectronDensity, ElectronDensityReciprocal, KohnShamPotential, EnergyComponents, SCFParameters, DFTSystem, AtomicSpecies, Atom, AtomicSystem, Wavefunction, KohnShamStates
using .PlaneWave
using .XCFunctionals
using .SelfConsistent
using .HSDParser
using .InputParser
using .CLI

export Lattice, SimpleCubic, FCC

export UniformElectronGas

export PlaneWaveBasis

export ElectronDensity, ElectronDensityReciprocal

export KohnShamPotential

export EnergyComponents

export SCFParameters

export DFTSystem

# NEW EXPORTS for Phase 2A
export AtomicSpecies, Atom, AtomicSystem
export Wavefunction, KohnShamStates
export cartesian_to_fractional, fractional_to_cartesian
export total_nuclear_charge, total_electrons
export get_atomic_number, ATOMIC_NUMBERS
export create_atomic_system

export lda_exchange_energy, lda_correlation_energy, lda_xc_energy

export lda_exchange_potential, lda_correlation_potential, lda_xc_potential

export compute_lda_xc

export compute_hartree_potential, compute_hartree_energy

export run_scf!, self_consistent_field

export initialize_uniform_density
export create_uniform_electron_gas, jellium_energy_per_electron, jellium_total_energy

# Export IO functions
export HSDNode, parse_hsd_file, parse_hsd_string
export InputConfig, parse_input_file, run_from_input, run_from_input_string

# Export CLI function
export main

end # module EVC_DFT
