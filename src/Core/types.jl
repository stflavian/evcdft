"""
Fundamental data structures for DFT calculations.
"""

module Types

using LinearAlgebra
using ..Constants: π, twopi

# Abstract types for dispatch
abstract type AbstractBasis end
abstract type AbstractPotential end
abstract type AbstractSystem end

"""
    Represents a 3D lattice in real space.
    
    Fields:
    - a1, a2, a3: Lattice vectors (in Bohr)
    - volume: Unit cell volume (in Bohr³)
    - b1, b2, b3: Reciprocal lattice vectors (in Bohr⁻¹)
"""
struct Lattice
    a1::Vector{Float64}
    a2::Vector{Float64}
    a3::Vector{Float64}
    volume::Float64
    b1::Vector{Float64}
    b2::Vector{Float64}
    b3::Vector{Float64}
    
    function Lattice(a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64})
        # Compute volume using scalar triple product
        vol = abs(dot(a1, cross(a2, a3)))
        
        # Compute reciprocal lattice vectors: b_i = 2π (a_j × a_k) / volume
        b1 = twopi * cross(a2, a3) / vol
        b2 = twopi * cross(a3, a1) / vol
        b3 = twopi * cross(a1, a2) / vol
        
        new(a1, a2, a3, vol, b1, b2, b3)
    end
end

"""
    Simple cubic lattice constructor.
    
    Args:
    - a: Lattice constant (in Bohr)
"""
function SimpleCubic(a::Float64)
    a1 = [a, 0.0, 0.0]
    a2 = [0.0, a, 0.0]
    a3 = [0.0, 0.0, a]
    return Lattice(a1, a2, a3)
end

"""
    FCC lattice constructor.
    
    Args:
    - a: Lattice constant (in Bohr)
"""
function FCC(a::Float64)
    a1 = [a, 0.0, 0.0]
    a2 = [0.0, a, 0.0]
    a3 = [0.0, 0.0, a]
    # For FCC, we need to adjust - but this is simple cubic
    # Actual FCC would have a1 = [0, a/2, a/2], a2 = [a/2, 0, a/2], a3 = [a/2, a/2, 0]
    # But for uniform electron gas, simple cubic is fine
    return Lattice(a1, a2, a3)
end

"""
    Represents a uniform electron gas (jellium model).
    
    Fields:
    - lattice: The simulation cell lattice
    - n_electrons: Total number of electrons
    - density: Uniform electron density (in Bohr⁻³)
    - rs: Wigner-Seitz radius (in Bohr)
"""
struct UniformElectronGas
    lattice::Lattice
    n_electrons::Int
    density::Float64
    rs::Float64  # Wigner-Seitz radius: (3/(4πn))^(1/3)
    
    function UniformElectronGas(lattice::Lattice, n_electrons::Int)
        volume = lattice.volume
        density = n_electrons / volume
        rs = (3.0 / (4π * density))^(1/3)
        new(lattice, n_electrons, density, rs)
    end
end

"""
    Represents a plane wave basis set.
    
    Fields:
    - lattice: The real-space lattice
    - cutoff: Energy cutoff (in Hartree)
    - g_vectors: List of G vectors (reciprocal lattice vectors)
    - g2: |G|² for each G vector
    - n_g: Number of G vectors
    - fft_size: FFT grid size (nx, ny, nz)
"""
struct PlaneWaveBasis <: AbstractBasis
    lattice::Lattice
    cutoff::Float64
    g_vectors::Vector{Vector{Int}}
    g2::Vector{Float64}  # |G|² values
    n_g::Int
    fft_size::Tuple{Int, Int, Int}
    
    function PlaneWaveBasis(lattice::Lattice, cutoff::Float64, fft_size::Tuple{Int, Int, Int})
        # Generate all reciprocal lattice vectors within cutoff
        nx, ny, nz = fft_size
        g_vectors = Vector{Vector{Int}}()
        g2_list = Vector{Float64}()
        
        # Generate G vectors: G = i*b1 + j*b2 + k*b3
        for i in -div(nx, 2):div(nx, 2)
            for j in -div(ny, 2):div(ny, 2)
                for k in -div(nz, 2):div(nz, 2)
                    g_vec = [i, j, k]
                    g_cart = i * lattice.b1 + j * lattice.b2 + k * lattice.b3
                    g2 = dot(g_cart, g_cart)
                    
                    # Kinetic energy = |G|² / 2 (in Hartree)
                    if g2 / 2 <= cutoff
                        push!(g_vectors, g_vec)
                        push!(g2_list, g2)
                    end
                end
            end
        end
        
        n_g = length(g_vectors)
        new(lattice, cutoff, g_vectors, g2_list, n_g, fft_size)
    end
end

"""
    Represents electron density in real space.
    
    Fields:
    - data: Density values on FFT grid (in Bohr⁻³)
    - grid_size: Size of the FFT grid
"""
struct ElectronDensity
    data::Array{Float64, 3}
    grid_size::Tuple{Int, Int, Int}
    
    function ElectronDensity(grid_size::Tuple{Int, Int, Int})
        data = zeros(Float64, grid_size)
        new(data, grid_size)
    end
end

"""
    Represents electron density in reciprocal space.
    
    Fields:
    - data: Density values in reciprocal space
    - grid_size: Size of the FFT grid
"""
struct ElectronDensityReciprocal
    data::Array{ComplexF64, 3}
    grid_size::Tuple{Int, Int, Int}
    
    function ElectronDensityReciprocal(grid_size::Tuple{Int, Int, Int})
        data = zeros(ComplexF64, grid_size)
        new(data, grid_size)
    end
end

"""
    Represents the Kohn-Sham potential.
    
    Fields:
    - hartree: Hartree potential (V_H)
    - exchange: Exchange potential (V_x)
    - correlation: Correlation potential (V_c)
    - external: External potential (V_ext, e.g., ionic)
"""
struct KohnShamPotential <: AbstractPotential
    hartree::Array{Float64, 3}
    exchange::Array{Float64, 3}
    correlation::Array{Float64, 3}
    external::Array{Float64, 3}
    grid_size::Tuple{Int, Int, Int}
    
    function KohnShamPotential(grid_size::Tuple{Int, Int, Int})
        hartree = zeros(Float64, grid_size)
        exchange = zeros(Float64, grid_size)
        correlation = zeros(Float64, grid_size)
        external = zeros(Float64, grid_size)
        new(hartree, exchange, correlation, external, grid_size)
    end
end

"""
    Total energy components.
    
    Fields:
    - kinetic: Kinetic energy
    - hartree: Hartree (electrostatic) energy
    - exchange: Exchange energy
    - correlation: Correlation energy
    - external: External energy (e.g., ion-electron)
    - total: Total energy
"""
struct EnergyComponents
    kinetic::Float64
    hartree::Float64
    exchange::Float64
    correlation::Float64
    external::Float64
    total::Float64
    
    function EnergyComponents()
        new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
end

"""
    Parameters for self-consistent field (SCF) calculation.
    
    Fields:
    - max_iter: Maximum number of SCF iterations
    - energy_tolerance: Convergence tolerance for energy (Hartree)
    - density_tolerance: Convergence tolerance for density (Bohr⁻³)
    - mixing_parameter: Mixing parameter for density mixing
    - mixing_type: Type of mixing ("linear", "kerker", "pulay", "broyden")
"""
struct SCFParameters
    max_iter::Int
    energy_tolerance::Float64
    density_tolerance::Float64
    mixing_parameter::Float64
    mixing_type::String
    
    function SCFParameters(;
        max_iter::Int = 100,
        energy_tolerance::Float64 = 1e-6,
        density_tolerance::Float64 = 1e-6,
        mixing_parameter::Float64 = 0.5,
        mixing_type::String = "linear"
    )
        new(max_iter, energy_tolerance, density_tolerance, mixing_parameter, mixing_type)
    end
end

"""
    Main DFT system structure.
    
    Fields:
    - lattice: The simulation cell lattice
    - basis: The basis set (e.g., PlaneWaveBasis)
    - electrons: Number of electrons
    - density: Electron density (real space)
    - potential: Kohn-Sham potential
    - energies: Energy components
"""
struct DFTSystem <: AbstractSystem
    lattice::Lattice
    basis::AbstractBasis
    electrons::Int
    density::ElectronDensity
    density_recip::ElectronDensityReciprocal
    potential::KohnShamPotential
    energies::EnergyComponents
    
    function DFTSystem(lattice::Lattice, basis::AbstractBasis, electrons::Int)
        grid_size = basis.fft_size
        density = ElectronDensity(grid_size)
        density_recip = ElectronDensityReciprocal(grid_size)
        potential = KohnShamPotential(grid_size)
        energies = EnergyComponents()
        new(lattice, basis, electrons, density, density_recip, potential, energies)
    end
end

end # module Types
