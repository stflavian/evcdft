"""
Plane wave basis set implementation for DFT calculations.
"""

module PlaneWave

using LinearAlgebra
using FFTW
using ..Types: Lattice, PlaneWaveBasis, ElectronDensity, ElectronDensityReciprocal
using ..Constants: twopi, π

import ..Types: AbstractBasis

# Export the main type and functions
export PlaneWaveBasis, compute_kinetic_energy, apply_laplacian, 
       fft_forward, fft_backward, generate_g_vectors

"""
    Generate all reciprocal lattice vectors within a given cutoff.
    
    Args:
    - lattice: The real-space lattice
    - cutoff: Energy cutoff (in Hartree). Kinetic energy = |G|² / 2
    - fft_size: Size of the FFT grid (nx, ny, nz)
    
    Returns:
    - g_vectors: List of G vector indices (i, j, k)
    - g2: |G|² for each G vector
    - g_cart: Cartesian coordinates of G vectors
"""
function generate_g_vectors(lattice::Lattice, cutoff::Float64, fft_size::Tuple{Int, Int, Int})
    nx, ny, nz = fft_size
    g_vectors = Vector{Vector{Int}}()
    g2_list = Vector{Float64}()
    g_cart_list = Vector{Vector{Float64}}()
    
    # Generate G vectors: G = i*b1 + j*b2 + k*b3
    # Note: FFTW uses a specific convention for negative frequencies
    # For even sizes, the negative frequencies are in the second half
    for i in 0:nx-1
        for j in 0:ny-1
            for k in 0:nz-1
                # Get the actual G vector index (FFTW convention)
                gi = i < nx ÷ 2 + 1 ? i : i - nx
                gj = j < ny ÷ 2 + 1 ? j : j - ny
                gk = k < nz ÷ 2 + 1 ? k : k - nz
                
                g_cart = gi * lattice.b1 + gj * lattice.b2 + gk * lattice.b3
                g2 = dot(g_cart, g_cart)
                
                # Kinetic energy = |G|² / 2 (in Hartree)
                if g2 / 2 <= cutoff
                    push!(g_vectors, [gi, gj, gk])
                    push!(g2_list, g2)
                    push!(g_cart_list, g_cart)
                end
            end
        end
    end
    
    return g_vectors, g2_list, g_cart_list
end

"""
    Compute the kinetic energy for each G vector.
    
    Args:
    - basis: PlaneWaveBasis
    
    Returns:
    - kinetic_energy: Kinetic energy (|G|² / 2) for each G vector
"""
function compute_kinetic_energy(basis::PlaneWaveBasis)
    return basis.g2 / 2.0
end

"""
    Apply the Laplacian operator in reciprocal space.
    
    The Laplacian in reciprocal space is -|G|².
    
    Args:
    - f_recip: Function in reciprocal space
    - basis: PlaneWaveBasis
    
    Returns:
    - laplacian: Laplacian of f in reciprocal space
"""
function apply_laplacian(f_recip::Array{ComplexF64, 3}, basis::PlaneWaveBasis)
    nx, ny, nz = basis.fft_size
    laplacian = similar(f_recip)
    
    # Compute G² for all points in the FFT grid
    for i in 1:nx, j in 1:ny, k in 1:nz
        # Get G vector index
        gi = i <= nx ÷ 2 + 1 ? i - 1 : i - 1 - nx
        gj = j <= ny ÷ 2 + 1 ? j - 1 : j - 1 - ny
        gk = k <= nz ÷ 2 + 1 ? k - 1 : k - 1 - nz
        
        g_cart = gi * basis.lattice.b1 + gj * basis.lattice.b2 + gk * basis.lattice.b3
        g2 = dot(g_cart, g_cart)
        
        laplacian[i, j, k] = -g2 * f_recip[i, j, k]
    end
    
    return laplacian
end

"""
    Forward FFT from real space to reciprocal space.
    
    Args:
    - f_real: Function in real space
    
    Returns:
    - f_recip: Function in reciprocal space
"""
function fft_forward(f_real::Array{Float64, 3})
    nx, ny, nz = size(f_real)
    f_recip = zeros(ComplexF64, nx, ny, nz)
    
    # Use FFTW for efficient FFT
    plan = FFTW.plan_fft(f_real)
    f_recip .= plan * f_real
    
    return f_recip
end

"""
    Backward FFT from reciprocal space to real space.
    
    Args:
    - f_recip: Function in reciprocal space
    
    Returns:
    - f_real: Function in real space
"""
function fft_backward(f_recip::Array{ComplexF64, 3})
    nx, ny, nz = size(f_recip)
    f_real = zeros(Float64, nx, ny, nz)
    
    # Use FFTW for efficient IFFT
    plan = FFTW.plan_ifft(f_recip)
    f_real .= real.(plan * f_recip)
    
    return f_real
end

"""
    Compute the Hartree potential from electron density using FFT.
    
    The Hartree potential satisfies Poisson's equation: ∇²V_H = -4πρ
    In reciprocal space: V_H(G) = -4πρ(G) / |G|²
    
    Args:
    - density: Electron density in real space
    - basis: PlaneWaveBasis
    
    Returns:
    - potential: Hartree potential in real space
"""
function compute_hartree_potential(density::ElectronDensity, basis::PlaneWaveBasis)
    nx, ny, nz = basis.fft_size
    lattice = basis.lattice
    
    # Forward FFT to get density in reciprocal space
    density_recip = fft_forward(density.data)
    
    # Compute Hartree potential in reciprocal space
    potential_recip = zeros(ComplexF64, nx, ny, nz)
    
    for i in 1:nx, j in 1:ny, k in 1:nz
        # Get G vector index
        gi = i <= nx ÷ 2 + 1 ? i - 1 : i - 1 - nx
        gj = j <= ny ÷ 2 + 1 ? j - 1 : j - 1 - ny
        gk = k <= nz ÷ 2 + 1 ? k - 1 : k - 1 - nz
        
        g_cart = gi * lattice.b1 + gj * lattice.b2 + gk * lattice.b3
        g2 = dot(g_cart, g_cart)
        
        # Handle G=0 term (uniform background charge)
        # For a neutral system, ρ(G=0) = 0, so V_H(G=0) = 0
        if g2 > 1e-12  # Avoid division by zero
            potential_recip[i, j, k] = -4π * density_recip[i, j, k] / g2
        else
            potential_recip[i, j, k] = 0.0
        end
    end
    
    # Backward FFT to get potential in real space
    potential_real = fft_backward(potential_recip)
    
    return potential_real
end

"""
    Compute the kinetic energy density from wavefunctions.
    
    For a single wavefunction ψ:
    Kinetic energy = (1/2) ∫ ψ* (-∇²) ψ d³r
    
    In reciprocal space: -∇²ψ(G) = |G|² ψ(G)
    
    Args:
    - wavefunction: Wavefunction in reciprocal space
    - basis: PlaneWaveBasis
    
    Returns:
    - kinetic_energy: Kinetic energy (Hartree)
"""
function compute_kinetic_energy(wavefunction::Array{ComplexF64, 3}, basis::PlaneWaveBasis)
    nx, ny, nz = size(wavefunction)
    lattice = basis.lattice
    volume = lattice.volume
    
    kinetic_energy = 0.0
    
    for i in 1:nx, j in 1:ny, k in 1:nz
        # Get G vector index
        gi = i <= nx ÷ 2 + 1 ? i - 1 : i - 1 - nx
        gj = j <= ny ÷ 2 + 1 ? j - 1 : j - 1 - ny
        gk = k <= nz ÷ 2 + 1 ? k - 1 : k - 1 - nz
        
        g_cart = gi * lattice.b1 + gj * lattice.b2 + gk * lattice.b3
        g2 = dot(g_cart, g_cart)
        
        # Kinetic energy contribution: (1/2) |G|² |ψ(G)|²
        kinetic_energy += 0.5 * g2 * abs2(wavefunction[i, j, k])
    end
    
    # Normalize by volume
    return kinetic_energy / volume
end

"""
    Compute the Hartree energy: (1/2) ∫ ρ V_H d³r
    
    Args:
    - density: Electron density in real space
    - potential: Hartree potential in real space
    - volume: Cell volume
    
    Returns:
    - hartree_energy: Hartree energy (Hartree)
"""
function compute_hartree_energy(density::Array{Float64, 3}, potential::Array{Float64, 3}, volume::Float64)
    # Hartree energy = (1/2) ∫ ρ V_H d³r
    hartree_energy = 0.5 * sum(density .* potential) * (volume / length(density))
    return hartree_energy
end

end # module PlaneWave
