"""
Self-consistent field (SCF) loop for DFT calculations.

This module implements the basic SCF loop for solving the Kohn-Sham equations.
"""

module SelfConsistent

using LinearAlgebra
using FFTW
using ..Types: DFTSystem, Lattice, PlaneWaveBasis, ElectronDensity, 
                ElectronDensityReciprocal, KohnShamPotential, EnergyComponents, SCFParameters
using ..PlaneWave: compute_hartree_potential, compute_hartree_energy, fft_forward, fft_backward
using ..XCFunctionals: compute_lda_xc

# Export the main function
export run_scf!, self_consistent_field, initialize_uniform_density

"""
    Mix two densities using linear mixing.
    
    Args:
    - density_new: New density
    - density_old: Old density
    - alpha: Mixing parameter (0 < alpha < 1)
    
    Returns:
    - density_mixed: Mixed density
"""
function linear_mixing(density_new::Array{Float64, 3}, density_old::Array{Float64, 3}, alpha::Float64)
    return alpha * density_new + (1.0 - alpha) * density_old
end

"""
    Simple density mixing (Kerker-style).
    
    Args:
    - density_new: New density
    - density_old: Old density
    - alpha: Mixing parameter
    
    Returns:
    - density_mixed: Mixed density
"""
function kerker_mixing(density_new::Array{Float64, 3}, density_old::Array{Float64, 3}, alpha::Float64)
    # For now, use linear mixing
    # Kerker mixing would involve screening in reciprocal space
    return linear_mixing(density_new, density_old, alpha)
end

"""
    Apply density mixing based on the mixing type.
    
    Args:
    - density_new: New density
    - density_old: Old density
    - params: SCFParameters
    
    Returns:
    - density_mixed: Mixed density
"""
function apply_mixing(density_new::Array{Float64, 3}, density_old::Array{Float64, 3}, params::SCFParameters)
    if params.mixing_type == "linear"
        return linear_mixing(density_new, density_old, params.mixing_parameter)
    elseif params.mixing_type == "kerker"
        return kerker_mixing(density_new, density_old, params.mixing_parameter)
    else
        # Default to linear mixing
        return linear_mixing(density_new, density_old, params.mixing_parameter)
    end
end

"""
    Compute the total electron density from Kohn-Sham orbitals.
    
    For a uniform electron gas or single orbital:
    ρ(r) = Σ |ψ_i(r)|² * f_i
    where f_i is the occupation number.
    
    Args:
    - wavefunctions: Array of wavefunctions in real space
    - occupations: Array of occupation numbers
    
    Returns:
    - density: Electron density in real space
"""
function compute_density(wavefunctions::Vector{Array{Float64, 3}}, occupations::Vector{Float64})
    nx, ny, nz = size(wavefunctions[1])
    density = zeros(Float64, nx, ny, nz)
    
    for (i, (psi, occ)) in enumerate(zip(wavefunctions, occupations))
        density .+= occ .* abs2.(psi)
    end
    
    return density
end

"""
    Initialize the electron density (e.g., uniform density for jellium).
    
    Args:
    - system: DFTSystem
    
    Returns:
    - density: Initial electron density
"""
function initialize_density(system::DFTSystem)
    nx, ny, nz = system.density.grid_size
    volume = system.lattice.volume
    n_electrons = system.electrons
    
    # Uniform density for jellium
    uniform_density = n_electrons / volume
    
    # Create uniform density array
    density_data = fill(uniform_density, nx, ny, nz)
    
    return ElectronDensity((nx, ny, nz), density_data)
end

"""
    Initialize the electron density for a uniform electron gas.
    
    Args:
    - system: DFTSystem
    
    Returns:
    - density: Uniform electron density
"""
function initialize_uniform_density(system::DFTSystem)
    nx, ny, nz = system.density.grid_size
    volume = system.lattice.volume
    n_electrons = system.electrons
    
    # Uniform density
    uniform_density = n_electrons / volume
    
    # Fill the density array
    system.density.data .= uniform_density
    
    return system.density
end

"""
    Compute the total energy from the current density and potential.
    
    Args:
    - system: DFTSystem
    
    Returns:
    - total_energy: Total energy (Hartree)
"""
function compute_total_energy(system::DFTSystem)
    volume = system.lattice.volume
    
    # For now, compute only Hartree and XC energies
    # Kinetic energy would require wavefunctions
    
    # Hartree energy
    hartree_potential = compute_hartree_potential(system.density, system.basis)
    hartree_energy = compute_hartree_energy(system.density.data, hartree_potential, volume)
    
    # XC energy
    xc_energy, xc_potential = compute_lda_xc(system.density, volume)
    
    # Total energy (simplified - missing kinetic energy)
    total_energy = hartree_energy + xc_energy
    
    return total_energy
end

"""
    Check if the SCF calculation has converged.
    
    Args:
    - system: DFTSystem
    - params: SCFParameters
    - old_energy: Energy from previous iteration
    - old_density: Density from previous iteration
    
    Returns:
    - is_converged: Boolean indicating convergence
"""
function check_convergence(system::DFTSystem, params::SCFParameters, 
                          old_energy::Float64, old_density::Array{Float64, 3})
    # Check energy convergence
    energy_diff = abs(system.energies.total - old_energy)
    if energy_diff > params.energy_tolerance
        return false
    end
    
    # Check density convergence
    density_diff = maximum(abs.(system.density.data - old_density))
    if density_diff > params.density_tolerance
        return false
    end
    
    return true
end

"""
    Run a single SCF iteration.
    
    For a uniform electron gas (jellium model), the density remains uniform.
    In this simplified implementation, we compute the potentials and energy
    but keep the density uniform, which is the correct solution for jellium.
    
    Args:
    - system: DFTSystem
    - params: SCFParameters
    
    Returns:
    - new_energy: Total energy after this iteration
    - new_density: New electron density (for jellium, same as input)
"""
function scf_iteration!(system::DFTSystem, params::SCFParameters)
    volume = system.lattice.volume
    n_electrons = system.electrons
    
    # Compute Hartree potential
    hartree_potential = compute_hartree_potential(system.density, system.basis)
    system.potential.hartree .= hartree_potential
    
    # Compute XC potential (LDA)
    xc_energy, xc_potential = compute_lda_xc(system.density, volume)
    system.potential.exchange .= xc_potential
    system.potential.correlation .= 0.0
    
    # For jellium model: add positive background potential to cancel Hartree
    # The positive background creates a potential V_b = -V_H
    # So total electrostatic potential = V_H + V_b = 0
    # This means the Kohn-Sham potential is just V_xc
    
    # Compute energy components
    system.energies.hartree = compute_hartree_energy(system.density.data, hartree_potential, volume)
    system.energies.exchange = xc_energy
    system.energies.correlation = 0.0
    
    # For uniform electron gas, add kinetic energy estimate
    # Kinetic energy from Thomas-Fermi: T = (3/10) * (3*pi^2)^(2/3) * integral(n^(5/3)) dr
    # For uniform density: T = (3/10) * (3*pi^2)^(2/3) * n^(2/3) * N
    uniform_density = n_electrons / volume
    if uniform_density > 0
        tf_kinetic = (3.0 / 10.0) * (3.0 * pi^2)^(2/3) * uniform_density^(2/3) * n_electrons
        system.energies.kinetic = tf_kinetic
    else
        system.energies.kinetic = 0.0
    end
    
    # Total energy
    system.energies.total = system.energies.kinetic + system.energies.hartree + system.energies.exchange
    
    # For uniform electron gas, the density remains uniform
    # This is the correct self-consistent solution for jellium
    new_density_data = copy(system.density.data)
    
    return system.energies.total, new_density_data
end

"""
    Run the self-consistent field calculation.
    
    Args:
    - system: DFTSystem
    - params: SCFParameters
    
    Returns:
    - system: Updated DFTSystem with converged density and energy
"""
function run_scf!(system::DFTSystem, params::SCFParameters)
    # Initialize density
    initialize_uniform_density(system)
    
    # Store old density for mixing
    old_density = copy(system.density.data)
    old_energy = 0.0
    
    # Print SCF header
    println("
***  SCF Iterations")
    println()
    println(" iSCF   Total Energy (Ha)    Energy Diff (Ha)    Density Diff")
    println("-"^78)
    
    # Print SCF header
    println("
***  SCF Iterations")
    println()
    println(" iSCF   Total Energy (Ha)    Energy Diff (Ha)    Density Diff")
    println("-"^78)
    
    # SCF loop
    for iteration in 1:params.max_iter
        # Store current density and energy
        old_density .= system.density.data
        old_energy = system.energies.total
        
        # Run SCF iteration
        new_energy, new_density = scf_iteration!(system, params)
        
        # Update density with mixing
        system.density.data .= apply_mixing(new_density, old_density, params)
        system.energies.total = new_energy
        
        # Calculate differences for output
        energy_diff = abs(system.energies.total - old_energy)
        density_diff = maximum(abs.(system.density.data - old_density))
        
        # Print iteration info
        iter_pad = lpad(string(iteration), 5)
        energy_pad = lpad(string(round(system.energies.total; digits=10)), 18)
        ediff_pad = lpad(string(round(energy_diff; digits=2)), 18)
        ddiff_pad = lpad(string(round(density_diff; digits=2)), 14)
        println("    ", iter_pad, "   ", energy_pad, "   ", ediff_pad, "   ", ddiff_pad)
", 
        
        # Check convergence
        if check_convergence(system, params, old_energy, old_density)
            println()
            println("SCF converged in ", iteration, " iterations")
            break
        end
    end
    
    return system
end

"""
    Self-consistent field calculation (non-mutating version).
    
    Args:
    - system: DFTSystem
    - params: SCFParameters
    
    Returns:
    - converged_system: New DFTSystem with converged results
"""
function self_consistent_field(system::DFTSystem, params::SCFParameters)
    # Create a copy of the system
    converged_system = deepcopy(system)
    
    # Run SCF
    run_scf!(converged_system, params)
    
    return converged_system
end

end # module SelfConsistent
