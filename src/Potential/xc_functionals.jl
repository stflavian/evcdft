"""
Exchange-correlation functionals for DFT calculations.

This module implements the Local Density Approximation (LDA) for exchange and correlation.
"""

module XCFunctionals

using ..Constants: π, twopi, sqrtpi, fourpi
using ..Types: ElectronDensity

# Export the main functions
export lda_exchange_energy, lda_correlation_energy, lda_xc_energy,
       lda_exchange_potential, lda_correlation_potential, lda_xc_potential,
       compute_lda_energy, compute_lda_potential

"""
    Exchange energy per electron in LDA (Dirac exchange).
    
    For a uniform electron gas:
    ε_x = - (3/4) * (3/π)^(1/3) * n^(1/3)
    
    Args:
    - density: Electron density (in Bohr⁻³)
    
    Returns:
    - exchange_energy_density: Exchange energy per electron (Hartree)
"""
function lda_exchange_energy(density::Float64)
    if density <= 0.0
        return 0.0
    end
    # Dirac exchange: ε_x = - (3/4π) * (3π²)^(1/3) * n^(1/3)
    # More standard form: ε_x = - (3/4) * (3/π)^(1/3) * n^(1/3)
    return - (3.0 / (4.0 * π)) * (3.0 * π^2)^(1/3) * density^(1/3)
end

"""
    Correlation energy per electron in LDA (Perdew-Zunger parameterization).
    
    For a uniform electron gas, the correlation energy is parameterized as:
    ε_c = -C * (1 + δ₁√r_s + δ₂r_s + δ₃r_s^(3/2) + ...) / (1 + β₁√r_s + β₂r_s)
    where r_s = (3/(4πn))^(1/3) is the Wigner-Seitz radius.
    
    Args:
    - density: Electron density (in Bohr⁻³)
    
    Returns:
    - correlation_energy_density: Correlation energy per electron (Hartree)
"""
function lda_correlation_energy(density::Float64)
    if density <= 0.0
        return 0.0
    end
    
    # Wigner-Seitz radius
    rs = (3.0 / (4π * density))^(1/3)
    
    # Perdew-Zunger parameterization (from Ceperley-Alder quantum Monte Carlo)
    # For rs >= 1
    if rs >= 1.0
        sqrt_rs = sqrt(rs)
        num = -0.1423 * (1.0 + 1.0529 * sqrt_rs + 0.3334 * rs)
        den = 1.0 + 1.0529 * sqrt_rs + 0.3334 * rs + 0.00264 * rs * sqrt_rs
        return num / den
    else
        # For rs < 1 (high density)
        log_rs = log(rs)
        num = -0.0480 * (1.0 + 1.1310 * sqrt(rs) + 0.5282 * rs + 0.0960 * rs * sqrt(rs))
        den = 1.0 + 1.6344 * sqrt(rs) + 0.4929 * rs + 0.1164 * rs * sqrt(rs)
        return num / den + (0.0311 * log_rs - 0.0470)
    end
end

"""
    Total LDA exchange-correlation energy per electron.
    
    Args:
    - density: Electron density (in Bohr⁻³)
    
    Returns:
    - xc_energy_density: XC energy per electron (Hartree)
"""
function lda_xc_energy(density::Float64)
    return lda_exchange_energy(density) + lda_correlation_energy(density)
end

"""
    Exchange potential in LDA.
    
    V_x = d(ε_x n) / dn = ε_x + n dε_x/dn
    For Dirac exchange: V_x = - (3/π)^(1/3) * n^(1/3)
    
    Args:
    - density: Electron density (in Bohr⁻³)
    
    Returns:
    - exchange_potential: Exchange potential (Hartree)
"""
function lda_exchange_potential(density::Float64)
    if density <= 0.0
        return 0.0
    end
    # V_x = - (3/π)^(1/3) * n^(1/3) * (4/3)
    # Actually: V_x = d(n ε_x)/dn = ε_x + n * dε_x/dn
    # ε_x = -C * n^(1/3), so dε_x/dn = -C/3 * n^(-2/3)
    # V_x = -C * n^(1/3) + n * (-C/3 * n^(-2/3)) = -C * n^(1/3) - C/3 * n^(1/3) = -4C/3 * n^(1/3)
    # where C = (3/4π) * (3π²)^(1/3)
    C = (3.0 / (4π)) * (3π^2)^(1/3)
    return - (4.0 / 3.0) * C * density^(1/3)
end

"""
    Correlation potential in LDA (Perdew-Zunger).
    
    V_c = d(ε_c n) / dn = ε_c + n dε_c/dn
    
    Args:
    - density: Electron density (in Bohr⁻³)
    
    Returns:
    - correlation_potential: Correlation potential (Hartree)
"""
function lda_correlation_potential(density::Float64)
    if density <= 0.0
        return 0.0
    end
    
    rs = (3.0 / (4π * density))^(1/3)
    
    if rs >= 1.0
        sqrt_rs = sqrt(rs)
        # Parameters for Perdew-Zunger
        A = -0.1423
        α1 = 1.0529
        α2 = 0.3334
        β1 = 1.0529
        β2 = 0.3334
        β3 = 0.00264
        
        denom = 1.0 + β1 * sqrt_rs + β2 * rs + β3 * rs * sqrt_rs
        ddenom_drs = 0.5 * (β1 / sqrt_rs + 2β2 + 1.5 * β3 * sqrt_rs)
        
        num = A * (1.0 + α1 * sqrt_rs + α2 * rs)
        dnum_drs = A * (0.5 * α1 / sqrt_rs + α2)
        
        ε_c = num / denom
        dεc_drs = (dnum_drs * denom - num * ddenom_drs) / (denom^2)
        
        # V_c = ε_c + n * dε_c/dn = ε_c - (rs/3) * dε_c/drs
        V_c = ε_c - (rs / 3.0) * dεc_drs
        return V_c
    else
        # For rs < 1
        sqrt_rs = sqrt(rs)
        A = -0.0480
        α1 = 1.1310
        α2 = 0.5282
        α3 = 0.0960
        β1 = 1.6344
        β2 = 0.4929
        β3 = 0.1164
        
        denom = 1.0 + β1 * sqrt_rs + β2 * rs + β3 * rs * sqrt_rs
        ddenom_drs = 0.5 * (β1 / sqrt_rs + 2β2 + 1.5 * β3 * sqrt_rs)
        
        num = A * (1.0 + α1 * sqrt_rs + α2 * rs + α3 * rs * sqrt_rs)
        dnum_drs = A * (0.5 * α1 / sqrt_rs + α2 + 1.5 * α3 * sqrt_rs)
        
        ε_c = num / denom + (0.0311 * log(rs) - 0.0470)
        dεc_drs = (dnum_drs * denom - num * ddenom_drs) / (denom^2) + 0.0311 / rs
        
        V_c = ε_c - (rs / 3.0) * dεc_drs
        return V_c
    end
end

"""
    Total LDA exchange-correlation potential.
    
    Args:
    - density: Electron density (in Bohr⁻³)
    
    Returns:
    - xc_potential: XC potential (Hartree)
"""
function lda_xc_potential(density::Float64)
    return lda_exchange_potential(density) + lda_correlation_potential(density)
end

"""
    Compute the total LDA exchange-correlation energy for a given density.
    
    Args:
    - density: ElectronDensity structure
    
    Returns:
    - xc_energy: Total XC energy (Hartree)
"""
function compute_lda_energy(density::ElectronDensity)
    data = density.data
    volume = prod(density.grid_size)  # Approximate volume (for uniform grid)
    
    xc_energy = 0.0
    for i in eachindex(data)
        xc_energy += lda_xc_energy(data[i])
    end
    
    # Multiply by volume element (assuming uniform grid)
    # For a proper calculation, we need the actual volume
    # This is a simplified version
    return xc_energy * (volume / length(data))
end

"""
    Compute the LDA exchange-correlation potential for a given density.
    
    Args:
    - density: ElectronDensity structure
    
    Returns:
    - xc_potential: Array of XC potential values (Hartree)
"""
function compute_lda_potential(density::ElectronDensity)
    data = density.data
    potential = similar(data)
    
    for i in eachindex(data)
        potential[i] = lda_xc_potential(data[i])
    end
    
    return potential
end

"""
    Compute XC energy and potential together (more efficient).
    
    Args:
    - density: Electron density in real space
    
    Returns:
    - xc_energy: Total XC energy (Hartree)
    - xc_potential: XC potential array (Hartree)
"""
function compute_lda_xc(density::ElectronDensity, volume::Float64)
    data = density.data
    nx, ny, nz = density.grid_size
    
    xc_energy = 0.0
    xc_potential = similar(data)
    
    dV = volume / (nx * ny * nz)  # Volume element
    
    for i in eachindex(data)
        dens = data[i]
        xc_potential[i] = lda_xc_potential(dens)
        xc_energy += lda_xc_energy(dens) * dens
    end
    
    return xc_energy * dV, xc_potential
end

end # module XCFunctionals
