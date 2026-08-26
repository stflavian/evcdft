"""
Physical constants for DFT calculations in atomic units (Hartree).

Atomic units:
- Length: Bohr (a0 = 4πϵ0 ħ² / m_e e²)
- Energy: Hartree (E_h = m_e e⁴ / (4πϵ0)² ħ²)
- Mass: Electron mass (m_e)
- Charge: Elementary charge (e)
"""

module Constants

# Fundamental constants in SI units
const ħ_SI = 1.0545718176461565e-34      # Reduced Planck constant (J·s)
const m_e_SI = 9.109383701528254e-31     # Electron mass (kg)
const e_SI = 1.602176634e-19             # Elementary charge (C)
const ε0_SI = 8.854187812890987e-12     # Vacuum permittivity (F/m)
const c_SI = 299792458.0                 # Speed of light (m/s)

# Derived atomic units
const a0 = 4π * ε0_SI * ħ_SI^2 / (m_e_SI * e_SI^2)  # Bohr radius (m)
const E_h = m_e_SI * e_SI^4 / (4π * ε0_SI)^2 / ħ_SI^2  # Hartree energy (J)

# Atomic unit values (all constants = 1 in atomic units)
const ħ = 1.0      # Reduced Planck constant (atomic units)
const m_e = 1.0    # Electron mass (atomic units)
const e = 1.0      # Elementary charge (atomic units)
const ε0 = 1.0/(4π) # Vacuum permittivity (atomic units)

# Conversion factors
const angstrom_to_bohr = 1.0 / 0.5291772109038427   # 1 Å = 1.8897259886 a0
const bohr_to_angstrom = 0.5291772109038427        # 1 a0 = 0.529177 Å

const hartree_to_ev = 27.211386245988     # 1 Ha = 27.2114 eV
const ev_to_hartree = 1.0 / 27.211386245988

const hartree_to_joule = E_h
const joule_to_hartree = 1.0 / E_h

const hartree_to_kcal_mol = hartree_to_ev * 23.060548   # 1 eV = 23.0605 kcal/mol
const kcal_mol_to_hartree = 1.0 / (hartree_to_ev * 23.060548)

# Mathematical constants
const π = Base.MathConstants.pi
const twopi = 2 * π
const sqrtpi = sqrt(π)
const sqrt2 = sqrt(2.0)

# Useful combinations
const fourpi = 4 * π
const twopi_sqrt = sqrt(twopi)

# Electron volt in atomic units
const ryberg = 0.5  # 1 Ry = 0.5 Ha

# Bohr magneton in atomic units
const μ_B = 0.5  # μ_B = e ħ / (2 m_e) in atomic units

# Fine structure constant
const α = e_SI^2 / (4π * ε0_SI * ħ_SI * c_SI)  # ≈ 1/137

end # module Constants
