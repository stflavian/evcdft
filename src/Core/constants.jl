"""
Physical constants for DFT calculations in atomic units (Hartree).

Atomic units:
- Length: Bohr (a0 = 4*pi*epsilon0*hbar^2 / m_e e^2)
- Energy: Hartree (E_h = m_e e^4 / (4*pi*epsilon0)^2 hbar^2)
- Mass: Electron mass (m_e)
- Charge: Elementary charge (e)
"""

module Constants

# Fundamental constants in SI units
const hbar_SI = 1.0545718176461565e-34      # Reduced Planck constant (J\u00b7s)
const m_e_SI = 9.109383701528254e-31     # Electron mass (kg)
const e_SI = 1.602176634e-19             # Elementary charge (C)
const epsilon0_SI = 8.854187812890987e-12     # Vacuum permittivity (F/m)
const c_SI = 299792458.0                 # Speed of light (m/s)

# Derived atomic units
const a0 = 4 * Base.MathConstants.pi * epsilon0_SI * hbar_SI^2 / (m_e_SI * e_SI^2)  # Bohr radius (m)
const E_h = m_e_SI * e_SI^4 / (4 * Base.MathConstants.pi * epsilon0_SI)^2 / hbar_SI^2  # Hartree energy (J)

# Atomic unit values (all constants = 1 in atomic units)
const hbar = 1.0      # Reduced Planck constant (atomic units)
const m_e = 1.0    # Electron mass (atomic units)
const e = 1.0      # Elementary charge (atomic units)
const epsilon0 = 1.0/(4 * Base.MathConstants.pi) # Vacuum permittivity (atomic units)

# Conversion factors
const angstrom_to_bohr = 1.0 / 0.5291772109038427   # 1 \u00c5 = 1.8897259886 a0
const bohr_to_angstrom = 0.5291772109038427        # 1 a0 = 0.529177 \u00c5

const hartree_to_ev = 27.211386245988     # 1 Ha = 27.2114 eV
const ev_to_hartree = 1.0 / 27.211386245988

const hartree_to_joule = E_h
const joule_to_hartree = 1.0 / E_h

const hartree_to_kcal_mol = hartree_to_ev * 23.060548   # 1 eV = 23.0605 kcal/mol
const kcal_mol_to_hartree = 1.0 / (hartree_to_ev * 23.060548)

# Mathematical constants
const pi = Base.MathConstants.pi
const twopi = 2 * pi
const sqrtpi = sqrt(pi)
const sqrt2 = sqrt(2.0)

# Useful combinations
const fourpi = 4 * pi
const twopi_sqrt = sqrt(twopi)

# Electron volt in atomic units
const ryberg = 0.5  # 1 Ry = 0.5 Ha

# Bohr magneton in atomic units
const mu_B = 0.5  # mu_B = e hbar / (2 m_e) in atomic units

# Fine structure constant
const alpha = e_SI^2 / (4 * pi * epsilon0_SI * hbar_SI * c_SI)  # \u2248 1/137

# Atomic numbers for common elements (hardcoded for convenience)
const ATOMIC_NUMBERS = Dict{String, Int}(
    # Period 1
    "H" => 1, "He" => 2,
    # Period 2
    "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7, "O" => 8, "F" => 9, "Ne" => 10,
    # Period 3
    "Na" => 11, "Mg" => 12, "Al" => 13, "Si" => 14, "P" => 15, "S" => 16, "Cl" => 17, "Ar" => 18,
    # Period 4
    "K" => 19, "Ca" => 20, "Sc" => 21, "Ti" => 22, "V" => 23, "Cr" => 24, "Mn" => 25,
    "Fe" => 26, "Co" => 27, "Ni" => 28, "Cu" => 29, "Zn" => 30,
    # Period 5
    "Ga" => 31, "Ge" => 32, "As" => 33, "Se" => 34, "Br" => 35, "Kr" => 36,
    "Rb" => 37, "Sr" => 38, "Y" => 39, "Zr" => 40, "Nb" => 41, "Mo" => 42,
    "Tc" => 43, "Ru" => 44, "Rh" => 45, "Pd" => 46, "Ag" => 47, "Cd" => 48,
    "In" => 49, "Sn" => 50, "Sb" => 51, "Te" => 52, "I" => 53, "Xe" => 54,
    # Period 6
    "Cs" => 55, "Ba" => 56, "La" => 57, "Ce" => 58, "Pr" => 59, "Nd" => 60,
    "Pm" => 61, "Sm" => 62, "Eu" => 63, "Gd" => 64, "Tb" => 65, "Dy" => 66,
    "Ho" => 67, "Er" => 68, "Tm" => 69, "Yb" => 70, "Lu" => 71, "Hf" => 72,
    "Ta" => 73, "W" => 74, "Re" => 75, "Os" => 76, "Ir" => 77, "Pt" => 78,
    "Au" => 79, "Hg" => 80, "Tl" => 81, "Pb" => 82, "Bi" => 83, "Po" => 84,
    "At" => 85, "Rn" => 86
)

"""
    Get the atomic number for a given element symbol.

    Args:
    - symbol: Element symbol (e.g., "H", "He", "Ga", "As")

    Returns:
    - Atomic number (Int)

    Throws:
    - ErrorException if element symbol is not found
"""
function get_atomic_number(symbol::String)::Int
    # Case-insensitive lookup: normalize to title case
    sym_normalized = uppercasefirst(lowercase(symbol))
    if haskey(ATOMIC_NUMBERS, sym_normalized)
        return ATOMIC_NUMBERS[sym_normalized]
    else
        error("Atomic number for element '" * symbol * "' not found in ATOMIC_NUMBERS dictionary")
    end
end

end # module Constants
