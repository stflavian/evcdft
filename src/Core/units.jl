"""
Unit conversions and utility functions for working with atomic units.
"""

module Units

using ..Constants

# Export commonly used conversions
export a0, E_h, angstrom_to_bohr, bohr_to_angstrom, hartree_to_ev, ev_to_hartree

"""
    Convert a length from angstroms to bohr.
"""
angstrom_to_bohr(x) = x * (1.0 / 0.5291772109038427)

"""
    Convert a length from bohr to angstroms.
"""
bohr_to_angstrom(x) = x * 0.5291772109038427

"""
    Convert energy from Hartree to eV.
"""
hartree_to_ev(x) = x * 27.211386245988

"""
    Convert energy from eV to Hartree.
"""
ev_to_hartree(x) = x * (1.0 / 27.211386245988)

"""
    Convert energy from Hartree to kcal/mol.
"""
hartree_to_kcal_mol(x) = x * 27.211386245988 * 23.060548

"""
    Convert energy from kcal/mol to Hartree.
"""
kcal_mol_to_hartree(x) = x / (27.211386245988 * 23.060548)

"""
    Convert a vector of lengths from angstroms to bohr.
"""
function angstrom_to_bohr(v::AbstractVector)
    return v .* (1.0 / 0.5291772109038427)
end

"""
    Convert a vector of lengths from bohr to angstroms.
"""
function bohr_to_angstrom(v::AbstractVector)
    return v .* 0.5291772109038427
end

"""
    Convert a matrix of lengths from angstroms to bohr.
"""
function angstrom_to_bohr(m::AbstractMatrix)
    return m .* (1.0 / 0.5291772109038427)
end

"""
    Convert a matrix of lengths from bohr to angstroms.
"""
function bohr_to_angstrom(m::AbstractMatrix)
    return m .* 0.5291772109038427
end

"""
    Convert a 3D array of lengths from angstroms to bohr.
"""
function angstrom_to_bohr(a::AbstractArray{T,3}) where T
    return a .* (1.0 / 0.5291772109038427)
end

"""
    Convert a 3D array of lengths from bohr to angstroms.
"""
function bohr_to_angstrom(a::AbstractArray{T,3}) where T
    return a .* 0.5291772109038427
end

# Type aliases for clarity
const Length = Float64  # In bohr
const Energy = Float64  # In Hartree
const Angle = Float64   # In radians

end # module Units
