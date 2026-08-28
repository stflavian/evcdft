"""
Unit tests for Phase 2A: Data Structures for Atomic Systems.

This file contains tests for:
- Atomic number lookup (Task 1.1)
- New types: AtomicSpecies, Atom, AtomicSystem, Wavefunction, KohnShamStates (Tasks 1.2-1.6)
- DFTSystem updates (Task 1.7)
- Utility functions (Tasks 2.1-2.2)
- Input parser updates (Tasks 4.1-4.6)

Run with: julia --project test/test_phase2a.jl
"""

using Test
using LinearAlgebra

# Include the main module
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "EVC_DFT.jl"))

# Use the module
using .EVC_DFT
using EVC_DFT.Constants
using EVC_DFT.Types

# =============================================================================
# TASK 1.1: Atomic Number Lookup Tests
# =============================================================================

@testset "Atomic Number Lookup - Basic Elements" begin
    # Test first few elements
    @test get_atomic_number("H") == 1
    @test get_atomic_number("He") == 2
    @test get_atomic_number("Li") == 3
    @test get_atomic_number("Be") == 4
    @test get_atomic_number("B") == 5
    @test get_atomic_number("C") == 6
    @test get_atomic_number("N") == 7
    @test get_atomic_number("O") == 8
    @test get_atomic_number("F") == 9
    @test get_atomic_number("Ne") == 10
end

@testset "Atomic Number Lookup - Second Period" begin
    @test get_atomic_number("Na") == 11
    @test get_atomic_number("Mg") == 12
    @test get_atomic_number("Al") == 13
    @test get_atomic_number("Si") == 14
    @test get_atomic_number("P") == 15
    @test get_atomic_number("S") == 16
    @test get_atomic_number("Cl") == 17
    @test get_atomic_number("Ar") == 18
end

@testset "Atomic Number Lookup - Third Period and Beyond" begin
    @test get_atomic_number("K") == 19
    @test get_atomic_number("Ca") == 20
    @test get_atomic_number("Sc") == 21
    @test get_atomic_number("Ti") == 22
    @test get_atomic_number("V") == 23
    @test get_atomic_number("Cr") == 24
    @test get_atomic_number("Mn") == 25
    @test get_atomic_number("Fe") == 26
    @test get_atomic_number("Co") == 27
    @test get_atomic_number("Ni") == 28
    @test get_atomic_number("Cu") == 29
    @test get_atomic_number("Zn") == 30
end

@testset "Atomic Number Lookup - GaAs Relevant Elements" begin
    # Elements needed for GaAs example
    @test get_atomic_number("Ga") == 31
    @test get_atomic_number("Ge") == 32
    @test get_atomic_number("As") == 33
    @test get_atomic_number("Se") == 34
    @test get_atomic_number("Br") == 35
    @test get_atomic_number("Kr") == 36
end

@testset "Atomic Number Lookup - Case Insensitivity" begin
    # Test case insensitivity
    @test get_atomic_number("h") == 1
    @test get_atomic_number("H") == 1
    @test get_atomic_number("He") == 2
    @test get_atomic_number("he") == 2
    @test get_atomic_number("HE") == 2
    @test get_atomic_number("ga") == 31
    @test get_atomic_number("Ga") == 31
    @test get_atomic_number("GA") == 31
    @test get_atomic_number("as") == 33
    @test get_atomic_number("As") == 33
    @test get_atomic_number("AS") == 33
end

@testset "Atomic Number Lookup - Dictionary Access" begin
    # Test direct dictionary access
    @test Constants.ATOMIC_NUMBERS["H"] == 1
    @test Constants.ATOMIC_NUMBERS["He"] == 2
    @test Constants.ATOMIC_NUMBERS["Ga"] == 31
    @test Constants.ATOMIC_NUMBERS["As"] == 33
end

@testset "Atomic Number Lookup - Error Handling" begin
    # Test that unknown elements throw an error
    @test_throws ErrorException get_atomic_number("Xx")
    @test_throws ErrorException get_atomic_number("Unknown")
    @test_throws ErrorException get_atomic_number("")
end

# =============================================================================
# TASK 1.2-1.7: Type Definitions (Placeholder tests - will be implemented later)
# =============================================================================

@testset "AtomicSpecies Type" begin
    # Test construction
    ga = AtomicSpecies("Ga", 31, mass=69.723)
    @test ga.name == "Ga"
    @test ga.atomic_number == 31
    @test ga.mass ≈ 69.723
    
    # Test default mass
    as = AtomicSpecies("As", 33)
    @test as.name == "As"
    @test as.atomic_number == 33
    @test as.mass == 0.0
end

@testset "Atom Type" begin
    lattice = SimpleCubic(5.0)
    
    # Test construction with origin
    atom1 = Atom(1, [0.0, 0.0, 0.0], lattice)
    @test atom1.species_index == 1
    @test atom1.position == [0.0, 0.0, 0.0]
    @test atom1.fractional_position == [0.0, 0.0, 0.0]
    
    # Test construction with non-origin position
    atom2 = Atom(1, [1.0, 2.0, 3.0], lattice)
    @test atom2.species_index == 1
    @test atom2.position == [1.0, 2.0, 3.0]
    @test atom2.fractional_position ≈ [0.2, 0.4, 0.6]
    
    # Test construction with position at lattice boundary
    atom3 = Atom(1, [5.0, 0.0, 0.0], lattice)
    @test atom3.position == [5.0, 0.0, 0.0]
    @test atom3.fractional_position == [1.0, 0.0, 0.0]
end

@testset "AtomicSystem Type" begin
    lattice = SimpleCubic(5.0)
    species_names = ["Ga", "As"]
    types_and_coords = [1 0.0 0.0 0.0; 2 1.0 1.0 1.0]  # Bohr coordinates
    
    # Test neutral system
    atomic_system = AtomicSystem(lattice, species_names, types_and_coords)
    @test length(atomic_system.atoms) == 2
    @test atomic_system.species_names == ["Ga", "As"]
    @test atomic_system.atomic_numbers == [31, 33]
    @test atomic_system.total_nuclear_charge == 64
    @test atomic_system.n_valence_electrons == 64
    @test atomic_system.net_charge == 0
    @test atomic_system.periodic == true
    
    # Test with net charge
    atomic_system_charged = AtomicSystem(lattice, species_names, types_and_coords, net_charge=2)
    @test atomic_system_charged.net_charge == 2
    @test atomic_system_charged.n_valence_electrons == 62
    
    # Test with non-periodic
    atomic_system_molecule = AtomicSystem(lattice, species_names, types_and_coords, periodic=false)
    @test atomic_system_molecule.periodic == false
end

@testset "Wavefunction Type" begin
    data = rand(ComplexF64, 8, 8, 8)
    wf = Wavefunction(data, -0.5, 1.0, 1)
    
    @test wf.energy == -0.5
    @test wf.occupation == 1.0
    @test wf.band_index == 1
    @test size(wf.data) == (8, 8, 8)
    @test wf.data === data  # Same object reference
end

@testset "KohnShamStates Type" begin
    states = KohnShamStates(10)
    
    @test states.n_bands == 10
    @test length(states.wavefunctions) == 10
    @test length(states.eigenvalues) == 10
    @test states.fermi_energy == 0.0
    
    # Test that wavefunctions are uninitialized
    @test states.wavefunctions[1] === nothing || typeof(states.wavefunctions[1]) <: Wavefunction
end

@testset "DFTSystem with AtomicSystem" begin
    lattice = SimpleCubic(5.0)
    species_names = ["H"]
    types_and_coords = [1 0.0 0.0 0.0]  # Bohr
    atomic_system = AtomicSystem(lattice, species_names, types_and_coords)
    basis = PlaneWaveBasis(lattice, 10.0, (8, 8, 8))
    system = DFTSystem(lattice, basis, atomic_system)
    
    @test system.atomic_system === atomic_system
    @test system.electrons == 1
    @test system.states !== nothing
    @test system.states.n_bands > 0
    @test system.temperature == 0.0
    @test system.smearing_type == "none"
end

# =============================================================================
# TASK 2.1-2.2: Utility Functions (Placeholder tests)
# =============================================================================

@testset "Coordinate Conversion" begin
    lattice = SimpleCubic(5.0)
    
    # Test origin
    @test cartesian_to_fractional([0.0, 0.0, 0.0], lattice) == [0.0, 0.0, 0.0]
    @test fractional_to_cartesian([0.0, 0.0, 0.0], lattice) == [0.0, 0.0, 0.0]
    
    # Test lattice vector
    @test cartesian_to_fractional([5.0, 0.0, 0.0], lattice) == [1.0, 0.0, 0.0]
    @test cartesian_to_fractional([0.0, 5.0, 0.0], lattice) == [0.0, 1.0, 0.0]
    @test cartesian_to_fractional([0.0, 0.0, 5.0], lattice) == [0.0, 0.0, 1.0]
    
    # Test fractional point
    @test cartesian_to_fractional([2.5, 2.5, 2.5], lattice) == [0.5, 0.5, 0.5]
    @test fractional_to_cartesian([0.5, 0.5, 0.5], lattice) == [2.5, 2.5, 2.5]
    
    # Test round-trip conversion
    for x in 0.0:0.1:5.0, y in 0.0:0.1:5.0, z in 0.0:0.1:5.0
        frac = cartesian_to_fractional([x, y, z], lattice)
        cart = fractional_to_cartesian(frac, lattice)
        @test cart ≈ [x, y, z] atol=1e-10
    end
end

@testset "Atomic System Utilities" begin
    lattice = SimpleCubic(5.0)
    species_names = ["H", "He"]
    types_and_coords = [1 0.0 0.0 0.0; 2 1.0 1.0 1.0]
    atomic_system = AtomicSystem(lattice, species_names, types_and_coords)
    
    @test total_nuclear_charge(atomic_system) == 3
    @test total_electrons(atomic_system) == 3
    
    positions = get_atom_positions(atomic_system)
    @test length(positions) == 2
    @test positions[1] == [0.0, 0.0, 0.0]
    @test positions[2] == [1.0, 1.0, 1.0]
    
    frac_positions = get_fractional_positions(atomic_system)
    @test length(frac_positions) == 2
    @test frac_positions[1] == [0.0, 0.0, 0.0]
    @test frac_positions[2] ≈ [0.2, 0.2, 0.2]
end

# =============================================================================
# TASK 4.1-4.6: Input Parser Tests (Placeholder tests)
# =============================================================================

@testset "Input Parser - GaAs Example" begin
    input_content = """
    Geometry = {
      TypeNames = { "Ga" "As" }
      TypesAndCoordinates [Angstrom] = {
        1 0.000000 0.000000 0.000000
        2 1.356773 1.356773 1.356773
      }
      Periodic = Yes
      LatticeVectors [Angstrom] = {
        2.713546 2.713546 0.0
        0.0 2.713546 2.713546
        2.713546 0.0 2.713546
      }
    }
    
    Options = {
      EnergyCutoff [Ha] = 10.0
      FFTGrid = (16, 16, 16)
    }
    
    Driver = {
      CalculationType = "Atomic"
    }
    """
    
    # Write to temp file
    temp_file = tempname() * ".hsd"
    write(temp_file, input_content)
    
    try
        # Parse the input
        system = parse_input_file(temp_file)
        
        # Check atomic system
        @test system.atomic_system !== nothing
        @test length(system.atomic_system.atoms) == 2
        @test system.atomic_system.species_names == ["Ga", "As"]
        @test system.atomic_system.atomic_numbers == [31, 33]
        @test system.atomic_system.total_nuclear_charge == 64
        @test system.electrons == 64  # Neutral system
        @test system.atomic_system.net_charge == 0
        @test system.atomic_system.periodic == true
        
        # Check lattice
        @test system.lattice.volume > 0
        
        # Check basis
        @test system.basis.cutoff ≈ 10.0
        @test system.basis.fft_size == (16, 16, 16)
        
    finally
        rm(temp_file, force=true)
    end
end

@testset "Input Parser - With Net Charge" begin
    input_content = """
    Geometry = {
      TypeNames = { "H" }
      TypesAndCoordinates [Angstrom] = {
        1 0.0 0.0 0.0
      }
      Periodic = No
    }
    
    Options = {
      EnergyCutoff [Ha] = 10.0
      FFTGrid = (8, 8, 8)
    }
    
    Driver = {
      CalculationType = "Atomic"
      NetCharge = 1
    }
    """
    
    temp_file = tempname() * ".hsd"
    write(temp_file, input_content)
    
    try
        system = parse_input_file(temp_file)
        
        @test system.atomic_system !== nothing
        @test system.atomic_system.net_charge == 1
        @test system.electrons == 0  # H+ ion has no electrons
        @test system.atomic_system.periodic == false
        
    finally
        rm(temp_file, force=true)
    end
end

@testset "Input Parser - Jellium Still Works" begin
    input_content = """
    Geometry = {
      LatticeVectors [Angstrom] {
        5.0  0.0  0.0
        0.0  5.0  0.0
        0.0  0.0  5.0
      }
      ElectronGas {
        NumElectrons = 2
      }
    }
    
    Options = {
      EnergyCutoff [Ha] = 10.0
      FFTGrid = (16, 16, 16)
    }
    
    Driver = {
      CalculationType = "Jellium"
    }
    """
    
    temp_file = tempname() * ".hsd"
    write(temp_file, input_content)
    
    try
        system = parse_input_file(temp_file)
        
        @test system.atomic_system === nothing
        @test system.electrons == 2
        @test system.lattice.volume > 0
        
    finally
        rm(temp_file, force=true)
    end
end

@testset "create_atomic_system Helper" begin
    # Test with simple cubic lattice
    lattice_constant = 5.0  # Bohr
    species_names = ["H", "He"]
    types_and_coords = [1 0.0 0.0 0.0; 2 1.0 1.0 1.0]  # Angstrom
    
    system = create_atomic_system(lattice_constant, species_names, types_and_coords)
    
    @test system.atomic_system !== nothing
    @test length(system.atomic_system.atoms) == 2
    @test system.atomic_system.species_names == ["H", "He"]
    @test system.electrons == 3
    @test system.basis.fft_size == (16, 16, 16)
    
    # Test with explicit lattice vectors
    lattice_vectors = [5.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 5.0]
    system2 = create_atomic_system(lattice_vectors, species_names, types_and_coords,
                                    cutoff=15.0, fft_size=(32, 32, 32))
    
    @test system2.basis.cutoff ≈ 15.0
    @test system2.basis.fft_size == (32, 32, 32)
end

# =============================================================================
# Run all tests and print summary
# =============================================================================

println("\n" * "="^60)
println("Phase 2A Tests Starting")
println("="^60)
