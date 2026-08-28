"""
Input parser for DFTB+ compatible HSD files.

This module provides functions to parse DFTB+ style input files
and convert them into EVC_DFT data structures.
"""

module InputParser

using ..Types: Lattice, DFTSystem, SCFParameters, PlaneWaveBasis, AtomicSystem
using ..HSDParser: HSDNode, parse_hsd_file, parse_hsd_string, get_node, get_value
using ..EVC_DFT: create_uniform_electron_gas, jellium_total_energy
using ..Constants
using ..Units

# Export main functions
export parse_input_file, run_from_input, InputConfig
    println("Driver keys: ", keys(config.driver))
    println("Net charge: ", config.net_charge)
    println("="^60 * "\n")
end

"""
    Configuration structure extracted from input file.
    
    Fields:
    - geometry: Geometry configuration
    - hamiltonian: Hamiltonian configuration
    - options: Options configuration
    - driver: Driver configuration
    - net_charge: Net charge of the system (default: 0)
"""
mutable struct InputConfig
    geometry::Dict{String, Any}
    hamiltonian::Dict{String, Any}
    options::Dict{String, Any}
    driver::Dict{String, Any}
    net_charge::Int
    
    function InputConfig()
        new(Dict(), Dict(), Dict(), Dict(), 0)
    end
    println("Net charge: ", config.net_charge)
    println("="^60 * "\n")
end
end

"""
    Extract configuration from HSD root node.
    
    Args:
    - root: Root HSDNode from parser
    
    Returns:
    - InputConfig structure
"""
function extract_config(root::HSDNode)::InputConfig
    config = InputConfig()
    
    # Debug: print parsed root structure
    println("
DEBUG: Root children: ", keys(root.children))
    
    # Extract Geometry block
    if haskey(root.children, "Geometry")
        config.geometry = extract_geometry(root.children["Geometry"])
    end
    
    # Extract Hamiltonian block
    if haskey(root.children, "Hamiltonian")
        config.hamiltonian = extract_hamiltonian(root.children["Hamiltonian"])
    end
    
    # Extract Options block
    if haskey(root.children, "Options")
        config.options = extract_options(root.children["Options"])
    end
    
    # Extract Driver block
    if haskey(root.children, "Driver")
        config.driver = extract_driver(root.children["Driver"])
    end
    
    # Extract net charge from Driver or Geometry
    if haskey(config.driver, "NetCharge")
        config.net_charge = Int(config.driver["NetCharge"])
    elseif haskey(config.geometry, "NetCharge")
        config.net_charge = Int(config.geometry["NetCharge"])
    end
    
    return config
end

"""
    Extract geometry configuration from Geometry block.
"""
function extract_geometry(geometry_node::HSDNode)::Dict{String, Any}
    config = Dict{String, Any}()
    
    # Check for LatticeVectors
    if haskey(geometry_node.children, "LatticeVectors")
        lattice_node = geometry_node.children["LatticeVectors"]
        lattice_config = extract_lattice_vectors(lattice_node)
        config["LatticeVectors"] = lattice_config
    end
    
    # Check for ElectronGas
    if haskey(geometry_node.children, "ElectronGas")
        eg_node = geometry_node.children["ElectronGas"]
        config["ElectronGas"] = extract_electron_gas(eg_node)
    end
    
    # Check for TypeNames (DFTB+ format for atomic systems)
    if haskey(geometry_node.children, "TypeNames")
        type_names_node = geometry_node.children["TypeNames"]
        config["TypeNames"] = extract_type_names(type_names_node)
    end
    
    # Check for TypesAndCoordinates (DFTB+ format for atomic systems)
    if haskey(geometry_node.children, "TypesAndCoordinates")
        tac_node = geometry_node.children["TypesAndCoordinates"]
        units = tac_node.units
        config["TypesAndCoordinates"] = extract_types_and_coordinates(tac_node, units)
    end
    
    # Check for Periodic flag
    if haskey(geometry_node.children, "Periodic")
        periodic_node = geometry_node.children["Periodic"]
        config["Periodic"] = get_value(periodic_node)
    else
        config["Periodic"] = true  # Default to periodic
    end
    
    # Check for NetCharge in Geometry block
    if haskey(geometry_node.children, "NetCharge")
        config["NetCharge"] = get_value(geometry_node.children["NetCharge"])
    end
    
    return config
end

"""
    Extract lattice vectors from LatticeVectors block.
"""
function extract_lattice_vectors(node::HSDNode)::Dict{String, Any}
    config = Dict{String, Any}()
    
    # Get units if specified
    units = node.units
    
    # Lattice vectors can be specified as:
    # 1. Three separate vectors
    # 2. A single matrix
    
    if haskey(node.children, "1") && haskey(node.children, "2") && haskey(node.children, "3")
        # Three separate vectors
        vec1 = get_value(node.children["1"])
        vec2 = get_value(node.children["2"])
        vec3 = get_value(node.children["3"])
        
        # Handle unit conversion if specified on the parent
        if units == "Angstrom" || units == "angstrom"
            vec1 = vec1 .* angstrom_to_bohr
            vec2 = vec2 .* angstrom_to_bohr
            vec3 = vec3 .* angstrom_to_bohr
        end
        
        config["a1"] = vec1
        config["a2"] = vec2
        config["a3"] = vec3
        config["units"] = units
    else
        # Try to get as a matrix or array
        if node.value isa Vector
            # Flatten the vector and reshape into 3x3
            flat = node.value
            if length(flat) == 9
                a1 = flat[1:3]
                a2 = flat[4:6]
                a3 = flat[7:9]
                config["a1"] = a1
                config["a2"] = a2
                config["a3"] = a3
            elseif length(flat) == 3
                # Single vector (cubic lattice)
                a = flat[1]
                config["a1"] = [a, 0.0, 0.0]
                config["a2"] = [0.0, a, 0.0]
                config["a3"] = [0.0, 0.0, a]
            end
        end
        config["units"] = units
    end
    
    return config
end

"""
    Extract electron gas configuration.
"""
function extract_electron_gas(node::HSDNode)::Dict{String, Any}
    config = Dict{String, Any}()
    
    for (key, child) in node.children
        config[key] = get_value(child)
    end
    
    return config
end

"""
    Extract type names from TypeNames block (DFTB+ format).
    
    Args:
    - node: HSDNode for TypeNames block
    
    Returns:
    - Vector of species names (Strings)
"""
function extract_type_names(node::HSDNode)::Vector{String}
    if node.value isa Vector
        return String.(node.value)
    elseif node.value isa String
        # Could be space-separated or comma-separated
        return split(node.value)
    else
        # Try to extract from children
        names = String[]
        for (key, child) in node.children
            val = get_value(child)
            if val isa String
                push!(names, val)
            end
        end
        return names
    end
end

"""
    Extract types and coordinates from TypesAndCoordinates block (DFTB+ format).
    
    Args:
    - node: HSDNode for TypesAndCoordinates block
    - units: Units specification (e.g., "Angstrom")
    
    Returns:
    - Matrix of type indices and coordinates (N_atoms x 4)
"""
function extract_types_and_coordinates(node::HSDNode, units::String)::Matrix{Float64}
    n_atoms = length(node.children)
    result = zeros(Float64, n_atoms, 4)
    
    for (idx, (key, child)) in enumerate(node.children)
        # Parse the value: could be "1 0.0 0.0 0.0" or [1, 0.0, 0.0, 0.0]
        val = get_value(child)
        if val isa Vector
            result[idx, :] = Float64.(val)
        elseif val isa String
            parts = parse.(Float64, split(val))
            result[idx, :] = parts
        end
    end
    
    # Convert coordinates from units to Bohr if needed
    if units == "Angstrom" || units == "angstrom"
        result[:, 2:4] .*= angstrom_to_bohr
    end
    
    return result
end

"""
    Extract atomic positions (for compatibility, deprecated in favor of TypesAndCoordinates).
"""
function extract_positions(node::HSDNode)::Vector{Dict{String, Any}}
    positions = Vector{Dict{String, Any}}()
    
    for (key, child) in node.children
        pos = Dict{String, Any}()
        pos["index"] = parse(Int, key)
        pos["coords"] = get_value(child)
        push!(positions, pos)
    end
    
    return positions
end

"""
    Extract Hamiltonian configuration.
"""
function extract_hamiltonian(hamiltonian_node::HSDNode)::Dict{String, Any}
    config = Dict{String, Any}()
    
    for (key, child) in hamiltonian_node.children
        config[key] = get_value(child)
    end
    
    return config
end

"""
    Extract options configuration.
"""
function extract_options(options_node::HSDNode)::Dict{String, Any}
    config = Dict{String, Any}()
    
    for (key, child) in options_node.children
        if key == "SCF"
            config[key] = extract_scf_options(child)
        elseif key == "FFTGrid"
            # FFTGrid can be a tuple or array
            value = get_value(child)
            if value isa Vector
                config[key] = Tuple(Int.(value))
            else
                config[key] = value
            end
        elseif key == "EnergyCutoff"
            config[key] = get_value(child)
        else
            config[key] = get_value(child)
        end
    end
    
    return config
end

"""
    Extract SCF options.
"""
function extract_scf_options(scf_node::HSDNode)::Dict{String, Any}
    config = Dict{String, Any}()
    
    for (key, child) in scf_node.children
        config[key] = get_value(child)
    end
    
    return config
end

"""
    Extract driver configuration.
"""
function extract_driver(driver_node::HSDNode)::Dict{String, Any}
    config = Dict{String, Any}()
    
    for (key, child) in driver_node.children
        config[key] = get_value(child)
    end
    
    return config
end

"""
    Validate input configuration.
    
    Args:
    - config: InputConfig to validate
    
    Throws:
    - Error if configuration is invalid
"""
function validate_input(config::InputConfig)
    # Check required blocks
    if isempty(config.geometry)
        error("Geometry block is required")
    end
    
    # For atomic systems, TypeNames + TypesAndCoordinates are valid
    # For jellium, ElectronGas is valid
    # For both, LatticeVectors can be present
    has_atomic = haskey(config.geometry, "TypeNames") && haskey(config.geometry, "TypesAndCoordinates")
    has_jellium = haskey(config.geometry, "ElectronGas")
    has_lattice = haskey(config.geometry, "LatticeVectors")
    
    if !has_atomic && !has_jellium
        error("Geometry block must contain ElectronGas or (TypeNames and TypesAndCoordinates)")
    end
    
    # For atomic systems
    if haskey(config.geometry, "TypeNames") && haskey(config.geometry, "TypesAndCoordinates")
        if !haskey(config.geometry, "LatticeVectors")
            error("LatticeVectors required for atomic systems")
        end
        if length(config.geometry["TypeNames"]) == 0
            error("TypeNames must contain at least one species")
        end
        if size(config.geometry["TypesAndCoordinates"], 1) == 0
            error("TypesAndCoordinates must contain at least one atom")
        end
    end
    
    # For jellium
    if haskey(config.geometry, "ElectronGas")
        if !haskey(config.geometry["ElectronGas"], "NumElectrons")
            error("ElectronGas block must contain NumElectrons")
        end
    end
    
    if isempty(config.options)
        error("Options block is required")
    end
    
    if !haskey(config.options, "EnergyCutoff")
        config.options["EnergyCutoff"] = 10.0  # Default
    end
    
    if !haskey(config.options, "FFTGrid")
        config.options["FFTGrid"] = (16, 16, 16)  # Default
    end
    
    if !haskey(config.options, "SCF")
        config.options["SCF"] = Dict()
    end
end

"""
    Build a Lattice from geometry configuration.
    
    Args:
    - config: InputConfig
    
    Returns:
    - Lattice object
"""
function build_lattice(config::InputConfig)::Lattice
    lattice_config = config.geometry["LatticeVectors"]
    
    a1 = lattice_config["a1"]
    a2 = lattice_config["a2"]
    a3 = lattice_config["a3"]
    
    return Lattice(a1, a2, a3)
end

"""
    Build SCFParameters from options configuration.
    
    Args:
    - config: InputConfig
    
    Returns:
    - SCFParameters object
"""
function build_scf_parameters(config::InputConfig)::SCFParameters
    scf_config = get(config.options, "SCF", Dict{String, Any}())
    
    return SCFParameters(
        max_iter = get(scf_config, "MaxIter", 100),
        energy_tolerance = get(scf_config, "EnergyTolerance", 1e-6),
        density_tolerance = get(scf_config, "DensityTolerance", 1e-6),
        mixing_parameter = get(scf_config, "MixingParameter", 0.5),
        mixing_type = get(scf_config, "Mixing", "Linear")
    )
end

"""
    Build an AtomicSystem from input configuration.
    
    Args:
    - config: InputConfig
    
    Returns:
    - AtomicSystem object
"""
function build_atomic_system(config::InputConfig)::AtomicSystem
    lattice = build_lattice(config)
    species_names = config.geometry["TypeNames"]
    types_and_coords = config.geometry["TypesAndCoordinates"]
    periodic = get(config.geometry, "Periodic", true)
    net_charge = config.net_charge
    
    return AtomicSystem(lattice, species_names, types_and_coords,
                       net_charge=net_charge, periodic=periodic)
end

"""
    Build a DFTSystem from input configuration.
    
    Args:
    - config: InputConfig
    
    Returns:
    - DFTSystem object
"""
function build_system(config::InputConfig)::DFTSystem
    # For jellium (uniform electron gas)
    if haskey(config.geometry, "ElectronGas")
        lattice = build_lattice(config)
        n_electrons = config.geometry["ElectronGas"]["NumElectrons"]
        cutoff = config.options["EnergyCutoff"]
        fft_size = config.options["FFTGrid"]
        
        return create_uniform_electron_gas(lattice.a1[1], n_electrons, cutoff, fft_size)
    # For atomic systems
    elseif haskey(config.geometry, "TypeNames") && haskey(config.geometry, "TypesAndCoordinates")
        atomic_system = build_atomic_system(config)
        lattice = atomic_system.lattice
        cutoff = config.options["EnergyCutoff"]
        fft_size = config.options["FFTGrid"]
        
        # Create plane wave basis
        basis = PlaneWaveBasis(lattice, cutoff, fft_size)
        
        # Create DFT system with atomic system
        return DFTSystem(lattice, basis, atomic_system)
    else
        error("Geometry must contain ElectronGas or (TypeNames and TypesAndCoordinates)")
    end
end

"""
    Parse an input file and return a DFTSystem.
    
    Args:
    - filename: Path to the input file
    
    Returns:
    - DFTSystem ready for calculation
"""
function parse_input_file(filename::String)::DFTSystem
    # Parse the HSD file
    root = parse_hsd_file(filename)
    
    # Extract configuration
    config = extract_config(root)
    
    # Validate input
    validate_input(config)
    
    # Build and return the system
    return build_system(config)
end

"""
    Parse an input string and return a DFTSystem.
    
    Args:
    - content: Input content as a string
    
    Returns:
    - DFTSystem ready for calculation
"""
function parse_input_string(content::String)::DFTSystem
    root = parse_hsd_string(content)
    config = extract_config(root)
    validate_input(config)
    return build_system(config)
end

"""
    Run a calculation from an input file.
    
    Args:
    - filename: Path to the input file
    
    Returns:
    - Converged DFTSystem with results
"""
function run_from_input(filename::String)::DFTSystem
    # Parse input
    system = parse_input_file(filename)
    
    # Get SCF parameters
    config = extract_config(parse_hsd_file(filename))
    params = build_scf_parameters(config)
    
    # Run SCF
    println("Running SCF calculation from input file: $filename")
    result = self_consistent_field(system, params)
    
    # Print summary
    print_summary(result, config)
    
    return result
end

"""
    Print a summary of the calculation results.
    
    Args:
    - result: Converged DFTSystem
    - config: InputConfig used for the calculation
"""
function print_summary(result::DFTSystem, config::InputConfig)
    println("\n" * "="^60)
    println("EVC_DFT - Calculation Summary")
    println("="^60)
    
    # Geometry
    println("\nGeometry:")
    println("  Lattice vectors:")
    println("    a1 = $(result.lattice.a1)")
    println("    a2 = $(result.lattice.a2)")
    println("    a3 = $(result.lattice.a3)")
    println("  Volume = $(result.lattice.volume) Bohr^3")
    
    # Atomic system info
    if result.atomic_system !== nothing
        println("\nAtomic System:")
        println("  Number of atoms = $(length(result.atomic_system.atoms))")
        println("  Species = $(result.atomic_system.species_names)")
        println("  Total nuclear charge = $(result.atomic_system.total_nuclear_charge)")
        println("  Net charge = $(result.atomic_system.net_charge)")
        println("  Periodic = $(result.atomic_system.periodic)")
    end
    
    # Electrons
    println("\nElectrons:")
    println("  Number = $(result.electrons)")
    
    # Basis
    println("\nBasis:")
    println("  Type = PlaneWave")
    println("  Energy cutoff = $(result.basis.cutoff) Hartree")
    println("  FFT grid = $(result.basis.fft_size)")
    println("  Number of G vectors = $(result.basis.n_g)")
    
    # SCF Parameters
    println("\nSCF Parameters:")
    scf_config = get(config.options, "SCF", Dict())
    println("  Max iterations = $(get(scf_config, "MaxIter", 100))")
    println("  Energy tolerance = $(get(scf_config, "EnergyTolerance", 1e-6)) Hartree")
    println("  Density tolerance = $(get(scf_config, "DensityTolerance", 1e-6)) Bohr^-3")
    println("  Mixing = $(get(scf_config, "Mixing", "Linear"))")
    println("  Mixing parameter = $(get(scf_config, "MixingParameter", 0.5))")
    
    # Results
    println("\nResults:")
    println("  Total energy = $(result.energies.total) Hartree")
    println("  Hartree energy = $(result.energies.hartree) Hartree")
    println("  Exchange energy = $(result.energies.exchange) Hartree")
    println("  Correlation energy = $(result.energies.correlation) Hartree")
    
    # Hamiltonian
    if haskey(config.hamiltonian, "XCFunctional")
        println("\nHamiltonian:")
        println("  XC Functional = $(config.hamiltonian["XCFunctional"])")
    end
    
    println("\n" * "="^60)
end

"""
    Run a calculation from an input string.
    
    Args:
    - content: Input content as a string
    
    Returns:
    - Converged DFTSystem with results
"""
function run_from_input_string(content::String)::DFTSystem
    root = parse_hsd_string(content)
    config = extract_config(root)
    validate_input(config)
    
    system = build_system(config)
    params = build_scf_parameters(config)
    
    result = self_consistent_field(system, params)
    print_summary(result, config)
    
    return result
end

end # module InputParser
