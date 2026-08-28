"""
Command-line interface for EVC_DFT.

This module provides a simple CLI to run DFT calculations from input files.
Usage:
    julia --project cli.jl <input.hsd>
    julia --project -e 'using EVC_DFT; main()' <input.hsd>
"""

module CLI

using ..InputParser: run_from_input, parse_input_file
using ..EVC_DFT: DFTSystem

# Export the main function
export main

"""
    Main entry point for the CLI.
    
    Usage:
        julia --project cli.jl <input.hsd>
        julia --project -e 'using EVC_DFT; main()' <input.hsd>
        julia --project -e 'using EVC_DFT; main()' --help
"""
function main()
    args = Base.ARGV
    
    # Check for help flag
    if "--help" in args || "-h" in args
        print_usage()
        return
    end
    
    if "--version" in args || "-v" in args
        print_version()
        return
    end
    
    # Check if input file is provided
    if length(args) == 0
        println("Error: No input file provided")
        print_usage()
        return
    end
    
    # Get input file (first non-flag argument)
    input_file = args[1]
    
    # Check if file exists
    if !isfile(input_file)
        println("Error: Input file '$input_file' not found")
        print_usage()
        return
    end
    
    # Run the calculation
    try
        result = run_from_input(input_file)
        println("\nCalculation completed successfully!")
        return result
    catch e
        println("\nError during calculation:")
        showerror(stdout, e)
        return nothing
    end
end

"""
    Print usage information.
"""
function print_usage()
    println("""
EVC_DFT - Density Functional Theory Calculator

Usage:
    julia --project cli.jl <input.hsd>
    julia --project -e 'using EVC_DFT; main()' <input.hsd>

Options:
    --help, -h    Show this help message
    --version, -v Show version information

Input File:
    The input file should be in DFTB+ compatible HSD format.
    
    Example (jellium.hsd):
    ---------------------
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
    
    Hamiltonian = DFT {
        XCFunctional = "LDA-PZ"
    }
    
    Options = {
        EnergyCutoff [Ha] = 10.0
        FFTGrid = (16, 16, 16)
        SCF {
            MaxIter = 50
            EnergyTolerance = 1e-8
            DensityTolerance = 1e-8
            Mixing = "Linear"
            MixingParameter = 0.5
        }
    }
    
    Driver = {
        CalculationType = "Jellium"
    }

    Example (atomic system, Phase 2A):
    --------------------------------
    Geometry = {
        TypeNames = { "Ga" "As" }
        TypesAndCoordinates [Angstrom] = {
            1 0.0 0.0 0.0
            2 1.356773 1.356773 1.356773
        }
        Periodic = Yes
        LatticeVectors [Angstrom] = {
            2.713546 2.713546 0.0
            0.0 2.713546 2.713546
            2.713556 0.0 2.713546
        }
    }
    
    Options = {
        EnergyCutoff [Ha] = 10.0
        FFTGrid = (16, 16, 16)
    }

Supported Features:
    Phase 1: Uniform electron gas (jellium) calculations
    Phase 2A: Atomic systems with periodic boundary conditions (Gamma-point only)
    - LDA exchange-correlation (Perdew-Zunger)
    - Plane wave basis sets
    - Linear and Kerker density mixing
    - Self-consistent field calculations

For more information, see the EVC_DFT documentation.
""")
end

"""
    Print version information.
"""
function print_version()
    println("EVC_DFT v0.2.0")
    println("Phase 2A: Atomic Systems with Periodic Boundary Conditions")
    println("GitHub: https://github.com/stflavian/evcdft")
end

end # module CLI
