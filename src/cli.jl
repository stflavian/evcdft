"""
Command-line interface for EVC_DFT.

This module provides a simple CLI to run DFT calculations from input files.
Usage:
    julia --project -e 'using EVC_DFT; main()' <input.hsd>
"""

module CLI

using .InputParser: run_from_input, parse_input_file
using .EVC_DFT: DFTSystem

# Export the main function
export main

"""
    Main entry point for the CLI.
    
    Usage:
        julia --project -e 'using EVC_DFT; main()' <input.hsd>
        julia --project -e 'using EVC_DFT; main()' --help
"""
function main()
    args = Base.ARGV
    
    if length(args) == 0
        print_usage()
        return
    end
    
    # Check for help flag
    if "--help" in args || "-h" in args
        print_usage()
        return
    end
    
    if "--version" in args || "-v" in args
        print_version()
        return
    end
    
    # Get input file
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
    catch e
        println("\nError during calculation:")
        showerror(stdout, e)
        return
    end
end

"""
    Print usage information.
"""
function print_usage()
    println("""
EVC_DFT - Density Functional Theory Calculator

Usage:
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

Supported Features (Phase 1):
    - Uniform electron gas (jellium) calculations
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
    println("EVC_DFT v0.1.0")
    println("Phase 1: Uniform Electron Gas (Jellium)")
    println("GitHub: https://github.com/stflavian/evcdft")
end

"""
    Interactive mode for testing.
    
    This allows running calculations interactively without command-line arguments.
"""
function interactive_mode()
    println("EVC_DFT Interactive Mode")
    println("Enter 'help' for commands, 'quit' to exit")
    
    while true
        print("evcdft> ")
        input = readline()
        
        if input == "quit" || input == "exit" || input == "q"
            break
        elseif input == "help" || input == "h"
            println("Commands:")
            println("  run <file>    - Run calculation from input file")
            println("  test          - Run test suite")
            println("  version       - Show version")
            println("  quit          - Exit")
        elseif startswith(input, "run ")
            filename = input[5:end]
            if isfile(filename)
                result = run_from_input(filename)
            else
                println("Error: File '$filename' not found")
            end
        elseif input == "test"
            println("Running tests...")
            include("test/test_phase1.jl")
        elseif input == "version"
            print_version()
        else
            println("Unknown command. Type 'help' for available commands.")
        end
    end
end

end # module CLI
