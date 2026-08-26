"""
Test runner for EVC_DFT package.

This script runs all tests for the EVC_DFT package.
Usage: julia runtests.jl
"""

# Include the main module
include("../src/evcdft.jl")

# Use the module
using .EVC_DFT

# Import all submodules for testing
using .EVC_DFT.Constants
using .EVC_DFT.Types
using .EVC_DFT.PlaneWave
using .EVC_DFT.XCFunctionals
using .EVC_DFT.SelfConsistent

using Test

println("Running EVC_DFT test suite...")
println("="^60)

# Include all test files
include("unit_tests.jl")
include("integration_tests.jl")
include("test_phase1.jl")

println("="^60)
println("\nAll tests completed successfully!")
