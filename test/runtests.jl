"""
Test runner for EVC_DFT package.

This script runs all tests for the EVC_DFT package.
Usage: julia runtests.jl
"""

# Add the parent directory to LOAD_PATH for proper module loading
push!(LOAD_PATH, joinpath(@__DIR__, ".."))

# Include the main module
include(joinpath(@__DIR__, "..", "src", "EVC_DFT.jl"))

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
include(joinpath(@__DIR__, "unit_tests.jl"))
include(joinpath(@__DIR__, "integration_tests.jl"))
include(joinpath(@__DIR__, "test_phase1.jl"))

println("="^60)
println("\nAll tests completed successfully!")
