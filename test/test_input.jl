"""
Tests for the HSD parser and input file handling.

This file contains tests for:
- HSD parser (lexer and parser)
- Input configuration extraction
- Unit conversion
- Input validation
"""

using Test
using EVC_DFT
using EVC_DFT.HSDParser
using EVC_DFT.InputParser
using EVC_DFT.Constants
using EVC_DFT.Units

@testset "HSD Lexer" begin
    # Test simple tokenization
    input = "Key = Value"
    tokens = HSDParser.tokenize(input)
    @test length(tokens) == 3  # Key, =, Value, EOF
    @test tokens[1].type == HSDParser.IDENTIFIER
    @test tokens[1].value == "Key"
    @test tokens[2].type == HSDParser.EQUALS
    @test tokens[3].type == HSDParser.IDENTIFIER
    @test tokens[3].value == "Value"
    
    # Test number tokenization
    input = "Number = 3.14159"
    tokens = HSDParser.tokenize(input)
    @test tokens[1].type == HSDParser.IDENTIFIER
    @test tokens[3].type == HSDParser.NUMBER
    @test tokens[3].value == "3.14159"
    
    # Test negative numbers
    input = "Negative = -5.0"
    tokens = HSDParser.tokenize(input)
    @test tokens[3].type == HSDParser.NUMBER
    @test tokens[3].value == "-5.0"
    
    # Test scientific notation
    input = "Scientific = 1.0e-10"
    tokens = HSDParser.tokenize(input)
    @test tokens[3].type == HSDParser.NUMBER
    @test tokens[3].value == "1.0e-10"
    
    # Test strings
    input = "String = \"hello\""
    tokens = HSDParser.tokenize(input)
    @test tokens[3].type == HSDParser.STRING
    @test tokens[3].value == "hello"
    
    # Test blocks
    input = "Block = { }"
    tokens = HSDParser.tokenize(input)
    @test tokens[3].type == HSDParser.BEGIN_BLOCK
    @test tokens[5].type == HSDParser.END_BLOCK
    
    # Test arrays
    input = "Array = (1, 2, 3)"
    tokens = HSDParser.tokenize(input)
    @test tokens[3].type == HSDParser.LPAREN
    @test tokens[4].type == HSDParser.NUMBER
    @test tokens[6].type == HSDParser.COMMA
    @test tokens[8].type == HSDParser.RPAREN
    
    # Test comments
    input = "# This is a comment"
    tokens = HSDParser.tokenize(input)
    @test length(tokens) == 1  # Only EOF
    
    input = "// This is also a comment"
    tokens = HSDParser.tokenize(input)
    @test length(tokens) == 1
    
    # Test units
    input = "Length [Angstrom] = 5.0"
    tokens = HSDParser.tokenize(input)
    @test tokens[3].type == HSDParser.LBRACKET
    @test tokens[4].type == HSDParser.IDENTIFIER
    @test tokens[4].value == "Angstrom"
    @test tokens[6].type == HSDParser.NUMBER
end

@testset "HSD Parser" begin
    # Test simple key-value
    input = "Key = Value"
    root = parse_hsd_string(input)
    @test haskey(root.children, "Key")
    @test root.children["Key"].value == "Value"
    
    # Test number values
    input = "Number = 42"
    root = parse_hsd_string(input)
    @test root.children["Number"].value == 42
    
    input = "Float = 3.14"
    root = parse_hsd_string(input)
    @test root.children["Float"].value ≈ 3.14
    
    # Test boolean values
    input = "Flag = Yes"
    root = parse_hsd_string(input)
    @test root.children["Flag"].value == true
    
    input = "Flag = No"
    root = parse_hsd_string(input)
    @test root.children["Flag"].value == false
    
    # Test blocks
    input = """
    Parent = {
        Child = Value
    }
    """
    root = parse_hsd_string(input)
    @test haskey(root.children, "Parent")
    @test haskey(root.children["Parent"].children, "Child")
    @test root.children["Parent"].children["Child"].value == "Value"
    
    # Test nested blocks
    input = """
    Level1 = {
        Level2 = {
            Level3 = Value
        }
    }
    """
    root = parse_hsd_string(input)
    @test root.children["Level1"].children["Level2"].children["Level3"].value == "Value"
    
    # Test arrays
    input = "Array = (1, 2, 3)"
    root = parse_hsd_string(input)
    @test root.children["Array"].value == [1, 2, 3]
    
    input = "Array2 = (1.0, 2.0, 3.0)"
    root = parse_hsd_string(input)
    @test root.children["Array2"].value == [1.0, 2.0, 3.0]
    
    # Test units
    input = "Length [Angstrom] = 5.0"
    root = parse_hsd_string(input)
    @test root.children["Length"].units == "Angstrom"
    @test root.children["Length"].value == 5.0
end

@testset "Unit Conversion" begin
    # Test Angstrom to Bohr conversion
    input = "Length [Angstrom] = 1.0"
    root = parse_hsd_string(input)
    node = root.children["Length"]
    value = get_value(node)
    @test value ≈ angstrom_to_bohr
    
    # Test eV to Hartree conversion
    input = "Energy [eV] = 1.0"
    root = parse_hsd_string(input)
    node = root.children["Energy"]
    value = get_value(node)
    @test value ≈ ev_to_hartree
    
    # Test array conversion
    input = "Vectors [Angstrom] = (1.0, 2.0, 3.0)"
    root = parse_hsd_string(input)
    node = root.children["Vectors"]
    value = get_value(node)
    expected = [1.0, 2.0, 3.0] .* angstrom_to_bohr
    @test value ≈ expected
end

@testset "Input Configuration Extraction" begin
    # Test parsing minimal input
    input = """
    Geometry = {
        LatticeVectors {
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
        EnergyCutoff = 10.0
        FFTGrid = (16, 16, 16)
        SCF {
            MaxIter = 50
        }
    }
    """
    
    root = parse_hsd_string(input)
    config = extract_config(root)
    
    # Check geometry
    @test haskey(config.geometry, "LatticeVectors")
    @test haskey(config.geometry, "ElectronGas")
    @test config.geometry["ElectronGas"]["NumElectrons"] == 2
    
    # Check hamiltonian
    @test config.hamiltonian["XCFunctional"] == "LDA-PZ"
    
    # Check options
    @test config.options["EnergyCutoff"] == 10.0
    @test config.options["FFTGrid"] == (16, 16, 16)
    @test config.options["SCF"]["MaxIter"] == 50
    
    # Test with units
    input = """
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
    }
    """
    
    root = parse_hsd_string(input)
    config = extract_config(root)
    lattice_config = config.geometry["LatticeVectors"]
    
    # Check that units were converted
    @test lattice_config["a1"] ≈ [5.0 * angstrom_to_bohr, 0.0, 0.0]
end

@testset "Input Validation" begin
    # Test missing Geometry block
    input = """
    Hamiltonian = DFT {
        XCFunctional = "LDA-PZ"
    }
    """
    
    root = parse_hsd_string(input)
    config = extract_config(root)
    @test_throws ErrorException("Geometry block is required") validate_input(config)
    
    # Test missing ElectronGas or LatticeVectors
    input = """
    Geometry = {
        SomeOtherBlock = {}
    }
    """
    
    root = parse_hsd_string(input)
    config = extract_config(root)
    @test_throws ErrorException("Geometry block must contain LatticeVectors or ElectronGas") validate_input(config)
    
    # Test missing NumElectrons
    input = """
    Geometry = {
        ElectronGas {
            SomeOtherKey = 5
        }
    }
    """
    
    root = parse_hsd_string(input)
    config = extract_config(root)
    @test_throws ErrorException("ElectronGas block must contain NumElectrons") validate_input(config)
    
    # Test valid input with defaults
    input = """
    Geometry = {
        LatticeVectors {
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
    """
    
    root = parse_hsd_string(input)
    config = extract_config(root)
    validate_input(config)
    
    # Check that defaults were added
    @test config.options["EnergyCutoff"] == 10.0
    @test config.options["FFTGrid"] == (16, 16, 16)
end

@testset "System Building" begin
    # Test building a system from input
    input = """
    Geometry = {
        LatticeVectors {
            5.0  0.0  0.0
            0.0  5.0  0.0
            0.0  0.0  5.0
        }
        ElectronGas {
            NumElectrons = 2
        }
    }
    Options = {
        EnergyCutoff = 10.0
        FFTGrid = (8, 8, 8)
    }
    """
    
    system = parse_input_string(input)
    
    @test system.lattice.volume ≈ 125.0
    @test system.electrons == 2
    @test system.basis.cutoff ≈ 10.0
    @test system.basis.fft_size == (8, 8, 8)
    
    # Test with Angstrom units
    input = """
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
        EnergyCutoff = 10.0
        FFTGrid = (8, 8, 8)
    }
    """
    
    system = parse_input_string(input)
    
    # 5 Angstrom = 5 * 1.8897259886 Bohr ≈ 9.4486 Bohr
    expected_volume = (5.0 * angstrom_to_bohr)^3
    @test system.lattice.volume ≈ expected_volume
end

@testset "SCF Parameters Building" begin
    input = """
    Options = {
        SCF {
            MaxIter = 100
            EnergyTolerance = 1e-8
            DensityTolerance = 1e-8
            Mixing = "Kerker"
            MixingParameter = 0.7
        }
    }
    """
    
    root = parse_hsd_string(input)
    config = extract_config(root)
    params = build_scf_parameters(config)
    
    @test params.max_iter == 100
    @test params.energy_tolerance ≈ 1e-8
    @test params.density_tolerance ≈ 1e-8
    @test params.mixing_type == "Kerker"
    @test params.mixing_parameter == 0.7
    
    # Test defaults
    input = """
    Options = {
        SCF {
            MaxIter = 50
        }
    }
    """
    
    root = parse_hsd_string(input)
    config = extract_config(root)
    params = build_scf_parameters(config)
    
    @test params.max_iter == 50
    @test params.energy_tolerance ≈ 1e-6  # Default
    @test params.density_tolerance ≈ 1e-6  # Default
    @test params.mixing_parameter == 0.5  # Default
end

@testset "File Parsing" begin
    # Test parsing from actual files
    # Note: These tests will only work if the files exist
    
    # Test jellium.hsd
    if isfile("examples/jellium.hsd")
        system = parse_input_file("examples/jellium.hsd")
        @test system.electrons == 2
        @test system.basis.cutoff ≈ 10.0
    end
    
    # Test minimal.hsd
    if isfile("examples/minimal.hsd")
        system = parse_input_file("examples/minimal.hsd")
        @test system.electrons == 2
    end
end

@testset "Run from Input" begin
    # Test running a calculation from input string
    input = """
    Geometry = {
        LatticeVectors {
            5.0  0.0  0.0
            0.0  5.0  0.0
            0.0  0.0  5.0
        }
        ElectronGas {
            NumElectrons = 2
        }
    }
    Options = {
        EnergyCutoff = 10.0
        FFTGrid = (8, 8, 8)
        SCF {
            MaxIter = 5
            EnergyTolerance = 1e-4
            DensityTolerance = 1e-4
        }
    }
    """
    
    # This should run without errors
    result = run_from_input_string(input)
    
    @test result.energies.total < 0.0  # Should have negative energy
    @test result.energies.hartree < 0.0
end

# Print test summary
println("\nInput Parser Tests Complete!")
println("All HSD parsing and input handling functionality is working.")
