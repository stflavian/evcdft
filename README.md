# EVC_DFT - Density Functional Theory Calculator in Julia

A modular, extensible DFT calculator implemented in Julia. Currently supports **uniform electron gas (jellium) calculations** with plane wave basis sets and LDA exchange-correlation functionals.

## Features

### Phase 1 (Completed) ✅
- **Uniform Electron Gas (Jellium)**: Full support for jellium model calculations
- **Plane Wave Basis**: FFT-based operations with energy cutoff
- **LDA Functionals**: Perdew-Zunger exchange-correlation
- **SCF Loop**: Self-consistent field with density mixing
- **DFTB+ Input Format**: Compatible with DFTB+ HSD file format

### Phase 2 (Planned)
- Pseudopotential support (UPF format)
- Atomic systems with periodic boundary conditions
- Monkhorst-Pack k-point sampling
- Local and non-local potential handling

### Phase 3 (Planned)
- GGA functionals (PBE, BLYP)
- Spin-polarized calculations
- Gaussian basis sets for molecules
- Forces and stress tensor

### Phase 4 (Planned)
- Hybrid functionals (PBE0, HSE)
- Meta-GGA functionals
- GPU acceleration (CUDA.jl)
- MPI parallelization

## Installation

```bash
# Clone the repository
git clone https://github.com/stflavian/evcdft.git
cd evcdft

# Instantiate the Julia environment
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Usage

### Command-Line Interface

```bash
# Run a calculation from an input file
julia --project -e 'using EVC_DFT; run_from_input("input.hsd")'

# Or use the main CLI
julia --project -e 'using EVC_DFT; main()' input.hsd

# Show help
julia --project -e 'using EVC_DFT; main()' --help
```

### Julia REPL

```julia
using EVC_DFT

# Create a uniform electron gas system
system = create_uniform_electron_gas(5.0, 2, 10.0, (16, 16, 16))

# Set up SCF parameters
params = SCFParameters(max_iter=50, energy_tolerance=1e-8)

# Run self-consistent calculation
result = self_consistent_field(system, params)

# Print results
println("Total energy: ", result.energies.total, " Hartree")
```

### Using Input Files

Create a file `jellium.hsd`:

```hsd
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
```

Then run:
```bash
julia --project -e 'using EVC_DFT; run_from_input("jellium.hsd")'
```

## Input File Format

EVC_DFT uses a **subset of the DFTB+ HSD format**. See `examples/` for working input files.

### Supported Blocks

| Block | Description | Required |
|-------|-------------|----------|
| `Geometry` | Lattice and electron specification | ✅ Yes |
| `Geometry.LatticeVectors` | Lattice vectors (3x3 matrix) | ✅ Yes |
| `Geometry.ElectronGas` | Electron count for jellium | ✅ Yes |
| `Hamiltonian` | Exchange-correlation functional | ✅ Yes |
| `Options` | Calculation parameters | ✅ Yes |
| `Options.EnergyCutoff` | Plane wave cutoff (Hartree) | ❌ No (default: 10.0) |
| `Options.FFTGrid` | FFT grid size (tuple) | ❌ No (default: (16,16,16)) |
| `Options.SCF` | SCF convergence parameters | ❌ No |
| `Driver` | Calculation type | ❌ No |

### Supported Units

| Unit | Description | Conversion |
|------|-------------|------------|
| `[Angstrom]` | Length in Angstrom | → Bohr |
| `[Bohr]` | Length in Bohr | No conversion |
| `[Ha]` or `[Hartree]` | Energy in Hartree | No conversion |
| `[eV]` | Energy in eV | → Hartree |
| `[nm]` | Length in nanometers | → Bohr |

### Example Input Files

- `examples/jellium.hsd` - Full jellium example with all options
- `examples/minimal.hsd` - Minimal input with defaults
- `examples/jellium_angstrom.hsd` - Example with Angstrom units

## Project Structure

```
evc_dft/
├── Project.toml              # Julia project file
├── README.md                 # This file
├── src/
│   ├── Core/
│   │   ├── constants.jl      # Physical constants
│   │   ├── units.jl          # Unit conversions
│   │   └── types.jl          # Data structures
│   ├── Basis/
│   │   └── plane_wave.jl     # Plane wave basis
│   ├── Potential/
│   │   └── xc_functionals.jl # XC functionals
│   ├── SCF/
│   │   └── self_consistent.jl # SCF loop
│   ├── IO/
│   │   ├── hsd_parser.jl     # HSD format parser
│   │   └── input_parser.jl   # Input file parser
│   ├── cli.jl                # Command-line interface
│   └── evcdft.jl             # Main module
├── test/
│   ├── test_phase1.jl        # Phase 1 tests
│   └── test_input.jl         # Input parser tests
└── examples/
    ├── jellium.hsd            # Jellium example
    ├── minimal.hsd            # Minimal example
    └── jellium_angstrom.hsd   # Angstrom units example
```

## API Documentation

### Core Types

```julia
# Lattice representation
struct Lattice
    a1::Vector{Float64}
    a2::Vector{Float64}
    a3::Vector{Float64}
    volume::Float64
    b1::Vector{Float64}  # Reciprocal lattice
    b2::Vector{Float64}
    b3::Vector{Float64}
end

# DFT system
struct DFTSystem
    lattice::Lattice
    basis::PlaneWaveBasis
    electrons::Int
    density::ElectronDensity
    potential::KohnShamPotential
    energies::EnergyComponents
end

# SCF parameters
struct SCFParameters
    max_iter::Int
    energy_tolerance::Float64
    density_tolerance::Float64
    mixing_parameter::Float64
    mixing_type::String
end
```

### Main Functions

```julia
# Create a uniform electron gas system
create_uniform_electron_gas(a::Float64, n_electrons::Int, cutoff::Float64, fft_size::Tuple)

# Run self-consistent field calculation
self_consistent_field(system::DFTSystem, params::SCFParameters)

# Parse input file
parse_input_file(filename::String)
run_from_input(filename::String)

# Parse input string
parse_input_string(content::String)
run_from_input_string(content::String)
```

### Exchange-Correlation Functionals

```julia
# LDA exchange energy per electron (Dirac)
lda_exchange_energy(density::Float64)

# LDA correlation energy per electron (Perdew-Zunger)
lda_correlation_energy(density::Float64)

# LDA exchange-correlation potential
lda_xc_potential(density::Float64)
```

## Testing

Run the test suite:

```bash
julia --project test/test_phase1.jl
julia --project test/test_input.jl
```

Or run all tests:
```bash
julia --project -e 'include("test/test_phase1.jl"); include("test/test_input.jl")'
```

## Performance Considerations

- **FFTW**: Pre-planning FFTs for better performance
- **Memory**: Pre-allocate arrays to reduce GC overhead
- **Parallelism**: Threads.jl support for multi-threaded FFTs
- **GPU**: CUDA.jl support planned for Phase 4

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/name`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin feature/name`)
5. Open a Pull Request

## License

MIT License - see LICENSE file for details.

## Acknowledgments

- Inspired by Quantum ESPRESSO, VASP, and DFTB+
- Built with Julia for performance and productivity
- Uses FFTW.jl for efficient FFT operations

## Contact

For questions or issues, please open a GitHub issue.
