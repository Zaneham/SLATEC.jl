# SLATEC.jl

A Julia port of SLATEC — the Sandia, Los Alamos, Air Force Weapons Laboratory Technical Exchange Committee library.

## What This Is

SLATEC is a collection of ~1,400 mathematical routines written in Fortran between 1982 and 1993 by scientists at US national laboratories. It contains some of the most battle-tested numerical code ever written.

This package ports those routines to Julia, one function at a time, with obsessive attention to correctness.

The code was written by people who couldn't afford to be wrong. When your calculations determine whether a missile hits its target or a reactor stays stable, you check your work. Then you check it again. Then you have someone else check it.

Modern software moves fast and breaks things. SLATEC moved carefully and broke nothing.

## Why Bother?

Julia already has numerical libraries. Many of them are quite good. But most modern implementations optimise for speed first and correctness second. SLATEC did it the other way around.

Modern numerical libraries are written by very clever people with PhDs and tenure. They optimise for speed, which is sensible. They assume IEEE 754 floating-point, which is reasonable. They trust the compiler, which is bold.

SLATEC was written by people who knew what hexadecimal floating-point did to your mantissa at 3am. They handled edge cases that modern developers don't know exist. They wrote code that would give the same answer on a CDC 6600 as on a Cray-1, which is a bit like writing code that runs identically on a calculator and a spacecraft.

The original SLATEC code has been validated against:
- IBM System/360 mainframes (hexadecimal floating-point)
- Abramowitz & Stegun reference tables
- NIST Digital Library of Mathematical Functions
- Decades of production use in weapons simulation, aerospace, and scientific computing

## On Testing

The original SLATEC authors validated their code against Abramowitz & Stegun and decades of production use. That was enough for the 1980s.

It's not enough now.

The four-level test framework in [slatec-modern](https://github.com/Zaneham/slatec) asks whether the code works, whether it matches mathematical truth, whether it matches IBM 360 historical output, and whether it survives hostile compiler flags. Most routines pass. Some don't. J₀(500) returns the wrong sign, which is the sort of thing you'd rather know before your spacecraft misses Mars.

If a routine doesn't pass all four levels, it doesn't get ported. There's no shame in admitting failure. The shame would be in shipping it anyway.

| Level | Question |
|-------|----------|
| **L1** | Does the code work? |
| **L2** | Does the output match mathematical truth? |
| **L3** | Does the output match IBM 360 historical values? |
| **L4** | Does the output survive hostile compiler flags? |

### How Julia Validation Works

The Fortran code in slatec-modern has been validated against IBM 360 hardware, mathematical references, and hostile compiler flags. Running Julia on a 1966 mainframe would be impressive but impractical.

Instead, the Julia tests verify that output matches the validated Fortran. If Fortran passed L1-L4, and Julia matches Fortran, then Julia inherits the validation. The hard work was done once; the Julia port piggybacks on it.

## Current Status

| Module | Routines | Status |
|--------|----------|--------|
| BLAS | DAXPY, DSCAL, DCOPY, DSWAP, DROT, DROTG | Pass |
| MINPACK | DENORM, ENORM | Pass |
| LINPACK | DGEFA, DGESL, DPOFA, DPOSL | Planned |
| Special Functions | Bessel, Gamma, Airy, Elliptic | Planned |

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Zaneham/SLATEC.jl")
```

## Usage

```julia
using SLATEC

# Euclidean norm with overflow/underflow protection
x = [1e-300, 1e-300, 1e-300]
denorm(x)  # Won't flush to zero like naive implementations

# BLAS Level 1
y = [1.0, 2.0, 3.0]
daxpy!(2.0, x, y)  # y = 2*x + y
```

## Philosophy

Every routine in this package exists because:

1. The original SLATEC implementation is mathematically sound
2. The algorithm has been validated against authoritative references
3. The code produces identical output to IBM 360/370 mainframes from the 1980s
4. The implementation survives adversarial compiler optimisations

Routines with known accuracy failures are documented but not ported. See [slatec-modern](https://github.com/Zaneham/slatec) for the four-level test framework and deviation analysis.

## Contact

- Issues and contributions: [GitHub Issues](https://github.com/Zaneham/SLATEC.jl/issues)
- Author: [Zane Hambly](https://github.com/Zaneham)

## License

Apache License 2.0. See [LICENSE](LICENSE) for details.

The original SLATEC library is public domain.
