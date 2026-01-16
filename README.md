# SLATEC.jl

A Julia port of SLATEC — the Sandia, Los Alamos, Air Force Weapons Laboratory Technical Exchange Committee library.

## What This Is

SLATEC is a collection of ~1,400 mathematical routines written in Fortran between 1982 and 1993 by scientists at US national laboratories. It contains some of the most battle-tested numerical code ever written.

This package ports those routines to Julia, one function at a time, with obsessive attention to correctness.

## Why Bother?

Julia already has numerical libraries. Many of them are quite good. But most modern implementations optimise for speed first and correctness second. SLATEC did it the other way around.

The original SLATEC code has been validated against:
- IBM System/360 mainframes (hexadecimal floating-point)
- Abramowitz & Stegun reference tables
- NIST Digital Library of Mathematical Functions
- Decades of production use in weapons simulation, aerospace, and scientific computing

We port only routines that pass a four-level test framework:

| Level | Question |
|-------|----------|
| **L1** | Does the code work? |
| **L2** | Does the output match mathematical truth? |
| **L3** | Does the output match IBM 360 historical values? |
| **L4** | Does the output survive hostile compiler flags? |

If a routine doesn't pass all four levels in the original Fortran, it doesn't get ported.

## Current Status

| Module | Routines | Ported | Tested |
|--------|----------|--------|--------|
| BLAS | DAXPY, DSCAL, DCOPY, DSWAP, DROT, DROTG | In Progress | — |
| MINPACK | DENORM, ENORM | In Progress | — |
| LINPACK | DGEFA, DGESL, DPOFA, DPOSL | Planned | — |
| Special Functions | Bessel, Gamma, Airy, Elliptic | Planned | — |

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
3. The code produces identical output to IBM 360 mainframes from the 1980s
4. The implementation survives adversarial compiler optimisations

We do not port routines with known accuracy failures. For example, SLATEC's Bessel J functions fail catastrophically for large arguments (J_0(500) returns the wrong sign). Those routines are documented but not ported.

See [slatec-modern](https://github.com/Zaneham/slatec-modern) for the four-level test framework and deviation analysis.

## License

Apache License 2.0. See [LICENSE](LICENSE) for details.

The original SLATEC library is public domain.
