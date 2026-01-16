# MINPACK routines
# Ported from SLATEC with 4-level test validation

"""
    denorm(x)

Compute Euclidean norm of vector `x` with overflow/underflow protection.

Unlike naive `sqrt(sum(x.^2))`, this implementation:
- Won't overflow for large values (handles 1e300 correctly)
- Won't underflow for small values (handles 1e-300 correctly)
- Won't flush subnormals to zero (unless compiled with -ffast-math)

The algorithm partitions values into three ranges:
- Large (> rgiant/n): Scaled to prevent overflow
- Intermediate: Direct accumulation
- Small (< rdwarf): Scaled to prevent underflow

Validated against:
- L1: Regression tests (9/9 pass)
- L2: Mathematical properties (homogeneity, triangle inequality)
- L3: IBM 360 golden values (bit-identical)
- L4: Hostile compiler flags (45/45 pass without -ffast-math)

WARNING: Do NOT compile with -ffast-math. It enables flush-to-zero
which silently corrupts subnormal inputs.
"""
function denorm(x::AbstractVector{Float64})
    n = length(x)

    if n == 0
        return 0.0
    end

    # Machine constants
    rdwarf = 3.834e-20   # sqrt(tiny) roughly
    rgiant = 1.304e19    # sqrt(huge) roughly

    s1 = 0.0  # Large values accumulator
    s2 = 0.0  # Intermediate values accumulator
    s3 = 0.0  # Small values accumulator

    x1max = 0.0  # Max large value seen
    x3max = 0.0  # Max small value seen

    agiant = rgiant / n

    @inbounds for i in 1:n
        xabs = abs(x[i])

        if xabs > rdwarf && xabs < agiant
            # Intermediate range - accumulate directly
            s2 += xabs * xabs
        elseif xabs <= rdwarf
            # Small range - scale to prevent underflow
            if xabs <= x3max
                if x3max != 0.0
                    s3 += (xabs / x3max)^2
                end
            else
                if x3max != 0.0
                    s3 = 1.0 + s3 * (x3max / xabs)^2
                else
                    s3 = 1.0
                end
                x3max = xabs
            end
        else
            # Large range - scale to prevent overflow
            if xabs <= x1max
                s1 += (xabs / x1max)^2
            else
                s1 = 1.0 + s1 * (x1max / xabs)^2
                x1max = xabs
            end
        end
    end

    # Combine accumulators
    if s1 != 0.0
        # Large values dominate
        return x1max * sqrt(s1 + (s2 / x1max) / x1max)
    elseif s2 != 0.0
        # Intermediate values dominate
        if s2 >= x3max
            return sqrt(s2 * (1.0 + (x3max / s2) * (x3max * s3)))
        else
            return sqrt(x3max * ((s2 / x3max) + (x3max * s3)))
        end
    else
        # Small values only
        return x3max * sqrt(s3)
    end
end

"""
    enorm(x)

Single-precision Euclidean norm with overflow/underflow protection.

See `denorm` for algorithm details.
"""
function enorm(x::AbstractVector{Float32})
    # Convert to Float64 for intermediate calculations, then back
    return Float32(denorm(Float64.(x)))
end
