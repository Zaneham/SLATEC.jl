module SLATEC

# BLAS Level 1
export daxpy!, dscal!, dcopy!, dswap!, drot!, drotg

# MINPACK
export denorm, enorm

# To be implemented:
# - LINPACK (DGEFA, DGESL, DPOFA, DPOSL)
# - Special functions (Bessel, Gamma, etc.)

include("blas.jl")
include("minpack.jl")

end # module
