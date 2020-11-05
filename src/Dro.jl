"""
Distributionally robust extension
"""

abstract type AmbiguitySet end
abstract type WassersteinSet <: AmbiguitySet end

"""
    set

Set the parameters for Wasserstein ambiguity set

# Arguments
- `p`: norm of the Wasserstein distance
- `ϵ`: the size of the Wasserstein ball
"""
function set(::Type{WassersteinSet}, p::Real, ϵ::Float64) end
