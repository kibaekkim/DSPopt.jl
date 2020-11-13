"""
Distributionally robust extension
"""

abstract type AmbiguitySet end
abstract type WassersteinSet <: AmbiguitySet end

mutable struct DroData{WassersteinSet}
    lp_norm::Float64
    ϵ::Float64
end

"""
    set

Set the parameters for Wasserstein ambiguity set

# Arguments
- `dsp`: DSPProblem struct
- `p`: norm of the Wasserstein distance
- `ϵ`: the size of the Wasserstein ball
"""
function set(::Type{WassersteinSet}, p, ϵ)
    dspenv.dro = DroData{WassersteinSet}(p, ϵ)
end

function loadDroData(dsp::DSPProblem, dro::Nothing) 
end

function loadDroData(dsp::DSPProblem, dro::DroData{WassersteinSet}) 
    @dsp_ccall("setWassersteinAmbiguitySet", Cvoid, 
        (Ptr{Cvoid}, Cdouble, Cdouble), 
        dsp.p, dro.lp_norm, dro.ϵ)
end
